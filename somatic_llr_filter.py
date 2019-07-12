import argparse
import sys
import os
import re
import vcfpy
import tempfile
import csv
import math
from collections import OrderedDict

def define_parser():
    parser = argparse.ArgumentParser('somatic-llr-filter')
    parser.add_argument(
        "input_vcf",
        help="A VCF file with at least two samples (tumor and normal) and readcount information"
    )
    parser.add_argument(
        "output_vcf",
        help="Path to write the output VCF file"
    )
    parser.add_argument(
        "--sequence-error-rate",
        help="expected sequencing error rate (range 0 to 1) - default 0.005",
        type=float,
        default=0.005
    )
    parser.add_argument(
        "--tumor-purity",
        help="tumor purity (range 0 to 1) - default 1",
        type=float,
        default=1
    )
    parser.add_argument(
        "--normal-contamination-rate",
        help="normal contamination rate (range 0 to 1) - default 0",
        type=float,
        default=0
    )
    parser.add_argument(
        "--allele-depth-field",
        help="field corresponding to allele depth - default AD",
        default="AD"
    )
    parser.add_argument(
        "--site-depth-field",
        help="field corresponding to site depth - default DP",
        default="DP"
    )
    parser.add_argument(
        "--tumor-sample-name",
        help="name of the vcf sample corresponding to the tumor - default TUMOR",
        default="TUMOR"
    )
    parser.add_argument(
        "--normal-sample-name",
        help="name of the vcf sample corresponding to the normal - default NORMAL",
        default="NORMAL"
    )
    parser.add_argument(
        "--llr-field",
        help="name of the vcf info field in which to store the log-likelihood ratio value",
        default="LLR"
    )
    parser.add_argument(
        "--somatic-field",
        help="name of the vcf info field in which to store the classification call from this script (somatic = 0 or 1)",
        default="SOMATIC"
    )
    parser.add_argument(
        "--llr-threshold",
        type=float,
        help="if set, variants that are not somatic or have an LLR value below this threshold will have their FILTER field set appropriately"
    )
    parser.add_argument(
        "--filter-field",
        help="if --llr-threshold is given, then failing variants will have this string added to their FILTER field",
        default="SOMATIC_LLR"
    )
   
    parser.add_argument('-w', "--overwrite", action='store_true',
        help="by default, this tool will raise an exception if the LLR or SOMATIC fields already exist in the VCF. This flag allows existing fields to be overwritten."
    )
    return parser


def create_vcf_reader(args):
    vcf_reader = vcfpy.Reader.from_path(args.input_vcf)
    #do some sanity checking
    sample_names = vcf_reader.header.samples.names
    is_multi_sample = len(sample_names) > 1
    if(not is_multi_sample):
        vcf_reader.close()
        raise Exception("A multisample VCF with both tumor and normal data is required")

    if not args.tumor_sample_name in sample_names:
        raise Exception("Could not find tumor sample name {} in sample names".format(args.tumor_sample_name))
    if not args.normal_sample_name in sample_names:
        raise Exception("Could not find normal sample name {} in sample names".format(args.normal_sample_name))

    #check for needed format fields
    if not args.allele_depth_field in vcf_reader.header.format_ids():
        vcf_reader.close()
        raise Exception("No " + args.allele_depth_field + " format field found. Annotate your VCF with readcounts first")
    if not args.site_depth_field in vcf_reader.header.format_ids():
        vcf_reader.close()
        raise Exception("No {} format field found. Annotate your VCF with readcounts first".format(args.site_depth_field))

    return vcf_reader


def create_vcf_writer(args, vcf_reader):
    output_file = args.output_vcf

    new_header = vcf_reader.header.copy()

    #check/add llr field in header
    if args.llr_field in vcf_reader.header.info_ids():
        if args.overwrite:  #verify compatibility
            if not vcf_reader.header.get_info_field_info(args.llr_field).type == "Float":
                vcf_reader.close()
                raise Exception("{} field to be overwritten must be of type 'Float'. Either modify this, or choose a new {} field".format(args.llr_field,args.llr_field))
            if not vcf_reader.header.get_info_field_info(args.llr_field).number == 1:
                vcf_reader.close()
                raise Exception("{} field to be overwritten must have Number '1'. Either modify this, or choose a new {} field".format(args.llr_field,args.llr_field))
        else: 
            vcf_reader.close()
            raise Exception("INFO already contains a {} field. Choose a different label, or use the --overwrite flag to retain this field description and overwrite values".format(args.llr_field))

    else:
        od = OrderedDict([('ID', args.llr_field), ('Number', '1'), ('Type', 'Float'), ('Description', 'log-likelihood ratio for the binomial filter call')])
        new_header.add_info_line(od)

    #check/add somatic field in header
    if args.somatic_field in vcf_reader.header.info_ids():
        if args.overwrite:
            if not vcf_reader.header.get_info_field_info(args.somatic_field).type == "Flag":
                vcf_reader.close()
                raise Exception("{} field to be overwritten must be of type 'Flag'. Either modify this, or choose a new {} field".format(args.somatic_field,args.somatic_field))
            if not vcf_reader.header.get_info_field_info(args.somatic_field).number == 0:
                vcf_reader.close()
                raise Exception("{} field to be overwritten must have Number '0'. Either modify this, or choose a new {} field".format(args.somatic_field,args.somatic_field))
        else:
            vcf_reader.close()
            raise Exception("INFO already contains a {} field. Choose a different label, or use the --overwrite flag to retain this field description and overwrite values".format(args.somatic_field))
    else:
        od = OrderedDict([('ID', args.somatic_field), ('Number', '0'), ('Type', 'Flag'), ('Description', 'Is a somatic mutation')])
        new_header.add_info_line(od)

    #check/add FILTER field in header
    if args.llr_threshold is not None: #filtering may not even be specified
        if args.filter_field in vcf_reader.header.filter_ids():
            if not args.overwrite:
                vcf_reader.close()
                raise Exception("FILTER {} already exists. Choose a different --filter-field, or use the --overwrite flag to retain this filter description and overwrite values".format(args.filter_field))
        else:
            od = OrderedDict([('ID', args.filter_field), ('Description', 'Is a somatic mutation with LLR greater than {}'.format(args.llr_threshold))])
            new_header.add_filter_line(od)


    return vcfpy.Writer.from_path(output_file, new_header)


#for each call possibility (ref, germ het, germ hom, etc)
#calculate the llr here
def calc_llr(normal_expect, tumor_expect, normal_ref, normal_var, tumor_ref, tumor_var):
    return(normal_ref * math.log(1-normal_expect) +
           normal_var * math.log(normal_expect) +
           tumor_ref * math.log(1-tumor_expect) +
           tumor_var * math.log(tumor_expect))

def get_llr(call, normal_ref, normal_var, tumor_ref, tumor_var, error_rate, heterozygous_expect,
            homozygous_expect, tumor_freq, tumor_purity, normal_contamination_rate, error_expect):
    if call == "Germline_het":
        return(calc_llr(heterozygous_expect, heterozygous_expect,
                       normal_ref, normal_var, tumor_ref, tumor_var))
    elif call == "Germline_hom":
        return(calc_llr(homozygous_expect, homozygous_expect,
                    normal_ref, normal_var, tumor_ref, tumor_var))
    elif call == "LOH_ref":
        return(calc_llr(heterozygous_expect, error_rate,
                       normal_ref, normal_var, tumor_ref, tumor_var))
    elif call == "LOH_variant":
        return(calc_llr(heterozygous_expect, homozygous_expect,
                       normal_ref, normal_var, tumor_ref, tumor_var))
    elif call == "NotSomatic":
        return(calc_llr(error_expect, error_expect,
                       normal_ref, normal_var, tumor_ref, tumor_var))
    elif call == "Reference":
        return(calc_llr(error_rate, error_rate,
                       normal_ref, normal_var, tumor_ref, tumor_var))
    elif call == "Somatic":
        return(calc_llr(error_rate + (tumor_freq / tumor_purity * normal_contamination_rate), tumor_freq,
                       normal_ref, normal_var, tumor_ref, tumor_var))
    raise Exception('Call "' + call + '" is not a recognized type')


def make_call(normal_ref, normal_var, tumor_ref, tumor_var, error_rate, heterozygous_expect,
              homozygous_expect, tumor_freq, tumor_purity, normal_contamination_rate, error_expect):
    #Want to test several models and pick the most likely one
    #Reference: normal expectation is 0.001 and tumor expectation is 0.01
    #Germline Het: normal expectation is 0.5 and tumor expectation is 0.5
    #Somatic Het: normal expectation is 0.001 and tumor expectation in 0.5
    #Germline Homozygote: normal expectation is 0.999 and tumor expectation is 0.999
    #LOH (variant): germline is 0.5 and tumor is 0.999
    #LOH (ref): germline is 0.5 and tumor is 0.001

    marginal_probability = None
    max_llr = None
    max2_llr = None
    max_call = None
    max2_call = None
    llr = 0

    if(tumor_var == 0): #special case that chokes some of the log code.  LLR is always zero if we have no supp reads
        return([0,"Reference"])

    call_types = ["Germline_het", "Germline_hom", "LOH_ref", "LOH_variant", "NotSomatic", "Reference", "Somatic"];
    for call in call_types:
        llr = get_llr(call, normal_ref, normal_var, tumor_ref, tumor_var, error_rate, heterozygous_expect,
                      homozygous_expect, tumor_freq, tumor_purity, normal_contamination_rate, error_expect)
        if marginal_probability is None:
            marginal_probability=llr
        elif llr > marginal_probability:
            marginal_probability = llr + math.log(1 + math.exp(marginal_probability-llr))
        else:
            marginal_probability = marginal_probability + math.log(1 + math.exp(llr-marginal_probability))

        max_llr_def = (max_llr is not None)
        max2_llr_def = (max2_llr is not None)

        if not (max_llr_def and (llr < max_llr)):
            #save the second-place, store the new best
            max2_call = max_call
            max2_llr = max_llr
            max_llr = llr
            max_call = call
        elif not (max2_llr_def and (llr < max2_llr)):
            max2_call = call
            max2_llr = llr
    #end for

    #get the final llr
    if max_llr is not None:
        llr = max_llr-max2_llr
    else:
        llr = 0

    #get the final call
    if max_call is None:
        max_call = "-"

    return (llr,max_call)



def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    #these aren't going to change per site
    homozygous_expect = 1 - args.sequence_error_rate;
    heterozygous_expect = 0.5;

    vcf_reader  = create_vcf_reader(args)
    vcf_writer = create_vcf_writer(args, vcf_reader)

    for entry in vcf_reader:
        #collect the needed info
        ref = entry.REF
        alts = entry.ALT

        #this code will mostly handle multiple alleles, as written, but the issue is that there is
        #no way to set the INFO (SOMATIC) field appropriately when there are multiple alleles per line,
        # so we're going to restrict it to decomposed VCFs - one variant per line
        if len(alts) > 1:
             raise Exception("site with multiple alleles detected. This tool requires a decomposed vcf (one allele per line) so that INFO fields can be set appropriately")

        def getFormatField(sample_name, field_name):
            if(sample_name in entry.call_for_sample and field_name in entry.call_for_sample[sample_name].data):
                return entry.call_for_sample[sample_name].data[field_name]
            return("NA")

        def missingVals(arr):
            for i in arr:
                if i == "NA":
                    return True
            return False

        ad_nrm = getFormatField(args.normal_sample_name,args.allele_depth_field)
        ad_tum = getFormatField(args.tumor_sample_name,args.allele_depth_field)
        normal_depth = getFormatField(args.normal_sample_name,args.site_depth_field)
        tumor_depth = getFormatField(args.tumor_sample_name,args.site_depth_field)
        normal_ref = ad_nrm[0]
        tumor_ref = ad_tum[0]

        call = ""
        llr = 0

        #TODO parse out per alt, retrieve calls
        for i in range(1,(len(alts)+1)):  #right now, this will only ever be one, due to above check.  Could be expanded to support multiple alleles - see above
            if(ad_nrm == "NA"):
                normal_var = "NA"
            else:
                normal_var = ad_nrm[i]
            if(ad_tum == "NA"):
                tumor_var = "NA"
            else:
                tumor_var = ad_tum[i]
            
            #if neither has any depth or vals or missing, then fail this up front
            if missingVals([normal_var,tumor_var,tumor_depth,normal_depth]) or (tumor_depth + normal_depth == 0):
                (llr,call) = (0,"Reference")
                continue
            
            #weighted average of the frequencies
            error_expect = (tumor_var + normal_var)/(tumor_depth + normal_depth)
            #dave had this line in there, but I don't know why...
            #error_expect ||= error_rate
            if error_expect == 1:
                error_expect = 1 - args.sequence_error_rate

            #handle case where depth is zero or missing
            if tumor_depth == 0:
                tumor_freq = 0
            else:
                tumor_freq = tumor_var/tumor_depth

            if tumor_freq == 0:
                tumor_freq = args.sequence_error_rate
            elif tumor_freq == 1:
                tumor_freq = tumor_freq - args.sequence_error_rate

            (llr,call) = make_call(normal_ref, normal_var, tumor_ref, tumor_var, args.sequence_error_rate, heterozygous_expect,
              homozygous_expect, tumor_freq, args.tumor_purity, args.normal_contamination_rate, error_expect)

        #Store it back in the entry before writing out.
        if call == "Somatic":
            entry.INFO[args.somatic_field] = 1
        else:
            entry.INFO.pop('SOMATIC',None)

        entry.INFO[args.llr_field] = llr

        #do filtering if specified
        if args.llr_threshold is not None:
            if (not call == "Somatic") or (llr < args.llr_threshold):
                entry.add_filter(args.filter_field);

        vcf_writer.write_record(entry)

    vcf_reader.close()
    vcf_writer.close()

if __name__ == '__main__':
    main()
