import os
import time
import argparse
import subprocess
from difflib import SequenceMatcher
import warnings

nextflow_path="/home/molecularvirology/miniconda2/envs/vdp_lrs/bin/"
nextflow_script_path=os.path.dirname(os.path.realpath(__file__))+"/"

host_ref_dir_path = "/media/molecularvirology/b6d973a5-06b6-4aae-9d77-bf28064b405e/Kijin/host_reference"
host_dict = { 
              'rhinolophus ferrumequinum' : host_ref_dir_path + "/GCF_004115265.1_mRhiFer1_v1.p_genomic.fna.gz",
              'pygoscelis antarctica' : host_ref_dir_path + "/GCA_010078415.1_BGI_Pant.V1_genomic.fna.gz",
              'haemaphysalis longicornis' : host_ref_dir_path + "/GCA_010078415.1_BGI_Pant.V1_genomic.fna.gz",
              'homo sapiens' : host_ref_dir_path + "/GCF_000001405.39_GRCh38.p13_genomic.fna.gz",
              'apodemus agrarius' : host_ref_dir_path + "/NC_016428.1_Apodemus_agrarius_mitochondrion.fasta",
              'tscherskia triton' : host_ref_dir_path + "/NC_013068.1_Tscherskia_triton_mitochondrion.fasta",
              'rattus norvegicus' : host_ref_dir_path + "/GCF_015227675.2_mRatBN7.2_genomic.fna.gz"
            }

def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exist!" % arg)
    return arg

class PrefixError(ValueError):
    def __str__(self):
        return "There is no prefix given. Please check the argument '--prefix'!"

class FastqError(ValueError):
    def __str__(self):
        return "Any required fastq is not given or does not exist. Please check the argument '--fastq' or '--fastq2'!"

class HostNameError(ValueError):
    def __str__(self):
        return "invalid host name. Please check the argument '--host'!"

def getHostRefPath(host):
    if host in host_dict:
        return host_dict[host]
    else:
        closest_cands = dict(filter(lambda elem: SequenceMatcher(None, host, elem[0]).ratio() > 0.75, host_dict.items()))
        for cand in closest_cands :
            ans = input("Did you mean " + cand + "?(y/n)")
            
            while ans != "y" and ans != "n" :
                ans = input("Please type either 'y' or 'n':")
        
            if ans == "y" :
                host = cand
            elif ans == "n":
                continue
    if not host in host_dict :
        print("There was no detected host name. Pipeline will run without host filtering.")
    else:
        return host_dict[host]


def parseArgsToArgsDictList(args):
    
    args_dict_list = []

    file_input = args.inputs_from_file
    # get inputs from file input
    if file_input:
        with open(file_input) as file:
            lines = file.readlines()
            header = [ x.strip() for x in lines[0].split('\t') ]
            for line in lines[1:]:
                if line == "\n":
                    continue
                li = [ x.strip() for x in line.split('\t') ]
                args_dict = {}
                args_dict["sequencing_type"] = args.sequencing_type
                host = ""
                for index in range(len(header)) :
                    col_name = header[index]
                    if col_name == "host":
                        host = li[index]
                    else:
                        args_dict[col_name] = '"' + li[index] + '"'

                if host:
                    args_dict["host"] = '"' + host.lower() + '"'
                    args_dict["host_ref_path"] = getHostRefPath(host.lower())
                # if host not given, host filtering is not executed

                if not "outdir" in args_dict and "prefix" in args_dict:
                    args_dict["outdir"] = args_dict["prefix"]
                                
                args_dict_list.append(args_dict)

    # one input analysis
    else:
        args_dict = {}
        args_dict["sequencing_type"] = args.sequencing_type
        args_dict["fastq"] = args.fastq
        if args_dict["sequencing_type"] == "illumina":
            args_dict["fastq2"] = args.fastq2
        args_dict["prefix"] = args.prefix
        if args.host:
            args_dict["host"] = '"' + args.host.lower() + '"'
            args_dict["host_ref_path"] = getHostRefPath(args.host.lower())
        # if host not given, host filtering is not executed
        if not args.outdir:
            args_dict["outdir"] = args_dict["prefix"]
        args_dict_list.append(args_dict)

    return args_dict_list

def args_dict_list_valid_check(args_dict_list):
    for args_dict in args_dict_list:
        if "fastq" not in args_dict :
            raise FastqError
        else:
            if not args_dict["fastq"] or not os.path.exists(args_dict["fastq"].strip('\"')):
                raise FastqError
        if args_dict["sequencing_type"] == "illumina":
            if "fastq2" not in args_dict:
                raise FastqError
            else:
                if not args_dict["fastq2"] or not os.path.exists(args_dict["fastq2"].strip('\"')):
                    raise FastqError
        else:
            if "fastq2" in args_dict:
                warnings.warn("fastq2 will not be used for nanoporeseq analysis.")

        if "prefix" not in args_dict or not args_dict["prefix"]:
            raise PrefixError

def parseArgsDictToCmd(args_dict):
    cmd = ""
    cmd = nextflow_path + "nextflow run " + nextflow_script_path 
    if args_dict["sequencing_type"] == "illumina":
        cmd = cmd + "main_nanopore.nf "
    else:
        cmd = cmd + "main_illumina.nf "
    for key in (k for k in args_dict.keys() if k != "seq_type"):
        cmd += "--" + key + " " + args_dict[key] + " "
    return cmd

def parseArgsDictListToCmdList(args_dict_list):
    cmd_list = []
    for args_dict in args_dict_list:
        cmd_list.append(parseArgsDictToCmd(args_dict))
    return cmd_list

def parseArgsToCmdList(args):
    args_dict_list = parseArgsToArgsDictList(args)
    args_dict_list_valid_check(args_dict_list)
    cmd_list = parseArgsDictListToCmdList(args_dict_list)
    return cmd_list
    

parser = argparse.ArgumentParser(description= 'Wrapper script for running pipeline')

parser.add_argument('sequencing_type', choices=['nanopore', 'illumina'])

parser.add_argument('--inputs_from_file', '-f', metavar='file.txt',
        help='Get inputs from this text file. Column header corresponding to the options should be provided')

parser.add_argument('--fastq', '-fq', metavar='fastq',
        help='First fastq file taken as input for the pipeline',
        type=lambda x: is_valid_file(parser, x))
    
parser.add_argument('--fastq2', '-fq2', metavar='fastq_2',
        help='Second fastq file taken as input for the pipeline in case of paired end sequencing',
        type=lambda x: is_valid_file(parser, x)) 

parser.add_argument('--prefix', '-p', metavar='name',
        help='Labels for each input fastq')

parser.add_argument('--outdir', '-d', metavar='dir',
        help='Name of output directory')

parser.add_argument('--host', metavar='host',
        help='Scientific name of host organism'),

parser.add_argument('--reference', metavar='ref',
        help='The location of the directory where reference virus sequences are located')

parser.add_argument('--background', '-bg', action='store_true',
        help='Whether nextflow run in background or not')

parser.add_argument('--resume', action='store_true',
        help='Whether nextflow resume the previous run')

parser.add_argument('--test', action='store_true',
        help='Whether it is test run')

args = parser.parse_args()

cmd_list = parseArgsToCmdList(args)

test = args.test
reference = args.reference
background = args.background
resume = args.resume

for cmd in cmd_list :
    t = time.localtime()
    curr_time = time.strftime("%Y_%m_%d_%H_%M_%S", t)
    cmd += "" if test else (" -with-report " + curr_time  + "_nextflow_running_report.html" + " -with-trace " + curr_time + "_nextflow_trace_report.txt" + " -with-timeline " + curr_time + "_nextflow_timeline_report.html")
    cmd += (" --reference_vir_path \"" + reference + "\"") if reference else ""
    cmd += " -bg " if background else ""
    cmd += " -resume " if resume else ""
    subprocess.run(cmd, shell=True, check=True)    
    #process = subprocess.Popen(cmd, shell=True)
    #time.sleep(1)


