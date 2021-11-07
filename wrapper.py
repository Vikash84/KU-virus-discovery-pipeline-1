import os
import time
import argparse
import subprocess
from difflib import SequenceMatcher
import warnings

nextflow_path=""
host_reference_dir_path = ""

nextflow_script_path=os.path.dirname(os.path.realpath(__file__))+"/"

host_reference_dict = { 
#              'rhinolophus ferrumequinum' : host_reference_dir_path + "/GCF_004115265.1_mRhiFer1_v1.p_genomic.fna.gz",  #bat
#              'pygoscelis antarctica' : host_reference_dir_path + "/GCA_010078415.1_BGI_Pant.V1_genomic.fna.gz",        #penguin
#              'haemaphysalis longicornis' : host_reference_dir_path + "/GCA_010078415.1_BGI_Pant.V1_genomic.fna.gz",    #tick
#              'homo sapiens' : host_reference_dir_path + "/GCF_000001405.39_GRCh38.p13_genomic.fna.gz",                 #human
#              'apodemus agrarius' : host_reference_dir_path + "/NC_016428.1_Apodemus_agrarius_mitochondrion.fasta",     
#              'apodemus peninsulae' : host_reference_dir_path + "/NC_016060.1_Apodemus_peninsulae_mitochondrion.fasta",
#              'apodemus chejuensis' : host_reference_dir_path + "/NC_016662.1_Apodemus_chejuensis_mitochondrion.fasta",
#              'tscherskia triton' : host_reference_dir_path + "/NC_013068.1_Tscherskia_triton_mitochondrion.fasta",
#              'rattus norvegicus' : host_reference_dir_path + "/GCF_015227675.2_mRatBN7.2_genomic.fna.gz",
            }

def is_valid_file(parser_, arg_):
    if not os.path.exists(arg_):
        parser_.error("The file %s does not exist!" % arg_)
    return arg_

def get_host_ref_path(host_):
    """
    Recommend the closest host reference
    when wrong host name is given
    """
    if host_ in host_reference_dict:
        return host_reference_dict[host_]
    else:
        closest_cands = dict(filter(lambda elem: SequenceMatcher(None, host_, elem[0]).ratio() > 0.75, host_reference_dict.items()))
        for cand in closest_cands :
            ans = input("Did you mean " + cand + "?(y/n)")
            
            while ans != "y" and ans != "n" :
                ans = input("Please type either 'y' or 'n':")
        
            if ans == "y" :
                host_ = cand
            elif ans == "n":
                continue
    if host_ not in host_reference_dict :
        print("There was no detected host name. Pipeline will run without host filtering step.")
    else:
        return host_reference_dict[host_]

class PrefixError(ValueError):
    def __str__(self):
        return "There is no prefix given. Please check the argument '--prefix'!"

class FastqError(ValueError):
    def __str__(self):
        return "Any required fastq is not given or does not exist. Please check the argument '--fastq' or '--fastq2'!"

class HostNameError(ValueError):
    def __str__(self):
        return "invalid host name. Please check the argument '--host'!"

class Arguments:
    def __init__(self, args_):
        self.args_dict_list = []

        file_input = args_.inputs_from_file
        # get multiple inputs from file with -f option
        if file_input:
            with open(file_input) as file:
                lines = file.readlines()
                header = [ x.strip() for x in lines[0].split('\t') ]
                for line in lines[1:]:
                    if line == "\n":
                        continue
                    li = [ x.strip() for x in line.split('\t') ]
                    args_dict = {}
                    args_dict["sequencing_type"] = args_.sequencing_type
                    host = ""
                    for index in range(len(header)) :
                        col_name = header[index]
                        if col_name == "host":
                            host = li[index]
                        else:
                            args_dict[col_name] = '"' + li[index] + '"'

                    if host and host.lower() != "none" and get_host_ref_path(host.lower()):
                        args_dict["host"] = '"' + host.lower() + '"'
                        args_dict["host_ref_path"] = get_host_ref_path(host.lower())
                    # if host not given, host filtering is not executed

                    if not "outdir" in args_dict and "prefix" in args_dict:
                        args_dict["outdir"] = args_dict["prefix"]
                                    
                    self.args_dict_list.append(args_dict)

        # one input
        else:
            args_dict = {}
            args_dict["sequencing_type"] = args_.sequencing_type
            args_dict["fastq"] = args_.fastq
            if args_dict["sequencing_type"] == "illumina":
                args_dict["fastq2"] = args_.fastq2
            args_dict["prefix"] = args_.prefix
            if args_.host and args_.host.lower() != "none" and get_host_ref_path(args_.host.lower()):
                args_dict["host"] = '"' + args_.host.lower() + '"'
                args_dict["host_ref_path"] = get_host_ref_path(args_.host.lower())
            # if host not given, host filtering is not executed
            if not args_.outdir:
                args_dict["outdir"] = args_dict["prefix"]
            self.args_dict_list.append(args_dict)

    def valid_check(self):
        for args_dict in self.args_dict_list:
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

    def parse_args_dict_to_cmd(self, args_dict_):
        cmd = ""
        cmd = nextflow_path + "nextflow run " + nextflow_script_path 
        if args_dict_["sequencing_type"] == "illumina":
            cmd = cmd + "main_illumina.nf "
        else:
            cmd = cmd + "main_nanopore.nf "
        for key in (k for k in args_dict_.keys() if k != "seq_type"):
            if not args_dict_[key]:
                continue
            cmd += "--" + key + " " + args_dict_[key] + " "
        return cmd
    
    def parse_args_dict_list_to_cmd_list(self):
        cmd_list = []
        for args_dict in self.args_dict_list:
            cmd_list.append(self.parse_args_dict_to_cmd(args_dict))
        return cmd_list
    

parser = argparse.ArgumentParser(description= 'Wrapper script for running pipeline')

parser.add_argument('sequencing_type', choices=['nanopore', 'illumina', 'host'])

parser.add_argument('--inputs_from_file', '-f', metavar='file.txt',
        help='Get inputs from this text file. Column header corresponding to the options should be provided')

parser.add_argument('--fastq', '-fq', metavar='fastq',
        help='First fastq file taken as input for the pipeline',
        type=lambda x: is_valid_file(parser, x))
    
parser.add_argument('--fastq2', '-fq2', metavar='fastq_2',
        help='Second fastq file taken as input for the pipeline in case of paired end sequencing',
        type=lambda x: is_valid_file(parser, x)) 

parser.add_argument('--prefix', '-p', metavar='name',
        help='Labels for input')

parser.add_argument('--outdir', '-d', metavar='dir',
        help='Name of output directory')

parser.add_argument('--host', metavar='host',
        help='Scientific name of host organism. Type "wrapper.py host" to get the list of avaible host names'),

parser.add_argument('--background', '-bg', action='store_true',
        help='Whether nextflow run in background or not')

parser.add_argument('--resume', action='store_true',
        help='Whether nextflow resume the previous run')

parser.add_argument('--test', action='store_true',
        help='Whether it is test run')

args = parser.parse_args()

if args.sequencing_type == "host":
    print('\n'.join(host_reference_dict.keys()))
    print('\nAbove are available host genomes.\nYou need to add to wrapper.py if you want to add new host genome.')
    exit()

args_obj = Arguments(args)
args_obj.valid_check()

cmd_list = args_obj.parse_args_dict_list_to_cmd_list()

test = args.test
background = args.background
resume = args.resume

for cmd in cmd_list :
    t = time.localtime()
    curr_time = time.strftime("%Y_%m_%d_%H_%M_%S", t)
    cmd += "" if test else (" -with-report " + curr_time  + "_nextflow_running_report.html" + " -with-trace " + curr_time + "_nextflow_trace_report.txt" + " -with-timeline " + curr_time + "_nextflow_timeline_report.html")
    cmd += " -bg " if background else ""
    cmd += " -resume " if resume else ""
    subprocess.run(cmd, shell=True, check=True)


