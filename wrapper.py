import os
import sys
import time
import argparse
import subprocess
from difflib import SequenceMatcher
import warnings
import pandas as pd

####################################################################
###########             path should be given             ###########
####################################################################
nextflow_path="" # nextflow binary path
host_reference_dir_path = "" # directory containing host reference sequences
####################################################################

nextflow_script_path=os.path.dirname(os.path.realpath(__file__))+"/"
scripts_path=nextflow_script_path+"/scripts"+"/"

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

def get_host_reference_path(host_):
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
        warnings.warn("There was no detected host genome. Pipeline will run without host filtering step.")
        return None
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
        
class IlluminaArguments:
    def __init__(self, prefix_, outdir_, host_, fastq_, fastq2_):
        self.platform = "illumina"
        self.prefix = prefix_
        self.outdir = outdir_
        
        self.host = None
        self.host_reference_path = None

        if host_:
            host_reference_path = get_host_reference_path(host_)
            if host_reference_path:
                self.host = host_
                self.host_reference_path = host_reference_path
        
        self.fastq = fastq_       
        self.fastq2 = fastq2_

class NanoporeArguments:
    def __init__(self, prefix_, outdir_, host_, fastq_):
        self.platform = "nanopore"
        self.prefix = prefix_
        self.outdir = outdir_

        self.host = None
        self.host_reference_path = None

        if host_:
            host_reference_path = get_host_reference_path(host_)
            if host_reference_path:
                self.host = host_
                self.host_reference_path = host_reference_path

        self.fastq = fastq_

class PseudoArgs:
    def __init__(self, prefix_, outdir_, host_, fastq_, fastq2_= None):
        self.prefix = prefix_
        self.outdir = outdir_
        self.host = host_
        self.fastq = fastq_
        self.fastq2 = fastq2_

def parse_one_input(platform, args_):

    if platform == "illumina":
        if args_.outdir:
            arguments = IlluminaArguments(args_.prefix, args_.outdir, args_.host, args_.fastq, args_.fastq2)
        else: # if outdir not provided, prefix will replace it
            arguments = IlluminaArguments(args_.prefix, args_.prefix, args_.host, args_.fastq, args_.fastq2)
    
    elif platform == "nanopore":
        if args_.outdir:
            arguments = NanoporeArguments(args_.prefix, args_.outdir, args_.host, args_.fastq)
        else: # if outdir not provided, prefix will replace it
            arguments = NanoporeArguments(args_.prefix, args_.prefix, args_.host, args_.fastq)
    
    else:
        sys.exit("This cannot happen. Ask to the developer.")

    return arguments

def parse_file_input(platform_, file_input_):

    df = pd.read_csv(file_input_, header=0, sep="\t")

    df_colnames = list(df)
    allowed_colnames = ["prefix", "outdir", "host", "fastq", "fastq2"]
    for colname in df_colnames:
        if colname not in allowed_colnames:
            sys.exit("Wrong column names are used. Check your input file.\nAllowed names are " + ",".join(allowed_colnames) + ".")

    if "prefix" not in df_colnames or "fastq" not in df_colnames:
        "At least prefix and fastq columns are required"
    
    if platform_ == "illumina" and "fastq2" not in df_colnames:
        sys.exit("Pipeline only supports paired-end Illumina sequencing. Input file doesn't have fastq2 column.")

    if platform_ == "nanopore" and "fastq2" in df_colnames:
        warnings.warn("fastq2 will not be used for Nanopore sequencing analysis.")

    if "host" not in df_colnames:
        df["host"] = ""
    if "outdir" not in df_colnames:
        df["outdir"] = ""

    arguments_list = []

    for index, row in df.iterrows():
        if platform_ == "illumina":
            new_args = PseudoArgs(row["prefix"], row["outdir"], row["host"], row["fastq"], row["fastq2"])

        else:
            new_args = PseudoArgs(row["prefix"], row["outdir"], row["host"], row["fastq"])
            
        arguments = parse_one_input(platform_, new_args)
        valid_check(arguments)
        arguments_list.append(arguments)

    return arguments_list


def valid_check(Arguments_):

    if not Arguments_.fastq :
        raise FastqError
    else:
        if not os.path.exists(Arguments_.fastq.strip('\"')):
            raise FastqError
    
    if Arguments_.platform == "illumina":
        if not Arguments_.fastq2:
            raise FastqError
        else:
            if not os.path.exists(Arguments_.fastq2.strip('\"')):
                raise FastqError

    if not Arguments_.prefix:
        raise PrefixError

def parse_arguments_to_cmd(Arguments_):
    cmd = ""
    cmd = nextflow_path + "nextflow run " + nextflow_script_path 
    if Arguments_.platform == "illumina":
        cmd = cmd + "main_illumina.nf "
    else:
        cmd = cmd + "main_nanopore.nf "

    cmd += "--prefix " + Arguments_.prefix + " "
    cmd += "--outdir " + Arguments_.outdir + " "
    if(Arguments_.host):
        cmd += "--host " + "\"" + Arguments_.host + "\" "
        cmd += "--host_reference_path " + Arguments_.host_reference_path + " "
    cmd += "--fastq " + Arguments_.fastq + " "

    if Arguments_.platform == "illumina":
        cmd += "--fastq2 " + Arguments_.fastq2

    
    return cmd
    
parser = argparse.ArgumentParser(description= 'Wrapper script for running pipeline')

parser.add_argument('platform', choices=['nanopore', 'illumina', 'host'])

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

parser.add_argument('--resume', action='store_true',
        help='Whether nextflow resume the previous run')

parser.add_argument('--test', action='store_true',
        help='Whether it is test run')

args = parser.parse_args()

platform = args.platform
if platform == "host":
    print('\n'.join(host_reference_dict.keys()))
    print('\nAbove are available host genomes.\nYou need to add to wrapper.py if you want to add new host genome.')
    exit()

file_input = args.inputs_from_file
if file_input:
    arguments_list = parse_file_input(platform, file_input)
else:
    if args.prefix is None or args.fastq is None:
        parser.error("At least --prefix and --fastq are required")
    if args.platform == "illumina" and not args.fastq2:
        sys.exit("Pipeline only supports paired-end Illumina sequencing. fastq2 is not given.")
    if args.platform == "nanopore" and args.fastq2:
        warnings.warn("fastq2 will not be used for Nanopore sequencing analysis.")
    arguments_list = [ parse_one_input(platform, args) ]

cmd_list = [ parse_arguments_to_cmd(x) for x in arguments_list ]

test = args.test
resume = args.resume

for i in range(len(cmd_list)) :
    cmd = cmd_list[i]
    cmd += " --nextflow_script_path " + nextflow_script_path
    t = time.localtime()
    curr_time = time.strftime("%Y_%m_%d_%H_%M_%S", t)
    cmd += "" if test else (" -with-report " + curr_time  + "_nextflow_running_report.html" + " -with-trace " + curr_time + "_nextflow_trace_report.txt" + " -with-timeline " + curr_time + "_nextflow_timeline_report.html")
    cmd += " -resume " if resume else ""
    subprocess.run(cmd, shell=True, check=True)

    prefix = arguments_list[i].prefix
    report_cmd = "Rscript " + scripts_path + "6_report_with_nozzle_" + platform + ".R --prefix " + prefix
    subprocess.run(report_cmd, shell=True, check=True)
    