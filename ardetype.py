"""
This is a wrapper script of ARDETYPE(?) pipeline.
Date: 2022-04-27
Version: 0.0
"""
import sys, argparse, yaml, subprocess

###Logic
"""
Pipeline can start from:
    a. bcl files - run all
    b. fastq files - run all, except demultiplexing module in core
    c. contigs.fasta - run only shell, relevant tip modules and shape
Pipeline output structure:
    bact_output/input_folder_name/bact_core/
        sample_id_raw fastq(excluding undetermined)
        sample_id_host_filtered fastq
        sample_id_contigs.fasta
        sample_id_contig_based_taxonomy
        sample_list (AQUAMIS format + majority genus from reads for each sample (kraken2))

from bcl:
    input - path to run folder
    check for sample_sheet
        if missing - exit with message indicating missing input
    qsub bcl2fastq
        if output files are missing after job is finished - exit with message indicating job error
from fastq:
    input - path to folder with raw fastq files
    generate sample_id_list
    qsub bact_core
        check output for each sample in id list >

        
    
"""

###Arguments template (https://docs.python.org/3/library/argparse.html)
# parser = argparse.ArgumentParser(description='This is a wrapper script of ARDETYPE(?) pipeline.')
# parser.add_argument('-p1', '--placeholder1', metavar='\b', help = 'Placeholder argument 1 - not required', default=3, required=False)
# req_arg_grp = parser.add_argument_group('required arguments') #to display argument under required header in help message
# req_arg_grp.add_argument('-p2', '--placeholder2', metavar='\b', help = 'Placeholder argument 2 - required', default=None, required=True)
# if len(sys.argv)==1: #if no command-line arguments provided - display help and stop script excecution
#     parser.print_help(sys.stderr)
#     sys.exit(1)
# args = parser.parse_args()

###Template to read config yaml file into dict
# with open("yaml_path", 'r') as yaml_handle:
#     config_dict=yaml.safe_load(yaml_handle)
#     print(config_dict)


###Templates to run shell command using subprocess
# subprocess.check_call(['qsub', '-F', f'{read_1} {read_2} {sample_id} {database}', 'job.sh'])
# subprocess.call(f'bash create_sampleSheet.sh --mode illumina --fastxDir {config["target_dir"]} --outDir {config["target_dir"]}'.split(' '))