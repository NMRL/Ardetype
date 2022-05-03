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

If starting from bcl:
    find sample sheet 
    qsub bcl2fastq - check for 
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