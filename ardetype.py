"""
This is a wrapper script of ARDETYPE pipeline.
Date: 2022-05-06
Version: 0.0
"""
import warnings, os, sys,  argparse
from ardetype_utilities import *
warnings.simplefilter(action='ignore', category=FutureWarning)


def parse_arguments():
    """
    Parse pre-defined set of arguments from the command line, returning a namespace (object),
    that allows accessing arguments as instance variables of namespace by their full name.
    """    
    ###Parsers
    parser = argparse.ArgumentParser(description='This is a wrapper script of ARDETYPE pipeline.', formatter_class=argparse.RawTextHelpFormatter)
    req_arg_grp = parser.add_argument_group('required arguments') #to display argument under required header in help message
    
    ###generic arguments
    #####Required
    req_arg_grp.add_argument('-m', '--mode',
        metavar='\b',
        help = """Selecting mode that allows to run specific modules of the pipeline:
        all - run all modules (starting from fastq files) (not active)
        core - run only bact_core module (starting from fastq files) 
        shell - run only bact_shell module (starting from fasta file) (not active)
        shell_tip - run bact_shell and bact_tip modules (starting from fasta file) (not active)
         """,
        default=None,
        required=True)

    #####Optional
    parser.add_argument('-c', '--config', metavar='\b', help = 'Path to the config file (if not supplied, the copy of the template config_modular.yaml will be used)', default="./config_modular.yaml", required=False)
    parser.add_argument('-r', '--skip_reporting', help = 'Use this flag to skip reporting trough bact_shape module (which will run by-default with any other option)', action='store_true')
    parser.add_argument('-o', '--output_dir', metavar='\b', help = 'Path to the output directory where the results will be stored (ardetype_output/ by-default).', default="ardetype_output/", required=False)
    parser.add_argument('-s', '--install_snakemake', help = 'Use this flag to install mamba and snakemake for the current HPC user, if it is not already installed.', action='store_true')
    parser.add_argument('-j', '--module_jobs', help='Use this flag to run modules as individual jobs on HPC (without submitting subjobs to computational nodes)', action='store_true')
    ###bact_core arguments

    #####Required
    req_arg_grp.add_argument('-f', '--fastq', metavar='\b', help = 'Path to directory that contains fastq files to be analysed (all files in subdirectories are included).', default=None, required=True)
    #####Optional

    ###If no command-line arguments provided - display help and stop script excecution
    if len(sys.argv)==1: 
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = parse_arguments()
    if args.install_snakemake:
        install_snakemake()
    if args.mode == "core":
        file_list = parse_folder(args.fastq,'.fastq.gz')
        try:
            fastq_formats = "(_S[0-9]*_R[1,2]_001.fastq.gz|_[1,2].fastq.gz|_S[0-9]*_R[1,2]_001_unclassified_out)"
            sample_sheet = create_sample_sheet(file_list, fastq_formats, mode=0)
            os.system(f"mkdir -p {args.output_dir}")
            sample_sheet.to_csv(f"{args.output_dir}sample_sheet.csv", header=True, index=False)
        except AssertionError as msg:
            print(f"Sample sheet generation error: {msg}")

        target_list = ['sample_sheet.csv']
        template_list = [
            "_contigs.fasta",
            "_bact_reads_classified_1.fastq.gz", 
            "_bact_reads_classified_2.fastq.gz",
            "_bact_reads_unclassified_1.fastq.gz",
            "_bact_reads_unclassified_2.fastq.gz",
            "_kraken2_contigs_report.txt",
            "_kraken2_host_filtering_report.txt"
        ]
        [target_list.append(f'{args.output_dir}{id}{tmpl}') for id in sample_sheet['sample_id'] for tmpl in template_list]
        target_list.remove("sample_sheet.csv")
        config_file = read_config(args.config)
        edit_config(config_file, "core_target_files", target_list)
        edit_config(config_file, "output_directory", args.output_dir)
        try:
            write_config(config_file, f'{args.output_dir}config_core.yaml')
        except AssertionError as msg:
            print(f"Configuration file manipulation error: {msg}")

        config_file_path = f'{os.path.abspath(args.output_dir)}/config_core.yaml'
        cluster_config_path = 'cluster.yaml'

        if args.module_jobs:
            job_id = submit_module_job('core',config_file_path, args.output_dir)
            check_job_completion(job_id,"bact_core",sleeping_time=5,output_dir=args.output_dir)
        else:
            run_module_cluster("core",config_file_path, cluster_config_path, 12)

        check_dict = check_module_output(file_list=target_list)
        id_check_dict = {id:"" for id in sample_sheet['sample_id']}

        for file in check_dict:
            split = file.split("/",1)[1]
            id_check_dict[split.split("_",1)[0]] += f"|{split}:{check_dict[file]}"
        
        sample_sheet = edit_sample_sheet(sample_sheet, id_check_dict, "check_note")
        sample_sheet.to_csv(f"{args.output_dir}sample_sheet.csv", header=True, index=False)
    else:
        print('Sorry, other options not supported yet.')
