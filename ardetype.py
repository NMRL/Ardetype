"""
This is a wrapper script of ARDETYPE(?) pipeline.
Date: 2022-04-27
Version: 0.0
"""
import os, sys, re, argparse, yaml, subprocess, pandas as pd

###Architecture
"""
Pipeline can start from:
    a. fastq files - run all, except demultiplexing module in core, run all or just perform assembly
    b. contigs.fasta - run only shell, relevant tip modules and shape
        shell + tip(+shape) - downstream agnostic + specific (option to report)
        shell(+shape) - downstream agnostic (option to report)
        tip(+shape) - downstream specific (option to report)
    
    *path to config file should be supplied as required argument
        template config file will be stored in github repository
        to use default settings - make a copy of the template (the file will be altered by the script)
            edit configurations in the copy of the template if customization required
            !As many options for as many tools should be accessible from config file

Pipeline output structure:
    bact_output/input_folder_name/bact_core/
        sample_id_raw fastq(excluding undetermined)
        sample_id_host_filtered fastq
        sample_id_contigs.fasta
        sample_id_contig_based_taxonomy
        sample_list (AQUAMIS format + majority genus from reads for each sample (kraken2))
        bact_shell/
            amr_pp/
            resfinder/
            kraken2_contigs/
            quast/
            mob_suite/
            rgi/
            mlst/
        bact_tip/
            based on kraken2_contigs output
            specific molecular typing for detected bacterial species
        bact_shape/
            hamronization/
            ...

from fastq:
    input - path to folder with raw fastq files
    generate sample_id_list (use aquamis script for sample id generation from different fastq formats)    qsub bact_core
        check output for each sample in id list (add status and check_note column to sample_list to indicate if check is succesful and add more information if it is not): 
            status:
                fail - cannot continue (check_note: missing <file_name> from <step_name> conducted by <module name>)
                warning - can continue, except for some steps downstream (check_note: missing <file_name> from <step_name> conducted by <module name>, affecting <affected_step_name_1>... steps)
                ok - can continue (check_note: empty)
        sample_id_raw fastq (excluding undetermined) - ok/fail
        sample_id_host_filtered fastq - ok/fail
        sample_id_contigs.fasta - ok/fail
        sample_id_contig_based_taxonomy - ok/warning
        sample_list (AQUAMIS format + majority genus from reads for each sample (kraken2)) - ok/fail
    bact_core itself will generate template string that contains all wildcards to be used by the bact_core rules
        example: '{sample_id_pattern}-{reference_sequence_pattern}-scaffolds/{sample_id_pattern}-{reference_sequence_pattern}-ragtag.scaffold.fasta' 

from contigs.fasta:
    input - path to contig file in fasta format named in format <sample_id>_contigs.fasta 
    OR 
    path to folder, containing multiple files named in format<sample_id>_contigs.fasta
    generate sample_id_list
        check file name format
            print warnings for files that does not match specified format
            if no file match format - exit with corresponding error message
            else create aquamis format table
                three columns with sample id and paths to 1st and 2nd fastq file:
                    sample	fq1	fq2
    bact_shell and other modules themselves will generate template strings that contains all wildcards to be used by the bact_shell, bact_tip & bact_shape rules
        example: '{sample_id_pattern}-{reference_sequence_pattern}-scaffolds/{sample_id_pattern}-{reference_sequence_pattern}-ragtag.scaffold.fasta' 
    if shell: qsub shell
    if tip: qsub tip
    check results
        if results contain minimum reportable amount of data and reporting option selected
            qsub shape
            check results
"""

def parse_arguments():
    """
    Parse pre-defined set of arguments from the command line, returning a namespace (object),
    that allows accessing arguments as instance variables of namespace by their full name.
    """    
    #Useful links:
    ###Argparse docs (https://docs.python.org/3/library/argparse.html)
    ###String editing within argparse: https://qa.ostack.cn/qa/?qa=985730/

    ###Parsers
    parser = argparse.ArgumentParser(description='This is a wrapper script of ARDETYPE pipeline.', formatter_class=argparse.RawTextHelpFormatter)
    req_arg_grp = parser.add_argument_group('required arguments') #to display argument under required header in help message
    
    ###generic arguments
    #####Required
    req_arg_grp.add_argument('-m', '--mode',
        metavar='\b',
        help = """Selecting mode that allows to run specific modules of the pipeline:
        all - run all modules (starting from fastq files)
        core - run only bact_core module (starting from fastq files)
        shell - run only bact_shell module (starting from fasta file)
        shell_tip - run bact_shell and bact_tip modules (starting from fasta file)
         """,
        default=None,
        required=True)

    #####Optional
    parser.add_argument('-c', '--config', metavar='\b', help = 'Path to the config file (if not supplied, the copy of the template config_modular.yaml will be used)', default=None, required=False)
    parser.add_argument('-r', '--skip_reporting', help = 'Use this flag to skip reporting trough bact_shape module (which will run by-default with any other option)', action='store_true')

    ###bact_core arguments
    #####Required
    req_arg_grp.add_argument('-f', '--fastq', metavar='\b', help = 'Placeholder argument 2 - required', default=None, required=True)
    
    #####Optional

    ###If no command-line arguments provided - display help and stop script excecution
    if len(sys.argv)==1: 
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    return args


def parse_folder(folder_pth_str, file_fmt_str, substr_lst=None, regstr_lst=None):
    '''
    Given path to the folder (str) and file format (str), returns a list, 
    containing paths of all files of specified format found in folder and subfolders,
    except for files that contain patterns to exclude (regex string list) or substrings to exclude (list of strings).    
    '''
    file_list = []
    for (root,dirs,files) in os.walk(folder_pth_str, topdown=True): #get list of file paths (from parent dir & subdir)
        for name in files:
            if (substr_lst is not None and regstr_lst is not None) and (len(substr_lst) + len(substr_lst) > 0):
                if file_fmt_str in name:
                    name_series = pd.Series([name]) #convert to pandas series for faster substring search
                    if any(name_series.str.contains('|'.join(substr_lst+regstr_lst))): continue #go to the next file name if it contains any of the substrings or matches one or more regex patterns
                    file_list.append(os.path.join(root, name)) 
            elif regstr_lst is not None and len(regstr_lst) > 0:
                if file_fmt_str in name: 
                    name_series = pd.Series([name]) #convert to pandas series for faster substring search
                    if len(regstr_lst) > 1:
                        if any(name_series.str.contains('|'.join(regstr_lst))): continue #go to the next file name if it contains any of the regex patterns
                    else:
                        if any(name_series.str.contains(regstr_lst[0])): continue #if there is only one regex pattern
                    file_list.append(os.path.join(root, name)) #if flow of excecution is here, patterns were not detected
            elif substr_lst is not None and len(substr_lst) > 0:
                if file_fmt_str in name: 
                    name_series = pd.Series([name]) #convert to pandas series for faster substring search
                    if len(substr_lst) > 1:
                        if any(name_series.str.contains('|'.join(substr_lst))): continue #go to the next file name if it contains any of the substrings
                    else:
                        if any(name_series.str.contains(substr_lst[0])): continue #if there is only one substring pattern
                    file_list.append(os.path.join(root, name)) #if flow of excecution is here, patterns were not detected
            else:
                if file_fmt_str in name: file_list.append(os.path.join(root, name))
    return file_list
                
            


def create_sample_sheet(file_lst, generic_str=None, regex_str=None):
    """
    Given (list) of paths to files and a generic part of the file name (e.g. _contigs.fasta or _R{1,2}_001.fastq.gz string)
    or sample_id regex pattern (regex string), returns pandas dataframe with sample_id column and one (fa for fasta) or two (fq1 fq2, for fastq) file path columns. 
    """

def edit_sample_sheet():
    """
    Given path to a path to the sample sheet file (str) or pandas dataframe (object), a dictionary where each sample id is matched with information to be added (dict),
    and a new column name (str), returns a pandas dataframe (object), that contains new column where new information is added to the corresponding sample id.
    """

def check_module_output():
    """
    Given (list) of paths to expected module output files, returns a dictionary where each file path is matched with the boolean (dict)
    indicating if it is present in the file system.
    """

def read_config():
    """
    Given path to a config.yaml file, return a dictionary (dict) form of the yaml file.
    """

def edit_config():
    """
    Given a dictionary (dict) that is generated from config yaml file, a (list) of parameter paths that lead to the value to be changed
    and a new value of the parameter (string), return edited dictionary were the value of specified parameter is changed.
    """

def write_config():
    """
    Given a dictionary (dict) and a path to the new config file, check if the structure of the dictionary corresponds to the config template structure
    (hardcoded or read from file), and if it fits, write the contents to the new config file.
    """

def submit_module_job():
    """
    Given snakemake module name (str) and a dictionary (dict) of module-specific parameters and job-specific parameters, 
    edit submition code string (bash template, hardcoded or read from file), 
    create temporary job script (removed after submission) and perform job submition to RTU HPC cluster.
    """
    

###Template to read config yaml file into dict
# with open("yaml_path", 'r') as yaml_handle:
#     config_dict=yaml.safe_load(yaml_handle)
#     print(config_dict)


###Templates to run shell command using subprocess
# subprocess.check_call(['qsub', '-F', f'{read_1} {read_2} {sample_id} {database}', 'job.sh'])
# subprocess.call(f'bash create_sampleSheet.sh --mode illumina --fastxDir {config["target_dir"]} --outDir {config["target_dir"]}'.split(' '))

if __name__ == "__main__":
    args = parse_arguments()
    if args.mode == "core":
        fastq_list = parse_folder(args.fastq,'.fastq.gz')
        print(fastq_list)
        
    else:
        print('Sorry, other options not supported yet.')
        
    
    """
    from fastq:
    input - path to folder with raw fastq files
    generate sample_id_list (use aquamis script for sample id generation from different fastq formats)    qsub bact_core
        check output for each sample in id list (add status and check_note column to sample_list to indicate if check is succesful and add more information if it is not): 
            status:
                fail - cannot continue (check_note: missing <file_name> from <step_name> conducted by <module name>)
                warning - can continue, except for some steps downstream (check_note: missing <file_name> from <step_name> conducted by <module name>, affecting <affected_step_name_1>... steps)
                ok - can continue (check_note: empty)
        sample_id_raw fastq (excluding undetermined) - ok/fail
        sample_id_host_filtered fastq - ok/fail
        sample_id_contigs.fasta - ok/fail
        sample_id_contig_based_taxonomy - ok/warning
        sample_list (AQUAMIS format + majority genus from reads for each sample (kraken2)) - ok/fail
    bact_core itself will generate template string that contains all wildcards to be used by the bact_core rules
        example: '{sample_id_pattern}-{reference_sequence_pattern}-scaffolds/{sample_id_pattern}-{reference_sequence_pattern}-ragtag.scaffold.fasta' 
    """