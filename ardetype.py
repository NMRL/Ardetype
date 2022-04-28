"""
This is a wrapper script of ARDETYPE(?) pipeline.
Date: 2022-04-27
Version: 0.0
"""
import sys, argparse, yaml, subprocess

###Architecture
"""
Pipeline can start from:
    a. bcl files - run all or just demultiplex
    b. fastq files - run all, except demultiplexing module in core, run all or just perform assembly
    c. contigs.fasta - run only shell, relevant tip modules and shape
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

from bcl:
    input - path to run folder
    check for sample_sheet
        if missing - exit with message indicating missing input
    qsub bcl2fastq
        if output files are missing after job is finished - exit with message indicating job error
        else proceed to fastq step

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
    generate template string that contains all wildcards to be used by the bact_core rules
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
    generate template string that contains all wildcards to be used by the bact_shell, bact_tip & bact_shape rules
        write to config_modular.yaml -> <module_name>_target/wildcard_string sections
        example: '{sample_id_pattern}-{reference_sequence_pattern}-scaffolds/{sample_id_pattern}-{reference_sequence_pattern}-ragtag.scaffold.fasta' 
    if shell: qsub shell
    if tip: qsub tip
    check results
        if results contain minimum reportable amount of data and reporting option selected
            qsub shape
            check results
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