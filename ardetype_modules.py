from ardetype_utilities import parse_folder, create_sample_sheet, read_config, edit_config, write_config, submit_module_job, run_module_cluster, check_job_completion, check_module_output, edit_sample_sheet, validate_config
import os, warnings, re, pandas as pd
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)


def run_core(args):
    """This is a function that runs bact_core module of the pipeline, after receiving a namespace (object) with command line arguments"""
    file_list = parse_folder(args.input,'.fastq.gz')
    fastq_patterns = "(_R[1,2]_001.fastq.gz|_[1,2].fastq.gz|_R[1,2]_001_unclassified_out)"
    try:
        sample_sheet = create_sample_sheet(file_list, fastq_patterns, mode=0)
        os.system(f"mkdir -p {args.output_dir}")
        sample_sheet.to_csv(f"{args.output_dir}sample_sheet.csv", header=True, index=False)
    except AssertionError as msg:
        print(f"bact_core : sample sheet generation error : {msg}")
    
    target_list, template_list  = [], [
        "_contigs.fasta",
        "_bact_reads_classified_1.fastq.gz", 
        "_bact_reads_classified_2.fastq.gz",
        "_bact_reads_unclassified_1.fastq.gz",
        "_bact_reads_unclassified_2.fastq.gz",
        "_kraken2_contigs_report.txt",
        "_kraken2_host_filtering_report.txt"
    ]
    [target_list.append(f'{args.output_dir}{id}{tmpl}') for id in sample_sheet['sample_id'] for tmpl in template_list]
    
    config_file = read_config(args.config)
    edit_result = edit_config(config_file, "core_target_files", target_list)
    assert edit_result == 0, f"bact_core: editing of the config file failed while trying to set 'core_target_files' value"
    edit_result = edit_config(config_file, "output_directory", args.output_dir)
    assert edit_result == 0, f"bact_core: editing of the config file failed while trying to set 'output_directory' value"
    validation_code = validate_config(config_file)
    assert validation_code == 0, f"bact_core: validation of the config file failed with code {validation_code}"
    write_config(config_file, f'{args.output_dir}config.yaml')
    
    config_file_path = f'{os.path.abspath(args.output_dir)}/config.yaml'
    cluster_config_path = 'cluster.yaml'

    if args.module_jobs:
        job_id = submit_module_job('core',config_file_path, args.output_dir)
        check_job_completion(job_id,"bact_core",sleeping_time=5,output_dir=args.output_dir)
    else:
        run_module_cluster("core",config_file_path, cluster_config_path, 12)
    
    check_dict = check_module_output(file_list=target_list)
    id_check_dict = {id:"" for id in sample_sheet['sample_id']}
    
    for file in check_dict:
        split = os.path.basename(file)
        for tmpl in template_list: split = re.sub(tmpl,"",split)
        id_check_dict[split] += f"|{file}:{check_dict[file]}"
    
    sample_sheet = edit_sample_sheet(sample_sheet, id_check_dict, "check_note_core")
    
    sample_sheet.to_csv(f"{args.output_dir}sample_sheet.csv", header=True, index=False)
    if args.mode == "all":
        return sample_sheet

def run_shell(args):
    """This is a function that runs bact_shell module of the pipeline, after receiving a namespace (object) with command line arguments"""
    fasta_patterns = '_contigs.fasta'
    fastq_patterns = "(_R[1,2]_001.fastq.gz|_[1,2].fastq.gz)"
    if args.mode == "all":
        file_list = parse_folder(args.output_dir, fasta_patterns)
        sample_sheet = pd.read_csv(f"{args.output_dir}/sample_sheet.csv")
        fasta_dict = {id:"" for id in sample_sheet["sample_id"]}
        id_extractor = lambda x: os.path.basename(x).replace(fasta_patterns, "")
        for file in file_list: fasta_dict[id_extractor(file)] = file
        edit_sample_sheet(sample_sheet, fasta_dict, "fa")
        sample_sheet.to_csv(f"{args.output_dir}sample_sheet.csv", header=True, index=False)
    else:
        fasta_list = parse_folder(args.input, fasta_patterns)
        fastq_list = parse_folder(args.input, fastq_patterns)

        try:
            os.system(f"mkdir -p {args.output_dir}")
            if len(fastq_list) >= 2:
                file_series = pd.Series(fastq_list, dtype="str") #to fascilitate filtering
                id_extractor = lambda x: re.sub(fastq_patterns,"",os.path.basename(x)) #extract id from string by using regex
                id_series = file_series.apply(id_extractor).drop_duplicates(keep = "first").sort_values().reset_index(drop=True)
                read_1_dict, read_2_dict = {}, {} #to use python mapping to ensure correspondance between id and path
                for id in id_series:
                    read_files = file_series[file_series.str.contains(id)].reset_index(drop=True).sort_values().reset_index(drop=True) #extract read paths
                    read_1_dict[id] = read_files[0]
                    read_2_dict[id] = read_files[1]
                sample_sheet = create_sample_sheet(fasta_list, fasta_patterns, mode=1)
                sample_sheet['fq1'] = sample_sheet['sample_id'].map(read_1_dict)
                sample_sheet['fq2'] = sample_sheet['sample_id'].map(read_2_dict)
            sample_sheet.to_csv(f"{args.output_dir}sample_sheet.csv", header=True, index=False)
        except AssertionError as msg:
            print(f"bact_shell: sample sheet generation error : {msg}")

    target_list, template_list = [], [
        ".rgi.txt",
        ".rgi.json",
        "_mlst_output.csv",
        "_resfinder/pheno_table.txt",
        "_resfinder/ResFinder_Hit_in_genome_seq.fsa",
        "_resfinder/ResFinder_Resistance_gene_seq.fsa",
        "_resfinder/ResFinder_results_tab.txt",
        "_resfinder/ResFinder_results.txt",
        "_amrpp/ResistomeResults/AMR_analytic_matrix.csv"
    ]
    [target_list.append(f'{args.output_dir}{id}{tmpl}') for id in sample_sheet['sample_id'] for tmpl in template_list]
    config_file = read_config(args.config)
    edit_result = edit_config(config_file, "shell_target_files", target_list)
    assert edit_result == 0, f"bact_shell: editing of the config file failed while trying to set 'shell_target_files' value"
    edit_result = edit_config(config_file, "output_directory", args.output_dir)
    assert edit_result == 0, f"bact_shell: editing of the config file failed while trying to set 'output_directory' value"
    validation_code = validate_config(config_file)
    assert validation_code == 0, f"bact_shell: validation of the config file failed with code {validation_code}"
    write_config(config_file, f'{args.output_dir}config.yaml')

    config_file_path = f'{os.path.abspath(args.output_dir)}/config.yaml'
    cluster_config_path = 'cluster.yaml'

    if args.module_jobs:
        job_id = submit_module_job('shell',config_file_path, args.output_dir)
        check_job_completion(job_id,"bact_shell",sleeping_time=5,output_dir=args.output_dir)
    else:
        run_module_cluster("shell",config_file_path, cluster_config_path, 12)
    
    check_dict = check_module_output(file_list=target_list)
    id_check_dict = {id:"" for id in sample_sheet['sample_id']}

    for file in check_dict:
        split = file.rsplit("/",3)
        if "resfinder" in split[-2]:
            split = f"{split[-2]}/{split[-1]}"
        elif "ResistomeResults" in split[-2]:
            split = f"{split[-3]}/{split[-2]}/{split[-1]}"
        else:
            split = split[-1]
        for tmpl in template_list: split = re.sub(tmpl,"",split)
        id_check_dict[split] += f"|{file}:{check_dict[file]}"
        
    sample_sheet = edit_sample_sheet(sample_sheet, id_check_dict, "check_note_shell")
    sample_sheet.to_csv(f"{args.output_dir}sample_sheet.csv", header=True, index=False)
    if args.mode == "all":
        return sample_sheet

def run_tip(args):
    """This is a function that runs bact_tip module of the pipeline, after receiving a namespace (object) with command line arguments"""


def run_shape(args):
    """This is a function that runs bact_shape module of the pipeline, after receiving a namespace (object) with command line arguments"""