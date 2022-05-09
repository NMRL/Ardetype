from ardetype_utilities import parse_folder, create_sample_sheet, read_config, edit_config, write_config, submit_module_job, run_module_cluster, check_job_completion, check_module_output, edit_sample_sheet
import os, warnings, re
warnings.simplefilter(action='ignore', category=FutureWarning)

def run_core(args):
    """This is a function that runs bact_core module of the pipeline, after receiving a namespace (object)"""
    file_list = parse_folder(args.fastq,'.fastq.gz')
    fastq_formats = "(_S[0-9]*_R[1,2]_001.fastq.gz|_[1,2].fastq.gz|_S[0-9]*_R[1,2]_001_unclassified_out)"
    try:
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
        for tmpl in template_list: split = re.sub(tmpl,"",split)
        id_check_dict[split] += f"|{split}:{check_dict[file]}"
        
    sample_sheet = edit_sample_sheet(sample_sheet, id_check_dict, "check_note")
    sample_sheet.to_csv(f"{args.output_dir}sample_sheet.csv", header=True, index=False)