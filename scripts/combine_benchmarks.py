import pandas as pd, sys, argparse, os, time

parser = argparse.ArgumentParser(description='A script to combine benchmark reports by snakemake.')
req_arg_grp = parser.add_argument_group('Required arguments.')
req_arg_grp.add_argument('-b', '--benchmark', metavar='\b', help='Path to the folder containing benchmark reports.')
req_arg_grp.add_argument('-s', '--sample_id', metavar='\b', help='Sample id based on which to filter the benchmark reports.')

args = parser.parse_args()
if args.benchmark is None or args.sample_id is None:
    parser.print_help(sys.stderr)
    sys.exit(1)

folder_path,sample_id = args.benchmark,args.sample_id
output_df_indexes = []
output_df_exist = []
for file in os.listdir(folder_path):
    if sample_id in file and 'bbmap_qc' in file:
        output_df_indexes.append('bbmap_df')
        bbmap_df = pd.read_csv(f'{folder_path}/{file}', delimiter='\t')
        output_df_exist.append( bbmap_df)
    elif sample_id in file and 'fastp' in file:
        output_df_indexes.append('fastp_df')
        fastp_df = pd.read_csv(f'{folder_path}/{file}', delimiter='\t')
        output_df_exist.append( fastp_df)
    elif sample_id in file and 'ragtag' in file:
        output_df_indexes.append(f'ragtag_df_{file.split("-",1)[1]}')
        ragtag_df = pd.read_csv(f'{folder_path}/{file}', delimiter='\t')
        output_df_exist.append(ragtag_df)
    elif sample_id in file and 'shovill' in file:
        output_df_indexes.append('shovill_df')
        shovill_df = pd.read_csv(f'{folder_path}/{file}', delimiter='\t')
        output_df_exist.append( shovill_df)
    elif sample_id in file and 'mlst' in file:
        output_df_indexes.append('mlst_df')
        mlst_df = pd.read_csv(f'{folder_path}/{file}', delimiter='\t')
        output_df_exist.append( mlst_df)
total_df = pd.concat(output_df_exist)


total_df.columns = ['Runtime_in_seconds', 'Runtime_in_h:m:s', 'RAM_size_MB', 'Virtual_memory_size_MB', 'Unique_Allocated_Memory_MB', 'Shared_Allocated_Memory_MB(Linux)', 'Number_of_MB_read', 'Number_of_MB_written', 'CPU_usage(%)_over time/Runtime_in_seconds', 'CPU_time(seconds)_summed_for_user&system']
sum_dict = {
    'Runtime_in_seconds': total_df['Runtime_in_seconds'].sum(), 
    'Runtime_in_h:m:s': time.strftime('%H:%M:%S', time.gmtime(total_df['Runtime_in_seconds'].sum())),
    'RAM_size_MB':total_df['RAM_size_MB'].sum(), 
    'Virtual_memory_size_MB':total_df['Virtual_memory_size_MB'].sum(), 
    'Unique_Allocated_Memory_MB':total_df['Unique_Allocated_Memory_MB'].sum(), 
    'Shared_Allocated_Memory_MB(Linux)':total_df['Shared_Allocated_Memory_MB(Linux)'].sum(), 
    'Number_of_MB_read':total_df['Number_of_MB_read'].sum(), 
    'Number_of_MB_written':total_df['Number_of_MB_written'].sum(), 
    'CPU_usage(%)_over time/Runtime_in_seconds':total_df['CPU_usage(%)_over time/Runtime_in_seconds'].sum(), 
    'CPU_time(seconds)_summed_for_user&system':total_df['CPU_time(seconds)_summed_for_user&system'].sum()
}
total_df = total_df.append([sum_dict], True)
output_df_indexes = list(set(output_df_indexes))
output_df_indexes.append('sum')
total_df.index = output_df_indexes
total_df.to_csv(f'{folder_path}/{sample_id}_combined_benchmark.csv', header=True, index=True)