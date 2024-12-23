import sys, os, subprocess as sp, json, time, pandas as pd
from datetime import datetime
from glob import glob

sys.path.insert(0, os.path.abspath('/home/group/pipelines/Ardetype'))
from subscripts.ardetype_utilities import Ardetype_housekeeper as hk


def update_cgmlst_files(batch_dir_path:str):
    '''Copy files for cgmlst-based cluster analysis 
    into specified location for a given processed batch folder path.
    '''
    # Statics
    history_path = f'/home/group/data_lake/analysis_history/01_Illumina/02_bact/01_current/'
    contig_path  = os.path.join(history_path,f'cgmlst_clusters/contigs_by_species')
    profile_path = os.path.join(history_path,f'cgmlst_clusters/profiles_by_species')
    tax_map_path = os.path.join(batch_dir_path, 'core_aggregated_taxonomy.json')

    # Path lists
    contig_plist  = list(glob(os.path.join(batch_dir_path, 'folded_*/*contigs.fasta')))
    profile_plist = list(glob(os.path.join(batch_dir_path, 'folded_*/*_chewbbaca/results_alleles.tsv')))

    # Taxonomy map
    with open(tax_map_path, "r+") as json_handle:             
        batch_tax_map = json.load(json_handle)
    species = set(batch_tax_map.values())


    # Creating copy of contigs in the species-specific folder
    start_time = datetime.now()
    start = time.time()

    # Copy contigs
    hk.copy_files_by_species(
        taxonomy_map = batch_tax_map,
        file_list = contig_plist,
        collection_path = contig_path,
        extension='contigs.fasta',
        batch_depth=2,
    )

    # Copy profiles
    hk.copy_files_by_species(
        taxonomy_map = batch_tax_map,
        file_list = profile_plist,
        collection_path = profile_path,
        extension='results_alleles.tsv',
        batch_depth=3,
    )
    end = time.time()
    end_time = datetime.now()

    species_joint = "\n    "+"\n    ".join(species)
    print(f'''-- Grouped contigs & cgMLST profiles by species --
Batch: {os.path.basename(batch_dir_path)}
    # of profiles: {len(profile_plist)}
    # of contigs: {len(contig_plist)}
    # of species: {len(species)}
Species: {species_joint}
Start time: {start_time}
End time: {end_time}
Finished in: {round(end-start, 4)} sec\n''')



def get_latest_file(path_list:list, suffix:str):
    cur_ts = None
    latest_file = ""
    for p in path_list:
        if p.endswith(suffix):
            fname = os.path.basename(p)
            ts = fname.replace(suffix, '')
            ts = datetime.strptime(ts, "%Y-%m-%d")
            update_cond = cur_ts is None or cur_ts < ts
            if update_cond:
                cur_ts = ts
                latest_file = p
    return latest_file

def update_seq_batches(batch_id:str):
    '''Wrapper to update key map file when new batch data is exported'''
    sbm_path = './seq_batch_map.csv'
    cmd = 'python ./update_seq_batches.py --mode update'

    # update the state of seq_batch_map based on other files

    start_time = datetime.now()
    start = time.time()
    
    sp.run(cmd, shell=True)
    
    end = time.time()
    end_time = datetime.now()
   

    # verify the new batch is present in seq_batch_map and give summary statistics
    sbm = pd.read_csv(sbm_path)
    sbm_batch = sbm[sbm['analysis_batch_id'] == batch_id]

    if sbm_batch.empty:
        print(f'''-- ERROR: Records for new batch not present in primary key table after update --
Index file: {sbm_path}
Batch: {batch_id}
Start time: {start_time}
End time: {end_time}
Finished in: {round(end-start, 4)} sec\n''')
    else:
        total_keys = sbm_batch.__len__()
        _reads_keys = sbm_batch[sbm_batch['sample_id'].str.endswith('_reads')].__len__()
        unq_keys = sbm_batch[~sbm_batch['sample_id'].str.endswith('_reads')]
        dup_keys = unq_keys['sample_id'].duplicated().sum()
        unq_keys = unq_keys[~unq_keys['sample_id'].duplicated()].__len__()
        
        print(f'''-- Index file succesfully updated --
Index file: {sbm_path}  
Batch: {batch_id}
Start time: {start_time}
End time: {end_time}
# of batch entries (total): {total_keys}
# of batch sample identifiers (unique): {unq_keys}
# of duplicate sample identifier in batch: {dup_keys}
Finished in: {round(end-start, 4)} sec\n''')



def main():
    batch_dir_path = sys.argv[1]

    # Test switches
    update_history = False
    update_cgmlst = False
    update_seq_batch = True

    if update_history:
        reports_path = os.path.join(batch_dir_path, 'reports')
        path_list = [os.path.join(reports_path, p) for p in os.listdir(reports_path)]
        latest_aquamis_report = get_latest_file(path_list, suffix="_aquamis_qc_report.csv")
        latest_ardetype_report = get_latest_file(path_list, suffix="_ardetype_report.csv")
        path_map = {
            "agnostic":{
                "k2c":os.path.join(reports_path,"kraken2contigs_report.csv"),
                "qst":os.path.join(reports_path,"quast_report.csv"),
                "aqc":latest_aquamis_report,
                "vir":os.path.join(reports_path,"virulencefinder_summary.csv")
            },
            "ardetype":{
                "ard":latest_ardetype_report
            },
            "plasmids":{
                "mbt":os.path.join(reports_path,"mobtyper_summary.csv"),
                "plf":os.path.join(reports_path,"plasmidfinder_summary.csv"),
            },
            "resistance":{
                "hamr":os.path.join(reports_path,"harmonized_resistance_profile.tsv"),
                "rf":os.path.join(reports_path,"resfinder_pheno_table_gathered.csv"),
                "pf":os.path.join(reports_path,"pointfinder_report.csv"),
                "af":os.path.join(reports_path,"amrfp_mutation_report.csv")
            },
            "software_logs":{
                "log":os.path.join(reports_path,"software_log.csv")
            },
            "specific":{
                "kbt":os.path.join(reports_path,"kleborate_report.csv"),
                "ect":os.path.join(reports_path,"ectyper_report.csv"),
                "stf":os.path.join(reports_path,"stectfinder_report.csv"),
                "agr":os.path.join(reports_path,"agrvate_report.csv"),
                "ss2":os.path.join(reports_path,"seqsero2_report.csv"),
                "sst":os.path.join(reports_path,"sistr_report.csv"),
                "lss":os.path.join(reports_path,"lissero_report.csv"),
                "mnt":os.path.join(reports_path,"meningotype_report.csv"),
                "lgt":os.path.join(reports_path,"legsta_report.csv"),
                "cbc":os.path.join(reports_path,"chewbbaca_qc_report.csv"),
                "lrf":os.path.join(reports_path,"lrefinder_report.csv"),
                "spa":os.path.join(reports_path,"spatyper_report.csv"),
                "sht":os.path.join(reports_path,"shigatyper_report.csv"),
                "sba":os.path.join(reports_path,"seroba_report.csv"),
                "emm":os.path.join(reports_path,"emmtyper_report.csv"),
                "hic":os.path.join(reports_path,"hicap_report.csv")
            }
        }

        cmd_list = []
        start_time = datetime.now()
        start = time.time()
        for dir in path_map:
            script_path = glob(f"./{dir}/*.py")[0]
            cmd = f'python {script_path}'
            for arg in path_map[dir]:
                cmd += f" --{arg} {path_map[dir][arg]}"
            sp.run(cmd, shell=True)
            cmd_list.append(cmd)

        end = time.time()
        end_time = datetime.now()
        
        cmd_joint = "\n    "+"\n    ".join(cmd_list)
        print(f'''-- Imported analysis results into history files --
Batch: {os.path.basename(batch_dir_path)}
Commands: {cmd_joint}
Start time: {start_time}
End time: {end_time}
Finished in: {round(end-start, 4)} sec\n''')

    if update_cgmlst:
        update_cgmlst_files(batch_dir_path=batch_dir_path)

    if update_seq_batch:
        batch_id = os.path.basename(batch_dir_path) if not batch_dir_path.endswith('/') else os.path.basename(os.path.dirname(batch_dir_path))
        update_seq_batches(batch_id=batch_id)

if __name__ == "__main__":
    main()