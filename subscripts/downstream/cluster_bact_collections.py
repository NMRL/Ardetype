import os
import sys
import pandas as pd
import shutil
import glob
import subprocess
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
sys.path.insert(0,'/mnt/beegfs2/home/groups/nmrl/bact_analysis/Ardetype/')
from subscripts.downstream import update_utilities as uu

# Constants
NMRL='/mnt/beegfs2/home/groups/nmrl/'
WD = os.path.join(NMRL,'bact_analysis/analysis_history/01_Illumina/01_current/')
PROFILE_COLLECTIONS = os.path.join(WD, 'cgmlst_clusters/profiles_by_species')
CONTIG_COLLECTIONS = os.path.join(WD, 'cgmlst_clusters/contigs_by_species')
CLUSTERING_SCRIPT_PATH = os.path.join(NMRL,'utils/phylogenetics_tools/process_chewbacca_profiles_edits.py')
CSP2_PATH = os.path.join(WD,'phylogenetics/CSP2')
REPORT_TIME = uu.get_current_timestamp()
OUTPUT_PATH = os.path.join(WD,'cgmlst_clusters/')
BACKUP_PATH = os.path.join(OUTPUT_PATH, 'backup')


# Functions
def get_collections(ctype:str):
    '''Get list of full paths for each subfolder in respective collection.'''
    if ctype=='profiles':
        collection_paths = [os.path.join(PROFILE_COLLECTIONS, dirname) for dirname in os.listdir(PROFILE_COLLECTIONS)]
    elif ctype=='contigs':
        collection_paths = [os.path.join(CONTIG_COLLECTIONS, dirname) for dirname in os.listdir(CONTIG_COLLECTIONS)]
    return collection_paths


def run_clustering_script(*args, **kwargs):
    # Extract the directory where the clustering script is located
    script_directory = os.path.dirname(CLUSTERING_SCRIPT_PATH)
    
    # Construct the command as a single string with 'cd' and then 'python script'
    cmd = f"cd {script_directory} && python {os.path.basename(CLUSTERING_SCRIPT_PATH)} " \
          f"-p {kwargs['collection_name']} -wc {kwargs['wildcard']} -s " \
          f"-t {kwargs['distance_threshold']} -st True -sc {kwargs['min_sample_count']} -sd True"
    
    # Execute the command using subprocess.run with shell=True
    try:
        result = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    except subprocess.CalledProcessError as e:
        print(f"Clustering failed with error:\n{e.stderr}")
    except Exception as e:
        print(f"Clustering - an error occurred: {e}")


def run_refchooser(species_name, clusters):
    '''Runs refchooser to generate top-1 reference for each species and each cluster (including outgroup)'''
    species_contig_path = os.path.join(CONTIG_COLLECTIONS, species_name)
    sketch_dir = os.path.join(CONTIG_COLLECTIONS, 'sketches')
    for label, group in clusters:
        target_dir = os.path.join(species_contig_path, label)
        num_files = len(os.listdir(target_dir))
        refch_rank_path = os.path.join(species_contig_path, f'{label}_refchooser.tsv')
        cmd = f"source activate nf && refchooser metrics --sort Score --top {num_files} {target_dir} {sketch_dir} > {refch_rank_path}"

        try:
            print(f"Running Refchooser for {target_dir}")
            result = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        except subprocess.CalledProcessError as e:
            print(f"Refchooser failed with error:\n{e.stderr}")
        except Exception as e:
            print(f"Refchooser - an error occurred: {e}")


def run_csp2(species_name, clusters):
    '''Runs CSP2 based on refchooser-selected reference to attempt constructing cluster-specific SNP-based distance matrix (excluding outgroup)'''
    species_contig_path = os.path.join(CONTIG_COLLECTIONS, species_name)
    for label, group in clusters:
        if label != 'outgroup':
            target_dir = os.path.join(species_contig_path, label)
            refch_rank_path = os.path.join(species_contig_path, f'{label}_refchooser.tsv')
            refch_rank = pd.read_csv(refch_rank_path, sep='\t')
            nf_log_path = os.path.join(species_contig_path, f'{label}_nf.log')
            config_path = os.path.join(species_contig_path, f'{label}_nf.config')
            nf_work_dir = os.path.join(species_contig_path, f'{label}_nf_wd')
            nfx_home = os.path.join(species_contig_path, f'{label}_nf_home')

            for f in [nf_work_dir, nfx_home]:
                if os.path.isdir(f):
                    shutil.rmtree(f)
                os.makedirs(f, mode=0o775)
            ref_fasta = list(refch_rank['Path'])[0]
            output_path = os.path.join(species_contig_path, f'{label}_CSP2')
            # Remove old results if present
            if os.path.isdir(output_path):
                shutil.rmtree(output_path)
            # Add config file for cluster if not already present
            if not os.path.isfile(config_path):
                shutil.copy(os.path.join(CSP2_PATH, 'nextflow.config'), config_path)
            # Remove old logs if present
            if os.path.isfile(nf_log_path):
                os.remove(nf_log_path)
            
            cmd = f"source activate nf && cd {CSP2_PATH} && export NXF_HOME={nfx_home} && nextflow -C {config_path} -log {nf_log_path} run -w {nf_work_dir} CSP2.nf --runmode snp --fasta {target_dir} --outroot {species_contig_path} --out {label}_CSP2 --ref_fasta {ref_fasta}"
            try:
                print(f"Running CSP2 for {target_dir}")
                result = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            except subprocess.CalledProcessError as e:
                print(f"CSP2 failed with error:\n{e.stderr}")
            except Exception as e:
                print(f"CSP2 - an error occurred: {e}")


def process_species_folder(species_folder, th=20):
    """Processes each species folder to handle clustering and file organization."""
    species_name = os.path.basename(species_folder)
    cluster_files = glob.glob(os.path.join(species_folder, 'clusters*.csv'))
    distance_files = glob.glob(os.path.join(species_folder, 'distances*.tsv'))

    if not cluster_files or not distance_files:
        print(f"No cluster or distance file found for {species_name}, skipping.")
        return

    # Assuming the first file is what we want for both clusters and distances
    clusters_df = pd.read_csv(cluster_files[0])
    distance_df = pd.read_csv(distance_files[0], index_col=0, sep='\t')
    # Assuming the cluster label is in the second column
    clusters_df[f'cluster_{th}_thr'] = clusters_df.iloc[:, 1].fillna('outgroup')
    clusters = clusters_df.groupby(f'cluster_{th}_thr')

    prepare_contig_collections(species_name, clusters.groups.keys())
    move_files_to_clusters(species_name, clusters)
    create_cluster_distance_matrices(species_name, clusters, distance_df)
    # run_refchooser(species_name, clusters)
    # run_csp2(species_name, clusters)


def get_species_clusters(species_folder, th=20):
    """Processes each species folder to handle clustering result extraction."""
    species_name = os.path.basename(species_folder)
    cluster_files = glob.glob(os.path.join(species_folder, 'clusters*.csv'))

    if not cluster_files:
        print(f"No cluster files found for {species_name}, skipping.")
        return

    # Assuming the first file is what we want for both clusters
    clusters_df = pd.read_csv(cluster_files[0])
    clusters_df[f'cluster_{th}_thr'] = clusters_df.iloc[:, 1].fillna('outgroup')
    clusters_df['species'] = species_name.replace("_", " ")
    return clusters_df


def prepare_contig_collections(species_name, cluster_labels):
    """Prepares the directories in CONTIG_COLLECTIONS based on cluster labels."""
    species_contig_path = os.path.join(CONTIG_COLLECTIONS, species_name)

    # Clean up and reorganize existing directories
    existing_folders = {folder for folder in os.listdir(species_contig_path) if os.path.isdir(os.path.join(species_contig_path, folder))}
    for folder in existing_folders:
        folder_path = os.path.join(species_contig_path, folder)
        
        for file in os.listdir(folder_path):
            if '_contigs.fasta' in file:
                shutil.move(os.path.join(folder_path, file), species_contig_path)
        shutil.rmtree(folder_path)
    
    #Remove old distance matrices

    # Create new directories based on cluster labels
    for label in cluster_labels:
        os.makedirs(os.path.join(species_contig_path, label), exist_ok=True)

def move_files_to_clusters(species_name, clusters):
    """Moves contig files into their respective cluster directories."""
    species_contig_path = os.path.join(CONTIG_COLLECTIONS, species_name)

    for label, group in clusters:
        target_dir = os.path.join(species_contig_path, label)
        for sample_id in group['sample_id']:
            # Assuming contigs follow a naming convention that includes the sample_id
            for contig_file in glob.glob(os.path.join(species_contig_path, f"{sample_id}*.fasta")):
                shutil.move(contig_file, target_dir)


def create_cluster_distance_matrices(species_name, clusters, distance_df):
    '''Create cgmlst distance matrix only for samples present in cluster'''
    species_contig_path = os.path.join(CONTIG_COLLECTIONS, species_name)
    for label, group in clusters:
        cluster_samples = group['sample_id']
        cluster_distance_matrix = distance_df.loc[cluster_samples, cluster_samples]
        cluster_distance_matrix.to_csv(os.path.join(species_contig_path, f"{label}_distance_matrix.tsv"), sep='\t')


def main():
    # Constant arguments
    kwargs = {
        'wildcard': '*_result*_alleles.tsv',
        'distance_threshold': 15,
        'min_sample_count': 3
    }
    
    # Get collection paths
    profile_paths = get_collections(ctype='profiles')
    
    # Create a task for each collection path
    tasks = [
        {**kwargs, 'collection_name': path}
        for path in profile_paths
    ]
    
    #Execute cgmlst clustering concurrently
    with ThreadPoolExecutor(max_workers=len(tasks)) as executor:
        future_to_task = {executor.submit(run_clustering_script, **task): task for task in tasks}
        
        for future in as_completed(future_to_task):
            task = future_to_task[future]
            try:
                future.result()  # Trigger any exceptions caught during task execution
                print(f"Task completed successfully for: {task}")
            except Exception as e:
                print(f"Task failed for {task}: {e}")


    # Prepare contigs collections and run refchooser concurrently
    profile_paths = get_collections(ctype='profiles')
    with ThreadPoolExecutor(max_workers=28) as executor:
        futures = [executor.submit(process_species_folder, profile_path, kwargs['distance_threshold']) for profile_path in profile_paths]
        for future in as_completed(futures):
            try:
                future.result()
            except Exception as e:
                print(f"An error occurred: {e}")

    #Aggregate clusters for all species
    profile_paths = get_collections(ctype='profiles')
    aggregated_clusters = pd.DataFrame()
    with ProcessPoolExecutor(max_workers=28) as executor:
        futures = [executor.submit(get_species_clusters, profile_path, kwargs['distance_threshold']) for profile_path in profile_paths]
        for future in as_completed(futures):
            try:
                result = future.result()
                if result is not None:
                    aggregated_clusters = pd.concat([aggregated_clusters, result])
            except Exception as e:
                print(f"An error occurred: {e}")

    #Find current aggregation file
    current_date_string = datetime.now().strftime('%Y-%m-%d')
    aggregated_clusters['date_cur'] = current_date_string
    current_cluster_file = glob.glob(os.path.join(OUTPUT_PATH, 'bact_clusters*.csv'))
    current_clustsize_file = glob.glob(os.path.join(OUTPUT_PATH, 'bclust_sizes*.csv'))
    #Moving old clustsize file to backup
    if current_clustsize_file:
        current_clustsize_file = current_clustsize_file[0]
        cur_sfile_name = os.path.basename(current_clustsize_file)
        shutil.move(current_clustsize_file, os.path.join(BACKUP_PATH, cur_sfile_name))
    #Combining old and new cluster registry files
    if current_cluster_file:
        current_cluster_file = current_cluster_file[0]
        cur_file_name = os.path.basename(current_cluster_file)
        old_clusters = pd.read_csv(current_cluster_file)
        #Move to backup
        shutil.move(current_cluster_file, os.path.join(BACKUP_PATH, cur_file_name))

        #Check if new clusters differ from old clusters
        check = old_clusters.merge(aggregated_clusters, on='sample_id', how='inner')
        # Get the current date as a string.
        
        check = check[check[f'cluster_{kwargs["distance_threshold"]}_thr_x'] != check[f'cluster_{kwargs["distance_threshold"]}_thr_y']]
        # check.to_csv('test.csv',header=True, index=False)
        if len(check) > 0:

            old_clusters.set_index('sample_id', inplace=True)
            check.set_index('sample_id', inplace=True)

            #in case cluster ids were updated during last iteration
            for index,row in check.iterrows():
                # Assuming 'cluster_{kwargs["distance_threshold"]}_thr' is the current cluster ID column
                new_cluster_id = row[f'cluster_{kwargs["distance_threshold"]}_thr_y']
                new_date = row['date_cur_y']  # Assuming this is correctly named and contains the current date string
                
                # Check if cluster_hist is NaN and update accordingly
                if pd.isna(old_clusters.at[index, 'cluster_hist']):
                    old_clusters.at[index, 'cluster_hist'] = new_cluster_id
                else:
                    old_id = old_clusters.at[index, 'cluster_hist'].split('_')[-1]
                    new_id = str(new_cluster_id).split('_')[-1]
                    print(index, old_id, new_id, new_cluster_id)
                    if old_id != new_id:
                        old_clusters.at[index, 'cluster_hist'] += "_" + new_id
                        old_clusters.at[index, f'cluster_{kwargs["distance_threshold"]}_thr'] = new_cluster_id

                # Check if date_hist is NaN and update accordingly
                if pd.isna(old_clusters.at[index, 'date_hist']):
                    old_clusters.at[index, 'date_hist'] = new_date
                else:
                    old_date = old_clusters.at[index, 'date_hist'].split('_')[-1]
                    if old_date != new_date:
                        old_clusters.at[index, 'date_hist'] += "_" + new_date
                        old_clusters.at[index, 'date_cur'] = new_date
        #Generate unique cluster identifier based on first assigned cluster# and species
        cluster_ids = []
        first_dates = []
        for row in old_clusters.iterrows():
            clust_hist = row[1]['cluster_hist'].split('_')
            first_dates.append(row[1]['date_hist'].split('_')[0])
            prefix = ''
            if clust_hist[0] == 'outgroup':
                prefix += 'outgroup'
            else:
                prefix += 'cluster'
                prefix += '_'+clust_hist[1]
            cluster_id = prefix + "_" + row[1]['species'].replace(" ", "-")
            cluster_ids.append(cluster_id)
        old_clusters['cluster_id'] = cluster_ids
        old_clusters['first_date'] = first_dates

        #Calculating cluster sizes
        clust_sizes = old_clusters.groupby([f'cluster_{kwargs["distance_threshold"]}_thr', 'species']).count()
        #Transforming to get sql-accepted shape
        clust_sizes = clust_sizes.reset_index()
        # clust_sizes.to_csv('test1.csv')
        #Keeping relevant columns
        clust_sizes = clust_sizes[[f'cluster_{kwargs["distance_threshold"]}_thr', 'species', 'date_cur']]
        #Applying informative column names
        clust_sizes = clust_sizes.rename(columns={'date_cur':'sample_count'})

        #Save updated history in case it was not a "cold start"
        # old_clusters.reset_index(inplace=True)
        old_clusters.to_csv(os.path.join(OUTPUT_PATH, f'bact_clusters_{REPORT_TIME}.csv'), header=True, index=False)
        #Save up-to-date cluster size table
        clust_sizes.to_csv(os.path.join(OUTPUT_PATH, f'bclust_sizes_{REPORT_TIME}.csv'), header=True, index=False)
        return


    #Calculating cluster sizes
    clust_sizes = aggregated_clusters.groupby([f'cluster_{kwargs["distance_threshold"]}_thr', 'species']).count()
    #Transforming to get sql-accepted shape
    clust_sizes = clust_sizes.reset_index()
    #Keeping relevant columns
    clust_sizes = clust_sizes[[f'cluster_{kwargs["distance_threshold"]}_thr', 'species', 'sample_id']]
    #Applying informative column names
    clust_sizes = clust_sizes.rename(columns={'sample_id':'sample_count'})

    #Save up-to-date clusters table
    aggregated_clusters.to_csv(os.path.join(OUTPUT_PATH, f'bact_clusters_{REPORT_TIME}.csv'), header=True, index=False)
    #Save up-to-date cluster size table
    clust_sizes.to_csv(os.path.join(OUTPUT_PATH, f'bclust_sizes_{REPORT_TIME}.csv'), header=True, index=False)

# Runtime
if __name__ == "__main__":
    main()
