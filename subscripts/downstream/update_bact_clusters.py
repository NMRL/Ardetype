import os
import sys
import pandas as pd
import glob
import shutil
from concurrent.futures import ProcessPoolExecutor, as_completed
sys.path.insert(0,'/mnt/beegfs2/home/groups/nmrl/bact_analysis/Ardetype/')
from subscripts.downstream import update_utilities as uu

# Constants
PROFILE_COLLECTIONS = '/mnt/beegfs2/home/groups/nmrl/bact_analysis/tax_profile_collections/'
CONTIG_COLLECTIONS = '/mnt/beegfs2/home/groups/nmrl/bact_analysis/tax_contig_collections/'
REPORT_TIME = uu.get_current_timestamp()
OUTPUT_PATH = '/mnt/beegfs2/home/groups/nmrl/bact_analysis/analysis_history/bact_clusters/'
BACKUP_PATH = os.path.join(OUTPUT_PATH, 'backup')

# Functions
def get_collections(ctype:str):
    '''Get list of full paths for each subfolder in respective collection.'''
    if ctype=='profiles':
        collection_paths = [os.path.join(PROFILE_COLLECTIONS, dirname) for dirname in os.listdir(PROFILE_COLLECTIONS)]
    elif ctype=='contigs':
        collection_paths = [os.path.join(CONTIG_COLLECTIONS, dirname) for dirname in os.listdir(CONTIG_COLLECTIONS)]
    return collection_paths


def process_species_folder(species_folder):
    """Processes each species folder to handle clustering result extraction."""
    species_name = os.path.basename(species_folder)
    cluster_files = glob.glob(os.path.join(species_folder, 'clusters*.csv'))

    if not cluster_files:
        print(f"No cluster files found for {species_name}, skipping.")
        return

    # Assuming the first file is what we want for both clusters
    clusters_df = pd.read_csv(cluster_files[0])
    clusters_df['cluster_20_thr'] = clusters_df.iloc[:, 1].fillna('outgroup')
    clusters_df['species'] = species_name.replace("_", " ")
    return clusters_df


# Runtime
def main():
    # Extract clustering results
    profile_paths = get_collections(ctype='profiles')
    aggregated_clusters = pd.DataFrame()
    with ProcessPoolExecutor(max_workers=28) as executor:
        futures = [executor.submit(process_species_folder, profile_path) for profile_path in profile_paths]
        for future in as_completed(futures):
            try:
                result = future.result()
                if result is not None:
                    aggregated_clusters = pd.concat([aggregated_clusters, result])
            except Exception as e:
                print(f"An error occurred: {e}")

    #Find current file
    current_cluster_file = glob.glob(os.path.join(OUTPUT_PATH, 'bact_clusters*.csv'))
    if current_cluster_file:
        current_cluster_file = glob.glob(os.path.join(OUTPUT_PATH, 'bact_clusters*.csv'))[0]
        cur_file_name = os.path.basename(current_cluster_file)
        #Move to backup
        shutil.move(current_cluster_file, os.path.join(BACKUP_PATH, cur_file_name))
    aggregated_clusters.to_csv(os.path.join(OUTPUT_PATH, f'bact_clusters_{REPORT_TIME}.csv'), header=True, index=False)

if __name__ == "__main__":
    main()