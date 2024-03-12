import os
import pandas as pd
import shutil
import glob
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed

# Constants
PROFILE_COLLECTIONS = '/mnt/beegfs2/home/groups/nmrl/bact_analysis/tax_profile_collections_copy/'
CONTIG_COLLECTIONS = '/mnt/beegfs2/home/groups/nmrl/bact_analysis/tax_contig_collections_copy/'
CLUSTERING_SCRIPT_PATH = '/mnt/beegfs2/home/groups/nmrl/utils/phylogenetics_tools/process_chewbacca_profiles_edits.py'


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
        print(f"Command failed with error:\n{e.stderr}")
    except Exception as e:
        print(f"An error occurred: {e}")


def process_species_folder(species_folder):
    """Processes each species folder to handle clustering and file organization."""
    species_name = os.path.basename(species_folder)
    cluster_files = glob.glob(os.path.join(species_folder, 'clusters*.csv'))
    if not cluster_files:
        print(f"No cluster file found for {species_name}, skipping.")
        return

    df = pd.read_csv(cluster_files[0])
    # Assuming the cluster label is in the second column
    df['cluster_20_thr'] = df.iloc[:, 1].fillna('outgroup')
    clusters = df.groupby('cluster_20_thr')

    prepare_contig_collections(species_name, clusters.groups.keys())
    move_files_to_clusters(species_name, clusters)

def prepare_contig_collections(species_name, cluster_labels):
    """Prepares the directories in CONTIG_COLLECTIONS based on cluster labels."""
    species_contig_path = os.path.join(CONTIG_COLLECTIONS, species_name)

    # Clean up and reorganize existing directories
    existing_folders = {folder for folder in os.listdir(species_contig_path) if os.path.isdir(os.path.join(species_contig_path, folder))}
    for folder in existing_folders:
        folder_path = os.path.join(species_contig_path, folder)
        for file in os.listdir(folder_path):
            shutil.move(os.path.join(folder_path, file), species_contig_path)
        os.rmdir(folder_path)

    # Create new directories based on cluster labels
    for label in cluster_labels:
        os.makedirs(os.path.join(species_contig_path, label), exist_ok=True)

def move_files_to_clusters(species_name, clusters):
    """Moves files into their respective cluster directories."""
    species_contig_path = os.path.join(CONTIG_COLLECTIONS, species_name)

    for label, group in clusters:
        target_dir = os.path.join(species_contig_path, label)
        for sample_id in group['sample_id']:
            # Assuming contigs follow a naming convention that includes the sample_id
            for contig_file in glob.glob(os.path.join(species_contig_path, f"{sample_id}*.fasta")):
                shutil.move(contig_file, target_dir)


def main():
    # Constant arguments
    kwargs = {
        'wildcard': '*_result*_alleles.tsv',
        'distance_threshold': 20,
        'min_sample_count': 3
    }
    
    # Get collection paths
    profile_paths = get_collections(ctype='profiles')
    
    # Create a task for each collection path
    tasks = [
        {**kwargs, 'collection_name': path}
        for path in profile_paths
    ]
    
    # Execute tasks concurrently
    with ProcessPoolExecutor(max_workers=len(tasks)) as executor:
        future_to_task = {executor.submit(run_clustering_script, **task): task for task in tasks}
        
        for future in as_completed(future_to_task):
            task = future_to_task[future]
            try:
                future.result()  # Trigger any exceptions caught during task execution
                print(f"Task completed successfully for: {task}")
            except Exception as e:
                print(f"Task failed for {task}: {e}")


    profile_paths = get_collections(ctype='profiles')
    with ProcessPoolExecutor(max_workers=os.cpu_count()) as executor:
        futures = [executor.submit(process_species_folder, profile_path) for profile_path in profile_paths]
        for future in as_completed(futures):
            try:
                future.result()
            except Exception as e:
                print(f"An error occurred: {e}")

# Ensure runtime block is correct
if __name__ == "__main__":
    main()