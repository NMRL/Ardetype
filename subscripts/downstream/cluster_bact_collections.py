import os
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed

# Constants
PROFILE_COLLECTIONS = '/mnt/beegfs2/home/groups/nmrl/bact_analysis/tax_profile_collections/'
CONTIG_COLLECTIONS = '/mnt/beegfs2/home/groups/nmrl/bact_analysis/tax_contig_collections/'
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



def main():
    # Constant arguments
    kwargs = {
        'wildcard': '*_result_alleles.tsv',
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

# Ensure runtime block is correct
if __name__ == "__main__":
    main()