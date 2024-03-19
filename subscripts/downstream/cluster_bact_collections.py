import os
import pandas as pd
import shutil
import glob
import subprocess
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed

# Constants
PROFILE_COLLECTIONS = '/mnt/beegfs2/home/groups/nmrl/bact_analysis/tax_profile_collections/'
CONTIG_COLLECTIONS = '/mnt/beegfs2/home/groups/nmrl/bact_analysis/tax_contig_collections/'
CLUSTERING_SCRIPT_PATH = '/mnt/beegfs2/home/groups/nmrl/utils/phylogenetics_tools/process_chewbacca_profiles_edits.py'
CSP2_PATH = '/mnt/beegfs2/home/groups/nmrl/bact_analysis/phylogenetics/CSP2'


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
            if os.path.isdir(nf_work_dir):
                shutil.rmtree(nf_work_dir)
            os.makedirs(nf_work_dir)
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
            
            cmd = f"source activate nf && cd {CSP2_PATH} && nextflow -C {config_path} -log {nf_log_path} run -w {nf_work_dir} CSP2.nf --runmode snp --fasta {target_dir} --outroot {species_contig_path} --out {label}_CSP2 --ref_fasta {ref_fasta}"
            print(cmd)
            try:
                print(f"Running CSP2 for {target_dir}")
                result = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            except subprocess.CalledProcessError as e:
                print(f"CSP2 failed with error:\n{e.stderr}")
            except Exception as e:
                print(f"CSP2 - an error occurred: {e}")


def process_species_folder(species_folder):
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
    clusters_df['cluster_20_thr'] = clusters_df.iloc[:, 1].fillna('outgroup')
    clusters = clusters_df.groupby('cluster_20_thr')

    prepare_contig_collections(species_name, clusters.groups.keys())
    move_files_to_clusters(species_name, clusters)
    create_cluster_distance_matrices(species_name, clusters, distance_df)
    run_refchooser(species_name, clusters)
    #run_csp2(species_name, clusters)


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
    """Moves files into their respective cluster directories."""
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
    
    # Execute cgmost clustering concurrently
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
        futures = [executor.submit(process_species_folder, profile_path) for profile_path in profile_paths]
        for future in as_completed(futures):
            try:
                future.result()
            except Exception as e:
                print(f"An error occurred: {e}")

# Ensure runtime block is correct
if __name__ == "__main__":
    main()
