import matplotlib.pyplot as plt
import os
import logging
import datetime
import argparse
import pandas as pd
import subprocess
import sys
import numpy as np
import re

from glob import glob
from collections import defaultdict
from scipy.cluster import hierarchy as sch
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform


'''
Development notes

08.09.2023 - JB
1. Moved iterative clustering logic from <gather_chewbacca_data> to <cluster_analysis> to keep the "one function - one task".
2. In <gather_chewbacca_data> replaced os.walk with glob + wildcard to avoid iterating over irrelevant files.
3. Added commented sections for readability.
4. Wrapped logging and argparsing into functions in helper section.
5. Added singularity image for cgmlst-dists (0.4.0) and an option to run the script using it.
6. Added __main__ logic to avoid running the script if something is imported from it.
7. Replace ClusterCounter class with global variable cluster_counter in perform_clustering.
8. Added docstrings.
9. Added squareform fix when there is only one sample per bin.

11.09.2023 - JB 
1. Finalized fix of 1-sample-per-bin case.
2. Added proper visualization for clustering results as dendrogram + thresholds.

18.09.2022 - JB
1. Added thresholds as command-line argument
2. Added filter_distance_outliers function to remove columns and rows of totally unrelated samples
3. Added option to filter_distance_outliers to keep the isolate closest to other isolates as outgroup.
'''


###############
#Global configs
###############

current_datetime = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
cluster_counter = 0

#######################
#Helpers & housekeeping
#######################

def parse_args():
    '''Parse arguments'''
    # define arguments and paths for running functions
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--path", required=True, help="Path to folder", type=str)
    parser.add_argument("-l", "--linkage", required=False, help="Linkage method: single, complete, average, ward", type=str, default='single')
    parser.add_argument("-ro", "--rem_outlier", required=False, help="Option to filter outliers from distance matrix; provide numeric distance threshold; all isolates that are further than threshold from nearest isolate will be excluded from clustering", type=int, default=None)
    parser.add_argument("-wc", "--wildcard", required=False, help="Option to specify wildcard to be used to identify allelic profiles in the <path> folder.", type=str, default="*_profile.tsv")
    parser.add_argument('-t', "--thresholds", help = "cgmlst distances for primary and secondary clustering; at least 1 value expected; separate by <space> if many", nargs='+', type=int, default=[5, 15])
    parser.add_argument('-st', "--skip_timestamp", help = "Set to True to avoid adding timestamp to output file names", type=bool, default=False)
    parser.add_argument('-sd', "--skip_dendrogram", help = "Set to True to avoid generating dendrogram image", type=bool, default=False)
    parser.add_argument("-sc", "--min_sample_count", help = 'Set to int value to avoid analyzing species with too few samples', type=int, default=1)
    parser.add_argument('-s', "--singularity", help = "Use singularity container for cgmlst-dists istead of local installation", action='store_true')
    parser.add_argument("-ls", "--log_scale", help="Log-scale the dendrogram image", action='store_true')  # New argument
    
    args = parser.parse_args()
    # show help text if run without required arguments
    if len(sys.argv)==1: 
        parser.print_help(sys.stderr)
        sys.exit(1)
    return args


def initialize_logging(path, current_datetime:datetime=current_datetime, use_timestamp=True) -> None:
    '''
    args - arparse namespace containing path argument
    current_datetime - datetime object
    '''
    if use_timestamp:
        log_file = os.path.join(path, f"aggregator_{current_datetime}.log")
    else:
        log_file = os.path.join(path, f"aggregator.log")
    logging.basicConfig(
        filename=log_file,
        level=logging.INFO,
        format="%(asctime)s:%(levelname)s:%(message)s"
        )


def plot_dendrogram(clusters, labels, thresholds, img_path, x_label, y_label):
    '''
    clusters - scipy linkage output, 
    labels - labels corresponding to distance matrix, 
    thresholds - thresholds to plot horizontal lines, 
    img_path - path to save the image, 
    x_label - x axis label,
    y_label - y axis label,
    '''
    
    sch.dendrogram(clusters, labels=labels)
    
    fig = plt.gcf()  # Get the current figure
    ax = plt.gca()
    xmin, xmax = ax.get_xlim()
    
    fig.set_size_inches(w=1.5*len(labels),h=1.5*len(labels))  # Set the width and height in centimeters
    for h in thresholds: plt.hlines(h, xmin=xmin, xmax=xmax, colors='r')
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    
    plt.savefig(img_path, dpi=500, bbox_inches='tight')
    plt.close()
    


def filter_distance_outliers(dist_matr_df, outlier_threshold=500, keep_outgroup=False) -> pd.DataFrame:
    '''
    Note - will not work correctly for negative thresholds!
    dist_matr_df - distance matrix as pandas dataframe
    outlier_threshold - distance value indicating that isolates are most probably unrelated
        if all distance values in a row/column are above this value, 
        the row and column for the corresponding isolate, 
        will be filtered from the distance matrix dataframe
        keep_outgroup - coupled together with outlier_threshold, keeps one outlier isolate closest to non-outliers as outgroup

    logic note:
        df_np > outlier_threshold - compare all values to threshold
        np.sum - count number of above threshold values for every column (axis = 1)
            max value for this count is number of rows in the column - 1, 
            since diagonal values in the matrix are always 0 and below threshold
        df_np.shape[0] - number of rows in each column
        if difference is 1, only the diagonal value is below the threshold and isolate is to be removed
        np.where returns indexes of rows and columns to be removed
    '''
    df_np = dist_matr_df.values
    #get indexes where all elements except diagonal are above threshold
    #since distance matrix is symmetric, row and columns indexes will be the same
    outlier_indexes = np.where(df_np.shape[0] - np.sum(df_np > outlier_threshold, axis=1) == 1) #see logic note in docstring

    #makes sense only if outliers exist
    if keep_outgroup and outlier_indexes[0].size != 0:
        #to select closest
        outliers = df_np[:, outlier_indexes[0]] 
        
        #absolute distance is appropriate for non-negative counts
        manhattan_distances = np.sum(np.abs(outliers), axis=0)
        
        #getting the index of the closest outlier
        min_distance_column_index = np.argmin(manhattan_distances)

        #removing closest outlier from filtering list
        outlier_indexes = (outlier_indexes[0][outlier_indexes[0] != outlier_indexes[0][min_distance_column_index]],)

    #remove outlier columns
    
    df_np_filtered = np.delete(df_np, outlier_indexes, axis=1)
    #remove outlier rows
    df_np_filtered = np.delete(df_np_filtered, outlier_indexes, axis=0)
    #update row indexes and column names
    rows = columns = [dist_matr_df.columns[i] for i in range(dist_matr_df.shape[1]) if i not in outlier_indexes[0]]
    
    #restore dataframe
    df_filtered = pd.DataFrame(df_np_filtered, columns = columns, index=rows)
    df_filtered.index.name = dist_matr_df.index.name

    #sort rows by index in incr. order to improve consistency
    df_filtered.sort_index(inplace=True)

    #reorder columns in incr. order (same)
    df_filtered = df_filtered[df_filtered.index]

    return df_filtered


###########
#Functions
###########

def gather_chewbacca_data(path:str, wildcard:str = "*_profile.tsv", min_sample_count:int = 1) -> defaultdict: #"folded*/*_chewbbaca/results_alleles.tsv"
    '''
    path - path to folder containing allelic profiles produced by chewbbaca
    wildcard - pattern used to locate chewbbaca allelic profiles in the <path>
    '''
    # create a dictionary of dataframes, one for each bin
    df_dict = defaultdict(lambda: pd.DataFrame())
    path_list = list(glob(os.path.join(path,wildcard)))
    if len(path_list) < min_sample_count:
        sys.exit(f'Less than {min_sample_count} files in {path} - Terminating.')
    
    # loop through all files matching the wildcard pattern
    for chew_path in path_list:
        try:
            df_temp = pd.read_csv(chew_path, sep="\t")
            if df_temp.shape[0] > 2:  # Ignore files with more than 2 rows
                logging.warning(f"File {chew_path} has more than 2 rows and will be ignored.")
                print(f"Warning: File {chew_path} has more than 2 rows and will be ignored.")
                continue
            
            # trim sample_id's in "FILE" column
            if 'FILE' in df_temp.columns:
                if '_results_alleles.tsv' in chew_path:
                    df_temp["FILE"] = os.path.basename(chew_path).replace('_results_alleles.tsv', '')
                else:
                    df_temp["FILE"] = os.path.basename(chew_path).replace(wildcard.replace('/','').replace('*',''), '') #df_temp["FILE"].str.split("_").str[0]
            elif "Genome" in df_temp.columns:
                if '_results_alleles.tsv' in chew_path:
                    df_temp["FILE"] = os.path.basename(chew_path).replace('_results_alleles.tsv', '')
                else:
                    df_temp["Genome"] = os.path.basename(chew_path).replace(wildcard.replace('/','').replace('*',''), '') #df_temp["Genome"].str.split("_").str[0]
            
            # get number of columns as bin number minux sample_id column
            bin_number = df_temp.shape[1] - 1

            df_dict[bin_number] = pd.concat([df_dict[bin_number], df_temp])
            logging.info(f"Processed file: {chew_path}")
        except Exception as e:
            logging.error(f"Error processing file {chew_path}: {e}")
    return df_dict


def perform_clustering(distance_file:str, thresholds:list, linkage_method:str, log_scale=False, outlier_filter:int=None, make_dendrogram:bool=True) -> pd.DataFrame:
    '''
    distance_file - path to the distance matrix in tsv format
    threshold - allelic distancee threshold to evaluate (list)
    hierarchical clustering method - see https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html for details
        requires globally-defined counter to be used accross all iterations <cluster_counter>
    log_scale - use log y-axis in dendrogram plot
    outlier_filter - perform outlier filtering based on provided value
    '''
    # Read the distance matrix
    distance_df = pd.read_csv(distance_file, sep="\t", index_col=0)
    if outlier_filter is not None:
        #update distance matrix
        distance_df = filter_distance_outliers(distance_df, outlier_threshold=outlier_filter)
        #save updated matrix
        outliers_rem_path = str(os.path.basename(distance_file)).replace(".tsv",f"_filtered_{outlier_filter}.tsv")
        outliers_rem_path = os.path.join(os.path.dirname(distance_file), outliers_rem_path)
        distance_df.to_csv(outliers_rem_path, header=True, index=True, sep='\t')

    # Convert the DataFrame to a distance matrix
    distance_matrix = squareform(distance_df.values)
    if distance_matrix.size == 0:
        if len(distance_df.columns) > 0:
            map_dict = {"sample_id":[distance_df.columns[0]]}
        else:
            map_dict = {"sample_id":['None']}
        for thr in thresholds:
            if len(distance_df.columns) == 0:
                map_dict[f"cluster_{thr}_thr"] = ["None"]
            else:    
                map_dict[f"cluster_{thr}_thr"] = ["single_per_bin"]
        return pd.DataFrame(map_dict)

    clusters_df = pd.DataFrame(distance_df.index)
    clusters_df.columns = ['sample_id']
    
    for thr in thresholds:
        #visualising heirclust results
        bin_stamp = str(os.path.basename(distance_file)).replace(".tsv","").replace("distance_","")
        img_path  = os.path.join(os.path.dirname(distance_file),f"dendrogram_{bin_stamp}.svg")
        
        if log_scale:
            clusters       = linkage(np.log1p(distance_matrix), method=linkage_method)
            cluster_labels = fcluster(clusters, t=np.log1p(thr), criterion='distance')
        else:
            clusters       = linkage(distance_matrix, method=linkage_method)
            cluster_labels = fcluster(clusters, t=thr, criterion='distance')

        if not os.path.isfile(img_path):
            y_label = f'Allelic distance defined by {linkage_method} linkage'
            if log_scale:
                y_label = f'log10(Allelic distance) defined by {linkage_method} linkage'
                thresholds = [np.log1p(threshold) for threshold in thresholds]
            
            if make_dendrogram:
                plot_dendrogram(
                    clusters=clusters, 
                    thresholds=thresholds,
                    labels=distance_df.index,
                    img_path=img_path,
                    x_label=f'{bin_stamp}',
                    y_label=y_label,
                    )

        # Create DataFrame from cluster labels
        temp_df = pd.DataFrame(cluster_labels, index=distance_df.index, columns=["cluster_id"])
        temp_df = temp_df[temp_df.duplicated(keep=False)]  # Only keep duplicated clusters (size > 1)

        # Prepend "cluster_" and ensure unique cluster_ids across files
        cluster_ids           = sorted(temp_df["cluster_id"].unique())
        unique_cluster_id_map = {}
        for old in cluster_ids:
            global cluster_counter
            cluster_counter += 1
            unique_cluster_id_map[old]= "cluster_" + str(cluster_counter)
        temp_df["cluster_id"] = temp_df["cluster_id"].map(unique_cluster_id_map)
        clusters_df           = pd.merge(clusters_df, temp_df, how="left", left_on="sample_id", right_index=True)
        clusters_df.rename(columns={"cluster_id": f"cluster_{thr}_thr"}, inplace=True)
    return clusters_df


def cluster_analysis(
        df_dict:dict, 
        path:str,  
        output_path:str=None,
        thresholds:list=[5, 15], 
        linkage_method:str='average', 
        singularity:bool=False, 
        singularity_path:str='/home/group/pipelines/Ardetype/subscripts/downstream/cgmlst-dists.sif',
        current_datetime:datetime=current_datetime,
        outlier_filter:int=None,
        log_scale=False,
        use_timestamp=True,
        skip_dendrogram=False) -> None:
    '''
    df_dict - profiles mapped to bin number
    threshold - allelic distancee threshold to evaluate (list)
    heirarchical clustering method - see https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html for details
    singularity - bool - if it is required to use singularity container instead of local installation for cgmlst-dists
    singularity_path - path to the singularity container
    log_scale - use log y scale in dendrogram plots
    outlier_filter - perform outlier filtering based on provided value
    '''
    
    if output_path is None:
        output_path = path
    
    # save dataframes as csv files and run cgmlst-dists for each
    for bin_number, df in df_dict.items():
        if use_timestamp:
            tsv_path = os.path.join(output_path, f"{bin_number}_{current_datetime}.tsv")
            dist_file = os.path.join(output_path, f"distances_{bin_number}_{current_datetime}.tsv")
        else:
            tsv_path = os.path.join(output_path, f"{bin_number}.tsv")
            dist_file = os.path.join(output_path, f"distances_{bin_number}.tsv")
        
        df = df.drop_duplicates(subset=['FILE'])
        df.to_csv(tsv_path, index=False, sep="\t")

       
        if singularity:
            cmd       = f"singularity --silent exec --bind {tsv_path},{dist_file}:{tsv_path},{dist_file} {singularity_path} cgmlst-dists {tsv_path} > {dist_file}"
        else:
            cmd       = f"cgmlst-dists {tsv_path} > {dist_file}"

        try:
            subprocess.run(cmd, shell=True, check=True)
            logging.info(f"Successfully ran cgmlst-dists for {tsv_path}. Distances written to {dist_file}")

            # Perform clustering on distances file
            cluster_df   = perform_clustering(
                dist_file, 
                thresholds=thresholds, 
                linkage_method=linkage_method, 
                log_scale=log_scale, 
                outlier_filter=outlier_filter, 
                make_dendrogram=not skip_dendrogram
                )

            if use_timestamp:
                cluster_path = os.path.join(output_path, f"clusters_{bin_number}_{current_datetime}.csv")
            else:
                cluster_path = os.path.join(output_path, f"clusters_{bin_number}.csv")
            cluster_df.to_csv(cluster_path, index=False)
        except subprocess.CalledProcessError as e:
            logging.error(f"Error occurred while running cgmlst-dists for {tsv_path}")
            logging.error(e)
            raise e


##############
#Runtime logic
##############

def main():
    args = parse_args()
    initialize_logging(args.path, current_datetime, use_timestamp = not args.skip_timestamp)
    chew_dict = gather_chewbacca_data(args.path, wildcard=args.wildcard, min_sample_count=args.min_sample_count)
    if args.singularity:
        cluster_analysis(
            path = args.path, 
            thresholds = args.thresholds, 
            df_dict=chew_dict, 
            singularity=True, 
            linkage_method=args.linkage, 
            log_scale=args.log_scale, 
            outlier_filter=args.rem_outlier,
            use_timestamp=not args.skip_timestamp,
            skip_dendrogram=args.skip_dendrogram)
    else:
        cluster_analysis(
            path = args.path, 
            thresholds = args.thresholds, 
            df_dict=chew_dict, 
            linkage_method=args.linkage, 
            log_scale=args.log_scale, 
            outlier_filter=args.rem_outlier,
            use_timestamp=not args.skip_timestamp,
            skip_dendrogram=args.skip_dendrogram)


if __name__ == "__main__":
    main()
