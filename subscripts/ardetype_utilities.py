import sys
import os
import time
import os
import re
import json
import pathlib
import shutil
import pandas as pd

from glob import glob
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor as ppe, as_completed, ThreadPoolExecutor
from bisect import bisect_left
sys.path.insert(0, os.path.dirname(os.path.dirname(Path(__file__).absolute())))

from subscripts.src.utilities import Housekeeper as hk


###############
#Global configs
###############

cluster_counter = 0
ardetype_path = os.path.dirname(Path(__file__).parents[0].absolute())

###############


class Ardetype_housekeeper(hk):
    '''Class extends the standard housekeeper class to implement functions required by specific pipeline'''

    @staticmethod
    def move_folder(src_path:str, tgt_path:str) -> str:
        '''Moves folder from src to tgt and returns tgt path.'''
        src_exists = os.path.isdir(src_path)

        if not src_exists:
            sys.exit('Source does not exist or is not a folder - teminating.')

        return shutil.move(src_path, tgt_path)

    @staticmethod
    def merge_paths(src_list, target_folder=None, ignore_collisions=True, exclude_files=[]):
        '''
        Combines contents from all paths supplied as src_list under target_folder.
        Excludes and deletes files from source paths specified in exclude_files.
        '''

        confirm = input(f'WARNING: merge WILL ATTEMPT TO MOVE everything FROM each --merge_from folder to --output_dir path and run the pipeline on it.\n Original folders will NOT be deleted, but logs will not be merged.\nPlease verify if you wish to proceed [yY/nN]')


        if confirm.lower() == "y":
            def move_item(src, dst):
                """Move item from src to dst."""
                if os.path.exists(dst):
                    if not dst.endswith("_mark"):
                        print(f"Collision detected at destination path '{dst}'", file=sys.stderr)
                        sys.exit(1)
                if not any(p in dst for p in ["data", "config", "reports", "scratch","logs","_mark"]):
                    shutil.move(src, dst)

            def delete_item(path):
                if os.path.isdir(path) and not os.path.islink(path):
                    shutil.rmtree(path)
                elif os.path.exists(path):
                    os.remove(path)

            if not target_folder:
                raise ValueError("target_folder argument is required")


            # Check for merge compatibility
            unmarked = []
            dtype_dict = {}
            for src in src_list:
                mark = glob(os.path.join(src,'*mark'))
                if mark:
                    dtype_dict[src] = os.path.basename(mark[0])
                else:
                    unmarked.append(src)
            
            if unmarked:
                nl = "\n"
                print(f'Was Unable to verify data type for the following batches:\n{nl.join(unmarked)}\nPlease check folder contents.', file=sys.stderr)
                sys.exit(1)
            
            dtype_check = [all(v == dtype for v in dtype_dict.values()) for dtype in ["ONT_mark","ILL_mark","HYB_mark", "FA_mark"]]
            if sum(dtype_check) != 1: # exatcly one is true
                print(f'Some batches in the --merged_from list are unmergable based on stated rules (ALLOWED MERGES: Illumina + Hybrid, Illumina + Illumina, Nanopore + Nanopore).\nPlease check data types.', file=sys.stderr)
                sys.exit(1)

            # Check for naming collisions among source directories and with the target directory
            src_name_to_paths = {}
            os.makedirs(target_folder, mode=0o775, exist_ok=True)


            for src in src_list:
                if os.path.abspath(src) == os.path.abspath(target_folder):
                    print("Source and target folders cannot be the same.", file=sys.stderr)
                    sys.exit(1)
                    
                for item in os.listdir(src):
                    dest_path = os.path.join(target_folder, item)
                    
                    if not ignore_collisions:
                        if item in src_name_to_paths:
                            print(f"Naming collision detected for '{item}' in sources '{src_name_to_paths[item]}' and '{src}'", file=sys.stderr)
                            sys.exit(1)
                        if os.path.exists(dest_path):
                            print(f"Naming collision detected for '{item}' between source '{src}' and target '{target_folder}'")
                            sys.exit(1)

                        src_name_to_paths[item] = src

            with ThreadPoolExecutor() as executor:
                futures = []
                for src in src_list:
                    for item in os.listdir(src):
                        src_item_path = os.path.join(src, item)
                        dest_item_path = os.path.join(target_folder, item)
                        if os.path.basename(item) in exclude_files:
                            futures.append(executor.submit(delete_item, src_item_path))
                        else:
                            futures.append(executor.submit(move_item, src_item_path, dest_item_path))

                for future in as_completed(futures):
                    # Handle any exceptions or errors that might have occurred
                    future.result()
        else:
            print('Batch merge aborted.')
            sys.exit(0)


    @staticmethod
    def type_fasta_scheme(contig_path: str, url:str) -> dict:
        '''
        Type contigs for single sample using pubmlst api. Returns a dictionary converted from json output of the HTTP request.
        If the HTTP request fails, raises and exception, indicating response code and response message (if any).
        The returned dictionary may not contain typing information, indicating that the typing has failed due to assembly quality/database issue/etc. 
        '''
        from rauth import OAuth1Session
        from dotenv import load_dotenv
        from subscripts.pubmlst_rest_auth import retrieve_token, get_session_token
        load_dotenv(os.path.join(ardetype_path,'config_files/.env'))
        CONSUMER_KEY = os.getenv('CONSUMER_KEY')
        CONSUMER_SECRET = os.getenv('CONSUMER_SECRET')
        SESSION_TOKEN = os.path.join(ardetype_path,'config_files/.session_token')
        route='/sequence'

        (token, secret) = None, None
        if os.path.isfile(SESSION_TOKEN):

            def older12(file_path):
                if not os.path.exists(file_path):
                    raise FileNotFoundError(f"File not found: {file_path}")
                last_modified_time = os.path.getmtime(file_path)  # Get last modification time (in seconds since epoch)
                time_diff = time.time() - last_modified_time  # Compute difference from current time
                return time_diff > 12 * 3600  # Check if older than 12 hours
            
            if older12(SESSION_TOKEN):
                print('Session token is more than 12 hours old - requesting fresh one.')
                (token, secret) = get_session_token(None, None)
            else:
                print('Obtaining session token from file.')
                (token, secret) = retrieve_token(SESSION_TOKEN)
        else:
            print('Session token file is not found - requesting fresh one.')
            (token, secret) = get_session_token(None, None)

        if route and not re.match(r"^/", route): route = "/" + route
        url = url + route

        print(f"Accessing authenticated resource ({url})...\n")
        session = OAuth1Session(CONSUMER_KEY, CONSUMER_SECRET, access_token=token, access_token_secret=secret)
        extra_params = {}
        
        if not os.path.exists(contig_path):
            print("Sequence file " + contig_path + " does not exist.")
            return

        with open(contig_path, "r") as seq_file:
            data = seq_file.read()
            extra_params["sequence"] = data
            extra_params["details"] = "true"
    
        response = session.post(
            url,
            params={},
            data=json.dumps(extra_params),
            headers={"Content-Type": "application/json"},
            header_auth=True,
        )
        
        if response.status_code == 200 or response.status_code == 201:
            if re.search("json", response.headers["content-type"], flags=0):
                return response.json()
        elif response.status_code == 400:
            print("Bad request")
            print(response.json()["message"])
        elif response.status_code == 401:
            if re.search("unauthorized", response.json()["message"]):
                print("Access denied - client is unauthorized")
                return
            else:
                if re.search("verification", response.json()["message"]):
                    print(response.json())
                    return
        else:
            print("Error:")
            print(response.text)

                    
#######################
# Quality control check
#######################
    @staticmethod
    def check_quality_thr(thresholds_db:dict, sample_id:str, batch_id:str, species:str, ardetype_report:pd.DataFrame, quast_report_path:str, thresholds_db_timestamp:str):
        '''Generate QC report using aquamis thresholds
        Applies species thresholds, where thresholds are available.
        Applies genus thresholds where genus match and thresholds are available and not covered by species-specific thr.
        Applies all Species thresholds, where specific thresholds are not available or metric not covered by specific thr.
        Returns pandas dataframe.
        '''
        #read quast report
        qst_rep = pd.read_csv(quast_report_path, sep="\t")

        #list of metrics to include
        quast_metrics = [
            'GC (%)', 
            'N50', 
            'Total length',
            '# contigs (>= 0 bp)', 
            '# contigs (>= 1000 bp)', 
            ]
        ardetype_metrics = [
            'assembly_coverageDepth',
            'contig_hit1_species_fraction', 
            'q30_rate_after'
        ]
        metrics = quast_metrics + ardetype_metrics

        #adding columns
        columns = []
        for metric in metrics:
            for suffix in ['value', 'range', 'status', 'color']:
                columns.append(f"{metric} {suffix}")

        #keeping genus from kraken2-assigned species
        genus = species.split(" ")[0]

        #get thresholds
        species_thr = thresholds_db['thresholds'].get(species,dict())
        genus_thr   = thresholds_db['thresholds'].get(genus, dict())
        all_thr     = thresholds_db['thresholds'].get('all Species')

        #threshold filters
        species_filter = ["GC (%)", "Total length", "N50", "# contigs (>= 0 bp)", "# contigs (>= 1000 bp) ", "assembly_coverageDepth"]
        genus_filter   = ["GC (%)", "Total length", "N50", "# contigs (>= 0 bp)", "# contigs (>= 1000 bp) ", "assembly_coverageDepth", "assembly_coverageDepth"]
        all_filter     = ["q30_rate_after", "assembly_coverageDepth", "contig_hit1_species_fraction"]

        #combine thresholds
        thr = {k:v for k,v in species_thr.items() if k in species_filter}
        thr.update({k:v for k,v in genus_thr.items() if k not in thr and k in genus_filter})
        thr.update({k:v for k,v in all_thr.items() if k not in thr and k in all_filter})

        #get ardetype data
        adt = ardetype_report[(ardetype_report['sample_id'] == sample_id) & (ardetype_report['analysis_batch_id'] == batch_id)][[
            'q30_raf', 
            "average_coverage", 
            "species_contig_%"
        ]]
        adt.rename(columns = {
            'q30_raf':'q30_rate_after', 
            "average_coverage": 'assembly_coverageDepth', 
            "species_contig_%":'contig_hit1_species_fraction'}, inplace=True)
        adt['contig_hit1_species_fraction'] = adt['contig_hit1_species_fraction']/100

        #get quast data
        qst = qst_rep[quast_metrics]

        #combine results
        data = pd.concat([df.reset_index(drop=True) for df in [adt, qst]], axis=1).to_dict(orient='list')

        #apply qc thresholds
        upd_data = {}
        for key in data:
            check                     = {}
            check[f"{key} value"]     = data[key][0]
            if key in thr:
                check[f"{key} range"]  = thr[key][0]['interval']

                placement = bisect_left(check[f"{key} range"], check[f"{key} value"])
                check[f"{key} color"] = [thr[key][0]["bincolor"][placement]]
                check[f"{key} value"] = [check[f"{key} value"]]

                max_score = len(check[f"{key} range"])
                score     = max_score - placement
                #if between to pass
                if thr[key][0]['binscore'][-1] == thr[key][0]['binscore'][0]:
                    if score > 0 and score < max_score:
                        check[f"{key} status"] = ['PASS']
                    else:
                        check[f"{key} status"] = ['FAIL']
                #if >= to pass
                if thr[key][0]['binscore'][-1] < thr[key][0]['binscore'][0]:
                    if score == 0:
                        check[f"{key} status"] = ['PASS']
                    else:
                        check[f"{key} status"] = ['FAIL']
                # if <= to pass
                if thr[key][0]['binscore'][-1] > thr[key][0]['binscore'][0]:
                    if score != max_score:
                        check[f"{key} status"] = ['FAIL']
                    else:
                        check[f"{key} status"] = ['PASS']
            else:
                check[f"{key} status"]   = 'N/A'
                check[f"{key} color"]    = 'N/A'
                check[f"{key} range"] = 'N/A'
            upd_data.update(check)
            if upd_data[f"{key} range"] != 'N/A': upd_data[f"{key} range"] = "-".join([str(v) for v in upd_data[f"{key} range"]])

        data.update(upd_data)
        data['sample_id']              = [sample_id]
        data['analysis_batch_id']      = [batch_id]
        data['species']                = [species]
        data['threshold_db_sync_date'] = [thresholds_db_timestamp]
        df = pd.DataFrame.from_dict(data)[['sample_id', 'analysis_batch_id', 'species']+columns+['threshold_db_sync_date']]
        return df






####################
# Agnostic extractors
####################

    @staticmethod
    def pointfinder_results(pf_results_path: str, batch: str):
        '''Returns dataframe containing sample_id and analysis_batch_id columns
        where the same sample id and batch id are matched to all findings for a given sample.'''
        df = pd.read_csv(pf_results_path, sep='\t')
        sample_id = re.sub(r'(_S[0-9]*)?_resfinder', '',
                           os.path.basename(os.path.dirname(pf_results_path)))
        df.insert(0, 'sample_id', [sample_id for _ in df.index])
        df.insert(1, 'analysis_batch_id',
                  os.path.basename(os.path.dirname(batch)))
        df = df.astype(str)
        if len(df[df['PMID'].str.contains(', ')]) > 0:
            df['PMID'] = df['PMID'].str.replace(', ', '')
        if len(df[df['Resistance'].str.contains(', ')]) > 0:
            df['Resistance'] = df['Resistance'].str.replace(', ', '')
        return df

    @staticmethod
    def respheno_results(rfp_result_path: str, batch: str) -> pd.DataFrame:
        '''Returns a dataframe - containing valid resistance genes.'''
        df = pd.read_csv(rfp_result_path, skiprows=16, sep="\t")
        sample_id = re.sub(r'(_S[0-9]*)?_resfinder_pheno.txt',
                           '', os.path.basename(rfp_result_path))

        # remove all rows that contain "No resistance" in "WGS-predicted phenotype" column
        # remove rows with missing phenotype
        df = df[df["WGS-predicted phenotype"] != "No resistance"]
        df = df[~df['Genetic background'].isna()]  # remove warnings
        # remove redundant headers
        df = df[df['# Antimicrobial'] != "# Feature_ID"]

        # Add identifiers
        df.insert(0, 'sample_id', [sample_id for _ in df.index])
        df.insert(1, 'analysis_batch_id', [os.path.basename(
            os.path.dirname(batch)) for _ in df.index])

        return df

    @staticmethod
    def plasmidfinder_results(plf_result_path: str, batch: str) -> pd.DataFrame:
        '''To combine plasmidfinder reports and map them to sample_id-batch pair.'''
        sample_id = re.sub(r'(_S[0-9]*)?_plasmidfinder', '',
                           os.path.basename(os.path.dirname(plf_result_path)))
        df = pd.read_csv(plf_result_path, sep='\t')
        df.insert(0, 'sample_id', [sample_id for _ in df.index])
        df.insert(1, 'analysis_batch_id', [os.path.basename(
            os.path.dirname(batch)) for _ in df.index])
        return df


    @staticmethod
    def virulencefinder_results(vir_result_path: str, batch: str) -> pd.DataFrame:
        '''To combine virulencefinder reports and map them to sample_id-batch pair.'''


        sample_id = re.sub(r'(_S[0-9]*)?_virulencefinder', '',
                           os.path.basename(os.path.dirname(vir_result_path)))
        try:
            df = pd.read_csv(vir_result_path, sep='\t')
        except:
            print(vir_result_path)
            df=pd.DataFrame(columns=["Database","Virulence factor","Identity","Query / Template length","Contig","Position in contig","Protein function","Accession number"])
            
        df.insert(0, 'sample_id', [sample_id for _ in df.index])
        df.insert(1, 'analysis_batch_id', [os.path.basename(
            os.path.dirname(batch)) for _ in df.index])
        return df


    @staticmethod
    def mobtyper_results(mbt_result_path: str, batch: str) -> pd.DataFrame:
        '''To combine mob_typer reports and map them to sample_id-batch pair.'''
        sample_id = re.sub(r'(_S[0-9]*)?_mob_recon',
                        '', os.path.basename(os.path.dirname(mbt_result_path)))
        df = pd.read_csv(mbt_result_path, sep='\t')
        df.rename(columns={'sample_id': 'genetic_element'}, inplace=True)
        df.insert(0, 'sample_id', [sample_id for _ in df.index])
        df.insert(1, 'analysis_batch_id', [os.path.basename(
            os.path.dirname(batch)) for _ in df.index])
        return df

    
    @staticmethod
    def mobtyper_contig_results(contig_result_path: str, batch:str) -> pd.DataFrame:
        '''To aggregate information about plasmid location in the assembly'''
        df = pd.read_csv(contig_result_path, sep='\t')
        df['sample_id'] = df['sample_id'].str.replace("_contigs","", regex=False)
        df['sample_id'] = df['sample_id'].str.replace(r"_S[0-9]*", "", regex=True)
        df = df[df['molecule_type'] == 'plasmid'] #keep only information about plasmids
        df.insert(1, 'analysis_batch_id', [os.path.basename(
            os.path.dirname(batch)) for _ in df.index])
        return df


    @staticmethod
    def kraken2contigs_results(k2c_report_path: str, batch: str) -> pd.DataFrame:
        '''To combine kraken2 contig reports and map them to sample_id-batch pair.'''
        sample_id = re.sub(
            r'(_S[0-9]*)?_kraken2_contigs_report.txt', '', os.path.basename(k2c_report_path))
        columns = [
            "cl_cov_frac",
            "cl_cov_abs",
            "tax_cov_abs",
            "rank_code",
            "taxid",
            "tax_name"
        ]
        df = pd.read_csv(k2c_report_path, sep='\t', header=None)
        df.columns = columns
        df.tax_name = df.tax_name.str.strip()
        df.insert(0, 'sample_id', [sample_id for _ in df.index])
        df.insert(1, 'analysis_batch_id', [os.path.basename(
            os.path.dirname(batch)) for _ in df.index])
        return df

    @staticmethod
    def spatyper_results(spa_result_path: str, batch: str) -> pd.DataFrame:
        '''To combine spatyper reports and map them to sample_id-batch pair.'''
        sample_id = re.sub(r'(_S[0-9]*)?_spatyper.txt', '',
                           os.path.basename(spa_result_path))
        df = pd.read_csv(spa_result_path, sep='\t')
        df.insert(0, 'sample_id', [sample_id for _ in df.index])
        df.insert(1, 'analysis_batch_id', [os.path.basename(
            os.path.dirname(batch)) for _ in df.index])
        return df


    @staticmethod
    def kraken2reads_results(k2r_report_path: str, batch: str) -> pd.DataFrame:
        '''To combine kraken2 read reports and map them to sample_id-batch pair.'''
        sample_id = re.sub(r'(_S[0-9]*)?_kraken2_reads_report.txt',
                           '', os.path.basename(k2r_report_path))
        columns = [
            "cl_cov_frac",
            "cl_cov_abs",
            "tax_cov_abs",
            "rank_code",
            "taxid",
            "tax_name"
        ]
        df = pd.read_csv(k2r_report_path, sep='\t', header=None)
        df.columns = columns
        df.tax_name = df.tax_name.str.strip()
        df.insert(0, 'sample_id', [sample_id for _ in df.index])
        df.insert(1, 'analysis_batch_id', [os.path.basename(
            os.path.dirname(batch)) for _ in df.index])
        return df


    @staticmethod
    def quast_results(qst_report_path: str, batch: str) -> pd.DataFrame:
        '''To combine quast reports and map them to sample_id-batch pair.'''
        sample_id = re.sub(r'(_S[0-9]*)?_quast', '', os.path.basename(os.path.dirname(qst_report_path)))
        df = pd.read_csv(qst_report_path, header=None, sep='\t')
        df = df.T
        df.columns = df.iloc[0]
        df = df.drop(0)
        df.drop('Assembly',axis=1, inplace=True)
        df.insert(0, 'sample_id', [sample_id for _ in df.index])
        df.insert(1, 'analysis_batch_id', [os.path.basename(
            os.path.dirname(batch)) for _ in df.index])
        return df

####################
# Specific extractors
####################

    @staticmethod
    def kleborate_results(klbt_result_path: str, batch: str) -> pd.DataFrame:
        '''To combine kleborate reports and map them to sample_id-batch pair.'''
        columns = [
            "strain", "species", "species_match", "contig_count", "N50",
            "largest_contig", "total_size", "ambiguous_bases", "QC_warnings", "ST", "virulence_score",
            "resistance_score", "num_resistance_classes", "num_resistance_genes", "Yersiniabactin",
            "YbST", "Colibactin", "CbST", "Aerobactin", "AbST", "Salmochelin", "SmST", "RmpADC",
            "RmST", "rmpA2", "wzi", "K_locus", "K_type", "K_locus_problems", "K_locus_confidence",
            "K_locus_identity", "K_locus_missing_genes", "O_locus", "O_type", "O_locus_problems",
            "O_locus_confidence", "O_locus_identity", "O_locus_missing_genes", "AGly_acquired",
            "Col_acquired", "Fcyn_acquired", "Flq_acquired", "Gly_acquired", "MLS_acquired",
            "Phe_acquired", "Rif_acquired", "Sul_acquired", "Tet_acquired", "Tgc_acquired",
            "Tmt_acquired", "Bla_acquired", "Bla_inhR_acquired", "Bla_ESBL_acquired",
            "Bla_ESBL_inhR_acquired", "Bla_Carb_acquired", "Bla_chr", "SHV_mutations",
            "Omp_mutations", "Col_mutations", "Flq_mutations", "truncated_resistance_hits",
            "spurious_resistance_hits", "Chr_ST", "gapA", "infB", "mdh", "pgi", "phoE",
            "rpoB", "tonB", "ybtS", "ybtX", "ybtQ", "ybtP", "ybtA", "irp2", "irp1", "ybtU",
            "ybtT", "ybtE", "fyuA", "clbA", "clbB", "clbC", "clbD", "clbE", "clbF", "clbG",
            "clbH", "clbI", "clbL", "clbM", "clbN", "clbO", "clbP", "clbQ", "iucA", "iucB",
            "iucC", "iucD", "iutA", "iroB", "iroC", "iroD", "iroN", "rmpA", "rmpD", "rmpC",
            "spurious_virulence_hits",
        ]
        try:
            df = pd.read_csv(klbt_result_path, sep='\t')
            df['strain'] = df['strain'].str.replace(
                r'(_S[0-9]*)?_contigs', '', regex=True)
            # reorder columns and aggregate results even if resistance scan was not performed
            df = df.reindex(columns=columns, fill_value=None)
            df.insert(1, 'analysis_batch_id', [os.path.basename(
                os.path.dirname(batch)) for _ in df.index])
            return df
        except Exception as e:
            print(e)
            return pd.DataFrame(columns=columns)

    @staticmethod
    def ectyper_results(ect_result_path: str, batch: str) -> pd.DataFrame:
        '''To combine ectyper reports and map them to sample_id-batch pair.'''
        df = pd.read_csv(ect_result_path, sep='\t')
        df.rename(columns={'Name': 'sample_id'}, inplace=True)
        df.sample_id = df.sample_id.str.replace(
            r'(_S[0-9]*)?_contigs', '', regex=True)
        df.insert(1, 'analysis_batch_id', [os.path.basename(
            os.path.dirname(batch)) for _ in df.index])
        return df


    @staticmethod
    def seroba_results(sba_result_path: str, batch: str):
        '''To combine seroba reports and map them to sample_id-batch pair.'''
        try:
            report = pd.read_csv(sba_result_path, sep="\t", header=None).astype(str)
            if len(report.columns) == 3:
                if report[2][0] != 'nan':
                    report[1] = report[1]+';'+report[2]
            else: 
                report[2] = ['nan']
            report.columns = ["sample_id","serotype", "analysis_batch_id"]
            report.analysis_batch_id = [os.path.basename(os.path.dirname(batch))]
            report["sample_id"] = report.sample_id.apply(lambda x:os.path.basename(x).split('_')[0])
            report = report[["sample_id","analysis_batch_id","serotype"]]
            return report
        except pd.errors.EmptyDataError:
            report = pd.DataFrame.from_dict({
                'sample_id':[os.path.basename(sba_result_path).replace('_seroba.tsv', '').split('_')[0]], 
                "analysis_batch_id":[os.path.basename(os.path.dirname(batch))],
                "serotype":[None],
                })
            return report


    @staticmethod
    def emmtyper_results(emm_result_path:str, batch:str):
        try:
            report = pd.read_csv(emm_result_path, sep="\t", header=None)
            report.columns = ["sample_id","num_clusters", "emm_type", "emm_like_alleles", "emm_cluster"]
            report['analysis_batch_id'] = [os.path.basename(os.path.dirname(batch))]
            report["sample_id"] = report.sample_id.apply(lambda x:x.replace('_contigs.tmp', '').split('_')[0])
            report = report[["sample_id", "analysis_batch_id", "num_clusters", "emm_type", "emm_like_alleles", "emm_cluster"]]
            return report
        except pd.errors.EmptyDataError:
            report = pd.DataFrame.from_dict({
                'sample_id':[os.path.basename(emm_result_path).replace('_strp_emmtyper.tsv', '').split('_')[0]], 
                "analysis_batch_id":[os.path.basename(os.path.dirname(batch))],
                "num_clusters":[None],
                "emm_type":[None],
                "emm_like_alleles":[None],
                "emm_cluster":[None]
                })
            return report


    @staticmethod
    def hicap_results(hi_result_path:str, batch:str):
        try:
            report = pd.read_csv(hi_result_path, sep="\t")
            report.columns = ['sample_id', 'predicted_serotype', 'attributes', 'genes_identified',
                    'locus_location', 'region_I_genes', 'region_II_genes',
                    'region_III_genes', 'IS1016_hits']
            report['analysis_batch_id'] = [os.path.basename(os.path.dirname(batch))]
            report["sample_id"] = report.sample_id.apply(lambda x:x.replace('_contigs', '').split('_')[0])
            report = report[
                ['sample_id', 'analysis_batch_id', 'predicted_serotype', 
                    'attributes', 'genes_identified','locus_location', 
                    'region_I_genes', 'region_II_genes', 'region_III_genes', 'IS1016_hits']
                ]
            return report
        except pd.errors.EmptyDataError:
            report = pd.DataFrame.from_dict({
                'sample_id':[os.path.basename(hi_result_path).replace('_hi_hicap.tsv', '').split('_')[0]], 
                'analysis_batch_id':[os.path.basename(os.path.dirname(batch))],
                'predicted_serotype':['NTHi'],
                'attributes':[None],
                'genes_identified':[None],
                'locus_location':[None],
                'region_I_genes':[None],
                'region_II_genes':[None],
                'region_III_genes':[None],
                'IS1016_hits':[None],
                })
            return report
        

    @staticmethod
    def pt_ab_results(pt_ab_result_path:str, batch:str):
        try:
            report = pd.read_csv(pt_ab_result_path, sep='\t', header=None)
            report.columns = [
                "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                "qstart", "qend", "sstart", "send", "evalue", "bitscore"
            ]
            sample_id = os.path.basename(pt_ab_result_path).replace('_abplasmidtype.txt', '')
            report["sample_id"] = [sample_id for _ in report.index]
            report["seq_batch"] = [os.path.basename(os.path.dirname(pt_ab_result_path)) for _ in report.index]
            report = report[[
                "sample_id","seq_batch",
                "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                "qstart", "qend", "sstart", "send", "evalue", "bitscore"
            ]]
            return report
        except pd.errors.EmptyDataError:
            report = pd.DataFrame.from_dict({
                'sample_id':[os.path.basename(pt_ab_result_path).replace('_abplasmidtype.txt', '')],
                'seq_batch':[os.path.basename(os.path.dirname(batch))],
                "qseqid":[None], 
                "sseqid":[None], 
                "pident":[None], 
                "length":[None], 
                "mismatch":[None], 
                "gapopen":[None],
                "qstart":[None], 
                "qend":[None], 
                "sstart":[None], 
                "send":[None], 
                "evalue":[None], 
                "bitscore":[None]
                })
            return report
    
    @staticmethod
    def capybara_results(cpb_result_path:str, batch:str):
        try:
            report = pd.read_csv(cpb_result_path, sep="\t")
            report.columns = [
                'query',
                "ESL",
                "Clade",
                "Coverage"
            ]
            report["sample_id"] = report["query"].apply(lambda x:os.path.basename(x).replace('_contigs.fasta', ''))
            report["seq_batch"] = report["query"].apply(lambda x:os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(x)))))
            report = report[[
                    "sample_id",
                    "seq_batch",
                    "ESL",
                    "Clade",
                    "Coverage"
                ]]
            return report
        except pd.errors.EmptyDataError:
            report = pd.DataFrame.from_dict({
                'sample_id':[os.path.basename(cpb_result_path).replace('_capybara.tsv', '').split('_')[0]],
                'seq_batch':[os.path.basename(os.path.dirname(batch))],
                "ESL":[None],
                "Clade":[None],
                "Coverage":[None]
                })
            return report
    

    @staticmethod
    def plassembler_results(plassembler_result_path:str, batch:str):
        columns = ["contig", "length", "mean_depth_long", "sd_depth_long", "q25_depth_long", "q75_depth_long",
                "plasmid_copy_number_long", "circularity", "PLSDB_hit", "NUCCORE_ACC", "mash_distance", "mash_pval",
                "mash_matching_hashes", "NUCCORE_UID", "NUCCORE_Description", "NUCCORE_CreateDate", "NUCCORE_Topology",
                "NUCCORE_Completeness", "NUCCORE_TaxonID", "NUCCORE_Genome", "NUCCORE_Length", "NUCCORE_DuplicatedEntry",
                "NUCCORE_Source", "NUCCORE_BiosampleID", "BIOSAMPLE_UID", "BIOSAMPLE_ACC", "BIOSAMPLE_Location",
                "BIOSAMPLE_Coordinates", "BIOSAMPLE_IsolationSource", "BIOSAMPLE_Host", "BIOSAMPLE_CollectionDate",
                "BIOSAMPLE_HostDisease", "BIOSAMPLE_SampleType", "ASSEMBLY_UID", "ASSEMBLY_ACC", "ASSEMBLY_Status",
                "ASSEMBLY_coverage", "ASSEMBLY_SeqReleaseDate", "ASSEMBLY_SubmissionDate", "ASSEMBLY_Lastest",
                "ASSEMBLY_BiosampleID", "TAXONOMY_superkingdom", "TAXONOMY_phylum", "TAXONOMY_class", "TAXONOMY_order",
                "TAXONOMY_family", "TAXONOMY_genus", "TAXONOMY_species", "TAXONOMY_strain", "TAXONOMY_UID",
                "TAXONOMY_taxon_rank", "TAXONOMY_taxon_name", "TAXONOMY_taxon_lineage", "TAXONOMY_superkingdom_id",
                "TAXONOMY_phylum_id", "TAXONOMY_class_id", "TAXONOMY_order_id", "TAXONOMY_family_id", "TAXONOMY_genus_id",
                "TAXONOMY_species_id", "TAXONOMY_strain_id", "has_biosample", "has_assembly", "has_location",
                "rMLST_hits", "rMLST_hitscount", "inclusions", "NUCCORE_GC", "Length", "BIOSAMPLE_Host_processed",
                "BIOSAMPLE_Host_processed_source", "BIOSAMPLE_Host_label", "BIOSAMPLE_HostDisease_processed",
                "loc_lat", "loc_lng", "loc_parsed", "D1", "D2", "plasmidfinder", "pmlst"
            ]
        try:
            report = pd.read_csv(plassembler_result_path, sep='\t')
            sample_id = os.path.basename(plassembler_result_path).replace('_with_sample.tsv', '')
            report["sample_id"] = [sample_id for _ in report.index]
            report["seq_batch"] = [os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(plassembler_result_path))))) for _ in report.index]
            report = report[columns]
            return report
        except pd.errors.EmptyDataError:
            report = pd.DataFrame.from_dict({
                'sample_id': [os.path.basename(plassembler_result_path).replace('_abplasmidtype.txt', '')],
                'seq_batch': [os.path.basename(os.path.dirname(batch))],
                **{col: [None] for col in columns}
            })
            return report

    @staticmethod
    def stecfinder_results(stf_result_path: str, batch: str) -> pd.DataFrame:
        '''To combine stecfinder reports and map them to sample_id-batch pair.'''
        df = pd.read_csv(stf_result_path, sep='\t')
        df.rename(columns={'Sample': 'sample_id'}, inplace=True)
        df.sample_id = df.sample_id.str.replace(
            r'(_S[0-9]*)?_bact_reads_classified', '', regex=True)
        df.insert(1, 'analysis_batch_id', [os.path.basename(
            os.path.dirname(batch)) for _ in df.index])
        return df

    @staticmethod
    def agrvate_results(agv_result_path: str, batch: str) -> pd.DataFrame:
        '''To combine agrvate reports and map them to sample_id-batch pair.'''
        df = pd.read_csv(agv_result_path, sep='\t')
        df.rename(columns={'#filename': 'sample_id'}, inplace=True)
        df.sample_id = df.sample_id.str.replace(
            r'(_S[0-9]*)?_contigs', '', regex=True)
        df.insert(1, 'analysis_batch_id', [os.path.basename(
            os.path.dirname(batch)) for _ in df.index])
        return df

    @staticmethod
    def legsta_results(lgs_result_path: str, batch: str) -> pd.DataFrame:
        '''To combine legsta reports and map them to sample_id-batch pair.'''
        df = pd.read_csv(lgs_result_path)
        df.rename(columns={'FILE': 'sample_id'}, inplace=True)
        df['sample_id'] = df['sample_id'].apply(
            lambda x: os.path.basename(x)).astype(str)
        df.sample_id = df.sample_id.str.replace(
            r'(_S[0-9]*)?_contigs.fasta', '', regex=True)
        df.insert(1, 'analysis_batch_id', [os.path.basename(
            os.path.dirname(batch)) for _ in df.index])
        return df

    @staticmethod
    def meningotype_results(mnt_result_path: str, batch: str) -> pd.DataFrame:
        '''To combine meningotype reports and map them to sample_id-batch pair.'''
        df = pd.read_csv(mnt_result_path, sep='\t')
        df.rename(columns={'SAMPLE_ID': 'sample_id'}, inplace=True)
        df['sample_id'] = df['sample_id'].apply(
            lambda x: os.path.basename(x)).astype(str)
        df.sample_id = df.sample_id.str.replace(
            r'(_S[0-9]*)?_contigs.fasta', '', regex=True)
        df.insert(1, 'analysis_batch_id', [os.path.basename(
            os.path.dirname(batch)) for _ in df.index])
        return df

    @staticmethod
    def lissero_results(lss_result_path: str, batch: str) -> pd.DataFrame:
        '''To combine lissero reports and map them to sample_id-batch pair.'''
        df = pd.read_csv(lss_result_path, sep='\t')
        df.rename(columns={'ID': 'sample_id'}, inplace=True)
        df['sample_id'] = df['sample_id'].apply(
            lambda x: os.path.basename(x)).astype(str)
        df.sample_id = df.sample_id.str.replace(
            r'(_S[0-9]*)?_contigs.fasta', '', regex=True)
        df.insert(1, 'analysis_batch_id', [os.path.basename(
            os.path.dirname(batch)) for _ in df.index])
        return df

    @staticmethod
    def sistr_results(ssr_result_path: str, batch: str):
        '''To combine sistr reports and map them to sample_id-batch pair.'''
        df = pd.read_csv(ssr_result_path)
        df.rename(columns={'genome': 'sample_id'}, inplace=True)
        df.drop('fasta_filepath', inplace=True, axis=1)
        df['sample_id'] = df['sample_id'].apply(
            lambda x: os.path.basename(x)).astype(str)
        df.sample_id = df.sample_id.str.replace(
            r'(_S[0-9]*)?_contigs', '', regex=True)
        df.insert(7, 'analysis_batch_id', [os.path.basename(
            os.path.dirname(batch)) for _ in df.index])
        return df

    @staticmethod
    def seqsero2_results(sqs2_result_path: str, batch: str):
        '''To combine seqsero2 reports and map them to sample_id-batch pair.'''
        df = pd.read_csv(sqs2_result_path, sep='\t')
        df.rename(columns={'Sample name': 'sample_id'}, inplace=True)
        df.drop('Output directory', inplace=True, axis=1)
        df.drop('Input files', inplace=True, axis=1)
        df['sample_id'] = df['sample_id'].astype(str)
        df['sample_id'] = df['sample_id'].apply(
            lambda x: os.path.basename(x)).astype(str)
        df.sample_id = df.sample_id.str.replace(
            r'_S[0-9]*', '', regex=True)
        df.insert(1, 'analysis_batch_id', [os.path.basename(
            os.path.dirname(batch)) for _ in df.index])
        return df

    @staticmethod
    def amrfpm_results(amrfpm_result_path: str, batch: str):
        '''To combine amrfinder+ mutation reports and map them to sample_id-batch pair.'''
        df = pd.read_csv(amrfpm_result_path, sep='\t')
        sample_id = re.sub(r'(_S[0-9]*)?_amrfinderplus_point.tab', '', os.path.basename(amrfpm_result_path))
        df.insert(0, 'sample_id', [sample_id for _ in df.index])
        df.insert(1, 'analysis_batch_id', [os.path.basename(os.path.dirname(batch)) for _ in df.index])
        return df
    
    @staticmethod
    def rmlst_results(rmlst_result_path: str, batch: str):
        '''To combine rmlst species reports and map them to sample_id'''
        # Parse json file
        try:
            with open(rmlst_result_path, 'r') as f:
                data = json.load(f)

            # Get species (if present) and sample identifier
            species = data.get('fields', {}).get('species', '')
            taxons = data.get('taxon_prediction',[])
            sample_id = os.path.basename(rmlst_result_path).replace('_rmlst.json','')

            # Get top species if mix inferred
            if species:
                splt = species.split(",")
                if len(splt) > 1:
                    species = splt[0]
            else:
                if taxons:
                    species = max(taxons, key=lambda x: x['support'])['taxon']

            # Convert to pandas df
            df = pd.DataFrame.from_dict({"sample_id":[sample_id], "species": species})
            return df
        except Exception as e:
            print(f'WARNING: rMLST failed with exception: {e}')
            sample_id = os.path.basename(rmlst_result_path).replace('_rmlst.json','')
            df = pd.DataFrame.from_dict({"sample_id":[sample_id], "species":[""]})
            return df
        
    @staticmethod
    def plsdb_results(plsdb_result_path: str, batch: str):
        try:
            df = pd.read_csv(plsdb_result_path)
            columns = ['sample_id', 'seq_batch', 'query'] + list(df.columns)[:-3]
            df = df[columns]
        except:
            df = pd.DataFrame(columns='NUCCORE_ACC,NUCCORE_Topology,Length,NUCCORE_GC,NUCCORE_Source,ASSEMBLY__ASSEMBLY_ACC,BIOSAMPLE__BIOSAMPLE_ACC,BIOSAMPLE__LOCATION_name,BIOSAMPLE__ECOSYSTEM_tags,BIOSAMPLE__ECOSYSTEM_taxid,BIOSAMPLE__DISEASE_ontid,TAXONOMY__TAXONOMY_UID,TAXONOMY__TAXONOMY_taxon_name,TAXONOMY__TAXONOMY_taxon_rank,BIOSAMPLE__ECOSYSTEM_query,NUCCORE_CreateDate,TAXONOMY__TAXONOMY_class,TAXONOMY__TAXONOMY_family,TAXONOMY__TAXONOMY_genus,TAXONOMY__TAXONOMY_order,TAXONOMY__TAXONOMY_phylum,TAXONOMY__TAXONOMY_species,TAXONOMY__TAXONOMY_strain,TAXONOMY__TAXONOMY_superkingdom,BIOSAMPLE__LOCATION_lat,BIOSAMPLE__LOCATION_lng,UMAP_D1,UMAP_D2,MOB__PMLST_scheme,MOB__rep_type,MOB__relaxase_type,MOB__mpf_type,qseqid,evalue,bitscore,qstart,qend,sstart,send,qcov,qcovhsp,pident,identity,shared_hashes,median_multiplicity,pvalue,query,seq_batch,sample_id'.split(','))
        return df

    @staticmethod
    def lrefinder_results(lrefinder_pos_path:str, batch: str):
        '''To combine lrefinder mutation reports and map them to sample_id-batch pair.'''
        df = pd.read_csv(lrefinder_pos_path, sep='\t')
        sample_id = re.sub(r'(_S[0-9]*)?.pos', '', os.path.basename(lrefinder_pos_path))
        df.insert(0, 'sample_id', [sample_id for _ in df.index])
        df.insert(1, 'analysis_batch_id', [os.path.basename(os.path.dirname(batch)) for _ in df.index])
        return df

    @staticmethod
    def shigatyper_results(shigatyper_result_path:str, batch: str):
        '''To combine shigatyper reports and map them to sample_id-batch pair.'''
        df = pd.read_csv(shigatyper_result_path, sep='\t')
        df.rename(columns={'sample':'sample_id'}, inplace=True)
        df.insert(1, 'analysis_batch_id', [os.path.basename(os.path.dirname(batch)) for _ in df.index])
        return df
    
    @staticmethod
    def cgmlst_quality_results(lrefinder_pos_path:str, batch: str):
        '''To combine chewbbaca allele calling reports and map them to sample_id-batch pair.'''
        df = pd.read_csv(lrefinder_pos_path, sep='\t')
        df.rename(columns={'FILE': 'sample_id'}, inplace=True)
        df.sample_id = df.sample_id.str.replace(r'(_S[0-9]*)?_contigs', '', regex=True)
        df.insert(1, 'analysis_batch_id', [os.path.basename(os.path.dirname(batch)) for _ in df.index])
        return df

####################
#Special aggregation
####################

    @staticmethod
    def hamr_plasmid_map(hamr_report_path:str, mbt_contig_report_path:str) -> pd.DataFrame:
        '''Matching plasmid reports with contig reports for Resfinder & mobtyper to search for potentially plasmid-localized AMRs'''
        hr = pd.read_csv(hamr_report_path, sep='\t')
        cr = pd.read_csv(mbt_contig_report_path)
        hr = hr[hr.reference_database_name == 'resfinder']
        hr.input_sequence_id = hr.input_sequence_id.str.replace(" ", "_")
        hr.input_sequence_id = hr.input_sequence_id.astype(str)
        hr.input_file_name = hr.input_file_name.astype(str)
        cr.sample_id = cr.sample_id.astype(str)
        cr.sample_id = cr.sample_id.str.replace(r"_S[0-9]*", '', regex=True)
        hamr_pls_merge = hr.merge(cr, 
                                  how='right', 
                                  left_on=['input_file_name', 'input_sequence_id'], 
                                  right_on=['sample_id', 'contig_id'])
        hamr_pls_merge.dropna(subset=['input_file_name'], inplace=True)
        hamr_pls_merge = hamr_pls_merge[[
            'sample_id', 
            'analysis_batch_id_x', 
            'gene_name', 
            'sequence_identity',
            'coverage_percentage',
            'analysis_software_version', 
            'reference_database_version',
            'rep_type(s)',
            'primary_cluster_id',
            'contig_id',
            "input_gene_start",
            "input_gene_stop",
            "reference_accession",
            "rep_type_accession(s)"
            ]]
        hamr_pls_merge.rename(columns={'analysis_batch_id_x': 'analysis_batch_id'}, inplace=True)
        return hamr_pls_merge


###########################
# Multiprocessing aggregator
###########################

    @staticmethod
    def aggregator(outfolder_path: str, proc_num: int, wildcard: str = None, extractor=pointfinder_results, pathlist: list = None):
        '''Runs result extractors in-parallel - used in snakefiles for aggregation of specific results.'''
        if wildcard is not None:
            pf_out_list = pathlib.Path(outfolder_path).rglob(wildcard)
        elif pathlist is None and wildcard is None:
            raise Exception(
                f'Either "pathlist" or "wildcard" should be not None.')
        else:
            pf_out_list = pathlist
        summary = pd.DataFrame()
        with ppe(max_workers=proc_num) as executor:
            results = [executor.submit(
                extractor, path, outfolder_path) for path in pf_out_list]
            for output in as_completed(results):
                summary = pd.concat([summary, output.result()])
        summary = summary.reset_index(drop=True)
        return summary
