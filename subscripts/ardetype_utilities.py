import sys
import os
import requests
import base64
import re
import pathlib
import pandas as pd
from pathlib import Path
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor as ppe, as_completed, ThreadPoolExecutor
from bisect import bisect_left
sys.path.insert(0, os.path.dirname(os.path.dirname(Path(__file__).absolute())))
sys.path.insert(0, '/mnt/beegfs2/home/groups/nmrl/utils/phylogenetics_tools/')
from subscripts.src.utilities import Housekeeper as hk
from process_chewbacca_profiles_edits import initialize_logging, gather_chewbacca_data, cluster_analysis, cluster_counter


###############
#Global configs
###############

cluster_counter = 0

###############


class Ardetype_housekeeper(hk):
    '''Class extends the standard housekeeper class to implement functions required by specific pipeline'''
        
    @staticmethod
    def merge_paths(src_list, target_folder=None, ignore_collisions=True, exclude_files=[]):
        '''
        Combines contents from all paths supplied as src_list under target_folder.
        Excludes and deletes files from source paths specified in exclude_files.
        '''

        def move_item(src, dst):
            """Move item from src to dst."""
            if os.path.exists(dst):
                raise ValueError(f"Collision detected at destination path '{dst}'")
            os.system(f'chmod -R 775 {src}')
            os.system(f'mv {src} {dst}')

        def delete_item(src):
            """Delete the specified item."""
            os.system(f'rm -rf {src}')

        if not target_folder:
            raise ValueError("target_folder argument is required")

        # Check for naming collisions among source directories and with the target directory
        src_name_to_paths = {}
        os.system(f'mkdir -m 775 -p {target_folder}')
        for src in src_list:
            if os.path.abspath(src) == os.path.abspath(target_folder):
                raise ValueError("Source and target folders cannot be the same.")
                
            for item in os.listdir(src):
                dest_path = os.path.join(target_folder, item)
                
                if not ignore_collisions:
                    if item in src_name_to_paths:
                        raise ValueError(f"Naming collision detected for '{item}' in sources '{src_name_to_paths[item]}' and '{src}'")
                    if os.path.exists(dest_path):
                        raise ValueError(f"Naming collision detected for '{item}' between source '{src}' and target '{target_folder}'")

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


    @staticmethod
    def type_fasta_scheme(contig_path: str, url: str) -> dict:
        '''
        Type contigs for single sample using pubmlst api. Returns a dictionary converted from json output of the HTTP request.
        If the HTTP request fails, raises and exception, indicating response code and response message (if any).
        The returned dictionary may not contain typing information, indicating that the typing has failed due to assembly quality/database issue/etc. 
        '''
        with open(contig_path, 'r') as x:  # read contigs from file
            fasta = x.read()
        payload = '{"base64":true,"details":true,"sequence":"' + base64.b64encode(
            fasta.encode()).decode() + '"}'  # construct payload of HTTP request
        # send the request and await the response
        response = requests.post(url, data=payload)
        if response.status_code == requests.codes.ok:  # response is valid if return code is valid
            data = response.json()  # get response payload into dictionary
            return data  # return path to fasta file and resulting dictionary
        else:  # if response is not valid
            raise Exception(f'''
            Typing has failed with the following status code: {response.status_code}.
            The following warning message was recieved: {response.text}
            ''')

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
        species_filter = ["GC (%)", "Total length", "N50", "contigs (>= 0 bp)", "assembly_coverageDepth"]
        genus_filter   = ["GC (%)", "Total length", "N50", "contigs (>= 0 bp)", "assembly_coverageDepth", "assembly_coverageDepth"]
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
                    if score != max_score:
                        check[f"{key} status"] = ['PASS']
                    else:
                        check[f"{key} status"] = ['FAIL']
                # if <= to pass
                if thr[key][0]['binscore'][-1] > thr[key][0]['binscore'][0]:
                    if score != 0:
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
        sample_id = re.sub(r'_S[0-9]*_resfinder', '',
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
        sample_id = re.sub(r'_S[0-9]*_resfinder_pheno.txt',
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
        sample_id = re.sub(r'_S[0-9]*_plasmidfinder', '',
                           os.path.basename(os.path.dirname(plf_result_path)))
        df = pd.read_csv(plf_result_path, sep='\t')
        df.insert(0, 'sample_id', [sample_id for _ in df.index])
        df.insert(1, 'analysis_batch_id', [os.path.basename(
            os.path.dirname(batch)) for _ in df.index])
        return df

    @staticmethod
    def mobtyper_results(mbt_result_path: str, batch: str) -> pd.DataFrame:
        '''To combine mob_typer reports and map them to sample_id-batch pair.'''
        sample_id = re.sub(r'(_S[0-9]*_mob_recon|_mob_recon)',
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
            r'_S[0-9]*_kraken2_contigs_report.txt', '', os.path.basename(k2c_report_path))
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
    def kraken2reads_results(k2r_report_path: str, batch: str) -> pd.DataFrame:
        '''To combine kraken2 read reports and map them to sample_id-batch pair.'''
        sample_id = re.sub(r'_S[0-9]*_kraken2_reads_report.txt',
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
        sample_id = re.sub(r'_S[0-9]*_quast', '', os.path.basename(os.path.dirname(qst_report_path)))
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

        df = pd.read_csv(klbt_result_path, sep='\t')
        df['strain'] = df['strain'].str.replace(
            r'_S[0-9]*_contigs', '', regex=True)
        # reorder columns and aggregate results even if resistance scan was not performed
        df = df.reindex(columns=columns, fill_value=None)
        df.insert(1, 'analysis_batch_id', [os.path.basename(
            os.path.dirname(batch)) for _ in df.index])
        return df

    @staticmethod
    def ectyper_results(ect_result_path: str, batch: str) -> pd.DataFrame:
        '''To combine ectyper reports and map them to sample_id-batch pair.'''
        df = pd.read_csv(ect_result_path, sep='\t')
        df.rename(columns={'Name': 'sample_id'}, inplace=True)
        df.sample_id = df.sample_id.str.replace(
            r'_S[0-9]*_contigs', '', regex=True)
        df.insert(1, 'analysis_batch_id', [os.path.basename(
            os.path.dirname(batch)) for _ in df.index])
        return df

    @staticmethod
    def stecfinder_results(stf_result_path: str, batch: str) -> pd.DataFrame:
        '''To combine stecfinder reports and map them to sample_id-batch pair.'''
        df = pd.read_csv(stf_result_path, sep='\t')
        df.rename(columns={'Sample': 'sample_id'}, inplace=True)
        df.sample_id = df.sample_id.str.replace(
            r'_S[0-9]*_bact_reads_classified', '', regex=True)
        df.insert(1, 'analysis_batch_id', [os.path.basename(
            os.path.dirname(batch)) for _ in df.index])
        return df

    @staticmethod
    def agrvate_results(agv_result_path: str, batch: str) -> pd.DataFrame:
        '''To combine agrvate reports and map them to sample_id-batch pair.'''
        df = pd.read_csv(agv_result_path, sep='\t')
        df.rename(columns={'#filename': 'sample_id'}, inplace=True)
        df.sample_id = df.sample_id.str.replace(
            r'_S[0-9]*_contigs', '', regex=True)
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
            r'_S[0-9]*_contigs.fasta', '', regex=True)
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
            r'_S[0-9]*_contigs.fasta', '', regex=True)
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
            r'_S[0-9]*_contigs.fasta', '', regex=True)
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
            r'_S[0-9]*_contigs', '', regex=True)
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
        sample_id = re.sub(r'_S[0-9]*_amrfinderplus_point.tab', '', os.path.basename(amrfpm_result_path))
        df.insert(0, 'sample_id', [sample_id for _ in df.index])
        df.insert(1, 'analysis_batch_id', [os.path.basename(os.path.dirname(batch)) for _ in df.index])
        return df

    @staticmethod
    def lrefinder_results(lrefinder_pos_path:str, batch: str):
        '''To combine lrefinder mutation reports and map them to sample_id-batch pair.'''
        df = pd.read_csv(lrefinder_pos_path, sep='\t')
        sample_id = re.sub(r'_S[0-9]*.pos', '', os.path.basename(lrefinder_pos_path))
        df.insert(0, 'sample_id', [sample_id for _ in df.index])
        df.insert(1, 'analysis_batch_id', [os.path.basename(os.path.dirname(batch)) for _ in df.index])
        return df
    
    @staticmethod
    def cgmlst_quality_results(lrefinder_pos_path:str, batch: str):
        '''To combine chewbbaca allele calling reports and map them to sample_id-batch pair.'''
        df = pd.read_csv(lrefinder_pos_path, sep='\t')
        df.rename(columns={'FILE': 'sample_id'}, inplace=True)
        df.sample_id = df.sample_id.str.replace(r'_S[0-9]*_contigs', '', regex=True)
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


    @staticmethod
    def cgmlst_cluster_analysis(path:str, wildcard:str="folded*/*_chewbbaca/results_alleles.tsv", thresholds:list=[], linkage:str='single', log_scale:bool=True, rem_outlier:int=None):
        current_datetime = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        output_path = os.path.join(path, 'cgmlst_clusters')
        os.makedirs(output_path, mode=0o775, exist_ok=True)
        initialize_logging(path = output_path, current_datetime = current_datetime)
        chew_dict = gather_chewbacca_data(path = path, wildcard = wildcard)
        cluster_analysis(
            path = path,
            output_path = output_path,
            thresholds = thresholds,
            df_dict=chew_dict,
            singularity=True, 
            linkage_method=linkage, 
            log_scale=log_scale, 
            outlier_filter=rem_outlier)


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
