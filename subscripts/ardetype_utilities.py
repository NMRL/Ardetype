import sys, os, requests, base64, re, pathlib, pandas as pd
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor as ppe, as_completed
sys.path.insert(0, os.path.dirname(os.path.dirname(Path(__file__).absolute())))
from subscripts.src.utilities import Housekeeper as hk


class Ardetype_housekeeper(hk):
    '''Class extends the standard housekeeper class to implement functions required by specific pipeline'''

    @staticmethod
    def type_fasta_scheme(contig_path:str, url:str)->dict:
        '''
        Type contigs for single sample using pubmlst api. Returns a dictionary converted from json output of the HTTP request.
        If the HTTP request fails, raises and exception, indicating response code and response message (if any).
        The returned dictionary may not contain typing information, indicating that the typing has failed due to assembly quality/database issue/etc. 
        '''
        with open(contig_path, 'r') as x: #read contigs from file
            fasta = x.read()
        payload  = '{"base64":true,"details":true,"sequence":"' + base64.b64encode(fasta.encode()).decode() + '"}' #construct payload of HTTP request
        response = requests.post(url, data=payload) #send the request and await the response
        if response.status_code == requests.codes.ok: #response is valid if return code is valid
            data = response.json() #get response payload into dictionary
            return data #return path to fasta file and resulting dictionary
        else: #if response is not valid
            raise Exception(f'''
            Typing has failed with the following status code: {response.status_code}.
            The following warning message was recieved: {response.text}
            ''')


    @staticmethod
    def pointfinder_results(pf_results_path:str, batch:str):
        '''Returns dataframe containing sample_id and analysis_batch_id columns
        where the same sample id and batch id are matched to all findings for a given sample.'''
        df = pd.read_csv(pf_results_path, sep='\t')
        sample_id = re.sub(r'_S[0-9]*_resfinder', '', os.path.basename(os.path.dirname(pf_results_path)))
        df.insert(0, 'sample_id', [sample_id for _ in df.index])
        df.insert(1, 'analysis_batch_id', os.path.basename(os.path.dirname(batch)))
        df = df.astype(str)
        if len(df[df['PMID'].str.contains(', ')]) > 0:
            df['PMID'] = df['PMID'].str.replace(', ','')
        if len(df[df['Resistance'].str.contains(', ')]) > 0:
            df['Resistance'] = df['Resistance'].str.replace(', ','')
        return df
    

    @staticmethod
    def kleborate_results(klbt_result_path:str, batch:str):
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
        
        df           = pd.read_csv(klbt_result_path, sep='\t')
        df['strain'] = df['strain'].str.replace(r'_S[0-9]*_contigs', '', regex=True)
        df           = df.reindex(columns=columns, fill_value=None) #reorder columns and aggregate results even if resistance scan was not performed
        df.insert(1, 'analysis_batch_id', [os.path.basename(os.path.dirname(batch)) for _ in df.index])
        return df


    @staticmethod
    def respheno_results(rfp_result_path:str, batch:str) -> pd.DataFrame:
        '''Returns a dataframe - containing valid resistance genes.'''
        df = pd.read_csv(rfp_result_path, skiprows=16,sep="\t")
        sample_id = re.sub(r'_S[0-9]*_resfinder_pheno.txt', '', os.path.basename(rfp_result_path))
        df.insert(0, 'sample_id', [sample_id for _ in df.index])
        df.insert(1, 'analysis_batch_id', [os.path.basename(os.path.dirname(batch)) for _ in df.index])
        
        #remove all rows that contain "No resistance" in "WGS-predicted phenotype" column
        df  = df[df["WGS-predicted phenotype"] != "No resistance"]
        df  = df[~df['Genetic background'].isna()] #remove warnings
        return df
    

    @staticmethod
    def plasmidfinder_results(plf_result_path:str, batch:str) -> pd.DataFrame:
        sample_id = re.sub(r'_S[0-9]*_plasmidfinder', '', os.path.basename(os.path.dirname(plf_result_path)))
        df = pd.read_csv(plf_result_path, sep='\t')
        df.insert(0, 'sample_id', [sample_id for _ in df.index])
        df.insert(1, 'analysis_batch_id', [os.path.basename(os.path.dirname(batch)) for _ in df.index])
        return df


    @staticmethod
    def mobtyper_results(mbt_result_path:str, batch:str) -> pd.DataFrame:
        sample_id = re.sub(r'_S[0-9]*_mob_typer.tab', '', os.path.basename(mbt_result_path))
        df = pd.read_csv(mbt_result_path, sep='\t')
        df.rename(columns = {'sample_id': 'genetic_element'}, inplace=True)
        df.insert(0, 'sample_id', [sample_id for _ in df.index])
        df.insert(1, 'analysis_batch_id', [os.path.basename(os.path.dirname(batch)) for _ in df.index])
        return df


    @staticmethod
    def aggregator(outfolder_path:str, proc_num:int, wildcard:str=None, extractor=pointfinder_results, pathlist:list=None):
        '''Runs result extractors in-parallel - used in snakefiles for aggregation of specific results.'''
        if wildcard is not None:
            pf_out_list = pathlib.Path(outfolder_path).rglob(wildcard)
        elif pathlist is None and wildcard is None:
            raise Exception(f'Either "pathlist" or "wildcard" should be not None.')
        else:
            pf_out_list = pathlist
        summary = pd.DataFrame()
        with ppe(max_workers=proc_num) as executor:
            results = [executor.submit(extractor, path, outfolder_path) for path in pf_out_list]
            for output in as_completed(results):
                summary = pd.concat([summary, output.result()])
        summary = summary.reset_index(drop=True)
        return summary
