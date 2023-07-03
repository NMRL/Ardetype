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
        '''Returns dataframe containing sample_id column where the same sample id is matched to all findings for a given sample.'''
        df = pd.read_csv(pf_results_path, sep='\t')
        sample_id = re.sub(r'_S[0-9]*_resfinder', '', os.path.basename(os.path.dirname(pf_results_path)))
        df.insert(0, 'sample_id', [sample_id for _ in df.index])
        df.insert(1, 'seq_batch', os.path.basename(os.path.dirname(batch)))
        df = df[df.columns[::-1]] #reverse column order
        df = df.astype(str)
        if len(df[df['PMID'].str.contains(', ')]) > 0:
            df['PMID'] = df['PMID'].str.replace(', ','')
        if len(df[df['Resistance'].str.contains(', ')]) > 0:
            df['Resistance'] = df['Resistance'].str.replace(', ','')
        return df


    @staticmethod 
    def aggregate_pointfinder(outfolder_path:str, proc_num:int):
        '''Runs pointfinder_results '''

        pf_out_list = pathlib.Path(outfolder_path).rglob("*/PointFinder_results.txt")
        summary = pd.DataFrame()
        with ppe(max_workers=proc_num) as executor:
            results = [executor.submit(Ardetype_housekeeper.pointfinder_results, path, outfolder_path) for path in pf_out_list]
            for output in as_completed(results):
                summary = pd.concat([summary, output.result()])
        summary = summary.reset_index(drop=True)
        return summary