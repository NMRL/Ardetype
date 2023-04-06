import sys, os, requests, base64
from pathlib import Path
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


    