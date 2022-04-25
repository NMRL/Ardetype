#!/usr/bin/env python3
# Modified from: https://pubmlst.org/species-id/species-identification-via-api
# Full Pasteur REST API: https://bigsdb.pasteur.fr/api
# API Docs: https://bigsdb.readthedocs.io/en/latest/rest.html#db-schemes-scheme-id-sequence

import requests, sys, argparse, base64, json

parser = argparse.ArgumentParser()
parser.add_argument('--file', '-f', type=str, default='contigs.fasta', help='assembly contig filename (FASTA format)')
args = parser.parse_args()

def main():
    uri = "https://bigsdb.pasteur.fr/api/db/pubmlst_listeria_seqdef/schemes/2/sequence"
    with open(args.file, 'r') as x: 
        fasta = x.read()
    req_start = '{"base64":true,"details":true,"sequence":"'
    req_stop = '"}'
    payload = req_start + base64.b64encode(fasta.encode()).decode() + req_stop
    response = requests.post(uri, data=payload)
    if response.status_code == requests.codes.ok:
        data = response.json()
        print(json.dumps(data, indent = 1)) #Can be redirected into file as needed
    else:
        print(response.status_code)
        print(response.text)

if __name__ == "__main__":
    main()
