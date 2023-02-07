import sys, os, re, pandas as pd, subprocess as sp
from pathlib import Path
sys.path.insert(0, os.path.dirname(os.path.dirname(Path(__file__).absolute())))
from subscripts.src.utilities import Housekeeper as hk


class Ardetype_housekeeper(hk):
    '''Class extends the standard housekeeper class to implement functions required by specific pipeline'''

    @staticmethod
    def snakemake_to_dict(path_to_smk:str) -> dict:
        '''Given a path to a snakefile, returns dictionary mapping each rule to the code defined in "run" or "shell" directives as defined in given snakefile.'''
        with open(path_to_smk, 'r') as f:
            contents = f.read()

        names = re.findall(r'rule .*\:', contents)
        idxes = [contents.index(r) for r in names]
        rules = [contents[idx:idxes[i+1]] for i, idx in enumerate(idxes) if not idx == idxes[-1]]
        templates = []
        for k in rules:
            if 'shell:' in k: 
                templates.append(k[k.index('shell:'):].replace('shell:', '').strip())
            else:
                templates.append(k[k.index('run:'):].replace('run:', '').strip())

        names = [s.replace(':','').replace('rule ', '') for s in re.findall(r'rule .*\:', contents)]
        data = dict(zip(names, templates))
        return data


class Wrapper():
    '''Toolkit class to store wrapper methods for different tools'''

    #pipeline configuration saved at module import
    _rule_dict_core = Ardetype_housekeeper.snakemake_to_dict('./snakefiles/bact_core')
    _rule_dict_shell = Ardetype_housekeeper.snakemake_to_dict('./snakefiles/bact_shell')
    _rule_dict_tip = Ardetype_housekeeper.snakemake_to_dict('./snakefiles/bact_tip')
    _config_dict    = hk.read_yaml("./config_files/yaml/config_modular.yaml")
    
    #tool versions saved at module import
    _fastp_version  = sp.run(
        f'module load singularity && singularity run {_config_dict["fastp_sif"]} fastp --version',
        stderr=sp.PIPE, shell=True).stderr.decode('utf-8').strip()

    #stdout version extractor example
    # _cutadapt_version  = sp.run(
    #     f'module load singularity && singularity run {_config_dict["fastq_sif"]} cutadapt --version',
    #     stdout=sp.PIPE, shell=True).stdout.decode('utf-8').strip()
    
    #more extractors to be added
    

    @staticmethod
    def parse_fastp(sample_id:str, path_to_report:str) -> None:
        '''Serializes report as json file'''
        
        #Defining paths to required fastp MRIs
        mri_map = {
            'before_filtering':         ['summary','before_filtering'],
            'after_filtering':          ['summary', 'after_filtering'],
            'filtering_results':        ['filtering_result'],
            'duplication':              ['duplication', 'rate'],
            'insert_size':              ['insert_size'],
            'adapter_cutting':          ['adapter_cutting'],
            'read1_before_filtering':   ['read1_before_filtering'],
            'read2_before_filtering':   ['read2_before_filtering'],
            'read1_after_filtering':    ['read1_after_filtering'],
            'read2_after_filtering':    ['read2_after_filtering'],
            'command':                  ["command"]
        }

        #extracting MRIs from fastp report
        report = hk.read_json_dict(path_to_report)
        data   = {mri:hk.find_in_nested_dict(report, mri_map[mri]) for mri in mri_map} 

        tags = [
            'total_reads',
            'total_bases',
            'q20_bases',
            'q30_bases',
            'total_cycles',
            'quality_curves',
            'content_curves'
            ]
        filters = [
            'read1_before_filtering', 
            'read2_before_filtering', 
            'read1_after_filtering', 
            'read2_after_filtering'
            ]
        
        #keeping only subfields defined in tags
        for mri in filters:
            record    = dict(zip(tags,[data[mri][tag] for tag in tags]))
            data[mri] = record

        #generating sequencing mode string based on real cycles
        forward_cycles          = data["read1_before_filtering"]["total_cycles"]
        reverse_cycles          = data["read2_before_filtering"]["total_cycles"]
        seq_mode                = f'paired end ({forward_cycles} cycles + {reverse_cycles})'

        data['sequencing_mode'] = seq_mode
        data['fastp_version']   = Wrapper._fastp_version
        data['options']         = Wrapper._config_dict['fastp']
    
        output_folder           = os.path.abspath(os.path.dirname(path_to_report))
        outpath                 = f'{output_folder}/{sample_id}_fastp_et.json'
        hk.write_json(data, outpath, indent=4)        
