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

    #stdout extractor example
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


    @staticmethod
    def combine_jsons(sample_id:str, output_path:str, report_paths:list) -> None:
        '''Per-sample aggregator resulting in a single json file containing information from all tools'''
        _template_path = Wrapper._config_dict['sample_template']                              #path to json template to use
        _template      = hk.read_json_dict(_template_path)                                    #reading in json template
        data = {}

        for path in report_paths:
            report_name = os.path.basename(path)

            #infer corresponding read number from serialized report naming convention
            if   '_1_et.json' in report_name:   read = 1                                      
            elif '_2_et.json' in report_name:   read = 2                                      
            else:                               read = ''                                     #non-read based report         

            report_type = report_name.replace(f'{sample_id}_','')
            key         = re.sub('(_[0-9])?_et.json', '', report_type)                        #tool name
            results     = hk.read_json_dict(path)                                             #serialized results

            conditions = [                                                                    #for processing separate reports for each read file
                not (key in data) and not (read == ''),                                       #read-based report first time
                key in data and read != ''                                                    #read-based report second time
            ]
 
            if   conditions[0]: data[key]                  = {f'read_{read}':results}
            elif conditions[1]: data[key][f'read_{read}']  = results
            else:               data[key]                  = results

        #add versions for all tools not explicitly serialized or resulting in separate report for each read
        tool_dict = {
            #    'cutadapt': [Wrapper._rule_dict['adapter_removal'].split('\n'),     Wrapper._cutadapt_version,  Wrapper._config_dict['cutadapt']],
            #     'bwa-mem': [Wrapper._rule_dict['read_alignment'].split('\n'),      Wrapper._bwamem_version,    Wrapper._config_dict['bwa-mem']],
            #      'picard': [Wrapper._rule_dict['indel_realignment'].split('\n'),   Wrapper._picard_version,    Wrapper._config_dict['picard']],
            #        'abra': [Wrapper._rule_dict['indel_realignment'].split('\n'),   Wrapper._abra_version,      Wrapper._config_dict['abra']],
            #      'snpEff': [Wrapper._rule_dict['variant_annotation'].split('\n'),  Wrapper._snpeff_version,    Wrapper._config_dict['snpEff']],
            #   'freebayes': [Wrapper._rule_dict['variant_calling'].split('\n'),     Wrapper._freebayes_version, Wrapper._config_dict['freebayes']],
            # 'fastqscreen': [Wrapper._rule_dict['fastq_screening'].split('\n'),     Wrapper._fastqc_version,    Wrapper._config_dict['fastqscreen']]
        }

        for tool in tool_dict:
            if tool in data:                                            #not explicitly serialized
                data[tool][f'{tool}_version'] =     tool_dict[tool][1]
                data[tool]['command']         =     tool_dict[tool][0]
                data[tool]['options']         =     tool_dict[tool][2]
            else:                                                       #resulting in separate report for each read
                data[tool] = {
                    f'{tool}_version' : tool_dict[tool][1], 
                    'command'         : tool_dict[tool][0], 
                    'options'         : tool_dict[tool][2]
                    }

        #populate the template
        for key in data: hk.edit_nested_dict(_template['data']['pipelines']['inhouse-ardetype'], key, data[key])

        outpath = f'{output_path}/{sample_id}_combined.json'
        hk.write_json(_template, outpath)