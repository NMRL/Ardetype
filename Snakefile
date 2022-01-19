configfile: 'config.yaml'
import os, pandas as pd

#GODLIKE EXPLAINATION HOW SNAKEMAKE WORKS
#https://vincebuffalo.com/blog/2020/03/04/understanding-snakemake.html

#READ KMERFINDER INPUT
kmerfinder_output = pd.read_csv('01_13_2022_assambled_kmerfinder_report.csv')
samples = list(kmerfinder_output['sample_id'])
reference_list = list(kmerfinder_output['accession'])

#GENERATE TARGET LIST
target_list, scaffold = [], '{sample_id_pattern}-{reference_sequence_pattern}-scaffolds/{sample_id_pattern}-{reference_sequence_pattern}-ragtag.scaffold.fasta'
for sample_id_pattern, reference_sequence_pattern in zip(samples, reference_list):
    benchmark = f'{sample_id_pattern}-{reference_sequence_pattern}-scaffolds/{sample_id_pattern}-{reference_sequence_pattern}-ragtag.scaffold.fasta'
    target_list+=[benchmark]

#GENERATING FOLDERS FOR EACH SAMPLE-REF COMBINATION SPECIFIED IN THE KMERFINDER OUTPUT
for i in range(len(samples)):
    os.system(f'mkdir -p {config["home_dir"]}{samples[i]}_output {config["home_dir"]}{samples[i]}_output/benchmarks')

#FINAL RULE
rule all:
    input: 
        target_list
    run:
        for i in range(len(samples)):
            os.system(f'python scripts/combine_benchmarks.py -b benchmarks -s {samples[i]}')
            os.system(f'mv {samples[i]}* {config["home_dir"]}{samples[i]}_output/')
            os.system(f'mv data/{samples[i]}_fastp* {config["home_dir"]}{samples[i]}_output/')
            os.system(f'mv benchmarks/{samples[i]}* {config["home_dir"]}{samples[i]}_output/benchmarks/')

#READ LENGTH & QUALITY TRIMMING
rule quality_control:
    input:
        read_1 = 'data/{sample_id_pattern}_R1_001.fastq.gz',
        read_2 = 'data/{sample_id_pattern}_R2_001.fastq.gz'
    threads: 4
    benchmark:
        temp('benchmarks/{sample_id_pattern}.fastp.benchmark.txt')
    output: 
        temp('{sample_id_pattern}.fastp.json'),
        '{sample_id_pattern}.fastp.html',
        read_1_tr = 'data/{sample_id_pattern}_fastp_R1_001.fastq.gz',
        read_2_tr = 'data/{sample_id_pattern}_fastp_R2_001.fastq.gz'
    shell:
        'fastp -j {wildcards.sample_id_pattern}.fastp.json -h {wildcards.sample_id_pattern}.fastp.html --in1 {input.read_1} --in2 {input.read_2} --out1 {output.read_1_tr} --out2 {output.read_2_tr} --thread {threads}'

#GENERATING CONTIGS FROM READS
rule contig_assembly:
    input:
        read_1 = 'data/{sample_id_pattern}_R1_001.fastq.gz',
        read_2 = 'data/{sample_id_pattern}_R2_001.fastq.gz'
    output:
        temp('{sample_id_pattern}_contigs/contigs.fa')
    threads: 4
    benchmark:
        temp('benchmarks/{sample_id_pattern}.shovill.benchmark.txt')
    shell:
        'shovill --depth {config[shovill_params][depth]} --ram {config[shovill_params][ram]} --minlen {config[shovill_params][minlen]} --force --outdir {wildcards.sample_id_pattern}_contigs --R1 {input.read_1} --R2 {input.read_2}'

#RENAMING CONTINGS
rule contig_id:
    input:
        'data/{sample_id_pattern}_fastp_R1_001.fastq.gz',
        cnt = '{sample_id_pattern}_contigs/contigs.fa'
    output:
        '{sample_id_pattern}_contigs/{sample_id_pattern}_contigs.fasta'
    shell:
        'cp {input.cnt} {output}'

#APPLYING CONTIG CORRECTIONS
rule contig_correction:
    input:
        reference = 'reference/{reference_sequence_pattern}.fasta',
        contigs = '{sample_id_pattern}_contigs/{sample_id_pattern}_contigs.fasta'
    output:
        '{sample_id_pattern}-{reference_sequence_pattern}-scaffolds/ragtag.correct.fasta'
    shell:
        "ragtag.py correct {input.reference} {input.contigs} -o {wildcards.sample_id_pattern}-{wildcards.reference_sequence_pattern}-scaffolds"

#GENERATING SCAFFOLDS USING EXISTING GENOME AS REFERENCE
rule scaffold_assembly:
    input:
        reference = 'reference/{reference_sequence_pattern}.fasta',
        contigs = '{sample_id_pattern}-{reference_sequence_pattern}-scaffolds/ragtag.correct.fasta'
    output:
        '{sample_id_pattern}-{reference_sequence_pattern}-scaffolds/ragtag.scaffold.fasta'
    benchmark:
        temp('benchmarks/{sample_id_pattern}-{reference_sequence_pattern}.ragtag.benchmark.txt')
    shell:
        "ragtag.py scaffold -o {wildcards.sample_id_pattern}-{wildcards.reference_sequence_pattern}-scaffolds -C {input.reference} {input.contigs}"

#RENAMING SCAFFOLDS
rule scaffold_id:
    input:
        scf = '{sample_id_pattern}-{reference_sequence_pattern}-scaffolds/ragtag.scaffold.fasta'
    output:
        '{sample_id_pattern}-{reference_sequence_pattern}-scaffolds/{sample_id_pattern}-{reference_sequence_pattern}-ragtag.scaffold.fasta'
    shell:
        'cp {input.scf} {output}'

#GENERATING QUALITY CONTROL METRICS FOR SCAFFOLDS
rule assembly_qc:
    input:
        # 'data/{sample_id}_fastp_R1_001.fastq.gz',
        scaffolds = '{sample_id_pattern}-{reference_sequence_pattern}-scaffolds/{sample_id_pattern}-{reference_sequence_pattern}-ragtag.scaffold.fasta'
    output:
        stats = '{sample_id_pattern}-{reference_sequence_pattern}.assembly_stats.txt'
    benchmark:
        temp('benchmarks/{sample_id_pattern}-{reference_sequence_pattern}.bbmap_qc.benchmark.txt')
    shell:
        'statswrapper.sh in={input.scaffolds} > {output.stats}'

# rule mlst:
#     input: 
#         'data/{sample_id}_fastp_R1_001.fastq.gz',
#         scaffolds = '{sample_id}-scaffolds/{sample_id}_ragtag.scaffold.fasta'
#     output:
#         mlst_output = '{sample_id}_mlst_output.csv',
#     threads: 4
#     benchmark:
#         temp('benchmarks/{sample_id}.mlst.benchmark.txt')
#     shell:
#         'mlst --csv {input} >> {output.mlst_output}'

   


# # on a cluster - https://carpentries-incubator.github.io/workflows-snakemake/09-cluster/index.html
