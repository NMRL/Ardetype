configfile: 'config.yaml'
import os, pandas as pd, shlex, subprocess

#GODLIKE EXPLAINATION HOW SNAKEMAKE WORKS
#https://vincebuffalo.com/blog/2020/03/04/understanding-snakemake.html


#CONFIGURATIONS
# os.system('module load singularity')

#READ KMERFINDER INPUT
kmerfinder_output = pd.read_csv('01_13_2022_assambled_kmerfinder_report.csv')
samples = list(kmerfinder_output['sample_id'])
reference_list = list(kmerfinder_output['accession'])
snakemake_sif = 'snakemake.sif'

#GENERATE TARGET LIST
target_list, scaffold = [], '{sample_id_pattern}-{reference_sequence_pattern}-scaffolds/{sample_id_pattern}-{reference_sequence_pattern}-ragtag.scaffold.fasta'
for sample_id_pattern, reference_sequence_pattern in zip(samples, reference_list):
    prokka = f'{sample_id_pattern}-prokka'
    rgi_txt = f'{sample_id_pattern}.rgi.txt'
    rgi_json = f'{sample_id_pattern}.rgi.json'
    mlst = f'{sample_id_pattern}-{reference_sequence_pattern}-scaffolds/{sample_id_pattern}-{reference_sequence_pattern}_mlst_output.csv'
    target_list+=[mlst]
    target_list+=[prokka]
    target_list+=[rgi_txt]
    target_list+=[rgi_json]

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
        sif_file = snakemake_sif,
        read_1 = 'data/{sample_id_pattern}_R1_001.fastq.gz',
        read_2 = 'data/{sample_id_pattern}_R2_001.fastq.gz'
    threads: 4
    envmodules:
        'singularity'
    benchmark:
        temp('benchmarks/{sample_id_pattern}.fastp.benchmark.txt')
    output: 
        temp('{sample_id_pattern}.fastp.json'),
        '{sample_id_pattern}.fastp.html',
        read_1_tr = 'data/{sample_id_pattern}_fastp_R1_001.fastq.gz',
        read_2_tr = 'data/{sample_id_pattern}_fastp_R2_001.fastq.gz'
    shell:
        'singularity run {input.sif_file} fastp -j {wildcards.sample_id_pattern}.fastp.json -h {wildcards.sample_id_pattern}.fastp.html --in1 {input.read_1} --in2 {input.read_2} --out1 {output.read_1_tr} --out2 {output.read_2_tr} --thread {threads}'

#GENERATING CONTIGS FROM READS
rule contig_assembly:
    input:
        sif_file = snakemake_sif,
        read_1 = 'data/{sample_id_pattern}_R1_001.fastq.gz',
        read_2 = 'data/{sample_id_pattern}_R2_001.fastq.gz'
    output:
        temp('{sample_id_pattern}_contigs/contigs.fa')
    envmodules:
        'singularity'
    threads: 4
    benchmark:
        temp('benchmarks/{sample_id_pattern}.shovill.benchmark.txt')
    shell:
        'singularity run {input.sif_file} shovill --depth {config[shovill_params][depth]} --ram {config[shovill_params][ram]} --minlen {config[shovill_params][minlen]} --force --outdir {wildcards.sample_id_pattern}_contigs --R1 {input.read_1} --R2 {input.read_2}'

#RENAMING CONTINGS
rule contig_id:
    input:
        'data/{sample_id_pattern}_fastp_R1_001.fastq.gz',
        cnt = '{sample_id_pattern}_contigs/contigs.fa'
    envmodules:
        'singularity'
    output:
        '{sample_id_pattern}_contigs/{sample_id_pattern}_contigs.fasta'
    shell:
        'cp {input.cnt} {output}'

#APPLYING CONTIG CORRECTIONS
rule contig_correction:
    input:
        sif_file = snakemake_sif,
        reference = 'reference/{reference_sequence_pattern}.fasta',
        contigs = '{sample_id_pattern}_contigs/{sample_id_pattern}_contigs.fasta'
    envmodules:
        'singularity'
    output:
        '{sample_id_pattern}-{reference_sequence_pattern}-scaffolds/ragtag.correct.fasta'
    shell:
        "singularity run {input.sif_file} ragtag.py correct {input.reference} {input.contigs} -o {wildcards.sample_id_pattern}-{wildcards.reference_sequence_pattern}-scaffolds"

#GENERATING SCAFFOLDS USING EXISTING GENOME AS REFERENCE 
#IF SCAFFOLDING FAILS, WHICH IS THE CASE FOR POOR REFERENCE MATCH - GENERATE EMPTY PLACEHOLDER FILE SO THAT PIPELINE IS NOT BROKEN
rule scaffold_assembly:
    input:
        sif_file = snakemake_sif,
        reference = 'reference/{reference_sequence_pattern}.fasta',
        contigs = '{sample_id_pattern}-{reference_sequence_pattern}-scaffolds/ragtag.correct.fasta'
    envmodules:
        'singularity'
    output:
        '{sample_id_pattern}-{reference_sequence_pattern}-scaffolds/ragtag.scaffold.fasta'
    benchmark:
        temp('benchmarks/{sample_id_pattern}-{reference_sequence_pattern}.ragtag.benchmark.txt')
    shell:
        """
        singularity run {input.sif_file} ragtag.py scaffold -o {wildcards.sample_id_pattern}-{wildcards.reference_sequence_pattern}-scaffolds -C {input.reference} {input.contigs}
        if [[ -f {sample_id_pattern}-{reference_sequence_pattern}-scaffolds/ragtag.scaffold.fasta ]] ; then echo 'Scaffolding succesful!' ; else touch {output} ; fi 
        """

#RENAMING SCAFFOLDS
rule scaffold_id:
    input:
        scf = '{sample_id_pattern}-{reference_sequence_pattern}-scaffolds/ragtag.scaffold.fasta'
    envmodules:
        'singularity'
    output:
        '{sample_id_pattern}-{reference_sequence_pattern}-scaffolds/{sample_id_pattern}-{reference_sequence_pattern}-ragtag.scaffold.fasta'
    shell:
        'cp {input.scf} {output}'

#GENERATING QUALITY CONTROL METRICS FOR SCAFFOLDS
rule assembly_qc:
    input:
        sif_file = snakemake_sif,
        scaffolds = '{sample_id_pattern}-{reference_sequence_pattern}-scaffolds/{sample_id_pattern}-{reference_sequence_pattern}-ragtag.scaffold.fasta'
    output:
        stats = '{sample_id_pattern}-{reference_sequence_pattern}.assembly_stats.txt'
    envmodules:
        'singularity'
    benchmark:
        temp('benchmarks/{sample_id_pattern}-{reference_sequence_pattern}.bbmap_qc.benchmark.txt')
    shell:
        'singularity run {input.sif_file} statswrapper.sh in={input.scaffolds} > {output.stats}'

#PERFORM SEROTYPING
rule mlst: 
    input: 
        mlst_quast_path = 'mlst_quast.sif',
        scaffolds = '{sample_id_pattern}-{reference_sequence_pattern}-scaffolds/{sample_id_pattern}-{reference_sequence_pattern}-ragtag.scaffold.fasta'
    output:
        mlst_output = '{sample_id_pattern}-{reference_sequence_pattern}-scaffolds/{sample_id_pattern}-{reference_sequence_pattern}_mlst_output.csv',
    threads: 4
    envmodules:
        'singularity'
    benchmark:
        temp('benchmarks/{sample_id_pattern}-{reference_sequence_pattern}.mlst.benchmark.txt')
    shell:
        """
        singularity run {input.mlst_quast_path} mlst --csv {input.scaffolds} >> {output.mlst_output}
        """

# #PERFORM HOST FILTERING
# rule kraken2: 
#     input: 
#         scaffolds = '{sample_id_pattern}-{reference_sequence_pattern}-scaffolds/{sample_id_pattern}-{reference_sequence_pattern}-ragtag.scaffold.fasta'
#     output:
#         mlst_output = '{sample_id_pattern}-{reference_sequence_pattern}-scaffolds/{sample_id_pattern}-{reference_sequence_pattern}_mlst_output.csv',
#     threads: 4
#     benchmark:
#         temp('benchmarks/{sample_id_pattern}-{reference_sequence_pattern}.mlst.benchmark.txt')
#     shell:
#         """
#         singularity run {input.mlst_quast_path} mlst --csv {input.scaffolds} >> {output.mlst_output}
#         """


#PERFORM GENE ANNOTATION
rule prokka: 
    input: 
        rgi_prokka_path = 'rgi_prokka.sif',
        contigs = '{sample_id_pattern}_contigs/{sample_id_pattern}_contigs.fasta'
    output:
        prokka_outdir = directory('{sample_id_pattern}-prokka')
    threads: 4
    envmodules:
        'singularity'
    benchmark:
        temp('benchmarks/{sample_id_pattern}.prokka.benchmark.txt')
    shell:
        """singularity run {input.rgi_prokka_path} prokka --outdir {output.prokka_outdir} {input.contigs}"""
        


#PERFORM RESISTANCE GENE EXTRACTION
rule res_gen_id: 
    input: 
        contigs = '{sample_id_pattern}_contigs/{sample_id_pattern}_contigs.fasta'
    output:
        rgi_txt_output = '{sample_id_pattern}.rgi.txt',
        rgi_json_output = '{sample_id_pattern}.rgi.json'
    threads: 4
    conda:
        'rgi_env.yaml'
    benchmark:
        temp('benchmarks/{sample_id_pattern}.rig.benchmark.txt')
    shell:
        """
        rgi main --input_sequence {input.contigs} --output_file {wildcards.sample_id_pattern}.rgi --input_type contig --clean
        """
# if [[ -f {sample_id_pattern}.rgi ]] ; then echo 'Gene search succesful!' ; else touch {output.rgi_output} ; fi 
# on a cluster - https://carpentries-incubator.github.io/workflows-snakemake/09-cluster/index.html
#Add KRAKEN2 HOST FILTERING (USE GH38 human assembly)
#Add prokka to extract genes
#Add quast to automate qc (for each reference + metaquast based on all references at once (using contigs))
#Fix benchmarking
#Add mlst assembly script