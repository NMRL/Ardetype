configfile: 'config.yaml'

import os
for id in config['samples']:
    os.system(f'mkdir -p {config["home_dir"]}Bact_assembly_output/{id}_output {config["home_dir"]}Bact_assembly_output/{id}_output/benchmarks')
    os.system(f'cp -u {config["input_dir"]}{id}* data/')

rule all:
    input: 
        #expand('{sample}_mlst_output.csv', sample=config['samples']),
        expand('{sample}.assembly_stats.txt', sample=config['samples']),
        expand('benchmarks/{sample}_combined_benchmark.csv', sample=config['samples'])
    run:
        for id in config['samples']:
            os.system(f'mv {id}* {config["home_dir"]}Bact_assembly_output/{id}_output/')
            os.system(f'mv data/{id}_fastp* {config["home_dir"]}Bact_assembly_output/{id}_output/')
            os.system(f'mv benchmarks/{id}* {config["home_dir"]}Bact_assembly_output/{id}_output/benchmarks')
            os.system(f'rm -r {id}*')
            os.system(f'rm data/{id}*')

rule quality_control:
    input:
        read_1 = 'data/{sample_id}_R1_001.fastq.gz',
        read_2 = 'data/{sample_id}_R2_001.fastq.gz'
    threads: 4
    benchmark:
        temp('benchmarks/{sample_id}.fastp.benchmark.txt')
    output: 
        temp('{sample_id}.fastp.json'),
        '{sample_id}.fastp.html',
        read_1_tr = 'data/{sample_id}_fastp_R1_001.fastq.gz',
        read_2_tr = 'data/{sample_id}_fastp_R2_001.fastq.gz'
    shell:
        'fastp -j {wildcards.sample_id}.fastp.json -h {wildcards.sample_id}.fastp.html --in1 {input.read_1} --in2 {input.read_2} --out1 {output.read_1_tr} --out2 {output.read_2_tr} --thread {threads}'

rule contig_assembly:
    input:
        read_1 = 'data/{sample_id}_fastp_R1_001.fastq.gz',
        read_2 = 'data/{sample_id}_fastp_R2_001.fastq.gz'
    output:
        temp('{sample_id}_contigs/contigs.fa')
    threads: 4
    benchmark:
        temp('benchmarks/{sample_id}.shovill.benchmark.txt')
    shell:
        'shovill --depth {config[shovill_params][depth]} --ram {config[shovill_params][ram]} --minlen {config[shovill_params][minlen]} --force --outdir {wildcards.sample_id}_contigs --R1 {input.read_1} --R2 {input.read_2}'

rule contig_id:
    input:
        'data/{sample_id}_fastp_R1_001.fastq.gz',
        cnt = '{sample_id}_contigs/contigs.fa'
    output:
        '{sample_id}_contigs/{sample_id}_contigs.fasta'
    shell:
        'cp {input.cnt} {output}'

rule scaffold_assembly:
    input:
        reference = config['reference'],
        contigs = '{sample_id}_contigs/{sample_id}_contigs.fasta'
    output:
        '{sample_id}_scaffolds/ragtag.scaffold.fasta'
    benchmark:
        temp('benchmarks/{sample_id}.ragtag.benchmark.txt')
    shell:
        'ragtag.py scaffold -o {wildcards.sample_id}_scaffolds -C {input.reference} {input.contigs}'

rule scaffold_id:
    input:
        'data/{sample_id}_fastp_R1_001.fastq.gz',
        scf = '{sample_id}_scaffolds/ragtag.scaffold.fasta'
    output:
        '{sample_id}_scaffolds/{sample_id}_ragtag.scaffold.fasta'
    shell:
        'cp {input.scf} {output}'

rule assembly_qc:
    input:
        'data/{sample_id}_fastp_R1_001.fastq.gz',
        scaffolds = '{sample_id}_scaffolds/{sample_id}_ragtag.scaffold.fasta'
    output:
        stats = '{sample_id}.assembly_stats.txt'
    benchmark:
        temp('benchmarks/{sample_id}.bbmap_qc.benchmark.txt')
    shell:
        'statswrapper.sh in={input.scaffolds} > {output.stats}'

# rule mlst:
#     input: 
#         'data/{sample_id}_fastp_R1_001.fastq.gz',
#         scaffolds = '{sample_id}_scaffolds/{sample_id}_ragtag.scaffold.fasta'
#     output:
#         mlst_output = '{sample_id}_mlst_output.csv',
#     threads: 4
#     benchmark:
#         temp('benchmarks/{sample_id}.mlst.benchmark.txt')
#     shell:
#         'mlst --csv {input} >> {output.mlst_output}'

rule combine_benchmark:
    input:
        'benchmarks/{sample_id}.bbmap_qc.benchmark.txt',
        'benchmarks/{sample_id}.fastp.benchmark.txt',
        #'benchmarks/{sample_id}.mlst.benchmark.txt',
        'benchmarks/{sample_id}.ragtag.benchmark.txt',
        'benchmarks/{sample_id}.shovill.benchmark.txt'
    output:
        'benchmarks/{sample_id}_combined_benchmark.csv'
    shell:
        'python scripts/combine_benchmarks.py -b benchmarks -s {wildcards.sample_id}'

# on a cluster - https://carpentries-incubator.github.io/workflows-snakemake/09-cluster/index.html