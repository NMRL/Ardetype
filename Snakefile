configfile: 'config.yaml'

rule all:
    input: 
        expand('{sample}_mlst_output.csv', sample=config['samples']),
        expand('{sample}.assembly_stats.txt', sample=config['samples'])

rule quality_control:
    input:
        read_1 = 'data/{sample_id}_R1_001.fastq.gz',
        read_2 = 'data/{sample_id}_R2_001.fastq.gz'
    threads: 4
    output: 
        read_1_tr = 'data/{sample_id}_fastp_R1_001.fastq.gz',
        read_2_tr = 'data/{sample_id}_fastp_R2_001.fastq.gz'
    shell:
        'fastp --in1 {input.read_1} --in2 {input.read_2} --out1 {output.read_1_tr} --out2 {output.read_2_tr} --thread {threads}'

rule contig_assembly:
    input:
        read_1 = 'data/{sample_id}_fastp_R1_001.fastq.gz',
        read_2 = 'data/{sample_id}_fastp_R2_001.fastq.gz'
    output:
        temp('{sample_id}_contigs/contigs.fa')
    threads: 4
    shell:
        'shovill --depth 100 --kmers 31,33,55,77,99,127 --ram 7 --minlen 500 --force --outdir {wildcards.sample_id}_contigs --R1 {input.read_1} --R2 {input.read_2}'

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
        reference = 'reference/GCF_900187225.1_51881_G01_genomic.fa',
        contigs = '{sample_id}_contigs/{sample_id}_contigs.fasta'
    output:
        '{sample_id}_scaffolds/ragtag.scaffold.fasta'
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
    shell:
        'statswrapper.sh in={input.scaffolds} > {output.stats}'

rule mlst:
    input: 
        'data/{sample_id}_fastp_R1_001.fastq.gz',
        scaffolds = '{sample_id}_scaffolds/{sample_id}_ragtag.scaffold.fasta'
    output:
        mlst_output = '{sample_id}_mlst_output.csv',
    threads: 4
    shell:
        'mlst --csv {input} >> {output.mlst_output}'

#using temp()
#adding sample id to contigs - bash script
#adding sample id to scaffolds - bash scripts
#finding resistance plasmids - Amr++