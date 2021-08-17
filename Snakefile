SAMPLES=['2104093601_S47_L001']

rule all:
    input: 
        expand('{sample}_mlst_output.csv', sample=SAMPLES),
        expand('{sample}.assembly_stats.txt', sample=SAMPLES)

rule quality_control:
    input:
        read_1 = expand('data/{sample}_R1_001.fastq.gz', sample=SAMPLES),
        read_2 = expand('data/{sample}_R2_001.fastq.gz', sample=SAMPLES)
    output: 
        read_1_tr = expand('data/{sample}_fastp_R1_001.fastq.gz', sample=SAMPLES),
        read_2_tr = expand('data/{sample}_fastp_R2_001.fastq.gz', sample=SAMPLES)
    threads: 4
    shell:
        'fastp --in1 {input.read_1} --in2 {input.read_2} --out1 {output.read_1_tr} --out2 {output.read_2_tr} --thread {threads}'

rule contig_assembly:
    input:
        read_1 = expand('data/{sample}_fastp_R1_001.fastq.gz', sample=SAMPLES),
        read_2 = expand('data/{sample}_fastp_R2_001.fastq.gz', sample=SAMPLES)
    output:
        contigs = expand('contigs/contigs.fa', sample=SAMPLES)
    threads: 4
    shell:
        'shovill --depth 100 --kmers 31,33,55,77,99,127 --ram 7 --minlen 500 --force --outdir contigs --R1 {input.read_1} --R2 {input.read_2}'

rule scaffold_assembly:
    input:
        reference = 'reference/GCF_900187225.1_51881_G01_genomic.fa',
        contigs =  expand("contigs/contigs.fa", sample=SAMPLES)
    output:
        scaffolds = '/ragtag_output/ragtag.scaffold.fasta'
    shell:
        'ragtag.py scaffold -C {input.reference} {input.contigs}'

rule assembly_qc:
    input:
        scaffolds = expand('ragtag_output/ragtag.scaffold.fasta', sample=SAMPLES)
    output:
        stats = expand('{sample}.assembly_stats.txt', sample=SAMPLES)
    shell:
        'statswrapper.sh in={input.scaffolds} > {output.stats}'

rule mlst:
    input: 
        scaffolds = 'ragtag_output/ragtag.scaffold.fasta'
    output:
        mlst_output = expand('{sample}_mlst_output.csv', sample=SAMPLES)
    threads: 4
    shell:
        'mlst --csv {input} >> {output.mlst_output}'

# How to add quast? - different python version required - explore qc options
# How to to clean-up? temp()
# How to input many files?