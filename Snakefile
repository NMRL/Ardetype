configfile: 'config.yaml'

rule all:
    input: 
        expand('{sample}_mlst_output.csv', sample=config['samples']),
        expand('{sample}.assembly_stats.txt', sample=config['samples'])

rule quality_control:
    input:
        read_1 = expand('data/{sample}_R1_001.fastq.gz', sample=config['samples']),
        read_2 = expand('data/{sample}_R2_001.fastq.gz', sample=config['samples'])
    threads: 4
    output: 
        read_1_tr = expand('data/{sample}_fastp_R1_001.fastq.gz', sample=config['samples']),
        read_2_tr = expand('data/{sample}_fastp_R2_001.fastq.gz', sample=config['samples'])
    shell:
        'fastp --in1 {input.read_1} --in2 {input.read_2} --out1 {output.read_1_tr} --out2 {output.read_2_tr} --thread {threads}'

rule contig_assembly:
    input:
        read_1 = expand('data/{sample}_fastp_R1_001.fastq.gz', sample=config['samples']),
        read_2 = expand('data/{sample}_fastp_R2_001.fastq.gz', sample=config['samples']),
    output:
        temp(expand('contigs/contigs.fa', sample=config['samples']))
    threads: 4
    shell:
        'shovill --depth 100 --kmers 31,33,55,77,99,127 --ram 7 --minlen 500 --force --outdir contigs --R1 {input.read_1} --R2 {input.read_2}'

rule scaffold_assembly:
    input:
        reference = 'reference/GCF_900187225.1_51881_G01_genomic.fa',
        contigs =  expand("contigs/contigs.fa", sample=config['samples']),
    output:
        expand('scaffolds/ragtag.scaffold.fasta', sample=config['samples'])
    shell:
        'ragtag.py scaffold -o scaffolds -C {input.reference} {input.contigs}'

rule assembly_qc:
    input:
        scaffolds = expand('scaffolds/ragtag.scaffold.fasta', sample=config['samples'])
    output:
        stats = expand('{sample}.assembly_stats.txt', sample=config['samples'])
    shell:
        'statswrapper.sh in={input.scaffolds} > {output.stats}'

rule mlst:
    input: 
        scaffolds = 'scaffolds/ragtag.scaffold.fasta'
    output:
        mlst_output = expand('{sample}_mlst_output.csv', sample=config['samples'])
    threads: 4
    shell:
        'mlst --csv {input} >> {output.mlst_output}'

#using temp()
#adding sample id to contigs - bash script
#adding sample id to scaffolds - bash scripts
#finding resistance plasmids - Amr++