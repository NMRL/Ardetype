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
    benchmark:
        'benchmarks/{sample_id}.fastp.benchmark.txt'
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
        'benchmarks/{sample_id}.shovill.benchmark.txt'
    shell:
        'shovill --depth {config[shovill_params][depth]} --kmers {config[shovill_params][k1]},{config[shovill_params][k2]},{config[shovill_params][k3]},{config[shovill_params][k4]},{config[shovill_params][k5]},{config[shovill_params][k6]} --ram {config[shovill_params][ram]} --minlen {config[shovill_params][minlen]} --force --outdir {wildcards.sample_id}_contigs --R1 {input.read_1} --R2 {input.read_2}'

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
    benchmark:
        'benchmarks/{sample_id}.ragtag.benchmark.txt'
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
        'benchmarks/{sample_id}.bbmap_qc.benchmark.txt'
    shell:
        'statswrapper.sh in={input.scaffolds} > {output.stats}'

rule mlst:
    input: 
        'data/{sample_id}_fastp_R1_001.fastq.gz',
        scaffolds = '{sample_id}_scaffolds/{sample_id}_ragtag.scaffold.fasta'
    output:
        mlst_output = '{sample_id}_mlst_output.csv',
    threads: 4
    benchmark:
        'benchmarks/{sample_id}.mlst.benchmark.txt'
    shell:
        'mlst --csv {input} >> {output.mlst_output}'

# on a cluster - https://carpentries-incubator.github.io/workflows-snakemake/09-cluster/index.html