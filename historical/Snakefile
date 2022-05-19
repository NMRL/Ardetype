configfile: 'config.yaml'
import os, pandas as pd, shlex, subprocess, shutil

#GODLIKE EXPLAINATION OF HOW SNAKEMAKE WORKS
#https://vincebuffalo.com/blog/2020/03/04/understanding-snakemake.html

#STATIC PATHS - TO BE MOVED TO CONFIG
bact_analysis_path = f'/home/groups/nmrl/bact_analysis/'
kfinder = f'{bact_analysis_path}kmerfinder_suite/'
kfinder_input = f'{kfinder}input_files/'
kfinder_output = f'{kfinder}output/'
kfinder_backup = f'{kfinder}backup/'
home_pipe = f'{bact_analysis_path}NMRL_Bact_Assembly_Inhouse/'
ref_db = f'/home/groups/nmrl/db/db-refseq/'
snakemake_sif = '/home/groups/nmrl/image_files/snakemake.sif'
mlst_quast_sif = '/home/groups/nmrl/image_files/mlst_quast.sif'
prokka_sif = '/home/groups/nmrl/image_files/prokka.sif'
aquamis = f'{bact_analysis_path}AQUAMIS/'
aquamis_scripts = f'{aquamis}scripts/'
kraken2_env = f"{home_pipe}/conda_defs/kraken2.yaml"
rgi_env = f'{home_pipe}/conda_defs/rgi_env.yaml'


#DEV FLAGS
run_kfinder = True

#RUN KMERFINDER
if run_kfinder:
    #Check kmerfinder folders - move files out to backup if any left from previous runs; purge output
    if any(os.scandir(kfinder_input)):
        os.system(f'mv {kfinder_input}* {kfinder_backup}')
    if any(os.scandir(kfinder_output)):
        os.system(f'rm -r {kfinder_output}*')
    # Move files to kmerfinder folder
    os.system(f'mv {config["target_dir"]}* {kfinder_input}')

    # run kmerfinder (abort if failed; wait until finished)
    try:
        os.chdir(kfinder)
        subprocess.check_call(f'python run_kmerfinder.py -n {config["kfinder_threads"]} -d {config["kfinder_db"]}'.split(" "))
        os.system(f'mv  {kfinder_input}* {config["target_dir"]}')
        kfinder_table = [f for f in os.listdir(kfinder_output) if 'assambled_kmerfinder_report.csv' in f][0]
    except Exception as error:
        os.system(f'mv  {kfinder_input}* {config["target_dir"]}')
        sys.exit(f'kmerfinder error: {error}')
else:
    kfinder_table = [f for f in os.listdir(kfinder_output) if 'assambled_kmerfinder_report.csv' in f][0]

#READ KMERFINDER INPUT
kmerfinder_output = pd.read_csv(f'{kfinder_output}{kfinder_table}')
samples = list(kmerfinder_output["sample_id"])
reference_list = list(kmerfinder_output["accession"])


#GET NEW REFERENCES
os.chdir(ref_db)
ref_list = os.listdir(f'{home_pipe}reference/')
for acc in kmerfinder_output['accession']:
    if f'{acc}.fasta' not in ref_list:
        df = pd.DataFrame({"accession":[acc]})
        df.to_csv(f'{ref_db}{acc}.csv', header=True, index=False)
        print(f'{acc}: Attempting download')
        subprocess.call(f'python manage_refseqdb.py --add-ncbi {acc}'.split(' ')) #download sequence 
        print(f'{acc}: Attempting export')
        subprocess.call(f'python manage_refseqdb.py --exp-fasta {acc}.csv'.split(' ')) #extract sequence
        try:
            shutil.move(f'{ref_db}exported_sequences.fasta', f'{home_pipe}reference/{acc}.fasta') # rename
        except Exception as e:
            print(f'{acc} export failed: {e}')
        os.remove(f'{ref_db}{acc}.csv')

#GENERATE AQUAMIS LIST
os.chdir(aquamis_scripts)
subprocess.call(f'bash create_sampleSheet.sh --mode illumina --fastxDir {config["target_dir"]} --outDir {config["target_dir"]}'.split(' '))
sample_list_path = f'{config["target_dir"]}samples.tsv'

#INIT SAMPLE ID & REFERENCE_SEQUENCE_PATTERN WILDCARDS TO BE USED IN RULES
scaffold = '{sample_id_pattern}-{reference_sequence_pattern}-scaffolds/{sample_id_pattern}-{reference_sequence_pattern}-ragtag.scaffold.fasta'

#INIT TARGET LIST
target_list = []

for sample_id_pattern, reference_sequence_pattern in zip(samples, reference_list):
    #TARGET FILE NAMES BASED ON SAMPLE IDS ARE CREATED HERE
    prokka = f'{sample_id_pattern}-prokka'
    rgi_txt = f'{sample_id_pattern}.rgi.txt'
    rgi_json = f'{sample_id_pattern}.rgi.json'
    mlst = f'{sample_id_pattern}-{reference_sequence_pattern}-scaffolds/{sample_id_pattern}-{reference_sequence_pattern}_mlst_output.csv'
    kraken2 = f'{sample_id_pattern}_contigs/{sample_id_pattern}_kraken2_report.txt'
    quast = f'{sample_id_pattern}-{reference_sequence_pattern}-scaffolds/icarus.html'

    #TARGET FILES BEING COMBINED IN A LIST
    target_list += [
        mlst,
        prokka,
        rgi_txt,
        rgi_json,
        kraken2,
        quast
        ]


#GENERATING FOLDERS FOR EACH SAMPLE-REF COMBINATION SPECIFIED IN THE KMERFINDER OUTPUT
os.system(f'mkdir -p {config["home_dir"]}{config["target_dir"].split("/")[-2]}_output')
os.system(f'mkdir -p {config["home_dir"]}{config["target_dir"].split("/")[-2]}_output/kmerfinder_output/')
os.system(f'mkdir -p {config["home_dir"]}{config["target_dir"].split("/")[-2]}_output/aquamis_output/')
for i in range(len(samples)):
    os.system(f'mkdir -p  {config["home_dir"]}{config["target_dir"].split("/")[-2]}_output/{samples[i]}/benchmarks')

#RUN AQUAMIS
#rewrite as rule - hard to lock files to use in further processing
#os.chdir(aquamis)   
#subprocess.check_call(['qsub', '-F', f'{sample_list_path} {config["home_dir"]}{config["target_dir"]}_output/aquamis_output/', 'run_aquamis.sh']) #replace output dir

os.chdir(home_pipe)
os.system(f'mv  {config["target_dir"]}* data/')


#FINAL RULE
rule all:
    input: 
        target_list
    run:
        for i in range(len(samples)):
            os.system(f'python scripts/combine_benchmarks.py -b benchmarks -s {samples[i]}')
            os.system(f'mv {samples[i]}* {config["home_dir"]}{config["target_dir"].split("/")[-1]}_output/{config["home_dir"]}{samples[i]}/')
            os.system(f'mv data/{samples[i]}_fastp* {config["home_dir"]}{config["target_dir"].split("/")[-1]}_output/{config["home_dir"]}{samples[i]}/')
            os.system(f'mv benchmarks/{samples[i]}* {config["home_dir"]}{config["target_dir"].split("/")[-1]}_output/{config["home_dir"]}{samples[i]}/benchmarks/')
            os.system(f'mv {kfinder_output}* {config["home_dir"]}{config["target_dir"].split("/")[-1]}_output/kmerfinder_output/')

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
        temp('{sample_id_pattern}.fastp.html'),
        read_1_tr = 'data/{sample_id_pattern}_fastp_R1_001.fastq.gz',
        read_2_tr = 'data/{sample_id_pattern}_fastp_R2_001.fastq.gz'
    shell:
        'singularity run {input.sif_file} fastp -j {wildcards.sample_id_pattern}.fastp.json -h {wildcards.sample_id_pattern}.fastp.html --in1 {input.read_1} --in2 {input.read_2} --out1 {output.read_1_tr} --out2 {output.read_2_tr} --thread {threads}'

#RUN KRAKEN2 TO FILTER OUT HOST CONTAMINATION
rule filter_host:
    input:
        #quality-trimmed reads
        read_1 = 'data/{sample_id_pattern}_fastp_R1_001.fastq.gz',
        read_2 = 'data/{sample_id_pattern}_fastp_R2_001.fastq.gz'
    output:
        #host reads (temp)
        temp('data/{sample_id_pattern}_host_1.fastq'),
        temp('data/{sample_id_pattern}_host_2.fastq'),
        #sample reads 
        sample_1 = 'data/{sample_id_pattern}_sample_1.fastq.gz',
        sample_2 = 'data/{sample_id_pattern}_sample_2.fastq.gz',
        report_name = '{sample_id_pattern}_contigs/{sample_id_pattern}_kraken2_report.txt'
    threads: 48
    conda:
        kraken2_env
    shell:
        """ 
        kraken2 --threads 48 --db /mnt/home/groups/nmrl/db/db-kraken2/human_reference/ --classified-out data/{wildcards.sample_id_pattern}_host#.fastq --unclassified-out data/{wildcards.sample_id_pattern}_sample#.fastq --report {output.report_name} --gzip-compressed --paired {input.read_1} {input.read_2}
        pigz data/{wildcards.sample_id_pattern}_sample_1.fastq
        pigz data/{wildcards.sample_id_pattern}_sample_2.fastq
        """

#GENERATING CONTIGS FROM READS
rule contig_assembly:
    input:
        sif_file = snakemake_sif,
        read_1 = 'data/{sample_id_pattern}_sample_1.fastq.gz',
        read_2 = 'data/{sample_id_pattern}_sample_2.fastq.gz'
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
        """ singularity run {input.sif_file} ragtag.py scaffold -o {wildcards.sample_id_pattern}-{wildcards.reference_sequence_pattern}-scaffolds -C {input.reference} {input.contigs}
        if [[ -f {sample_id_pattern}-{reference_sequence_pattern}-scaffolds/ragtag.scaffold.fasta ]] ; then echo 'Scaffolding succesful!' ; else touch {output} ; fi """

#RENAMING SCAFFOLDS
rule scaffold_id:
    input:
        scf = '{sample_id_pattern}-{reference_sequence_pattern}-scaffolds/ragtag.scaffold.fasta'
    output:
        '{sample_id_pattern}-{reference_sequence_pattern}-scaffolds/{sample_id_pattern}-{reference_sequence_pattern}-ragtag.scaffold.fasta'
    shell:
        'cp {input.scf} {output}'

#RUN QUAST
rule quast_scaffolds:
    input:
        sif_file = mlst_quast_sif,
        reference = 'reference/{reference_sequence_pattern}.fasta',
        scaffold = '{sample_id_pattern}-{reference_sequence_pattern}-scaffolds/{sample_id_pattern}-{reference_sequence_pattern}-ragtag.scaffold.fasta'
    envmodules:
        'singularity'
    output:
        '{sample_id_pattern}-{reference_sequence_pattern}-scaffolds/icarus.html'
    shell:
        'singularity run {input.sif_file} quast -r {input.reference} -o {wildcards.sample_id_pattern}-{wildcards.reference_sequence_pattern}-scaffolds {input.scaffold}'

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
        mlst_quast_path = mlst_quast_sif,
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


#PERFORM GENE ANNOTATION
rule prokka: 
    input: 
        prokka_path = prokka_sif,
        contigs = '{sample_id_pattern}_contigs/{sample_id_pattern}_contigs.fasta'
    output:
        prokka_outdir = directory('{sample_id_pattern}-prokka')
    threads: 4
    envmodules:
        'singularity'
    benchmark:
        temp('benchmarks/{sample_id_pattern}.prokka.benchmark.txt')
    shell:
        """singularity run {input.prokka_path} prokka --outdir {output.prokka_outdir} {input.contigs}"""
        

#PERFORM RESISTANCE GENE EXTRACTION
rule res_gen_id: 
    input: 
        contigs = '{sample_id_pattern}_contigs/{sample_id_pattern}_contigs.fasta'
    output:
        temp('{sample_id_pattern}.rgi.txt'),
        temp('{sample_id_pattern}.rgi.json')
    threads: 4
    conda:
        rgi_env
    benchmark:
        temp('benchmarks/{sample_id_pattern}.rig.benchmark.txt')
    shell:
        """rgi main --input_sequence {input.contigs} --output_file {wildcards.sample_id_pattern}.rgi --input_type contig --clean"""


# on a cluster - https://carpentries-incubator.github.io/workflows-snakemake/09-cluster/index.html
# Add KRAKEN2 HOST FILTERING (USE GH38 human assembly)
# Add quast to automate qc (for each reference + metaquast based on all references at once (using contigs))
# Fix benchmarking