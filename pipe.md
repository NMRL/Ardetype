## tools 
   - [snakemake](https://snakemake.readthedocs.io/en/stable/) - controlling flow of excecution
   - [bactpipe](https://bactpipe.readthedocs.io/en/latest/running.html) - read assembly into contigs
   - [ragtag](https://github.com/malonge/RagTag) - contig assembly into genome
   - [quast](https://github.com/ablab/quast) - quality control
## input: 
   - species of bacteria (to download reference genome)
   - two fastq files of 150bp reads (should contain sample id in format [0-9]{10})
## output:
   - fasta file with longest assembled scaffold produced by ragstag
   - quast report on contigs produced by bactpipe
   - quast report on scaffold produced by ragtag
   - [agp](https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/) file with information about the process of scaffold assembly
## workflow:
   - [download](https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/) reference genome
   - run bactpipe on fastq files > contig fasta file
   - run quast on contig fasta file > contig quality report
   - run ragtag on contig fasta file > scaffolds fasta file
   - run quast on scaffolds fasta file > scaffold quality report
## further options:
   - https://github.com/tseemann/mlst - database access
