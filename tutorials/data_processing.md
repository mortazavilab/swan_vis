# Data processing and TALON

For the Swan manuscript, these are the data processing steps we took to process our data before using Swan. Obtaining a transcriptome via TALON is not necessary and you can use other tools that yield transcriptomes as input to Swan.

If you don't want to process the raw data from square one and want to get started using Swan please see the [Getting started](getting_started.md) section.

You can download and view documentation for TALON [here](https://github.com/mortazavilab/TALON).

Each of the below steps should be run in your Bash terminal

## Download reference genome and annotation
```bash
mkdir data
mkdir figures
cd data/

# download reference files
wget https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz
gunzip GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz
awk '{print $1}' GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta > hg38.fa

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz
gunzip gencode.v29.annotation.gtf.gz
```

## Download mapped PacBio RNA-seq data for each replicate
```bash
# hepg2 rep 1
wget https://www.encodeproject.org/files/ENCFF575ABA/@@download/ENCFF575ABA.bam
mv ENCFF575ABA.bam hepg2_1.bam

# hepg2 rep 2
wget wget https://www.encodeproject.org/files/ENCFF463JKW/@@download/ENCFF463JKW.bam
mv ENCFF463JKW.bam hepg2_2.bam

# hffc6 rep 1
wget https://www.encodeproject.org/files/ENCFF677TYB/@@download/ENCFF677TYB.bam
mv ENCFF677TYB.bam hffc6_1.bam 

# hffc6 rep 2
wget https://www.encodeproject.org/files/ENCFF584UQW/@@download/ENCFF584UQW.bam
mv ENCFF584UQW.bam hffc6_2.bam

# hffc6 rep 3
wget https://www.encodeproject.org/files/ENCFF964PGS/@@download/ENCFF964PGS.bam
mv ENCFF964PGS.bam hffc6_3.bam
```

## Convert each bam to sam 
```bash
module load samtools
for file in *.bam;
do
	sam_file=`basename $file .bam`
	sam_file=${sam_file}.sam
	samtools view -h $file > $sam_file
done
```

## Label reads with TALON for detection of internal priming
```bash
FASTA=hg38.fa
mkdir labelled
for file in *.sam;
do
	base=`basename $file .sam`
	talon_label_reads \
		--f $file \
		--g $FASTA \
		--t 8 \
		--tmpDir ${base}_temp \
		--deleteTmp \
		--o labelled/${base}
done
```

## Initialize TALON database with GENCODE v29 annotation
```bash 
cd ../

GTF=data/gencode.v29.annotation.gtf
talon_initialize_database \
	--f $GTF \
	--g hg38 \
	--a gencode_v29 \
	--l 0 \
	--idprefix TALON \
	--5p 500 \
	--3p 300 \
	--o talon
```

## Create config file to run TALON 
```bash
plat="Sequel"
touch config.csv
for file in data/labelled/*.sam;
do
	base=`basename $file .sam`
	base="${base::${#base}-8}"
	sample="${base::${#base}-2}"
	printf "${base},${sample},${plat},${file}\n" >> config.csv
done
```

## Run TALON
```bash
talon \
	--f config.csv \
	--db talon.db \
	--t 16 \
	--build hg38 \
	--cov 0.9 \
	--identity 0.8 \
	--o talon_swan
```

## Filter annotated transcripts for internal priming and reproducibility for each cell line
```bash
talon_filter_transcripts \
	--db talon.db \
	-a gencode_v29 \
	--datasets hepg2_1,hepg2_2 \
	--maxFracA 0.5 \
	--minCount 5 \
	--minDatasets 2 \
	--o hepg2_pass_list.csv

talon_filter_transcripts \
	--db talon.db \
	-a gencode_v29 \
	--datasets hffc6_1,hffc6_2,hffc6_3 \
	--maxFracA 0.5 \
	--minCount 5 \
	--minDatasets 2 \
	--o hffc6_pass_list.csv
    
cat hepg2_pass_list.csv > temp.csv
cat hffc6_pass_list.csv >> temp.csv 
cat temp.csv | uniq > all_pass_list.csv
rm temp.csv
```

## Create a GTF of annotated transcripts for all observed transcripts that passed filtering
```bash
talon_create_GTF \
	--db talon.db \
	-b hg38 \
	-a gencode_v29 \
	--whitelist all_pass_list.csv \
	--observed \
	--o all
```

## Obtain a filtered abundance file of each transcript that passed filtering
```bash
talon_abundance \
	--db talon.db \
	-a gencode_v29 \
	-b hg38 \
	--whitelist all_pass_list.csv \
	--o all
```