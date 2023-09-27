#Date of analysis: Nov. 12, 2021
#Doing the analysis on the Graham server, starting from Final RNA Seq folder, so after fastp processing


#Trimming reads for quality, had to restart with trimming, since there was some weird artifact with gzipped reads
#This is exectued in the Final RNA seq folder
	nano fastp2.sh #create a script to trim the low-quality reads and adaptors from each run 
	#!/bin/bash
	#SBATCH --time=10:00:00 	
	#SBATCH --cpus-per-task=8
	#SBATCH --mem=16000
	#SBATCH --account=def-bfinlay
	set -e 
	set -u
	set -o pipefail
	for i in "4UIN" "4UIWT" "490WT" "490N" "4180N" "4180WT" "5UIN" "5UIWT" "545WT" "545N" "590WT" "590N" "5180N" "5180WT" "6UIN" "6UIWT" "645WT" "645N" "690WT" "690N" "6180N" "6180WT" #inititate a for loop for each accession number
	do
		./fastp -i ${i}*L001_R1_001.fastq.gz -I ${i}*L001_R2_001.fastq.gz -o ${i}_L001_1_out.fastq -O ${i}_L001_2_out.fastq -h ${i}_L001.html -3 -M 28 --trim_poly_x #use fastp program to trim and filter forward and reverse reads, output files into
		#"accession number"_1_out.fq or "accession number"_2_out.fq and create an html report named after the accession number. Trim low-qulaity 3' basepairs, trim basepairs with 
		#quality lower than 28 and trim polyA tails, if they are present
		./fastp -i ${i}*L002_R1_001.fastq.gz -I ${i}*L002_R2_001.fastq.gz -o ${i}_L002_1_out.fastq -O ${i}_L002_2_out.fastq -h ${i}_L002.html -3 -M 28 --trim_poly_x
		./fastp -i ${i}*L003_R1_001.fastq.gz -I ${i}*L003_R2_001.fastq.gz -o ${i}_L003_1_out.fastq -O ${i}_L003_2_out.fastq -h ${i}_L003.html -3 -M 28 --trim_poly_x
		./fastp -i ${i}*L004_R1_001.fastq.gz -I ${i}*L004_R2_001.fastq.gz -o ${i}_L004_1_out.fastq -O ${i}_L004_2_out.fastq -h ${i}_L004.html -3 -M 28 --trim_poly_x
		cat ${i}*1_out.fastq > merged_${i}_1_out.fastq #merge the fastq files
		cat ${i}*2_out.fastq > merged_${i}_2_out.fastq
	done
	mv *out.fastq /home/zakhar/projects/def-bfinlay/zakhar/Final_RNA_seq/Trimmed_reads #move all filtered reads into Trimmed_reads directory
	 #exit nano
	chmod +x fastp2.sh #give permissions to run
	sbatch fastp2.sh #submit a job

# link out towards the human genome index
ln -s /project/6003396/zakhar/Polar_RNA_seq/Stranded/Trimmed_reads/gencode38



#Run a test align script from the Trimmed Reads directory
nano test_align.sh
#!/bin/bash
#SBATCH --account=def-bfinlay
#SBATCH --time=6:00:00 	
#SBATCH --cpus-per-task=16
#SBATCH --mem=64000
module load StdEnv/2020 nixpkgs/16.09 gcc/7.3.0
module load star/2.7.5a
STAR --runThreadN 16 --genomeDir ./gencode38/STAR_index --readFilesIn merged_590WT_1_out.fastq merged_590WT_2_out.fastq --outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx --quantMode GeneCounts #exit nano
#Ran well, now can submit the rest of mappingg

nano 4_align.sh
	#!/bin/bash
	#SBATCH --account=def-bfinlay
	#SBATCH --time=6:00:00 	
	#SBATCH --cpus-per-task=16
	#SBATCH --mem=64000
	module load StdEnv/2020 nixpkgs/16.09 gcc/7.3.0
	module load star/2.7.5a
	for i in "4UIN" "4UIWT" "490WT" "490N" "4180N" "4180WT"
	do
	mkdir ./UI9${i}_STAR_out
	STAR --runThreadN 16 --genomeDir ./gencode38/STAR_index --readFilesIn merged_${i}_1_out.fastq merged_${i}_2_out.fastq --outSAMtype BAM Unsorted SortedByCoordinate \
	--outReadsUnmapped Fastx --quantMode GeneCounts --outFileNamePrefix ./${i}_STAR_out
	done #need redo, because of lack of size. Failed at 4180WT
nano 5_align.sh
	#!/bin/bash
	#SBATCH --account=def-bfinlay
	#SBATCH --time=6:00:00 	
	#SBATCH --cpus-per-task=16
	#SBATCH --mem=64000
	module load StdEnv/2020 nixpkgs/16.09 gcc/7.3.0
	module load star/2.7.5a
	for i in "5UIN" "5UIWT" "545WT" "545N" "590WT" "590N" "5180N" "5180WT"
	do
	mkdir ./UI9${i}_STAR_out
	STAR --runThreadN 16 --genomeDir ./gencode38/STAR_index --readFilesIn merged_${i}_1_out.fastq merged_${i}_2_out.fastq --outSAMtype BAM Unsorted SortedByCoordinate \
	--outReadsUnmapped Fastx --quantMode GeneCounts --outFileNamePrefix ./${i}_STAR_out
	done #need redo, failed at 5180WT

nano 6_align.sh
#!/bin/bash
#SBATCH --account=def-bfinlay
#SBATCH --time=6:00:00 	
#SBATCH --cpus-per-task=16
#SBATCH --mem=64000
module load StdEnv/2020 nixpkgs/16.09 gcc/7.3.0
module load star/2.7.5a
for i in "6UIN" "6UIWT" "645WT" "645N" "690WT" "690N" "6180N" "6180WT"
do
STAR --runThreadN 16 --genomeDir ./gencode38/STAR_index --readFilesIn merged_${i}_1_out.fastq merged_${i}_2_out.fastq --outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx --quantMode GeneCounts --outFileNamePrefix ./${i}_STAR_out
done #worked well

nano leftover_align.sh
#!/bin/bash
#SBATCH --account=def-bfinlay
#SBATCH --time=6:00:00 	
#SBATCH --cpus-per-task=16
#SBATCH --mem=64000
module load StdEnv/2020 nixpkgs/16.09 gcc/7.3.0
module load star/2.7.5a
for i in "4180WT" "5180WT"
do
STAR --runThreadN 16 --genomeDir ./gencode38/STAR_index --readFilesIn merged_${i}_1_out.fastq merged_${i}_2_out.fastq --outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx --quantMode GeneCounts --outFileNamePrefix ./${i}_STAR_out
done 


#It seems like they got like ~30-70% of reads mapped nicely, likely due to rRNA contamination (high number of multi-mappers), anyway, combine the gene counts to a nice table and take it out
paste *4UIN*ReadsPerGene* *4UIWT*ReadsPerGene* *490WT*ReadsPerGene* *490N*ReadsPerGene* *4180N*ReadsPerGene* *4180WT*ReadsPerGene* *5UIN*ReadsPerGene* *5UIWT*ReadsPerGene* *545WT*ReadsPerGene* *545N*ReadsPerGene* *590WT*ReadsPerGene* *590N*ReadsPerGene* *5180N*ReadsPerGene* *5180WT*ReadsPerGene* *6UIN*ReadsPerGene* *6UIWT*ReadsPerGene* *645WT*ReadsPerGene* *645N*ReadsPerGene* *690WT*ReadsPerGene* *690N*ReadsPerGene* *6180N*ReadsPerGene* *6180WT*ReadsPerGene*  > STAR_reads_all.txt
#continue this in R, the same way as for the polar RNA-seq analysis

##########################################################
#This step is done in April 28, 2022
#Map and count unmapped reads to EPEC genome
#Generate genome index for EPEC using ENSEMBL genome
#Download EPEC genome and gtf annotation
wget -O EPEC_genome.rel-41.fa.gz http://ftp.ensemblgenomes.org/pub/bacteria/release-41/fasta/bacteria_10_collection/escherichia_coli_o127_h6_str_e2348_69/dna/Escherichia_coli_o127_h6_str_e2348_69.ASM2654v1.dna.toplevel.fa.gz
wget -O EPEC_annot.rel-41.gff3.gz http://ftp.ensemblgenomes.org/pub/bacteria/release-41/gff3/bacteria_10_collection/escherichia_coli_o127_h6_str_e2348_69/Escherichia_coli_o127_h6_str_e2348_69.ASM2654v1.37.gff3.gz
gunzip EPEC_genome.rel-41.fa.gz
gunzip EPEC_annot.rel-41.gff3.gz
#Convert the gff3 format into gtf format
module load StdEnv/2020 gcc/9.3.0 cufflinks/2.2.1
gffread EPEC_annot.rel-41.gff3 -T -o EPEC_annot.rel-41.gtf


mkdir STAR_EPEC_index
nano EPEC_idx_build.sh
#!/bin/bash
#SBATCH --account=def-bfinlay
#SBATCH --time=3:00:00 	
#SBATCH --cpus-per-task=16
#SBATCH --mem=32000
module load StdEnv/2020 nixpkgs/16.09 gcc/7.3.0
module load star/2.7.5a
STAR --runThreadN 16 --runMode genomeGenerate \
--genomeDir ./STAR_EPEC_index --genomeFastaFiles ./EPEC_genome.rel-41.fa \
--sjdbGTFfile ./EPEC_annot.rel-41.gtf --sjdbOverhang 149 --genomeSAindexNbases 8
#exit nano
sbatch EPEC_idx_build.sh

#align and count EPEC reads
nano EPEC_align_4_5.sh
#!/bin/bash
#SBATCH --account=def-bfinlay
#SBATCH --time=6:00:00 	
#SBATCH --cpus-per-task=16
#SBATCH --mem=64000
module load StdEnv/2020 star/2.7.5a
for i in "4UIN" "4UIWT" "490WT" "490N" "4180N" "4180WT" "5UIN" "5UIWT" "545WT" "545N" "590WT" "590N" "5180N" "5180WT"
do
STAR --runThreadN 16 --genomeDir ./STAR_EPEC_index --readFilesIn ${i}_STAR_outUnmapped.out.mate1 ${i}_STAR_outUnmapped.out.mate2 --outSAMtype BAM Unsorted \
--alignIntronMax 1 --quantMode GeneCounts --outFileNamePrefix ./EPEC_mapping/${i}_STAR_out
done
#exit nano
sbatch EPEC_align_4_5.sh


nano EPEC_align_6.sh
#!/bin/bash
#SBATCH --account=def-bfinlay
#SBATCH --time=6:00:00 	
#SBATCH --cpus-per-task=16
#SBATCH --mem=64000
module load StdEnv/2020 star/2.7.5a
for i in "6UIN" "6UIWT" "645WT" "645N" "690WT" "690N" "6180N" "6180WT"
do
STAR --runThreadN 16 --genomeDir ./STAR_EPEC_index --readFilesIn ${i}_STAR_outUnmapped.out.mate1 ${i}_STAR_outUnmapped.out.mate2 --outSAMtype BAM Unsorted \
--alignIntronMax 1 --quantMode GeneCounts --outFileNamePrefix ./EPEC_mapping/${i}_STAR_out
done
#exit nano

#It seems like they got like ~30-70% of reads mapped nicely, likely due to rRNA contamination (high number of multi-mappers), anyway, combine the gene counts to a nice table and tak it out
paste *4UIN*ReadsPerGene* *4UIWT*ReadsPerGene* *490WT*ReadsPerGene* *490N*ReadsPerGene* *4180N*ReadsPerGene* *4180WT*ReadsPerGene* *5UIN*ReadsPerGene* *5UIWT*ReadsPerGene* *545WT*ReadsPerGene* *545N*ReadsPerGene* *590WT*ReadsPerGene* *590N*ReadsPerGene* *5180N*ReadsPerGene* *5180WT*ReadsPerGene* *6UIN*ReadsPerGene* *6UIWT*ReadsPerGene* *645WT*ReadsPerGene* *645N*ReadsPerGene* *690WT*ReadsPerGene* *690N*ReadsPerGene* *6180N*ReadsPerGene* *6180WT*ReadsPerGene*  > STAR_EPEC_reads_all.txt
#continue this in R, the same way as for the polar RNA-seq analysis

######################################
######################################
#Date of analysis: Nov. 27, 2022
#Found the 445 samples - so should add them to analysis
#Doing the analysis on the Graham server, starting from Final RNA Seq folder, so after fastp processing
#Trimming reads for quality, had to restart with trimming, since there was some weird artifact with gzipped reads
#This is exectued in the Final RNA seq folder
screen -S 445
nano fastp445.sh #create a script to trim the low-quality reads and adaptors from each run 
	#!/bin/bash
	#SBATCH --time=1:00:00 	
	#SBATCH --cpus-per-task=8
	#SBATCH --mem=16000
	#SBATCH --account=def-bfinlay
	set -e 
	set -u
	set -o pipefail
	for i in "445WT" "445N" #inititate a for loop for each accession number
	do
		./fastp -i ${i}*L001_R1_001.fastq.gz -I ${i}*L001_R2_001.fastq.gz -o ${i}_L001_1_out.fastq -O ${i}_L001_2_out.fastq -h ${i}_L001.html -3 -M 28 --trim_poly_x #use fastp program to trim and filter forward and reverse reads, output files into
		#"accession number"_1_out.fq or "accession number"_2_out.fq and create an html report named after the accession number. Trim low-qulaity 3' basepairs, trim basepairs with 
		#quality lower than 28 and trim polyA tails, if they are present
		./fastp -i ${i}*L002_R1_001.fastq.gz -I ${i}*L002_R2_001.fastq.gz -o ${i}_L002_1_out.fastq -O ${i}_L002_2_out.fastq -h ${i}_L002.html -3 -M 28 --trim_poly_x
		./fastp -i ${i}*L003_R1_001.fastq.gz -I ${i}*L003_R2_001.fastq.gz -o ${i}_L003_1_out.fastq -O ${i}_L003_2_out.fastq -h ${i}_L003.html -3 -M 28 --trim_poly_x
		./fastp -i ${i}*L004_R1_001.fastq.gz -I ${i}*L004_R2_001.fastq.gz -o ${i}_L004_1_out.fastq -O ${i}_L004_2_out.fastq -h ${i}_L004.html -3 -M 28 --trim_poly_x
		cat ${i}*1_out.fastq > merged_${i}_1_out.fastq #merge the fastq files
		cat ${i}*2_out.fastq > merged_${i}_2_out.fastq
	done
	mv *out.fastq /home/zakhar/projects/def-bfinlay/zakhar/Final_RNA_seq/Trimmed_reads #move all filtered reads into Trimmed_reads directory
#exit nano
chmod +x fastp445.sh #give permissions to run
sbatch fastp445.sh #submit a job
##Run align script from the Trimmed Reads directory
nano 445_align.sh
	#!/bin/bash
	#SBATCH --account=def-bfinlay
	#SBATCH --time=2:00:00 	
	#SBATCH --cpus-per-task=16
	#SBATCH --mem=64000
	module load StdEnv/2020 nixpkgs/16.09 gcc/7.3.0
	module load star/2.7.5a
	for i in "445WT" "445N"
	do
	STAR --runThreadN 16 --genomeDir ./gencode38/STAR_index --readFilesIn merged_${i}_1_out.fastq merged_${i}_2_out.fastq --outSAMtype BAM Unsorted SortedByCoordinate \
	--outReadsUnmapped Fastx --quantMode GeneCounts --outFileNamePrefix ./${i}_STAR_out
	done 
#exit nano	
chmod +x 445_align.sh
sbatch 445_align.sh
paste *445N*ReadsPerGene* *445WT*ReadsPerGene* > STAR_reads_445.txt
paste STAR_reads_all.txt STAR_reads_445.txt > STAR_reads_updated.txt
#continue this in R, the same way as for the polar RNA-seq analysis
#Do the analysis for EPEC counts
#align and count EPEC reads
nano EPEC_align_445.sh
#!/bin/bash
#SBATCH --account=def-bfinlay
#SBATCH --time=6:00:00 	
#SBATCH --cpus-per-task=16
#SBATCH --mem=64000
module load StdEnv/2020 star/2.7.5a
for i in "445WT" "445N"
do
STAR --runThreadN 16 --genomeDir ./STAR_EPEC_index --readFilesIn ${i}_STAR_outUnmapped.out.mate1 ${i}_STAR_outUnmapped.out.mate2 --outSAMtype BAM Unsorted \
--alignIntronMax 1 --quantMode GeneCounts --outFileNamePrefix ./EPEC_mapping/${i}_STAR_out
done
#exit nano
chmod +x EPEC_align_445.sh
sbatch EPEC_align_445.sh
#Combine the EPEC reads in the EPEC_mapping folder
paste *445N*ReadsPerGene* *445WT*ReadsPerGene* > STAR_EPEC_reads_445.txt
paste STAR_EPEC_reads_all.txt STAR_EPEC_reads_445.txt > STAT_EPEC_reads_updated.txt