#Date of analysis is starting Nov. 12, 2020
#The samples are from Genewiz, stranded library prep of my infection with EPEC for 9hrs of polarized TC7 cells
#The analysis is going to be as follows: confirm transfer (md5sum) -> QC the reads (fastp) -> align the reads (STAR) -> count the reads (feautreCounts) 
#-> analyze with EDAseq, DESeq2

#In the Stranded directory in Graham server confirm all the md5 sums
nano sumcheck.sh #create a script
#!/bin/bash
#SBATCH --account=def-bfinlay
#SBATCH --time=4:00:00 	
#SBATCH --cpus-per-task=16
#SBATCH --mem=32000
set -e 
set -u
set -o pipefail
for i in "*" #inititate a for loop for each accession number
do
	md5sum ${i} >> local_check.txt
done #exit nano
chmod +x sumcheck.sh
sbatch sumcheck.sh

#Trimming reads for quality
	nano fastp1.sh #create a script to trim the low-quality reads and adaptors from each run 
#!/bin/bash
#SBATCH --account=def-bfinlay
#SBATCH --time=4:00:00 	
#SBATCH --cpus-per-task=8
#SBATCH --mem=16000
set -e 
set -u
set -o pipefail
module load fastp
for i in "UI9" "N9" 
do 
for j in "2" "3" "6"
do
	fastp -i ${i}${j}_R1_001.fastq.gz -I ${i}${j}_R2_001.fastq.gz -o ${i}${j}_1_out.fastq -O ${i}${j}_2_out.fastq -h ${i}${j}.html -3 -M 20 --trim_poly_x & #use fastp program to trim and filter forward and reverse reads, output files into
	#"accession number"_1_out.fq or "accession number"_2_out.fq and create an html report named after the accession number. Trim low-qulaity 3' basepairs, trim basepairs with 
	#quality lower than 28 and trim polyA tails, if they are present	
done
done
for i in "WT9"
do
for j in "1" "3" "5"
do
	fastp -i ${i}${j}_R1_001.fastq.gz -I ${i}${j}_R2_001.fastq.gz -o ${i}${j}_1_out.fastq -O ${i}${j}_2_out.fastq -h ${i}${j}.html -3 -M 20 --trim_poly_x #use fastp program to trim and filter forward and reverse reads, output files into
	#"accession number"_1_out.fq or "accession number"_2_out.fq and create an html report named after the accession number. Trim low-qulaity 3' basepairs, trim basepairs with 
	#quality lower than 28 and trim polyA tails, if they are present	
done
done
mv *out.fastq ./Trimmed_reads #move all filtered reads into Trimmed_reads directory
mv *.html ./html_qc_files
 #exit nano
chmod +x fastp1.sh #give permissions to run
sbatch fastp1.sh #submit a job
#Prepare STAR to map the reads onto the human genome
module load StdEnv/2020 nixpkgs/16.09 gcc/7.3.0
module load star/2.7.5a
mkdir gencode38
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.primary_assembly.annotation.gtf.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/GRCh38.primary_assembly.genome.fa.gz
gunzip gencode.v35.primary_assembly.annotation.gtf.gz
gunzip GRCh38.primary_assembly.genome.fa.gz

nano idx_gen.sh #create a script to generate genome index
#!/bin/bash
#SBATCH --account=def-bfinlay
#SBATCH --time=6:00:00 	
#SBATCH --cpus-per-task=16
#SBATCH --mem=64000
module load StdEnv/2020 nixpkgs/16.09 gcc/7.3.0
module load star/2.7.5a
STAR --runThreadN 16 --runMode genomeGenerate \
--genomeDir ./STAR_index --genomeFastaFiles ./GRCh38.primary_assembly.genome.fa \
--sjdbGTFfile ./gencode.v35.primary_assembly.annotation.gtf --sjdbOverhang 149
#exit nano

#Next step is to map the reads onto this annotation
nano test_align.sh #create a test script to align the reads to the genome
#!/bin/bash
#SBATCH --account=def-bfinlay
#SBATCH --time=6:00:00 	
#SBATCH --cpus-per-task=16
#SBATCH --mem=64000
module load StdEnv/2020 nixpkgs/16.09 gcc/7.3.0
module load star/2.7.5a
STAR --runThreadN 16 --genomeDir ./gencode38/STAR_index --readFilesIn UI92_1_out.fastq UI92_2_out.fastq --outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx --quantMode GeneCounts #exit nano
#Okay this all went well, now create three scripts to create the 
nano align_UI.sh #create a test script to align the reads to the genome
	#!/bin/bash
	#SBATCH --account=def-bfinlay
	#SBATCH --time=3:00:00 	
	#SBATCH --cpus-per-task=16
	#SBATCH --mem=64000
	module load StdEnv/2020 nixpkgs/16.09 gcc/7.3.0
	module load star/2.7.5a
	for i in "2" "3" "6"
	do
	mkdir ./UI9${i}_STAR_out
	STAR --runThreadN 16 --genomeDir ./gencode38/STAR_index --readFilesIn UI9${i}_1_out.fastq UI9${i}_2_out.fastq --outSAMtype BAM Unsorted SortedByCoordinate \
	--outReadsUnmapped Fastx --quantMode GeneCounts --outFileNamePrefix ./UI9${i}_STAR_out
	done
#exit nano
nano allign_N.sh #create a script to align dEscN reads to the genome
	#!/bin/bash
	#SBATCH --account=def-bfinlay
	#SBATCH --time=3:00:00 	
	#SBATCH --cpus-per-task=16
	#SBATCH --mem=64000
	module load StdEnv/2020 nixpkgs/16.09 gcc/7.3.0
	module load star/2.7.5a
	for i in "2" "3" "6"
	do
	mkdir ./N9${i}_STAR_out
	STAR --runThreadN 16 --genomeDir ./gencode38/STAR_index --readFilesIn N9${i}_1_out.fastq N9${i}_2_out.fastq --outSAMtype BAM Unsorted SortedByCoordinate \
	--outReadsUnmapped Fastx --quantMode GeneCounts --outFileNamePrefix ./N9${i}_STAR_out
	done
nano align_WT.sh #now aling WT reads
	#!/bin/bash
	#SBATCH --account=def-bfinlay
	#SBATCH --time=3:00:00 	
	#SBATCH --cpus-per-task=16
	#SBATCH --mem=64000
	module load StdEnv/2020 nixpkgs/16.09 gcc/7.3.0
	module load star/2.7.5a
	for i in "1" "3" "5"
	do
	mkdir ./WT9${i}_STAR_out
	STAR --runThreadN 16 --genomeDir ./gencode38/STAR_index --readFilesIn WT9${i}_1_out.fastq WT9${i}_2_out.fastq --outSAMtype BAM Unsorted SortedByCoordinate \
	--outReadsUnmapped Fastx --quantMode GeneCounts --outFileNamePrefix ./WT9${i}_STAR_out
	done
#The reads are counted as htseq count as part of the STAR command. Can just export the STAR_gene_counts for further analysis in R 
#The last step before exporting to R is to combine all files into one count table
paste UI92*Reads* UI93*Reads* UI96*Reads* N92*Reads* N93*Reads* N96*Reads* WT91*Reads* WT93*Reads* WT95*Reads* > all_reads_STAR


#Try and align/count EPEC genes

##########################################################
#This step is done in April 28, 2022
#Map and count unmapped reads to EPEC genome
#Generate genome index for EPEC using ENSEMBL genome
#Download EPEC genome and gtf annotation
#Use index built for short ifnection
#Do this in the Trimmed reads directory in Polar RNA seq Stranded folde
ln -s /home/zakhar/Trimmed_reads/STAR_EPEC_index/
#align and count EPEC reads
nano EPEC_align_9.sh
#!/bin/bash
#SBATCH --account=def-bfinlay
#SBATCH --time=6:00:00 	
#SBATCH --cpus-per-task=16
#SBATCH --mem=64000
module load StdEnv/2020 star/2.7.5a
for i in "UI92" "UI93" "UI96" "N92" "N93" "N96" "WT91" "WT93" "WT95"
do
STAR --runThreadN 16 --genomeDir ./STAR_EPEC_index --readFilesIn ${i}_STAR_outUnmapped.out.mate1 ${i}_STAR_outUnmapped.out.mate2 --outSAMtype BAM Unsorted \
--alignIntronMax 1 --quantMode GeneCounts --outFileNamePrefix ./EPEC_polar_mapping/${i}_STAR_out
done
#exit nano
sbatch EPEC_align_9.sh

#The reads are counted as htseq count as part of the STAR command. Can just export the STAR_gene_counts for further analysis in R 
#The last step before exporting to R is to combine all files into one count table
paste UI92*Reads* UI93*Reads* UI96*Reads* N92*Reads* N93*Reads* N96*Reads* WT91*Reads* WT93*Reads* WT95*Reads* > all_EPEC_reads_STAR.txt




















#Now we can repeat the analysis for the miRNA library using miRDeep
#First check that all the files transferred properly
sbatch ../Stranded/sumcheck.sh
#All samples were fine, now try running fastp analysis to see how the data looks
sbatch ../Stranded/fastp1.sh
#Probably from here we can move onto mirDeep analysis
#In the Trimmed_reads directory:
#Clone and install git directory for mirDeep2
git clone https://github.com/rajewsky-lab/mirdeep2.git
#In the mirDeep2 directory clone the patch for mirdeep2
git clone https://github.com/Drmirdeep/mirdeep2_patch.git
######################################################## MIRDEEP2 TIME########################## version 2.0.1.2
#########################################
##############################
#Install mirdeep2
module load StdEnv/2020
module load perl/5.30.2
perl install.pl #make sure it says install successful
#patch the mirdeep2 as well
cd mirdeep2_patch
bash patchme.sh
#Now we wait for the 3' adapter sequence from genewiz to perform analysis
#Get the human mature and hairpin miRNA sequences
wget ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.gz
perl mirdeep2/src/extract_miRNAs.pl hairpin.fa hsa > hairpin_ref.fa
wget ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz
perl mirdeep2/src/extract_miRNAs.pl mature.fa hsa > mature_ref.fa
#Get the human genome and build an index
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz
nano index_build.sh
#!/bin/bash
#SBATCH --account=def-bfinlay
#SBATCH --time=3:00:00 	
#SBATCH --cpus-per-task=16
#SBATCH --mem=64000
module load StdEnv/2020
module load perl/5.30.2
perl mirdeep2/bin/bowtie-build GRCh38.primary_assembly.genome.fa Gencode38_human
#Trying mapping, and quantifying one of the files
nano map-test.sh
#!/bin/bash
#SBATCH --account=def-bfinlay
#SBATCH --time=3:00:00 	
#SBATCH --cpus-per-task=16
#SBATCH --mem=64000
module load StdEnv/2020
module load perl/5.30.2
perl mirdeep2/bin/mapper.pl  N92_1_out.fastq -e -j -h -k TGGAATTCTCGGGTGCCAAGGC  -l 18 -m -p Gencode38_human -s N92_collapsed.fa -t N92_collapsed_vs_genome.arf -v
perl mirdeep2/bin/quantifier.pl -p hairpin_ref.fa -m mature_ref.fa -r N92_collapsed.fa -t hsa -y 16_19

nano quant-test.sh
#!/bin/bash
#SBATCH --account=def-bfinlay
#SBATCH --time=3:00:00 	
#SBATCH --cpus-per-task=16
#SBATCH --mem=64000
module load StdEnv/2020 perl/5.30.2
module load perl/5.30.2
perl mirdeep2/bin/quantifier.pl -p hairpin_ref.fa -m mature_ref.fa -r N92_collapsed.fa -t hsa -y 16_19 -n -x

perl mirdeep2/bin/quantifier.pl -p hairpin_ref.fa -m mature_ref.fa -r N92_collapsed.fa -t hsa -y 16_19 -n -x 
#Try a workaround to get that dumb PDF/API2 into the server
wget https://fastapi.metacpan.org/source/SSIMMS/PDF-API2-2.038/lib/PDF/API2.pm  #Work on that later
#Need to merge the reads with BBmerge for proper processing
nano merge-test.sh
#!/bin/bash
#SBATCH --account=def-bfinlay
#SBATCH --time=3:00:00 	
#SBATCH --cpus-per-task=16
#SBATCH --mem=32000
module load StdEnv/2020 bbmap/38.86
bbmerge.sh in1=N92_1_out.fastq in2=N92_2_out.fastq out=N92_merged.fastq ihist=N92_merging.txt nzo=t showhiststats=t mininsert=12 minoverlap=12 
nano merge_all.sh
#!/bin/bash
#SBATCH --account=def-bfinlay
#SBATCH --time=3:00:00 	
#SBATCH --cpus-per-task=16
#SBATCH --mem=32000
module load StdEnv/2020 bbmap/38.86
for i in "UI9" "N9" 
do 
for j in "2" "3" "6"
do
	bbmerge.sh in1=${i}${j}_1_out.fastq in2=${i}${j}_2_out.fastq out=${i}${j}_merged.fastq ihist=${i}${j}_merging.txt nzo=t showhiststats=t mininsert=12 minoverlap=12 	
done
done
for i in "WT9"
do
for j in "1" "3" "5"
do
	bbmerge.sh in1=${i}${j}_1_out.fastq in2=${i}${j}_2_out.fastq out=${i}${j}_merged.fastq ihist=${i}${j}_merging.txt nzo=t showhiststats=t mininsert=12 minoverlap=12 	
done
done
#try installing the PDF::API2 module manually
cpan
install PDF::API2
nano map_and_quant_all.sh
#!/bin/bash
#SBATCH --account=def-bfinlay
#SBATCH --time=6:00:00 	
#SBATCH --cpus-per-task=16
#SBATCH --mem=64000
set -e 
set -u
set -o pipefail
module load StdEnv/2020 perl/5.30.2
module load perl/5.30.2
for i in "UI9" "N9" 
do 
for j in "2" "3" "6"
do
	perl mirdeep2/bin/mapper.pl  ${i}${j}_merged.fastq -e -j -h -k TGGAATTCTCGGGTGCCAAGGC  -l 18 -m -p Gencode38_human -s ${i}${j}_collapsed.fa -t ${i}${j}_collapsed_vs_genome.arf -v
	perl mirdeep2/bin/quantifier.pl -p hairpin_ref.fa -m mature_ref.fa -r ${i}${j}_collapsed.fa -t hsa -y ${i}${j}
done
done
for i in "WT9"
do
for j in "1" "3" "5"
do
	perl mirdeep2/bin/mapper.pl  ${i}${j}_merged.fastq -e -j -h -k TGGAATTCTCGGGTGCCAAGGC  -l 18 -m -p Gencode38_human -s ${i}${j}_collapsed.fa -t ${i}${j}_collapsed_vs_genome.arf -v
	perl mirdeep2/bin/quantifier.pl -p hairpin_ref.fa -m mature_ref.fa -r ${i}${j}_collapsed.fa -t hsa -y ${i}${j} 
done
done

nano map_and_quant_after_N92.sh
#!/bin/bash
#SBATCH --account=def-bfinlay
#SBATCH --time=6:00:00 	
#SBATCH --cpus-per-task=16
#SBATCH --mem=64000
set -e 
set -u
set -o pipefail
module load StdEnv/2020 perl/5.30.2
module load perl/5.30.2
for i in "N9" 
do 
for j in "2" "3" "6"
do
	perl mirdeep2/bin/mapper.pl  ${i}${j}_merged.fastq -e -j -h -k TGGAATTCTCGGGTGCCAAGGC  -l 18 -m -p Gencode38_human -s ${i}${j}_collapsed.fa -t ${i}${j}_collapsed_vs_genome.arf -v
	perl mirdeep2/bin/quantifier.pl -p hairpin_ref.fa -m mature_ref.fa -r ${i}${j}_collapsed.fa -t hsa -y ${i}${j}
done
done
for i in "WT9"
do
for j in "1" "3" "5"
do
	perl mirdeep2/bin/mapper.pl  ${i}${j}_merged.fastq -e -j -h -k TGGAATTCTCGGGTGCCAAGGC  -l 18 -m -p Gencode38_human -s ${i}${j}_collapsed.fa -t ${i}${j}_collapsed_vs_genome.arf -v
	perl mirdeep2/bin/quantifier.pl -p hairpin_ref.fa -m mature_ref.fa -r ${i}${j}_collapsed.fa -t hsa -y ${i}${j} 
done
done #sbatch 41718862
#The last step is to merge the expression files into one big table
#The last step before exporting to R is to combine all files into one count table
paste miRNAs_*UI92* miRNAs_*UI93* miRNAs_*UI96* miRNAs_*N92* miRNAs_*N93* miRNAs_*N96* miRNAs_*WT91* miRNAs_*WT93* miRNAs_*WT95* > all_reads_mirDEEP2


