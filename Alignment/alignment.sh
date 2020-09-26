#############################
#Title:STAR Alignment
#Written by: Conor Cremin


mkdir out
fastqDIR=$PWD  #Specify directory were fastq files are located.

for file in $(ls $fastqDIR/*.fastq | cut -c98- | rev | cut -c10- | rev | uniq); 
do
	/software/STAR/2.7.2a/bin/STAR --readFilesCommand gunzip -c \
				--outFilterScoreMinOverLread 0.3 \
				--outFilterMatchNminOverLread 0.3 \
				--outFileNamePrefix $out/{file} \
				--runThreadN 2 --genomeDir $PWD/ensembl_STAR_Hamester_COVID-19 \
				--outSAMtype BAM SortedByCoordinate --quantMode GeneCounts \
				--readFilesIn ${file}.1.fastq.gz ${file}.2.fastq.gz
done
