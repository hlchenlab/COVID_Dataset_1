#############################
#Title:STAR Index Build
#Written by: Conor Cremin


mkdir STAR_Hamester_COVID-19
STAR --runThreadN 11 \
     --runMode genomeGenerate \
     --genomeDir STAR_Hamester_COVID-19 \
     --genomeFastaFiles STAR_Hamester_COVID-19.fa \
     --sjdbGTFfile STAR_Hamester_COVID-19.gtf \
     --sjdbOverhang 150 
