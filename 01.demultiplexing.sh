#Demultiplexing and primer trimming
#input files - primers.fasta and matrix.txt


ADAPT_PATH=$(echo /media/kreising/DATA/data/RUN_16_3_2019_areny/primers.fasta)
MATRIX_PATH=$(echo /media/kreising/DATA/data/RUN_16_3_2019_areny/matrix.txt)


cd /media/kreising/DATA/data/RUN_16_3_2019_areny/00.RAW_DATA
mkdir ../demultiplexed
rm ../demultiplexed/*

SAMPLES=$(ls -1| grep "R[12].fastq.gz"|sed -r 's/R[12].fastq.gz//'|sort| uniq)
for i in $SAMPLES
 do 
  skewer -x $ADAPT_PATH -M $MATRIX_PATH -b -m head -k 35 -d 0 -t 8 "$i"R1.fastq.gz "$i"R2.fastq.gz -o ../demultiplexed/"$i"
done >> ./log


cd ../demultiplexed/
gzip *fastq

cd ..
mkdir 02A.DEMULTI.16s
cd ./demultiplexed

SAMPLES_16S=$(ls -1| grep "assigned-F_"|sed -r 's/-pair[12].fastq.gz//'|sort| uniq)
for i in $SAMPLES_16S
 do 
   skewer -x NNNNNTACGGNNGGCWGCAG -y GACTACHVGGGTATCTAATCC  -m head -k 35 -d 0 -t 8 "$i"-pair1.fastq.gz "$i"-pair2.fastq.gz -o ../02A.DEMULTI.16s/"$i"
done >> ../02A.DEMULTI.16s/log.trimmed.head

cd ../02A.DEMULTI.16s/

gzip ./*fastq


