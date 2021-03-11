###Bash script to process TT-Seq data:

STAR --runThreadN 64 --genomeDir /path/to/index/STAR-hg19/ --readFilesIn MINT-Seq_R1.fastq.gz MINT-Seq_R2.fastq.gz --outFileNamePrefix MINT-Seq --outSAMtype BAM Unsorted --readFilesCommand zcat
makeTagDirectory MINT-Seq-homer MINT-SeqAligned.out.bam  -tbp 1 -sspe -flip
makeBigWig.pl MINT-Seq-homer hg19 -strand -o auto -webdir ./ -url ./

STAR --runThreadN 64 --genomeDir /path/to/index/STAR-ERCC --readFilesIn MINT-Seq_R1.fastq.gz MINT-Seq_R2.fastq.gz --outFileNamePrefix MINT-Seq-ERCC --outSAMtype BAM Unsorted --readFilesCommand zcat
samtools sort -@ 24 MINT-Seq-ERCCAligned.out.bam -o MINT-Seq-ERCC.sorted.bam
samtools index MINT-Seq-ERCC.sorted.bam
samtools idxstats MINT-Seq-ERCC.sorted.bam > MINT-Seq-ERCC.out

STAR --runThreadN 64 --genomeDir /path/to/index/STAR-m6Aspikein --readFilesIn MINT-Seq_R1.fastq.gz MINT-Seq_R2.fastq.gz --outFileNamePrefix MINT-Seq-m6Aspike --outSAMtype BAM Unsorted --readFilesCommand zcat
samtools sort -@ 24 MINT-Seq-m6AspikeAligned.out.bam -o MINT-Seq-m6Aspike.sorted.bam
samtools index MINT-Seq-m6Aspike.sorted.bam
samtools idxstats MINT-Seq-m6Aspike.sorted.bam > MINT-Seq-m6Aspike.out

###Bash script to prepocess MINT-Seq data
STAR --runThreadN 64 --genomeDir /path/to/index/STAR-hg19/ --readFilesIn TT-Seq_R1.fastq.gz TT-Seq_R2.fastq.gz --outFileNamePrefix TT-Seq --outSAMtype BAM Unsorted --readFilesCommand zcat
makeTagDirectory TT-Seq-homer TT-SeqAligned.out.bam  -tbp 1 -sspe -flip
makeBigWig.pl TT-Seq-homer hg19 -strand -o auto -webdir ./ -url ./

STAR --runThreadN 64 --genomeDir /path/to/index/STAR-ERCC --readFilesIn TT-Seq_R1.fastq.gz TT-Seq_R2.fastq.gz --outFileNamePrefix MINT-Seq-ERCC --outSAMtype BAM Unsorted --readFilesCommand zcat
samtools sort -@ 24 TT-Seq-ERCCAligned.out.bam -o TT-Seq-ERCC.sorted.bam
samtools index TT-Seq-ERCC.sorted.bam
samtools idxstats TT-Seq-ERCC.sorted.bam > TT-Seq-ERCC.out

STAR --runThreadN 64 --genomeDir /path/to/index/STAR-m6Aspikein --readFilesIn TT-Seq_R1.fastq.gz TT-Seq_R2.fastq.gz --outFileNamePrefix MINT-Seq-m6Aspike --outSAMtype BAM Unsorted --readFilesCommand zcat
samtools sort -@ 24 TT-Seq-m6AspikeAligned.out.bam -o TT-Seq-m6Aspike.sorted.bam
samtools index TT-Seq-m6Aspike.sorted.bam
samtools idxstats TT-Seq-m6Aspike.sorted.bam > TT-Seq-m6Aspike.out

###L1 quantification
analyzeRepeats.pl repeats hg19 -L3 L1 -fpkm -strand + -noCondensing -d TT-Seq-homer MINT-Seq-homer > MINT-TT.L1.rpkm.txt

###MACS2 peak calling
macs2 callpeak  -t MINT-SeqAligned.out.bam  -c TT-SeqAligned.out.bam  -f BAM -g hs -n MINT-Seq-macs2 -B -q 0.05