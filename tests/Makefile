all: ex1.glf ex1.pileup.gz ex1.bam.bai ex1.glfview.gz \
	ex2.sam.gz ex2.sam ex1.sam \
	ex2.bam \
	ex3.bam ex3.bam.bai \
	ex4.bam ex4.bam.bai \
	ex5.bam ex5.bam.bai \
	ex6.bam \
	ex8.bam

ex2.sam.gz: ex1.bam ex1.bam.bai
		samtools view -h ex1.bam | gzip > ex2.sam.gz

%.bam: %.sam ex1.fa.fai
	samtools import ex1.fa.fai $< $@

%.sam: %.sam.gz
	gunzip < $< > $@

ex1.fa.fai:ex1.fa
		samtools faidx ex1.fa
ex1.bam:ex1.sam.gz ex1.fa.fai
		samtools import ex1.fa.fai ex1.sam.gz ex1.bam
%.bam.bai:%.bam
		samtools index $<
ex1.pileup.gz:ex1.bam ex1.fa
		samtools pileup -cf ex1.fa ex1.bam | gzip > ex1.pileup.gz
ex1.glf:ex1.bam ex1.fa	
		samtools pileup -gf ex1.fa ex1.bam > ex1.glf
ex1.glfview.gz:ex1.glf
		samtools glfview ex1.glf | gzip > ex1.glfview.gz

clean:
		rm -fr *.bam *.bai *.glf* *.fai *.pileup* *~ calDepth *.dSYM pysam_*.sam ex2.sam ex2.sam.gz ex1.sam
