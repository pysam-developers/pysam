all: ex1.glf ex1.pileup.gz ex1.bam.bai ex1.glfview.gz ex2.sam.gz ex2.sam ex1.sam ex3.bam ex3.bam.bai ex4.bam ex4.bam.bai
		@echo; echo \# You can now launch the viewer with: \'samtools tview ex1.bam ex1.fa\'; echo;

ex2.sam.gz: ex1.bam ex1.bam.bai
		samtools view -h ex1.bam | gzip > ex2.sam.gz

ex3.bam: ex3.sam ex1.fa.fai
	samtools import ex1.fa.fai ex3.sam ex3.bam

ex4.bam: ex4.sam ex1.fa.fai
	samtools import ex1.fa.fai ex4.sam ex4.bam

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
