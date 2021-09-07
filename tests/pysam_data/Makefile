SAM=$(wildcard *.sam)
BAM=$(SAM:%.sam=%.bam)
BAI=$(BAM:%.bam=%.bam.bai)
CRAM=ex1.cram ex2.cram ex3.cram
CRAI=$(CRAM:%.cram=%.cram.crai)
NO_PG:=$(findstring --no-PG,$(shell samtools view '-?'))

# ex2.bam - bam file without index

all: all.stamp

all.stamp: ex1.pileup.gz \
	ex1.sam ex1.bam \
	ex2.sam.gz ex2.sam ex2.bam ex2.bam.bai \
	with_md.sam.gz with_md.bam with_md.bam.bai \
        uncompressed.bam \
	$(BAM) $(BAI) \
	$(CRAM) $(CRAI) \
	example_bai.bam \
        rg_with_tab.bam \
	ex2_truncated.bam \
	ex2_corrupt.bam \
	empty.bam empty.bam.bai \
	explicit_index.bam explicit_index.cram \
	faidx_empty_seq.fq.gz \
	ex1.fa.gz ex1.fa.gz.csi \
	ex1_csi.bam \
	example_reverse_complement.bam \
	example_dash_in_chr.bam
	touch $@

# ex2.sam - as ex1.sam, but with header
ex2.sam.gz: ex1.bam ex1.bam.bai
	samtools view $(NO_PG) -h ex1.bam | gzip > ex2.sam.gz

with_md.sam.gz: ex2.bam ex1.fa
	samtools calmd $(NO_PG) --output-fmt BAM $^ > $@

#%.bam: %.sam ex1.fa.fai
#	samtools view $(NO_PG) -bo $@ -t ex1.fa.fai $<

uncompressed.bam: ex2.sam
	samtools view $(NO_PG) -bu -o $@ $<

%.bam: %.sam
	samtools view $(NO_PG) -bo $@ $<

%.cram: %.sam
	samtools view $(NO_PG) -Co $@ -T ex1.fa $<

%.cram.crai: %.cram
	samtools index $<

%.sam: %.sam.gz
	gunzip < $< > $@

ex1.fa.fai:ex1.fa
		samtools faidx ex1.fa

ex1.bam:ex1.sam.gz ex1.fa.fai
		samtools view $(NO_PG) -bo ex1.bam -t ex1.fa.fai ex1.sam.gz

%.bam.bai:%.bam
		samtools index $<

ex1.pileup.gz:ex1.bam ex1.fa
		samtools mpileup -f ex1.fa ex1.bam | gzip > ex1.pileup.gz

ex2_truncated.bam: ex2.bam
	head -c 124000 ex2.bam > ex2_truncated.bam

# Append a corrupt read with block_size < sizeof(bam_core_t fields)
ex2_corrupt.bam: ex2.bam
	(bgzip -d < $<; printf '\37\0\0\0\1\0\0\0') | bgzip > $@

ex1_csi.bam: ex1.bam
	cp ex1.bam ex1_csi.bam
	samtools index -c ex1_csi.bam

empty.bam: ex2.sam
	grep "^@" $< | samtools view $(NO_PG) -bo $@ -

example_unmapped_reads_no_sq.bam: example_unmapped_reads_no_sq.sam
	touch tmp.list
	samtools view $(NO_PG) -bo $@ -t tmp.list $<
	rm -f tmp.list

example_bai.bam: ex1.bam
	cp ex1.bam $@
	samtools index $@
	mv $@.bai example_bai.bai

explicit_index.bam: ex1.bam
	cp ex1.bam $@

explicit_index.cram: ex1.cram
	cp ex1.cram $@

clean:
	rm -fr [a-z]*.bam *.bai *.csi *.fai *.gzi *.pileup* [a-z]*.cram *.crai \
	all.stamp *~ calDepth *.dSYM pysam_*.sam \
	ex2.sam ex2.sam.gz ex1.sam ex1.fa.gz \
	with_md.sam.gz \
	*.fq.gz

%.fq.gz: %.fq
	gzip < $< > $@

%.fa.gz: %.fa
	bgzip < $< > $@

%.fa.gz.csi: %.fa.gz
	samtools faidx $<
