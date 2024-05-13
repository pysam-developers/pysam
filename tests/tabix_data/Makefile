all: all.stamp

DERIVED_FILES = \
	empty.bed.gz.tbi \
	example.bed.gz example.bed.gz.tbi \
	example.gff2.gz.tbi \
	example.gff3.gz.tbi \
	example.sam.gz.tbi \
	example.vcf.gz.tbi \
	example_badcomments.bed.gz.tbi \
	example_badcomments.gtf.gz.tbi \
	example_badcomments.sam.gz.tbi \
	example_badcomments.vcf.gz.tbi \
	example_comments.bed.gz.tbi \
	example_comments.gtf.gz.tbi \
	example_comments.sam.gz.tbi \
	example_comments.vcf.gz.tbi \
	example_large.bed.gz.tbi \
	fivecolumns.bed.gz fivecolumns.bed.gz.tbi

all.stamp: $(DERIVED_FILES)
	touch $@

%.gz: %
	bgzip -c $< > $@

%.gff2.gz.tbi: %.gff2.gz
	tabix -f -p gff $<

%.gff3.gz.tbi: %.gff3.gz
	tabix -f -p gff $<

%.gtf.gz.tbi: %.gtf.gz
	tabix -f -p gff $<

%.gz.tbi: %.gz
	tabix -f -p $(subst .,,$(suffix $*)) $<

clean:
	-rm -f all.stamp $(DERIVED_FILES)
