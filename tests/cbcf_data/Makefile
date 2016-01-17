VCF=$(wildcard *.vcf)

VCFGZ=$(VCF:%.vcf=%.vcf.gz)

all: $(VCFGZ)

%.vcf.gz: %.vcf
	bgzip < $< > $@
	tabix -p vcf $@

clean:
	rm -f *.gz *.tbi
