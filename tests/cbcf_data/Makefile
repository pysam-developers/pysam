VCF=$(wildcard *.vcf)

VCFGZ=$(VCF:%.vcf=%.vcf.gz)
BCF=$(VCF:%.vcf=%.bcf)

all: $(VCFGZ) $(BCF)

%.vcf.gz: %.vcf
	bgzip < $< > $@
	tabix -p vcf $@    # create tbi index
	bcftools index $@  # create csi index

%.bcf: %.vcf.gz
	bcftools view -O b $< -o $@
	bcftools index $@

example_empty.bcf: example_empty.vcf.gz
	touch $@

clean:
	rm -f *.gz *.tbi *.csi *.bcf

