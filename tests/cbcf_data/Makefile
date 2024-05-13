ALL_VCF=$(wildcard *.vcf)
VCF=$(filter-out example_empty.vcf,$(ALL_VCF))

VCFGZ=$(VCF:%.vcf=%.vcf.gz)
BCF=$(VCF:%.vcf=%.bcf)

all: all.stamp

all.stamp: $(VCFGZ) $(BCF)
	touch $@

%.vcf.gz: %.vcf
	bgzip < $< > $@
	tabix -f -p vcf $@    # create tbi index
	bcftools index -f $@  # create csi index

%.bcf: %.vcf.gz
	bcftools view -O b $< -o $@
	bcftools index -f $@

example_empty.bcf: example_empty.vcf.gz
	touch $@

clean:
	-rm -f all.stamp *.gz *.tbi *.csi *.bcf
