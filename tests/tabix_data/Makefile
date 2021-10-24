all: all.stamp

DERIVED_FILES = \
	fivecolumns.bed.gz fivecolumns.bed.gz.tbi

all.stamp: $(DERIVED_FILES)
	touch $@

%.gz: %
	bgzip -c $< > $@

%.gz.tbi: %.gz
	tabix -p $(subst .,,$(suffix $*)) $<

clean:
	-rm -f all.stamp $(DERIVED_FILES)
