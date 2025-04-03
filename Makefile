SHELL := /bin/zsh

##############PARAMETERS##############
pwd=$(shell pwd)
base ?= /Users/shehzailabbas/Desktop/FUSE
file ?= test.tar.gz
out ?= output
snp ?= "Marker Name"
a1 ?= A1
a2 ?= A2
p ?= "P-value"
frq ?= "Allele Frequency in meta-analysis"
signedSumstats ?= "Log(Odds Ratio)"
N?= 387856
nullV?=0

Rscript ?= $(base)/fusion_twas-master/FUSION.assoc_test.R
overlapScript ?= $(base)/fusion_twas-master/overlap.R
sumstats ?= output.sumstats 
weights ?= $(base)/WEIGHTS/GTEx.Whole_Blood.pos 
weights_dir ?= $(base)/WEIGHTS/ 
ref_ld_chr ?= $(base)/fusion_twas-master/LDREF/1000G.EUR. 
chr ?= -1
outF ?= out.dat
SupMultiStudy ?= 0
spacer=-----------------------------------------------------------------
outFolder=OUTPUT
Makeloc=$(base)/Makefile

##############PARAMETERS##############

CONDA_PREFIX = $(shell conda info --base)
.SILENT:
.ONESHELL:
.EXPORT_ALL_VARIABLES:

master:
	@if [ "$(chr)" = "-1" ]; then \
		mkdir -p $(outFolder)/temp $(outFolder)/chromosomes ;\
		echo "Running for all chromosomes..."; \
		for chro in $(shell seq 1 22); do \
			outFile="chr$${chro}_$(outF)"; \
			echo $(spacer)"Processing chr $$chro"$(spacer); \
			if [ -e $(outFolder)/$(out).sumstats ]; then \
    			$(MAKE) -f $(base)/Makefile compSumstats file=$(outFolder)/$(out).sumstats chr=$$chro outF=$(outFolder)/temp/$$outFile; \
				rm -r temp; \
			else \
    			$(MAKE) -f $(base)/Makefile compSumstats chr=$$chro outF=$(outFolder)/temp/$$outFile; \
				mv  temp/$(out).sumstats $(outFolder); \
				rm -r temp; \
			fi; \
		done; \
	else \
		mkdir -p $(outFolder)/chromosomes $(outFolder)/temp; \
		echo $(spacer)"Running for chr=$(chr)"$(spacer); \
		$(MAKE) -f $(base)/Makefile compSumstats outF=$(outFolder)/temp/$(outF); \
		mv temp/$(out).sumstats $(outFolder); \
		rm -r temp; \
	fi

compSumstats:
	set -e
	source $(CONDA_PREFIX)/etc/profile.d/conda.sh && conda activate fiz && echo "starting...." && \
	python $(base)/build.py --file $(file) --out $(out) --snp $(snp) \
	--a1 $(a1) --a2 $(a2) --p $(p) --frq $(frq) --signedSumstats \
	$(signedSumstats) --N $(N) --nullV $(nullV) --Rscript $(Rscript) --sumstats $(sumstats) --weights $(pwd)/$(weights) \
	--weights_dir $(pwd)/$(weights_dir) --ref_ld_chr $(ref_ld_chr) --chr $(chr) \
	--outF $(outF) --SupMultiStudy $(SupMultiStudy)
	make -C temp 


help:
	python $(base)/build.py -h

overlap:
	source $(CONDA_PREFIX)/etc/profile.d/conda.sh && conda activate ldsc & Rscript $(overlapScript) $(weights_dir) $(file) $(weights)

gui:
	echo "edit the parameters in this file $(Makeloc)"

option:
	echo "overlap gui help compSumstats master"
