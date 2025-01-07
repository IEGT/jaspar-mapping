
CXX=g++
CXXFLAGS=-std=c++23
CXXFLAGS += -g
LDFLAGS=-lz -lbz2 -lm 

SCRATCHDIR=/tmp
SRCS=$(wildcard *.cpp)

JASPAR=JASPAR2022_CORE_non-redundant_pfms_jaspar.txt
GENOME=Homo_sapiens.GRCh38.dna.primary_assembly.fasta
GENOMEGZ=Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

.SUFFIXES: .gz .bed.gz .cpp .o .fasta .fa.gz _positive_1.bed _positive_1.bed.gz _negative_1.bed _negative_1.bed.gz _bidirect_1.bed.gz

all: pssm_scan gtf_file_region_retrieval context

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $<

context: context.o compressed_file_reader.o
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDFLAGS)

pssm_scan: pssm_scan.cpp progress.o pssm.o compressed_file_reader.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

gtf_file_region_retrieval: gtf_file_region_retrieval.cpp progress.o gtf_file_region.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

clean:
	$(RM) gtf_file_region_retrieval pssm_scan

$(JASPAR):
	wget https://jaspar2022.genereg.net/download/data/2022/CORE/JASPAR2022_CORE_non-redundant_pfms_jaspar.txt
jaspar: $(JASPAR)

$(GENOMEGZ):
	wget https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/$(GENOMEGZ)

$(GENOME): $(GENOMEGZ)
	gunzip -c $< > $@

genome: $(GENOME)
genomegz: $(GENOMEGZ)

# Define the pattern rule for generating .bed files
output_Chr1/%_negative_1.bed output_Chr1/%_positive_1.bed: NAME=$(shell echo $* | sed -e 's/_MA.*$$//')
output_Chr1/%_negative_1.bed output_Chr1/%_positive_1.bed: ACC=$(shell echo $* | tr "_" "\n" |grep ^MA|head -n 1)
output_Chr1/%_negative_1.bed output_Chr1/%_positive_1.bed:
	@echo "NAME=$(NAME) ACC=$(ACC)"
	if [ ! -f output_Chr1/$(NAME)_$(ACC)_positive_1.bed ] ; then \
	    echo "Missing: output_Chr1/$(NAME)_$(ACC)_positive_1.bed" ; \
	    ./pssm_scan --outdir output_Chr1 --genome $(GENOME) -l 0 -m $(ACC) --chr 1 ; \
	elif [ ! -f output_Chr1/$(NAME)_$(ACC)_negative_1.bed ]; then \
	    echo "Missing: ort -k 1,1 -k2,2n" ; \
		./pssm_scan --outdir output_Chr1 --genome $(GENOME) -l 0 -m $(ACC) --chr 1 ; \
	fi

%_bidirect_1.bed.gz: %_negative_1.bed.gz %_positive_1.bed.gz
	zcat $^ | sort -k 1,1 -k2,2n | gzip -c > $@

%.bed.gz: %.bed
	gzip $<

# Generate the list of targets
SHELL=bash
BED_FILES := $(shell grep "^>" JASPAR2022_CORE_non-redundant_pfms_jaspar.txt | sed -e 's%[/:()]%-%g' | awk '{print $$NF "_" $$1 "_positive_1.bed.gz"}' | sed -e 's/[>]//')
BIDIRECT_FILES := $(shell grep "^>" JASPAR2022_CORE_non-redundant_pfms_jaspar.txt | sed -e 's%[/:()]%-%g' | awk '{print $$NF "_" $$1 "_bidirect_1.bed.gz"}' | sed -e 's/[>]//')
echo_bed:
	@echo $(BED_FILES)|sort
echo_bidirect:
	@echo $(BIDIRECT_FILES)|sort

$(shell basename $(GENOME) .fasta )_top500000.fasta: $(GENOME)
	head -n 500000 $< > Homo_sapiens.GRCh38.dna.primary_assembly_top500000.fasta
$(shell basename $(GENOME) .fasta )_bottom500000.fasta: $(GENOME)
	( echo ">44 nonsense" ; tail -n 500000 $< ) > Homo_sapiens.GRCh38.dna.primary_assembly_bottom500000.fasta

genome_testdata: $(shell basename $(GENOME) .fasta )_bottom500000.fasta $(shell basename $(GENOME) .fasta )_top500000.fasta
genome_testdata_gz: genome_testdata
	gzip -k $(shell basename $(GENOME) .fasta )_bottom500000.fasta $(shell basename $(GENOME) .fasta )_top500000.fasta

testGTF: gtf_file_region_retrieval
	echo "TP73" |  ./gtf_file_region_retrieval

test: pssm_scan Homo_sapiens.GRCh38.dna.primary_assembly_top500000.fasta Homo_sapiens.GRCh38.dna.primary_assembly_bottom500000.fasta
	#./pssm_scan --genome Homo_sapiens.GRCh38.dna.primary_assembly_top500000.fasta -l -5400 --verbose -m MA0861.1 --chr 1 --from 100000 --to 103000 --help
	./pssm_scan --genome Homo_sapiens.GRCh38.dna.primary_assembly_top500000.fasta -l 0 --verbose -m MA0861.1 --chr 1 --from 100000 --to 130000
	./pssm_scan --genome Homo_sapiens.GRCh38.dna.primary_assembly_top500000.fasta -l -500 --verbose -m MA0861.1 --chr 1 --from 100000 --to 130000
	./pssm_scan --genome Homo_sapiens.GRCh38.dna.primary_assembly_top500000.fasta -l 0 --verbose -m MA0861.1 --chr 1 --from 100001 --to 103001 -s
	./pssm_scan --genome Homo_sapiens.GRCh38.dna.primary_assembly_top500000.fasta -l -500 --verbose -m MA0861.1 --chr 1 --from 100001 --to 103001 -s
	#./pssm_scan --genome Homo_sapiens.GRCh38.dna.primary_assembly_bottom500000.fasta -l -500 --verbose -o output_bottom --chr 44 --from 100000 --to 103000
	#./pssm_scan --genome Homo_sapiens.GRCh38.dna.primary_assembly_bottom500000.fasta -l -500 --verbose -o output_bottom --from 100000 --to 103000

.PHONY: test all output_Chr1 jaspar genome genomegz genome_testdata count datatable
.PRECIOUS: $(GENOME) $(GENOMEGZ)

#output_Chr1: $(addprefix output_Chr1/,$(BED_FILES))
output_Chr1: $(addprefix output_Chr1/,$(BIDIRECT_FILES))

TP73_datatable.bed.gz: context output_Chr1/TP73_MA0861.1_bidirect_1.combined.bed.gz
	./context output_Chr1/TP73_MA0861.1_bidirect_1.combined.bed.gz output_Chr1/*bidirect*.bed.gz | gzip -c > $@ || echo "I: Check ulimit -n 3000 if failing to open files"

datatable: TP73_datatable.bed.gz

count:
	find output_Chr1 -name "*_bidirect_*.bed.gz" | wc -l

#depend: .depend

#.depend: $(SRCS)
#	rm -f "$@"
#	$(CC) $(CFLAGS) -MM $^ -MF "$@"

#include .depend
