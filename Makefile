
CXX=g++
CXXFLAGS=-std=c++23
CXXFLAGS += -g
LDFLAGS=-lm

SRCS=$(wildcard *.cpp)

.SUFFIXES: .cpp .o .fasta .fa.gz

all: depend pssm_scan gtf_file_region_retrieval context

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $<

context: context.o
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

pssm_scan: pssm_scan.cpp progress.o pssm.o
	$(CXX) $(CXXFLAGS) -o $@ $< progress.o pssm.o $(LDFLAGS)

gtf_file_region_retrieval: gtf_file_region_retrieval.cpp progress.o gtf_file_region.o
	$(CXX) $(CXXFLAGS) -o $@ $< progress.o gtf_file_region.o $(LDFLAGS)

clean:
	$(RM) gtf_file_region_retrieval pssm_scan

Homo_sapiens.GRCh38.dna.primary_assembly.fasta: Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
	gunzip -c $< > $@


# Define the pattern rule for generating .bed files
output_Chr1/%_negative_1.bed output_Chr1/%_positive_1.bed: ACC=$(word 1,$(subst _, ,$*))
output_Chr1/%_negative_1.bed output_Chr1/%_positive_1.bed: NAME=$(word 2,$(subst _, ,$*))
output_Chr1/%_negative_1.bed output_Chr1/%_positive_1.bed:
	@echo "NAME=$(NAME) ACC=$(ACC)"
	@if [ ! -f output_Chr1/$(NAME)_$(ACC)_positive_1.bed ] ; then \
	 	echo "Missing: output_Chr1/$(NAME)_$(ACC)_positive_1.bed" ; \
		./pssm_scan --outdir output_Chr1 --genome Homo_sapiens.GRCh38.dna.primary_assembly.fasta -l 0 -m $(ACC) --chr 1 ; \
	elif [ ! -f output_Chr1/$(NAME)_$(ACC)_negative_1.bed ]; then \
	 	echo "Missing: output_Chr1/$(NAME)_$(ACC)_negative_1.bed" ; \
		./pssm_scan --outdir output_Chr1 --genome Homo_sapiens.GRCh38.dna.primary_assembly.fasta -l 0 -m $(ACC) --chr 1 ; \
	fi

# Generate the list of targets
SHELL=bash
BED_FILES := $(shell grep "^>" JASPAR2022_CORE_non-redundant_pfms_jaspar.txt | sed -e 's%[/:()]%-%g' | awk '{print $$1 "_" $$NF "_positive_1.bed"}' | sed -e 's/^>//')
echo:
	echo $(BED_FILES)

Homo_sapiens.GRCh38.dna.primary_assembly_top500000.fasta: Homo_sapiens.GRCh38.dna.primary_assembly.fasta
	head -n 500000 $< > Homo_sapiens.GRCh38.dna.primary_assembly_top500000.fasta
Homo_sapiens.GRCh38.dna.primary_assembly_bottom500000.fasta: Homo_sapiens.GRCh38.dna.primary_assembly.fasta
	( echo ">44 nonsense" ; tail -n 500000 $< ) > Homo_sapiens.GRCh38.dna.primary_assembly_bottom500000.fasta

testGTF: gtf_file_region_retrieval
	echo "TP73" |  ./gtf_file_region_retrieval

test: pssm_scan Homo_sapiens.GRCh38.dna.primary_assembly_top500000.fasta Homo_sapiens.GRCh38.dna.primary_assembly_bottom500000.fasta
	#./pssm_scan --genome Homo_sapiens.GRCh38.dna.primary_assembly_top500000.fasta -l -5400 --verbose -m MA0861.1 --chr 1 --from 100000 --to 103000 --help
	./pssm_scan --genome Homo_sapiens.GRCh38.dna.primary_assembly_top500000.fasta -l 0 --verbose -m MA0861.1 --chr 1 --from 100000 --to 103000
	./pssm_scan --genome Homo_sapiens.GRCh38.dna.primary_assembly_top500000.fasta -l -500 --verbose -m MA0861.1 --chr 1 --from 100000 --to 103000
	./pssm_scan --genome Homo_sapiens.GRCh38.dna.primary_assembly_top500000.fasta -l 0 --verbose -m MA0861.1 --chr 1 --from 100001 --to 103001 -s
	./pssm_scan --genome Homo_sapiens.GRCh38.dna.primary_assembly_top500000.fasta -l -500 --verbose -m MA0861.1 --chr 1 --from 100001 --to 103001 -s
	#./pssm_scan --genome Homo_sapiens.GRCh38.dna.primary_assembly_bottom500000.fasta -l -500 --verbose -o output_bottom --chr 44 --from 100000 --to 103000
	#./pssm_scan --genome Homo_sapiens.GRCh38.dna.primary_assembly_bottom500000.fasta -l -500 --verbose -o output_bottom --from 100000 --to 103000

.PHONY: test all depend output_Chr1

output_Chr1: $(addprefix output_Chr1/,$(BED_FILES))

depend: .depend

.depend: $(SRCS)
	rm -f "$@"
	$(CC) $(CFLAGS) -MM $^ -MF "$@"

include .depend
