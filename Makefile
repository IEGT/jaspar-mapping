
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

.PHONY: test all depend

depend: .depend

.depend: $(SRCS)
	rm -f "$@"
	$(CC) $(CFLAGS) -MM $^ -MF "$@"

include .depend
