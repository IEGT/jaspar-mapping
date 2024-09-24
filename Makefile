
CXX=g++
CXXFLAGS=-std=c++23
LDFLAGS=-lm

pssm_scan: pssm_scan.cpp
	$(CXX) $(CXXFLAGS) -o pssm_scan pssm_scan.cpp $(LDFLAGS)

Homo_sapiens.GRCh38.dna.primary_assembly_top500000.fasta: Homo_sapiens.GRCh38.dna.primary_assembly.fasta
	head -n 500000 Homo_sapiens.GRCh38.dna.primary_assembly.fasta > Homo_sapiens.GRCh38.dna.primary_assembly_top500000.fasta
Homo_sapiens.GRCh38.dna.primary_assembly_bottom500000.fasta: Homo_sapiens.GRCh38.dna.primary_assembly.fasta
	( echo ">44 nonsense" ; tail -n 500000 Homo_sapiens.GRCh38.dna.primary_assembly.fasta ) > Homo_sapiens.GRCh38.dna.primary_assembly_bottom500000.fasta

test: pssm_scan Homo_sapiens.GRCh38.dna.primary_assembly_top500000.fasta Homo_sapiens.GRCh38.dna.primary_assembly_bottom500000.fasta
	./pssm_scan --genome Homo_sapiens.GRCh38.dna.primary_assembly_top500000.fasta -l 0 --verbose -m MA0861.1 --chr 1 --from 100000 --to 103000
	./pssm_scan --genome Homo_sapiens.GRCh38.dna.primary_assembly_bottom500000.fasta -l -500 --verbose -o output_bottom --chr 44 --from 100000 --to 103000
	./pssm_scan --genome Homo_sapiens.GRCh38.dna.primary_assembly_bottom500000.fasta -l -500 --verbose -o output_bottom --from 100000 --to 103000

.PHONE: test
