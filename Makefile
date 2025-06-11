
CXX=g++
CXXFLAGS= -std=c++23
CXXFLAGS += -I/home/sm718/miniconda3/include
CXXFLAGS += -g
CXXFLAGS += -O3
#LDFLAGS=-lz -lbz2
#LDFLAGS=/home/sm718/miniconda3/pkgs/zlib-1.3.1-h4ab18f5_1/lib/libz.a /home/sm718/miniconda3/pkgs/bzip2-1.0.8-h4bc722e_7/lib/libbz2.a
LDFLAGS += -lz -lbz2
LDFLAGS += -lm

SCRATCHDIR=/tmp
SRCS=$(wildcard *.cpp)
CHR=unset
OUTPUTDIR=output_RelativeRisk_20250217

JASPAR=JASPAR2022_CORE_non-redundant_pfms_jaspar.txt
GENOME=Homo_sapiens.GRCh38.dna.primary_assembly.fasta
GENOMEGZ=Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

.SUFFIXES: .gz .bed.gz .cpp .o .fasta .fa.gz _positive_$(CHR).bed _positive_$(CHR).bed.gz _negative_$(CHR).bed _negative_$(CHR).bed.gz _bidirect_$(CHR).bed.gz .bed .bedGraph .combined.bed

BINARIES=pssm_scan gtf_file_region_retrieval context

all: $(BINARIES)

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $<

context: context.o compressed_file_reader.o
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDFLAGS)

pssm_scan: pssm_scan.cpp progress.o pssm.o compressed_file_reader.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

gtf_file_region_retrieval: gtf_file_region_retrieval.cpp progress.o gtf_file_region.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

clean:
	$(RM) -f *.o

distclean:
	$(RM) -f $(BINARIES) *.o

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
$(OUTPUTDIR)/$(CHR)/%_negative_$(CHR).bed $(OUTPUTDIR)/$(CHR)/%_positive_$(CHR).bed:
	echo $*
	NAME=$(shell echo $* | sed -e 's/_MA.*$$//') ; \
	ACC=$(shell echo $* | tr "_" "\n" |grep -E "^MA[0-9][0-9][0-9][0-9]" |head -n 1) ; \
	echo "NAME=$$NAME ACC=$$ACC" ; \
	if [ ! -f $(OUTPUTDIR)/$(CHR)/$${NAME}_$${ACC}_positive_$(CHR).bed ] ; then \
	    echo "Missing: $(OUTPUTDIR)/$(CHR)/$${NAME}_$${ACC}_positive_$(CHR).bed" ; \
	    ./pssm_scan --outdir $(OUTPUTDIR)/$(CHR) --genome $(GENOME) -l 0 -m $$ACC --chr $(CHR) ; \
	elif [ ! -f $(OUTPUTDIR)/$(CHR)/$${NAME}_$${ACC}_negative_$(CHR).bed ]; then \
	    echo "Missing: sort -k 1,1 -k2,2n" ; \
	   ./pssm_scan --outdir $(OUTPUTDIR)/$(CHR) --genome $(GENOME) -l 0 -m $$ACC --chr $(CHR) ; \
	fi

#$(OUTPUTDIR)/%_bidirect_$(CHR).bed.gz: $(OUTPUTDIR)/%_negative_$(CHR).bed.gz $(OUTPUTDIR)/%_positive_$(CHR).bed.gz
#	zcat $^ | sort -S 2G -k 1,1 -k2,2n | gzip -9 -n -c > $@

%.bed.gz: %.bed
	gzip -9 -n $<

# Generate the list of targets
SHELL=bash
BED_FILES := $(shell grep "^>" JASPAR2022_CORE_non-redundant_pfms_jaspar.txt | sed -e 's%[/:()]%-%g' | awk '{print $$NF "_" $$1 "_positive_$(CHR).bed.gz"}' | sed -e 's/[>]//')
BIDIRECT_FILES := $(shell grep "^>" JASPAR2022_CORE_non-redundant_pfms_jaspar.txt | sed -e 's%[/:()]%-%g' | awk '{print $$NF "_" $$1 "_bidirect_$(CHR).bed.gz"}' | sed -e 's/[>]//')
echo_bed:
	@echo $(BED_FILES)|sort -S 2G
echo_bidirect:
	@echo $(BIDIRECT_FILES)|sort -S 2G

$(shell basename $(GENOME) .fasta )_top500000.fasta: $(GENOME)
	head -n 500000 $< > Homo_sapiens.GRCh38.dna.primary_assembly_top500000.fasta
$(shell basename $(GENOME) .fasta )_bottom500000.fasta: $(GENOME)
	( echo ">44 nonsense" ; tail -n 500000 $< ) > Homo_sapiens.GRCh38.dna.primary_assembly_bottom500000.fasta

genome_testdata: $(shell basename $(GENOME) .fasta )_bottom500000.fasta $(shell basename $(GENOME) .fasta )_top500000.fasta
genome_testdata_gz: genome_testdata
	gzip -9 -n -k $(shell basename $(GENOME) .fasta )_bottom500000.fasta $(shell basename $(GENOME) .fasta )_top500000.fasta

testGTF: gtf_file_region_retrieval
	echo "TP73" |  ./gtf_file_region_retrieval

test: pssm_scan Homo_sapiens.GRCh38.dna.primary_assembly_top500000.fasta Homo_sapiens.GRCh38.dna.primary_assembly_bottom500000.fasta
	#./pssm_scan --genome Homo_sapiens.GRCh38.dna.primary_assembly_top500000.fasta -l -5400 --verbose -m MA0861.1 --chr $(CHR) --from 100000 --to 103000 --help
	./pssm_scan --genome Homo_sapiens.GRCh38.dna.primary_assembly_top500000.fasta -l 0 --verbose -m MA0059.1 --chr $(CHR) --from 100001 --to 103001 -s
	#./pssm_scan --genome Homo_sapiens.GRCh38.dna.primary_assembly_top500000.fasta -l 0 --verbose -m MA0019.1 --chr $(CHR) --from 100000 --to 130000
	#./pssm_scan --genome Homo_sapiens.GRCh38.dna.primary_assembly_top500000.fasta -l -500 --verbose -m MA1001.3 --chr $(CHR) --from 100000 --to 130000
	#./pssm_scan --genome Homo_sapiens.GRCh38.dna.primary_assembly_top500000.fasta -l -500 --verbose -m MA0861.1 --chr $(CHR) --from 100001 --to 103001 -s
	#./pssm_scan --genome Homo_sapiens.GRCh38.dna.primary_assembly_bottom500000.fasta -l -500 --verbose -o output_bottom --chr 44 --from 100000 --to 103000
	#./pssm_scan --genome Homo_sapiens.GRCh38.dna.primary_assembly_bottom500000.fasta -l -500 --verbose -o output_bottom --from 100000 --to 103000

.PHONY: test all $(OUTPUTDIR)/$(CHR) jaspar genome genomegz genome_testdata count datatables files_cutandrun_clean TP73_datatable
.PRECIOUS: $(GENOME) $(GENOMEGZ)

#PATH_CUTNRUN=cutandrun_20240313_nodupes
#PATH_CUTNRUN=cutandrun_20250516_withDuplicates
PATH_CUTNRUN=cutandrun_20250602_noDuplicates
FILES_CUTNRUN= $(PATH_CUTNRUN)/pos_saos2_DN_R1.clipped.clean.bed \
	$(PATH_CUTNRUN)/pos_saos2_GFP_R1.clipped.clean.bed \
	$(PATH_CUTNRUN)/pos_saos2_TA_R1.clipped.clean.bed \
	$(PATH_CUTNRUN)/pos_skmel29_2_DN_R1.clipped.clean.bed \
	$(PATH_CUTNRUN)/pos_skmel29_2_GFP_R1.clipped.clean.bed \
	$(PATH_CUTNRUN)/pos_skmel29_2_TA_R1.clipped.clean.bed \
	$(PATH_CUTNRUN)/tp73_saos2_DN_R1.clipped.clean.bed \
	$(PATH_CUTNRUN)/tp73_saos2_GFP_R1.clipped.clean.bed \
	$(PATH_CUTNRUN)/tp73_saos2_TA_R1.clipped.clean.bed \
	$(PATH_CUTNRUN)/tp73_skmel29_2_DN_R1.clipped.clean.bed \
	$(PATH_CUTNRUN)/tp73_skmel29_2_GFP_R1.clipped.clean.bed \
	$(PATH_CUTNRUN)/tp73_skmel29_2_TA_R1.clipped.clean.bed

# Derives .bed files from the bedGraphs
files_cutandrun_clean: $(FILES_CUTNRUN)

#test2: $(OUTPUTDIR)/$(CHR)/TP73_MA0861.1_bidirect_$(CHR).combined.bed.gz
#$(OUTPUTDIR)/$(CHR)/TP73_MA0861.1_bidirect_$(CHR).combined.bed.gz: $(OUTPUTDIR)/$(CHR)/TP73_MA0861.1_bidirect_$(CHR).bed.gz
#	ls $(OUTPUTDIR)/$(CHR)/TP73_MA0861.1_bidirect_$(CHR).bed.gz

%_$(CHR).combined.bed: $(OUTPUTDIR)/$(CHR)/TP73_MA0861.1_bidirect_$(CHR).bed.gz $(FILES_CUTNRUN)
	if ! which bedtools; then echo "E: Need bedtools in path."; exit 1; fi

	#echo -n "Chr\tFrom\tTo\tName\tScore\tStrand" > "$$a_tmp"
	#cp $$i $$i_tmp

	a_tmp=$(shell mktemp -p . -u --suffix="_a_tmp_Chr_$(CHR).bed") ; \
	b_tmp=$(shell mktemp -p . -u --suffix="_b_tmp_Chr_$(CHR).bed") ; \
	i_tmp=$(shell mktemp -p . -u --suffix="_i_tmp_Chr_$(CHR).bed") ; \
	TP73_refgz="$(OUTPUTDIR)/$(CHR)/TP73_MA0861.1_bidirect_$(CHR).bed.gz" ; \
	zcat "$$TP73_refgz" | grep -vi ^Chr >> "$$a_tmp" ; \
	outputfile="TP73_MA0861.1_bidirect_$(CHR).combined.bed" ; \
	echo -e -n "Chr\tFrom\tTo\tName\tScore\tStrand" > "$$outputfile" ; \
	\
	for i in $(FILES_CUTNRUN); do \
		echo "I: Working on '$$i' for '$$outputfile'" ; \
		echo -n "	" >> "$$outputfile" ; \
		echo "I:    Cleaned scientific notations" ; \
		echo -n "$$i" | sed -e 's/_R1.clipped.clean.bed$$//' >> "$$outputfile" ; \
		if [ -f "$$a_tmp" ]; then \
			if ! LANG=C awk 'BEGIN {OFS="\t"} {for (i=1; i<=NF; i++) if ($$i ~ /^[0-9.eE+-]+$$/) $$i = sprintf("%.0f", $$i); print}' < $$i > "$$i_tmp" ; then \
				echo "awk failed" ; \
				exit 1 ; \
			fi ; \
			if ! LANG=C bedtools map -null 0 -a "$$a_tmp" -b "$$i_tmp" -o min > "$$b_tmp" ; then \
				echo "bedtools failed" ; \
				exit 1 ; \
			fi ; \
			mv "$$b_tmp" "$$a_tmp" ; \
		else \
	        echo "E: File '$$a_tmp' should be existing" ; \
	        exit ; \
		fi ; \
	done ; \
	echo >> "$$outputfile" ; \
	\
	cat "$$a_tmp" >> "$$outputfile" ; \
	rm "$$a_tmp" "$$i_tmp"

# Rule to generate clean.bed files from bedGraph files in PATH_CUTNRUN
$(PATH_CUTNRUN)/%.clean.bed: $(PATH_CUTNRUN)/%.bedGraph
	@echo "$< -> $@"
	grep -v hromosome $< | sed -e 's%^chr%%' | awk '{print $$1"\t"$$2"\t"$$3"\tcnr\t"$$4}' > $@


#$(CHR): $(addprefix $(CHR)/,$(BED_FILES))
$(OUTPUTDIR)/$(CHR): $(addprefix $(OUTPUTDIR)/$(CHR)/,$(BIDIRECT_FILES))

datatables: context
	#for chr in 2; do
	@for chr in $(shell seq 1 22) X Y; do \
	    echo "Chr $$chr" ; \
	    $(MAKE) CHR=$$chr TP73_datatable.bed.gz ; \
	done

TP73_datatable.bed.gz: TP73_datatable_$(CHR).bed.gz

TP73_datatable_$(CHR).bed.gz: TP73_MA0861.1_bidirect_$(CHR).combined.bed.gz
	./context $< $(OUTPUTDIR)/$(CHR)/*bidirect*.bed.gz | gzip -9 -n -c > $@ || echo "I: Check ulimit -n 3000 if failing to open files"


count:
	@if [ -z "$(CHR)" ] || [ "unset" = "$(CHR)" ]; then \
	    echo "E: Define CHR to be on of 1 .. 22 X Y" ; \
	else \
	    find $(OUTPUTDIR)/$(CHR) -name "*_bidirect_*.bed.gz" | wc -l ; \
	fi

#depend: .depend

#.depend: $(SRCS)
#	rm -f "$@"
#	$(CC) $(CFLAGS) -MM $^ -MF "$@"

#include .depend

