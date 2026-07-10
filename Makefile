
CXX=g++
CXXFLAGS?= -std=c++23
CXXFLAGS += -g
CXXFLAGS += -O3
CXXFLAGS += -DNDEBUG
CXXOPTFLAGS?=
CXXFLAGS += $(CXXOPTFLAGS)
TEST_CXXFLAGS?=$(filter-out -std=c++23 -O3 -g,$(CXXFLAGS)) -std=c++17 -O0
LTO?=0
ifeq ($(LTO),1)
CXXFLAGS += -flto
endif
PKG_CONFIG?=pkg-config
PARQUET?=0
PARQUET_CXXFLAGS?=$(shell $(PKG_CONFIG) --cflags arrow parquet 2>/dev/null)
PARQUET_LDFLAGS?=$(shell $(PKG_CONFIG) --libs arrow parquet 2>/dev/null)
ifeq ($(PARQUET),1)
CXXFLAGS += -DPSSM_SCAN_WITH_PARQUET $(PARQUET_CXXFLAGS)
endif
#LDFLAGS=-lz -lbz2
#LDFLAGS=/home/sm718/miniconda3/pkgs/zlib-1.3.1-h4ab18f5_1/lib/libz.a /home/sm718/miniconda3/pkgs/bzip2-1.0.8-h4bc722e_7/lib/libbz2.a
LDFLAGS += -lz -lbz2
LDFLAGS += -lm
LDOPTFLAGS?=
LDFLAGS += $(LDOPTFLAGS)
ifeq ($(LTO),1)
LDFLAGS += -flto
endif
ifeq ($(PARQUET),1)
LDFLAGS += $(PARQUET_LDFLAGS)
endif

SCRATCHDIR=/tmp
SRCS=$(wildcard *.cpp)
CHR=unset
SCORE_MODE?=log2_relative_risk
OUTPUTDIR?=$(if $(filter log2_relative_risk,$(SCORE_MODE)),output_RelativeRisk_20250217,output_$(SCORE_MODE)_20250217)
PSEUDOCOUNT?=
PSEUDOCOUNT_FLAG=$(if $(PSEUDOCOUNT),--pseudocount $(PSEUDOCOUNT),)
SCAN_CHR_FLAGS?=--coordinate-mode bed --show-sequence $(PSEUDOCOUNT_FLAG)

JASPAR_VERSION=2026
JASPAR_DIR=.
JASPAR_BASENAME=JASPAR$(JASPAR_VERSION)_CORE_non-redundant_pfms_jaspar.txt
JASPAR=$(JASPAR_DIR)/$(JASPAR_BASENAME)
JASPAR_URL=https://jaspar.elixir.no/download/data/$(JASPAR_VERSION)/CORE/$(JASPAR_BASENAME)
GENOME=Homo_sapiens.GRCh38.dna.primary_assembly.fasta
GENOMEGZ=Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
GENOME_INDEX=$(GENOME).fai
TP73_MOTIF_ID?=$(if $(filter 2022,$(JASPAR_VERSION)),MA0861.1,MA0861.2)
TP73_MOTIF=TP73_$(TP73_MOTIF_ID)
TP73_BIDIRECT=$(TP73_MOTIF)_bidirect_$(CHR)
TP73_BIDIRECT_GZ=$(OUTPUTDIR)/$(CHR)/$(TP73_BIDIRECT).bed.gz
TP73_COMBINED_BED=$(TP73_BIDIRECT).combined.bed
TEST_REFERENCE_OUTPUTDIR=output_test_reference_chr1
TEST_REFERENCE_CHR=1
TEST_REFERENCE_FROM=3651800
TEST_REFERENCE_TO=3652600
TEST_REFERENCE_THRESHOLD=0
TEST_REFERENCE_MOTIFS=MA0024.3 MA0106.3 MA0525.2 $(TP73_MOTIF_ID)
TEST_REFERENCE_E2F1_POSITIVE=$(TEST_REFERENCE_OUTPUTDIR)/E2F1_MA0024.3_positive_$(TEST_REFERENCE_CHR)_$(TEST_REFERENCE_FROM)-$(TEST_REFERENCE_TO).bed
TEST_REFERENCE_E2F1_NEGATIVE=$(TEST_REFERENCE_OUTPUTDIR)/E2F1_MA0024.3_negative_$(TEST_REFERENCE_CHR)_$(TEST_REFERENCE_FROM)-$(TEST_REFERENCE_TO).bed
DISTRIBUTION_CHR=1
DISTRIBUTION_BIN_WIDTH?=adaptive
DISTRIBUTION_PSEUDOCOUNT_LABEL=$(if $(PSEUDOCOUNT),_pseudocount_$(PSEUDOCOUNT),)
DISTRIBUTIONDIR?=score_distributions_$(SCORE_MODE)_bins_$(DISTRIBUTION_BIN_WIDTH)_JASPAR$(JASPAR_VERSION)_chr$(DISTRIBUTION_CHR)$(DISTRIBUTION_PSEUDOCOUNT_LABEL)

.SUFFIXES: .gz .bed.gz .cpp .o .fasta .fa.gz _positive_$(CHR).bed _positive_$(CHR).bed.gz _negative_$(CHR).bed _negative_$(CHR).bed.gz _bidirect_$(CHR).bed.gz .bed .bedGraph .combined.bed

BINARIES=pssm_scan gtf_file_region_retrieval context
PARQUET_BINARIES=pssm_scan_parquet
TEST_BINARIES=tests/test_pssm_scan tests/test_gtf_file_region

all: $(BINARIES)

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $<

compressed_file_reader.o: compressed_file_reader.cpp compressed_file_reader.h

progress.o: progress.cpp progress.h

pssm.o: pssm.cpp pssm.h

pssm_scan_core.o: pssm_scan_core.cpp pssm_scan_core.h pssm.h

gtf_file_region.o: gtf_file_region.cpp gtf_file_region.h progress.h

context: context.o compressed_file_reader.o
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDFLAGS)

pssm_scan: pssm_scan.cpp pssm.h pssm_scan_core.h progress.o pssm.o pssm_scan_core.o compressed_file_reader.o
	$(CXX) $(CXXFLAGS) -o $@ pssm_scan.cpp progress.o pssm.o pssm_scan_core.o compressed_file_reader.o $(LDFLAGS)

pssm_scan_parquet: pssm_scan.cpp pssm.h pssm_scan_core.h progress.o pssm.o pssm_scan_core.o compressed_file_reader.o
	$(CXX) $(CXXFLAGS) -DPSSM_SCAN_WITH_PARQUET $(PARQUET_CXXFLAGS) -o $@ pssm_scan.cpp progress.o pssm.o pssm_scan_core.o compressed_file_reader.o $(LDFLAGS) $(PARQUET_LDFLAGS)

tests/test_pssm_scan: tests/test_pssm_scan.cpp pssm_scan_core.h pssm.h pssm.o pssm_scan_core.o
	$(CXX) $(TEST_CXXFLAGS) -o $@ tests/test_pssm_scan.cpp pssm.o pssm_scan_core.o

tests/test_gtf_file_region: tests/test_gtf_file_region.cpp gtf_file_region.h gtf_file_region.o progress.o
	$(CXX) $(TEST_CXXFLAGS) -o $@ tests/test_gtf_file_region.cpp gtf_file_region.o progress.o

check: $(TEST_BINARIES)
	./tests/test_pssm_scan
	./tests/test_gtf_file_region

check-r:
	Rscript tests/test_analyze_bed_cutandrun.R

gtf_file_region_retrieval: gtf_file_region_retrieval.cpp progress.o gtf_file_region.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

clean:
	$(RM) -f *.o

distclean:
	$(RM) -f $(BINARIES) $(PARQUET_BINARIES) $(TEST_BINARIES) *.o

$(JASPAR):
	mkdir -p $(dir $@)
	wget -O $@ $(JASPAR_URL)
jaspar: $(JASPAR)

$(GENOMEGZ):
	wget https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/$(GENOMEGZ)

$(GENOME): $(GENOMEGZ)
	gunzip -c $< > $@

genome: $(GENOME)
genomegz: $(GENOMEGZ)

$(GENOME_INDEX): $(GENOME)
	samtools faidx $<

genome_index: $(GENOME_INDEX)

scan_chr_all_motifs: pssm_scan $(JASPAR) $(GENOME) $(GENOME_INDEX)
	mkdir -p $(OUTPUTDIR)/$(CHR)
	./pssm_scan --outdir $(OUTPUTDIR)/$(CHR) --genome $(GENOME) --pssm $(JASPAR) --score-mode $(SCORE_MODE) -l 0 --chr $(CHR) $(SCAN_CHR_FLAGS)

# Define the pattern rule for generating .bed files
$(OUTPUTDIR)/$(CHR)/%_negative_$(CHR).bed $(OUTPUTDIR)/$(CHR)/%_positive_$(CHR).bed:
	mkdir -p $(OUTPUTDIR)/$(CHR)
	echo $*
	NAME=$(shell echo $* | sed -e 's/_MA.*$$//') ; \
	ACC=$(shell echo $* | tr "_" "\n" |grep -E "^MA[0-9][0-9][0-9][0-9]" |head -n 1) ; \
	echo "NAME=$$NAME ACC=$$ACC" ; \
	if [ ! -f $(OUTPUTDIR)/$(CHR)/$${NAME}_$${ACC}_positive_$(CHR).bed ] ; then \
	    echo "Missing: $(OUTPUTDIR)/$(CHR)/$${NAME}_$${ACC}_positive_$(CHR).bed" ; \
	    ./pssm_scan --outdir $(OUTPUTDIR)/$(CHR) --genome $(GENOME) --pssm $(JASPAR) --score-mode $(SCORE_MODE) -l 0 -m $$ACC --chr $(CHR) ; \
	elif [ ! -f $(OUTPUTDIR)/$(CHR)/$${NAME}_$${ACC}_negative_$(CHR).bed ]; then \
	    echo "Missing: sort -k 1,1 -k2,2n" ; \
	    ./pssm_scan --outdir $(OUTPUTDIR)/$(CHR) --genome $(GENOME) --pssm $(JASPAR) --score-mode $(SCORE_MODE) -l 0 -m $$ACC --chr $(CHR) ; \
	fi

$(OUTPUTDIR)/%_bidirect_$(CHR).bed.gz: $(OUTPUTDIR)/%_negative_$(CHR).bed.gz $(OUTPUTDIR)/%_positive_$(CHR).bed.gz
	zcat $^ | sort -S 2G -k 1,1 -k2,2n | gzip -9 -n -c > $@

%.bed.gz: %.bed
	gzip -9 -n $<

# Generate the list of targets
SHELL=bash
BED_FILES := $(shell if [ -f "$(JASPAR)" ]; then grep "^>" "$(JASPAR)" | sed -e 's%[/:()]%-%g' | awk '{print $$NF "_" $$1 "_positive_$(CHR).bed.gz"}' | sed -e 's/[>]//'; fi)
BIDIRECT_FILES := $(shell if [ -f "$(JASPAR)" ]; then grep "^>" "$(JASPAR)" | sed -e 's%[/:()]%-%g' | awk '{print $$NF "_" $$1 "_bidirect_$(CHR).bed.gz"}' | sed -e 's/[>]//'; fi)
SCORE_DISTRIBUTION_FILES := $(shell if [ -f "$(JASPAR)" ]; then grep "^>" "$(JASPAR)" | sed -e 's%[/:()]%-%g' | awk '{print "$(DISTRIBUTIONDIR)/" $$NF "_" $$1 "_score_distribution_$(SCORE_MODE)_bins_$(DISTRIBUTION_BIN_WIDTH)_positive_$(DISTRIBUTION_CHR)$(DISTRIBUTION_PSEUDOCOUNT_LABEL).tsv"}' | sed -e 's/[>]//'; fi)
echo_bed:
	@echo $(BED_FILES)|sort -S 2G
echo_bidirect:
	@echo $(BIDIRECT_FILES)|sort -S 2G
echo_score_distributions:
	@echo $(SCORE_DISTRIBUTION_FILES)|sort -S 2G

score_distributions_chr1: $(SCORE_DISTRIBUTION_FILES)

$(DISTRIBUTIONDIR)/%_score_distribution_$(SCORE_MODE)_bins_$(DISTRIBUTION_BIN_WIDTH)_positive_$(DISTRIBUTION_CHR)$(DISTRIBUTION_PSEUDOCOUNT_LABEL).tsv: pssm_scan $(JASPAR) $(GENOME)
	mkdir -p $(DISTRIBUTIONDIR)
	echo $*
	NAME=$(shell echo $* | sed -e 's/_MA.*$$//') ; \
	ACC=$(shell echo $* | tr "_" "\n" |grep -E "^MA[0-9][0-9][0-9][0-9]" |head -n 1) ; \
	echo "NAME=$$NAME ACC=$$ACC" ; \
	./pssm_scan --score-distribution --distribution-bin-width $(DISTRIBUTION_BIN_WIDTH) --outdir $(DISTRIBUTIONDIR) --genome $(GENOME) --pssm $(JASPAR) --score-mode $(SCORE_MODE) $(PSEUDOCOUNT_FLAG) --motif $$ACC --chr $(DISTRIBUTION_CHR)

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
	#./pssm_scan --genome Homo_sapiens.GRCh38.dna.primary_assembly_top500000.fasta -l -5400 --verbose -m $(TP73_MOTIF_ID) --chr $(CHR) --from 100000 --to 103000 --help
	./pssm_scan --pssm $(JASPAR) --genome Homo_sapiens.GRCh38.dna.primary_assembly_top500000.fasta --score-mode $(SCORE_MODE) -l 0 --verbose -m MA0059.1 --chr $(CHR) --from 100001 --to 103001 -s
	#./pssm_scan --genome Homo_sapiens.GRCh38.dna.primary_assembly_top500000.fasta -l 0 --verbose -m MA0019.1 --chr $(CHR) --from 100000 --to 130000
	#./pssm_scan --genome Homo_sapiens.GRCh38.dna.primary_assembly_top500000.fasta -l -500 --verbose -m MA1001.3 --chr $(CHR) --from 100000 --to 130000
	#./pssm_scan --genome Homo_sapiens.GRCh38.dna.primary_assembly_top500000.fasta -l -500 --verbose -m $(TP73_MOTIF_ID) --chr $(CHR) --from 100001 --to 103001 -s
	#./pssm_scan --genome Homo_sapiens.GRCh38.dna.primary_assembly_bottom500000.fasta -l -500 --verbose -o output_bottom --chr 44 --from 100000 --to 103000
	#./pssm_scan --genome Homo_sapiens.GRCh38.dna.primary_assembly_bottom500000.fasta -l -500 --verbose -o output_bottom --from 100000 --to 103000

test_reference_tp73_promoter_chr1: pssm_scan $(JASPAR) $(GENOME)
	@if [ "$(TEST_REFERENCE_THRESHOLD)" != "0" ]; then \
	    echo "E: test_reference_tp73_promoter_chr1 expects TEST_REFERENCE_THRESHOLD=0 because pssm_scan changes output filenames for non-zero thresholds." ; \
	    exit 1 ; \
	fi
	mkdir -p $(TEST_REFERENCE_OUTPUTDIR)
	for motif in $(TEST_REFERENCE_MOTIFS); do \
	    ./pssm_scan --outdir $(TEST_REFERENCE_OUTPUTDIR) --genome $(GENOME) --pssm $(JASPAR) --score-mode $(SCORE_MODE) --motif $$motif --chr $(TEST_REFERENCE_CHR) --from $(TEST_REFERENCE_FROM) --to $(TEST_REFERENCE_TO) --threshold $(TEST_REFERENCE_THRESHOLD) --show-sequence ; \
	done
	for file in "$(TEST_REFERENCE_E2F1_POSITIVE)" "$(TEST_REFERENCE_E2F1_NEGATIVE)"; do \
	    if [ ! -f "$$file" ]; then \
	        echo "E: Expected E2F1 reference output '$$file' was not created." ; \
	        exit 1 ; \
	    fi ; \
	done ; \
	hits=$$(awk 'FNR > 1 { count++ } END { print count + 0 }' "$(TEST_REFERENCE_E2F1_POSITIVE)" "$(TEST_REFERENCE_E2F1_NEGATIVE)") ; \
	if [ "$$hits" -lt 1 ]; then \
	    echo "E: Expected E2F1 hits in TP73 promoter reference window; see PMC1531693." ; \
	    exit 1 ; \
	fi ; \
	echo "I: E2F1 hits in TP73 promoter reference window: $$hits"

.PHONY: test check check-r all $(OUTPUTDIR)/$(CHR) jaspar genome genomegz genome_index scan_chr_all_motifs genome_testdata count datatables files_cutandrun_clean TP73_datatable test_reference_tp73_promoter_chr1 score_distributions_chr1
.PRECIOUS: $(GENOME) $(GENOMEGZ) %.bed %.bed.gz
.SECONDARY:

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

#test2: $(TP73_COMBINED_BED).gz
#$(TP73_COMBINED_BED).gz: $(TP73_BIDIRECT_GZ)
#	ls $(TP73_BIDIRECT_GZ)

%_$(CHR).combined.bed: $(TP73_BIDIRECT_GZ) $(FILES_CUTNRUN)
	if ! which bedtools; then echo "E: Need bedtools in path."; exit 1; fi

	#echo -n "Chr\tFrom\tTo\tName\tScore\tStrand" > "$$a_tmp"
	#cp $$i $$i_tmp

	a_tmp=$(shell mktemp -p . -u --suffix="_a_tmp_Chr_$(CHR).bed") ; \
	b_tmp=$(shell mktemp -p . -u --suffix="_b_tmp_Chr_$(CHR).bed") ; \
	i_tmp=$(shell mktemp -p . -u --suffix="_i_tmp_Chr_$(CHR).bed") ; \
	TP73_refgz="$(TP73_BIDIRECT_GZ)" ; \
	zcat "$$TP73_refgz" | grep -vi ^Chr >> "$$a_tmp" ; \
	outputfile="$(TP73_COMBINED_BED)" ; \
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

TP73_datatable_$(CHR).bed.gz: $(TP73_COMBINED_BED).gz
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
