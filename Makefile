
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
file_number_label = $(shell printf '%.12g' '$(1)' | tr '+-' 'pm')
EFFECTIVE_PSEUDOCOUNT_VALUE=$(if $(PSEUDOCOUNT),$(PSEUDOCOUNT),$(if $(filter log_odds,$(SCORE_MODE)),1,0))
EFFECTIVE_PSEUDOCOUNT_LABEL=$(call file_number_label,$(EFFECTIVE_PSEUDOCOUNT_VALUE))
SCAN_THRESHOLD?=0
SCAN_THRESHOLD_LABEL=$(call file_number_label,$(SCAN_THRESHOLD))
SCAN_COORDINATE_MODE?=bed
SCAN_SHOW_SEQUENCE?=1
SCAN_N_POLICY?=skip
MIN_PWM_RELATIVE_SCORE?=
MAX_PWM_RELATIVE_SCORE?=
PWM_RELATIVE_FLAGS=$(if $(MIN_PWM_RELATIVE_SCORE),--min-pwm-relative-score $(MIN_PWM_RELATIVE_SCORE),) $(if $(MAX_PWM_RELATIVE_SCORE),--max-pwm-relative-score $(MAX_PWM_RELATIVE_SCORE),)
SCAN_N_FLAG=$(if $(filter neutral,$(SCAN_N_POLICY)),--neutral-N,--skip-N)
HIT_N_POLICY_LABEL=$(if $(filter neutral,$(SCAN_N_POLICY)),neutral,skip)
SCAN_CHR_FLAGS?=--coordinate-mode $(SCAN_COORDINATE_MODE) $(if $(filter 1,$(SCAN_SHOW_SEQUENCE)),--show-sequence,) $(SCAN_N_FLAG) $(PSEUDOCOUNT_FLAG) $(PWM_RELATIVE_FLAGS)
HIT_SEQUENCE_LABEL=$(if $(filter 1,$(SCAN_SHOW_SEQUENCE)),included,omitted)
HIT_PWM_RELATIVE_MIN_LABEL=$(if $(MIN_PWM_RELATIVE_SCORE),$(call file_number_label,$(MIN_PWM_RELATIVE_SCORE)),none)
HIT_PWM_RELATIVE_MAX_LABEL=$(if $(MAX_PWM_RELATIVE_SCORE),$(call file_number_label,$(MAX_PWM_RELATIVE_SCORE)),none)
HIT_CONFIG_SUBDIR=hits/score_mode-$(SCORE_MODE)/pseudocount-$(EFFECTIVE_PSEUDOCOUNT_LABEL)/threshold-$(SCAN_THRESHOLD_LABEL)/pwm_relative_min-$(HIT_PWM_RELATIVE_MIN_LABEL)/pwm_relative_max-$(HIT_PWM_RELATIVE_MAX_LABEL)/coordinate_mode-$(SCAN_COORDINATE_MODE)/sequence-$(HIT_SEQUENCE_LABEL)/n_policy-$(HIT_N_POLICY_LABEL)
HIT_OUTPUTDIR=$(OUTPUTDIR)/$(CHR)/$(HIT_CONFIG_SUBDIR)

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
TP73_BIDIRECT_GZ=$(HIT_OUTPUTDIR)/$(TP73_BIDIRECT).bed.gz
TP73_COMBINED_BED=$(TP73_BIDIRECT).combined.bed
TEST_REFERENCE_OUTPUTDIR=output_test_reference_chr1
TEST_REFERENCE_CHR=1
TEST_REFERENCE_FROM=3651800
TEST_REFERENCE_TO=3652600
TEST_REFERENCE_THRESHOLD=0
TEST_REFERENCE_THRESHOLD_LABEL=$(call file_number_label,$(TEST_REFERENCE_THRESHOLD))
TEST_REFERENCE_MOTIFS=MA0024.3 MA0106.3 MA0525.2 $(TP73_MOTIF_ID)
TEST_REFERENCE_CONFIG_SUBDIR=hits/score_mode-$(SCORE_MODE)/pseudocount-$(EFFECTIVE_PSEUDOCOUNT_LABEL)/threshold-$(TEST_REFERENCE_THRESHOLD_LABEL)/pwm_relative_min-$(HIT_PWM_RELATIVE_MIN_LABEL)/pwm_relative_max-$(HIT_PWM_RELATIVE_MAX_LABEL)/coordinate_mode-$(SCAN_COORDINATE_MODE)/sequence-$(HIT_SEQUENCE_LABEL)/n_policy-$(HIT_N_POLICY_LABEL)
TEST_REFERENCE_HITDIR=$(TEST_REFERENCE_OUTPUTDIR)/$(TEST_REFERENCE_CONFIG_SUBDIR)
TEST_REFERENCE_E2F1_POSITIVE=$(TEST_REFERENCE_HITDIR)/E2F1_MA0024.3_positive_$(TEST_REFERENCE_CHR)_$(TEST_REFERENCE_FROM)-$(TEST_REFERENCE_TO).bed
TEST_REFERENCE_E2F1_NEGATIVE=$(TEST_REFERENCE_HITDIR)/E2F1_MA0024.3_negative_$(TEST_REFERENCE_CHR)_$(TEST_REFERENCE_FROM)-$(TEST_REFERENCE_TO).bed
DISTRIBUTION_CHR=1
DISTRIBUTION_BIN_WIDTH?=adaptive
DISTRIBUTION_PSEUDOCOUNT_LABEL=$(if $(PSEUDOCOUNT),_pseudocount_$(PSEUDOCOUNT),)
DISTRIBUTIONDIR?=score_distributions_$(SCORE_MODE)_bins_$(DISTRIBUTION_BIN_WIDTH)_JASPAR$(JASPAR_VERSION)_chr$(DISTRIBUTION_CHR)$(DISTRIBUTION_PSEUDOCOUNT_LABEL)
DRY_RUN_FROM?=3600000
DRY_RUN_TO?=3700000
DRY_RUN_FULL_CHR1?=0
DRY_RUN_OUTPUT?=dry_runs/chr1_patz1_tp73_$(if $(filter 1,$(DRY_RUN_FULL_CHR1)),full,from-$(DRY_RUN_FROM)-to-$(DRY_RUN_TO))
DRY_RUN_RANGE_FLAGS=$(if $(filter 1,$(DRY_RUN_FULL_CHR1)),--full-chr1,--from $(DRY_RUN_FROM) --to $(DRY_RUN_TO))
CONTEXT_FLANK?=150
CONTEXT_BATCH_SIZE?=32
CONTEXT_TMPDIR?=$(SCRATCHDIR)
CONTEXT_DATATABLE=TP73_datatable_$(CHR)_flank-$(CONTEXT_FLANK)_motif-center.bed.gz

.SUFFIXES: .gz .bed.gz .cpp .o .fasta .fa.gz _positive_$(CHR).bed _positive_$(CHR).bed.gz _negative_$(CHR).bed _negative_$(CHR).bed.gz _bidirect_$(CHR).bed.gz .bed .bedGraph .combined.bed

BINARIES=pssm_scan gtf_file_region_retrieval context
PARQUET_BINARIES=pssm_scan_parquet
TEST_BINARIES=tests/test_pssm_scan tests/test_gtf_file_region tests/test_compressed_file_reader tests/test_context

all: $(BINARIES)

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $<

compressed_file_reader.o: compressed_file_reader.cpp compressed_file_reader.h

context_core.o: context_core.cpp context_core.h compressed_file_reader.h

context.o: context.cpp context_core.h

progress.o: progress.cpp progress.h

pssm.o: pssm.cpp pssm.h

pssm_scan_core.o: pssm_scan_core.cpp pssm_scan_core.h pssm.h

gtf_file_region.o: gtf_file_region.cpp gtf_file_region.h progress.h

context: context.o context_core.o compressed_file_reader.o
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LDFLAGS)

pssm_scan: pssm_scan.cpp pssm.h pssm_scan_core.h progress.o pssm.o pssm_scan_core.o compressed_file_reader.o
	$(CXX) $(CXXFLAGS) -o $@ pssm_scan.cpp progress.o pssm.o pssm_scan_core.o compressed_file_reader.o $(LDFLAGS)

pssm_scan_parquet: pssm_scan.cpp pssm.h pssm_scan_core.h progress.o pssm.o pssm_scan_core.o compressed_file_reader.o
	$(CXX) $(CXXFLAGS) -DPSSM_SCAN_WITH_PARQUET $(PARQUET_CXXFLAGS) -o $@ pssm_scan.cpp progress.o pssm.o pssm_scan_core.o compressed_file_reader.o $(LDFLAGS) $(PARQUET_LDFLAGS)

tests/test_pssm_scan: tests/test_pssm_scan.cpp pssm_scan_core.h pssm.h pssm.o pssm_scan_core.o
	$(CXX) $(TEST_CXXFLAGS) -o $@ tests/test_pssm_scan.cpp pssm.o pssm_scan_core.o

tests/test_gtf_file_region: tests/test_gtf_file_region.cpp gtf_file_region.h gtf_file_region.o progress.o
	$(CXX) $(TEST_CXXFLAGS) -o $@ tests/test_gtf_file_region.cpp gtf_file_region.o progress.o

tests/test_compressed_file_reader: tests/test_compressed_file_reader.cpp compressed_file_reader.h compressed_file_reader.o
	$(CXX) $(TEST_CXXFLAGS) -o $@ tests/test_compressed_file_reader.cpp compressed_file_reader.o $(LDFLAGS)

tests/test_context: tests/test_context.cpp context_core.h context_core.o compressed_file_reader.o
	$(CXX) $(TEST_CXXFLAGS) -o $@ tests/test_context.cpp context_core.o compressed_file_reader.o $(LDFLAGS)

check: pssm_scan $(TEST_BINARIES)
	./tests/test_pssm_scan
	./tests/test_gtf_file_region
	./tests/test_compressed_file_reader
	./tests/test_context
	bash tests/test_context_cli.sh
	bash tests/test_fix_missing_bidirect.sh
	bash tests/test_localMaxSkmelTADN.sh
	bash tests/test_indexed_genome_scan.sh
	bash tests/test_build_fasta_index.sh
	bash tests/test_cutandrun_containment.sh
	bash tests/test_script_help.sh

check-r:
	Rscript tests/test_analyze_bed_cutandrun.R

check-duckdb: pssm_scan_parquet
	bash tests/test_duckdb_contract.sh
	bash tests/test_chr1_dense_duckdb.sh
	bash tests/test_export_dense_bed.sh
	bash tests/test_dense_cutandrun_calibration.sh
	bash tests/test_motif_context.sh
	bash tests/test_sparse_parquet.sh

check_synthetic_dense: pssm_scan_parquet
	bash tests/test_synthetic_dense_dataset.sh

check_synthetic_sparse: pssm_scan_parquet
	bash tests/test_sparse_parquet.sh

synthetic_dense_example: pssm_scan_parquet
	bash tests/test_synthetic_dense_dataset.sh dry_runs/synthetic_dense

check_cutandrun_containment:
	bash tests/test_cutandrun_containment.sh

synthetic_cutandrun_example:
	python3 scripts/analyze_cutandrun_containment.py \
		--motifs test_files/synthetic_cutandrun/tp73_motifs.bed \
		--coverage-bed test_files/synthetic_cutandrun/tp73_fragments.bed \
		--output-dir dry_runs/synthetic_cutandrun_coverage_union \
		--chrom 1 --sample-id synthetic_tp73 --score-mode synthetic_score

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
	@if command -v samtools >/dev/null 2>&1; then \
		samtools faidx $<; \
	else \
		python3 scripts/build_fasta_index.py $< --output $@; \
	fi

%.fasta.fai: %.fasta
	@if command -v samtools >/dev/null 2>&1; then \
		samtools faidx $<; \
	else \
		python3 scripts/build_fasta_index.py $< --output $@; \
	fi

genome_index: $(GENOME_INDEX)

scan_chr_all_motifs: pssm_scan $(JASPAR) $(GENOME) $(GENOME_INDEX)
	mkdir -p $(OUTPUTDIR)/$(CHR)
	./pssm_scan --outdir $(OUTPUTDIR)/$(CHR) --genome $(GENOME) --pssm $(JASPAR) --score-mode $(SCORE_MODE) --threshold $(SCAN_THRESHOLD) --chr $(CHR) $(SCAN_CHR_FLAGS)

dry_run_chr1_patz1_tp73: pssm_scan_parquet $(JASPAR) $(GENOME) $(GENOME_INDEX)
	bash scripts/run_chr1_patz1_tp73_dry_run.sh --output $(DRY_RUN_OUTPUT) $(DRY_RUN_RANGE_FLAGS)

inspect_chr1_patz1_tp73: dry_run_chr1_patz1_tp73
	bash scripts/inspect_chr1_dense_dry_run.sh --package $(DRY_RUN_OUTPUT) overview

# Define the pattern rule for generating .bed files
$(HIT_OUTPUTDIR)/%_negative_$(CHR).bed $(HIT_OUTPUTDIR)/%_positive_$(CHR).bed: pssm_scan $(JASPAR) $(GENOME) $(GENOME_INDEX)
	mkdir -p $(HIT_OUTPUTDIR)
	echo $*
	NAME=$(shell echo $* | sed -e 's/_MA.*$$//') ; \
	ACC=$(shell echo $* | tr "_" "\n" |grep -E "^MA[0-9][0-9][0-9][0-9]" |head -n 1) ; \
	echo "NAME=$$NAME ACC=$$ACC" ; \
	if [ ! -f $(HIT_OUTPUTDIR)/$${NAME}_$${ACC}_positive_$(CHR).bed ] ; then \
	    echo "Missing: $(HIT_OUTPUTDIR)/$${NAME}_$${ACC}_positive_$(CHR).bed" ; \
	    ./pssm_scan --outdir $(OUTPUTDIR)/$(CHR) --genome $(GENOME) --pssm $(JASPAR) --score-mode $(SCORE_MODE) --threshold $(SCAN_THRESHOLD) --motif $$ACC --chr $(CHR) $(SCAN_CHR_FLAGS) ; \
	elif [ ! -f $(HIT_OUTPUTDIR)/$${NAME}_$${ACC}_negative_$(CHR).bed ]; then \
	    echo "Missing: sort -k 1,1 -k2,2n" ; \
	    ./pssm_scan --outdir $(OUTPUTDIR)/$(CHR) --genome $(GENOME) --pssm $(JASPAR) --score-mode $(SCORE_MODE) --threshold $(SCAN_THRESHOLD) --motif $$ACC --chr $(CHR) $(SCAN_CHR_FLAGS) ; \
	fi

$(HIT_OUTPUTDIR)/%_bidirect_$(CHR).bed.gz: $(HIT_OUTPUTDIR)/%_negative_$(CHR).bed.gz $(HIT_OUTPUTDIR)/%_positive_$(CHR).bed.gz
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

$(DISTRIBUTIONDIR)/%_score_distribution_$(SCORE_MODE)_bins_$(DISTRIBUTION_BIN_WIDTH)_positive_$(DISTRIBUTION_CHR)$(DISTRIBUTION_PSEUDOCOUNT_LABEL).tsv: pssm_scan $(JASPAR) $(GENOME) $(GENOME_INDEX)
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

test: pssm_scan $(JASPAR) Homo_sapiens.GRCh38.dna.primary_assembly_top500000.fasta Homo_sapiens.GRCh38.dna.primary_assembly_top500000.fasta.fai Homo_sapiens.GRCh38.dna.primary_assembly_bottom500000.fasta
	#./pssm_scan --genome Homo_sapiens.GRCh38.dna.primary_assembly_top500000.fasta -l -5400 --verbose -m $(TP73_MOTIF_ID) --chr $(CHR) --from 100000 --to 103000 --help
	./pssm_scan --pssm $(JASPAR) --genome Homo_sapiens.GRCh38.dna.primary_assembly_top500000.fasta --score-mode $(SCORE_MODE) -l 0 --verbose -m MA0059.1 --chr $(CHR) --from 100001 --to 103001 -s
	#./pssm_scan --genome Homo_sapiens.GRCh38.dna.primary_assembly_top500000.fasta -l 0 --verbose -m MA0019.1 --chr $(CHR) --from 100000 --to 130000
	#./pssm_scan --genome Homo_sapiens.GRCh38.dna.primary_assembly_top500000.fasta -l -500 --verbose -m MA1001.3 --chr $(CHR) --from 100000 --to 130000
	#./pssm_scan --genome Homo_sapiens.GRCh38.dna.primary_assembly_top500000.fasta -l -500 --verbose -m $(TP73_MOTIF_ID) --chr $(CHR) --from 100001 --to 103001 -s
	#./pssm_scan --genome Homo_sapiens.GRCh38.dna.primary_assembly_bottom500000.fasta -l -500 --verbose -o output_bottom --chr 44 --from 100000 --to 103000
	#./pssm_scan --genome Homo_sapiens.GRCh38.dna.primary_assembly_bottom500000.fasta -l -500 --verbose -o output_bottom --from 100000 --to 103000

test_reference_tp73_promoter_chr1: pssm_scan $(JASPAR) $(GENOME) $(GENOME_INDEX)
	@if [ "$(TEST_REFERENCE_THRESHOLD)" != "0" ]; then \
	    echo "E: test_reference_tp73_promoter_chr1 is calibrated for TEST_REFERENCE_THRESHOLD=0." ; \
	    exit 1 ; \
	fi
	mkdir -p $(TEST_REFERENCE_OUTPUTDIR)
	for motif in $(TEST_REFERENCE_MOTIFS); do \
	    ./pssm_scan --outdir $(TEST_REFERENCE_OUTPUTDIR) --genome $(GENOME) --pssm $(JASPAR) --score-mode $(SCORE_MODE) --motif $$motif --chr $(TEST_REFERENCE_CHR) --from $(TEST_REFERENCE_FROM) --to $(TEST_REFERENCE_TO) --threshold $(TEST_REFERENCE_THRESHOLD) $(SCAN_CHR_FLAGS) ; \
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

.PHONY: test check check-r check-duckdb check_synthetic_dense check_synthetic_sparse synthetic_dense_example check_cutandrun_containment synthetic_cutandrun_example all $(OUTPUTDIR)/$(CHR) jaspar genome genomegz genome_index scan_chr_all_motifs dry_run_chr1_patz1_tp73 inspect_chr1_patz1_tp73 genome_testdata count datatables files_cutandrun_clean TP73_datatable test_reference_tp73_promoter_chr1 score_distributions_chr1
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
$(OUTPUTDIR)/$(CHR): $(addprefix $(HIT_OUTPUTDIR)/,$(BIDIRECT_FILES))

datatables: context
	#for chr in 2; do
	@for chr in $(shell seq 1 22) X Y; do \
	    echo "Chr $$chr" ; \
	    $(MAKE) CHR=$$chr TP73_datatable ; \
	done

TP73_datatable: $(CONTEXT_DATATABLE)

$(CONTEXT_DATATABLE): $(TP73_COMBINED_BED).gz context
	@mkdir -p "$(CONTEXT_TMPDIR)"
	@set -euo pipefail; \
	    temporary=$$(mktemp "$@.tmp.XXXXXX"); \
	    trap 'rm -f "$$temporary"' EXIT HUP INT TERM; \
	    ./context --flank "$(CONTEXT_FLANK)" \
	        --batch-size "$(CONTEXT_BATCH_SIZE)" \
	        --temp-dir "$(CONTEXT_TMPDIR)" \
	        "$<" $(HIT_OUTPUTDIR)/*bidirect*.bed.gz \
	      | gzip -9 -n -c > "$$temporary"; \
	    mv "$$temporary" "$@"; \
	    trap - EXIT HUP INT TERM


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
