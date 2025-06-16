The motivation behind this repository was to get some simple means to compare the peaks of a cut'n'run data set with the expected binding sites of a transcription factor as described by the JASPAR database (or another source of PSSM). Data provided by web services or readily usable .BED files only provide matches above a threshold. "Weak binding sites" are now shown, so a local peak cannot be recoginized as such.

Data provided by this script shall then be used augment data shown in IGV, so novel binding sites determined from the cut'n'run experiment can be distinguished from the expected established ones.

To execute this program
 1. Compile pssm_scan.cpp
    - Install dependencies:
      UNIX: sudo apt install libbz2-dev make
    - Compile:
      UNiX: make
 2. Download the exact same genome that your cut'n'run experiments were run against, like +
    `wget https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz`
    and
 3. unpack it by +
   `gunzip -c Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz > Homo_sapiens.GRCh38.dna.primary_assembly.fasta`
 4. Download JASPAR, like +
    `wget https://jaspar2022.genereg.net/download/data/2022/CORE/JASPAR2022_CORE_non-redundant_pfms_jaspar.txt`
 5. Execute `./pssm_scan`, expect separate outputs for every motif and the respective positive and negative strands. The option "--help" provides extra details.

Unconstrained, for the full genome, expect an output of ~30 GB per motif. Options are provided to retrieve matches above a given threshold or for particular chromosomal regions.

For an integration of data from a CUT&RUN experiment without referal to peaks, the bedGraph data may be transformed to .bed files and subsequently be used with bedtools to combined these data with TFBS data.
 
 . copy or link "*.clipped.clean.bedGraph" files to a local directory
 . inform Makefile about paths and transform .bedGraph to .bed by `make files_cutandrun_clean`

A helper routine was created to filter GTF annotations of genome. This can be used to retrieve genomic coordinates for a set of genes of interest, much like a "quick local BioMart". To execute that program
 1. Compile gtf_file_region_retrieval, again by invocating make.
 2. Download the GTF file accompanying the FASTA file, i.e. https://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz
 3. Execute it ./gtf_file_region_retrieval

--
  Steffen Möller, IEGT, 9/2024
