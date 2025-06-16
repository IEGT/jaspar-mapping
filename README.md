The motivation behind this repository was to get some simple means to compare the peaks of a CUT&RUN data set with the expected binding sites of a transcription factor as described by the JASPAR database (or another source of PSSM). Data provided by web services or readily usable .BED files only provide matches above a threshold.

Data provided by this script shall then be used augment data shown in IGV, so novel binding sites determined from the CUT&RUN experiment can be matched against the expected established ones. Not all predicted TFBS will be confirmed by CUT&RUN and vice versa. It may be worthwhile to also generate new descriptions of the binding sites if the transcription factor binding site is not yet described in a database like JASPAR. This software will also identify those transcription factor binding sites in JASPAR that are co-located with the prime TFBS-motif investigated and that are more frequent in those predicted binding sites that have a confirmation by the CUT&RUN data.

To execute this program
 1. Compile binaries
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

Unconstrained, for the full genome, expect an output of ~30 GB for each of the ~2000 motifs of JASPAR. Options are provided to retrieve matches above a given threshold or for particular chromosomal regions.

For an integration of data from a CUT&RUN experiment without referal to peaks, the bedGraph data may be transformed to .bed files and subsequently be used with bedtools to combined these data with TFBS data.
 
 . copy or link "*.clipped.clean.bedGraph" files to a local directory
 . inform Makefile about paths and transform .bedGraph to .bed by `make files_cutandrun_clean`

A helper routine was created to filter GTF annotations of genome. This can be used to retrieve genomic coordinates for a set of genes of interest, much like a "quick local BioMart". To execute that program
 1. Compile gtf_file_region_retrieval, again by invocating make.
 2. Download the GTF file accompanying the FASTA file, i.e. https://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz
 3. Execute it ./gtf_file_region_retrieval

Associate all filtered matches of your transcription factor's predicted binding sites with CUT&RUN matches.

 1. Create a folder with a copy of the data you aim to use as a reference of your CUT&RUN data to a dedicated folder. This could be the .bed files or .bigGraph files showing the detected peaks or the coverage genomic regions.
 2. Adjust the variable "PATH_CUTNRUN" in the Makefile.
 3. Adjust the list of .bed files that shall be merged with the predicted p73 binding sites.
 4. Auto-transform the .bigGraph files to .bed files:
    make files_cutandrun_clean
    and inspect the created files.
 5. Map the CUT&RUN results to the predicted TFBS with 
    for i in $(seq 1 22) X Y; do echo $i; make CHR=$i TP73_MA0861.1_bidirect_$i.combined.bed.gz; done
    which internally invokes "bedtools map" with the output 0 when no map is possible.

--
  Steffen Möller, IEGT, 9/2024
