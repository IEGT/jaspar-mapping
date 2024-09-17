The motivation behind this repository was to get some simple means to compare the peaks of a cut'n'run data set with the expected binding sites of a transcription factor as described by the JASPAR database (or another source of PSSM). Data provided by web services or readily usable .BED files only provide matches above a threshold. "Weak binding sites" are now shown, so a local peak cannot be recoginized as such.

Data provided by this script shall then be used augment data shown in IGV, so novel binding sites determined from the cut'n'run experiment can be distinguished from the expected established ones.

To execute this program
 1. Compile pssm_scan.cpp (on UNIX a mere `make pssm_scan` should suffice)
 2. Download the exact same genome that your cut'n'run experiments were run against, like `wget https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz`
 3. Download JASPAR, like https://jaspar2022.genereg.net/download/data/2022/CORE/JASPAR2022_CORE_non-redundant_pfms_jaspar.txt
 4. Execute `./pssm_scan`, expect separate outputs for every motif and the respective positive and negative strands.

For the full genome, expect an output of ~22 GB per motif.

--
  Steffen Möller, IEGT, 9/2024
