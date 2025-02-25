#!/bin/bash -e

if [ ! -d "scripts" ]; then
        echo "E: Execute the script from the root of the Jaspar folder"
        exit 1
fi

ofolder="GeneLists"
mkdir -p "$ofolder"

URL_prefix="https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp?geneSetName="
URL_postfix="&fileType=TSV"
lists="HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION HALLMARK_WNT_BETA_CATENIN_SIGNALING HALLMARK_ANGIOGENESIS KEGG_MELANOMA REACTOME_NEURONAL_SYSTEM"
for list in $lists; do
    echo "$list"
    URL="$URL_prefix$list$URL_postfix"
    echo "$URL"

    fname="$ofolder/$list.tsv"
    if [ -f "$fname" ]; then
        echo "I: $fname already exists - skipping download"
    else
        echo "I: Downloading $list to $fname"
        wget -O "$fname" "$URL"
    fi
    if [ -f "$ofolder/$list.promoter.bed" ] && [ -f "$ofolder/$list.transcript.bed" ]; then
        echo "I: $list.promoter.bed and $list.transcript.bed already exist - skipping"
        continue
    else
        genes=$(grep GENE_SYMBOLS "$fname" | cut -f2 | grep -v "^GENE_SYMBOL")
        ./gtf_file_region_retrieval -n "$genes" -f promoter > "$ofolder/$list.promoter.bed"
        ./gtf_file_region_retrieval -n "$genes" > "$ofolder/$list.transcript.bed"
    fi
done
