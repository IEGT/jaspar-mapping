#!/bin/bash -e

if [ ! -d "scripts" ]; then
        echo "E: Execute the script from the root of the Jaspar folder"
        exit 1
fi

ofolder="GeneLists"
mkdir -p "$ofolder"

URL_prefix="https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp?geneSetName="
URL_postfix="&fileType=TSV"
lists="HALLMARK_ANGIOGENESIS HALLMARK_APOPTOSIS HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION HALLMARK_INFLAMMATORY_RESPONSE HALLMARK_P53_PATHWAY HALLMARK_WNT_BETA_CATENIN_SIGNALING KEGG_MELANOMA REACTOME_NEURONAL_SYSTEM Nico_Analysis_DN_20250508 Nico_Analysis_DN_20250618 Nico_Analysis_DN_Methup_20250618 Nico_Analysis_TA_20250508 Nico_Analysis_TA_20250618 Nico_Analysis_TAandDN_20250508 Nico_Analysis_TA_MethUp_20250618"

for list in $lists; do
    echo "$list"

    fname="$ofolder/$list.tsv"
    if [ -f "$fname" ]; then
        echo "I: $fname already exists - skipping download"
    else
        URL="$URL_prefix$list$URL_postfix"
        echo "$URL"
        echo "I: Downloading $list to $fname"
        wget -O "$fname" "$URL"
    fi
    if [ -f "$ofolder/$list.promoter.bed" ] && [ -f "$ofolder/$list.transcript.bed" ]&& [ -f "$ofolder/$list.utr.bed" ]; then
        echo "I: $list.promoter.bed and $list.transcript.bed already exist - skipping"
        continue
    else
        #echo "D: Grepping genes"
        if ! genes=$(grep GENE_SYMBOLS "$fname" | cut -f2 | grep -v "^GENE_SYMBOL"); then
            echo "E: Could not grep for GENE_SYMBOLS, check if a true tab is separating the genes from the row name."
            exit 1
        fi
        #echo "D: Grepping genes found $genes"
        ./gtf_file_region_retrieval -c 1 -n "$genes" -f promoter > "$ofolder/$list.promoter.bed"
        ./gtf_file_region_retrieval -c 1 -n "$genes" -f utr > "$ofolder/$list.utr.bed"
        ./gtf_file_region_retrieval -c 1 -n "$genes" > "$ofolder/$list.transcript.bed"
    fi
done
