#!/bin/bash
# Step 2 Master Script - Updated for Clean 6-Column Intron Intersections

# Case Statement to test arguments
until [ -z "${1:-}" ]; do
    case $1 in
        --dataset)
            shift
            dataset=$1;;
        --output)
            shift
            output=$1;;
        --strict)
            shift
            strict=$1;;
        --strict_num)
            shift
            strict_num=$1;;
        --spliced_beds)
            shift
            spliced_beds=$1;;
        --exon_GFF)
            shift
            exon_GFF=$1;;
        --intron_BED)
            shift
            intron_BED=$1;;
        -* )
            echo "unrecognised argument: $1"
            exit 1;;
    esac
    shift
    if [ "$#" = "0" ]; then break; fi
done

echo "Output Prefix: $output"
echo "Strict Mode:   $strict"
echo "Strict Window: $strict_num"

## Concatenate all spliced intron bed files together and sort by start and coordinate
cat ${spliced_beds}/*spliced.introns.bed | sort -S 50% -k1,1 -k2,2n > ${output}.sorted.bed

# Strict mode - Merge cryptic tags if they are within configured windowbp of each other.
if [[ "$strict" == "yes" ]]; then
    bedtools merge -i ${output}.sorted.bed -d ${strict_num} -c 1 -o count > ${output}.strict.${strict_num}.overlap.merged.bed
    bedtools subtract -a ${output}.strict.${strict_num}.overlap.merged.bed -b ${exon_GFF} > ${output}.strict.${strict_num}.merged.bed
    output=${output}.strict.${strict_num}
elif [[ "$strict" == "no" ]]; then
    bedtools merge -i ${output}.sorted.bed -c 1 -o count > ${output}.merged.bed
fi

## Intersect against 6-column intron BED reference to drop intergenic artifacts.
## BED-4 (from -a) + BED-6 (from -b) sets the reference token name to column 8 and strand to column 10.
bedtools intersect -a ${output}.merged.bed -b ${intron_BED} -wb | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$10,$8}' | sort -k1,1V -k2,2n > ${output}.cryptics.merged.bed

## Create unique list of gene/strand tags from column 6 to iterate through
cat ${output}.cryptics.merged.bed | awk '{print $6}' | sort -V | uniq > ${output}.unique_gene_introns.tab

# Create clean temporary working directory for thread-safe output
TMP_DIR="${output}_tmp_parts"
rm -rf "$TMP_DIR"
mkdir -p "$TMP_DIR"

date
# N is the number of allowed concurrently running forks
N=8
i=0

for entry in $(cat ${output}.unique_gene_introns.tab); do 
    ((i=i%N)); ((i++==0)) && wait
    awk -v target="$entry" 'BEGIN{OFS="\t"; s=1} $6 == target {print $1,$2,$3,$4,$5,$6"__i"s; s+=1}' ${output}.cryptics.merged.bed > "${TMP_DIR}/${entry}.bed" &
done
wait
date

# Combine thread-safe temp files into annotated bed
cat "${TMP_DIR}"/*.bed > ${output}.merged.annotated.bed
rm -rf "$TMP_DIR"

## Convert into a valid DEXSeq-ready GFF file, splitting tokens on double-underscores
cat ${output}.merged.annotated.bed | awk 'BEGIN{OFS="\t"}{
    split($6, a, "__");
    gene = a[1];
    part = (a[2] != "" ? a[2] : "i1");
    print $1, "STAR_hg38", "exonic_part", $2, $3, ".", $5, ".", "transcripts \"cryptic_exon\"; exonic_part_number \"" part "\"; gene_id \"" gene "\""
}' | sort -k1,1V -k4,4n > ${output}.cryptics.gff

## Merge cryptic GFF with background exon GFF safely by checking the gene_id attribute specifically
cat ${output}.cryptics.gff ${exon_GFF} | sort -k1,1V -k4,4n -k5,5n | awk 'BEGIN{FS=OFS="\t"} $9 ~ /gene_id "ENSG/ || $9 ~ /gene_id "NOVEL/' > ${output}.total.cryptics.gff

echo "✅ Step 2 Complete. Output saved to: ${output}.total.cryptics.gff"
exit 0
