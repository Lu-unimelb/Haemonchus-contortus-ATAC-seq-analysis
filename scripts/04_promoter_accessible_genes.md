# Extract promoter windows from transcription start site (TSS) BED file
# For + strand: [-1000 bp upstream, +500 bp downstream]
# For – strand: [-500 bp upstream, +1000 bp downstream]

awk 'BEGIN{OFS="\t"}
{
  if ($6 == "+") {
    start = $2 - 1000; if (start < 0) start = 0;
    end = $3 + 500;
  } else {
    start = $2 - 500; if (start < 0) start = 0;
    end = $3 + 1000;
  }
  print $1, start, end, $4, $5, $6;
}' /results/02_ATACseq_data_analysis/genome/heacon5_annotation.tss.bed > /results/04_promoter_accessible_genes/promoter_window.bed


# filter sig.value($5) >= 540(= IDR < 0.05)
awk '$5 >= 540 && $10 != -1 {
  start = $2 + $10;
  end = start + 1;
  print $1, start, end;
}' OFS="\t" /results/03_peaks_qc_and_filtering/IDR_HCA1_HCA2 > /results/04_promoter_accessible_genes/IDR_filtered_peaks_summit.bed


# Construct gene and their transcripts mapping file
awk '$3 == "transcript" && $0 ~ /gene_id/ && $0 ~ /transcript_id/ {
    match($0, /gene_id "([^"]+)"/, g);
    match($0, /transcript_id "([^"]+)"/, t);
    if (g[1] && t[1]) print g[1], t[1];
}' /data/heacon5_annotation.gtf > /results/04_promoter_accessible_genes/gene_transcript_mapping.txt

# Extract genes whose promoters overlap peaks summit
# The output is a two-column table: the first column contains promoter-accessible gene IDs(nor unique), and the second column contains their transcript IDs.
bedtools intersect -a /results/04_promoter_accessible_genes/promoter_window.bed -b /results/04_promoter_accessible_genes/IDR_filtered_peaks_summit.bed -wa -u | \
cut -f4 | \
awk 'BEGIN{OFS="\t"} 
     FNR==NR {tx2gene[$2]=$1; next} 
     {gene = tx2gene[$1] ? tx2gene[$1] : "NA"; print gene, $1}' \
     /results/04_promoter_accessible_genes/gene_transcript_mapping.txt - > /results/04_promoter_accessible_genes/promoter_accessible_genes.txt