
findMotifsGenome.pl /results/03_peaks_qc_and_filtering/IDR_HCA1_HCA2 /data/heacon5_genome.fna /results/08_motif_analysis/findMotifsGenome/ -size 200 -mask

# 1. ATACorrect
TOBIAS ATACorrect --bam /results/02_ATACseq_data_analysis/bwa/merged_replicates/HCA.mLb.clN.sorted.bam --genome /data/heacon5_genome.fna --peaks /results/02_ATACseq_data_analysis/bwa/merged_replicate/macs2/narrow_peak/HCA.mRp.clN_peaks.narrowPeak --outdir /results/08_motif_analysis/ATACorrect_results --cores 20

# 2. FootprintScores
TOBIAS FootprintScores --signal /results/08_motif_analysis/ATACorrect_results/HCA.mRp.clN.sorted_corrected.bw --regions /results/04_promoter_accessible_genes/heacon5_IDR_filtered_peaks.bed --output /results/08_motif_analysis/FootprintScores_results/footprints.bw --cores 30

# 3. BINDetect
TOBIAS BINDetect --motifs /results/08_motif_analysis/findMotifsGenome/<motif_name>.meme --signals /results/08_motif_analysis/FootprintScores_results/footprints.bw --genome /data/heacon5_genome.fna --peaks /results/04_promoter_accessible_genes/heacon5_IDR_filtered_peaks.bed --outdir /results/08_motif_analysis/BINDetect_results --cond_names hc --cores 30