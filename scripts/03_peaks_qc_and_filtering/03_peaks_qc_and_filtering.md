
idr --samples /results/02_ATACseq_data_analysis/bwa/merged_library/macs2/narrow_peak/HCA_REP1.mLb.clN_peaks.narrowPeak /results/02_ATACseq_data_analysis/bwa/merged_library/macs2/narrow_peak/HCA_REP2.mLb.clN_peaks.narrowPeak --input-file-type narrowPeak --rank p.value --output-file /results/03_peaks_qc_and_filtering/IDR_HCA1_HCA2 --plot

awk '$5 > 540' /results/03_peaks_qc_and_filtering/IDR_HCA1_HCA2 > IDR_filtered_peaks.narrowPeak

/scripts/03_peaks_qc_and_filtering/split_bam_and_callpeaks.sh /results/02_ATACseq_data_analysis/bwa/merged_library/HCA_REP1.mLb.clN.sorted.bam /results/03_peaks_qc_and_filtering/HCA1_pseudo

idr --samples /results/03_peaks_qc_and_filtering/HCA1_pseudo/HCA_REP1.mRp.clN.sorted_pr1_peaks.narrowPeak /results/03_peaks_qc_and_filtering/HCA1_pseudo/HCA_REP1.mRp.clN.sorted_pr2_peaks.narrowPeak --input-file-type narrowPeak --rank p.value --output-file /results/03_peaks_qc_and_filtering/IDR_HCA1 --plot

/scripts/03_peaks_qc_and_filtering/split_bam_and_callpeaks.sh /results/02_ATACseq_data_analysis/bwa/merged_library/HCA_REP2.mLb.clN.sorted.bam /results/03_peaks_qc_and_filtering/HCA2_pseudo

idr --samples /results/03_peaks_qc_and_filtering/HCA2_pseudo/HCA_REP2.mRp.clN.sorted_pr1_peaks.narrowPeak /results/03_peaks_qc_and_filtering/HCA2_pseudo/HCA_REP2.mRp.clN.sorted_pr2_peaks.narrowPeak --input-file-type narrowPeak --rank p.value --output-file /results/03_peaks_qc_and_filtering/IDR_HCA2 --plot

/scripts/03_peaks_qc_and_filtering/split_bam_and_callpeaks.sh /results/02_ATACseq_data_analysis/bwa/merged_replicates/HCA.mLb.clN.sorted.bam /results/03_peaks_qc_and_filtering/HCA_pseudo

idr --samples /results/03_peaks_qc_and_filtering/HCA_pseudo/HCA.mRp.clN.sorted_pr1_peaks.narrowPeak /results/03_peaks_qc_and_filtering/HCA_pseudo/HCA.mRp.clN.sorted_pr2_peaks.narrowPeak --input-file-type narrowPeak --rank p.value --output-file /results/03_peaks_qc_and_filtering/IDR_HCA --plot --idr-threshold 0.01