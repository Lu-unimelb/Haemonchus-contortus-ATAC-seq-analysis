
nextflow run nf-core/atacseq -profile conda --input /scripts/02_ATACseq_data_analysis/samplesheet.csv --outdir /results/02_ATACseq_data_analysis/ --fasta /data/heacon5_genome.fa --gtf /data/heacon5_annotation.gtf.gz --save_reference --read_length 150 --mito_name MT --narrow_peak
