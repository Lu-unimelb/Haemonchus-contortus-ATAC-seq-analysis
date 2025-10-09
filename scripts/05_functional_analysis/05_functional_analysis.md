# extract the longest protein isoform of each gene
awk '/^>/ {
	if (seq) {
        len = length(seq); # Calculate the length of the previous sequence
        if (len > maxlen[gene]) {
            maxlen[gene] = len; # Updated to new maximum length
            header[gene] = name;
            seqs[gene] = seq;
        }
    }
    name = $0; # Save the new/current header row
	match($0, /V3C99_[0-9]+/, a);
    gene = a[0];  # Save gene variable
	seq = ""; 
	next;
}
{
    seq = seq $0; # Continuously splice the current sequence before reading the next header (>) line
}
END {
    for (g in seqs) {
        print header[g];
        print seqs[g];
    }
}' /data/heacon5_proteins.faa > /results/05_functional_analysis/heacon5_proteins_longest.faa


download_eggnog_data.py -d 6231 --data_dir /results/05_functional_analysis

emapper.py -i /results/05_functional_analysis/heacon5_proteins_longest.faa --data_dir /results/05_functional_analysis --dmnd_db /results/05_functional_analysis/6231.dmnd -o heacon5_gene_functional_annotation --cpu 96 --itype proteins


Rscript enrichment_and_plot.R --outdir=/results/05_functional_analysis/
