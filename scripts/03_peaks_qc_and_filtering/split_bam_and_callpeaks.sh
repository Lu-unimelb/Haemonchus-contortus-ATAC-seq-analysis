#!/bin/bash
# Usage: ./split_bam_and_callpeaks.sh input.bam output_folder

# Input BAM file
BAM=$1
# Output folder
OUTDIR=$2

if [ -z "$BAM" ] || [ -z "$OUTDIR" ]; then
    echo "[ERROR] Missing arguments."
    echo "Usage: $0 input.bam output_folder"
    exit 1
fi

# Extract prefix (filename without path and .bam extension)
PREFIX=$(basename "$BAM" .bam)

# Create output folder
mkdir -p "$OUTDIR"

echo "[INFO] Sorting BAM by read name..."
samtools sort -n -o "$OUTDIR/${PREFIX}_sorted.bam" "$BAM"

echo "[INFO] Extracting properly paired reads..."
samtools view -f 2 "$OUTDIR/${PREFIX}_sorted.bam" > "$OUTDIR/paired_reads.sam"

echo "[INFO] Randomly shuffling read pairs..."
awk '{printf("%s%s",$0,(NR%2==0)?"\n":"\0")}' "$OUTDIR/paired_reads.sam" | \
shuf | \
tr "\0" "\n" > "$OUTDIR/shuffled.sam"

echo "[INFO] Splitting shuffled read pairs into two pseudo-replicates..."
split -d -l $(($(wc -l < "$OUTDIR/shuffled.sam") / 2)) "$OUTDIR/shuffled.sam" "$OUTDIR/pseudo_rep_"

echo "[INFO] Extracting header from original BAM..."
samtools view -H "$BAM" > "$OUTDIR/header.sam"

echo "[INFO] Creating pseudo-replicate BAM files..."
cat "$OUTDIR/header.sam" "$OUTDIR/pseudo_rep_00" | samtools view -b -o "$OUTDIR/${PREFIX}_pr1.bam"
cat "$OUTDIR/header.sam" "$OUTDIR/pseudo_rep_01" | samtools view -b -o "$OUTDIR/${PREFIX}_pr2.bam"

# -------------------------------
# MACS2 Peak Calling
# -------------------------------
echo "[INFO] Calling peaks using MACS2..."

# Genome size (C. elegans = 267,733,863 bp)
GSIZE=267733863

# Call peaks for pseudo-replicate 1
macs2 callpeak \
  --keep-dup all --nomodel \
  --gsize $GSIZE \
  --format BAMPE \
  --name ${PREFIX}_pr1 \
  --treatment "$OUTDIR/${PREFIX}_pr1.bam" \
  --outdir "$OUTDIR"

# Call peaks for pseudo-replicate 2
macs2 callpeak \
  --keep-dup all --nomodel \
  --gsize $GSIZE \
  --format BAMPE \
  --name ${PREFIX}_pr2 \
  --treatment "$OUTDIR/${PREFIX}_pr2.bam" \
  --outdir "$OUTDIR"

echo "[DONE] Peaks generated in $OUTDIR:"
echo "    $OUTDIR/${PREFIX}_pr1_peaks.narrowPeak"
echo "    $OUTDIR/${PREFIX}_pr2_peaks.narrowPeak"

