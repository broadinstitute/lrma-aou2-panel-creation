#!/bin/bash

if [ "$#" -ne 6 ]; then
    echo "Usage: $0 <panel.bcf> <input.vcf.gz> <output.bcf> <num_sites> <rate_ref_missing> <rate_alt>"
    echo "Example: $0 panel.bcf cohort.vcf.gz output.bcf 100 0.0005 0.1"
    exit 1
fi

PANEL=$1
INPUT_VCF=$2
OUTPUT=$3
NUM_SITES=$4
RATE_REF=$5
RATE_ALT=$6

# How far upstream/downstream to look in the input to catch shifted indels
VIEW_WINDOW=50 

echo "[INFO] Streaming a genotype-aware sample of $NUM_SITES target sites from $OUTPUT..."
echo "[INFO] Hom-Ref/Missing Rate: $RATE_REF | Alt Rate: $RATE_ALT"

# 1. Query the OUTPUT file. [%GT,] comma-separates the cohort's genotypes.
# 2. awk checks if ANY sample in the cohort possesses an alternate allele (1).
# 3. awk strips the GT column back off so the downstream loop reads the expected 4 columns.
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT,]\n' "$OUTPUT" \
    | awk -v r_ref="$RATE_REF" -v r_alt="$RATE_ALT" '
        BEGIN { srand() }
        {
            # If the GT string contains a "1", at least one sample has the variant
            if ($5 ~ /1/) {
                rate = r_alt
            } else {
                rate = r_ref
            }
            
            # Roll the dice
            if (rand() < rate) {
                print $1 "\t" $2 "\t" $3 "\t" $4
            }
        }' \
    | head -n "$NUM_SITES" \
    | sort -V -k1,1 -k2,2n > target_sites.txt

echo "[INFO] Generating comparison..."
echo "--------------------------------------------------------"

while read -r CHROM POS REF ALT; do
    echo -e "\n=== Target: $CHROM:$POS ==="
    echo -e "Panel: \t\t$REF\t$ALT"
    
    START=$((POS - VIEW_WINDOW))
    END=$((POS + VIEW_WINDOW))
    if [ "$START" -lt 1 ]; then START=1; fi

    # Query the Input with a wide window to catch shifted variants
    # Note: \t inside the bracket cleanly spaces out multi-sample joint cohorts
    echo "Input Window ($START - $END):"
    bcftools query -r "$CHROM:$START-$END" -f '\t%POS\t%REF\t%ALT\t%END\t[%GT:%GQ:%PL\t]\n' "$INPUT_VCF" | sed 's/^/       /'
    
    # Query the output exactly at the coordinate, STRICTLY filtered for the sampled alleles
    echo "Output:"
    bcftools query -r "$CHROM:$POS" -i "REF=\"$REF\" && ALT=\"$ALT\"" -f '\t%POS\t%REF\t%ALT\t[%GT:%PL\t]\n' "$OUTPUT" | sed 's/^/       /'

done < target_sites.txt

echo -e "\n--------------------------------------------------------"
echo "[INFO] Verification complete. Cleaning up..."
rm target_sites.txt
