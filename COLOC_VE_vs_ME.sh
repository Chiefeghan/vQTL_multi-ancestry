#!/bin/bash

# Load the required modules
module load plink/2.00-alpha
module load r-3.6.1-gcc-5.4.0-zrytncq

######## Set base directory for bfiles and phenotype files
bfile_base="/rds/user/cb2253/hpc-work/plink_file/chr"
pheno_base="/rds/user/cb2253/hpc-work/plink_file/vQTL_3000/"
out_base="protein_linear_EUR_test_results"
PROTEIN_CHR_FILE="protein_mvmr.txt"

##### Create dirs
mkdir -p snp_files
TEMP_DIR=$(mktemp -d)

#Initialize counter for the SNP files
counter=1

#Create temp result file in the temp dir
TEMP_RESULT_FILE="${TEMP_DIR}/temp_results.txt"
touch "$TEMP_RESULT_FILE"

# Read input 
tail -n +2 "$PROTEIN_CHR_FILE" | while read -r line; do
    olink_tag=$(echo "$line" | awk '{print $1}')
    chr=$(echo "$line" | awk '{print $2}')
    gene_start=$(echo "$line" | awk '{print $3}')

    echo "Processing: olink_tag=$olink_tag, chr=$chr, gene_start=$gene_start"
    
    #Dfine cis-boundaries
    lower_bound=$((TSS - 1e6))
    upper_bound=$((TSS + 1e6))

    echo "Boundaries: lower_bound=$lower_bound, upper_bound=$upper_bound"

    # Determine the corresponding .vqtl file directory
    PLINK_FILE_DIR="${pheno_base}${olink_tag}"
    
    # Determine the corresponding .vqtl file (the only file with .vqtl suffix in the directory)
    vqtl_file=$(find "$PLINK_FILE_DIR" -type f -name "*.vqtl" | head -n 1)
    
    if [ -z "$vqtl_file" ]; then
        echo "No .vqtl file found for $olink_tag in $PLINK_FILE_DIR"
        continue
    fi

    echo "Found .vqtl file: $vqtl_file"

    # Extract the SNPs within the defined boundaries from the .vqtl file based on the 6th column (bp)
    awk -v lb="$lower_bound" -v ub="$upper_bound" '$6 >= lb && $6 <= ub' "$vqtl_file" > "$TEMP_DIR/${olink_tag}_snps.txt"
    
    # Create SNP file for PLINK
    snp_file="snp_files/file${counter}.txt"
    awk '{print $2}' "$TEMP_DIR/${olink_tag}_snps.txt" > "$snp_file"

    # Construct the bfile path
    bfile="${bfile_base}${chr}_imputed_filtered_v3"
    
    # Run PLINK command for the current SNP
    pheno="${pheno_base}${olink_tag}.txt"
    output_file="${out_base}_snp${counter}"
    plink2 --bfile "$bfile" --glm --pheno "$pheno" --out "$output_file" --extract "$snp_file"
    
    # Check if PLINK output file exists
    plink_output_file=$(ls ${output_file}*.glm.linear | head -n 1)
    if [ ! -f "$plink_output_file" ]; then
        echo "PLINK output file not found for $olink_tag. Skipping..."
        continue
    fi

    echo "PLINK output file found: $plink_output_file"
    
    # Extract BETA and SE from PLINK results
    awk 'NR > 1 {print $3, $9, $10}' "$plink_output_file" > "${TEMP_DIR}/${olink_tag}_plink_results.txt"
    
    # Run the R script to perform colocalization analysis
    Rscript -e "
    library(data.table)
    library(dplyr)
    library(coloc)

    # Load the input files
    p1 <- fread('$TEMP_DIR/${olink_tag}_snps.txt', col.names = c('chr', 'SNP', 'A1', 'A2', 'freq', 'bp', 'F-statistic', 'df1', 'df2', 'beta', 'se', 'P', 'NMISS'))
    p2 <- fread('${TEMP_DIR}/${olink_tag}_plink_results.txt', col.names = c('SNP', 'BETA', 'SE'))

    # Subset p1 based on lower and upper bounds
    p1 <- p1 %>% filter(bp >= $lower_bound & bp <= $upper_bound)

    # Select relevant columns for p1
    p1 <- p1 %>% select(SNP, beta, se, freq)
    colnames(p1) <- c('SNP', 'BETA', 'SE', 'FREQ')

    # Select columns for p2 and rename
    p2 <- p2 %>% select(SNP, BETA, SE)

    # Print the heads of p1 and p2 for debugging
    cat('p1 head:\n')
    print(head(p1))
    cat('p2 head:\n')
    print(head(p2))

    # Ensure columns are numeric where necessary
    p1\$varbeta <- as.numeric(p1\$SE)^2
    p2\$varbeta <- as.numeric(p2\$SE)^2
    p1\$FREQ <- as.numeric(p1\$FREQ)
    
    # Filter out invalid MAF values
    p1 <- p1[!is.na(p1\$FREQ) & p1\$FREQ > 0 & p1\$FREQ < 1]

    # Ensure both datasets have the same number of rows
    common_snps <- intersect(p1\$SNP, p2\$SNP)
    p1 <- p1 %>% filter(SNP %in% common_snps)
    p2 <- p2 %>% filter(SNP %in% common_snps)

    # Print the cleaned p1 and p2 for debugging
    cat('Cleaned p1 head:\n')
    print(head(p1))
    cat('Cleaned p2 head:\n')
    print(head(p2))

    # Perform colocalization analysis
    coloc_result2 <- coloc.abf(
        dataset1 = list(snp = p1\$SNP, beta = p1\$BETA, varbeta = p1\$varbeta, N = 45486, type = 'quant', MAF = p1\$FREQ),
        dataset2 = list(snp = p2\$SNP, beta = p2\$BETA, varbeta = p2\$varbeta, N = 45486, type = 'quant', MAF = p1\$FREQ)
    )

    # Extract the results
    nsnps <- coloc_result2\$summary[1]
    h3pp <- coloc_result2\$summary[5]
    h4pp <- coloc_result2\$summary[6]
  
    # Append the result to the temporary result file
    cat('$olink_tag $chr $gene_start', nsnps, h3pp, h4pp, '\n', file = '$TEMP_RESULT_FILE', append = TRUE)
    "
    
    # Remove temporary files for the current loop iteration
    rm "$TEMP_DIR/${olink_tag}_snps.txt" "$snp_file" "$plink_output_file"
    
    # Increment the counter
    counter=$((counter + 1))
done

# Merge colocalization results into the final output file
{
  echo "olink_tag chr gene_start nsnps H3PP_coloc H4PP_coloc"
  cat "$TEMP_RESULT_FILE"
} > "${PROTEIN_CHR_FILE%.txt}_with_results.txt"

echo "Colocalization results merged into ${PROTEIN_CHR_FILE%.txt}_with_results.txt"

# Remove the temp files and dir
rm "$TEMP_RESULT_FILE"
rm -r "$TEMP_DIR"
