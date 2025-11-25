#!/bin/bash

module load  R/3.6

###Directory containing the .txt files
TXT_DIR="/rds/user/cb2253/hpc-work/plink_file/MVMR/MVMR_EUR/MVMR_LDstore/test"

##Loop through all .txt files in the directory
for txt_file in $TXT_DIR/*.txt; do
  # Extract the file base name without extension
  base_name=$(basename "$txt_file" .txt)

  ##Create a corresponding .z file name
  z_file="$TXT_DIR/$base_name.z"

  #remove the 6th column only if it exists, and save as .z file
  Rscript -e "
    library(data.table)
    df1 <- fread('$txt_file')
    fwrite(df1[,-c(6)], file='$z_file', sep=' ', row.names=FALSE, quote=FALSE)
  "
done

echo "Conversion complete for all .txt files."


