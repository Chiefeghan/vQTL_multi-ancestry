########create ld matrices

#!/bin/bash

#SBATCH --job-name=test_run
#SBATCH -n 32
#SBATCH -N 1
#SBATCH --mem=180G
#SBATCH -t 4:45:00
#SBATCH -A INOUYE-SL3-CPU
#SBATCH --no-requeue
#SBATCH -p icelake-himem
#SBATCH -e test_run.error
#SBATCH -o test_run.out



##Set base dir
WORK_DIR="/rds/user/cb2253/hpc-work/plink_file/MVMR"
LDSTORE="/home/cb2253/chief/chief_apps/ldstore_v2.0_x86_64/ldstore_v2.0_x86_64"

##Directory with z files (for LDstorev2)
TXT_DIR="/rds/user/cb2253/hpc-work/plink_file/MVMR/MVMR_EUR/MVMR_LDstore"

####Log outputs and errors
LOG_FILE="$TXT_DIR/ldstore_run.log"
echo "LDstore run started at $(date)" > "$LOG_FILE"

####Loop throug all .z files in the directory
for z_file in $TXT_DIR/*.z; do
  #Extract the chr and prot name
  base_name=$(basename "$z_file" .z)
  chr=$(echo "$base_name" | grep -o "chr[0-9]\+")

  #create the config file in TXT_DIR ( needed to run LDsrtorev2!)
  config_file="$TXT_DIR/$base_name"
  echo "z;bgen;bgi;sample;bdose;bcor;ld;n_samples" > "$config_file"
  echo "$z_file;$WORK_DIR/${chr}_new.bgen;$WORK_DIR/${chr}_new.bgen.bgi;$WORK_DIR/${chr}_new.sample;$base_name.bdose;$base_name.bcor;$base_name.ld;45486" >> "$config_file"

  ##Log
  echo "Processing $base_name" >> "$LOG_FILE"

  ##Run ldstore 1: Generate .bcor file 
  $LDSTORE --in-files "$config_file" --write-bcor --read-only-bgen >> "$LOG_FILE" 2>&1

  ##Run ldstore 2: Convert .bcor to .ld (text format) and log output
  $LDSTORE --in-files "$config_file" --bcor-to-text >> "$LOG_FILE" 2>&1

  ##Keep the intermediate files, don't delete any .bcor or .bdose files
  echo "Keeping .bcor, .bdose, and .z files for $base_name" >> "$LOG_FILE"

  #message for completion
  echo "Processing complete for $base_name, kept .ld and .bcor files." >> "$LOG_FILE"
done

######Log the end of the run
echo "LDstore run finished at $(date)" >> "$LOG_FILE"
