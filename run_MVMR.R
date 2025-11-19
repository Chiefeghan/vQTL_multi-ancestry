library(data.table)
library(MendelianRandomization)

###paths
#/Users/cb2253/OneDrive - University of Cambridge/vQTL_analysis/MVMR_EUR/summary_data_main_vqtl/MVMV_EUR_TXT/test/
data_dir <- "/rds/user/cb2253/hpc-work/plink_file/MVMR/MVMR_EUR/MVMR_LDstore//test/"
#/Users/cb2253/OneDrive - University of Cambridge/vQTL_analysis/MVMR_EUR/summary_data_main_vqtl/MVMV_EUR_TXT/LDmatrices/first_set
###list of df_* files in the working directory
df_files <- list.files(data_dir, pattern = "^df_.*\\.txt$", full.names = TRUE)

#Initialize an empty list to store the results for all proteins
all_results <- list()

###Loop through each df_* file
for (df_file in df_files) {
  
  ##extract the protein name from the df_* file name
  file_name <- basename(df_file)
  protein_name <- sub("^df_(.*)\\.txt$", "\\1", file_name)
  
  #find any .ld file with a matching protein name and chromosome prefix (chr1_ to chr22_)
  ld_file_pattern <- paste0("chr[1-9][0-9]?_", protein_name, "\\.ld$")
  ld_file <- list.files(data_dir, pattern = ld_file_pattern, full.names = TRUE)
  
  # skip to the next protein if no corresponding .ld file is found,
  if (length(ld_file) == 0) {
    message(paste("LD file not found for:", protein_name))
    next
  }
  
  ####Read the df_ file
  df1 <- fread(df_file)
  
  #Ensure the .ld file exists and can be read
  ld_file_path <- ld_file[1]
  if (file.exists(ld_file_path)) {
    r.mat <- fread(ld_file_path)
    r.mat <- as.matrix(r.mat)
  } else {
    message(paste("LD file path is invalid:", ld_file_path))
    next
  }
  
  #prep input data for MR analysis
  input_dat <- mr_mvinput(
    bx = cbind(df1$beta.exposure, df1$beta.exposure.vqtl),
    bxse = cbind(df1$se.exposure, df1$se.exposure.vqtl),
    by = df1$beta.outcome,
    byse = df1$se.outcome,
    correlation = r.mat
  )
  
  # Define sample sizes
  nx.main <- 46485
  nx.vqtl <- 46485
  nx <- c(nx.main, nx.vqtl)
  ny <- 296525
  
  #Perform MR analysis prop variance 99.9
  res <- mr_mvpcgmm(input_dat, nx = nx, ny = ny, thres = 0.999)
  
  # Capture the results
  result <- data.frame(
    Exposure1 = res@Exposure[1],
    Exposure2 = res@Exposure[2],
    Estimate1 = res@Estimate[1],
    Estimate2 = res@Estimate[2],
    StdError1 = res@StdError[1],
    StdError2 = res@StdError[2],
    CILower1 = res@CILower[1],
    CILower2 = res@CILower[2],
    CIUpper1 = res@CIUpper[1],
    CIUpper2 = res@CIUpper[2],
    Pvalue1 = res@Pvalue[1],
    Pvalue2 = res@Pvalue[2],
    CondFstat1 = res@CondFstat[1],
    CondFstat1 = res@CondFstat[2],
    Overdispersion1 = res@Overdispersion[1],
    HeterStat1 = res@Heter.Stat[1],
    PCs1 = res@PCs[1],
    Protein = protein_name
  )
  
  # Append the result to the list
  all_results[[protein_name]] <- result
}

# Combine all results into a single data frame
final_results <- do.call(rbind, all_results)

#Write the final result to a file
write.table(final_results, file = "combined_results.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)

