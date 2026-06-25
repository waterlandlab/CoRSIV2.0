library(IlluminaHumanMethylationEPICanno.ilm10b5.hg38) # or use EPIC equivalent if needed
library(IlluminaHumanMethylationEPICmanifest)         # or EPIC equivalent
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(DMRcate)
library(limma)
library(minfi)
CoRSIV_control = "CoRSIV"
case_control_stats <- read.csv("./metadata.csv")
case_control_stats <- case_control_stats[case_control_stats$schizophrenia==" TRUE",]


beta <- read.csv(paste0("./",CoRSIV_control,"_cleaned_data.csv"),row.names = 1)  # Make sure probe IDs are rownames , row.names = 1 
beta <- beta[case_control_stats$ID_REF]

beta <- as.matrix(beta)
stopifnot(all(beta >= 0 & beta <= 1))
Mval <- log2(beta / (1 - beta))
# Example metadata
group <- factor(case_control_stats$case_control)  # Adjust to your data
gender <- factor(case_control_stats$gender) 
SZ <- factor(case_control_stats$schizophrenia) 
design <- model.matrix(~ case_control,data = case_control_stats) #+gender+age+cd8t+cd4t+nk+bcell+mono
colnames(case_control_stats)
# Get annotation data
annotation <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b5.hg38)

# Match probes
annotation <- annotation[rownames(Mval), ]

# Create GenomicRatioSet object
grSet <- makeGenomicRatioSetFromMatrix(beta, 
                                       pData = data.frame(group = group),
                                       array = "IlluminaHumanMethylationEPIC",
                                       annotation = "ilm10b5.hg38")



# Fit model and compute moderated t-statistics
fit <- lmFit(Mval, design)
fit <- eBayes(fit)

# Annotate for DMRcate
myAnnotation <- cpg.annotate(object = Mval, datatype = "array", what = "M",
                             analysis.type = "differential", design = design, 
                             coef = 2, arraytype = "EPIC",pcutoff=1)

# Run DMRcate
dmrcoutput <- dmrcate(myAnnotation, lambda = 1000, C = 2,pcutoff=1)
# Extract significant DMRs
results.ranges <- extractRanges(dmrcoutput, genome = "hg38")


library(GenomicRanges)
cpg_gr <- myAnnotation@ranges
hits <- findOverlaps(results.ranges, cpg_gr)

# Group probe IDs by DMR
dmr_probes <- split(names(cpg_gr)[subjectHits(hits)], queryHits(hits))
names(dmr_probes) <- paste0("DMR_", names(dmr_probes))
mcols(results.ranges)$probes <- sapply(dmr_probes, function(x) paste(x, collapse = ";"))
# View top regions
output_results <- data.frame(results.ranges)

write.table(output_results,file = "CoRSIV_DMRcate_output.txt",quote = F,row.names = F,sep = "\t")

probe_list <- strsplit(output_results$probes, ";")
all_selected_probes <- unlist(strsplit(output_results$probes, ";"))
write.table(x = all_selected_probes,file = "selected_CoRSIV_probes.txt",sep = "\t",quote = F,col.names = F,row.names = F)

head(output_results)
bed_file <- output_results[c(1,2,3)]
write.table(bed_file,file = paste0(adult_neonate,"_",CoRSIV_control,".bed.txt"),quote = F,col.names = F, row.names = F,sep = "\t")













library(ggpubr)

# Create an empty list to store plots
plots <- list()

for (i in 1:dim(output_results)[1]) {
  beta <- read.csv(paste0("../cleaned_filtered_",CoRSIV_control,".csv"))
  
  # Extract probe IDs from output_results column 14
  probe_ids <- strsplit(output_results[i, 14][1], ";")[[1]]
  
  # Compute mean across selected probes
  selected_region <- data.frame(
    colMeans(beta[beta$ID_REF %in% probe_ids, -1])
  )
  colnames(selected_region) <- "beta_values"
  
  # Add metadata
  case_control_stats <- read.csv("./case_control_status.csv")
  selected_region$adult_neonate <- case_control_stats$adult_neonate
  selected_region$case_control <- case_control_stats$case_control
  
  # Filter for neonates
  selected_region <- subset(selected_region, adult_neonate == adult_neonate)
  
  # Create plot
  p <- ggdensity(xlim = c(0, 1),
    selected_region,
    x = "beta_values",
    add = "mean",
    rug = TRUE,
    color = "case_control",
    fill = "case_control",
    palette = c("red", "blue")
  ) + ggtitle(paste0(output_results[i,c(1,2,3)],collapse = "_"))
  
  # Save plot in list
  plots[[i]] <- p
}

# Arrange and save all plots into one PDF file
pdf(paste0(adult_neonate,"_",CoRSIV_control,"_ART_not_ART.pdf"), width = 12, height = 10)
ggarrange(
  plotlist = plots,
  ncol = 5, nrow = 4
)
dev.off()


data <- read.csv("./output/adult_controls.fdr0.1_maxdiff0.05.csv",sep = ",",header = TRUE,flush = T)
data <- data[c(1,2,3,7,11,13)]
# Load required libraries
library(dplyr)
library(stringr)

# Set working directory to the folder containing your files
setwd("./output")  # Replace with your actual path

# List the 4 specific files
files <- c(
  "adult_controls.fdr0.1_maxdiff0.05.csv",
  "adult_corsivs.fdr0.1_maxdiff0.05.csv",
  "neonate_controls.fdr0.1_maxdiff0.05.csv",
  "neonate_corsivs.fdr0.1_maxdiff0.05.csv"
)

# Read and combine all files with a new column from filename prefix
combined_df <- lapply(files, function(file) {
  data <- read.csv(file,sep = ",",header = TRUE,flush = T)
  data <- data[c(1,2,3,7,11,13)]
  data$Source <- str_extract(file, "^[^\\.]+")  # Get text before first "."
  return(data)
}) %>%
  bind_rows()

# View result
head(combined_df)

write.csv(combined_df,file = "./combined_pval_diff.csv",row.names = F)

combined_df <- combined_df[combined_df$Source %in% c("neonate_controls","neonate_corsivs"),]

# Transform p-value
combined_df$log10_pval <- -log10(combined_df$min_smoothed_fdr)

# Plot
ggplot(combined_df, aes(x = maxdiff, y = log10_pval, color = Source)) +
  geom_point(alpha = 1, size = 2.5) +
  theme_minimal(base_size = 16) +
  labs(
    title = "Volcano Plot CoRSIV and Controls",
    x = "Mean Difference",
    y = "-log10(FDR)"
  ) +geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") + 
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "red")+
  geom_vline(xintercept = -0.05, linetype = "dashed", color = "red")+
  theme(legend.title = element_blank())+
  scale_color_manual(values = c(
    "adult_controls" = "skyblue",
    "adult_corsivs" = "coral",
    "neonate_controls" = "skyblue",
    "neonate_corsivs" = "coral"
  ))


CoRSIV_data_ML <- read.csv("../cleaned_datasets/corsiv_probes_ART_nonART_data.csv")

CoRSIV_data_ML_neonate <- CoRSIV_data_ML[CoRSIV_data_ML$Donor_time=="neonate",]
CoRSIV_data_ML_adult <- CoRSIV_data_ML[CoRSIV_data_ML$Donor_time=="adult",]
CoRSIV_data_ML_neonate <- CoRSIV_data_ML_neonate[,-c(1,2)]
CoRSIV_data_ML_adult <- CoRSIV_data_ML_adult[,-c(1,2)]

data <- merge(CoRSIV_data_ML_neonate,CoRSIV_data_ML_adult,by="Donor_ID")



