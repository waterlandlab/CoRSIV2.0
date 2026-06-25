library(IlluminaHumanMethylationEPICanno.ilm10b5.hg38) # or use EPIC equivalent if needed
library(IlluminaHumanMethylationEPICmanifest)         # or EPIC equivalent
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(DMRcate)
library(limma)
library(minfi)

adult_neonate="adult"
CoRSIV_control = "controls"

case_control_stats <- read.csv("./case_control_status.csv")
cellcounts <- read.csv("./Novakovic_EpiDISH_blood_fractions_RPC.csv")
case_control_stats <- merge(case_control_stats,cellcounts,by.x="ID_REF",by.y="Sample_ID")
case_control_stats <- case_control_stats[case_control_stats$adult_neonate==adult_neonate,]

beta <- read.csv(paste0("./cleaned_filtered_",CoRSIV_control,".csv"), row.names = 1)  # Make sure probe IDs are rownames 
beta <- beta[case_control_stats$ID_REF]

beta <- as.matrix(beta)
stopifnot(all(beta >= 0 & beta <= 1))
Mval <- log2(beta / (1 - beta))
# Example metadata
group <- factor(case_control_stats$case_control)  # Adjust to your data
cellcounts <- case_control_stats[, c("B","NK","CD4T","CD8T","Mono","Neutro","Eosino")]
design <- model.matrix(~ group + B + NK + CD4T + CD8T + Mono + Neutro + Eosino,data = case_control_stats)

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
                             coef = 2, arraytype = "EPICv1",pcutoff=1)

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

with(output_results, plot(meandiff, -log10(min_smoothed_fdr), pch = 20,
                          main = "Volcano Plot Adult data - Control Probes", xlab = "Mean Diff", ylab = "-log10(min_smoothed_fdr)"))
abline(h = -log10(0.05), col = "red", lty = 2)


#write.table(output_results,file = paste0("./",adult_neonate,"_",CoRSIV_control,"_dmrcoutput.csv"),quote = F,row.names = F,sep = ",")

output_results <- output_results[output_results$min_smoothed_fdr<0.1 & abs(output_results$meandiff)>0.05,]###change this for more dmrs

probe_list <- strsplit(output_results$probes, ";")
all_selected_probes <- unlist(strsplit(output_results$probes, ";"))
write.table(x = all_selected_probes,file = paste0("./",adult_neonate,"_",CoRSIV_control,"probes.txt"),sep = "\t",quote = F,col.names = F,row.names = F)

write.table(output_results[c(1,2,3)],file =paste0("./",adult_neonate,"_",CoRSIV_control,".BED"), col.names = FALSE,row.names = FALSE,sep = "\t",quote = F)
