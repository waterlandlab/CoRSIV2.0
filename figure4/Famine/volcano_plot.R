# Load required libraries
library(dplyr)
library(stringr)

# List the 4 specific files
files <- c(
  "CoRSIV.DMRcate_output.txt",
  "Control.DMRcate_output.txt"

)
file="CoRSIV.DMRcate_output.txt"
# Read and combine all files with a new column from filename prefix
combined_df <- lapply(files, function(file) {
  data <- read.csv(file,sep = "\t",header = TRUE)
  data <- data[c(1,2,3,7,12,13,14)]
  data$Source <- str_extract(file, "^[^\\.]+")  # Get text before first "."
  return(data)
}) %>%
  bind_rows()

# View result
head(combined_df)

write.csv(combined_df,file = "./combined_pval_diff.csv",row.names = F)

# Transform p-value
combined_df$log10_pval <- -log10(combined_df$min_smoothed_fdr)

# Plot
ggplot(combined_df, aes(x = meandiff, y = log10_pval, color = Source)) +
  geom_point(alpha = 1, size = 2.5) +
  theme_minimal(base_size = 16) +
  labs(
    x = "Mean Difference",
    y = "-log10(FDR)"
  ) + geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "red") + geom_vline(xintercept = 0.05, linetype = "dashed", color = "red")+
  geom_vline(xintercept = -0.05, linetype = "dashed", color = "red") + theme(legend.title = element_blank()) + scale_color_manual(values = c(
    "Control" = "gray",
    "CoRSIV" = "darkgreen"

  )) + theme(legend.position = "none") + xlim(-0.12, 0.12)


test <- combined_df[abs(combined_df$meandiff)>0.05 & combined_df$min_smoothed_fdr<0.05,]

library(ggplot2)
ggplot(
  combined_df,
  aes(
    x = meandiff,
    y = log10_pval,
    color = Source,
    shape = Source
  )
) +
  geom_point(size = 2.5) +
  labs(
    title = "",
    x = "Mean Difference",
    y = "-log10(FDR)"
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(-0.05, 0.05), linetype = "dashed", color = "red") +
  scale_color_manual(values = c(
    "Control" = "gray",
    "CoRSIV" = "darkgreen"
  )) +
  scale_shape_manual(values = c(
    "Control" = 17,
    "CoRSIV" = 16
  )) +
  coord_cartesian(
    xlim = c(-0.12, 0.12),
    ylim = c(0, 7.5)
  ) +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18)
  )

