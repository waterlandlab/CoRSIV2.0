# Transcription Factor Overlap Analysis Methodology

To evaluate the enrichment of transcription factor (TF) binding sites within our query genomic regions, we performed an overlap analysis using the ReMap 2022 database. The methodology consisted of the following steps:

## 1. Input Data
- **Query Regions:** Genomic intervals were provided in standard BED format (located within the `query/` directory). The total number of valid regions ($N$) within each query file was quantified, excluding any empty lines.
- **TF Binding Data:** We utilized MACS2-called peak files from the ReMap 2022 database (`remap2022_all_macs2_hg38_v1_0_*.bed`), which catalog experimentally determined TF binding sites across multiple cell lines (e.g., HepG2, K562, MCF7).

## 2. Genomic Overlap Calculation
We used `bedtools intersect` (with `-wa -wb` parameters) to identify overlaps between the query regions and the ReMap 2022 TF peaks. 

For each resulting overlap:
- The specific TF name was parsed from the annotation field of the ReMap interval.
- To account for redundant overlaps, we recorded whether a unique query interval (defined by chromosome, start, and end coordinates) overlapped a given TF's peak. 
- The raw count (`query_count`) was defined as the total number of unique query intervals overlapping at least one peak for a specific TF.

## 3. Normalization
To enable fair comparison across query files with varying numbers of regions, the raw TF overlap counts were normalized against the total number of intervals in the corresponding query file. 

The normalized count metric was calculated as **overlaps per 1,000 query regions**:
$$ \text{Normalized Count} = \left( \frac{\text{query\_count}}{N} \right) \times 1000 $$
where $N$ is the total number of regions in the query BED file.

**Example Calculation (Black CoRSIVs and CTCF in HepG2 cells):**
To illustrate, consider the evaluation of CTCF binding within Black CoRSIV regions in the HepG2 cell line:
1. **Total Query Regions ($N$):** The `Black_CpGge5_Rangegt20.bed` file contains exactly 4,354 regions ($N = 4,354$).
2. **Raw Overlaps (`query_count`):** Overlapping those 4,354 regions against CTCF peaks sequenced in HepG2 cells yielded 341 unique query intervals with at least one overlap (`query_count` = 341).
3. **Normalized Count:** Integrating these values into the formula yields $\left( \frac{341}{4354} \right) \times 1000 \approx 78.32$. 

This means that for every 1,000 Black CoRSIV regions, approximately 78 contain a CTCF binding site in HepG2 cells. This normalization ensures that this metric can be accurately and fairly compared against other datasets (e.g., White CoRSIV regions) regardless of differing initial sample sizes.

## 4. Output and Visualization
The overlap counts and normalized values for every TF across all query files and cell lines were aggregated into a tabular format (`tf_overlap_normalized_counts.tsv`). 

Subsequent exploratory visualizations included:
- **Bar Plots:** Displaying the top 10 most enriched TFs (by normalized count) for each query file and cell line combination.
- **Heatmaps:** Illustrating the differential enrichment profiles of the top 50 TFs per cell line across different query files. To construct the heatmaps, normalized counts were aggregated by TF and query file, followed by a row-wise min-max scaling to represent relative enrichment on a 0 to 1 scale.
- **Scatter Plots:** Highlighting the comparative differential enrichment of the top 100 TFs (ranked by total normalized counts across White, Black, and Random regions) per cell line. In these plots, the normalized enrichment in White CoRSIVs (x-axis) is plotted simultaneously against the enrichment in Black CoRSIVs and the Random background (y-axis). A $y=x$ reference diagonal is overlaid to identify TFs with equal enrichment between the query sets. To facilitate interpretation of the strongest signals, the top 20 TFs most enriched in White CoRSIVs are explicitly labeled using text-repelling algorithms to prevent label overlap.
## 5. Manuscript Contribution

To analyze the localized enrichment of transcription factor (TF) binding at our target loci, we assessed the genomic intersection of region sets (e.g., Black and White CoRSIVs) against experimentally determined TF peaks across multiple cell lines derived from the ReMap 2022 database (MACS2, hg38). Overlaps were identified using BEDTools, quantifying the number of unique target regions intersecting at least one peak for each individual TF. To mitigate inherent sample-size bias and allow for robust cross-comparisons between differently sized region sets, the raw overlapping counts for each TF were normalized relative to the total number of regions in the query set, yielding a standardized metric of unique overlaps per 1,000 query regions. Using this pipeline, we demonstrated differential TF occupancy, enabling side-by-side evaluations across variable genomic sequence classes regardless of the overall input library size (e.g., resulting in an observed ~78.3 CTCF binding overlaps per 1,000 Black CoRSIV regions modeled in HepG2 cells).
