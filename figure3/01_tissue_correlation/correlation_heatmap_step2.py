import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.colors as mcolors

# Load data
df = pd.read_csv("correlation_matrix.csv")

# Prepare numeric data
df_numeric = df.set_index('region_id')

# Remove rows with more than 50% missing values
threshold = int(df_numeric.shape[1] * 1)
df_filtered = df_numeric.dropna(thresh=threshold)

df_filtered.to_csv("filtered_region_tissue_paircorr.csv", index=True, header=True)

# Fill remaining missing values with zero
df_filtered_filled = df_filtered.fillna(0)

# Define custom colormap
custom_colors = ["white","darkblue", "lightblue","yellow", "lightsalmon", "darkorange"]
custom_cmap = mcolors.LinearSegmentedColormap.from_list("custom_cmap", custom_colors)

# Create clustered heatmap
g = sns.clustermap(
    df_filtered_filled,
    metric="euclidean",
    method="average",
    cmap=custom_cmap,
    vmin=0,
    vmax=1,
    figsize=(15, 12),
    row_cluster=False
)
g.savefig("correlation_heatmap_2025_Aug14.pdf")
#plt.title("Clustered Heatmap (Custom Color Scale 0–1)", pad=100)
#plt.show()