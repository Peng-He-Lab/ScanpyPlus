# ScanpyPlus
A set of files to do single-cell analysis in python

To use any of these script collections, just run these two lines in your *python kernel / Jupyter notebook*:
```python
sys.path.append('/home/ubuntu/tools/ScanpyPlus')
import Scanpyplus
```
## Citation:
[He, Lim and Sun et al.](https://www.cell.com/cell/fulltext/S0092-8674(22)01415-5)
A human fetal lung cell atlas uncovers proximal-distal gradients of differentiation and key regulators of epithelial fates

## DeepTree feature selection
![image](https://user-images.githubusercontent.com/4110443/146441826-a4079e4c-c9de-4d93-9ebe-3e1c07227eb1.png)


Among the functions in *Scanpyplus*, there's also a function to do feature gene selection (*DeepTree* algorithm). It removes garbage among highly variable genes, mitigate batch effect if you remove garbage batch by batch, and increases signal-to-noise ratio of the top PCs to promote rare cell type discovery.

[Here](https://nbviewer.jupyter.org/github/Peng-He-Lab/ScanpyPlus/tree/master/DeepTree_algorithm_demo.ipynb) is a [notebook](https://github.com/Peng-He-Lab/ScanpyPlus/tree/master/DeepTree_algorithm_demo.ipynb) to use DeepTree algorithm to "de-noise" highly-variable genes and improve initial clustering. 

A *MATLAB* implementation can be found [here](https://github.com/brianpenghe/Matlab-genomics).

This algorithm can be potentially used to reduce batch effect when fearing overcorrection, especially comparing conditions or time points. Two notebooks are provided showing "soft integration" of [embryonic limb](https://nbviewer.jupyter.org/github/Peng-He-Lab/ScanpyPlus/tree/master/Soft_integration_limb.ipynb) and [pancreas](https://nbviewer.jupyter.org/github/Peng-He-Lab/ScanpyPlus/tree/master/Soft_integration_pancreas.ipynb) data.

More conventionally, `HVGbyBatch` calls `sc.pp.highly_variable_genes` per batch, stores per-batch boolean flags and an aggregate count `highly_variable_n`
```python
Scanpyplus.HVGbyBatch(adata, batch_key='sample', min_disp=0.5)
```

`HVG_cutoff(adata, range_int=10, cutoff=5000, HVG_var='highly_variable_n', fig_size=(8,6))` Plots HVG count vs. threshold and prints the minimal threshold meeting your target count. This can also be used to determine the n for DeepTree-derived genes across multiple samples.
```python
i = Scanpyplus.HVG_cutoff(adata, cutoff=5000)
```

`HVG_Venn_Upset` compares HVG sets across provided flags (e.g., per-batch).
```python
Scanpyplus.HVG_Venn_Upset(adata, ['highly_variable_A','highly_variable_B'])
```

## Doublet Cluster Labeling (DouCLing)
![unnamed](https://user-images.githubusercontent.com/4110443/146441371-e7b4bec2-9e87-4a9d-98ad-3f3401ce13ed.jpg)

There are 4 types of doublets:

![image](https://user-images.githubusercontent.com/4110443/146040113-1c1b27e6-453e-48fa-a4e8-786ff8c759ec.png)

Cross-sample doublets can usually be identified by hastags or genetic backgrounds. Theoretically, (n-1)/n of doublets can be identified as cross-sample doublets when n samples with different hashtags/genetics are pooled equally.

Heterotypic doublets can sometimes trick data scientists into thinking they are a new type of dual-feature cell type like NKT cells etc. 
Heterotypic doublets are usually identified by matching individual cells to synthetic doublets regardless of manually curated clusters. Algorithms like Scrublet can remove a substantial part of doublet cells but not all of them. The survivor doublets can still aggregate into tiny clusters picked up by the annotaters when doing subclustering. Doublets of rarer cell types are also often missed, which obscures the discoveries of new cell types and states.
To leverage the input from biologists' manual parsing and the increased sensitivity of cluster-average signatures, I introduce here an alternative approach to facilitate heterotypic doublet cluster identification. This approach scans through individual tiny clusters and look for its "Parent 2" that gives it a unique feature that's different from its sibling subclusters sharing the same "Parent 1". 
[A notebook using published PBMC data](https://nbviewer.jupyter.org/github/Peng-He-Lab/ScanpyPlus/tree/master/DOUblet_Cluster_Labeling.ipynb) is provided.

<details>
  <summary><b>Other functions in Scanpyplus:</b></summary>

### An alternative way to call doublet subclusters based on *Scrublet* and [the gastrulation paper](https://www.nature.com/articles/s41586-019-0933-9)
`Bertie(adata,Resln=1,batch_key='batch')` was written with the help from [K. Polanski](https://github.com/ktpolanski). This script aggregates *Scrublet* scores from subclusters and makes threshold cuts based on subcluster p-values. And this is done batch by batch.
```python
Scanpyplus.Bertie(adata, Resln=1, batch_key='sample')
```
A variant version `Bertie_preclustered` applies the same subcluster-enrichment logic using an existing cluster assignment (e.g., Leiden).
```python
Scanpyplus.Bertie_preclustered(adata, cluster_key='leiden', batch_key='donor')
```

### Color utilities:
You can extract the color dict of a variable from an anndata object using `ExtractColor(adata,obsKey='louvain',keytype=int)`, 

and manually edit the color dict, and then use it to update colors in `adata.uns` using `UpdateUnsColor`.

```python
Scanpyplus.UpdateUnsColor(adata, 'cell_type', {'T cell':'#1f77b4', 'B cell':'#ff7f0e'})
```

You can also force a category to render white using `MakeWhite`
```python
Scanpyplus.MakeWhite(adata, 'condition', value='background')
```

### Metadata / obs (observation) helpers:
You can plot sankey graph between two variables of an anndata object using `ScanpySankey`. 

```python
Scanpyplus.ScanpySankey(adata, 'louvain', 'cell_type')
```

`orderGroups(adata, obsKey, new_order)` → reorder category levels.
```python
Scanpyplus.orderGroups(adata, 'louvain', ['0','2','1','3'])
```

`remove_barcode_suffix` removes the suffix after the '-' in the cell (barcode) name.
```python
Scanpyplus.remove_barcode_suffix(adata)
```

`CopyMeta` copies the metadata (both obs and var) from one object to another.
```python
Scanpyplus.CopyMeta(ref, query, obs_keys=['cell_type', 'donor'])
```

`AddMeta` stores a dataframe of obs values per each cell into an object.
```python
Scanpyplus.AddMeta(adata, df)
```

`AddMetaBatch` reads a dataframe of obs values per batch into an object. This format of metadata (rows are batch names, columns are obs categories) is more common, compact and human readable that is usually stored in *Excel* spreadsheets.
```python
Scanpyplus.AddMetaBatch(adata, batch_table, batch_key='donor', on='cell_type')
```

### Gene / var (variable)  metadata:
`OrthoTranslate` translates mouse genes to human orthologs and filter out poorly conserved genes, based on ortholog table that can be derived from Biomart etc.

### Converting file types:
`file2gz` creates .gz files which is useful for creating artificial 10X files.
```python3
Scanpyplus.file2gz('matrix.mtx')
```

`Scanpy2MM` saves an *anndata* into *MatrixMarket* form.
```python3
Scanpyplus.Scanpy2MM(adata, 'mm_out', layer='counts')
```

`mtx2df` reads *MatrixMarket* files into a dataframe.
```python
Scanpyplus.mtx2df('matrix.mtx', 'features.tsv', 'barcodes.tsv')
```

### Matrix utilities:
`GetRaw` copy .raw (or layer) to .X 
```python
Scanpyplus.GetRaw(adata)
```

 `CalculateRaw` → rebuild a raw count matrix from log-transformed counts based on `n_counts` 
```python
Scanpyplus.CalculateRaw(adata, counts_key='n_counts', log1p=True)
```

`CheckGAPDH` Quickly inspects first few entries for a given gene, optionally from a sparse matrix.
```python
Scanpyplus.CheckGAPDH(adata, sparse=True, gene='GAPDH')
```

`FindSimilarGenes` Computes gene–gene correlations and returns a sorted similarity series.
```python
Scanpyplus.FindSimilarGenes(adata, genename='XIST')
```

For large matrices, cells can be `DownSample`d based on labels such as cell types.
```python
Scanpyplus.DownSample(adata, 'cell_type', n_per_group=500)
```

Sometimes `PseudoBulk` profiles are also useful to generate, whether it's the mean, median or max.
```python
Scanpyplus.PseudoBulk(adata, groupby='sample', layer='counts')
```

### Embedding utilities:
`ShiftEmbedding` creates a platter that juxtaposes subsets of the data (batches, stages etc.) to visualize side by side.
```python
Scanpyplus.ShiftEmbedding(adata, 'batch', obsm_key='X_umap')
```

`CopyEmbedding` copies the embedding of one object to another.
```python
Scanpyplus.CopyEmbedding(ref, query, obsm_key='X_umap')
```

### Plotting stacked barplots of cell-type/condition proportions:
`celltype_per_stage_plot` and `stage_per_celltype_plot` plot horizontal and vertical bar plots respectively based on two metadata variables (cell type and stage, for example).

### UMAP plotting utilities:
`Plot3DimUMAP` generates a 3D plot (by *plotly*) of the UMAP after sc.tl.umap produces the 3D coordinates.

`plot_umap_with_labels` Plots UMAP with categorical labels drawn on top; auto-adjusts text to avoid overlap and saves PNG/PDF.
```python
Scanpyplus.plot_umap_with_labels(adata, 'cell_type', 'plot.png', 'plot.pdf')
```

### Gene-level calculation and plotting:
`DEmarkers` calculates, filters and plots differentially expressed genes between two populations.

`GlobalMarkers` calculates marker genes for every cell cluster and filters them.


`GPT_annotation_genes` Ranks genes per cluster (Wilcoxon) and returns a formatted text prompt listing top markers per cluster for quick human/GPT annotation.  
```python
Scanpyplus.GPT_annotation_genes(adata, leiden_key='leiden', top_n=10)
```

`ClusterGenes` transposes a log-transformed *adata* object and performs clustering and dimension reduction to classify genes.

`Dotplot2D` plots the expression levels of a gene across two metadata categories (e.g. samples and cell types). It can be used to trace maternal contimation by plotting XIST and check a key gene's expression patterns against cell types and age etc.

### Seaborn utilities:
`snsSplitViolin` plots splitviolin plots for two populations.
```python
Scanpyplus.snsSplitViolin(adata, 'GATA3', 'cell_type', 'condition')
```

`snsCluster` plots clustermaps using an *anndata object* as input. This has been helped by Bao Zhang from [Zhang lab](https://github.com/ZhangHongbo-Lab)

`markSeaborn` marks specific genes on a *Seaborn* plot.

`extractSeabornRows` extracts the rowlabels of a *Seaborn object* and saves into a *Series*.

### Plotting Venn / UpSet diagram:
`Venn_Upset` can be used to directly plot upset plots (bar plots of each category of intersections).

### Treemap:
`Treemap(adata, output='temp', branchlist=['project','batch'], width=1000, height=700, title='title')`
Treemap of obs counts by hierarchical categories; saves PDF and CSV.
```python
Scanpyplus.Treemap(adata, output='out/treemap', branchlist=['donor','cell_type'])
```

### Label transfer:
`LogisticRegressionCellType` can learn the defining features of a variable (such as cell type) of the reference object and predict the corresponding labels of a query object. 

The saved model files and also be re-used to predict a new query object in future by `LogisticPrediction`.
</details>

<details>
  <summary><b>Functions in pandasPlus:</b></summary>

`DF2Ann` converts a dataframe into an *anndata* object.

`UpSetFromLists` plots an upset plot (barplot of Venn diagram intersections) based on lists of lists.

`show_graph_with_labels` plots an interaction graph using edges to represent connection strength (max at 1, at least 0.9 to be shown).

Dataframe values can also be used to calculate `zscore` and `Ginni` coefficients.

`cellphonedb_n_interaction_Mat` and `cellphonedb_mat_per_interaction` are useful to reformat cellphonedb outputs.
</details>
