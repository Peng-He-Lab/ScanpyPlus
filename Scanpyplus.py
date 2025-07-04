import gc
#import scrublet as scr
import scipy.io
from scipy import sparse
import matplotlib.pyplot as plt
from matplotlib import rcParams
import seaborn as sns
import numpy as np
import random
import sys
from sklearn.manifold import TSNE
from sklearn.preprocessing import scale
from sklearn.decomposition import TruncatedSVD
from sklearn.cluster import SpectralClustering
from sklearn.model_selection import cross_val_score, StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.feature_extraction.text import TfidfTransformer
import joblib

import numpy as np
import pandas as pd
import scanpy as sc
import anndata
import os
from scipy import sparse
from scipy import cluster
from glob import iglob
import gzip



def GPT_annotation_genes(adata, leiden_key="leiden", top_n=10):
    
    sc.tl.rank_genes_groups(adata, groupby=leiden_key, method="wilcoxon")

    marker_df = pd.DataFrame(
        {cluster: adata.uns["rank_genes_groups"]["names"][cluster][:top_n]
         for cluster in adata.uns["rank_genes_groups"]["names"].dtype.names}
    )

    output_lines = [
        "Identify cell types of human periodontal ligament cells using the following markers separately for each row.",
        "Only provide the cell type name. Do not show numbers before the name.",
        "Some can be a mixture of multiple cell types."
    ]

    for cluster in marker_df.columns:
        genes = ",".join(marker_df[cluster].tolist())
        output_lines.append(f"{cluster}:{genes}")

    formatted_output = "\n".join(output_lines)
    return formatted_output




def plot_umap_with_labels(
    adata,
    color_column,
    output_png_path,
    output_pdf_path,
    figsize=(12, 10),
    point_size=5,
    label_fontsize=8,
    label_weight='bold',
    arrow_color='red',
    dpi=300,
    palette=None
):
    """
    Plot a UMAP scatter plot with labels and save it as PNG and PDF files.

    Parameters:
        adata: AnnData object containing UMAP coordinates and grouping information.
        color_column: str, column in adata.obs used for coloring points.
        output_png_path: str, path to save the PNG file.
        output_pdf_path: str, path to save the PDF file.
        figsize: tuple, size of the figure (width, height).
        point_size: int, size of the scatter plot points.
        label_fontsize: int, font size of the labels.
        label_weight: str, font weight of the labels (e.g., 'bold').
        arrow_color: str, color of the arrows used to adjust text labels.
        dpi: int, resolution of the saved images.
        palette: list, optional, color palette for the scatter plot.

    Returns:
        None
    """
    # Import necessary packages inside the function
    import matplotlib.pyplot as plt
    from adjustText import adjust_text
    import matplotlib.colors
    import scanpy as sc

    # Create a figure and axis object
    fig, ax = plt.subplots(figsize=figsize)

    # Plot UMAP scatter plot
    sc.pl.umap(
        adata,
        color=[color_column],
        ax=ax,
        size=point_size,
        legend_loc=None,  # Disable default legend
        palette=palette if palette else list(matplotlib.colors.CSS4_COLORS.values()),
        show=False
    )

    # Add custom text labels
    texts = []
    for label in adata.obs[color_column].unique():
        subset = adata.obs[adata.obs[color_column] == label]
        idx = adata.obs.index.get_loc(subset.index[0])  # Get the integer index of the first occurrence
        x = adata.obsm['X_umap'][idx, 0]
        y = adata.obsm['X_umap'][idx, 1]
        texts.append(ax.text(x, y, label, fontsize=label_fontsize, weight=label_weight))

    # Adjust text to avoid overlaps
    adjust_text(texts, arrowprops=dict(arrowstyle='->', color=arrow_color))

    # Adjust figure margins
    plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)

    # Save the plot
    if output_png_path:
        plt.savefig(output_png_path, format='png', dpi=dpi, bbox_inches='tight')
    if output_pdf_path:
        plt.savefig(output_pdf_path, format='pdf', dpi=dpi, bbox_inches='tight')

    # Show the plot
    plt.show()
# Example usage
# from anndata import AnnData
# import scanpy as sc
#
# # Assume `adata_highQ_epi` is a preloaded AnnData object with UMAP coordinates
# adata_highQ_epi = AnnData()  # Replace this with actual AnnData loading
# adata_highQ_epi.obsm['X_umap'] = [[1, 2], [3, 4], [5, 6]]  # Mock UMAP coordinates
# adata_highQ_epi.obs['lev3_celltype'] = ['Type1', 'Type2', 'Type3']  # Mock cell types
#
# # Call the function with detailed parameters
# plot_umap_with_labels(
#     adata=adata_highQ_epi,
#     color_column='lev3_celltype',
#     output_png_path='/path/to/Epi_umap_plot.png',
#     output_pdf_path='/path/to/Epi_umap_plot.pdf',
#     figsize=(10, 8),
#     point_size=10,
#     label_fontsize=12,
#     label_weight='bold',
#     arrow_color='blue',
#     dpi=300
# )

def ExtractColor(adata,obsKey='louvain',keytype=int):
#   labels=sorted(adata.obs[obsKey].unique().to_list(),key=keytype)
    labels=adata.obs[obsKey].cat.categories
    colors=adata.uns[obsKey+'_colors']
    return dict(zip(labels,colors))

def UpdateUnsColor(adata,ColorDict,obsKey='louvain'):
    #ColorDict is like {'Secretory 3 C1': '#c0c0c5','Secretory 4 C1': '#87ceee'}
    ColorUns=ExtractColor(adata,obsKey,keytype=str)
    ColorUns.update(ColorDict)
    adata.uns[obsKey+'_colors']=list(ColorUns.values())
    return adata

def Plot3DimUMAP(adata,obsKey='leiden',obsmKey='X_umap'):
    #Make sure adata.obsm['X_umap'] contains three columns
    #The obsKey must points to a str/categorical variable
    import plotly.express as px
    ThreeDdata=pd.DataFrame(adata.obsm[obsmKey],index=adata.obs_names,columns=['x','y','z'])
    ThreeDdata[obsKey]=adata.obs[obsKey]
    fig=px.scatter_3d(ThreeDdata,x='x',y='y',z='z',color=obsKey,opacity=1,
                 color_discrete_map=ExtractColor(adata,obsKey,str))
    fig.update_traces(marker=dict(size=2, 
                              line=dict(width=0,
                                        color='black')))
    return fig
    
def ScanpySankey(adata,var1,var2,aspect=20,fontsize=12, figureName="cell type",
                 leftLabels=['Default'], rightLabels=['Default']):
    
    from pySankey.sankey import sankey #correct import command
    
    df = adata.obs[[var1, var2]].dropna()
    left = df[var1]
    right = df[var2]

    colordict={**ExtractColor(adata,var1,str),
               **ExtractColor(adata,var2,str)}  
    #Ensure type conforming
    if leftLabels == ['Default']:
        leftLabels = sorted(df[var1].unique().tolist())
    else:
        leftLabels = sorted(df[var1].dtype.type(x) for x in leftLabels)
    
    if rightLabels == ['Default']:
        rightLabels = sorted(df[var2].unique().tolist())
    else:
        rightLabels = sorted(df[var2].dtype.type(x) for x in rightLabels)
  
    return sankey(left,right,aspect=aspect,fontsize=fontsize,figure_name=figureName) # removed contradictory argument for left right/label which caused pysankey function error due to reindexing

def iRODS_stats_starsolo(samples):
    #samples should be a list of library IDs
    qc = pd.DataFrame(0, index=samples, columns=['n_cells', 'median_n_counts'])
    for sample in samples:
        #download and import data
        os.system('iget -Kr /archive/HCA/10X/'+sample+'/starsolo/counts/Gene/cr3')
        adata = sc.read_10x_mtx('cr3')
        #this gets .obs['n_counts'] computed
        sc.pp.filter_cells(adata, min_counts=1)
        #compute the qc metrics and store them
        qc.loc[sample, 'n_cells'] = adata.shape[0]
        qc.loc[sample, 'median_n_counts'] = np.median(adata.obs['n_counts']).astype(float)
        #delete downloaded count matrix
        os.system('rm -r cr3')
    return qc

def orderGroups(adata,groupby='leiden'):
    #this returns a list of group names
    sc.tl.dendrogram(adata,groupby=groupby)
    return adata.uns[f'dendrogram_'+groupby]['dendrogram_info']['ivl']

def MakeWhite(adata,obsKey='louvain',whiteCat='nan',type=str):
    temp=ExtractColor(adata,obsKey,type)
    temp[whiteCat]='#FFFFFF'
    UpdateUnsColor(adata,temp,obsKey)
    return adata

def GetRaw(adata_all):
    adata=anndata.AnnData(X=adata_all.raw.X,obs=adata_all.obs,var=adata_all.raw.var,\
obsm=adata_all.obsm,uns=adata_all.uns,obsp=adata_all.obsp)
    adata.raw=adata
    return adata

def CalculateRaw(adata,scaling_factor=10000):
    #update by Polanski in Feb 2022
    #The object must contain a log-transformed matrix
    #This function returns an integer-count object
    #The normalization constant is assumed to be 10000
    #return anndata.AnnData(X=sparse.csr_matrix(np.rint(np.array(np.expm1(adata.X).todense().transpose())*(adata.obs['n_counts'].values).transpose() / scaling_factor).transpose()),\
    #              obs=adata.obs,var=adata.var,obsm=adata.obsm,varm=adata.varm)
    adata.X=adata.X.tocsr() #this step makes sure the datamatrix is in csr not csc
    X = np.expm1(adata.X)
    scaling_vector = adata.obs['n_counts'].values / scaling_factor
    #.indptr[i]:.indptr[i+1] provides the .data coordinates where the i'th row of the data resides in CSR
    #which happens to be a cell, which happens to be what we have a unique entry in scaling_vector for
    for i in np.arange(X.shape[0]):
        X.data[X.indptr[i]:X.indptr[i+1]] = X.data[X.indptr[i]:X.indptr[i+1]] * scaling_vector[i]
    return anndata.AnnData(X=np.rint(X),obs=adata.obs,var=adata.var,obsm=adata.obsm,varm=adata.varm)


def CalculateRawAuto(adata):
    X = np.expm1(adata.X)
    #.indptr[i]:.indptr[i+1] provides the .data coordinates where the i'th row of the data resides in CSR
    #which happens to be a cell, which happens to be what we need to reverse
    for i in np.arange(X.shape[0]):
        #the object is cursed, locate lowest count for each cell. treat that as 1
        #divide other counts by it. don't round for post-fact checks
        norm_one = np.min(X.data[X.indptr[i]:X.indptr[i+1]])
        X.data[X.indptr[i]:X.indptr[i+1]] = X.data[X.indptr[i]:X.indptr[i+1]] / norm_one
    #originally this had X=np.rint(X) but we actually want the full value space here
    return anndata.AnnData(X=X,obs=adata.obs,var=adata.var,obsm=adata.obsm,varm=adata.varm)

def CheckGAPDH(adata,sparse=True,gene='GAPDH'):
    if sparse==True:
        return adata[:,gene].X[0:5].todense()
    else:
        return adata[:,gene].X[0:5]

def FindSimilarGenes(adata,genename='GAPDH'):
    #This function finds the most correlated genes for a given gene
    temp=adata.to_df()
    corr_temp=np.corrcoef(temp,rowvar=False)
    corr_temp_series=pd.Series(corr_temp[:,temp.columns.get_loc(genename)],
index=temp.columns)
    return corr_temp_series.sort_values(ascending=False)

def OrthoTranslate(adata,\
oTable='~/refseq/Mouse-Human_orthologs_only.csv'):
    adata.var_names_make_unique(join='-')
    OrthologTable = pd.read_csv(oTable).dropna()
    MouseGenes=OrthologTable.loc[:,'Gene name'].drop_duplicates(keep=False)
    HumanGenes=OrthologTable.loc[:,'Human gene name'].drop_duplicates(keep=False)
    FilteredTable=OrthologTable.loc[((OrthologTable.loc[:,'Gene name'].isin(MouseGenes)) &\
                            (OrthologTable.loc[:,'Human gene name'].isin(HumanGenes))),:]
    bdata=adata[:,adata.var_names.isin(FilteredTable.loc[:,'Gene name'])]
    FilteredTable.set_index('Gene name',inplace=True,drop=False)
    bdata.var_names=FilteredTable.loc[bdata.var_names,'Human gene name']
    return bdata

def remove_barcode_suffix(adata):
    bdata=adata.copy()
    bdata.obs_names=pd.Index([i[0] for i in bdata.obs_names.str.split('-',expand=True)])
    return bdata

def file2gz(file,delete_original=True):
    with open(file,'rb') as src, gzip.open(file+'.gz','wb') as dst:
        dst.writelines(src)
    if delete_original==True:
        os.remove(file)

def Scanpy2MM(adata,prefix='temp',write2Dobsm=['No']):
    #Scanpy2MM(adata,"./")
    #please make sure the object contains raw counts (using our CalculateRaw function)
    adata.var['feature_types']='Gene Expression'
    scipy.io.mmwrite(prefix+'matrix.mtx',adata.X.transpose(),field='integer')
    if 'gene_ids' not in adata.var.columns.unique():
        adata.var['gene_ids']=adata.var_names
    adata.var[['gene_ids','feature_types']].reset_index().set_index(keys='gene_ids').to_csv(prefix+"features.tsv", \
            sep = "\t", index= True,header=False)
    if 'No' in write2Dobsm:
        print('No embeddings written')
    else:
        for basis in write2Dobsm: #save embeddings in obs
            adata.obs[basis+'_x']=adata.obsm[basis][:,0]
            adata.obs[basis+'_y']=adata.obsm[basis][:,1]
    adata.obs.to_csv(prefix+"barcodes.tsv", sep = "\t", columns=[],header= False)
    adata.obs.to_csv(prefix+"metadata.tsv", sep = "\t", index= True)
    if 'No' in write2Dobsm:
        print('No embeddings written')
    else:
        for basis in write2Dobsm: #delete obsm->obs
            del adata.obs[basis+'_x']
            del adata.obs[basis+'_y']
    file2gz(prefix+"matrix.mtx")
    file2gz(prefix+"barcodes.tsv")
    ##file2gz(prefix+"metadata.tsv")
    file2gz(prefix+"features.tsv")

def ShiftEmbedding(adata,domain_key='batch',embedding='X_umap',nrows=3,alpha=0.9):
    from sklearn import preprocessing
    scaler = preprocessing.MinMaxScaler()
    adata.obs[embedding+'0']=adata.obsm[embedding][:,0]
    adata.obs[embedding+'1']=adata.obsm[embedding][:,1]
    X=adata.obs[embedding+'0']
    Y=adata.obs[embedding+'1']
    batch_categories=adata.obs[domain_key].unique()
    for i in list(range(len(batch_categories))):
        temp=adata[adata.obs[domain_key]==batch_categories[i]].obsm[embedding]
        scaler.fit(temp)
        X.loc[adata.obs[domain_key]==batch_categories[i]]=(scaler.transform(temp)*alpha+ [int(i/nrows),i%nrows])[:,0]
        Y.loc[adata.obs[domain_key]==batch_categories[i]]=(scaler.transform(temp)*alpha+ [i/nrows,i%nrows])[:,1]
    adata.obsm[embedding]=np.vstack((X.values,Y.values)).T
    del adata.obs[embedding+'0']
    del adata.obs[embedding+'1']
    return adata

def CopyEmbedding(aFrom,aTo,embedding='X_umap'):
    aFrom.obs['temp0']=aFrom.obsm[embedding][:,0]
    aFrom.obs['temp1']=aFrom.obsm[embedding][:,1]
    aTo.obs['temp0']=0
    aTo.obs['temp1']=0
    aTo.obs.loc[aFrom.obs_names,'temp0']=aFrom.obs['temp0']
    aTo.obs.loc[aFrom.obs_names,'temp1']=aFrom.obs['temp1']
    aTo.obsm[embedding]=np.vstack((aTo.obs['temp0'],aTo.obs['temp1'])).T
    del aFrom.obs['temp0']
    del aFrom.obs['temp1']
    del aTo.obs['temp0']
    del aTo.obs['temp1']
    return aTo

def CopyMeta(aFro,aTo,overwrite=False):
    #This function copies the metadata of one object to another
    aFrom=aFro[aFro.obs_names.isin(aTo.obs_names)][:,
               aFro.var_names.isin(aTo.var_names)]
    aFrom
    if overwrite==True:
        obs_items=aFrom.obs.columns
        var_items=aFrom.var.columns
    else:
        obs_items=aFrom.obs.columns[~aFrom.obs.columns.isin(aTo.obs.columns)]
        var_items=aFrom.var.columns[~aFrom.var.columns.isin(aTo.var.columns)]
    aTo.obs[obs_items]=np.nan
    aTo.var[var_items]=np.nan
    aTo.obs.loc[aFrom.obs_names,obs_items]=aFrom.obs.loc[:,obs_items]
    aTo.var.loc[aFrom.var_names,var_items]=aFrom.var.loc[:,var_items]
    return aTo

def AddMeta(adata,meta):
    meta_df=meta.loc[meta.index.isin(adata.obs_names),:]
    meta_df=meta_df.loc[meta_df.index.drop_duplicates(keep=False),:]
    temp=adata.copy()
#    temp.obs=temp.obs.combine_first(meta_df) #it has barcode sliding problem!
    for i in meta_df.columns:
        print("copying "+i+"\n")
        temp.obs[i]=np.nan
        temp.obs.loc[meta_df.index,i]=meta_df.loc[:,i]
    return temp

def AddMetaBatch(adata,meta_compact,batch_key='batch'):
    #import your csv file into a df. The index should be batch IDs
    temp=adata.copy()
    for i in meta_compact.columns:
        temp.obs[i]=temp.obs[batch_key].replace(to_replace=meta_compact.loc[:,i].to_dict())
    return temp

def ExtractMetaBatch(adata,batch_key='batch'):
    #return a dataframe of the most frequent value for each variable per batch key
    #This can be regarded as the reverse of AddMetaBatch except for numeric variables
    return adata.obs.groupby(batch_key).agg(pd.Series.mode)

def celltype_per_stage_plot(adata,celltypekey='louvain',stagekey='batch',plotlabel=True,\
    celltypelist=['default'],stagelist=['default'],celltypekeytype=int,stagekeytype=str,
    fontsize='x-small',yfontsize='x-small',legend_pos=(1,0.5),savefig=None):
    # this is a function for horizonal bar plots
    if 'default' in celltypelist:
        celltypelist = sorted(adata.obs[celltypekey].unique().tolist(),key=celltypekeytype)
    if 'default' in stagelist:
        stagelist = sorted(adata.obs[stagekey].unique().tolist(),key=stagekeytype)
    celltypelist=[i for i in celltypelist if i in adata.obs[celltypekey].unique()]
    stagelist=[i for i in stagelist if i in adata.obs[stagekey].unique()]
    colors=ExtractColor(adata,celltypekey,keytype=str)
    count_array=np.array(pd.crosstab(adata.obs[celltypekey],adata.obs[stagekey]).loc[celltypelist,stagelist])
    count_ratio_array=count_array / np.sum(count_array,axis=0)
    for i in range(len(celltypelist)):
        plt.barh(stagelist[::-1],count_ratio_array[i,::-1],
            left=np.sum(count_ratio_array[0:i,::-1],axis=0),color=colors[celltypelist[i]],label=celltypelist[i])
        plt.yticks(fontsize=yfontsize)
    plt.grid(visible=False)
    if plotlabel:
        plt.legend(celltypelist,fontsize=fontsize,bbox_to_anchor=legend_pos)
    if savefig is not None:
        plt.savefig(savefig+'.pdf',bbox_inches='tight')

def stage_per_celltype_plot(adata,celltypekey='louvain',stagekey='batch',plotlabel=True,\
    # this is a function for vertical bar plots
    # please remember to run pl.umap to assign colors
    celltypelist=['default'],stagelist=['default'],celltypekeytype=int,stagekeytype=str,
    fontsize='x-small',xfontsize='x-small',legend_pos=(1,1),savefig=None):
    if 'default' in celltypelist:
        celltypelist = sorted(adata.obs[celltypekey].unique().tolist(),key=celltypekeytype)
    if 'default' in stagelist:
        stagelist = sorted(adata.obs[stagekey].unique().tolist(),key=stagekeytype)
    celltypelist=[i for i in celltypelist if i in adata.obs[celltypekey].unique()]
    stagelist=[i for i in stagelist if i in adata.obs[stagekey].unique()]
    colors=ExtractColor(adata,stagekey,keytype=str)
    count_array=np.array(pd.crosstab(adata.obs[celltypekey],adata.obs[stagekey]).loc[celltypelist,stagelist])
    count_ratio_array=count_array.transpose() / np.sum(count_array,axis=1)
    for i in range(len(stagelist)):
        plt.bar(celltypelist,count_ratio_array[i,:],
            bottom=1-np.sum(count_ratio_array[0:i+1,:],axis=0),
           color=colors[stagelist[i]],label=stagelist[i])
        plt.xticks(fontsize=xfontsize)
    plt.grid(visible=False)
    plt.legend(stagelist,fontsize=fontsize,bbox_to_anchor=legend_pos)
    plt.xticks(rotation=90)
    if savefig is not None:
        plt.savefig(savefig+'.pdf',bbox_inches='tight')

def mtx2df(mtx,idx,col):
    #mtx is the name/location of the matrix.mtx file
    #idx is the index file (rownames)
    #col is the colnames file
    count = scipy.io.mmread(mtx)
    idxs = [i.strip() for i in open(idx)]
    cols = [i.strip() for i in open(col)]
    sc_count = pd.DataFrame(data=count.toarray(),
                        index=idxs,
                        columns=cols)
    return sc_count

def returnDEres(adata, column = None, key= None, remove_mito_ribo = True):
    import functools
    if key is None:
        key = 'rank_genes_groups'
    else:
        key = key

    if column is None:
        column = list(adata.uns[key]['scores'].dtype.fields.keys())[0]
    else:
        column = column

    scores = pd.DataFrame(data = adata.uns[key]['scores'][column], index = adata.uns[key]['names'][column])
    lfc = pd.DataFrame(data = adata.uns[key]['logfoldchanges'][column], index = adata.uns[key]['names'][column])
    pvals = pd.DataFrame(data = adata.uns[key]['pvals'][column], index = adata.uns[key]['names'][column])
    padj = pd.DataFrame(data = adata.uns[key]['pvals_adj'][column], index = adata.uns[key]['names'][column])
    try:
        pts = pd.DataFrame(data = adata.uns[key]['pts'][column], index = adata.uns[key]['names'][column])       
    except:
        pass
    scores = scores.loc[scores.index.dropna()]
    lfc = lfc.loc[lfc.index.dropna()]
    pvals = pvals.loc[pvals.index.dropna()]
    padj = padj.loc[padj.index.dropna()]
    try:
        pts = pts.loc[pts.index.dropna()]
    except:
        pass
    try:
        dfs = [scores, lfc, pvals, padj, pts]
    except:
        dfs = [scores, lfc, pvals, padj]
    df_final = functools.reduce(lambda left, right: pd.merge(left, right, left_index = True, right_index = True), dfs)
    try:
        df_final.columns = ['scores', 'logfoldchanges', 'pvals', 'pvals_adj', 'pts']
    except:
        df_final.columns = ['scores', 'logfoldchanges', 'pvals', 'pvals_adj']
    if remove_mito_ribo:
        df_final = df_final[~df_final.index.isin(list(df_final.filter(regex='^RPL|^RPS|^MRPS|^MRPL|^MT-', axis = 0).index))]
        df_final = df_final[~df_final.index.isin(list(df_final.filter(regex='^Rpl|^Rps|^Mrps|^Mrpl|^mt-', axis = 0).index))]
    return(df_final)

def DEmarkers(adata,celltype,reference,obs,max_out_group_fraction=0.25,\
use_raw=False,length=100,obslist=['percent_mito','n_genes','batch'],\
min_fold_change=2,min_in_group_fraction=0.25,log=True,method='wilcoxon',
embedding='X_umap'):
    celltype=celltype
    sc.tl.rank_genes_groups(adata, obs, groups=[celltype],n_genes=length,
                        reference=reference,method=method,log=log,pts=True)
#    temp=returnDEres(adata,key='rank_genes_groups',column=celltype)
    sc.tl.filter_rank_genes_groups(adata, groupby=obs,\
                    max_out_group_fraction=max_out_group_fraction,
                    min_fold_change=min_fold_change,use_raw=use_raw,
                    min_in_group_fraction=min_in_group_fraction)
    GeneList=pd.DataFrame(adata.uns['rank_genes_groups_filtered']['names']).loc[:,celltype].dropna().head(length).transpose().tolist()
#    temp1=pd.concat([temp,
#        (adata[:,temp.index][adata.obs[obs]==celltype].to_df()>0).mean(axis=0).rename('pct1'),
#        (adata[:,temp.index][adata.obs[obs]==reference].to_df()>0).mean(axis=0).rename('pct2')],
#        axis=1)
    import math
#    GeneList=temp1.loc[(temp1.pvals < 0.05) & (temp1.pct1 >= min_in_group_fraction) & \
#(temp1.logfoldchanges > math.log(min_fold_change)) & (temp1.pct2 <= max_out_group_fraction),:].index.tolist()
    sc.pl.embedding(adata,basis=embedding,color=GeneList+obslist,
           color_map = 'jet',use_raw=use_raw)
    sc.tl.dendrogram(adata, groupby=obs, var_names=GeneList)

    sc.pl.dotplot(adata,var_names=GeneList,
             groupby=obs,use_raw=use_raw,standard_scale='var', dendrogram=True)
    sc.pl.stacked_violin(adata[adata.obs[obs].isin([celltype,reference]),:],var_names=GeneList,groupby=obs,
           swap_axes=True)
    del adata.uns['rank_genes_groups']
    del adata.uns['rank_genes_groups_filtered']
    return GeneList

def GlobalMarkers(adata,obs,max_out_group_fraction=0.25,min_fold_change=2,\
min_in_group_fraction=0.25,use_raw=False,method='wilcoxon'):
    sc.tl.rank_genes_groups(adata,groupby=obs,n_genes=len(adata.var_names),
                  method=method)
    sc.tl.filter_rank_genes_groups(adata,groupby=obs,
                    max_out_group_fraction=max_out_group_fraction,
                    min_fold_change=min_fold_change,use_raw=use_raw,
                    min_in_group_fraction=min_in_group_fraction)
    Markers=pd.DataFrame(adata.uns['rank_genes_groups_filtered']['names'])
    return Markers.apply(lambda x: pd.Series(x.dropna().values))

def HVGbyBatch(adata,batch_key='batch',min_mean=0.0125, max_mean=3, min_disp=0.5,\
min_clustersize=100,genenames=['default']):
    if 'default' in genenames:
        genenames = adata.var_names
    sc.settings.verbosity=0
    batchlist=adata.obs[batch_key].value_counts()
    for key in batchlist[batchlist>min_clustersize].index:
        adata_sample = adata[adata.obs[batch_key]==key,:][:,genenames]
        print(key)
        sc.pp.highly_variable_genes(adata_sample, min_mean=min_mean, max_mean=max_mean, min_disp=min_disp)
        adata.var['highly_variable'+'_'+key]=pd.Series(adata.var_names,\
            index=adata.var_names).isin(adata_sample.var_names[adata_sample.var['highly_variable']])
    sc.settings.verbosity=3
    adata.var['highly_variable_n']=0
    temp=adata.var['highly_variable_n'].astype('int32')
    for key in batchlist[batchlist>min_clustersize].index:
        temp=temp+adata.var['highly_variable_'+key].astype('int32')
    adata.var['highly_variable_n']=temp
    return adata

def HVG_cutoff(adata,range_int=10,cutoff=5000,HVG_var='highly_variable_n',fig_size=(8,6)):
    if range_int>max(adata.var[HVG_var]):
        range_int=max(adata.var[HVG_var])
    HVG_list=[]
    for i in list(range(range_int)):
        HVG_list.append((adata.var[HVG_var]>i).value_counts()[True])
    plt.figure(figsize=fig_size)
    plt.plot(list(range(range_int)),HVG_list)
    plt.axhline(y=cutoff, color='r', linestyle='--')
    HVG_n=next((x for x in reversed(HVG_list) if x >= cutoff), HVG_list[0])
    HVG_i=len(HVG_list) - 1 - next((i for i, x in enumerate(reversed(HVG_list)) if x >= cutoff), len(HVG_list)-1)
    print("The smallest i for intersection to achieve more than "+str(cutoff)+" HVGs is "+str(HVG_i)+" , yielding "+str(HVG_n)+" genes")
    return HVG_i

def Bertie(adata,Resln=1,batch_key='batch'):
    import scrublet as scr
    scorenames = ['scrublet_score','scrublet_cluster_score','bh_pval']
    adata.obs['doublet_scores']=0
    def bh(pvalues):
        '''
        Computes the Benjamini-Hochberg FDR correction.

        Input:
            * pvals - vector of p-values to correct
        '''
        n = int(pvalues.shape[0])
        new_pvalues = np.empty(n)
        values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]
        values.sort()
        values.reverse()
        new_values = []
        for i, vals in enumerate(values):
            rank = n - i
            pvalue, index = vals
            new_values.append((n/rank) * pvalue)
        for i in range(0, int(n)-1):
            if new_values[i] < new_values[i+1]:
                new_values[i+1] = new_values[i]
        for i, vals in enumerate(values):
            pvalue, index = vals
            new_pvalues[index] = new_values[i]
        return new_pvalues

    for i in np.unique(adata.obs[batch_key]):
        print(i)
        adata_sample = adata[adata.obs[batch_key]==i,:]
        scrub = scr.Scrublet(adata_sample.X)
        doublet_scores, predicted_doublets = scrub.scrub_doublets(verbose=False)
        adata_sample.obs['scrublet_score'] = doublet_scores
        adata_sample=adata_sample.copy()
        sc.pp.filter_genes(adata_sample, min_cells=3)
        sc.pp.normalize_per_cell(adata_sample, counts_per_cell_after=1e4)
        sc.pp.log1p(adata_sample)
        sc.pp.highly_variable_genes(adata_sample, min_mean=0.0125, max_mean=3, min_disp=0.5)
        adata_sample = adata_sample[:, adata_sample.var['highly_variable']]
        sc.pp.scale(adata_sample, max_value=10)
        sc.tl.pca(adata_sample, svd_solver='arpack')
        adata_sample = adata_sample.copy()
        # del adata_sample.obsm['X_diffmap']
        sc.pp.neighbors(adata_sample)
        #eoverclustering proper - do basic clustering first, then cluster each cluster
        sc.tl.louvain(adata_sample)
        for clus in np.unique(adata_sample.obs['louvain']):
            sc.tl.louvain(adata_sample, restrict_to=('louvain',[clus]),resolution=Resln)
            adata_sample.obs['louvain'] = adata_sample.obs['louvain_R']
        #compute the cluster scores - the median of Scrublet scores per overclustered cluster
        for clus in np.unique(adata_sample.obs['louvain']):
            adata_sample.obs.loc[adata_sample.obs['louvain']==clus, 'scrublet_cluster_score'] = \
                np.median(adata_sample.obs.loc[adata_sample.obs['louvain']==clus, 'scrublet_score'])
        #now compute doublet p-values. figure out the median and mad (from above-median values) for the distribution
        med = np.median(adata_sample.obs['scrublet_cluster_score'])
        mask = adata_sample.obs['scrublet_cluster_score']>med
        mad = np.median(adata_sample.obs['scrublet_cluster_score'][mask]-med)
        #let's do a one-sided test. the Bertie write-up does not address this but it makes sense
        pvals = 1-scipy.stats.norm.cdf(adata_sample.obs['scrublet_cluster_score'], loc=med, scale=1.4826*mad)
        adata_sample.obs['bh_pval'] = bh(pvals)
        #create results data frame for single sample and copy stuff over from the adata object
        scrublet_sample = pd.DataFrame(0, index=adata_sample.obs_names, columns=scorenames)
        for meta in scorenames:
            scrublet_sample[meta] = adata_sample.obs[meta]
        #write out complete sample scores
        #scrublet_sample.to_csv('scrublet-scores/'+i+'.csv')

        #scrub.plot_histogram();
        #plt.savefig('limb/sample_'+i+'_doulet_histogram.pdf')
        adata.obs.loc[adata.obs[batch_key]==i,'doublet_scores']=doublet_scores
        adata.obs.loc[adata.obs[batch_key]==i,'bh_pval'] = bh(pvals)
        del adata_sample
    return adata

def Bertie_preclustered(adata,batch_key='batch',cluster_key='louvain'):
    import scrublet as scr
    scorenames = ['scrublet_score','scrublet_cluster_score','bh_pval']
    adata.obs['doublet_scores']=0
    def bh(pvalues):
        '''
        Computes the Benjamini-Hochberg FDR correction.

        Input:
            * pvals - vector of p-values to correct
        '''
        n = int(pvalues.shape[0])
        new_pvalues = np.empty(n)
        values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]
        values.sort()
        values.reverse()
        new_values = []
        for i, vals in enumerate(values):
            rank = n - i
            pvalue, index = vals
            new_values.append((n/rank) * pvalue)
        for i in range(0, int(n)-1):
            if new_values[i] < new_values[i+1]:
                new_values[i+1] = new_values[i]
        for i, vals in enumerate(values):
            pvalue, index = vals
            new_pvalues[index] = new_values[i]
        return new_pvalues

    for i in np.unique(adata.obs[batch_key]):
        adata_sample = adata[adata.obs[batch_key]==i,:]
        scrub = scr.Scrublet(adata_sample.X)
        doublet_scores, predicted_doublets = scrub.scrub_doublets(verbose=False)
        adata_sample.obs['scrublet_score'] = doublet_scores
        adata_sample=adata_sample.copy()

        for clus in np.unique(adata_sample.obs[cluster_key]):
            adata_sample.obs.loc[adata_sample.obs[cluster_key]==clus, 'scrublet_cluster_score'] = \
                np.median(adata_sample.obs.loc[adata_sample.obs[cluster_key]==clus, 'scrublet_score'])

        med = np.median(adata_sample.obs['scrublet_cluster_score'])
        mask = adata_sample.obs['scrublet_cluster_score']>med
        mad = np.median(adata_sample.obs['scrublet_cluster_score'][mask]-med)
        #let's do a one-sided test. the Bertie write-up does not address this but it makes sense
        pvals = 1-scipy.stats.norm.cdf(adata_sample.obs['scrublet_cluster_score'], loc=med, scale=1.4826*mad)
        adata_sample.obs['bh_pval'] = bh(pvals)
        #create results data frame for single sample and copy stuff over from the adata object
        scrublet_sample = pd.DataFrame(0, index=adata_sample.obs_names, columns=scorenames)
        for meta in scorenames:
            scrublet_sample[meta] = adata_sample.obs[meta]
        #write out complete sample scores
        #scrublet_sample.to_csv('scrublet-scores/'+i+'.csv')

        #scrub.plot_histogram();
        #plt.savefig('limb/sample_'+i+'_doulet_histogram.pdf')
        adata.obs.loc[adata.obs[batch_key]==i,'doublet_scores']=doublet_scores
        adata.obs.loc[adata.obs[batch_key]==i,'bh_pval'] = bh(pvals)
        del adata_sample
    return adata


def snsSplitViolin(adata,genelist,celltype='leiden',celltypelist=['0','1']):
    df=sc.get.obs_df(adata[adata.obs[celltype].isin(celltypelist)],genelist+[celltype])
    df = df.set_index(celltype).stack().reset_index()
    df.columns=[celltype,'gene','value']
    sns.violinplot(data=df, x='gene', y='value', hue=celltype,
                split=True, inner="quart", linewidth=1)

def DownSample(MouseC1data,cell_type='leiden',downsampleTo=10):
    NewIndex3=[]
    if ( downsampleTo > 0 ) & (isinstance(downsampleTo, int)):
        for i in MouseC1data.obs[cell_type].sort_values().unique():
            NewIndex3=NewIndex3+random.sample(\
                   population=MouseC1data[MouseC1data.obs[cell_type]==i].obs_names.tolist(),
            k=min(downsampleTo,len(MouseC1data[MouseC1data.obs[cell_type]==i\
                         ].obs_names.tolist())))
    return MouseC1data[NewIndex3]

def snsCluster(MouseC1data,MouseC1ColorDict2={False:'#000000',True:'#00FFFF'},cell_type='louvain',gene_type='highly_variable',\
            cellnames=['default'],genenames=['default'],figsize=(10,7),row_cluster=False,col_cluster=False,\
            robust=True,xticklabels=False,yticklabels=False,method='complete',metric='correlation',cmap='jet',\
            downsampleTo=0):
    if 'default' in cellnames:
        cellnames = MouseC1data.obs_names
        if ( downsampleTo > 0 ) & (isinstance(downsampleTo, int)):
            cellnames = DownSample(MouseC1data,cell_type,downsampleTo).obs_names
    if 'default' in genenames:
        genenames = MouseC1data.var_names
    genenames = [i for i in genenames if i in MouseC1data.var_names]
    cellnames = [i for i in cellnames if i in MouseC1data.obs_names]
    cell_types=cell_type
    gene_types=gene_type
    if type(cell_type) == str:
        cell_types=[cell_type]
    if type(gene_type) == str:
        gene_types=[gene_type]
    louvain_col_colors=[]
    for key in cell_types: 
        MouseC1data_df = MouseC1data[MouseC1data.obs_names][:,genenames].to_df()
        MouseC1data_df[key] = MouseC1data[MouseC1data.obs_names].obs[key]
        MouseC1data_df = MouseC1data_df.sort_values(by=key)
        MouseC1data_df3 = MouseC1data_df.loc[pd.Series(cellnames,index=cellnames).index,:]
        cluster_names=MouseC1data_df3.pop(key)
        louvain_col_colors.append(cluster_names.map(ExtractColor(MouseC1data,obsKey=key,keytype=str)).astype(str))
    adata_for_plotting = MouseC1data_df.loc[cellnames,MouseC1data_df.columns.isin(genenames)]
    adata_for_plotting = adata_for_plotting.reindex(columns=genenames)
    if len(louvain_col_colors) > 1:
        louvain_col_colors=pd.concat(louvain_col_colors,axis=1)
    else:
        louvain_col_colors=louvain_col_colors[0]
    if 'null' in gene_types:
        cg1_0point2=sns.clustermap(adata_for_plotting.transpose(),metric=metric,cmap=cmap,\
                 figsize=figsize,row_cluster=row_cluster,col_cluster=col_cluster,robust=robust,xticklabels=xticklabels,\
                 yticklabels=yticklabels,z_score=0,vmin=-2.5,vmax=2.5,col_colors=louvain_col_colors,method=method)
    else:
        celltype_row_colors=[]
        for key in gene_types:
            genegroup_names=MouseC1data[:,genenames].var[key]
            celltype_row_colors.append(genegroup_names.map(MouseC1ColorDict2).astype(str))
        if len(celltype_row_colors) > 1:
            celltype_row_colors=pd.concat(celltype_row_colors,axis=1)
        else:
            celltype_row_colors=celltype_row_colors[0]
        cg1_0point2=sns.clustermap(adata_for_plotting.transpose(),metric=metric,cmap=cmap,\
                 figsize=figsize,row_cluster=row_cluster,col_cluster=col_cluster,robust=robust,xticklabels=xticklabels,\
                 yticklabels=yticklabels,z_score=0,vmin=-2.5,vmax=2.5,col_colors=louvain_col_colors,row_colors=celltype_row_colors,method=method)

    return cg1_0point2

def extractSeabornRows(snsObj):
    #This function returns a Series containing row labels
    NewIndex=pd.DataFrame(np.asarray([snsObj.data.index[i] for i in snsObj.dendrogram_row.reordered_ind])).iloc[:,0]
    return NewIndex

def markSeaborn(snsObj,genes,clustermap=True):
    if clustermap == True:
        NewIndex=pd.DataFrame(np.asarray([snsObj.data.index[i] for i in snsObj.dendrogram_row.reordered_ind])).iloc[:,0]
        NewIndex2=NewIndex.isin(genes)
        snsObj.ax_heatmap.set_yticks(NewIndex[NewIndex2].index.values.tolist())
        snsObj.ax_heatmap.set_yticklabels(NewIndex[NewIndex2].values.tolist())
    else:
        NewIndex=pd.DataFrame(np.asarray(snsObj.data.index))
        NewIndex2=snsObj.data.index.isin(genes)
        snsObj.ax_heatmap.set_yticks(NewIndex[NewIndex2].index.values)
        snsObj.ax_heatmap.set_yticklabels(NewIndex[NewIndex2].values[:,0])
    #snsObj.fig
    return snsObj.fig

def PseudoBulk(adata, group_key, layer=None, gene_symbols=None):
#This function was written by ivirshup
#https://github.com/scverse/scanpy/issues/181#issuecomment-534867254
    if layer is not None:
        getX = lambda x: x.layers[layer]
    else:
        getX = lambda x: x.X
    if gene_symbols is not None:
        new_idx = adata.var[idx]
    else:
        new_idx = adata.var_names

    grouped = adata.obs.groupby(group_key)
    out = pd.DataFrame(
        np.zeros((adata.shape[1], len(grouped)), dtype=np.float64),
        columns=list(grouped.groups.keys()),
        index=adata.var_names
    )

    for group, idx in grouped.indices.items():
        X = getX(adata[idx])
        out[group] = np.ravel(X.mean(axis=0, dtype=np.float64))
    return out

def Dotplot2D(adata,obs1,obs2,gene,cmap='OrRd', min_count=1):
    #This function was modified from K Polanski's codes. It can plot a gene such as XIST across samples and cell types
    #require at least these many cells in a batch+celltype intersection to process it

    #extract a simpler form of all the needed data - the gene's expression and the two obs columns
    #this way things run way quicker downstream
    expression = np.array(adata[:,gene].X)
    batches = adata.obs[obs1].values
    celltypes = adata.obs[obs2].values

    dot_size_df = pd.DataFrame(0.0, index=np.unique(batches), columns=np.unique(celltypes))
    dot_color_df = pd.DataFrame(0.0, index=np.unique(batches), columns=np.unique(celltypes))

    for batch in np.unique(batches):
        mask_batch = (batches == batch)
        for celltype in np.unique(celltypes):
            mask_celltype = (celltypes == celltype)
            #skip if there's not enough data for spot
            if np.sum(mask_batch & mask_celltype) >= min_count:
                sub = expression[mask_batch & mask_celltype]
                #color is mean expression
                dot_color_df.loc[batch, celltype] = np.mean(sub)
                #fraction expressed can be easily computed
                #by making all expressed cells be 1, and then doing a mean again
                sub[sub>0] = 1
                dot_size_df.loc[batch, celltype] = np.mean(sub)

    #reduce dimensions - no need for all-zero rows/cols
    dot_size_df = dot_size_df.loc[(dot_size_df.sum(axis=1) != 0), (dot_size_df.sum(axis=0) != 0)]
    dot_color_df = dot_color_df.loc[(dot_color_df.sum(axis=1) != 0), (dot_color_df.sum(axis=0) != 0)]

    import anndata
    from scanpy.pl import DotPlot

    bdata = anndata.AnnData(np.zeros(dot_size_df.shape))
    bdata.var_names = dot_size_df.columns
    bdata.obs_names = list(dot_size_df.index)
    bdata.obs[obs1] = dot_size_df.index
    bdp = DotPlot(bdata, dot_size_df.columns, obs1, dot_size_df=dot_size_df, dot_color_df=dot_color_df)
    bdp = bdp.style(cmap=cmap)
    bdp.make_figure()


def DeepTree(adata,MouseC1ColorDict2,cell_type='louvain',gene_type='highly_variable',\
            cellnames=['default'],genenames=['default'],figsize=(10,7),row_cluster=True,col_cluster=True,\
            method='complete',metric='correlation',Cutoff=0.8,CladeSize=2):
    if 'default' in cellnames:
        cellnames = adata.obs_names
    if 'default' in genenames:
        genenames = adata.var_names
    test=snsCluster(adata,\
                           MouseC1ColorDict2=MouseC1ColorDict2,\
                           genenames=genenames, cellnames=cellnames,\
                           gene_type=gene_type, cell_type=cell_type,method=method,\
                           figsize=figsize,row_cluster=row_cluster,col_cluster=col_cluster,metric=metric)
    cutree = cluster.hierarchy.cut_tree(test.dendrogram_row.linkage,height=Cutoff)
    TreeDict=dict(zip(*np.unique(cutree, return_counts=True)))
    TreeDF=pd.DataFrame(TreeDict,index=[0])
    DeepIndex=[i in TreeDF.loc[:,TreeDF.iloc[0,:] > CladeSize].columns.values for i in cutree]
    bdata=adata[:,test.data.index][cellnames]
    bdata.var['Deep']=DeepIndex
    test1=snsCluster(bdata,\
                           MouseC1ColorDict2=MouseC1ColorDict2,\
                           cellnames=cellnames,gene_type='Deep',cell_type=cell_type,method=method,\
                           figsize=figsize,row_cluster=True,col_cluster=True,metric=metric)
    test2=snsCluster(bdata[:,DeepIndex],\
                           MouseC1ColorDict2=MouseC1ColorDict2,\
                           cellnames=cellnames,gene_type='null',cell_type=cell_type,method=method,\
                           figsize=figsize,row_cluster=True,col_cluster=True,metric=metric)
    return [bdata,test,test1,test2]

def DeepTree_per_batch(adata,batch_key='batch',obslist=['batch'],min_clustersize=100,Cutoff=0.8,CladeSize=2):
    batchlist=adata.obs[batch_key].value_counts()
    for key in batchlist[batchlist>min_clustersize].index:
        print(key)
        bdata=adata[:,adata.var['highly_variable_'+key]][adata.obs[batch_key]==key,:]
        sc.pp.filter_genes(bdata,min_cells=3)
        sc.pl.umap(bdata,color=obslist)
        [bdata,test, test1, test2]=DeepTree(bdata,
                        MouseC1ColorDict2={False:'#000000',True:'#00FFFF'},
                        cell_type=obslist,
                        cellnames=adata[adata.obs[batch_key]==key,:].obs_names.tolist(),
                        genenames=adata[:,adata.var['highly_variable_'+key]].var_names.tolist(),
                         row_cluster=True,col_cluster=True,Cutoff=Cutoff,CladeSize=CladeSize)
        adata.var['Deep_'+key]=pd.Series(adata.var_names,index=adata.var_names).isin((bdata)[:,bdata.var['Deep']].var_names)
#        sc.pl.umap(adata,color=obslist)
    adata.var['Deep_n']=0
    temp=adata.var['Deep_n'].astype('int32')
    for key in batchlist[batchlist>min_clustersize].index:
        temp=temp+adata.var['Deep_'+key].astype('int32')
    adata.var['Deep_n']=temp
    return adata

def HVG_Venn_Upset(adata,genelists,size_height=3):
    from upsetplot import UpSet
    from upsetplot import plot
    #gene lists can be ['Deep_1','Deep_2']
    deepgenes=pd.DataFrame(adata.var[['highly_variable']+genelists])
    deepgenes=deepgenes.set_index(genelists)
    upset = UpSet(deepgenes, subset_size='count', intersection_plot_elements=size_height)
    return upset

def Treemap(adata,output="temp",branchlist=['project','batch'],width=1000,height=700,title='title'):
    import pandas as pd
    import numpy as np
    temp=adata.obs.groupby(by=branchlist).size()
    temp=temp[temp>0]
    import plotly.express as px
    fig = px.treemap(temp.reset_index(),
                 path=branchlist,values=0)
    fig.update_layout(title=title,
                  width=width, height=height)
    fig.write_image(output+'.pdf')
    temp.to_csv(output+'.csv')
    return fig

def DeepTree2(adata,method='complete',metric='correlation',cellnames=['default'],genenames=['default'],\
               Cutoff=0.8,CladeSize=2):
    if 'default' in cellnames:
        cellnames = adata.obs_names
    if 'default' in genenames:
        genenames = adata.var_names
    adata_df=adata[cellnames,:][:,genenames].to_df()
    testscipy=scipy.cluster.hierarchy.fclusterdata(adata_df.transpose(),\
                    metric=metric,method=method,t=Cutoff,criterion="distance")
    TreeDict=dict(zip(*np.unique(testscipy, return_counts=True)))
    TreeDF=pd.DataFrame(TreeDict,index=[0])
    DeepIndex=[i in TreeDF.loc[:,TreeDF.iloc[0,:] > CladeSize].columns.values for i in testscipy]
    bdata=adata[cellnames,:][:,genenames]
    bdata.var['Deep']=DeepIndex
    return bdata                 

def LoadLogitModel(model_addr):
    import joblib
    lr = joblib.load(open(model_addr,'rb'))
    return lr

def LoadLogitGenes(genecsv):
    CT_genes=pd.read_csv(genecsv,header=None)
    CT_genes['idx'] = CT_genes.index
    CT_genes.columns = ['symbol', 'idx']
    CT_genes = np.array(CT_genes['symbol'])
    return CT_genes

def ExtractLogitScores(adata,model,CT_genes):
    P=model.predict_proba(adata[:,CT_genes].X)
    df=pd.DataFrame(P,index=adata.obs_names,columns=model.classes_+'_prob')
    df['lr_score']=df.max(axis=1)
    return df

def Model2Coeff(model_pkl,genelistcsv):
    lr = LoadLogitModel(model_pkl)
    CT_genes = LoadLogitGenes(genelistcsv)
    return pd.DataFrame(lr.coef_,columns=CT_genes,index=lr.classes_)

from datetime import date
def LogisticRegressionCellType(Reference, Query, Category = 'louvain', DoValidate = False,\
    multi_class='ovr',n_jobs=15,max_iter=1000,C=0.2,tol=1e-4,keyword=''):
    #This function doesn't do normalization or scaling
    #The logistic regression function returns the updated Query object with predicted info stored
    Reference.var_names_make_unique()
    Query.var_names_make_unique()
    IntersectGenes = np.intersect1d(Reference.var_names,Query.var_names)
    Reference2 = Reference[:,IntersectGenes]
    Query2 = Query[:,IntersectGenes]
    X = Reference2.X
    y = Reference2.obs[Category].replace(np.nan,'None',regex=True)
    x = Query2.X
    cv = StratifiedKFold(n_splits=5, random_state=None, shuffle=False)
    logit = LogisticRegression(penalty='l2',
                           random_state=42,
                           C=C,
                           tol=tol,
                           n_jobs=n_jobs,
                           solver='sag',
                           multi_class=multi_class,
                           max_iter=max_iter,
                           verbose=10
                          )
    result=logit.fit(X, y)
    y_predict=result.predict(x)
    today = date.today()
    if DoValidate is True:
        scores = cross_val_score(logit, X, y, cv=cv)
        print(scores)
    _ = joblib.dump(result,str(today)+'Sklearn.result.'+keyword+'.pkl',compress=9)
    np.savetxt(str(today)+'Sklearn.result.'+keyword+'.csv',IntersectGenes,fmt='%s',delimiter=',')
    Query.obs['Predicted'] = y_predict
    return Query

def LogisticPrediction(adata, model_pkl, genelistcsv, scores=True, compute_umap=False,
                      n_neighbors=15, n_pcs=40, metric='euclidean', 
                      min_dist=0.5, spread=1.0, random_state=0):
    """
    This function imports saved logistic model and gene list csv to predict cell types for adata
    adata has better been scaled if you trained a model using scaled AnnData
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data object
    model_pkl : str
        Path to saved logistic regression model
    genelistcsv : str
        Path to gene list CSV file
    scores : bool, default=True
        Whether to return logistic scores
    compute_umap : bool, default=False
        Whether to compute UMAP on logistic scores
    n_neighbors : int, default=15
        Number of neighbors for neighbor graph computation
    n_pcs : int, default=40
        Number of principal components to use
    metric : str, default='euclidean'
        Distance metric for neighbor computation
    min_dist : float, default=0.5
        Minimum distance for UMAP
    spread : float, default=1.0
        Spread parameter for UMAP
    random_state : int, default=0
        Random state for reproducibility
    
    Returns:
    --------
    AnnData object with predictions and optionally logistic scores and UMAP coordinates

    UMAP feature written by Alex :)
    """
    
    CT_genes = LoadLogitGenes(genelistcsv)
    lr = LoadLogitModel(model_pkl)
    lr.features = CT_genes
    features = adata.var_names
    k_x = features.isin(list(CT_genes))
    print(f'{k_x.sum()} features used for prediction', file=sys.stderr)
    k_x_idx = np.where(k_x)[0]
    temp = adata
    X = temp.X[:, k_x_idx]
    features = features[k_x]
    ad_ft = pd.DataFrame(features.values, columns=['ad_features']).reset_index().rename(columns={'index': 'ad_idx'})
    lr_ft = pd.DataFrame(lr.features, columns=['lr_features']).reset_index().rename(columns={'index': 'lr_idx'})
    lr_idx = lr_ft.merge(ad_ft, left_on='lr_features', right_on='ad_features').sort_values(by='ad_idx').lr_idx.values
    lr.n_features_in_ = lr_idx.size
    lr.features = lr.features[lr_idx]
    lr.coef_ = lr.coef_[:, lr_idx]
    predicted_hi = lr.predict(X)
    adata.obs['Predicted'] = predicted_hi
    
    # Extract logistic scores for UMAP computation
    if compute_umap or scores:
        logit_scores = ExtractLogitScores(adata, lr, CT_genes)
    
    # Compute UMAP on logistic scores if requested
    if compute_umap:
        print("Computing UMAP on logistic regression scores...", file=sys.stderr)
        # Create temporary AnnData object with logit scores as X
        temp_adata = sc.AnnData(X=logit_scores)
        temp_adata.obs_names = adata.obs_names
        
        # Compute neighbors on logit scores with custom parameters
        sc.pp.neighbors(temp_adata, use_rep='X', n_neighbors=n_neighbors, 
                       n_pcs=n_pcs, metric=metric, random_state=random_state)
        
        # Compute UMAP on logit scores with custom parameters
        sc.tl.umap(temp_adata, min_dist=min_dist, spread=spread, random_state=random_state)
        
        # Save UMAP coordinates back to original adata
        adata.obsm['X_umap'] = temp_adata.obsm['X_umap']
        print(f"UMAP coordinates (from logit scores) saved to adata.obsm['X_umap'] with shape {adata.obsm['X_umap'].shape}", 
              file=sys.stderr)
    
    # Return with or without scores
    if scores:
        if compute_umap:
            return AddMeta(adata, logit_scores)  # Reuse already computed scores
        else:
            return AddMeta(adata, ExtractLogitScores(adata, lr, CT_genes))
    else:
        return adata

def DouCLing(adata,hi_type,lo_type,rm_genes=[],print_marker_genes=False, fraction_threshold=0.6):
    DoubletScores=pd.DataFrame(0,index=adata.obs[hi_type].unique(),
                          columns=['Parent1','Parent2','Parent1_count','Parent2_count','All_count','p_value'])
    #Dou(blet)C(luster)L(abe)ling method
    #hi_type='leiden_R' is the key for high-resolution cell types that main include doublet clusters
    #lo_type='leiden' is the key for low-resolution cell types that represent compartments
    #this function aims to identify cross-compartment doublet clusters. Same-compartment doublet clusters 
    # look more like transitional cell types closer to homotypic doublets which are difficult to catch
    #rm_genes=['TYMS','MKI67'] is a list of genes you don't want to include, such as cell-cycle genes
    import time
    import scipy.stats as ss
    start_time = time.time()
    alpha=adata.obs[lo_type].value_counts()
    for j in adata.obs[lo_type].sort_values().unique():
        temp=adata[adata.obs[lo_type]==j][:,~adata.var_names.isin(rm_genes)]
        if len(temp.obs[hi_type].value_counts())==1:
            continue
        sc.tl.rank_genes_groups(temp,groupby=hi_type,n_genes=50)
        Markers=pd.DataFrame(temp.uns['rank_genes_groups']['names'])
        for i in pd.DataFrame(temp.uns['rank_genes_groups']['names']).columns:
            scoring_genes=Markers.loc[~(Markers.loc[:,i].isin(rm_genes)),i]
            if print_marker_genes:
                print(scoring_genes)
            if len(scoring_genes)<20:
                continue
            sc.tl.score_genes(adata,gene_list=scoring_genes,score_name=i+'_score')
            DoubletScores.loc[i,'Parent1']=adata[adata.obs[hi_type]==i\
                        ].obs[lo_type].value_counts().index[0]
            cutoff=adata[adata.obs[hi_type]==i].obs[i+'_score'].quantile(q=0.75)

            beta=adata.obs.loc[adata.obs.loc[:,i+'_score'\
                    ]>cutoff,lo_type].value_counts()
            DoubletScores.loc[i,'Parent2']=beta.index[0]
#        if DoubletScores.loc[i,'Parent1']==DoubletScores.loc[i,'Parent2']:
#            pass
#        else:
            DoubletScores.loc[i,'Parent2_count']=beta[0]
            DoubletScores.loc[i,'Parent1_count']=beta.loc[DoubletScores.loc[i,'Parent1']]
            DoubletScores.loc[i,'All_count']=beta.sum()
            hpd=ss.hypergeom(alpha.sum()-alpha.loc[DoubletScores.loc[i,'Parent1']],
                         alpha.loc[DoubletScores.loc[i,'Parent2']],
                        beta.sum()-beta.loc[DoubletScores.loc[i,'Parent1']])
            DoubletScores.loc[i,'p_value']=hpd.pmf(DoubletScores.loc[i,'Parent2_count'])
    print("--- %s seconds ---" % (time.time() - start_time))
    DoubletScores.loc[:,'Is_doublet_cluster']=(DoubletScores.loc[:,'Parent2_count'] / DoubletScores.loc[:,'All_count'] > fraction_threshold) & \
~(DoubletScores.loc[:,'Parent1'] == DoubletScores.loc[:,'Parent2'])
    return DoubletScores

def ClusterGenes(adata,num_pcs=50,embedding='tsne'):
    #adata is already log-transformed
    bdata = adata.copy()
    sc.pp.scale(bdata)
    bdata = bdata.T
    sc.tl.pca(bdata,n_comps=num_pcs)
    sc.pl.pca_variance_ratio(bdata, log=True,n_pcs=50)
    sc.pp.neighbors(bdata,n_pcs=num_pcs)
    if embedding=='umap':
        sc.tl.umap(bdata)
    if embedding=='tsne':
        sc.tl.tsne(bdata)
    sc.tl.leiden(bdata,resolution=0.5)
    #bdata.obs['s_genes'] = [i in s_genes for i in bdata.obs_names]
    #bdata.obs['g2m_genes'] = [i in g2m_genes for i in bdata.obs_names]
    return bdata

def ConvertString(adata, columns_to_convert):
# Convert defined columns in anndata.obs to string, provide a list of columns' names

    for col in columns_to_convert:
        if col in adata.obs.columns:
            adata.obs[col] = adata.obs[col].astype(str)

def PlotCrosstab(
    adata, row_key='celltype', col_key='fine_annotation', 
    title="Confusion Matrix", outfile="confusion_matrix.png"):
#Creates a confusion matrix reindexed and containing all the categories even if they do not exist for consistency.

    row_vals = adata.obs[row_key].astype(str)
    col_vals = adata.obs[col_key].astype(str)

    row_categories = row_vals.unique()
    col_categories = col_vals.unique()

    conf_mat = pd.crosstab(
        row_vals,
        col_vals,
        rownames=[row_key],
        colnames=[col_key]
    ).reindex(index=row_categories, columns=col_categories, fill_value=0)

    fig_height = max(10, 0.5 * len(conf_mat.index))
    fig_width = max(12, 0.3 * len(conf_mat.columns))

    plt.figure(figsize=(fig_width, fig_height))
    sns.heatmap(conf_mat, annot=True, fmt="d", cmap="Blues")
    plt.xlabel(col_key)
    plt.ylabel(row_key)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outfile, format="png")
    plt.show()

def QC(adata, species="human", mt_prefix=None,
    log1p=True, inplace=True):
# Annotates mt, ribo and hb genes based on species name: "human" or "mouse", and calculated QC metrics while storing them in th anndata.obs    
    if mt_prefix is None:
        mt_prefix = "MT-" if species.lower() == "human" else "mt-"

    # Gene category annotations
    adata.var["mt"] = adata.var_names.str.upper().str.startswith(mt_prefix.upper())
    adata.var["ribo"] = adata.var_names.str.upper().str.startswith(("RPS", "RPL"))
    adata.var["hb"] = adata.var_names.str.upper().str.contains(r"^HB[^P]")

    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"], log1p=log1p,inplace=inplace)

def MapCategories(adata, key, corrections_dict):
# Mapping old categories to new names based on corrections_dict and the adata.obs['key'] column. Column needs to be category. Also, adds new categories to category's keys. Requires a mapping dictionary.
 
    adata.obs[key] = adata.obs[key].cat.rename_categories(corrections_dict)

def dropmeta(adata, columns_to_drop):
#Wrapper around adata.obs.drop for deleting specified columns from anndata. Takes a list of columns names.
    
    existing_cols = [col for col in columns_to_drop if col in adata.obs.columns]
    if existing_cols:
        adata.obs.drop(columns=existing_cols, inplace=True)

def SubclusterAll(adata, parent_key="leiden", resolution=0.3, result_key="leiden_R", clusters_to_use=None):
#Creates subclusters of defined clusters or all clusters and provides hierarchical labeling based on the parent's cluster name.
    
    adata.obs[result_key] = adata.obs[parent_key].astype(str)
    all_clusters = adata.obs[parent_key].unique().tolist()
    clusters = clusters_to_use if clusters_to_use is not None else all_clusters
    
# Iterate over each cluster
    for cluster in clusters:
        # Subset the data for the current cluster
        mask = adata.obs[parent_key] == cluster
        sub_adata = adata[mask].copy()

        sc.tl.leiden(sub_adata, resolution=resolution, key_added="leiden_subcluster")

        subcluster_labels = [
            f"{cluster}_{label}" for label in sub_adata.obs["leiden_subcluster"]
        ]

        # Assign hierarchical labels to original adata
        adata.obs.loc[mask, result_key] = subcluster_labels
def OverlapBarcode(adatas, samplename):
    #Identify overlapping barcodes between samples
    import pandas as pd
    
    barcode_dict = {}
    for adata, sample in zip(adatas, samplename):
        barcode_dict[sample] = pd.Index([i.split('-')[0] for i in adata.obs_names])
    
    overlap_matrix = pd.DataFrame(0, index=samplename, columns=samplename)
    for i, sample1 in enumerate(samplename):
        for j, sample2 in enumerate(samplename):
            if i != j:
                overlap = len(set(barcode_dict[sample1]).intersection(barcode_dict[sample2]))
                overlap_matrix.loc[sample1, sample2] = overlap
    
    return {
        'barcode_dict': barcode_dict,
        'overlap_matrix': overlap_matrix,
        'overlap_percentage': overlap_matrix.div(overlap_matrix.sum(axis=1), axis=0) * 100
    }

def show_graph_with_labels(adjacency_matrix, threshold=0.9):
    #Display network graph from adjacency matrix
    import networkx as nx
    import numpy as np
    import matplotlib.pyplot as plt
    
    rows, cols = np.where(adjacency_matrix.values >= threshold)
    gr = nx.Graph()
    gr.add_edges_from(zip(rows.tolist(), cols.tolist()))
    nx.draw(gr, node_size=100, 
            labels={i: adjacency_matrix.index[i] for i in range(len(adjacency_matrix.index))}, 
            with_labels=True)
    plt.show()

def DF2Ann(DF):
    #Convert DataFrame to AnnData
    import anndata
    
    return anndata.AnnData(DF)

def UpSetFromLists(listOflist, labels, size_height=3, showplot=True):
    #Create UpSet plot from list of lists
    from upsetplot import UpSet
    import pandas as pd
    
    listall = list(set([j for i in listOflist for j in i]))
    temp = pd.Series(listall, index=listall)
    temp2 = pd.concat([temp.isin(i) for i in listOflist], axis=1)
    temp2.columns = labels
    temp2 = temp2.set_index(labels)
    upset = UpSet(temp2, subset_size='count', intersection_plot_elements=size_height)
    if showplot:
        upset.plot()
    return upset

def zscore(DF, dropna=True, axis=1):
    #Calculate z-scores for DataFrame
    import scipy.stats
    import pandas as pd
    
    output = pd.DataFrame(scipy.stats.zscore(DF.values, axis=axis),
                         index=DF.index, columns=DF.columns)
    return output.dropna(axis=1-axis) if dropna else output

def Ginni(beta):
    #Calculate Gini coefficient for distributions
    import numpy as np
    
    beta.loc['Ginni'] = [np.abs(np.subtract.outer(beta.loc[:,i].values, beta.loc[:,i].values)).mean() / 
                         np.mean(beta.loc[:,i].values) * 0.5 for i in beta.columns]
    return beta

def cellphonedb_p2adjMat(cpdb_p_loc='pvalues.txt', pval=0.05):
    #Convert CellPhoneDB p-values to adjacency matrix with significant pairs
    import pandas as pd
    
    cpdb_p = pd.read_csv(cpdb_p_loc, sep='\t')
    cellpairs = pd.Series(cpdb_p.iloc[0, 11:].index).str.split('|', expand=True)
    cell0, cell1 = cellpairs.loc[:,0].unique(), cellpairs.loc[:,1].unique()
    
    InterMat = pd.DataFrame('', index=cell0, columns=cell1)
    for i in cell0:
        for j in cell1:
            InterMat.loc[i,j] = '|'.join(cpdb_p.loc[cpdb_p.loc[:,i+'|'+j]<pval, 'interacting_pair'].tolist())
    
    InterMat.to_csv(cpdb_p_loc+'.pairs'+str(pval)+'.csv')
    return InterMat

def cellphonedb_n_interaction_Mat(cpdb_p_loc='pvalues.txt', pval=0.05):
    #Count significant interactions between cell types
    import pandas as pd
    
    cpdb_p = pd.read_csv(cpdb_p_loc, sep='\t')
    cellpairs = pd.Series(cpdb_p.iloc[0, 11:].index).str.split('|', expand=True)
    cell0, cell1 = cellpairs.loc[:,0].unique(), cellpairs.loc[:,1].unique()
    
    InterMat = pd.DataFrame(0, index=cell0, columns=cell1)
    for i in cell0:
        for j in cell1:
            InterMat.loc[i,j] = len(cpdb_p.loc[cpdb_p.loc[:,i+'|'+j]<pval, 'interacting_pair'])
    
    InterMat.to_csv(cpdb_p_loc+'.n_pairs'+str(pval)+'.csv')
    return InterMat

def cellphonedb_mat_per_interaction(interacting_pair, cpdb_p_loc='pvalues.txt'):
    #Get p-value matrix for specific interaction pair
    import pandas as pd
    import numpy as np
    
    cpdb_p = pd.read_csv(cpdb_p_loc, sep='\t')
    cellpairs = pd.Series(cpdb_p.iloc[0, 11:].index).str.split('|', expand=True)
    cell0, cell1 = cellpairs.loc[:,0].unique(), cellpairs.loc[:,1].unique()
    
    InterMat = pd.DataFrame(0, index=cell0, columns=cell1)
    for i in cell0:
        for j in cell1:
            values = cpdb_p.loc[cpdb_p['interacting_pair']==interacting_pair, i+'|'+j].values
            InterMat.loc[i,j] = values[0] if len(values) > 0 else np.nan
    
    InterMat.to_csv(cpdb_p_loc+'.pairs'+interacting_pair.replace('/', '_')+'.csv')
    return InterMat

def ListsOverlap(A, B):
    #Calculate overlap counts between lists in A and B
    import pandas as pd
    
    counts = [[len(set(a_list).intersection(b_list)) for b_list in B] for a_list in A]
    return pd.DataFrame(counts, 
                       columns=[f"B{i+1}" for i in range(len(B))], 
                       index=[f"A{i+1}" for i in range(len(A))])

def read_folder(base_dir, data_pattern="filtered_feature_bc_matrix.h5", dirs=None,
                dir_pattern=None, dir_filter=None, batch_key='sample', concat=False,
                save_counts=True, verbose=True, **kwargs):
    #Load data from multiple directories
    import fnmatch
    import os
    import scanpy as sc
    import gc
    
    # Get subdirectories
    all_subdirs = [d for d in os.listdir(base_dir) 
                  if os.path.isdir(os.path.join(base_dir, d))]
    
    # Filter directories
    if dirs is not None:
        filtered_dirs = [d for d in all_subdirs if d in dirs]
    elif dir_pattern is not None:
        filtered_dirs = fnmatch.filter(all_subdirs, dir_pattern)
    elif dir_filter is not None:
        filtered_dirs = [d for d in all_subdirs if dir_filter(d)]
    else:
        filtered_dirs = all_subdirs
    
    if verbose:
        print(f"Found {len(filtered_dirs)} directories to process")
    
    adatas = []
    for idx, dir_name in enumerate(filtered_dirs):
        if verbose:
            print(f"Processing {idx+1}/{len(filtered_dirs)}: {dir_name}")
        
        dir_path = os.path.join(base_dir, dir_name)
        data_file = None
        
        # Find data file
        potential_file = os.path.join(dir_path, data_pattern)
        if os.path.exists(potential_file):
            data_file = potential_file
        else:
            for file in os.listdir(dir_path):
                if fnmatch.fnmatch(file, data_pattern):
                    data_file = os.path.join(dir_path, file)
                    break
        
        if data_file is None:
            if verbose:
                print(f"Warning: No data file found in {dir_path}. Skipping.")
            continue
        
        # Load data
        try:
            if data_file.endswith('.h5'):
                adata = sc.read_10x_h5(data_file, **kwargs)
            elif data_file.endswith(('.mtx', '.mtx.gz')):
                adata = sc.read_10x_mtx(os.path.dirname(data_file), **kwargs)
            else:
                adata = sc.read(data_file, **kwargs)
        except Exception as e:
            if verbose:
                print(f"Error loading {data_file}: {e}. Skipping.")
            continue
        
        adata.var_names_make_unique()
        if save_counts:
            adata.layers["counts"] = adata.X.copy()
        adata.obs[batch_key] = dir_name
        adatas.append(adata)
        
        if verbose:
            print(f"Added {dir_name}: {adata.n_obs} cells, {adata.n_vars} genes")
        gc.collect()
    
    if concat and adatas:
        if verbose:
            print("Concatenating datasets...")
        try:
            return sc.concat(adatas, join='outer', index_unique='-')
        except Exception as e:
            print(f"Concatenation failed: {e}. Returning list.")
            return adatas
    return adatas

def read_numbered_folders(base_dir, min_dir=1, max_dir=100, exclude_dirs=None, **kwargs):
    #Load data from numbered directories
    all_dirs = [str(i) for i in range(min_dir, max_dir + 1)]
    
    if exclude_dirs is not None:
        exclude_dirs_str = [str(d) for d in exclude_dirs]
        filtered_dirs = [d for d in all_dirs if d not in exclude_dirs_str]
    else:
        filtered_dirs = all_dirs
    
    return read_folder(base_dir, dirs=filtered_dirs, **kwargs)

def run_celltypist(adata, model_name, output_prefix=None, majority_voting=True, 
                   plot_umap=True, plot_size=10, figsize=(16, 16), legend_loc='on data',
                   color_map=None, save_plot=True, plot_file_prefix=None, dpi=300, **kwargs):
    #Run CellTypist annotation
    try:
        import celltypist
        from celltypist import models
    except ImportError:
        raise ImportError("Install celltypist: pip install celltypist")
    
    import matplotlib.pyplot as plt
    import matplotlib.colors
    import scanpy as sc
    
    prefix = output_prefix if output_prefix else ""
    
    if plot_umap and 'X_umap' not in adata.obsm_keys():
        print("Warning: No UMAP found. Skipping plots.")
        plot_umap = False
    
    # Load model
    try:
        model = models.Model.load(model=model_name)
    except:
        try:
            model = models.Model.load(model_file=model_name)
        except:
            raise ValueError(f"Cannot load model '{model_name}'")
    
    # Run annotation
    result = celltypist.annotate(adata, model=model, majority_voting=majority_voting, **kwargs)
    adata_annotated = result.to_adata()
    
    # Copy results
    adata_annotated.obs[f'{prefix}predicted_labels'] = adata_annotated.obs['predicted_labels'].copy()
    if majority_voting:
        adata_annotated.obs[f'{prefix}majority_voting'] = adata_annotated.obs['majority_voting'].copy()
    
    # Plot if requested
    if plot_umap:
        if color_map is None:
            color_map = list(matplotlib.colors.CSS4_COLORS.values())
        
        plot_prefix = plot_file_prefix or prefix or "celltypist"
        
        for col, suffix in [(f'{prefix}predicted_labels', 'predicted'),
                           (f'{prefix}majority_voting', 'majority')] if majority_voting else \
                          [(f'{prefix}predicted_labels', 'predicted')]:
            fig, ax = plt.subplots(figsize=figsize)
            sc.pl.umap(adata_annotated, color=[col], ax=ax, size=plot_size, 
                      legend_loc=legend_loc, palette=color_map, show=False)
            plt.tight_layout()
            if save_plot:
                plt.savefig(f'{plot_prefix}_{suffix}_umap.png', dpi=dpi, bbox_inches='tight')
            plt.show()
    
    return adata_annotated

def extract_barcode_prefix(barcode):
    #Extract barcode prefix before '-'
    return barcode.split('-')[0] if '-' in barcode else barcode

def import_souporcell(adata, clusters_file, genome_file=None, counts_matrix_file=None,
                     prefix='', barcode_col=0, assignment_col=1, status_col=2,
                     add_status=True, doublet_status='doublet', unassigned_status='unassigned',
                     verbose=True):
    #Import Souporcell results into AnnData
    import os
    import pandas as pd
    
    # Check files
    if not os.path.exists(clusters_file):
        raise FileNotFoundError(f"Clusters file not found: {clusters_file}")
    
    if verbose:
        print(f"Reading Souporcell clusters: {clusters_file}")
    
    # Read clusters
    try:
        clusters_df = pd.read_csv(clusters_file, sep='\t', header=None)
        barcodes = clusters_df.iloc[:, barcode_col].values
        assignments = clusters_df.iloc[:, assignment_col].values
        status = clusters_df.iloc[:, status_col].values if add_status and status_col < clusters_df.shape[1] else None
    except Exception as e:
        if verbose:
            print(f"Error reading file: {e}")
        return adata
    
    # Process barcodes
    processed_barcodes = [extract_barcode_prefix(bc) for bc in barcodes]
    barcode_to_assignment = dict(zip(processed_barcodes, assignments))
    
    # Add assignments
    cluster_col = f"{prefix}souporcell_cluster"
    adata.obs[cluster_col] = "unassigned"
    
    cells_assigned = 0
    for barcode in adata.obs_names:
        bc_key = extract_barcode_prefix(barcode)
        if bc_key in barcode_to_assignment:
            adata.obs.loc[barcode, cluster_col] = barcode_to_assignment[bc_key]
            cells_assigned += 1
    
    # Add status if requested
    if add_status and status is not None:
        barcode_to_status = dict(zip(processed_barcodes, status))
        status_col_name = f"{prefix}souporcell_status"
        adata.obs[status_col_name] = unassigned_status
        
        for barcode in adata.obs_names:
            bc_key = extract_barcode_prefix(barcode)
            if bc_key in barcode_to_status:
                adata.obs.loc[barcode, status_col_name] = barcode_to_status[bc_key]
        
        adata.obs[f"{prefix}is_doublet"] = (adata.obs[status_col_name] == doublet_status)
        
        if verbose:
            doublet_count = adata.obs[f"{prefix}is_doublet"].sum()
            print(f"Found {doublet_count} doublets")
    
    if verbose:
        print(f"Assigned {cells_assigned} cells to clusters")
    
    return adata
def iStar_extract_locs_raw(adata, output_path="locs-raw.tsv"):
    """
    Extracts spatial coordinates and saves as `locs-raw.tsv`.

    Parameters
    ----------
    adata : AnnData
        Annotated data object with spatial info in `adata.obsm["spatial"]`.
    output_path : str
        Path to save the output TSV file.

    Returns
    -------
    str : Output path
    """
    import pandas as pd

    df = pd.DataFrame(adata.obsm["spatial"], index=adata.obs.index, columns=["x", "y"])
    df.insert(0, "Spot_ID", df.index)
    df.to_csv(output_path, sep="\t", index=False)
    return output_path


def iStar_extract_counts_raw(adata, output_path="cnts.tsv"):
    """
    Extracts expression matrix and saves as `cnts.tsv`.

    Parameters
    ----------
    adata : AnnData
        AnnData object with `.X` as count matrix.
    output_path : str
        Path to save the output TSV file.

    Returns
    -------
    str : Output path
    """
    import pandas as pd

    matrix = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X
    df = pd.DataFrame(matrix, index=adata.obs.index, columns=adata.var.index)
    df.insert(0, "Spot_ID", df.index)
    df.to_csv(output_path, sep="\t", index=False)
    return output_path


def iStar_extract_scalefactors(adata, formatted_name, pixel_file="pixel-size-raw.txt", radius_file="radius-raw.txt"):
    """
    Extracts pixel size and radius from spatial scalefactors.

    Parameters
    ----------
    adata : AnnData
        AnnData object with `adata.uns['spatial'][formatted_name]['scalefactors']`.
    formatted_name : str
        Sample key used in `adata.uns['spatial']`.
    pixel_file : str
        Output file for pixel size.
    radius_file : str
        Output file for radius.

    Returns
    -------
    tuple[str, str] : Paths to the saved files
    """
    scalefactors = adata.uns['spatial'][formatted_name]['scalefactors']
    pixel_size_raw = (8000 / 2000) * scalefactors["tissue_hires_scalef"]
    radius_raw = scalefactors["spot_diameter_fullres"] * 0.5

    with open(pixel_file, "w") as f:
        f.write(str(pixel_size_raw))
    with open(radius_file, "w") as f:
        f.write(str(radius_raw))

    return pixel_file, radius_file


def iStar_generate_all_inputs(h5ad_path, formatted_name, output_dir="."):
    """
    Main function to extract all iStar input files from a `.h5ad` file.

    Parameters
    ----------
    h5ad_path : str
        Path to the input `.h5ad` file.
    formatted_name : str
        Sample name key used inside `adata.uns["spatial"]`.
    output_dir : str
        Output directory for the generated files.

    Returns
    -------
    dict : File paths for each output.
    """
    import scanpy as sc
    import os

    if not h5ad_path.endswith(".h5ad"):
        h5ad_path += ".h5ad"
    adata = sc.read_h5ad(h5ad_path)

    os.makedirs(output_dir, exist_ok=True)

    locs_path = os.path.join(output_dir, "locs-raw.tsv")
    cnts_path = os.path.join(output_dir, "cnts.tsv")
    pixel_path = os.path.join(output_dir, "pixel-size-raw.txt")
    radius_path = os.path.join(output_dir, "radius-raw.txt")

    iStar_extract_locs_raw(adata, locs_path)
    iStar_extract_counts_raw(adata, cnts_path)
    iStar_extract_scalefactors(adata, formatted_name, pixel_path, radius_path)

    return {
        "locs-raw.tsv": locs_path,
        "cnts.tsv": cnts_path,
        "pixel-size-raw.txt": pixel_path,
        "radius-raw.txt": radius_path,
    }
