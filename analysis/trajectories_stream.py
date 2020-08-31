#!/usr/bin/env python3

### Here we perform Clustering Analysis and find Markers of each cluster by using Seurat V.3.0.1 R package 
### Download counts matrix here: 


SINGLE CELL TRAJECTORIES WITH STREAM - THIS VERSION USES THE ALGORITHM 'SE' IN DIMENSION REDUCTION (FASTER THAN MLLE AND PRESERVES THE STRUCTURE)

##### THIS VERSION USES THE NORMALIZED DATA FROM SEURAT (LOG NORMALIZED)


import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()
import stream as st
import os.path
import pickle



#### Read Counts table
adata=st.read(file_name='./counts.tsv', workdir='./')
# Read Cell labels table
st.add_cell_labels(adata,file_name='./cell_label.tsv')
# Add random colors to each sample
st.add_cell_colors(adata, file_name='./cell_color.tsv')


### CHECK FOR VARIABLE GENES
# Check if the blue (loess) curve fits the points well
st.select_variable_genes(adata)
# Open plot file
plt.savefig('loess.png')
# Close Plot 
plt.close('loess.png')


# Adjust the blue curve to fits better
st.select_variable_genes(adata,loess_frac=0.01)
plt.savefig('adjust_loess.png')
plt.close('adjust_loess.png')


# Dimension reduction
st.dimension_reduction(adata, method = 'se', n_components=3, nb_pct =0.01)

#st.plot_dimension_reduction(adata)
plt.savefig('se_reduction.png')
plt.close('se_reduction.png')


# Plot UMAP
st.plot_visualization_2D(adata)
plt.savefig('umap.png')
plt.close('umap.png')


# Clustering to get branches
st.seed_elastic_principal_graph(adata, clustering = "kmeans")

st.plot_branches(adata)
plt.savefig('branches.png')
plt.close('branches.png')

st.plot_branches_with_cells(adata)
plt.savefig('branches_and_cells.png')
plt.close('branches_and_cells.png')


# Elastic principal graph
st.elastic_principal_graph(adata)

st.plot_branches(adata)
plt.savefig('epg_branches.png')
plt.close('epg_branches.png')

st.plot_branches_with_cells(adata)
plt.savefig('epg_branches_and_cells.png')
plt.close('epg_branches_and_cells.png')

# Optimizing step
st.optimize_branching(adata)

st.plot_branches(adata)
plt.savefig('OPTepg_branches.png')
plt.close('OPTepg_branches.png')

st.plot_branches_with_cells(adata)
plt.savefig('OPTepg_branches_and_cells.png')
plt.close('OPTepg_branches_and_cells.png')


# Extend leaf branch to reach further cells
st.extend_elastic_principal_graph(adata,epg_trimmingradius=0.1)
st.plot_branches(adata)
plt.savefig('EXTepg_branches.png')
plt.close('EXTepg_branches.png')

st.plot_branches_with_cells(adata)
plt.savefig('EXTepg_branches_and_cells.png')
plt.close('EXTepg_branches_and_cells.png')


# Plot flat tree
st.plot_flat_tree(adata,fig_legend_ncol=6,fig_size=(15,15))
plt.savefig('flat_tree.png')
plt.close('flat_tree.png')

# Validate the learned structure by visualizing the branch assignment
st.plot_visualization_2D(adata,fig_legend_ncol=3, fig_size=(15,15))
plt.savefig('branch_assign.png')
plt.close('branch_assign.png')

st.plot_visualization_2D(adata,color_by='branch',fig_legend_ncol=3, fig_size=(15,15))
plt.savefig('branch_assign_color.png')
plt.close('branch_assign_color.png')

##### AT THIS STEP WE MUST CHECK WHICH BRANCH CONTAINS THE iPS CELLS IN ORDER TO CHOOSE THE INITIAL STATE - SEE THE FIGURE 'flat_tree.png'
##### LET'S SUPPOSE THE BRANCH 'S1' HOSTS THE iPS CELLS

#### Subway map
st.subwaymap_plot(adata,root='S1',fig_legend_ncol=6, fig_size=(15,15)) 
plt.savefig('subwaymap.png')
plt.close('subwaymap.png')

#### Subway map - increase distance
st.subwaymap_plot(adata,root='S1',fig_legend_ncol=6,percentile_dist=100, fig_size=(15,15)) 
plt.savefig('subwaymap_pd.png')
plt.close('subwaymap_pd.png')

#### Smooth pseudotime
st.stream_plot(adata,root='S1',fig_legend_ncol=3,fig_size=(15,15),flag_log_view=True, factor_zoomin=200)
plt.savefig('smooth_pseudo.png')
plt.close('smooth_pseudo.png')

### Detect transition genes
st.detect_transistion_genes(adata,root='S1', n_jobs=16)

### Detect leaf genes (marker genes for each branch)
st.detect_leaf_genes(adata,root='S1')

### Detect DE genes among branches
st.detect_de_genes(adata,root='S1')
st.plot_de_genes(adata)
plt.savefig('deg_branches.png', dpi=400)
plt.close('deg_branches.png')
