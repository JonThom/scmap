## scmap
script for annotating clusters in a test set with an annotated reference set using scmap

## dependencies 
devtools, here, optparse, parallel, dplyr, ggplot2, Seurat 3, SingleCellExperiment, scmap

## return 
* list of single cell experiment objects (RDS.gz)
* full scmap results (RDS.gz)
* mapping results as table (.csv)
* test dataset as Seurat Object with new cell labels added to metadata (RDS.gz)
* heatmap(s) of features selected for mapping
* UMAP plot of test dataset with new cell labels

## usage

time Rscript ./scmap_pipeline.R --path_datExpr_test 'c("perslab"="/projects/jonatan/pub-perslab/18-liver-fred/output/liver_perslab_int_seurat_7_SCTint_finalLabels_seuratObj.RDS.gz")'  --metadata_test_id_col cluster_perslab  --path_datExpr_ref 'c("macparland" = "/projects/jonatan/pub-perslab/18-liver-fred/data/macparland_seurat_obj3_Samples245.RDS.gz")' --metadata_ref_id_col Cluster_annot  --threshold 0.5  --prefix_run scmap_perslab_macparland_test_1 --dir_out /projects/jonatan/190910_scmapTest/ --RAM_Gb_max 250

## info on parameters
Rscript ./scmap_pipeline.R --help 
