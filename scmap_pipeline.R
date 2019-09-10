# Script to map unannotated clusters in test dataset to annotated clusters in one or more reference sets using scmap

#' @usage time Rscript /projects/jonatan/tools/scmap/scmap_pipeline.R --path_datExpr_test 'c("perslab"="/projects/jonatan/pub-perslab/18-liver-fred/output/liver_perslab_int_seurat_7_SCTint_finalLabels_seuratObj.RDS.gz")'  --metadata_test_id_col cluster_perslab  --path_datExpr_ref 'c("macparland" = "/projects/jonatan/pub-perslab/18-liver-fred/data/macparland_seurat_obj3_Samples245.RDS.gz")' --metadata_ref_id_col Cluster_annot  --threshold 0.5  --prefix_run scmap_perslab_macparland_test_1 --dir_out /projects/jonatan/190910_scmapTest/ --RAM_Gb_max 250
#' @depends devtools, here, optparse, parallel, dplyr, ggplot2, Seurat 3, SingleCellExperiment, scmap
#' @return list of single cell experiment objects (RDS.gz)
#' @return full scmap results (RDS.gz)
#' @return mapping results as table (.csv)
#' @return test dataset as Seurat Object with new cell labels added to metadata (RDS.gz)
#' @return heatmap(s) of features selected for mapping
#' @return UMAP plot of test dataset with new cell labels
#' @author Jonatan Thompson, rkm at ku dot dk
#' 
#' TODO 
#' * Get Sankey plot working (googleVis)

# References
# * SingleCellExperiment bioconductor: https://bioconductor.org/packages/release/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html
# * scmap vignette: http://bioconductor.org/packages/release/bioc/vignettes/scmap/inst/doc/scmap.html
# * scmap article: https://www.nature.com/articles/nmeth.4644

######################################################################
########################### OptParse #################################
######################################################################
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(optparse))

option_list <- list(

  make_option("--path_datExpr_test", type="character",
              help = "Quoted named vector of length 1 with full path to test data in .RData or .RDS format, optionally gzip compressed, containing gene*cell matrix, data.frame or seurat object with @data"), #TODO:do we need e.g. var.genes?
  make_option("--path_metadata_test", type="character", default=NULL,
            help = "character with full path to test metadata in .RData or .RDS format, optionally gzip compressed, data.frame. Not required if datExpr_test is a seurat object with metadata set as Idents"), #TODO:do we need e.g. var.genes?
  make_option("--metadata_test_id_col", type="character", default= NULL,
              help = "If path_metadata_test is provided, also provide the name of the column containing cell identities"), #TODO:do we need e.g. var.genes?
  make_option("--path_datExpr_ref", type="character",
              help = "Quoted named vector of length 1 with full path to reference expression data in .RData or .RDS format, optionally gzip compressed, containing gene*cell matrix, data.frame or seurat 3 object with logNormalized data"),
  make_option("--path_metadata_ref", type="character", default=NULL,
              help = "If cell identities are not included in the expression data file, a character with the full path to reference set metadata in .RData or .RDS format, optionally gzip compressed, data.frame. Not required if datExpr_ref is a Seurat 3 object with cell identities in the Idents slot"), #TODO:do we need e.g. var.genes?
  make_option("--metadata_ref_id_col", type="character", default=NULL,
              help = "If cell identities are provided in path_metadata_ref, character giving the column name"), #TODO:do we need e.g. var.genes?
  make_option("--vec_n_features", type="character", default='c(500)',
              help = "Quoted vector. How many features to use in mapping cells to clusters. Can give several values to try, e.g. 'c(500,1000,1500)' 500 recommended for best sensitivity/specificity balance but do inspect plots https://www.nature.com/articles/nmeth.4644, [default %default]"), #TODO:do we need e.g. var.genes?
  make_option("--threshold", type="numeric", default=0.5,
              help = "Threshold parameter to pass to scmapCluster(), [default %default]"), #TODO:do we need e.g. var.genes?
  make_option("--prefix_run", type="character", default = "1",
              help = "Run prefix for output files, [default %default]"),
  make_option("--dir_out", type="character",
              help = "Project directory to which to write files. Should include subdirectories /tables, /RObjects, /plots, /log, else they will be created"),
  make_option("--RAM_Gb_max", type="integer", default=250,
              help = "Upper limit on Gb RAM available. Taken into account when setting up parallel processes. [default %default]")
)


######################################################################
######################### UTILITY FUNCTIONS ##########################
######################################################################

source(file=here("perslab-sc-library","utility_functions.R"))
source(file=here("perslab-sc-library","functions_sc.R"))

######################################################################
########################### PACKAGES #################################
######################################################################

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(scmap))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(SingleCellExperiment))
#library(googleVis) #required for Sankey plot. Currently not implemented

######################################################################
########################### GET OPTIONS ##############################
######################################################################

opt <- parse_args(OptionParser(option_list=option_list))

path_datExpr_test <- eval(parse(text=opt$path_datExpr_test))
path_metadata_test <- opt$path_metadata_test
#if (!is.null(path_metadata_test)) path_metadata_test <- eval(parse(text=path_metadata_test))
metadata_test_id_col <- opt$metadata_test_id_col
path_datExpr_ref <- eval(parse(text=opt$path_datExpr_ref))
path_metadata_ref <- opt$path_metadata_ref
metadata_ref_id_col <- opt$metadata_ref_id_col
#if (!is.null(metadata_ref_id_cols)) metadata_ref_id_cols <- eval(parse(text=metadata_ref_id_cols))
vec_n_features <- eval(parse(text=opt$vec_n_features))
threshold <- opt$threshold
prefix_run <- opt$prefix_run
dir_out <- opt$dir_out
RAM_Gb_max <- opt$RAM_Gb_max

######################################################################
########################### SET OPTIONS ##############################
######################################################################

options(stringsAsFactors = F, na.action="na.omit", warn=1)

######################################################################
############################ CONSTANTS ###############################
######################################################################

if (!is.null(names(path_datExpr_test))) prefix_test <- names(path_datExpr_test) else stop("path_datExpr_test must be a named vector of length one")
if (!is.null(names(path_datExpr_ref))) prefix_ref <- names(path_datExpr_ref) else stop("path_datExpr_ref must be a named vector of length one")

flag_date = substr(gsub("-","",as.character(Sys.Date())),3,1000)

# if specified output directory doesn't exist, create it
if (!file.exists(dir_out)) {
  dir.create(dir_out)
  message("Project directory not found, new one created")
}

dir_plots = paste0(dir_out,"plots/")
if (!file.exists(dir_plots)) dir.create(dir_plots)

dir_tables = paste0(dir_out,"tables/")
if (!file.exists(dir_tables)) dir.create(dir_tables)

dir_RObjects = paste0(dir_out,"RObjects/")
if (!file.exists(dir_RObjects)) dir.create(dir_RObjects)

dir_log = paste0(dir_out,"log/")
if (!file.exists(dir_log)) dir.create(dir_log)

randomSeed = 12345
set.seed(randomSeed)
######################################################################
########################## LOAD METADATA ############################
######################################################################

ident_test <- NULL
ident_ref <- NULL

if (!is.null(path_metadata_test)) {
  metadata_test <- load_obj(path_metadata_test)
  # if (any(grepl(pattern="X", x = colnames(metadata_test), fixed = T))){
  #   rownames(metadata_test) <- metadata_test[["X"]]
  # }
  ident_test <- metadata_test[[metadata_test_id_col]]
}

if (!is.null(path_metadata_ref)) {
  ident_ref <- #mapply(function(path_metadata_ref, id_col) {
    #if (!is.na(path_metadata_ref)) {
      tryCatch({
        metadata_ref_tmp <- load_obj(path_metadata_ref)}, error= function(err){stop(paste0(path_metadata_ref, " not found"))})
      if (any(grepl(pattern="X", x = colnames(metadata_ref_tmp), fixed = T))) {
        rownames(metadata_ref_tmp) <- metadata_ref_tmp[["X"]]  # get default rowname column
      }
      ident_ref <- metadata_ref_tmp[[metatadata_ref_id_col]]
      names(ident_ref) <- rownames(metadata_ref_tmp)
      #return(ident_ref)
    #} else {
    #  NA_character_
    #}
  #},path_metadata_ref = paths_metadata_ref, id_col=metadata_ref_id_cols, SIMPLIFY=F)
} else {
  ident_ref <- NULL#as.list(rep(NA_character_, times=length(paths_datExpr_ref)))
}

#names(ident_ref) <- names(paths_datExpr_ref)
######################################################################
####################### LOAD EXPRESSION MATRICES #####################
######################################################################

message("Loading test and reference expression matrices")

fun = function(path_datExpr) {
  datExpr_tmp <- load_obj(path_datExpr)
  if ("Seurat" %in% class(datExpr_tmp)) {
    #metadata <- datExpr_tmp@meta.data
    datExpr <- GetAssayData(datExpr_tmp, slot="data", assay="RNA")
    ident <- Idents(datExpr_tmp) %>% as.character
    names(ident) <- colnames(datExpr_tmp)
    #datExpr <- datExpr_tmp@data
  } else {
    datExpr <- datExpr_tmp
    if (any(grepl(pattern="X", x = colnames(datExpr), fixed = T))) {
      rownames(datExpr) <- datExpr[["X"]]
      datExpr[["X"]] <- NULL
    }
  }
  rm(datExpr_tmp)
  return(list("datExpr"=datExpr, "ident"=ident))
}

### Test data ###

list_data_test <- fun(path_datExpr_test)
datExpr_test <- list_data_test[["datExpr"]]

if (all(is.null(ident_test))) {
  ident_test <- list_data_test[["ident"]]
  if (all(is.null(ident_test)) | all(is.na(ident_test))) stop("provide cell identity metadata for test data set either in seurat object or separate file")
}

rm(list_data_test)

### Ref data ###

list_data_ref <- fun(path_datExpr_ref)
datExpr_ref <- list_data_ref[["datExpr"]]

if (all(is.null(ident_ref))) {
  ident_ref <- list_data_ref[["ident"]]
  if (all(is.null(ident_ref)) | all(is.na(ident_ref))) stop("provide cell identity metadata for ref data set either in seurat object or separate file")
}

rm(list_data_ref)

# list_iterable=list("X"=paths_datExpr_ref)
# list_list_data_ref <- safeParallel(fun=fun, list_iterable=list_iterable)
#
# list_datExpr_ref <- lapply(list_list_data_ref, function(list_data_ref) list_data_ref[["datExpr"]] %>% as.matrix)
#
# list_ident_ref <- mapply(function(list_data_ref, ident_ref) {
#     if (all(is.na(ident_ref))) {
#       ident_ref <- list_data_ref[["ident"]]
#     }
#     return(ident_ref)
#   },
#   list_data_ref=list_list_data_ref,
#   ident_ref = list_ident_ref,
#   SIMPLIFY=F)
#
#
# if (any(sapply(list_ident_ref, is.na))) stop("provide metadata for reference sets either in seurat objects or separate file")
# rm(list_list_data_ref)

######################################################################
#################### MATCH DATEXPR AND METADATA ######################
######################################################################

fun = function(colnames_datExpr, ident) {
  ident <- ident[as.integer(na.omit(match(names(ident), colnames_datExpr)))]
  return(ident)
}

# test
ident_test <- fun(colnames_datExpr=colnames(datExpr_test), ident=ident_test)
# ref
ident_ref <- fun(colnames_datExpr = colnames(datExpr_ref), ident=ident_ref)
# list_iterable=list("colnames_datExpr"=lapply(list_datExpr_ref, colnames),"ident"=list_ident_ref)
# list_ident_ref <- safeParallel(fun=fun, list_iterable=list_iterable)

######################################################################
########## MAKE CELL CLUSTER SCE OBJECTS, SELECT FEATURES ############
######################################################################

message("Subsetting test data")

list_datExpr_test_sub <- lapply(names(table(ident_test)), function(id){
  datExpr_test[,ident_test==id]
  #SubsetData(datExpr_test, ident.use = id, do.clean = T)
})

names(list_datExpr_test_sub) <- names(table(ident_test))

message("Making single cell experiment objects and selecting features")

outfile = paste0(dir_log, prefix_test, "_", prefix_run, "_subset_make_sce.txt")

#For each subset, try a range of values for n features (genes) to use in aligning with clusters in the reference dataset

fun <- function(datExpr_test_sub, vec_n_features) {
  sce_test_sub <- SingleCellExperiment(assays=list(logcounts = as.matrix(datExpr_test_sub)),
                                  colData = list(colnames = colnames(datExpr_test_sub)),
                                  rowData = list(rownames = rownames(datExpr_test_sub)))
  rm(datExpr_test_sub)

  rowData(sce_test_sub)$feature_symbol <- rowData(sce_test_sub)$rownames

  list_sce_test_sub = lapply(vec_n_features, function(n_features) selectFeatures(sce_test_sub,
                                                               suppress_plot=T,
                                                               n_features = n_features))
  for (i in 1:length(vec_n_features)) {
    name <- paste0("scmap_features", i)
    rowData(sce_test_sub)[name] <- rowData(list_sce_test_sub[[i]])$scmap_features
  }

  rm(list_sce_test_sub)
  return(sce_test_sub)
}

list_iterable = list("X"=list_datExpr_test_sub)

list_sce_test_sub <- safeParallel(fun=fun, list_iterable=list_iterable, outfile=outfile, vec_n_features=vec_n_features)

rm(list_datExpr_test_sub)

######################################################################
######## MAKE REF SET SCE OBJECT AND SELECT ALIGNMENT FEATURES #######
######################################################################

outfile = paste0(dir_log, prefix_run, "_", prefix_test, "_make_sce_ref_findfeatures.txt")

fun <- function(datExpr_ref, ident_ref) {

  #rownames(datExpr_ref) <- NULL #what is this for? Seems we need it immediately below

  sce_ref <- SingleCellExperiment(assays=list(logcounts = as.matrix(datExpr_ref)),
                       colData = list(colNames = as.character(colnames(datExpr_ref))),
                       rowData = list(rownames = as.character(rownames(datExpr_ref))))
  #Add column (cellwise) metadata
  colData(sce_ref)$cell_type1 <- gsub(" ", "_", ident_ref)

  #Add rowdata - feature_symbol is needed for alignment
  rowData(sce_ref)$feature_symbol <- rowData(sce_ref)$rownames

  message("Selecting reference set features")

  fun1  = function(n_features) selectFeatures(sce_ref,
                                             suppress_plot=T,
                                             n_features = n_features)

  list_sce_ref_tmp <- lapply("X"=vec_n_features, FUN = fun1)

  for (i in 1:length(vec_n_features)) {
    name <- paste0("scmap_features", i)
    rowData(sce_ref)[name] <- rowData(list_sce_ref_tmp[[i]])$scmap_features
  }

  rm(list_sce_ref_tmp)

  return(sce_ref)

}

# list_iterable=list(ident_ref = list_ident_ref,
#           datExpr_ref = list_datExpr_ref)

sce_ref <- fun(ident_ref=ident_ref, datExpr_ref=datExpr_ref)

#safeParallel(fun=fun,list_iterable=list_iterable, outfile=outfile, vec_n_features=vec_n_features)

######################################################################
###################### FIND SCMAP CLUSTER INDEX ######################
#####################################################################
# The scmap-cluster index of a reference dataset is created by finding the median gene expression for each cluster.
# By default scmap uses the cell_type1 column of the colData slot in the reference to identify clusters.
# (Other columns can be manually selected by adjusting indexCluster's cluster_col parameter)

message(paste0(c("Computing scmap cluster index for", vec_n_features, "features"), collapse = " "))

fun <- function(i) {
  index <- paste0("scmap_features",i)
  rowData(sce_ref)$scmap_features <- rowData(sce_ref)[[index]]
  return(indexCluster(sce_ref))
}

list_iterable = list("X"=1:length(vec_n_features))
# test sets
#list_list_sce_test_sub <- lapply("X"=list_sce_test_sub, FUN= function(sce_test_sub) safeParallel(fun=fun, list_iterable=list_iterable, sce_obj=sce_test_sub))

# ref sets
# list_list_sce_ref <- list()
# for (name in names(list_sce_ref)) {
#sce_obj=list_sce_ref[[name]]
list_sce_ref <- safeParallel(fun=fun, list_iterable=list_iterable)
names(list_sce_ref) <- vec_n_features


# Why does the vectorised version not work??
#list_list_sce_ref <- lapply("X"=list_sce_ref, FUN = function(sce_obj) safeParallel(fun=fun, list_iterable=list_iterable, sce_obj=sce_obj))

# For each n_features,plot heatmaps of cluster index
message(paste0(c("Plotting scmap cluster index heatmaps for", vec_n_features, "features"), collapse = " "))

try({
  par(mfrow=c(1, length(vec_n_features)))
  pdf(file = paste0(dir_plots, prefix_test, "_", prefix_run, "_sce_ref_cluster_index.pdf"))
  for (sce_ref in list_sce_ref) {
    heatmap(as.matrix(metadata(sce_ref)$scmap_cluster_index))
  }
  dev.off()
})
######################################################################
############################ PROJECTION ##############################
######################################################################
# Once the scmap-cluster index has been generated we can use it to project our dataset to another.
# This can be done with one index at a time, but scmap also allows for simultaneous projection to multiple indexes if they are provided as a list:

message("scmapCluster projection")
outfile = paste0(dir_log, prefix_test, "_", prefix_run, "_scmap_align.txt")

listClust_listnFeat_scmapOuts <- list()

#for (refset in names(list_list_sce_ref)){  # loop over reference set names
# list_sce_ref <- list_list_sce_ref[[refset]]
  for (cluster in names(list_sce_test_sub)) { # loop over test sets, i.e. clusters
    sce_test_sub <- list_sce_test_sub[[cluster]]
    fun = function(i) {
      feat <- paste0("scmap_features", i)
      # for hs.liver, select the set of features stores in rowData and set it in position to be found by scmapCluster
      rowData(sce_test_sub)$scmap_features <- rowData(sce_test_sub)[[feat]]

      result <- tryCatch({scmapCluster(projection = sce_test_sub,
                             index_list = list(ref = metadata(list_sce_ref[[i]])$scmap_cluster_index),
                             threshold=threshold)},
                         error = function(err) {
                           message(paste0(ref_set, ": ", cluster, ": scmapCluster failed with error: ", err))
                           return(NA_character_)
                         })
      return(result)
    }
    list_iterable = list("X" = 1:length(vec_n_features))
    listClust_listnFeat_scmapOuts[[cluster]] <- safeParallel(fun=fun, 
                                                             list_iterable=list_iterable, 
                                                             #sce_test_sub=sce_test_sub, 
                                                             #list_sce_ref=list_sce_ref, 
                                                             outfile=outfile)
    names(listClust_listnFeat_scmapOuts[[cluster]]) <- vec_n_features
  }
#}

try(warnings())

saveMeta(savefnc=saveRDS, object = listClust_listnFeat_scmapOuts, file=paste0(dir_RObjects, prefix_run, "_", prefix_test, "_full_scmap_results.RDS.gz"), compress="gzip")
saveMeta(savefnc=saveRDS, object=list_sce_test_sub, file=paste0(dir_RObjects, prefix_run, "_", prefix_test, "_list_sce_test_sub.RDS.gz"), compress="gzip")
######################################################################
########## GET ALIGN STATS FOR EACH TEST CLUSTER, N_FEATURES #########
######################################################################
# scmap-cluster projects the query dataset to all projections defined in the index_list.
# The results of cell label assignements are merged into one matrix
message("Preparing alignment summary dataframe")

#list_results_df <- list()
#k=1
#for (refset in names(list_list_sce_ref)) {
  # one dataframe per reference set
 # list_results_df[[refset]] <- if (k==1) {
  results_df <-
    data.frame(test_cluster = rep(names(table(ident_test)),
                 times = rep(length(vec_n_features),
                             times= length(names(table(ident_test))))),
                 n_features = vec_n_features, # recycles
    prop.assign = 0,
    ref_dataset = NA_character_,
    ref_celltype_1 = NA_character_,
    n_cells_assigned_1 = NA_integer_,
    ref_celltype_2 = NA_character_,
    n_cells_assigned_2 = NA_integer_,
    ref_celltype_3 = NA_character_,
    n_cells_assigned_3 = NA_integer_,
    ref_celltype_4 = NA_character_,
    n_cells_assigned_4 = NA_integer_,
    ref_celltype_5 = NA_character_,
    n_cells_assigned_5 = NA_integer_,
    remainder=NA_integer_,
    stringsAsFactors = F)
  # } else {
  #   data.frame(prop.assign = rep(0, nrow(list_results_df[[1]])),
  #              ref_dataset = NA_character_,
  #              ref_celltype_1 = NA_character_,
  #              n_cells_assigned_1 = NA_integer_,
  #              ref_celltype_2 = NA_character_,
  #              n_cells_assigned_2 = NA_integer_,
  #              ref_celltype_3 = NA_character_,
  #              n_cells_assigned_3 = NA_integer_,
  #              ref_celltype_4 = NA_character_,
  #              n_cells_assigned_4 = NA_integer_,
  #              ref_celltype_5 = NA_character_,
  #              n_cells_assigned_5 = NA_integer_,
  #              stringsAsFactors = F)
  # }
  refset <- names(path_datExpr_ref)
  for (i in 1:length(listClust_listnFeat_scmapOuts)) { # loop over cell clusters in test set to align
    for (j  in 1:length(listClust_listnFeat_scmapOuts[[i]])) { # loop over n features
      idx_row <- (i-1)*length(vec_n_features)+j
      results_df[[idx_row,'prop.assign']] <- if(!is.null(listClust_listnFeat_scmapOuts[[i]][[j]]$scmap_cluster_labs)) 1-sum(listClust_listnFeat_scmapOuts[[i]][[j]]$scmap_cluster_labs[,'ref']=="unassigned")/nrow(listClust_listnFeat_scmapOuts[[i]][[j]]$scmap_cluster_labs) else NA_real_
      results_df[[idx_row,'ref_dataset']] <- refset
      if (!is.null(listClust_listnFeat_scmapOuts[[i]][[j]]$scmap_cluster_labs)) {
        results_df[[idx_row,'ref_celltype_1']] <- if (!is.null(listClust_listnFeat_scmapOuts[[i]][[j]]$scmap_cluster_labs[,'ref'])) {
          tryCatch({listClust_listnFeat_scmapOuts[[i]][[j]]$scmap_cluster_labs[,'ref'] %>% table %>% sort(.,decreasing = T) %>% names %>% '[['(1)
            }, error = function(err){NA_character_})
          } else {NA_character_}
        results_df[[idx_row,'n_cells_assigned_1']] <- if(!is.null(listClust_listnFeat_scmapOuts[[i]][[j]]$scmap_cluster_labs[,'ref'])) {
          tryCatch({listClust_listnFeat_scmapOuts[[i]][[j]]$scmap_cluster_labs[,'ref'] %>% table %>% sort(.,decreasing = T) %>% '[['(1)
            }, error = function(err){NA_integer_})
        } else {NA_integer_}
        results_df[[idx_row,'ref_celltype_2']] <- if (!is.null(listClust_listnFeat_scmapOuts[[i]][[j]]$scmap_cluster_labs[,'ref'])) {
          tryCatch({listClust_listnFeat_scmapOuts[[i]][[j]]$scmap_cluster_labs[,'ref'] %>% table %>% sort(.,decreasing = T) %>% names %>% '[['(2)
            }, error = function(err){NA_character_})
        } else {NA_character_}
        results_df[[idx_row,'n_cells_assigned_2']] <- if (!is.null(listClust_listnFeat_scmapOuts[[i]][[j]]$scmap_cluster_labs[,'ref'])) {
          tryCatch({listClust_listnFeat_scmapOuts[[i]][[j]]$scmap_cluster_labs[,'ref'] %>% table %>% sort(.,decreasing = T) %>% '[['(2)
            }, error = function(err){NA_integer_})
        } else {NA_integer_}
        results_df[[idx_row,'ref_celltype_3']] <- if (!is.null(listClust_listnFeat_scmapOuts[[i]][[j]]$scmap_cluster_labs[,'ref'])) {
          tryCatch({
            listClust_listnFeat_scmapOuts[[i]][[j]]$scmap_cluster_labs[,'ref'] %>% table %>% sort(.,decreasing = T) %>% names %>% '[['(3)
            }, error = function(err){NA_character_})
        } else {NA_character_}
        results_df[[idx_row,'n_cells_assigned_3']] <- if(!is.null(listClust_listnFeat_scmapOuts[[i]][[j]]$scmap_cluster_labs[,'ref'])) {
          tryCatch({listClust_listnFeat_scmapOuts[[i]][[j]]$scmap_cluster_labs[,'ref']  %>% table %>% sort(.,decreasing = T) %>% '[['(3)
            }, error = function(err){NA_integer_})
        } else {NA_integer_}
        results_df[[idx_row,'ref_celltype_4']] <- if (!is.null(listClust_listnFeat_scmapOuts[[i]][[j]]$scmap_cluster_labs[,'ref'])) {
          tryCatch({
            listClust_listnFeat_scmapOuts[[i]][[j]]$scmap_cluster_labs[,'ref'] %>% table %>% sort(.,decreasing = T) %>% names %>% '[['(4)
          }, error = function(err){NA_character_})
        } else {NA_character_}
        results_df[[idx_row,'n_cells_assigned_4']] <- if(!is.null(listClust_listnFeat_scmapOuts[[i]][[j]]$scmap_cluster_labs[,'ref'])) {
          tryCatch({listClust_listnFeat_scmapOuts[[i]][[j]]$scmap_cluster_labs[,'ref']  %>% table %>% sort(.,decreasing = T) %>% '[['(4)
          }, error = function(err){NA_integer_})
        } else {NA_integer_}
        results_df[[idx_row,'ref_celltype_5']] <- if (!is.null(listClust_listnFeat_scmapOuts[[i]][[j]]$scmap_cluster_labs[,'ref'])) {
          tryCatch({
            listClust_listnFeat_scmapOuts[[i]][[j]]$scmap_cluster_labs[,'ref'] %>% table %>% sort(.,decreasing = T) %>% names %>% '[['(5)
          }, error = function(err){NA_character_})
        } else {NA_character_}
        results_df[[idx_row,'n_cells_assigned_5']] <- if(!is.null(listClust_listnFeat_scmapOuts[[i]][[j]]$scmap_cluster_labs[,'ref'])) {
          tryCatch({listClust_listnFeat_scmapOuts[[i]][[j]]$scmap_cluster_labs[,'ref']  %>% table %>% sort(.,decreasing = T) %>% '[['(5)
          }, error = function(err){NA_integer_})
        } else {NA_integer_}
        results_df[[idx_row,'remainder']] <- if(!is.null(listClust_listnFeat_scmapOuts[[i]][[j]]$scmap_cluster_labs[,'ref'])) {
          tryCatch({
            length(listClust_listnFeat_scmapOuts[[i]][[j]]$scmap_cluster_labs[,'ref']) -sum(table(listClust_listnFeat_scmapOuts[[i]][[j]]$scmap_cluster_labs[,'ref']) %>% sort(., decreasing=T) %>% '['(1:5))
          }, error = function(err){NA_integer_})
        } else {NA_integer_}
      }
    }
  }
  #k=k+1
#}

#results_df <- if (length(list_results_df)==1) list_results_df[[1]] else Reduce(f=cbind, x=list_results_df)

# write to disk
saveMeta(savefnc=write.csv, x=results_df, file=sprintf("%s%s_%s_scmap_results_df.csv", dir_tables, prefix_test, prefix_run), row.names=F)

######################################################################
############################# PLOTS ##################################
######################################################################

message("Plotting")
# Corresponding similarities are stored in the scmap_cluster_siml item:

# for (refset in names(listRef_listClust_listnFeat_scmapOuts)) {
#   for (i in 1:length(listRef_listClust_listnFeat_scmapOuts[[refset]])) {
#     pdf(sprintf("%s%s_%s_scmapCluster_similarities_hist_cluster_%s.pdf", dir_plots, prefix_test, refset, names(table(as.character(metadata_test[[metadata_test_id_col]])))[i]),h=10,w=16)
#     for (j  in 1:length(listRef_listClust_listnFeat_scmapOuts[[refset]][[i]])) {
#        hist(listRef_listClust_listnFeat_scmapOuts[[refset]][[i]][[j]]$scmap_cluster_siml[,'ref'],breaks = 100, main = paste0("test - ", refset, "  similarities, ",vec_n_features[j], " features"))
#     }
#     dev.off()
#   }
# }

# The results of scmap-cluster can be visualized as a Sankey diagram to show how cell-clusters are matched (getSankey() function).
# Note that the Sankey diagram will only be informative if both the query and the reference datasets have been clustered,
# but it is not necessary to have meaningful labels assigned to the query (cluster1, cluster2 etc. is sufficient)

# for (name in names(list_sce_test_sub)) {
#   colData(list_sce_test_sub[[i]])$cell_type1 <- name
# }

# for (i in 1:length(listClust_listnFeat_scmapOuts)) {
#   for (j  in 1:length(listClust_listnFeat_scmapOuts[[i]])) {
#     pdf(sprintf("%s%s_%s_scmapCluster_Sankey_cluster_%s_n_features_%s.pdf", dir_plots, prefix_test, prefix_run, names(table(as.character(metadata_test[[metadata_test_id_col]])))[i], vec_n_features[j]),h=10,w=16)
#     # TODO: need to address error: GDK_BACKEND does not match available displays
#     # https://askubuntu.com/questions/359753/gtk-warning-locale-not-supported-by-c-library-when-starting-apps-from-th
#     p <- getSankey(reference = colData(list_sce_test_sub[[j]])$cell_type1,
#                 clusters = listClust_listnFeat_scmapOuts[[i]][[j]]$scmap_cluster_labs[,'ref'],
#                 plot_height = 400)

#     #plot(getSankey(ann[,1], ann[,1]))
#     #dev.off()
#   }
# }

######################################################################
#################### FOR EACH CLUSTER GET BEST BET ###################
######################################################################
#
# if (FALSE) {
#   sapply(sort(unique(results_df[["test_cluster"]])), function(cluster) {
#     results_df_clust <- results_df[results_df[["test_cluster"]]==cluster,]
#     idx_col_prop_assign <- seq(from=3, by=3, to=ncol(results_df_clust))
#     idx_idx_rowmax_col <- apply(X = results_df_clust[,idx_col_prop_assign, drop=F], MARGIN=1, FUN = which.max)
#     idx_col_max <- sapply(idx_idx_rowmax_col, function(idx) idx_col_prop_assign[idx])
#     idx_row_max <- which.max(sapply(idx_col_max, function(idx) results_df_clust[,idx, drop=F], simplify = T))
#     results_df_clust[idx_row_max, c(idx_col_max[idx_row_max]:(idx_col_max[idx_row_max]+2))]
#   }, simplify = T) %>% t -> bestFits
#
#   # write to disk
#
#   write.csv(bestFits, file = paste0(dir_tables, prefix_test, "_", prefix_run, "_scmap_bestfits.csv"), quote = F,row.names = T)
# }
######################################################################
############### PLOT SCMAP CLUSTER ANNOTATION ON TSNE ################
######################################################################

message("Plotting new cluster assignment in test dataset")

message("loading dataset")

datExpr_test <- load_obj(path_datExpr_test)

if (!"Seurat" %in% class(datExpr_test)) {
  datExpr_test <- CreateSeuratObject(counts=datExpr_test,
  project=names(path_datExpr_test))
}

if (is.null(datExpr_test@reductions$pca)) {
  message("Running PCA")
  datExpr_test <- RunPCA(datExpr_test,
                         npcs = 30, 
                         seed.use = randomSeed)
}
if (is.null(datExpr_test@reductions$tsne)) {
  message("Running UMAP")
  datExpr_test <- RunUMAP(datExpr_test, 
                          dims=1:30,
                          seed.use= randomSeed)
}
  #datExpr_test[[metadata_test_id_col]] else
  #Idents(datExpr_test) <- ident_ref datExpr_test[[metadata_test_id_col]]
scmap_annot_labels <- paste0("(",results_df[["test_cluster"]], ") ", results_df[["ref_celltype_1"]], ":", results_df[["n_cells_assigned_1"]], ", ", results_df[["ref_celltype_2"]], ":", results_df[["n_cells_assigned_2"]], ", ",results_df[["ref_celltype_3"]], ":", results_df[["n_cells_assigned_3"]])
names(scmap_annot_labels) <- names(listClust_listnFeat_scmapOuts)
datExpr_test$scmap_annot <- scmap_annot_labels[match(ident_test,  names(scmap_annot_labels))]

#Idents(datExpr_test) <- datExpr_test$scmap_annot
DimPlot(object = datExpr_test,
        reduction="tsne",
        group.by = "scmap_annot",
        label.size= 6 ,
        label=T)

message("Saving test dataset with new annotations added to metadata")

saveMeta(savefnc=ggsave, filename = paste0(dir_plots, prefix_test, "_", prefix_run, "_tSNE_scmap_annot.pdf"), width = 40, height=30)

saveMeta(savefnc=saveRDS, object=datExpr_test, file=gsub("\\.RDS\\..*", paste0("_", prefix_run, ".RDS.gz"), path_datExpr_test), compress = "gzip")

######################################################################
############################ WRAP UP #################################
######################################################################

message("Script done!")
