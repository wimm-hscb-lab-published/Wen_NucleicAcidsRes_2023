# Set working directory
#setwd("/Users/seanwen/Documents/MARVEL/Github_WIMM/Wen_NucleicAcidsRes_2023/scripts/scripts_02_prepare_marvel_input/")

# Load packages
library(SingCellaR)
library(data.table)
library(plyr)
library(Matrix)
library(SingleCellExperiment)
    
# Load integrated R object
path <- "../../data/SingCellaR_output/integrated/"
file <- "Integrated_iPSC_CardioDay2_4_10.rdata"
load(file=paste(path, file, sep=""))

# Retrive metadata
df.pheno <- object@meta.data

######################################################################
####################### SRR9008754 (IPSC) ############################
######################################################################

# Read files
path <- "../../data/STARsolo/SRR9008754/Solo_out/Gene/"
folder <- "filtered/"
data_matrices_dir <- paste(path, folder, sep="")
object.2 <- new("SingCellaR")
object.2@dir_path_10x_matrix <- data_matrices_dir
object.2@sample_uniq_id <- ""
load_matrices_from_cellranger(object.2, cellranger.version=3)

# Retrieve cell ids
cell.ids <- colnames(object.2)
df.pheno.2 <- data.frame("cell.id."=cell.ids,
                         "cell.id"=cell.ids,
                         stringsAsFactors=FALSE
                         )

# Annotate donor id
df.pheno.2$cell.id <- paste(df.pheno.2$cell.id, "SRR9008754", sep="")

# Subset overlapping donor cells
overlap <- intersect(df.pheno$Cell, df.pheno.2$cell.id)
print(paste(length(df.pheno$Cell), " found by SingCellaR", sep=""))
print(paste(length(overlap), " found by STARsolo", sep=""))

df.pheno.2.small <- df.pheno.2[which(df.pheno.2$cell.id %in% overlap), ]
    
# Match matrix to phenoData
    # Retrieve matrix
    df <- counts(object.2)
    
    # Subset overlapping cells
    df.small <- df[, df.pheno.2.small$cell.id.]
    
    # Update cell ids
    colnames(df.small) <- df.pheno.2.small$cell.id
    
# Save as new object
df.SRR9008754 <- df.small

######################################################################
################## SRR9008755 (CARDIO DAY 2) #########################
######################################################################

# Read files
path <- "../../data/STARsolo/SRR9008755/Solo_out/Gene/"
folder <- "filtered/"
data_matrices_dir <- paste(path, folder, sep="")
object.2 <- new("SingCellaR")
object.2@dir_path_10x_matrix <- data_matrices_dir
object.2@sample_uniq_id <- ""
load_matrices_from_cellranger(object.2, cellranger.version=3)

# Retrieve cell ids
cell.ids <- colnames(object.2)
df.pheno.2 <- data.frame("cell.id."=cell.ids,
                         "cell.id"=cell.ids,
                         stringsAsFactors=FALSE
                         )

# Annotate donor id
df.pheno.2$cell.id <- paste(df.pheno.2$cell.id, "SRR9008755", sep="")

# Subset overlapping donor cells
overlap <- intersect(df.pheno$Cell, df.pheno.2$cell.id)
print(paste(length(df.pheno$Cell), " found by SingCellaR", sep=""))
print(paste(length(overlap), " found by STARsolo", sep=""))

df.pheno.2.small <- df.pheno.2[which(df.pheno.2$cell.id %in% overlap), ]
    
# Match matrix to phenoData
    # Retrieve matrix
    df <- counts(object.2)
    
    # Subset overlapping cells
    df.small <- df[, df.pheno.2.small$cell.id.]
    
    # Update cell ids
    colnames(df.small) <- df.pheno.2.small$cell.id
    
# Save as new object
df.SRR9008755 <- df.small

######################################################################
################## SRR9008752 (CARDIO DAY 4) #########################
######################################################################

# Read files
path <- "../../data/STARsolo/SRR9008752/Solo_out/Gene/"
folder <- "filtered/"
data_matrices_dir <- paste(path, folder, sep="")
object.2 <- new("SingCellaR")
object.2@dir_path_10x_matrix <- data_matrices_dir
object.2@sample_uniq_id <- ""
load_matrices_from_cellranger(object.2, cellranger.version=3)

# Retrieve cell ids
cell.ids <- colnames(object.2)
df.pheno.2 <- data.frame("cell.id."=cell.ids,
                         "cell.id"=cell.ids,
                         stringsAsFactors=FALSE
                         )

# Annotate donor id
df.pheno.2$cell.id <- paste(df.pheno.2$cell.id, "SRR9008752", sep="")

# Subset overlapping donor cells
overlap <- intersect(df.pheno$Cell, df.pheno.2$cell.id)
print(paste(length(df.pheno$Cell), " found by SingCellaR", sep=""))
print(paste(length(overlap), " found by STARsolo", sep=""))

df.pheno.2.small <- df.pheno.2[which(df.pheno.2$cell.id %in% overlap), ]
    
# Match matrix to phenoData
    # Retrieve matrix
    df <- counts(object.2)
    
    # Subset overlapping cells
    df.small <- df[, df.pheno.2.small$cell.id.]
    
    # Update cell ids
    colnames(df.small) <- df.pheno.2.small$cell.id
    
# Save as new object
df.SRR9008752 <- df.small

######################################################################
################# SRR9008753 (CARDIO DAY 10) #########################
######################################################################

# Read files
path <- "../../data/STARsolo/SRR9008753/Solo_out/Gene/"
folder <- "filtered/"
data_matrices_dir <- paste(path, folder, sep="")
object.2 <- new("SingCellaR")
object.2@dir_path_10x_matrix <- data_matrices_dir
object.2@sample_uniq_id <- ""
load_matrices_from_cellranger(object.2, cellranger.version=3)

# Retrieve cell ids
cell.ids <- colnames(object.2)
df.pheno.2 <- data.frame("cell.id."=cell.ids,
                         "cell.id"=cell.ids,
                         stringsAsFactors=FALSE
                         )

# Annotate donor id
df.pheno.2$cell.id <- paste(df.pheno.2$cell.id, "SRR9008753", sep="")

# Subset overlapping donor cells
overlap <- intersect(df.pheno$Cell, df.pheno.2$cell.id)
print(paste(length(df.pheno$Cell), " found by SingCellaR", sep=""))
print(paste(length(overlap), " found by STARsolo", sep=""))

df.pheno.2.small <- df.pheno.2[which(df.pheno.2$cell.id %in% overlap), ]
    
# Match matrix to phenoData
    # Retrieve matrix
    df <- counts(object.2)
    
    # Subset overlapping cells
    df.small <- df[, df.pheno.2.small$cell.id.]
    
    # Update cell ids
    colnames(df.small) <- df.pheno.2.small$cell.id
    
# Save as new object
df.SRR9008753 <- df.small

######################################################################
############################# MERGE ##################################
######################################################################

# Merge
df.merged <- cbind(df.SRR9008754, df.SRR9008755, df.SRR9008752, df.SRR9008753)
dim(df.merged); class(df.merged)

# Save file
    # Matrix
    path <- "../../data/MARVEL_input/Gene_STARsolo/"
    file <- "matrix_counts.mtx"
    writeMM(df.merged, file=paste(path, file, sep=""))
    
    # phenoData
    df.pheno.merged <- data.frame("cell.id"=colnames(df.merged))
    
    path <- "../../data/MARVEL_input/Gene_STARsolo/"
    file <- "phenoData.txt"
    write.table(df.pheno.merged, paste(path, file, sep=""), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
    
    # featureData
    df.feature.merged <- data.frame("gene_short_name"=rownames(df.merged))
    
    path <- "../../data/MARVEL_input/Gene_STARsolo/"
    file <- "featureData.txt"
    write.table(df.feature.merged, paste(path, file, sep=""), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
    
    
    
