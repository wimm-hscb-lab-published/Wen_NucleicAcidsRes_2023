# Set working directory
#setwd("/Users/seanwen/Documents/MARVEL/Github_WIMM/Wen_NucleicAcidsRes_2023/scripts/scripts_02_prepare_marvel_input/")

# Load packages
library(SingCellaR)
library(plyr)
library(Matrix)

# Load integrated R object
path <- "../../data/SingCellaR_output/integrated/"
file <- "Integrated_iPSC_CardioDay2_4_10.rdata"
load(file=paste(path, file, sep=""))

#########################################################################

# Prepare phenoData
    # Merge sample,  cluster, metadata
    df.pheno <- object@meta.data

    # Rename cell, donor id column
    names(df.pheno)[which(names(df.pheno)=="Cell")] <- "cell.id"
    names(df.pheno)[which(names(df.pheno)=="sampleID")] <- "donor.id"
    
    # Indicate cell type
    df.pheno$cell.type <- NA
    df.pheno$cell.type[grep("SRR9008754", df.pheno$donor.id)] <- "iPSC"
    df.pheno$cell.type[grep("SRR9008755", df.pheno$donor.id)] <- "Cardio day 2"
    df.pheno$cell.type[grep("SRR9008752", df.pheno$donor.id)] <- "Cardio day 4"
    df.pheno$cell.type[grep("SRR9008753", df.pheno$donor.id)] <- "Cardio day 10"
    table(df.pheno$cell.type, df.pheno$donor.id)
    
    # Subset relevant columns
    cols <- c("cell.id", "donor.id", "cell.type")
    df.pheno <- df.pheno[, cols]
    
    # Save file
    path <- "../../data/MARVEL_input/Gene_SingCellaR/"
    file <- "phenoData.txt"
    write.table(df.pheno, paste(path, file, sep=""), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

# Prepare gene metadata
    # Retrieve
    df.feature <- get_genes_metadata(object)
    
    # Create column for gene name
    df.feature$gene_short_name <- row.names(df.feature)
    df.feature <- df.feature[, "gene_short_name", drop=FALSE]
    row.names(df.feature) <- NULL
    
    # Save file
    path <- "../../data/MARVEL_input/Gene_SingCellaR/"
    file <- "featureData.txt"
    write.table(df.feature, paste(path, file, sep=""), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
    
# Prepare matrix
    # Retrieve matrix
    df <- get_normalized_umi(object)

    # Check alignmnet
    table(colnames(df)==df.pheno$cell.id)
    table(rownames(df)==df.feature$gene_short_name)

    # Save file
    path <- "../../data/MARVEL_input/Gene_SingCellaR/"
    file <- "matrix_normalised.mtx"
    writeMM(df, file=paste(path, file, sep=""))
