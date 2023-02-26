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

# Read barcode file
    # Read file
    path <- "../../data/STARsolo/SRR9008754/Solo_out/SJ/raw/"
    file <- "barcodes.tsv"
    df.pheno.2 <- as.data.frame(fread(paste(path, file, sep=""), sep="\t", header=FALSE, stringsAsFactors=FALSE))
    
    # Retrieve cell ids
    names(df.pheno.2) <- "cell.id."
    df.pheno.2$cell.id <- df.pheno.2$cell.id.

    # Indicate donor id
    df.pheno.2$cell.id <- paste(df.pheno.2$cell.id, "_SRR9008754", sep="")

# Read feature file
    # Read file
    path <- "../../data/STARsolo/SRR9008754/Solo_out/SJ/raw/"
    file <- "features.tsv"
    df.feature <- as.data.frame(fread(paste(path, file, sep=""), sep="\t", header=FALSE, stringsAsFactors=FALSE))
    
    # Subset relevant columns
    cols <- c("V1", "V2", "V3")
    df.feature <- df.feature[, cols]
    names(df.feature) <- c("chr", "start", "end")
    df.feature$chr <- paste("chr", df.feature$chr, sep="")
    
    # Create unique junction id
    . <- data.frame("coord.intron"=paste(df.feature$chr, df.feature$start, df.feature$end, sep=":"))
    df.feature <- cbind.data.frame(., df.feature)
    
# Read matrix
    # Read file
    path <- "../../data/STARsolo/SRR9008754/Solo_out/SJ/raw/"
    file <- "matrix.mtx"
    df <- readMM(paste(path, file, sep=""))
    
    # Check alignment
    ncol(df)==nrow(df.pheno.2)
    nrow(df)==nrow(df.feature)

    # Annotate columns
    colnames(df) <- df.pheno.2$cell.id.
    rownames(df) <- df.feature$coord.intron

# Subset overlapping cells
overlap <- intersect(df.pheno$Cell, df.pheno.2$cell.id)
print(paste(length(df.pheno$Cell), " found by SingCellaR", sep=""))
print(paste(length(overlap), " found by STARsolo", sep=""))

df.pheno.2.small <- df.pheno.2[which(df.pheno.2$cell.id %in% overlap), , drop=FALSE]

# Match matrix to phenoData
    # Subset overlapping cells
    df.small <- df[, df.pheno.2.small$cell.id.]
    
    # Update cell ids
    colnames(df.small) <- df.pheno.2.small$cell.id

# Save as new object
df.SRR9008754 <- df.small
df.pheno.SRR9008754 <- df.pheno.2.small
df.feature.SRR9008754 <- df.feature

# Check alignment
table(colnames(df.SRR9008754)==df.pheno.SRR9008754$cell.id)
table(rownames(df.SRR9008754)==df.feature.SRR9008754$coord.intron)

######################################################################
################## SRR9008755 (CARDIO DAY 2) #########################
######################################################################

# Read barcode file
    # Read file
    path <- "../../data/STARsolo/SRR9008755/Solo_out/SJ/raw/"
    file <- "barcodes.tsv"
    df.pheno.2 <- as.data.frame(fread(paste(path, file, sep=""), sep="\t", header=FALSE, stringsAsFactors=FALSE))
    
    # Retrieve cell ids
    names(df.pheno.2) <- "cell.id."
    df.pheno.2$cell.id <- df.pheno.2$cell.id.

    # Indicate donor id
    df.pheno.2$cell.id <- paste(df.pheno.2$cell.id, "_SRR9008755", sep="")

# Read feature file
    # Read file
    path <- "../../data/STARsolo/SRR9008755/Solo_out/SJ/raw/"
    file <- "features.tsv"
    df.feature <- as.data.frame(fread(paste(path, file, sep=""), sep="\t", header=FALSE, stringsAsFactors=FALSE))
    
    # Subset relevant columns
    cols <- c("V1", "V2", "V3")
    df.feature <- df.feature[, cols]
    names(df.feature) <- c("chr", "start", "end")
    df.feature$chr <- paste("chr", df.feature$chr, sep="")
    
    # Create unique junction id
    . <- data.frame("coord.intron"=paste(df.feature$chr, df.feature$start, df.feature$end, sep=":"))
    df.feature <- cbind.data.frame(., df.feature)
    
# Read matrix
    # Read file
    path <- "../../data/STARsolo/SRR9008755/Solo_out/SJ/raw/"
    file <- "matrix.mtx"
    df <- readMM(paste(path, file, sep=""))
    
    # Check alignment
    ncol(df)==nrow(df.pheno.2)
    nrow(df)==nrow(df.feature)

    # Annotate columns
    colnames(df) <- df.pheno.2$cell.id.
    rownames(df) <- df.feature$coord.intron

# Subset overlapping cells
overlap <- intersect(df.pheno$Cell, df.pheno.2$cell.id)
print(paste(length(df.pheno$Cell), " found by SingCellaR", sep=""))
print(paste(length(overlap), " found by STARsolo", sep=""))

df.pheno.2.small <- df.pheno.2[which(df.pheno.2$cell.id %in% overlap), , drop=FALSE]

# Match matrix to phenoData
    # Subset overlapping cells
    df.small <- df[, df.pheno.2.small$cell.id.]
    
    # Update cell ids
    colnames(df.small) <- df.pheno.2.small$cell.id

# Save as new object
df.SRR9008755 <- df.small
df.pheno.SRR9008755 <- df.pheno.2.small
df.feature.SRR9008755 <- df.feature

# Check alignment
table(colnames(df.SRR9008755)==df.pheno.SRR9008755$cell.id)
table(rownames(df.SRR9008755)==df.feature.SRR9008755$coord.intron)

######################################################################
################## SRR9008752 (CARDIO DAY 4) #########################
######################################################################

# Read barcode file
    # Read file
    path <- "../../data/STARsolo/SRR9008752/Solo_out/SJ/raw/"
    file <- "barcodes.tsv"
    df.pheno.2 <- as.data.frame(fread(paste(path, file, sep=""), sep="\t", header=FALSE, stringsAsFactors=FALSE))
    
    # Retrieve cell ids
    names(df.pheno.2) <- "cell.id."
    df.pheno.2$cell.id <- df.pheno.2$cell.id.

    # Indicate donor id
    df.pheno.2$cell.id <- paste(df.pheno.2$cell.id, "_SRR9008752", sep="")

# Read feature file
    # Read file
    path <- "../../data/STARsolo/SRR9008752/Solo_out/SJ/raw/"
    file <- "features.tsv"
    df.feature <- as.data.frame(fread(paste(path, file, sep=""), sep="\t", header=FALSE, stringsAsFactors=FALSE))
    
    # Subset relevant columns
    cols <- c("V1", "V2", "V3")
    df.feature <- df.feature[, cols]
    names(df.feature) <- c("chr", "start", "end")
    df.feature$chr <- paste("chr", df.feature$chr, sep="")
    
    # Create unique junction id
    . <- data.frame("coord.intron"=paste(df.feature$chr, df.feature$start, df.feature$end, sep=":"))
    df.feature <- cbind.data.frame(., df.feature)
    
# Read matrix
    # Read file
    path <- "../../data/STARsolo/SRR9008752/Solo_out/SJ/raw/"
    file <- "matrix.mtx"
    df <- readMM(paste(path, file, sep=""))
    
    # Check alignment
    ncol(df)==nrow(df.pheno.2)
    nrow(df)==nrow(df.feature)

    # Annotate columns
    colnames(df) <- df.pheno.2$cell.id.
    rownames(df) <- df.feature$coord.intron

# Subset overlapping cells
overlap <- intersect(df.pheno$Cell, df.pheno.2$cell.id)
print(paste(length(df.pheno$Cell), " found by SingCellaR", sep=""))
print(paste(length(overlap), " found by STARsolo", sep=""))

df.pheno.2.small <- df.pheno.2[which(df.pheno.2$cell.id %in% overlap), , drop=FALSE]

# Match matrix to phenoData
    # Subset overlapping cells
    df.small <- df[, df.pheno.2.small$cell.id.]
    
    # Update cell ids
    colnames(df.small) <- df.pheno.2.small$cell.id

# Save as new object
df.SRR9008752 <- df.small
df.pheno.SRR9008752 <- df.pheno.2.small
df.feature.SRR9008752 <- df.feature

# Check alignment
table(colnames(df.SRR9008752)==df.pheno.SRR9008752$cell.id)
table(rownames(df.SRR9008752)==df.feature.SRR9008752$coord.intron)

######################################################################
################# SRR9008753 (CARDIO DAY 10) #########################
######################################################################

# Read barcode file
    # Read file
    path <- "../../data/STARsolo/SRR9008753/Solo_out/SJ/raw/"
    file <- "barcodes.tsv"
    df.pheno.2 <- as.data.frame(fread(paste(path, file, sep=""), sep="\t", header=FALSE, stringsAsFactors=FALSE))
    
    # Retrieve cell ids
    names(df.pheno.2) <- "cell.id."
    df.pheno.2$cell.id <- df.pheno.2$cell.id.

    # Indicate donor id
    df.pheno.2$cell.id <- paste(df.pheno.2$cell.id, "_SRR9008753", sep="")

# Read feature file
    # Read file
    path <- "../../data/STARsolo/SRR9008753/Solo_out/SJ/raw/"
    file <- "features.tsv"
    df.feature <- as.data.frame(fread(paste(path, file, sep=""), sep="\t", header=FALSE, stringsAsFactors=FALSE))
    
    # Subset relevant columns
    cols <- c("V1", "V2", "V3")
    df.feature <- df.feature[, cols]
    names(df.feature) <- c("chr", "start", "end")
    df.feature$chr <- paste("chr", df.feature$chr, sep="")
    
    # Create unique junction id
    . <- data.frame("coord.intron"=paste(df.feature$chr, df.feature$start, df.feature$end, sep=":"))
    df.feature <- cbind.data.frame(., df.feature)
    
# Read matrix
    # Read file
    path <- "../../data/STARsolo/SRR9008753/Solo_out/SJ/raw/"
    file <- "matrix.mtx"
    df <- readMM(paste(path, file, sep=""))
    
    # Check alignment
    ncol(df)==nrow(df.pheno.2)
    nrow(df)==nrow(df.feature)

    # Annotate columns
    colnames(df) <- df.pheno.2$cell.id.
    rownames(df) <- df.feature$coord.intron

# Subset overlapping cells
overlap <- intersect(df.pheno$Cell, df.pheno.2$cell.id)
print(paste(length(df.pheno$Cell), " found by SingCellaR", sep=""))
print(paste(length(overlap), " found by STARsolo", sep=""))

df.pheno.2.small <- df.pheno.2[which(df.pheno.2$cell.id %in% overlap), , drop=FALSE]

# Match matrix to phenoData
    # Subset overlapping cells
    df.small <- df[, df.pheno.2.small$cell.id.]
    
    # Update cell ids
    colnames(df.small) <- df.pheno.2.small$cell.id

# Save as new object
df.SRR9008753 <- df.small
df.pheno.SRR9008753 <- df.pheno.2.small
df.feature.SRR9008753 <- df.feature

# Check alignment
table(colnames(df.SRR9008753)==df.pheno.SRR9008753$cell.id)
table(rownames(df.SRR9008753)==df.feature.SRR9008753$coord.intron)

######################################################################
############################# MERGE ##################################
######################################################################

# Merge phenoData
df.pheno.merged <- rbind.data.frame(df.pheno.SRR9008754,
                                    df.pheno.SRR9008755,
                                    df.pheno.SRR9008752,
                                    df.pheno.SRR9008753
                                    )
                                    
# Merge featureData
df.feature.merged <- rbind.data.frame(df.feature.SRR9008754,
                                      df.feature.SRR9008755,
                                      df.feature.SRR9008752,
                                      df.feature.SRR9008753
                                      )

df.feature.merged <- unique(df.feature.merged)

# Merge matrix
    # Merge into list
    df.list <- list(df.SRR9008754,
                    df.SRR9008755,
                    df.SRR9008752,
                    df.SRR9008753
                    )

    # Merge
    .list <- list()
    
    for(i in 1:length(df.list)) {
    
        if(i == 1){
        
            # Retrieve base matrix
            df <- df.list[[1]]
            
            # Save into list
            .list[[1]] <- df
        
        } else if(i != 1) {
        
            # Retrieve base matrix
            df.1 <- .list[[1]]
            
            # Retrieve matrix to add
            df.2 <- df.list[[i]]
           
            # Common junctions
                # Retrieve junction ids
                coord.introns <- intersect(rownames(df.1), rownames(df.2))
                length(coord.introns)
                
                # Subset
                    # Matrix 1
                    df.1.common <- df.1[coord.introns, ]
                    
                    # Matrix 2
                    df.2.common <- df.2[coord.introns, ]
            
            # Unique to Matrix 1
                # Retrieve junction ids
                coord.introns <- setdiff(rownames(df.1), rownames(df.2))
                length(coord.introns)
                
                # Subset Matrix 1
                df.1.small.1 <- df.1[coord.introns, ]
                
                # Create dummy matrix for Matrix 2
                . <- matrix(0, nrow=length(coord.introns), ncol=ncol(df.2))
                . <- Matrix(., sparse=TRUE)
                colnames(.) <- colnames(df.2)
                rownames(.) <- coord.introns
                df.2.small.1 <- .
        
            # Unique to Matrix 2
                # Retrieve junction ids
                coord.introns <- setdiff(rownames(df.2), rownames(df.1))
                length(coord.introns)
                
                # Create dummy matrix for Matrix 1
                . <- matrix(0, nrow=length(coord.introns), ncol=ncol(df.1))
                . <- Matrix(., sparse=TRUE)
                colnames(.) <- colnames(df.1)
                rownames(.) <- coord.introns
                df.1.small.2 <- .
                
                # Subset Matrix 2
                df.2.small.2 <- df.2[coord.introns, ]

            # Merge individual matrix
            df.1 <- rbind(df.1.common,
                          df.1.small.1,
                          df.1.small.2
                          )
          
            df.2 <- rbind(df.2.common,
                          df.2.small.1,
                          df.2.small.2
                          )
    
            # Check alignment
            print(table(rownames(df.1)==rownames(df.2)))
                    
            # Merge all matrices
            df <- cbind(df.1, df.2)
            
            # Replace superceded matrix
            .list[[1]] <- df
            
        }
        
        # Track progress
        print(paste("Matrix ", i, " done", sep=""))
    
    }

df.merged <- .list[[1]]

# Match featureData rows
nrow(df.merged) ; nrow(df.feature.merged)
df.merged <- df.merged[df.feature.merged$coord.intron, ]

# Check alignment
table(colnames(df.merged)==df.pheno.merged$cell.id)
table(rownames(df.merged)==df.feature.merged$coord.intron)

# Save files
    # Matrix
    path <- "../../data/MARVEL_input/SJ_STARsolo/"
    file <- "matrix_counts.mtx"
    writeMM(df.merged, file=paste(path, file, sep=""))
    
    # phenoData
    path <- "../../data/MARVEL_input/SJ_STARsolo/"
    file <- "phenoData.txt"
    write.table(df.pheno.merged[,"cell.id",drop=FALSE], paste(path, file, sep=""), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
    
    # featureData
    path <- "../../data/MARVEL_input/SJ_STARsolo/"
    file <- "featureData.txt"
    write.table(df.feature.merged[,"coord.intron", drop=FALSE], paste(path, file, sep=""), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
