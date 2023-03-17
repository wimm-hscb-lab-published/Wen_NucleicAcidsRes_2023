# Set working directory
#setwd("/Users/seanwen/Documents/MARVEL/Github_WIMM/Wen_NucleicAcidsRes_2023/scripts/scripts_01_preprocess_singcellar/")

# Load packages
library(SingCellaR)
library(data.table)

# Read SingCellaR files
sample.1 <- "SRR9008753.SingCellaR.rdata"
sample.2 <- "SRR9008754.SingCellaR.rdata"

object <- new("SingCellaR_Int")
object@dir_path_SingCellR_object_files <- "../../data/SingCellaR_output/individual/"
object@SingCellR_object_files=c(sample.1, sample.2)
preprocess_integration(object)

# Add relevant information to object
df.pheno <- get_cells_annotation(object)

filter_cells_and_genes(object,
                       min_UMIs=min(df.pheno$UMI_count),
                       max_UMIs=max(df.pheno$UMI_count),
                       min_detected_genes=min(df.pheno$detectedGenesPerCell),
                       max_detected_genes=max(df.pheno$detectedGenesPerCell),
                       max_percent_mito=15,
                       genes_with_expressing_cells=10,
                       isRemovedDoublets=FALSE
                       )

# Annotate sample ID in metadata
    # Retrieve metadata
    df.pheno <- object@meta.data
        
    # Annotate
    df.pheno$sampleID[grep("_SRR9008753", df.pheno$Cell, fixed=TRUE)] <- "SRR9008753"
    df.pheno$sampleID[grep("_SRR9008754", df.pheno$Cell, fixed=TRUE)] <- "SRR9008754"
    table(df.pheno$sampleID)

    # Save metadata
    object@meta.data <- df.pheno
    
# Normalise counts
normalize_UMIs(object, use.scaled.factor = T)

# Retrieve highly variable genes
get_variable_genes_by_fitting_GLM_model(object,
                                        mean_expr_cutoff=0.1,
                                        disp_zscore_cutoff=0.1
                                        )
                                        
plot_variable_genes(object)

# Reduce dimensions
    # PCA
    runPCA(object,
          use.components=50,
          use.regressout.data=FALSE
          )
    
    # Determine PCs for UMAP
    plot_PCA_Elbowplot(object)
    
    # UMAP
    #runUMAP(object,
            #dim_reduction_method="pca",
            #n.dims.use=10,
            #n.neighbors=30,
            #uwot.metric="euclidean"
            #)
                        
    # Plot UMAP
    #plot_umap_label_by_a_feature_of_interest(object,
                                             #feature="sampleID",
                                             #point.size = 0.5
                                             #)
                                             
    # tSNE
    runTSNE(object,
            dim_reduction_method="pca",
            n.dims.use=10,
            )
                        
    # Plot tSNE
    plot_tsne_label_by_a_feature_of_interest(object,
                                             feature="sampleID",
                                             point.size = 0.5
                                             )

# Save integrated R object
#path <- "/Users/seanwen/Documents/U2AF1_2019/Ou/Gene/Integrated/"
#file <- "Integrated_iPSC_CardioDay10.rdata"
#save(object, file=paste(path, file, sep=""))

# Retrieve coordinates: tSNE
    # Retrieve UMAP coordinates
    df.coord <- get_tsne.result(object)
    df.coord <- df.coord[,c("Cell", "TSNE1", "TSNE2")]
    names(df.coord) <- c("cell.id", "x", "y")
    
    # Save file for MARVEL input
    path <- "../../data/SingCellaR_output/integrated/"
    file <- "dim_red_coordinates_iPSC_CardioDay10.txt"
    write.table(df.coord, paste(path, file, sep=""), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
