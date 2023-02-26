# Set working directory
#setwd("/Users/seanwen/Documents/MARVEL/Github_WIMM/Wen_NucleicAcidsRes_2023/scripts/scripts_01_preprocess_singcellar/")

# Load packages
library(SingCellaR)
library(data.table)

# Read SingCellaR files
sample.1 <- "SRR9008752.SingCellaR.rdata"
sample.2 <- "SRR9008753.SingCellaR.rdata"
sample.3 <- "SRR9008754.SingCellaR.rdata"
sample.4 <- "SRR9008755.SingCellaR.rdata"

object <- new("SingCellaR_Int")
object@dir_path_SingCellR_object_files <- "../../data/SingCellaR_output/individual/"
object@SingCellR_object_files=c(sample.1, sample.2, sample.3, sample.4)
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

# Add Sample ID
    # Retrieve metadata
    df.pheno <- object@meta.data
        
    # Annotate
    df.pheno$sampleID[grep("_SRR9008752", df.pheno$Cell, fixed=TRUE)] <- "SRR9008752"
    df.pheno$sampleID[grep("_SRR9008753", df.pheno$Cell, fixed=TRUE)] <- "SRR9008753"
    df.pheno$sampleID[grep("_SRR9008754", df.pheno$Cell, fixed=TRUE)] <- "SRR9008754"
    df.pheno$sampleID[grep("_SRR9008755", df.pheno$Cell, fixed=TRUE)] <- "SRR9008755"
    table(df.pheno$sampleID)

    # Save metadata
    object@meta.data <- df.pheno

# Normalise counts (Checking w/ GW if this is neccessary)
normalize_UMIs(object,use.scaled.factor = T)

# Retrieve highly variable genes
get_variable_genes_by_fitting_GLM_model(object,
                                        mean_expr_cutoff=0.1,
                                        disp_zscore_cutoff=0.1
                                        )
                                        
plot_variable_genes(object)

# Save integrated R object
path <- "../../data/SingCellaR_output/integrated/"
file <- "Integrated_iPSC_CardioDay2_4_10.rdata"
save(object, file=paste(path, file, sep=""))
