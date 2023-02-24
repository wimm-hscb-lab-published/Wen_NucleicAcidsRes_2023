# Set working directory
#setwd("/Users/seanwen/Documents/MARVEL/Github_WIMM/Wen_NucleicAcidsRes_2023/scripts/scripts_01_preprocess_singcellar/")

# Load packages
library(SingCellaR)
library(ggplot2)
library(data.table)
library(Matrix)

# Read STARsolo files
data_matrices_dir <- "../../data/STARsolo/SRR9008754/Solo_out/Gene/filtered/"
Sample <- new("SingCellaR")
Sample@dir_path_10x_matrix <- data_matrices_dir
Sample@sample_uniq_id <- "SRR9008754"

load_matrices_from_cellranger(Sample, cellranger.version=3)

# Check
Sample

# Compute % mt
process_cells_annotation(Sample, mito_genes_start_with="MT-")

# QC plots
    # Histogram
    plot_cells_annotation(Sample, type="histogram")

    # Boxplot
    plot_cells_annotation(Sample, type="boxplot")
    
    # UMIs vs. the number of detected genes/cell
        # Draft
        plot_UMIs_vs_Detected_genes(Sample)
        
        # Definition
        data <- get_cells_annotation(Sample)
        x <- data$UMI_count
        y <- data$detectedGenesPerCell
        maintitle <- ""
        ytitle <- "Genes Detected per Cell"
        xtitle <- "UMI Count per Cell"
        xmin <- 0 ; xmax <- max(x) ; xinterval <- 10000
        ymin <- signif(min(y)-250, digits=1) ; ymax <- max(y) ; yinterval <- 1000
    
        # Plot
        ggplot() +
            geom_point(data, mapping=aes(x=x, y=y), size=0.1) +
            scale_x_continuous(breaks=seq(xmin, xmax, by=xinterval), limits=c(xmin, xmax)) +
            scale_y_continuous(breaks=seq(ymin, ymax, by=yinterval), limits=c(ymin, ymax)) +
            labs(title=maintitle, x=xtitle, y=ytitle) +
            theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                panel.border=element_blank(),
                plot.title=element_text(hjust = 0.5, size=15),
                plot.subtitle=element_text(hjust = 0.5, size=15),
                axis.line.y.left = element_line(color="black"),
                axis.line.x = element_line(color="black"),
                axis.title=element_text(size=12),
                axis.text=element_text(size=12),
                axis.text.x=element_text(size=10, colour="black", angle=45,  vjust=1.0, hjust=1.0),
                axis.text.y=element_text(size=10, colour="black"),
                #legend.position="none",
                legend.title=element_text(size=8),
                legend.text=element_text(size=8)
                )
                
        # Indicate outlier cells (check plot 1st)
        data$outlier <- ifelse(data$UMI_count > 45000, "yes", "no")
        data$outlier <- factor(data$outlier, levels=c("no", "yes"))
        z <- data$outlier
        
        # Plot
        plot <- ggplot() +
            geom_point(data, mapping=aes(x=x, y=y, color=z), size=0.1) +
            scale_x_continuous(breaks=seq(xmin, xmax, by=xinterval), limits=c(xmin, xmax)) +
            scale_y_continuous(breaks=seq(ymin, ymax, by=yinterval), limits=c(ymin, ymax)) +
            scale_color_manual(values=c("black", "red")) +
            labs(title=maintitle, x=xtitle, y=ytitle) +
            theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                panel.border=element_blank(),
                plot.title=element_text(hjust = 0.5, size=15),
                plot.subtitle=element_text(hjust = 0.5, size=15),
                axis.line.y.left = element_line(color="black"),
                axis.line.x = element_line(color="black"),
                axis.title=element_text(size=12),
                axis.text=element_text(size=12),
                axis.text.x=element_text(size=10, colour="black", angle=45,  vjust=1.0, hjust=1.0),
                axis.text.y=element_text(size=10, colour="black"),
                legend.position="none",
                legend.title=element_text(size=8),
                legend.text=element_text(size=8)
                )
        
        # Save file
        #path <- "/Users/seanwen/Documents/U2AF1_2019/Ou/Gene/SRR9008754/QC/"
        #dir.create(path, showWarnings=FALSE)
        #file <- "Genes Detected vs UMI Counts.pdf"
        #ggsave(paste(path, file, sep=""), plot, width=5, height=5)
        
# Filter
filter_cells_and_genes(Sample,
                       min_UMIs=1000,
                       max_UMIs=45000,
                       min_detected_genes=500,
                       max_detected_genes=max(y),
                       max_percent_mito=15,
                       genes_with_expressing_cells=10,
                       isRemovedDoublets=FALSE
                       )
                       
# Normalise UMI counts
normalize_UMIs(Sample, use.scaled.factor=TRUE)

# Regress out unwanted source of variations
remove_unwanted_confounders(Sample, residualModelFormulaStr="~UMI_count+percent_mito")

# Retrieve highly-variable genes (needed to load object for integration)
get_variable_genes_by_fitting_GLM_model(Sample,
                                        mean_expr_cutoff=0.05,
                                        disp_zscore_cutoff = 0.05
                                        )
                                        
# Save object
path <- "../../data/SingCellaR_output/individual/"
dir.create(path, showWarnings=FALSE)
file <- "SRR9008754.SingCellaR.rdata"
save(Sample, file=paste(path, file, sep=""))
