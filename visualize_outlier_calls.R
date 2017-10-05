args = commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(ggthemes)



# Loads in pvalue matrices
load_in_outlier_file <- function(file_name) {
    full_data <- read.table(file_name, header=TRUE, sep="\t")  
    cluster_names <- full_data[,1]  # First column is row labels (cluster names)
    pvalues <- full_data[,2:dim(full_data)[2]]  # 2nd column to end is pvalues
    rownames(pvalues) <- cluster_names  # assign row labels to be actually row labels
    return(-log10(pvalues + .0000001))
}


number_of_outlier_calls_per_sample_boxplot <- function(sample_info, pvalue_matrix, output_file, threshold, bgrd_distribution_string) {
    number_of_outlier_calls_per_sample <- colSums(pvalue_matrix > threshold, na.rm=TRUE)

    df <- data.frame(number_of_outlier_calls=number_of_outlier_calls_per_sample, sample_id=factor(sample_info$sample_id), status=factor(sample_info$status))
    # PLOT
    p <-ggplot(df, aes(reorder(sample_id,-number_of_outlier_calls), number_of_outlier_calls))
    p <- p + geom_bar(stat = "identity", aes(fill = status)) 
    p <- p + theme(text = element_text(size=15),axis.text.x = element_text(angle = 90, size=6,vjust=.41), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
    p <- p + labs(x = "Sample ID", y = paste0("# outlier calls (-log10(pvalue) > ", threshold,")"), title= bgrd_distribution_string)
    ggsave(p, file=output_file,width = 15,height=10.5,units="cm")

}

distribution_of_outlier_calls_boxplot <- function(sample_info, pvalue_matrix, output_file, bgrd_distribution_string) {
    num_samples <- dim(pvalue_matrix)[2]
    pvalues <- c()
    status <- c() 
    for (sample_num in 1:num_samples) {
        if (sample_num != 18 & sample_num != 22) {
        full_pval_vec <- pvalue_matrix[,sample_num]
        pval_vec_no_nan <- full_pval_vec[!is.na(full_pval_vec)]
        sample_status <- sample_info$status[sample_num]
        status <- c(status, as.character(rep(sample_status,length(pval_vec_no_nan))))
        pvalues <- c(pvalues, pval_vec_no_nan)
        }
    }
    df <- data.frame(status = factor(status), pvalues = pvalues)
    box_plot <- ggplot(df, aes(x=status, y=pvalues, fill=status)) + geom_boxplot(notch=TRUE)
    box_plot <- box_plot + theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    box_plot <- box_plot + labs(y = "-log10(pvalue)",title = bgrd_distribution_string)
    ggsave(box_plot, file=output_file,width = 15,height=10.5,units="cm")

}

outlier_vs_pli_score_scatter <- function(sample_info, pvalue_matrix, pli_info, output_file, bgrd_distribution_string) {
    num_samples <- dim(pvalue_matrix)[2]
    num_clusters <- dim(pvalue_matrix)[1]
    max_pli_score <- as.numeric(pli_info$max_pli_score)


    pvalues <- c()
    status <- c()
    pli_scores <- c()

    for (sample_num in 1:num_samples) {
        if (sample_num != 18 & sample_num != 22) {
        full_pval_vec <- pvalue_matrix[,sample_num]
        valid_indices <- !is.na(full_pval_vec) & !is.na(max_pli_score) & full_pval_vec > 1.5
        sample_status <- sample_info$status[sample_num]

        pval_vec_valid <- full_pval_vec[valid_indices]
        pli_scores_valid <- max_pli_score[valid_indices]

        status <- c(status, as.character(rep(sample_status,length(pval_vec_valid))))
        pvalues <- c(pvalues, pval_vec_valid)
        pli_scores <- c(pli_scores, pli_scores_valid)
        }
    }
    df <- data.frame(status = factor(status), pvalues = pvalues, pli_scores = pli_scores)
    # PLOT AWAY

    scatter <- ggplot(df, aes(x = pvalues, y = pli_scores, colour = status)) + geom_point(size=.3) 
    scatter <- scatter + theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    scatter <- scatter + labs(colour="status",x = "-log10(pvalue)", y = "pLI score", title = bgrd_distribution_string)
    ggsave(scatter, file=output_file,width = 15,height=10.5,units="cm")
}

pli_scores_of_outliers_boxplot <- function(sample_info, pvalue_matrix, pli_info, output_file, bgrd_distribution_string) {
    num_samples <- dim(pvalue_matrix)[2]
    num_clusters <- dim(pvalue_matrix)[1]
    max_pli_score <- as.numeric(pli_info$max_pli_score)


    pvalues <- c()
    status <- c()
    pli_scores <- c()

    for (sample_num in 1:num_samples) {
        if (sample_num != 18 & sample_num != 22) {
        full_pval_vec <- pvalue_matrix[,sample_num]
        valid_indices <- !is.na(full_pval_vec) & !is.na(max_pli_score) & full_pval_vec > 5.0
        sample_status <- sample_info$status[sample_num]

        pval_vec_valid <- full_pval_vec[valid_indices]
        pli_scores_valid <- max_pli_score[valid_indices]

        status <- c(status, as.character(rep(sample_status,length(pval_vec_valid))))
        pvalues <- c(pvalues, pval_vec_valid)
        pli_scores <- c(pli_scores, pli_scores_valid)
        }
    }
    df <- data.frame(status = factor(status), pvalues = pvalues, pli_scores = pli_scores)
    # PLOT AWAY
    box_plot <- ggplot(df, aes(x=status, y=pli_scores, fill=status)) + geom_violin()
    box_plot <- box_plot + theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    box_plot <- box_plot + labs(y = "pLI score",title = bgrd_distribution_string)
    ggsave(box_plot, file=output_file,width = 15,height=10.5,units="cm")

}

pvalues_of_high_pli_scores_boxplot <- function(sample_info, pvalue_matrix, pli_info, output_file, bgrd_distribution_string) {
    num_samples <- dim(pvalue_matrix)[2]
    num_clusters <- dim(pvalue_matrix)[1]
    max_pli_score <- as.numeric(pli_info$max_pli_score)


    pvalues <- c()
    status <- c()
    pli_scores <- c()

    for (sample_num in 1:num_samples) {
        if (sample_num != 18 & sample_num != 22) {
        full_pval_vec <- pvalue_matrix[,sample_num]
        valid_indices <- !is.na(full_pval_vec) & !is.na(max_pli_score) & max_pli_score > .9 & full_pval_vec > 3.5
        sample_status <- sample_info$status[sample_num]

        pval_vec_valid <- full_pval_vec[valid_indices]
        pli_scores_valid <- max_pli_score[valid_indices]

        status <- c(status, as.character(rep(sample_status,length(pval_vec_valid))))
        pvalues <- c(pvalues, pval_vec_valid)
        pli_scores <- c(pli_scores, pli_scores_valid)
        }
    }
    df <- data.frame(status = factor(status), pvalues = pvalues, pli_scores = pli_scores)
    # PLOT AWAY
    box_plot <- ggplot(df, aes(x=status, y=pvalues, fill=status)) + geom_violin()
    box_plot <- box_plot + theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    box_plot <- box_plot + labs(y = "-log10(pvalue)",title = bgrd_distribution_string)
    ggsave(box_plot, file=output_file,width = 15,height=10.5,units="cm")

}


# Wrapper to make a bunch of plots related to outlier levels
outlier_plotting_wrapper <- function(sample_info, pvalue_matrix, pli_info, output_root, bgrd_distribution_string) {
    # Make plot showing number of outlier calls greater than a certain shreshold
    threshold <- 3.3
    number_of_outlier_calls_per_sample_boxplot(sample_info, pvalue_matrix, paste0(output_root, "number_of_outlier_calls_per_sample_boxplot_",threshold,".png"), threshold, bgrd_distribution_string)
    threshold <- 5.0
    number_of_outlier_calls_per_sample_boxplot(sample_info, pvalue_matrix, paste0(output_root, "number_of_outlier_calls_per_sample_boxplot_",threshold,".png"), threshold, bgrd_distribution_string)

    distribution_of_outlier_calls_boxplot(sample_info, pvalue_matrix, paste0(output_root, "distribution_of_outlier_calls_boxplot.png"), bgrd_distribution_string)

    outlier_vs_pli_score_scatter(sample_info, pvalue_matrix, pli_info, paste0(output_root,"outlier_vs_pli_score_scatter.png"), bgrd_distribution_string)

    pli_scores_of_outliers_boxplot(sample_info, pvalue_matrix, pli_info, paste0(output_root,"pli_scores_of_outliers_boxplot.png"), bgrd_distribution_string)

    pvalues_of_high_pli_scores_boxplot(sample_info,pvalue_matrix,pli_info, paste0(output_root,"pvalues_of_pli_scores_boxplot.png"), bgrd_distribution_string)
}





options(warn=1)
#####################################
# Load in data
#####################################
output_dir = args[1]  # Directory to save plots to
gtex_rare_combined_outlier_file = args[2]  # Outlier calls for rare samples when background is gtex whole blood and rare samples
rare_only_outlier_file = args[3]  # Outlier calls for rare samples when background is rare samples
sra_gtex_rare_combined_outlier_file = args[4]  # Outlier calls for rare samples when background is gtex whole blood, sra, and rare samples
sample_info_file = args[5]  # File containing case control status of samples
gtex_rare_combined_pli_score_file  = args[6]  # File containing pli scores for each cluster
rare_only_pli_score_file = args[7]
sra_gtex_rare_combined_pli_score_file = args[8]


# Load in sample info
sample_info <- read.table(sample_info_file)
colnames(sample_info) <- c("sample_id","institution_id","status")


# Load in pvalue matrices
rare_only_pvalue <- load_in_outlier_file(rare_only_outlier_file)
gtex_rare_pvalue <- load_in_outlier_file(gtex_rare_combined_outlier_file)
sra_gtex_rare_pvalue <- load_in_outlier_file(sra_gtex_rare_combined_outlier_file)

# Load in pli score files
gtex_rare_pli <- read.table(gtex_rare_combined_pli_score_file, header=TRUE)
rare_pli <- read.table(rare_only_pli_score_file, header=TRUE)
sra_gtex_rare_pli <- read.table(sra_gtex_rare_combined_pli_score_file,header=TRUE)

# For rare only
outlier_plotting_wrapper(sample_info, rare_only_pvalue, rare_pli, paste0(output_dir,"bgrd_rare_only_"), "background distr. = rare samples")
outlier_plotting_wrapper(sample_info, gtex_rare_pvalue, gtex_rare_pli, paste0(output_dir,"bgrd_gtex_rare_"), "background distr. = gtex & rare samples")
outlier_plotting_wrapper(sample_info, sra_gtex_rare_pvalue, sra_gtex_rare_pli, paste0(output_dir,"bgrd_sra_gtex_rare_"), "background distr. = sra & gtex & rare samples")






