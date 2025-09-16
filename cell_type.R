library(EWCE)
library(ewceData)
library(ggplot2)
library(cowplot)
library(limma)
library(readxl)

home <- "PLS genes results"
setwd(home)

genes <- read.csv("subtype_1_PLS1_geneWeights.csv", header = FALSE)
genes_1 <- genes[which(genes$V3 > 5), 1]

genes <- read.csv("subtype_2_PLS1_geneWeights.csv", header = FALSE)
genes_2 <- genes[which(genes$V3 > 5), 1]

genes <- read.csv("subtype_3_PLS1_geneWeights.csv", header = FALSE)
genes_3 <- genes[which(genes$V3 > 5), 1]

# http://127.0.0.1:14261/library/EWCE/doc/EWCE.html#preparing-gene-lists 

# Use 100 bootstrap lists for speed, for publishable analysis use >10000
reps=10000

# Bootstrap significance test, no control for transcript length and GC content
set.seed(12345678)
# ctd <- ctd()

# ctd <- ewceData::ctd()
load('cell-type/ctd_aibsMultipleCrtxSmrtSeq_2019.rda')

# https://github.com/mgrantpeters/RBP_EWCE_analysis/blob/main/code/00-run_ewce.R
# https://github.com/RHReynolds/MarkerGenes/tree/master/specificity_matrices_new

full_results_1 = bootstrap_enrichment_test(sct_data=ctd,
                                         sctSpecies = 'human',
                                         genelistSpecies = 'human',
                                         hits=genes_1,
                                         reps=`reps`,
                                         annotLevel=1,
                                         geneSizeControl = TRUE)

full_results_2 = bootstrap_enrichment_test(sct_data=ctd,
                                           sctSpecies = 'human',
                                           genelistSpecies = 'human',
                                           hits=genes_2,
                                           reps=`reps`,
                                           annotLevel=1,
                                           geneSizeControl = TRUE)


full_results_3 = bootstrap_enrichment_test(sct_data=ctd,
                                           sctSpecies = 'human',
                                           genelistSpecies = 'human',
                                           hits=genes_3,
                                           reps=`reps`,
                                           annotLevel=1,
                                           geneSizeControl = TRUE)

# ewce_plot(full_results$results,mtc_method="BH")

results_1 <- full_results_1$results
results_2 <- full_results_2$results
results_3 <- full_results_3$results

upperLim <- max(max(abs(results_1$sd_from_mean), na.rm = TRUE),
                max(abs(results_2$sd_from_mean), na.rm = TRUE),
                max(abs(results_3$sd_from_mean), na.rm = TRUE))

ewce_plot_1 <- function (total_res, upperLim, mtc_method = "bonferroni", ctd = NULL, 
          align = "v", rel_heights = c(0.3, 1), axis = "lr") 
{
  requireNamespace("cowplot")
  requireNamespace("gridExtra")
  requireNamespace("grid")
  err_msg <- paste0("ERROR: Invalid mtc_method argument. Please see", 
                    " '?p.adjust' for valid methods.")
  if (!mtc_method %in% c("holm", "hochberg", "hommel", "bonferroni", 
                         "BH", "BY", "fdr", "none")) {
    stop(err_msg)
  }
  multiList <- TRUE
  if (is.null(total_res$list)) {
    multiList <- FALSE
  }
  make_dendro <- FALSE
  if (!is.null(ctd)) {
    make_dendro <- TRUE
    cells_in_ctd <- function(ctdIN, cells) {
      if (sum(!cells %in% colnames(ctdIN$specificity) == 
              0)) {
        return(1)
      }
      else {
        return(0)
      }
    }
    if (length(ctd[[1]]$plotting) > 0) {
      annotLevel <- which(unlist(lapply(ctd, FUN = cells_in_ctd, 
                                        cells = as.character(total_res$CellType))) == 
                            1)
      err_msg2 <- paste0("All of the cells within total_res should come", 
                         " from a single annotation layer of the CTD")
      if (length(annotLevel) == 0) {
        stop(err_msg2)
      }
    }
    if (length(ctd[[annotLevel]]$plotting) > 0) {
      total_res$CellType <- factor(total_res$CellType, 
                                   levels = ctd[[annotLevel]]$plotting$cell_ordering)
    }
  }
  if ("q" %in% colnames(total_res)) {
    total_res$q <- stats::p.adjust(total_res$p, method = mtc_method)
  }
  ast_q <- rep("", dim(total_res)[1])
  ast_q[total_res$q < 0.05] <- "*"
  total_res$ast_q <- ast_q
  total_res$sd_from_mean[total_res$sd_from_mean < 0] <- 0
  graph_theme <- theme_bw(base_size = 12, base_family = "Helvetica") + 
    theme(text = element_text(size = 14), axis.title.y = element_text(vjust = 0.6), 
          strip.background = element_rect(fill = "white"), 
          strip.text = element_text(color = "black"))
  # upperLim <- max(abs(total_res$sd_from_mean), na.rm = TRUE)
  total_res$y_ast <- total_res$sd_from_mean * 1.05
  total_res$abs_sd <- abs(total_res$sd_from_mean)
  if ("Direction" %in% colnames(total_res)) {
    # the_plot <- ggplot(total_res) + geom_bar(aes_string(x = "CellType", 
    #                                                     y = "abs_sd", fill = "Direction"), position = "dodge", 
    #                                          stat = "identity") + graph_theme
    the_plot <- ggplot(total_res) + geom_bar(aes_string(x = "CellType", 
                                                        y = "abs_sd"), position = "dodge", 
                                             stat = "identity", color='#1f77b4', fill='#1f77b4') + graph_theme
  }
  else {
    # the_plot <- ggplot(total_res) + geom_bar(aes_string(x = "CellType", 
    #                                                     y = "abs_sd", fill = "abs_sd"), stat = "identity") + 
    the_plot <- ggplot(total_res) + geom_bar(aes_string(x = "CellType", 
                                                        y = "abs_sd"), stat = "identity", color='#1f77b4', fill='#1f77b4') + 
      # scale_fill_gradient(low = "blue", high = "red") + 
      graph_theme + theme(legend.position = "none")
  }
  the_plot <- the_plot + theme(plot.margin = grid::unit(c(1, 0, 0, 0), "mm"), axis.text.x = element_text(angle = 55, hjust = 1, size=22, color = 'black'), 
                               axis.text.y= element_text(size=22), axis.title.y = element_text(size = 22)) + 
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) + xlab("") + theme(strip.text.y = element_text(angle = 0)) + 
    coord_cartesian(ylim = c(0, 1.1 * upperLim)) + 
    ylab("Std.Devs. from the mean")  +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  the_plot <- the_plot + scale_y_continuous(breaks = c(0, ceiling(upperLim * 0.66))) + geom_text(aes_string(label = "ast_q", x = "CellType", y = "y_ast"), size = 10)
  if (multiList) {
    the_plot <- the_plot + facet_grid("list ~ .", scales = "free", 
                                      space = "free_x")
  }
  output <- list()
  output$plain <- the_plot
  if (make_dendro) {
    the_dendrogram <- ctd[[annotLevel]]$plotting$ggdendro_horizontal + 
      theme(plot.margin = unit(c(0, 0, 0, 0), units = "cm")) + 
      scale_x_discrete(breaks = total_res$CellType)
    combined_plot <- cowplot::plot_grid(the_dendrogram, 
                                        the_plot, align = align, rel_heights = rel_heights, 
                                        axis = axis, ncol = 1)
    output$withDendro <- combined_plot
  }
  return(output)
}

ewce_plot_1(results_1,upperLim,mtc_method="BH")
ewce_plot_1(results_2,upperLim,mtc_method="BH")
ewce_plot_1(results_3,upperLim,mtc_method="BH")
