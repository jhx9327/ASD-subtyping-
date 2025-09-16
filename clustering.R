library(NbClust)
library(corrplot)
library(ggplot2)
library(factoextra)


data <- read.csv('meanMS_regional_corrected_neuroCombat_lm_asd.csv', header = FALSE)


sub_info <- read.csv('meanMS_regional_df.csv')
asd_info <- sub_info[which(sub_info$diagnosis=='ASD'),]

fviz_nbclust(data, hcut, method = "silhouette")+labs(subtitle = "Silhouette method")


#########
res <- fviz_nbclust(data, hcut, method = "silhouette")

df_sil <- res$data
df_sil$clusters <- as.numeric(as.character(df_sil$clusters))  


optimal_cluster <- df_sil$clusters[which.max(df_sil$y)]


ggplot(df_sil, aes(x = clusters, y = y)) +
  geom_line(size = 1.5) +
  geom_point(size = 4, color = "steelblue") +
  geom_vline(xintercept = optimal_cluster, linetype = "dashed", color = "red", size = 1) +  
  labs(title = "Optimal number of clusters", 
       subtitle = "Silhouette method",
       x = "Number of clusters", 
       y = "Average silhouette width") +
  scale_x_continuous(breaks = seq(min(df_sil$clusters), max(df_sil$clusters), by = 1)) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 16),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18),
    axis.line = element_line(size = 0.8, colour = "black"),
    panel.grid = element_blank(),
    panel.border = element_blank()
  )






d <- dist(data, method = "euclidean")  
hc <- hclust(d, method = "ward.D")


k <- 3
cluster_cut <- cutree(hc, k = k)


plot(hc, labels = FALSE, hang = -1, main = "Hierarchical Clustering Dendrogram")


rect.hclust(hc, k = k, border = 2:(k + 1))  # border = colors


n_iter <- 20000
n <- nrow(data)
stability_matrix <- matrix(0, n, n)
count_matrix <- matrix(0, n, n)

for (i in 1:n_iter) {
  sample_idx <- sample(1:n, size = floor(0.8 * n), replace = FALSE)
  data_sub <- data[sample_idx, ]
  
  d <- dist(data_sub, method = "euclidean")
  hc <- hclust(d, method = "ward.D")
  clusCut <- cutree(hc, k = 3)
  

  for (j in 1:length(sample_idx)) {
    for (l in 1:length(sample_idx)) {
      idx_j <- sample_idx[j]
      idx_l <- sample_idx[l]
      count_matrix[idx_j, idx_l] <- count_matrix[idx_j, idx_l] + 1
      if (clusCut[j] == clusCut[l]) {
        stability_matrix[idx_j, idx_l] <- stability_matrix[idx_j, idx_l] + 1
      }
    }
  }
}


consensus_matrix <- stability_matrix / count_matrix

save(consensus_matrix,file = paste0(n_iter, "_hclust_bootstrap_k_3.RData"))

load(file = paste0(n_iter, "_hclust_bootstrap_k_3.RData"))


colC <- colorRampPalette(c("white","royalblue4","seagreen","darkgoldenrod2","gold"))
# Plot Stability matrix
par(mfrow=c(1,1))
title <- "Stability Matrix"
try(corrplot(consensus_matrix,order="hclust",tl.pos='n', method="color",addgrid.col=NA,col=colC(100),is.corr = FALSE,cl.lim = c(0,1), mar=c(0,0,1,0),
             cl.cex = 2))


d <- dist(1-consensus_matrix, method = "euclidean")  
hc <- hclust(d, method = 'ward.D')
plot(hc, labels = FALSE, hang = -1, main = "Hierarchical Clustering Dendrogram")



fviz_nbclust(consensus_matrix, hcut, method = "silhouette")+labs(subtitle = "Silhouette method")

res <- fviz_nbclust(consensus_matrix, hcut, method = "silhouette")

df_sil <- res$data
df_sil$clusters <- as.numeric(as.character(df_sil$clusters))  


optimal_cluster <- df_sil$clusters[which.max(df_sil$y)]

ggplot(df_sil, aes(x = clusters, y = y)) +
  geom_line(size = 1.5) +
  geom_point(size = 4, color = "steelblue") +
  geom_vline(xintercept = optimal_cluster, linetype = "dashed", color = "red", size = 1) + 
  labs(title = "Optimal number of clusters", 
       subtitle = "Silhouette method",
       x = "Number of clusters", 
       y = "Average silhouette width") +
  scale_x_continuous(breaks = seq(min(df_sil$clusters), max(df_sil$clusters), by = 1)) +  
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 16),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18),
    axis.line = element_line(size = 0.8, colour = "black"),
    panel.grid = element_blank(),
    panel.border = element_blank()
  )




clusCut <- cutree(hc, k = 3)

asd_info$class_id_3 <- clusCut
asd_info$class_id_3 <- factor(asd_info$class_id_3)

write.csv(asd_info, "asd_info_class_id_3.csv", row.names = TRUE)


asd_info <- read.csv("asd_info_class_id_3.csv")
asd_info$class_id_3 <- factor(asd_info$class_id_3, levels = c(3, 2, 1))

##############
# plot stability matrix with class_id_3
sorted_idx <- order(asd_info$class_id_3)

consensus_sorted <- consensus_matrix[sorted_idx, sorted_idx]

class_sorted <- asd_info$class_id_3[sorted_idx]


colC <- colorRampPalette(c("white", "royalblue"))  

class_table <- table(class_sorted)
boundaries <- cumsum(class_table)  
starts <- c(1, head(boundaries, -1) + 1)  

corrplot(
  consensus_sorted,
  order = "original",
  tl.pos = "n",
  method = "color",
  col = colC(100),         
  is.corr = FALSE,
  cl.lim = c(0, 1),
  cl.cex = 1.5,             
  col.row = NULL,
  col.col = NULL,
  addrect = FALSE,        
  tl.cex = 0.6,
  mar = c(0, 0, 1, 0),
  main = "Consensus Matrix Sorted by class_id_3"
)


n <- nrow(consensus_sorted)  

for (i in seq_along(starts)) {
  from <- starts[i] - 0.5
  to <- boundaries[i] + 0.5
  
  rect(
    xleft   = from,
    ybottom = n - to + 1,
    xright  = to,
    ytop    = n - from + 1,
    border  = "red", lwd = 2
  )
}

###########

# Compute the analysis of variance

res.aov <- aov(ADI.R.RRB ~ class_id_3, data = asd_info)
summary(res.aov)
pairwise.t.test(asd_info$ADI.R.RRB, asd_info$class_id_3,
                p.adjust.method = "fdr")
t.test(asd_info[which(asd_info$class_id_3==1),]$ADI.R.RRB,
       asd_info[which(asd_info$class_id_3==2),]$ADI.R.RRB)

t.test(asd_info[which(asd_info$class_id_3==1),]$ADI.R.RRB,
       asd_info[which(asd_info$class_id_3==3),]$ADI.R.RRB)

t.test(asd_info[which(asd_info$class_id_3==2),]$ADI.R.RRB,
       asd_info[which(asd_info$class_id_3==3),]$ADI.R.RRB)

ggplot(asd_info, aes(x = class_id_3, y = ADI.R.RRB)) +
  geom_boxplot(fill = "lightblue", outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.3, size = 4) +
  theme_minimal() +
  labs(x = "Subtypes", y = "ADI-R RRB Score") +
  scale_x_discrete(labels = c("3" = "subtype-1", 
                              "2" = "subtype-2", 
                              "1" = "subtype-3")) +
  geom_segment(aes(x = 1, xend = 3, y = 15, yend = 15), size=1.5) +  # 1 vs 3
  annotate("text", x = 2, y = 15.2, label = "*", size = 8, fontface = 'bold') +
  theme_minimal(base_size = 18) +
  theme(
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 18, face = "bold", colour = 'black'),
    axis.line.y = element_line(color = "black"),    
    axis.ticks.y = element_line(color = "black"),    
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold")
  )