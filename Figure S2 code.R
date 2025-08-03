Immune <- readRDS("F:/mouse paper object/Immune.rds")
MB <- subset(Immune, idents="M/B")
DefaultAssay(MB) <- "integrated"
MB <- FindVariableFeatures(MB)
MB <- ScaleData(MB)
MB <- RunPCA(MB)
########Decide the PC numbers
# Determine percent of variation associated with each PC
pct <- MB[["pca"]]@stdev / sum(MB[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1
 
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2
 
pcs <- min(co1, co2)
pcs
 
plot_df <- data.frame(pct = pct, cumu = cumu, rank = 1:length(pct))

# Elbow plot to visualize 
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()

################################ runUMAP
MB <- RunUMAP(MB, dims = 1:11)
MB <- RunTSNE(MB, dims = 1:11)
MB <- FindNeighbors(MB)
MB <- FindClusters(MB, resolution = 0.05)
DimPlot(MB, label=T)+NoLegend()

clustree(MB)
MB <- subset(MB, exp=="OXA")


DimPlot(MB, reduction="umap", label=T, pt.size = 1,  label.size=5, cols = c("slateblue", "chocolate"))+NoLegend()+
  theme(text=element_text(family="Arial",size=15)) +
  theme(panel.background = element_rect(fill='white', colour='black'), 
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',family="Arial",size=18),axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'), 
        axis.ticks.margin = unit(0.6,"lines"),
        axis.line = element_line(colour  = "black"), 
        axis.title.x=element_text(colour='black', size=12),
        axis.title.y=element_text(colour='black', size=12),
        axis.text=element_text(colour='black',size=12),
        legend.title=element_blank(),
        legend.text=element_text(family="Arial", size=15),
        legend.key=element_blank())+
  theme(plot.title = element_text(size=22,colour = "black",face = "bold"))+ 
  guides(colour = guide_legend(override.aes = list(size=5)))

DefaultAssay(MB) <- "RNA"
p1 <- FeaturePlot(MB, reduction="umap", features = "Mcpt8", pt.size=1, order=T)+
  theme(panel.background = element_rect(fill='white', colour='black'), 
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',family="Arial",size=18),axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'), 
        axis.ticks.margin = unit(0.6,"lines"),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=12),
        axis.title.y=element_text(colour='black', size=12),
        axis.text=element_text(colour='black',size=12),
        plot.title=element_text(color="black", size=14, face="bold.italic"),
        legend.text=element_text(family="Arial", size=12),
        legend.key=element_blank())

p2 <- FeaturePlot(MB, reduction="umap", features = "Itgam", pt.size=1, order=T)+
  theme(panel.background = element_rect(fill='white', colour='black'), 
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',family="Arial",size=18),axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'), 
        axis.ticks.margin = unit(0.6,"lines"),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=12),
        axis.title.y=element_text(colour='black', size=12),
        axis.text=element_text(colour='black',size=12),
        plot.title=element_text(color="black", size=14, face="bold.italic"),
        legend.text=element_text(family="Arial", size=12),
        legend.key=element_blank())

p3 <- FeaturePlot(MB, reduction="umap", features = "Mcpt4", pt.size=1, order=T)+
  theme(panel.background = element_rect(fill='white', colour='black'), 
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',family="Arial",size=18),axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'), 
        axis.ticks.margin = unit(0.6,"lines"),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=12),
        axis.title.y=element_text(colour='black', size=12),
        axis.text=element_text(colour='black',size=12),
        plot.title=element_text(color="black", size=14, face="bold.italic"),
        legend.text=element_text(family="Arial", size=12),
        legend.key=element_blank())

p4 <- FeaturePlot(MB, reduction="umap", features = "Kit", pt.size=1, order=T)+
  theme(panel.background = element_rect(fill='white', colour='black'), 
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',family="Arial",size=18),axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'), 
        axis.ticks.margin = unit(0.6,"lines"),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=12),
        axis.title.y=element_text(colour='black', size=12),
        axis.text=element_text(colour='black',size=12),
        plot.title=element_text(color="black", size=14, face="bold.italic"),
        legend.text=element_text(family="Arial", size=12),
        legend.key=element_blank())

plot_grid(p1, p2, p3, p4,ncol=2)

Idents(MB) <- MB$seurat_clusters
MB <- RenameIdents(MB, `0` = "Basophils", `1` = "Mast")
MB[["Ident"]] <- Idents(object = MB)

DimPlot(MB, label=T)+NoLegend()


saveRDS(MB, "F:/wound healing flow data/single cell data/MB_oxa.rds")

dittoBarPlot(MB, "Ident", "stim", color.panel = c("slateblue", "chocolate"), x.labels.rotate=FALSE)+ scale_x_discrete(labels=c("OXA-C", "OXA"))+
  theme(panel.background = element_rect(fill='white', colour='black'), 
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',family="Arial",size=18),axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'), 
        axis.ticks.margin = unit(0.6,"lines"),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=12),
        axis.title.y=element_text(colour='black', size=12),
        axis.text=element_text(colour='black',size=12),
        plot.title=element_text(color="black", size=14, face="bold.italic"),
        legend.text=element_text(family="Arial", size=12),
        legend.key=element_blank())

##################DEGs
marker <- FindMarkers(MB, ident.1 = "Basophils")
write.xlsx(marker, "F:/wound healing flow data/single cell data/Baso_Cluster_marker.xlsx", rowNames=T)

##################DEGs
marker <- FindMarkers(MB, ident.1 = "OXA", ident.2 = "ETOH", group.by = "stim", subset.ident = "0")
write.xlsx(marker, "G:/wound healing/Baso_OXA_ETOH.xlsx", rowNames=T)

####Volcona Plot for DEG##################
library(openxlsx)
BasoDEG <- read.xlsx("F:/wound healing flow data/single cell data/Baso_OXA_EtOH.xlsx", sheet = 1, rowNames = T)

BasoDEG$thershold <- ifelse(BasoDEG$avg_log2FC > 1 & BasoDEG$p_val_adj < 0.05, "Up", 
                          ifelse(BasoDEG$avg_log2FC < -1 & BasoDEG$p_val_adj < 0.05, "Down", "Nonsig"))
BasoDEG$genelabels <- ""
BasoDEG$genelabels <- ifelse(BasoDEG$X1=="Mcpt4"
                           |BasoDEG$X1=="Ccl7"
                           |BasoDEG$X1=="Cpa3"
                           |BasoDEG$X1=="Ccl2"
                           |BasoDEG$X1=="Mcpt8"
                           |BasoDEG$X1=="Ccl3"
                           |BasoDEG$X1=="Ccl4"
                           |BasoDEG$X1=="Ccl6"
                           |BasoDEG$X1=="Ccl9"
                           |BasoDEG$X1=="Il1b"
                           |BasoDEG$X1=="Ifitm1"
                           |BasoDEG$X1=="Cxcl2"
                           |BasoDEG$X1=="Il4"
                           |BasoDEG$X1=="Il13"
                           |BasoDEG$X1=="Il6"
                           |BasoDEG$X1=="Cxcr4"
                           |BasoDEG$X1=="Gata2"
                           |BasoDEG$X1=="Kit"
                           |BasoDEG$X1=="Il1b"
                           |BasoDEG$X1=="Cxcl2"
                           |BasoDEG$X1=="Lilr4b"
                           |BasoDEG$X1=="Mcpt8", TRUE, FALSE)

cols <- c("Nonsig" = "darkgrey", "black" ="black", "Down" = "blue",  "Up" = "red")
ggplot(BasoDEG)+
  geom_point(aes(avg_log2FC, -log10(p_val_adj), col=thershold))+
  geom_text_repel(aes(avg_log2FC, -log10(p_val_adj)), label= ifelse(BasoDEG$genelabels, as.character(BasoDEG$X1), ""), size = 6, vjust = 2, nudge_y = 0)+
  theme(legend.title = element_blank(), text=element_text(size=20))+
  scale_color_manual(values = cols)+
  theme_bw()

####Heatmap to explore the difference between OXA vs IMQ