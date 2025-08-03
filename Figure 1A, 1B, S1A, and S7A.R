library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(openxlsx)
library(ggthemes)

#################data integration and non-immune remove
immune.combined <- readRDS("H:/Subjects/育飞-Baso_wound/Basophils_Wound healing/Wound data analysis/immune_combined_cellrep_data.rds")

##############Figure 1A
DimPlot(immune.combined, reduction = "umap",  pt.size = .8, group.by="celltype", label = T) + NoLegend()+
  scale_color_manual(values = colorRampPalette(brewer.pal(8, "Set1"))(15))+
  theme(text=element_text(family="Arial",size=15)) +
  theme(panel.background = element_rect(fill='white', colour='black'), 
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',family="Arial",size=18),axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'), 
        axis.ticks.margin = unit(0.6,"lines"),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=12),
        axis.title.y=element_text(colour='black', size=12),
        axis.text=element_text(colour='black',size=12),
        legend.title=element_blank(),
        legend.text=element_text(family="Arial", size=12),
        legend.key=element_blank())+
  theme(plot.title = element_text(size=15,colour = "black",face = "bold"))+ 
  guides(colour = guide_legend(override.aes = list(size=5)))

#########Figure S1B
rashx_markers <- c("Mafb",	"Msr1",		"Adgre1",	"Fcgr1", "Fcgr2b","Mrc1",
                   "Cd14",	"Ly6c2", "Plac8", 
                   "H2-Ab1", "H2-Aa", "Cd68",
                   "Aif1",	"C1qa",	"C1qb",	"C1qc","Arg1", "Chil3",	"Folr2","Cd163", 
                   "S100a8",	"S100a9",	"S100a11","Cxcr2",	"Csf3r",	"Lcn2",	
                   "Cd3d",	"Trdv4",	"Tcrg-C1", "Thy1", "Trdc",	
                   "Il17a", "Il17f",	"Il23r","Eomes",	"Tbx21",	"Ncr1",	"Gzma",	
                   "Rorc",	"Rora",	"Ccr6",	"Il18r1",
                   "Trbc2", "Trac",	"Cd4", "Icos",
                   "Foxp3",	"Il2ra",	"Ctla4",
                   "Cd80",	"Cd209a",	 "Mgl2",	"Clec10a",	"Itgax",
                   "Irf8", "Batf3",#cDC1
                   "Fscn1", "Cacnb3", "Ccr7",#mature DC
                   "Pclaf",	"Birc5",	"Tpx2", "Mki67",
                   "Cd207", "Epcam",
                   "Mcpt8",	"Fcer1a", "Il4",	"Cd200r3",	
                   "Siglech",	"Spib",	"Mzb1",	"Fcrla",
                   "Mcpt4", "Mrgprb1",	"Mrgprb2",	"Tpsab1")

DotPlot(immune.combined, features = rashx_markers, cols=c("#DDA0DD", "#6A0DAD", "#502380", "#290916"), assay = "RNA", col.min = 0.1, col.max = 1, dot.scale = 1,
        cluster.idents=F, group.by = "celltype")+
  scale_size(range = c(0, 5))+ 
  scale_size_area(max_size = 5)+ 
  #scale_color_viridis_c(name = 'log2 (count + 1)') + 
  cowplot::theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8, family="TT Times New Roman", face = "italic"),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, size = 9, family="TT Times New Roman"),
        legend.text = element_text(size=8),
        legend.title = element_text(size = 9)) +
  labs(x="", y="")


#########Figure 1B
Idents(immune.combined) <- "celltype"   # or use cluster identities
table(immune.combined$celltype, immune.combined$group)
library(dplyr)

# Tabulate counts
cell_counts <- table(immune.combined$celltype, immune.combined$group)
cell_counts_df <- as.data.frame(cell_counts)
colnames(cell_counts_df) <- c("CellType", "Group", "Count")

# Calculate proportion
cell_counts_df <- cell_counts_df %>%
  group_by(Group) %>%
  mutate(Proportion = Count / sum(Count))

cell_counts_df$Group <-factor(cell_counts_df$Group, levels = c("UnWounded", "Wounded"))

ggplot(cell_counts_df, aes(x = CellType, y = Proportion, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  ylab("Proportion") +
  theme_minimal() +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Cell Type Proportions Between Groups") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_manual(values = c("steelblue", "tomato")) +  # customize your colors here
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
    axis.text.y = element_text(angle = 0, hjust = 1, colour = "black"),
    legend.title = element_blank()
  )


##########################Figure S7A
rashx_markers <- c("Tnf", "Il1b", "Il6", "Il4", "Il13", "Il10", "Areg")
DotPlot(immune.combined, features = rashx_markers, cols=c("#DDA0DD", "#6A0DAD", "#502380", "#290916"), assay = "RNA", col.min = 0.1, col.max = 1, dot.min=0.1, dot.scale = 1,
        cluster.idents=F, group.by = "celltype")+
  scale_size(range = c(0, 5))+ 
  scale_size_area(max_size = 5)+ 
  #scale_color_viridis_c(name = 'log2 (count + 1)') + 
  cowplot::theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10, family="TT Times New Roman"),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, size = 10, family="TT Times New Roman", face = "italic"),
        legend.text = element_text(size=8),
        legend.title = element_text(size = 9)) +
  labs(x="", y="")+
  coord_flip()
