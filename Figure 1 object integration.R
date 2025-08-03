###############load original data

UnWounded1.data <-Read10X("E:/wound and un-wound data/Un-Wounded_1")
UnWounded2.data <-Read10X("E:/wound and un-wound data/Un-Wounded_2")

Wounded1.data <-Read10X("E:/wound and un-wound data/Wounded_1")
Wounded2.data <-Read10X("E:/wound and un-wound data/Wounded_2")
Wounded3.data <-Read10X("E:/wound and un-wound data/Wounded_3")

UnWounded1<- CreateSeuratObject(counts = UnWounded1.data, project = "UnWounded1")
UnWounded2<- CreateSeuratObject(counts = UnWounded2.data, project = "UnWounded2")
Wounded1<- CreateSeuratObject(counts = Wounded1.data, project = "Wounded1")
Wounded2<- CreateSeuratObject(counts = Wounded2.data, project = "Wounded2")
Wounded3<- CreateSeuratObject(counts = Wounded3.data, project = "Wounded3")


# QC FILTERING AND DATA EXAMINATION
UnWounded1[["percent.mt"]] <- PercentageFeatureSet(UnWounded1, pattern = "^mt-")
UnWounded2[["percent.mt"]] <- PercentageFeatureSet(UnWounded2, pattern = "^mt-")
Wounded1[["percent.mt"]] <- PercentageFeatureSet(Wounded1, pattern = "^mt-")
Wounded2[["percent.mt"]] <- PercentageFeatureSet(Wounded2, pattern = "^mt-")
Wounded3[["percent.mt"]] <- PercentageFeatureSet(Wounded3, pattern = "^mt-")

VlnPlot(UnWounded1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(UnWounded2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(Wounded1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(Wounded2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(Wounded3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


####subset the data
UnWounded1 <- subset(UnWounded1, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & percent.mt < 10)
UnWounded2 <- subset(UnWounded2, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & percent.mt < 12.5)
Wounded1 <- subset(Wounded1, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & percent.mt < 10)
Wounded2 <- subset(Wounded2, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & percent.mt < 10)
Wounded3 <- subset(Wounded3, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & percent.mt < 8)


###set up the condition
UnWounded1$cond <- "UnWounded"
UnWounded2$cond <- "UnWounded"
Wounded1$cond <- "Wounded"
Wounded2$cond <- "Wounded"
Wounded3$cond <- "Wounded"


############Subset CD45+ cells
DefaultAssay(Wound_Integrated_Log_LogNormalize) <- "RNA"
Immune <- subset(Wound_Integrated_Log_LogNormalize, Ptprc > 0)
DimPlot(Immune, reduction = "tsne")
Idents(Immune) <- Immune$orig.ident
UnWounded1 <- subset(Immune, idents = "UnWounded1")
UnWounded2 <- subset(Immune, idents = "UnWounded2")
Wounded1 <- subset(Immune, idents = "Wounded1")
Wounded2 <- subset(Immune, idents = "Wounded2")
Wounded3 <- subset(Immune, idents = "Wounded3")



#merge 2 UnWounded data#
UnWounded.big <- merge(UnWounded1, y = UnWounded2, add.cell.ids = c("UnWounded1", "UnWounded2"), project = "UnWounded")
UnWounded.big

#merge 3 Wounded data#
Wounded.big <- merge(Wounded1, y = c(Wounded2, Wounded3), add.cell.ids = c("Wounded1", "Wounded2", "Wounded3"), project = "Wounded")
Wounded.big


#Combie object into a list, with each dataset as an element
Wound.list <- list(UnWounded.big, Wounded.big)
#Normalize and FindVariableFeatures in all datasets
for (i in 1:length(Wound.list )) {
  Wound.list [[i]] <- NormalizeData(Wound.list [[i]], verbose = TRUE, normalization.method = "RC")
  Wound.list [[i]] <- FindVariableFeatures(Wound.list [[i]], selection.method = "vst", 
                                           nfeatures = 2000, verbose = TRUE)
}


#Integration steps
Wound.anchors <- FindIntegrationAnchors(Wound.list, anchor.features = 2000, dims = 1:20, k.filter = 30, k.score = 5)
Wound.integrated <- IntegrateData(anchorset = Wound.anchors, dims = 1:20, k.weight = 10)

Wound.integrated <- ScaleData(Wound.integrated, verbose = TRUE)
Wound.integrated <- RunPCA(Wound.integrated, npcs = 30, verbose = TRUE)
ElbowPlot(Wound.integrated, ndims = 30)
Wound.integrated <- RunUMAP(Wound.integrated, reduction = "pca", dims = 1:15)
Wound.integrated <- RunTSNE(Wound.integrated, reduction = "pca", dims = 1:15, check_duplicates = FALSE)
Wound.integrated <- FindNeighbors(Wound.integrated, reduction = "pca", dims = 1:15)
Wound.integrated <- FindClusters(Wound.integrated, resolution = 0.3, algorithm = 4, verbose = T)
clustree(Wound.integrated, prefix = "integrated_snn_res.")

saveRDS(Wound.integrated, "E:/wound and un-wound data/Results figures-Seurat/CRM, algorithm=/Wound_Integrated_RC_Leiden.rds")

DimPlot(Wound.integrated, reduction = "tsne", split.by = "cond", pt.size = 1, label = T, label.size = 6,
        cols = c("Indianred4", "slateblue", "darkviolet", "deeppink3", "deepskyblue4",
                 "dodgerblue2", "firebrick3", "gold4", "green4", "hotpink4", "khaki4", "tan4", "purple3", "darkorange3"))+NoLegend()
DimPlot(Wound.integrated, reduction = "tsne", pt.size = 1, label = T)
DimPlot(Wound.integrated, reduction = "tsne", pt.size = 1, group.by = "cond", split.by = "cond", cols = c("deeppink3", "slateblue"))
DimPlot(Wound.integrated, reduction = "umap", pt.size = .85, label = TRUE, 
        cols = c("Indianred4", "slateblue", "darkviolet", "deeppink3", "deepskyblue4",
                 "dodgerblue2", "firebrick3", "gold4", "green4", "hotpink4", "khaki4", "tan4", "purple3", "darkorange3"))

FeaturePlot(Wound.integrated, reduction = "tsne", features = "Igha", label = T)

DoHeatmap(Wound.integrated, features = top10$gene, slot = "data",
          group.colors = c("Indianred4", "slateblue", "darkviolet", "deeppink3", "deepskyblue4",
                           "dodgerblue2", "firebrick3", "gold4", "green4", "hotpink4", "khaki4"),
          size = 4, angle = 0, label = F)+
  scale_fill_gradient(low = "gray98", high = "hotpink3")


##Tabulate cells by cluster ID, replicate, or both#
Idents(Wound.integrated) <- "seurat_clusters"
table(Idents(Wound.integrated))
table(Wound.integrated$orig.ident)
prop.table(table(Idents(Wound.integrated)))
table(Idents(Wound.integrated), Wound.integrated$orig.ident)
prop.table(table(Idents(Wound.integrated), Wound.integrated$orig.ident), margin = 2)


###Find Cluster markers and DEG
DefaultAssay(Wound.integrated) <-"RNA"
Wound.integrated <- NormalizeData(Wound.integrated, normalization.method = "LogNormalize", scale.factor = 10000)

for (i in 11:13){
  marker_i <- FindMarkers(Wound.integrated, ident.1 ="Wounded", ident.2 ="UnWounded", group.by = "cond", verbose = TRUE, subset.ident = i, test.use = "MAST")
  tab  <- paste0("Cluster", " ", i)
  write.xlsx(marker_i, file = "E:/wound and un-wound data/Results figures-Seurat/New Immune integrate/Wou-UnWou-DEG.xlsx", sheetName = tab, append = TRUE)
}


DefaultAssay(Wound.integrated) <- "RNA"

for (i in 1:13){
  marker_i <- FindMarkers(Wound.integrated, ident.1 = i, verbose =TRUE, min.cells.group = 0, slot = "data", assay = "RNA", test.use = "MAST") 
  tab  <- paste0("Cluster", " ", i)
  write.xlsx(marker_i, file = "E:/wound and un-wound data/Results figures-Seurat/New Immune integrate/Wou-UnWou-FindMarkers.xlsx", sheetName = tab, append = TRUE)
}


