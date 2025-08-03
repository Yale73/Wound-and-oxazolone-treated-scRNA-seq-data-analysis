library("Rcpp")
library("RcppEigen")
library("ggsci")
library("viridis")
library("tidyverse")
library(Seurat)
library(scater)
library(scran)
library(dplyr)
library(Matrix)
library(cowplot)
library(ggplot2)
library(CellChat)
#devtools::install_github("sqjin/CellChat")
library(patchwork)
library(harmony)
library(scibetR)
library(CellChat)
library(patchwork)
library(ggplot2)
library(Seurat)
library(ggalluvial)#绘制桑基图
library(expm)
library(sna)
library(NMF)
library(ComplexHeatmap)
options(stringsAsFactors = FALSE)


immune.combined <- readRDS("H:/Subjects/育飞-Baso_wound/Basophils_Wound healing/Wound data analysis/immune_combined_cellrep_data.rds")

###########Merge interest cell types
Idents(immune.combined) <- immune.combined$celltype
Interg <- immune.combined

###############For CellChat
Idents(Interg) <- Interg@active.ident
Interg$Ident <- Interg@active.ident
levels(Interg$Ident)

###make meta
data.input = Matrix(Interg@assays$RNA@data)#需标准化的基因表达量矩阵和细胞分组信息文件
#dat <- as.data.frame(data.input)#可以查看表达矩阵 [1] 17328  7563
meta = Interg@meta.data


#分别构建cellchat文件
##构建BP组cellchat数据
cellchat <- createCellChat(object = data.input, #支持normalized表达矩阵，Seurat对象，和SingleCellExperiment对象
                           meta = meta, #meta文件
                           group.by = 'celltype') #meta中的细胞分类列


#重新加载数据进行cellchat分析#################################################
##Wounded group
cell.use = colnames(Interg)[Interg$group == 'Wounded'] #提取LS的细胞名称
data.input = data.input[, cell.use]#提取LS表达矩阵
meta = meta[cell.use, ]#提取LS细胞信息
identical(rownames(meta),colnames(data.input)) #检查矩阵列名和分组文件行名是否一致


##构建BP组cellchat数据
cellchat.Wounded <- createCellChat(object = data.input, #支持normalized表达矩阵，Seurat对象，和SingleCellExperiment对象
                                   meta = meta, #meta文件
                                   group.by = 'celltype') #meta中的细胞分类列

saveRDS(cellchat.Wounded, "H:/Subjects/育飞-Baso_wound/Basophils_Wound healing/Analysis_Final/cellchat__CR_Wounded.rds")
#cellchatHC组的cellchat分析
cellchat <- cellchat.Wounded
cellchat <- setIdent(cellchat, ident.use = 'celltype') #将label设置为显示的默认顺序
levels(cellchat@idents) #查看celltype和factor顺序
table(cellchat@idents) #每个celltype中的细胞数
cellchat@DB <- CellChatDB.mouse##设置配受体数据库(CellChatDB):
cellchat <- subsetData(cellchat)##信号基因的表达矩阵子集 赋值到cellchat@data.Signaling
cellchat <- identifyOverExpressedGenes(cellchat)##鉴定与每个细胞亚群相关的过表达信号基因
cellchat <- identifyOverExpressedInteractions(cellchat)##识别过表达基因配体-受体互作
cellchat <- projectData(cellchat, PPI.mouse)#将基因表达数据映射到PPI网络(可跳过)
cellchat@idents = droplevels(cellchat@idents, exclude = setdiff(levels(cellchat@idents),unique(cellchat@idents)))#将空的cluster过滤掉，以防影响后续分析
cellchat <- computeCommunProb(cellchat, raw.use = TRUE) #计算细胞通讯概率
cellchat <- filterCommunication(cellchat, min.cells = 10)##细胞通讯过滤
cellchat <- computeCommunProbPathway(cellchat)#计算信号通路水平上的通讯概率
cellchat <- aggregateNet(cellchat)#计算细胞对间通讯的数量和概率强度
cellchat <- netAnalysis_computeCentrality(cellchat,#计算网络中心性权重：识别每类细胞在信号通路中的角色/作用
                                          slot.name = "netP")

saveRDS(cellchat, "H:/Subjects/育飞-Baso_wound/Basophils_Wound healing/Analysis_Final/cellchat_CR_Wounded_analysis.rds")

##cellchat <- readRDS("H:/Subjects/育飞-Baso_wound/Basophils_Wound healing/Analysis_Final/cellchat_CR_Wounded_analysis.rds")
#######################################################
levels(cellchat@idents)

netVisual_heatmap(cellchat, measure = "weight", color.heatmap = c("gray90", "slategray3","#b2182b"), title.name = "Interaction strength")
netVisual_heatmap(cellchat, measure = "count", color.heatmap = c("gray90", "slategray3","#b2182b"), title.name = "Interaction number")


dev.off()

#3.信号通路水平的通讯差异

# 传出信号通路水平热图
pathway.union <- cellchat@netP$pathways

netAnalysis_signalingRole_heatmap(cellchat,
                                  pattern = "outgoing", #传出
                                  signaling = pathway.union,
                                  title = names(cellchat),
                                  width = 5,
                                  height = 20)

#传入信号通路水平热图
netAnalysis_signalingRole_heatmap(cellchat,
                                  pattern = "incoming", #传入
                                  signaling = pathway.union,
                                  title = names(cellchat),
                                  width = 5, height = 20,
                                  color.heatmap = "GnBu")


#总体信号通路水平热图
netAnalysis_signalingRole_heatmap(cellchat,
                                  pattern = "all", #总体
                                  signaling = pathway.union,
                                  title = names(cellchat),
                                  width = 6, height = 20,
                                  color.heatmap = "OrRd")




##Figure S6A
netVisual_bubble(cellchat,
                 sources.use = c("Baso"),
                 targets.use = c("M/MdM","cDC" ,"M2 Mac","Neu","DETC","gdT17", "NK" ,"CD4+T","Treg","mDC" ,"Mono","LC" ,"Baso" ,"pDC","Mast" ),
                 max.dataset = 2,
                 title.name = "Baso targeting signals after skin wound",
                 angle.x = 45,
                 remove.isolate = T) #Increased为比较组通讯概率更强的配受体对信息

###############################################
pathways.show <- c("IL4")
par(mfrow = c(1,2), xpd = TRUE)
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds", title.name = paste(pathways.show, "signaling ", names(cellchat)))

netAnalysis_contribution(cellchat, signaling = pathways.show, title = "Contribution of each L-R pair")

##########Figure S6B
pairLR.IL4 <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE) ########Three in total
LR.show <- pairLR.IL4[1,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)

netVisual_aggregate(cellchat,signaling = pathways.show, vertex.receiver =vertex.receiver)
