library(linkET)
library(tidymodels)
library(caret)
library(e1071)
library(glmnetUtils)
library(fmsb)
library(Seurat) 
library(dplyr) 
library(ggplot2)
library(patchwork)
library(harmony)
library(NMF)
library(pheatmap)
library(RColorBrewer)
library(igraph)
library(circlize)
library(RColorBrewer)
library(igraph)
library(ggsignif)
library(sjPlot)
library(clipr)
library(enrichplot)
library(DOSE)
#library("clusterProfiler")

library("org.Hs.eg.db")
#library("enrichplot")
library("ggplot2")
library("pathview")
library("ggnewscale")
library("DOSE")
library(stringr)
library(AnnotationHub)    #library导入需要使用的数据包
library('KEGG.db')

library(Seurat)
library(SCopeLoomR)
library(AUCell)
library(pbmcNIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(scRNAseq)
library(patchwork)
library(ggplot2) 
library(stringr)
library(circlize)
library(RcisTarget)

library(tidyr)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(vegan)

library(circlize)

library(tidyverse)

library(msigdbr)  
library(GSVA) 
library(GSEABase)
library(pheatmap)
library(limma)
library(BiocParallel)
library(msigdbr)
library(Seurat) 
library(dplyr) 
library(ggplot2)
library(patchwork)
library(monocle)
library(igraph)
library(knitr)
library(kableExtra)
library(table1) 
library(boot)
library(cellcall)

library(clusterProfiler)
library(rTensor)
library(scRank)
library(rjson)
library(easyTCGA)
library(glmnet)
library(Seurat)
library(dplyr)
library(viridis)









# === Step 1. 读取并构建表达矩阵 expr_data ===
# 读取第一行获取 cell barcodes
lines <- readLines("GSE141017_ALL.csv")
barcodes <- strsplit(lines[1], "\t")[[1]]  # 去掉第一列 header

# 读取剩余部分（跳过第一行）
expr_raw <- read.table("GSE141017_ALL.csv", header = FALSE, sep = "\t", skip = 1, stringsAsFactors = FALSE)

# 提取 gene symbol 与表达量矩阵
genes <- expr_raw[, 1]
expr_data <- as.matrix(expr_raw[, -1])
rownames(expr_data) <- genes
colnames(expr_data) <- barcodes

# 转换为数值型
mode(expr_data) <- "numeric"








# === Step 2. 读取 metadata 文件 ===
metadata_raw <- read.table("GSE141017_ALL_barcode_ident.csv",
                           header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 清洗 cell_barcode 字段（去除前缀 + 转换 .1 → -1）
metadata_raw$cell_barcode <- gsub("^X[_]*", "", metadata_raw$cell_barcode)
metadata_raw$cell_barcode <- gsub("\\.(\\d+)$", "-\\1", metadata_raw$cell_barcode)

# === Step 3. 特殊处理唯一重复 barcode（WT5M2_AGCTTGAGTACACCGC-1）
# 让两条记录都保留，并人为加后缀使其唯一（如 _1 和 _2）
metadata_raw <- metadata_raw %>%
  group_by(cell_barcode) %>%
  mutate(cell_barcode_unique = if (n() == 1) cell_barcode else paste0(cell_barcode, "_", row_number())) %>%
  ungroup()

# === Step 4. 构建 metadata_clean，设置行为 cell_barcode_unique
metadata_clean <- metadata_raw[, c("cell_barcode_unique", "ident", "time_point")]
rownames(metadata_clean) <- metadata_clean$cell_barcode_unique
metadata_clean$cell_barcode_unique <- NULL

# === Step 5. 同步修改 expr_data 的 column 名（确保唯一）
# 替换重复的 column barcode（这里只有 1 个重复）
colnames(expr_data)[colnames(expr_data) == "WT5M2_AGCTTGAGTACACCGC-1"] <- "WT5M2_AGCTTGAGTACACCGC-1_1"

# === Step 6. 构建 Seurat v5 对象，禁止默认过滤
seurat_obj <- CreateSeuratObject(
  counts = expr_data,
  assay = "RNA",
  meta.data = metadata_clean[colnames(expr_data), ])

# === Step 7. 添加 orig.ident 信息（从 cell barcode 中提取）

# 提取 orig.ident 为第一个非 X/下划线开头、直到第二个 "_" 之前的部分
seurat_obj$orig.ident <- gsub("^X*_*([A-Z0-9]+(?:_[A-Z0-9]+)?)_.*", "\\1", colnames(seurat_obj))



# === Step 8. 验证完整性
print(seurat_obj)
dim(LayerData(seurat_obj[["RNA"]], layer = "counts"))  # 应为 41140 x 41139
dim(seurat_obj@meta.data)                              # 应为 41139 x N
table(seurat_obj$orig.ident)
table(seurat_obj$time_point)


save(seurat_obj,file = '250408.Rdata')




pbmc <- NormalizeData(object = seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)


pbmc <- FindVariableFeatures(object = pbmc,selection.method = "vst", nfeatures = 2000)

all.gene <- rownames(pbmc)
pbmc <- ScaleData(pbmc,features = all.gene)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))


ElbowPlot(pbmc,ndim=50)
pbmc <- FindNeighbors(pbmc, dims = 1:20)#pa1-pc20
pbmc <- FindClusters(pbmc, resolution = 0.5)




pbmc <- RunUMAP(pbmc, dims = 1:20) 
DimPlot(pbmc, reduction = "umap",group.by = "seurat_clusters") 



save(pbmc,file = '250410.Rdata')










# === Step 1. 读取并构建表达矩阵 expr_data ===
# 读取第一行获取 cell barcodes
lines <- readLines("GSM4293555_Human.csv")
barcodes <- strsplit(lines[1], "\t")[[1]]  # 去掉第一列 header

# 读取剩余部分（跳过第一行）
expr_raw <- read.table("GSM4293555_Human.csv", header = FALSE, sep = "\t", skip = 1, stringsAsFactors = FALSE)

# 提取 gene symbol 与表达量矩阵
genes <- expr_raw[, 1]
expr_data <- as.matrix(expr_raw[, -1])
rownames(expr_data) <- genes
colnames(expr_data) <- barcodes

# 转换为数值型
mode(expr_data) <- "numeric"








# === Step 2. 读取 metadata 文件 ===
metadata_raw <- read.table("GSE141017_ALL_barcode_ident.csv",
                           header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 清洗 cell_barcode 字段（去除前缀 + 转换 .1 → -1）
metadata_raw$cell_barcode <- gsub("^X[_]*", "", metadata_raw$cell_barcode)
metadata_raw$cell_barcode <- gsub("\\.(\\d+)$", "-\\1", metadata_raw$cell_barcode)

# === Step 3. 特殊处理唯一重复 barcode（WT5M2_AGCTTGAGTACACCGC-1）
# 让两条记录都保留，并人为加后缀使其唯一（如 _1 和 _2）
metadata_raw <- metadata_raw %>%
  group_by(cell_barcode) %>%
  mutate(cell_barcode_unique = if (n() == 1) cell_barcode else paste0(cell_barcode, "_", row_number())) %>%
  ungroup()

# === Step 4. 构建 metadata_clean，设置行为 cell_barcode_unique
metadata_clean <- metadata_raw[, c("cell_barcode_unique", "ident", "time_point")]
rownames(metadata_clean) <- metadata_clean$cell_barcode_unique
metadata_clean$cell_barcode_unique <- NULL

# === Step 5. 同步修改 expr_data 的 column 名（确保唯一）
# 替换重复的 column barcode（这里只有 1 个重复）
colnames(expr_data)[colnames(expr_data) == "WT5M2_AGCTTGAGTACACCGC-1"] <- "WT5M2_AGCTTGAGTACACCGC-1_1"

# === Step 6. 构建 Seurat v5 对象，禁止默认过滤
seurat_obj <- CreateSeuratObject(
  counts = expr_data,
  assay = "RNA",
  meta.data = metadata_clean[colnames(expr_data), ])

# === Step 7. 添加 orig.ident 信息（从 cell barcode 中提取）

# 提取 orig.ident 为第一个非 X/下划线开头、直到第二个 "_" 之前的部分
seurat_obj$orig.ident <- gsub("^X*_*([A-Z0-9]+(?:_[A-Z0-9]+)?)_.*", "\\1", colnames(seurat_obj))



# === Step 8. 验证完整性
print(seurat_obj)
dim(LayerData(seurat_obj[["RNA"]], layer = "counts"))  # 应为 41140 x 41139
dim(seurat_obj@meta.data)                              # 应为 41139 x N
table(seurat_obj$orig.ident)
table(seurat_obj$time_point)


seurat_obj$orig.ident<-'Hm_PDAC'


save(pbmc,file = '250410.hm.Rdata')




pbmc <- NormalizeData(object = seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)


pbmc <- FindVariableFeatures(object = pbmc,selection.method = "vst", nfeatures = 2000)

all.gene <- rownames(pbmc)
pbmc <- ScaleData(pbmc,features = all.gene)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))


ElbowPlot(pbmc,ndim=50)
pbmc <- FindNeighbors(pbmc, dims = 1:20)#pa1-pc20
pbmc <- FindClusters(pbmc, resolution = 0.5)




pbmc <- RunUMAP(pbmc, dims = 1:20) 
DimPlot(pbmc, reduction = "umap",group.by = "seurat_clusters",label = T) 




#人类数据中所给基因能找到的：

features= c('SSTR1', 'SSTR2','VIPR1', 'VIPR2',  'NMUR1')
features= c('CALCRL', 'RAMP1')
DotPlot(pbmc, features = unique(features)) + RotatedAxis()+ggplot2:::coord_flip()


#小鼠数据中所给基因能找到的：

features= c('Tacr1', 'Sstr1','Sstr2','Sstr3','Vipr1', 'Vipr2')
features= c("Calcrl", "Ramp1")


DotPlot(pbmc, features = unique(features)) + RotatedAxis()+ggplot2:::coord_flip()



Idents(pbmc) <- 'seurat_clusters'
table(pbmc@meta.data$seurat_clusters)
new.cluster.ids <- c("FIB","epi", "FIB","endo", "4","5", "epi",
                     "FIB","epi", "endo","10", "peri","12", "13",
                     "epi","15", "epi","17", "FIB","19", "20",
                     "epi","22", "23","peri", "25","26")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5 )

pbmc$kinds = as.character(Idents(pbmc))
table(pbmc$kinds )
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5,group.by ='kinds' )




















load("E:/R对象/淑婷师姐/250410.mouse.Rdata")

Idents(pbmc) <- 'orig.ident'


pbmc = pbmc[, Idents(pbmc) %in% c( "CTRL2" ,  "WT15M")]


pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(object = pbmc,selection.method = "vst", nfeatures = 2000)
all.gene <- rownames(pbmc)
pbmc <- ScaleData(pbmc,features = all.gene)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
ElbowPlot(pbmc,ndim=50)
pbmc <- FindNeighbors(pbmc, dims = 1:20)#pa1-pc20
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:20) 
DimPlot(pbmc, reduction = "umap",group.by = "seurat_clusters",label = T) 
DimPlot(pbmc, reduction = "umap",group.by = "kinds",label = T) 
DimPlot(pbmc, reduction = "umap",group.by = "orig.ident",label = T) 
save(pbmc,file = '250411.俩老鼠.Rdata')








#结合http://117.50.127.228/CellMarker/

VlnPlot(pbmc, features = c("Dcn",'Fap'),slot = "counts", log = TRUE)
VlnPlot(pbmc, features = c('Pecam1','Flt1'),slot = "counts", log = TRUE)

VlnPlot(pbmc, features = c('Epcam'),slot = "counts", log = TRUE)
VlnPlot(pbmc, features = c('Krt8'),slot = "counts", log = TRUE)
VlnPlot(pbmc, features = c('Krt18'),slot = "counts", log = TRUE)


#T细胞
VlnPlot(pbmc, features = c('Cd3d'),slot = "counts", log = TRUE)


VlnPlot(pbmc, features = c('Ncam1'),slot = "counts", log = TRUE)
VlnPlot(pbmc, features = c('Fcgr3a'),slot = "counts", log = TRUE)
VlnPlot(pbmc, features = c('Gzmb'),slot = "counts", log = TRUE)


VlnPlot(pbmc, features = c('Cd14'),slot = "counts", log = TRUE)


VlnPlot(pbmc, features = c('Cd80'),slot = "counts", log = TRUE)
VlnPlot(pbmc, features = c('Cd86'),slot = "counts", log = TRUE)
VlnPlot(pbmc, features = c('Cd163'),slot = "counts", log = TRUE)


VlnPlot(pbmc, features = c('CD1C'),slot = "counts", log = TRUE)
VlnPlot(pbmc, features = c('Clec9a'),slot = "counts", log = TRUE)
VlnPlot(pbmc, features = c('Lilra4'),slot = "counts", log = TRUE)


VlnPlot(pbmc, features = c('Pdgfra'),slot = "counts", log = TRUE)

VlnPlot(pbmc, features = c('Pdgfra'),slot = "counts", log = TRUE)



VlnPlot(pbmc, features = c('Cd79a'),slot = "counts", log = TRUE)#浆细胞16
VlnPlot(pbmc, features = c('Ms4a1'),slot = "counts", log = TRUE)#B细胞12



VlnPlot(pbmc, features = c('Col3a1'),slot = "counts", log = TRUE)#成纤维3,8,9,14,19


VlnPlot(pbmc, features = c('tdTomato'),slot = "counts", log = TRUE)
VlnPlot(pbmc, features = c('Cpa1'),slot = "counts", log = TRUE)
VlnPlot(pbmc, features = c('Krt19'),slot = "counts", log = TRUE)#导管细胞2,6,7,15

#11,18，Acinar
#6,7,metaplastic,
#2,15,ductal

VlnPlot(pbmc, features = c('Chga'),slot = "counts", log = TRUE)#Endocrine，19

#巨噬marker
VlnPlot(pbmc, features = c('Ccl5'),slot = "counts", log = TRUE)
VlnPlot(pbmc, features = c('Mrc1'),slot = "counts", log = TRUE)#M2，13
VlnPlot(pbmc, features = c('Il12b'),slot = "counts", log = TRUE)#M2

VlnPlot(pbmc, features = c('Xcr1'),slot = "counts", log = TRUE)

VlnPlot(pbmc, features = c('Ramp1'),slot = "counts", log = TRUE)


clusters17.markers <- FindMarkers(pbmc, ident.1=17, logfc.threshold= 0.25,min.pct = 0.25)

top10 <- clusters17.markers  %>% top_n(n=10,wt=avg_log2FC)

#17可能是中性粒细胞



Idents(pbmc) <- 'seurat_clusters'
table(pbmc@meta.data$seurat_clusters)
new.cluster.ids <- c("Fib","endo_1", "ductal","Fib", "mono","T", "metaplastic",
                     "metaplastic","Fib", "peri","endo_2", "Acinar","B", "M2",
                     "Fib","ductal", "Plasma","Neutrophils", "Acinar","Endocrine")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5 )

pbmc$celltype = as.character(Idents(pbmc))
table(pbmc$celltype )
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5,group.by ='celltype' )
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5,group.by ='kinds' )

features= c("Calcrl", "Ramp1")


DotPlot(pbmc, features = unique(features)) + RotatedAxis()+ggplot2:::coord_flip()




table(pbmc@meta.data$orig.ident)

table(pbmc@meta.data$celltype)



# 指定目标基因
genes <- c("Calcrl", "Ramp1")
genes= c('Tacr1', 'Sstr1','Sstr2','Sstr3','Vipr1', 'Vipr2')
# 画点图，分组依据 orig.ident
DotPlot(pbmc, features = genes, group.by = "celltype") +
  RotatedAxis() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "")

DotPlot(pbmc, features = genes, group.by = "celltype") +
  RotatedAxis() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "") +
  scale_color_gradient(limits = c(-0.5, 2.5))

DotPlot(pbmc, features = genes, group.by = "celltype") +
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "") +
  scale_color_gradientn(colors = c("blue", "white", "red"),
                        limits = c(-0.5, 2.5),
                        oob = scales::squish) +
  scale_size(range = c(1, 6), 
             breaks = c(0, 20, 40, 60),
             name = "% expressed")




# 设定基因与样本
genes <- c("Calcrl", "Ramp1")
samples <- c("CTRL2", "WT15M")

# 提取表达矩阵（log-normalized）
expr <- GetAssayData(pbmc, slot = "data")[genes, ]

# 提取元信息
meta <- pbmc@meta.data[, c("orig.ident", "celltype")]

# 合并表达矩阵和元信息
df <- as.data.frame(t(expr))
df$cell_barcode <- rownames(df)
df <- cbind(df, meta)

# 转换为long格式，便于ggplot绘图
df_long <- df %>%
  pivot_longer(cols = all_of(genes), names_to = "gene", values_to = "expression") %>%
  mutate(group = paste(gene, orig.ident, sep = "_"))

# 计算每个 celltype-gene_sample 对的平均表达和表达比例
summary_df <- df_long %>%
  group_by(celltype, group) %>%
  summarise(avg_expr = mean(expression),
            pct_expr = mean(expression > 0),
            .groups = "drop")

# 画图
ggplot(summary_df, aes(x = group, y = celltype)) +
  geom_point(aes(size = pct_expr, color = avg_expr)) +
  scale_color_gradient(low = "gray90", high = "red") +
  theme_minimal() +
  labs(
    x = "Gene_Condition",
    y = "Cell Type",
    color = "Average Expression",
    size = "Percent Expressed"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))





# 基础FeaturePlot绘图
p <- FeaturePlot(pbmc, features = "Ramp1", reduction = "umap") 

# 使用主题和颜色调整美化
p + 
  scale_color_viridis_c(option = "viridis") +  # 使用viridis色系
  theme_void() +  # 去掉背景轴
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "italic")  # 标题居中斜体
  )

# 提取坐标和表达数据
umap_coords <- Embeddings(pbmc, reduction = "umap")
colnames(umap_coords) <- c("UMAP_1", "UMAP_2")  # ✅ 关键修复行！

expr_data <- FetchData(pbmc, vars = "Ramp1")

# 合并
plot_data <- cbind(umap_coords, expr_data)

# 画图
ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, fill = Ramp1)) +
  geom_point(shape = 21, size = 1.6, color = "black", stroke = 0.1) +
  scale_fill_viridis(option = "viridis", name = "Ramp1") +
  theme_void() +
  ggtitle("Ramp1") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "italic", size = 14),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )



# 基础FeaturePlot绘图
p <- FeaturePlot(pbmc, features = "Calcrl", reduction = "umap") 

# 使用主题和颜色调整美化
p + 
  scale_color_viridis_c(option = "viridis") +  # 使用viridis色系
  theme_void() +  # 去掉背景轴
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "italic")  # 标题居中斜体
  )




# 提取坐标和表达数据
umap_coords <- Embeddings(pbmc, reduction = "umap")
colnames(umap_coords) <- c("UMAP_1", "UMAP_2")

expr_data <- FetchData(pbmc, vars = "Calcrl")
plot_data <- cbind(umap_coords, expr_data)

# ⭐ 关键步骤：排序，确保高表达的点画在上方
plot_data <- plot_data[order(plot_data$Calcrl), ]

# 绘图
ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, fill = Calcrl)) +
  geom_point(shape = 21, size = 1.6, color = "black", stroke = 0.1) +
  scale_fill_viridis_c(option = "viridis", name = "Calcrl") +
  theme_void() +
  ggtitle("Calcrl") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "italic", size = 14),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )



Idents(pbmc)<-'celltype'
table(pbmc@meta.data$kinds)

VlnPlot(pbmc, features = c('Calcrl'),slot = "counts", log = TRUE)


Idents(pbmc) <- 'celltype'



#数据框操作：

pbmc = pbmc[, Idents(pbmc) %in% c( "T")]







# 提取坐标和表达数据
umap_coords <- Embeddings(pbmc, reduction = "umap")
colnames(umap_coords) <- c("UMAP_1", "UMAP_2")

expr_data <- FetchData(pbmc, vars = "Calcrl")
plot_data <- cbind(umap_coords, expr_data)

# ⭐ 关键步骤：排序，确保高表达的点画在上方
plot_data <- plot_data[order(plot_data$Calcrl), ]

# 绘图
ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, fill = Calcrl)) +
  geom_point(shape = 21, size = 1.6, color = "black", stroke = 0.1) +
  scale_fill_viridis_c(option = "viridis", name = "Calcrl") +
  theme_void() +
  ggtitle("Calcrl") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "italic", size = 14),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )











table(pbmc@meta.data$orig.ident)




library(Seurat)
library(Matrix)

# 创建 Seurat 对象（min.features 过滤表达基因数，min.cells 过滤在少于 3 个细胞中表达的基因）
pbmc <- CreateSeuratObject(counts = pbmc@assays$RNA$counts, min.cells = 3, min.features = 200, project = "PBMC")


pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")


pbmc <- subset(pbmc, 
               subset = nFeature_RNA >= 200 & nCount_RNA >= 1000 & percent.mt <= 10)


pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)


pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)


all.gene <- rownames(pbmc)

pbmc <- ScaleData(pbmc,features = all.gene)

# 后续 PCA 分析可继续进行
pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc))



ElbowPlot(pbmc,ndim=50)
pbmc <- FindNeighbors(pbmc, dims = 1:20)#pa1-pc20
pbmc <- FindClusters(pbmc, resolution = 0.5)#分辨率越小，聚类数越少；数字越大，聚类越多




pbmc <- RunUMAP(pbmc, dims = 1:20) 
DimPlot(pbmc, reduction = "umap",label = T,group.by = 'seurat_clusters') #umap更注重全局的关系


VlnPlot(pbmc, features = c("ONECUT2", "FOXQ1",'EPCAM','KRT19','MUC1'),slot = "counts", log = TRUE)



Idents(pbmc)<-'seurat_clusters'


VlnPlot(pbmc, features = c("DCN",'FAP'),slot = "counts", log = TRUE)
VlnPlot(pbmc, features = c('Pecam1','Flt1'),slot = "counts", log = TRUE)
VlnPlot(pbmc, features = c('PECAM1'),slot = "counts", log = TRUE)
VlnPlot(pbmc, features = c('EPCAM'),slot = "counts", log = TRUE)
VlnPlot(pbmc, features = c('Krt8'),slot = "counts", log = TRUE)
VlnPlot(pbmc, features = c('Krt18'),slot = "counts", log = TRUE)


#T细胞
VlnPlot(pbmc, features = c('CD3D'),slot = "counts", log = TRUE)


VlnPlot(pbmc, features = c('Ncam1'),slot = "counts", log = TRUE)
VlnPlot(pbmc, features = c('Fcgr3a'),slot = "counts", log = TRUE)
VlnPlot(pbmc, features = c('Gzmb'),slot = "counts", log = TRUE)


VlnPlot(pbmc, features = c('Cd14'),slot = "counts", log = TRUE)


VlnPlot(pbmc, features = c('Cd80'),slot = "counts", log = TRUE)
VlnPlot(pbmc, features = c('Cd86'),slot = "counts", log = TRUE)
VlnPlot(pbmc, features = c('Cd163'),slot = "counts", log = TRUE)


VlnPlot(pbmc, features = c('CD1C'),slot = "counts", log = TRUE)
VlnPlot(pbmc, features = c('Clec9a'),slot = "counts", log = TRUE)
VlnPlot(pbmc, features = c('Lilra4'),slot = "counts", log = TRUE)


VlnPlot(pbmc, features = c('Pdgfra'),slot = "counts", log = TRUE)

VlnPlot(pbmc, features = c('Pdgfra'),slot = "counts", log = TRUE)



VlnPlot(pbmc, features = c('CD79A'),slot = "counts", log = TRUE)#浆细胞16
VlnPlot(pbmc, features = c('MS4A1'),slot = "counts", log = TRUE)#B细胞12



VlnPlot(pbmc, features = c('Col3a1'),slot = "counts", log = TRUE)#成纤维3,8,9,14,19


VlnPlot(pbmc, features = c('tdTomato'),slot = "counts", log = TRUE)
VlnPlot(pbmc, features = c('CPA1'),slot = "counts", log = TRUE)
VlnPlot(pbmc, features = c('KRT19'),slot = "counts", log = TRUE)#导管细胞2,6,7,15


VlnPlot(pbmc, features = c('Chga'),slot = "counts", log = TRUE)#Endocrine，19

#巨噬marker
VlnPlot(pbmc, features = c('CCL5'),slot = "counts", log = TRUE)
VlnPlot(pbmc, features = c('MRC1'),slot = "counts", log = TRUE)#M2，13
VlnPlot(pbmc, features = c('IL12b'),slot = "counts", log = TRUE)#M2

VlnPlot(pbmc, features = c('CD86'),slot = "counts", log = TRUE)
VlnPlot(pbmc, features = c('CD163'),slot = "counts", log = TRUE)



VlnPlot(pbmc, features = c('CD4'),slot = "counts", log = TRUE)


VlnPlot(pbmc, features = c('CPA3'),slot = "counts", log = TRUE)



#13内皮，3,7,10成纤维，2,14,巨噬  1,9,12免疫    2,M1/M2 （CD86,CD163）  14巨噬？    





#人类数据中所给基因能找到的：

features= c('SSTR1', 'SSTR2','VIPR1', 'VIPR2',  'NMUR1')
features= c('CALCRL', 'RAMP1')
DotPlot(pbmc, features = unique(features)) + RotatedAxis()+ggplot2:::coord_flip()


#小鼠数据中所给基因能找到的：

features= c('Tacr1', 'Sstr1','Sstr2','Sstr3','Vipr1', 'Vipr2')
features= c("Calcrl", "Ramp1")


DotPlot(pbmc, features = unique(features)) + RotatedAxis()+ggplot2:::coord_flip()



clusters2.markers <- FindMarkers(pbmc, ident.1=2, logfc.threshold= 0.25,min.pct = 0.25)

top10 <- clusters2.markers  %>% top_n(n=10,wt=avg_log2FC)


clusters11.markers <- FindMarkers(pbmc, ident.1=11, logfc.threshold= 0.25,min.pct = 0.25)

Z.top10 <- clusters11.markers  %>% top_n(n=10,wt=avg_log2FC)


clusters13.markers <- FindMarkers(pbmc, ident.1=13, logfc.threshold= 0.25,min.pct = 0.25)

ZZ.top10 <- clusters13.markers  %>% top_n(n=10,wt=avg_log2FC)



Idents(pbmc) <- 'seurat_clusters'
table(pbmc@meta.data$seurat_clusters)
new.cluster.ids <- c("T",
                     "PDAC", "Mac","Fib", "PDAC","PDAC",
                     "T","PDAC","PDAC", "Fib","PDAC", 
                     "Mac","PDAC", "Fib","endo")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5 )

pbmc$celltype = as.character(Idents(pbmc))
table(pbmc$celltype )
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5,group.by ='celltype' )

save(pbmc,file = '人类注释好.Rdata')

features= c("CALCRL", "RAMP1")


DotPlot(pbmc, features = features, group.by = "celltype") +
  RotatedAxis() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "")


DotPlot(pbmc, features = features, group.by = "celltype") +
  RotatedAxis() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "") +
  scale_color_gradient(limits = c(-0.5, 2.5))





# 基础FeaturePlot绘图
p <- FeaturePlot(pbmc, features = "CALCRL", reduction = "umap") 

# 使用主题和颜色调整美化
p + 
  scale_color_viridis_c(option = "viridis") +  # 使用viridis色系
  theme_void() +  # 去掉背景轴
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "italic")  # 标题居中斜体
  )





# 提取坐标和表达数据
umap_coords <- Embeddings(pbmc, reduction = "umap")
colnames(umap_coords) <- c("UMAP_1", "UMAP_2")

expr_data <- FetchData(pbmc, vars = "CALCRL")
plot_data <- cbind(umap_coords, expr_data)

# ⭐ 关键步骤：排序，确保高表达的点画在上方
plot_data <- plot_data[order(plot_data$CALCRL), ]

# 绘图
ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, fill = CALCRL)) +
  geom_point(shape = 21, size = 1.6, color = "black", stroke = 0.1) +
  scale_fill_viridis_c(option = "viridis", name = "CALCRL") +
  theme_void() +
  ggtitle("CALCRL") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "italic", size = 14),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )





library(devtools)
install_github("YosefLab/VISION")

devtools::install_github("746443qjb/scKidney")





#GSE278694的h5对象的处理


library(Seurat)
library(scKidney)
library(dplyr)


# 设置工作路径
setwd("E:/R对象/淑婷师姐/GSE278694")


# 获取所有 h5 文件列表
h5_files <- list.files(pattern = "\\.h5$")

# 打印文件列表确认
print(h5_files)

# 批量读取并创建 Seurat 对象
seurat_list <- lapply(h5_files, function(f) {
  message("正在读取: ", f)
  
  # 1. 读取h5文件（10x格式）
  counts <- Read10X_h5(f)
  
  # 2. 创建Seurat对象
  obj <- CreateSeuratObject(
    counts = counts,
    project = gsub("_filtered_feature_bc_matrix.h5", "", f),
    min.cells = 3,
    min.features = 200
  )
  
  # 3. 添加样本标签
  obj$orig.ident <- gsub("_filtered_feature_bc_matrix.h5", "", f)
  
  return(obj)
})

# 合并所有Seurat对象为一个
seurat_all <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = gsub("_filtered_feature_bc_matrix.h5", "", h5_files),
  project = "GSE278694"
)

# 查看结果
seurat_all




seurat_all <- JoinLayers(seurat_all)
table(seurat_all@meta.data$orig.ident)

pbmc<-seurat_all
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)


all_genes <- rownames(pbmc)
remove_pattern <- "^MT-|^RPS|^RPL|^HB"
keep_genes <- all_genes[!grepl(remove_pattern, all_genes, ignore.case = TRUE)]
pbmc <- subset(pbmc, features = keep_genes)
dim(pbmc)

pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(object = pbmc,selection.method = "vst", nfeatures = 2000)
all.gene <- rownames(pbmc)
pbmc <- ScaleData(pbmc,features = all.gene)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
DimPlot(pbmc,reduction = "pca")
ElbowPlot(pbmc,ndim=50)
pbmc <- FindNeighbors(pbmc, dims = 1:30)#pa1-pc20
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:50) 
DimPlot(pbmc, reduction = "umap",label = T) #umap
pbmc <- RunTSNE(pbmc,dims = 1:30)
DimPlot(pbmc, reduction = "tsne") 





pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n=10,wt=avg_log2FC)

VlnPlot(pbmc, features = c("CD79A"),slot = "counts", log = TRUE)#0是B细胞

VlnPlot(pbmc, features = c("CD68"),slot = "counts", log = TRUE)

VlnPlot(pbmc, features = c("POSTN"),slot = "counts", log = TRUE)
VlnPlot(pbmc, features = c("EPCAM"),slot = "counts", log = TRUE)
Idents(pbmc)<-'seurat_clusters'
new.cluster.ids <- c("B",
                     "endocrine cells", "ductal cells","fibroblast", "myeloid cells","T/NK", 
                     "PDAC",#
                     "endocrine cells",
                     "PDAC",#
                     "endothelial cells",
                     "PDAC",#
                     "stellate cells","acinar cells", "Pericyte",
                     "PDAC")
names(new.cluster.ids) <- levels(pbmc)# 注意levels(pbmc)是啥顺序就怎么命名，和table(pbmc)顺序不一样
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

pbmc$celltype = as.character(Idents(pbmc))
table(pbmc$celltype )







features= c('POSTN')
DotPlot(pbmc, features = unique(features)) + RotatedAxis()+ggplot2:::coord_flip()







#去掉2簇

Idents(pbmc)<-'seurat_clusters'
pbmc = pbmc[, Idents(pbmc) %in% c( "0","1","3","4","5","6" ,"7","8","9","10","11" ,"12","13","14")]

pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(object = pbmc,selection.method = "vst", nfeatures = 2000)
all.gene <- rownames(pbmc)
pbmc <- ScaleData(pbmc,features = all.gene)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
DimPlot(pbmc,reduction = "pca")
ElbowPlot(pbmc,ndim=50)
pbmc <- FindNeighbors(pbmc, dims = 1:30)#pa1-pc20
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:50) 
DimPlot(pbmc, reduction = "umap",label = T) 


pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n=10,wt=avg_log2FC)

top10$gene


new.cluster.ids <- c("B",
                     "Stellate", "Myeloid","CAF", "T/NK","Ductal", 
                     "Endocrine","Endocrine","Ductal","Endothelial","Stellate",
                     "Ductal","Acinar", "Endothelial","Ductal")

names(new.cluster.ids) <- levels(pbmc)# 注意levels(pbmc)是啥顺序就怎么命名，和table(pbmc)顺序不一样
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

pbmc$celltype = as.character(Idents(pbmc))
table(pbmc$celltype )

Idents(pbmc)<-'celltype'
features= c('POSTN')
DotPlot(pbmc, features = unique(features)) + RotatedAxis()+ggplot2:::coord_flip()

features= c('ACAT2')
DotPlot(pbmc, features = unique(features)) + RotatedAxis()+ggplot2:::coord_flip()



Idents(pbmc)<-'orig.ident'
new.cluster.ids <- c("high", "high","low", "high","low", 
                     "low","high","high","low","low",
                     "high","high")

names(new.cluster.ids) <- levels(pbmc)# 注意levels(pbmc)是啥顺序就怎么命名，和table(pbmc)顺序不一样
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

pbmc$NI = as.character(Idents(pbmc))
table(pbmc$NI)






Idents(pbmc)<-'NI'
table(Idents(pbmc))

features= c('ACAT2')
DotPlot(pbmc, features = unique(features)) + RotatedAxis()+ggplot2:::coord_flip()


features= c('POSTN')
DotPlot(pbmc, features = unique(features)) + RotatedAxis()+ggplot2:::coord_flip()

features= c('POSTN','ACAT2')

DotPlot(pbmc, features = unique(features), scale = TRUE) + 
  RotatedAxis() + 
  coord_flip() +
  scale_size_continuous(range = c(4, 6)) +  # 设置圆圈的大小范围
  scale_color_gradient(
    low = "gray",  # 设置低表达的颜色
    high = "red",  # 设置高表达的颜色
    limits = c(0, 1),   # 设置色阶的范围
    breaks = seq(0, 1, by = 0.1)  # 设置色阶间隔
  )




DotPlot(pbmc, features = unique(features), scale = TRUE, group.by = 'NI') + 
  RotatedAxis() + 
  coord_flip() +
  scale_size_continuous(range = c(4, 6)) +  # 设置圆圈的大小范围
  scale_color_gradient(
    low = "yellow",  # 设置低表达的颜色
    high = "darkblue",  # 设置高表达的颜色
    limits = c(0, 1),   # 设置色阶的范围
    breaks = seq(0, 1, by = 0.1)  # 设置色阶间隔
  ) 


DotPlot(pbmc, features = unique(features), scale = TRUE, group.by = 'NI') + 
  RotatedAxis() + 
  coord_flip() +
  scale_size_continuous(range = c(4, 6)) +  # 设置圆圈的大小范围
  scale_color_viridis(option = "D", direction = 1)  # 使用 viridis 渐变色




library(MySeuratWrappers)

plot_markers <- sort(c('ACAT2', 'POSTN'))
VlnPlot(pbmc, features = plot_markers,stacked=T,pt.size=0.5)





DotPlot(pbmc, features = unique(features), scale = TRUE) + 
  RotatedAxis() + 
  coord_flip() +
  scale_size_continuous(range = c(4, 6)) +  # 设置圆圈的大小范围
  scale_color_gradient(
    low = "yellow",  # 设置低表达的颜色
    high = "darkblue",  # 设置高表达的颜色
    limits = c(0, 1),   # 设置色阶的范围
    breaks = seq(0, 1, by = 0.1)  # 设置色阶间隔
  )



DotPlot(pbmc, features = unique(features), scale = TRUE) + 
  RotatedAxis() + 
  coord_flip() +
  scale_size_continuous(range = c(4, 6)) +  # 设置圆圈的大小范围
  scale_color_gradient(
    low = "yellow",  # 设置低表达的颜色
    high = "darkblue",  # 设置高表达的颜色
    limits = c(0, 1),   # 设置色阶的范围
    breaks = seq(0, 1, by = 0.1)  # 设置色阶间隔
  ) +
  guides(color = guide_colorbar(title = "Average Expression"))




DotPlot(pbmc, features = unique(features), scale = TRUE) + 
  RotatedAxis() + 
  coord_flip() +
  guides(
    size = guide_legend(title = "Percent Expressed"),  # 调整圆圈大小的图例
    color = guide_colorbar(title = "Average Expression")  # 调整颜色的图例
  )








































































































































































































































































































































































































































