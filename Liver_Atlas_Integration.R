#All data for this project can be found at the following links
# Aizarani - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE124395
# Macparland - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115469
# Ramachandran - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE136103
# Segal - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130473
# Tamburini - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129933


#Libraries
library(ggplot2)
library(dplyr)
library(Seurat)
library(patchwork)

###Read in data from each of the primary papers and do initial QC and filtering
##Aizarani
data<-readRDS('GSE124395_Normalhumanliverdata.RData')
aizarani_healthy<-CreateSeuratObject(data, project="Aizarani_Healthy", min.cells=3, min.features = 300)
VlnPlot(aizarani_healthy, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
aizarani_healthy <- subset(aizarani_healthy, subset = nFeature_RNA > 300 & nFeature_RNA < 2500 & percent.mt < 30)

##MacParland
healthy1_data<-Read10X(data.dir = 'Healthy1/')
healthy1<-CreateSeuratObject(counts=healthy1_data,project = "Healthy1", min.cells = 3, min.features = 300)
healthy2_data<-Read10X(data.dir = 'Healthy2/')
healthy2<-CreateSeuratObject(counts=healthy2_data,project = "Healthy2", min.cells = 3, min.features = 300)
healthy3_data<-Read10X(data.dir = 'Healthy3/')
healthy3<-CreateSeuratObject(counts=healthy3_data,project = "Healthy3", min.cells = 3, min.features = 300)
healthy4_data<-Read10X(data.dir = 'Healthy4/')
healthy4<-CreateSeuratObject(counts=healthy4_data,project = "Healthy4", min.cells = 3, min.features = 300)
healthy5_data<-Read10X(data.dir = 'Healthy5/')
healthy5<-CreateSeuratObject(counts=healthy5_data,project = "Healthy5", min.cells = 3, min.features = 300)
macparland_healthy <- merge(healthy1, y=c(healthy2, healthy3, healthy4, healthy5), add.cell.ids = c("healthy1","healthy2","healthy3","healthy4","healthy5"),project="Macparland_Healthy")
macparland_healthy[["percent.mt"]] <- PercentageFeatureSet(macparland_healthy, pattern = "^MT-")
VlnPlot(macparland_healthy, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
macparland_healthy <- subset(macparland_healthy, subset = nFeature_RNA > 300 & nFeature_RNA < 2500 & percent.mt < 30)

##Ramachandran
healthy1a_data<-Read10X(data.dir = 'Healthy1A/')
healthy1a<-CreateSeuratObject(counts=healthy1a_data,project = "Healthy1a", min.cells = 3, min.features = 300)
healthy1b_data<-Read10X(data.dir = 'Healthy1B/')
healthy1b<-CreateSeuratObject(counts=healthy1b_data,project = "Healthy1b", min.cells = 3, min.features = 300)
healthy2a_data<-Read10X(data.dir = 'Healthy2/')
healthy2a<-CreateSeuratObject(counts=healthy2a_data,project = "Healthy2a", min.cells = 3, min.features = 300)
healthy3a_data<-Read10X(data.dir = 'Healthy3A/')
healthy3a<-CreateSeuratObject(counts=healthy3a_data,project = "Healthy3a", min.cells = 3, min.features = 300)
healthy3b_data<-Read10X(data.dir = 'Healthy3B/')
healthy3b<-CreateSeuratObject(counts=healthy3b_data,project = "Healthy3b", min.cells = 3, min.features = 300)
healthy4a_data<-Read10X(data.dir = 'Healthy4/')
healthy4a<-CreateSeuratObject(counts=healthy4a_data,project = "Healthy4a", min.cells = 3, min.features = 300)
ramachandran_healthy <- merge(healthy1a, y=c(healthy1b, healthy2a, healthy3a, healthy3b, healthy4a), add.cell.ids = c("healthy1a","healthy1b","healthy2","healthy3a","healthy3b","healthy4"),project="Ramachandran_Healthy")
ramachandran_healthy[["percent.mt"]] <- PercentageFeatureSet(ramachandran_healthy, pattern = "^MT-")
VlnPlot(ramachandran_healthy, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ramachandran_healthy <- subset(ramachandran_healthy, subset = nFeature_RNA > 300 & nFeature_RNA < 2500 & percent.mt < 30)

##Segal
#Import in counts and metadata
segal_count<-read.csv("GSE130473_Series_count_matrix.csv", header=T, row.names = 1)
segal_meta<-read.csv("GSE130473_Series_feature_matrix.csv", header = T)
#Filter out the adult samples
segal_adult<-subset(segal_meta, segal_meta$TISSUE=="ADULT")
segal_count_adult<-segal_count[segal_adult$X]
#Change row names from ensembl gene ids to gene symbols (Note this takes 15 min)
encode_ids<-row.names(segal_count_adult)
encode_ids <- gsub('\\..+$', '', encode_ids)
library(EnsDb.Hsapiens.v79)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= encode_ids, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
encode_ids<-as.data.frame(encode_ids)
geneIDs2<-merge(encode_ids,geneIDs1, by.y="GENEID", by.x="encode_ids", all=T)
geneIDs2$SYMBOL <- ifelse(is.na(geneIDs2$SYMBOL), geneIDs2$encode_ids, geneIDs2$SYMBOL)
segal_count_adult<-as.matrix(segal_count_adult)
row.names(segal_count_adult)<-geneIDs2$SYMBOL
segal_healthy<-CreateSeuratObject(segal_count_adult, project="segal_Healthy", min.cells=3, min.features = 300)
segal_healthy[["percent.mt"]] <- PercentageFeatureSet(segal_healthy, pattern = "^MT-")
VlnPlot(segal_healthy, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
segal_healthy <- subset(segal_healthy, subset = nFeature_RNA > 300 & nFeature_RNA < 2500 & percent.mt < 30)

##Tamburini
data<-read.table("GSE129933_count_matrix.tsv", header=T)
row.names(data)<-data[,1]
data<-data[,-1]
tamburini_healthy<-CreateSeuratObject(data, project="Tamburini_Healthy", min.cells=3, min.features = 300)
tamburini_healthy <-subset(tamburini_healthy, orig.ident=="nd1"|orig.ident=="nd2")
tamburini_healthy[["percent.mt"]] <- PercentageFeatureSet(tamburini_healthy, pattern = "^MT-")
VlnPlot(tamburini_healthy, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
tamburini_healthy <- subset(tamburini_healthy, subset = nFeature_RNA > 300 & nFeature_RNA < 2500 & percent.mt < 30)

###Normalize all data
aizarani_healthy <- NormalizeData(aizarani_healthy)
macparland_healthy <- NormalizeData(macparland_healthy)
ramachandran_healthy <- NormalizeData(ramachandran_healthy)
segal_healthy <- NormalizeData(segal_healthy)
tamburini_healthy <- NormalizeData(tamburini_healthy)

###Find Variable Features
aizarani_healthy <- FindVariableFeatures(aizarani_healthy, selection.method = "vst", nfeatures = 2000)
macparland_healthy <- FindVariableFeatures(macparland_healthy, selection.method = "vst", nfeatures = 2000)
ramachandran_healthy <- FindVariableFeatures(ramachandran_healthy, selection.method = "vst", nfeatures = 2000)
segal_healthy <- FindVariableFeatures(segal_healthy, selection.method = "vst", nfeatures = 2000)
tamburini_healthy <- FindVariableFeatures(tamburini_healthy, selection.method = "vst", nfeatures = 2000)

###Integration
liver_list=c(aizarani_healthy,macparland_healthy,ramachandran_healthy,segal_healthy,tamburini_healthy)
k.filter <- min(200, min(sapply(liver_list, ncol)))
liver.anchors <- FindIntegrationAnchors(object.list = liver_list, dims = 1:30, k.filter = k.filter)
liver.integrated <- IntegrateData(anchorset = liver.anchors, dims = 1:30)

# switch to integrated assay. 
DefaultAssay(object = liver.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
liver.integrated <- ScaleData(object = liver.integrated, verbose = FALSE)
liver.integrated <- RunPCA(object = liver.integrated, npcs = 50, verbose = FALSE)
liver.integrated <- RunUMAP(object = liver.integrated, reduction = "pca", 
                            dims = 1:50)
#Cluster the cells
liver.integrated <- FindNeighbors(liver.integrated, dims = 1:50)
liver.integrated <- FindClusters(liver.integrated, resolution = 0.8)
head(Idents(liver.integrated), 5)
DimPlot(liver.integrated, reduction = "umap",label=T, pt.size = 0.5) + NoLegend()
paper<-c(rep("Aizarani",11868), rep("Macparland",12476), rep("Ramachandran",11570),rep("Segal",94),rep("Tamburini",180))
liver.integrated$paper<-paper
DimPlot(object = liver.integrated, reduction = "umap", group.by = "paper")

cluster.markers <- FindMarkers(liver.integrated, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
head(cluster.markers, n = 10)
new.cluster.ids <- c("NK,NKT,T Cells","NK,NKT,T Cells","Hepatocytes","Endothelial Cells","Mononuclear Phagocytes",
                     "Innate Lymphoid Cells","Cholangiocytes","Mononuclear Phagocytes","Hepatocytes",
                     "Endothelial Cells","Innate Lymphoid Cells","Endothelial Cells","Hepatocytes",
                     "Endothelial Cells","Plasma Cells","Hepatocytes","Mesenchymal Cells","Mesenchymal Cells",
                     "Hepatocytes","Endothelial Cells","Undetermined","B Cells","Cycling","Cholangiocytes","pDCs")

names(new.cluster.ids) <- levels(liver.integrated)
liver.integrated <- RenameIdents(liver.integrated, new.cluster.ids)
DimPlot(liver.integrated, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
