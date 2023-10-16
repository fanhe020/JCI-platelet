library("Seurat");
library("sctransform");
library("dplyr");
library("RColorBrewer");
library("ggthemes");
library("ggplot2");
library("cowplot");
library("data.table");
library(patchwork)
#basic processing
{samples = c("HI105","HI108","HI111","PT130","PT139","PT153","PT158","PT164")
samples_dataset = c("HI","HI","HI","ET","ET","ET","ET","ET")
samples_mut = c("N","N","N","J","C","C","J","J")
samples_Throm = c("N","N","N","Y","N","N","N","Y")
samples_PLTcount = c("L","L","L","H","H","H","L","L")
data.10x = list(); # first declare an empty list in which to hold the feature-barcode matrices
data.10x[[1]] <- Read10X(data.dir ="D:/PBMC_scRNA-seq/rawdata/HI105");
data.10x[[2]] <- Read10X(data.dir = "D:/PBMC_scRNA-seq/rawdata/HI108");
data.10x[[3]] <- Read10X(data.dir = "D:/PBMC_scRNA-seq/rawdata/HI111");
data.10x[[4]] <- Read10X(data.dir ="D:/PBMC_scRNA-seq/rawdata/PT130");
data.10x[[5]] <- Read10X(data.dir = "D:/PBMC_scRNA-seq/rawdata/PT139");
data.10x[[6]] <- Read10X(data.dir = "D:/PBMC_scRNA-seq/rawdata/PT153");
data.10x[[7]] <- Read10X(data.dir = "D:/PBMC_scRNA-seq/rawdata/PT158");
data.10x[[8]] <- Read10X(data.dir = "D:/PBMC_scRNA-seq/rawdata/PT164");

scrna.list = list(); # First create an empty list to hold the Seurat objects
for (i in 1:length(data.10x)) {
  scrna.list[[i]] = CreateSeuratObject(counts = data.10x[[i]], min.cells=100, min.features=700, project=samples[i]);
  scrna.list[[i]][["DataSet"]] = samples_dataset[i];
  scrna.list[[i]][["Mutation"]] = samples_mut[i];
  scrna.list[[i]][["Thrombosis"]] = samples_Throm[i];
  scrna.list[[i]][["PLTcount"]] = samples_PLTcount[i];
  
}

rm(data.10x);

scrna <- merge(x=scrna.list[[1]], y=c(scrna.list[[2]],scrna.list[[3]],scrna.list[[4]],scrna.list[[5]],scrna.list[[6]],scrna.list[[7]],scrna.list[[8]]), project="ET_PBMC");
rm(scrna.list); # save some space
str(scrna@meta.data) # examine the structure of the Seurat object meta data

scrna[[]];
scrna@meta.data;
str(scrna@meta.data); # Examine structure and contents of meta data
head(scrna@meta.data$nFeature_RNA); # Access genes ("Features") for each cell
head(scrna@meta.data$nCount_RNA); # Access number of UMIs for each cell:
levels(x=scrna); # List the items in the current default cell identity class
length(unique(scrna@meta.data$seurat_clusters)); # How many clusters are there? Note that there will not be any clusters in the meta data until you perform clustering.
unique(scrna@meta.data$PLTcount); # What  are included in this data set?
scrna$NewIdentity <- vector_of_annotations; # Assign new cell annotations to a new "identity class" in the meta data


jpeg("VlnPlot_total.jpg", width = 20, height = 8, units="in", res=300);
scrna[["percent.mt"]] <- PercentageFeatureSet(scrna, pattern = "^MT-")
scrna[["percent.ribo"]] <- PercentageFeatureSet(scrna, pattern = "^RP[SL][[:digit:]]")
vln <-VlnPlot(scrna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo"), ncol = 4,pt.size=0)
print(vln)
dev.off();

# can draw the VlnPlot in different grouping: Mutation, DataSet, Thrombosis, PLTcount
VlnPlot(scrna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo"), group.by = 'DataSet', ncol = 4,pt.size=0)
# check how many cells before filtering
table(scrna@meta.data$orig.ident)
nrow(scrna@meta.data)
Count95 <- quantile(scrna@meta.data$nCount_RNA, 0.95) # calculate value in the 95rd percentile
Count95


# filter cells with 3 parameters setting, I changed the mito ratio to <0.2 instead of 0.1 in the protocol
scrna.filt <- subset(x = scrna, subset = nFeature_RNA > 700  & nCount_RNA < Count95 & percent.mt < 20)
nrow(scrna.filt@meta.data)
table(scrna.filt@meta.data$orig.ident)
save(scrna.filt, file = "scrna.filt.RData")
load("scrna.filt.RData")

# plot vlnplot after filtering
jpeg("VlnPlot_filt.jpg", width = 20, height = 8, units="in", res=300);
vln <-VlnPlot(scrna.filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo"), ncol = 4,pt.size=0)
print(vln)
dev.off();


# my RAM is not large enough to hold 7.5GB vector during scaling, so I can only sample first
subcells <- sample(Cells(scrna.filt), size=40000, replace=F)
scrna.sub <- subset(scrna.filt, cells=subcells)
table(scrna.sub@meta.data$orig.ident)

# normalize data
scrna.norm <- NormalizeData(object = scrna.sub, normalization.method = "LogNormalize", scale.factor = 1e5);

# Identification of highly variable features (feature selection)
scrna.norm <- FindVariableFeatures(scrna.norm, selection.method = "vst",nfeatures = 2000)

saveRDS(scrna.norm,"norm_ET_PBMC.rds")

scrna.norm<-readRDS("norm_ET_PBMC.rds")

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(scrna.norm), 10)


# plot variable features with and without labels, nice results: PPBP, PF4 are top VGs
jpeg("Variable_features.jpg", width = 20, height = 8, units="in", res=300);
plot1 <- VariableFeaturePlot(scrna.norm)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE) # add label
print(plot1 + plot2)
dev.off()

# scaling data: a standard pre-processing step prior to dimensional reduction techniques like PCA
scrna.norm.scaled <- ScaleData(object = scrna.norm, features = rownames(x = scrna.norm))

# check scaled data: GetAssayData(scrna.norm.scaled,slot="scale.data",assay="RNA")[1:8,1:4] 
# Perform linear dimensional reduction
pbmc <- RunPCA(scrna.norm.scaled, features = VariableFeatures(object = scrna.norm.scaled))

jpeg("Dimplot.jpg", width = 8, height = 8, units="in", res=300)
Dim<-DimPlot(pbmc, reduction = "pca")
print(Dim)
dev.off()
DimPlot(pbmc, reduction = "pca", group.by = 'DataSet')
DimPlot(pbmc, reduction = "pca", group.by = 'Thrombosis')


jpeg("PCA.heatmap.jpg", width = 8.5, height = 24, units="in", res=300);
hm <- DimHeatmap(object = pbmc, dims = 1:12, cells = 500, balanced = TRUE);
print(hm);
dev.off();

elbow <- ElbowPlot(object = pbmc)
jpeg(sprintf("PCA.elbow.jpg"), width = 6, height = 8, units="in", res=300);
print(elbow);
dev.off();
saveRDS(pbmc,"pbmc.rds")

pbmc<-readRDS("pbmc.rds")

pbmc <- JackStraw(object = pbmc, dims = 20, num.replicate = 100); # takes around 4 minutes
pbmc <- ScoreJackStraw(object = pbmc, dims = 1:20)

jpeg("PCA.jackstraw.jpg", width = 10, height = 10, units="in", res=300);
js <- JackStrawPlot(object = pbmc, dims = 1:20)
print(js);
dev.off();

pc.pval <- pbmc@reductions$pca@jackstraw@overall.p.values; # get p-value for each PC
write.csv(pc.pval, file="PCA.jackstraw.scores.csv")

# Clustering cells, dimension and resolution setting determines the number of clusters
pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Look at cluster IDs of the first 5 cells, it seems 28 clusters were created
table(Idents(pbmc))

# Run non-linear dimensional reduction (UMAP/tSNE)
nPC = 20;
pbmc <- RunTSNE(object = pbmc, dims = 1:nPC);
pbmc <- RunUMAP(object = pbmc, dims = 1:nPC);

saveRDS(pbmc, "clustered_ET_PBMC.rds")
rm(list=ls())
# start from here
pbmc<-readRDS("clustered_ET_PBMC.rds")
str(pbmc@meta.data)

# Finding differentially expressed features (cluster biomarkers)
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) %>%
write.csv(file="cluster.markers.csv")

# Assigning cell type identity to clusters
new.cluster.ids <- c("CD4", "Mono", "CD4", "NK", "CD8", "CD4","CD8","B",
                     "B","Mono","PLT","Mono","B","Mono","B","PLT",
                     "CD8", "PLT","HSC", "Mono", "HSC")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
str(pbmc@meta.data)
pbmc[["celltype"]] <- Idents(pbmc)
str(pbmc@meta.data)
saveRDS(pbmc,"renamed_clustered_ET_PBMC.rds")
}

rm(list=ls())
pbmc<-readRDS("renamed_clustered_PBMC.rds")

# general plots
{tiff("UMAP.tiff", units="in", width=4, height=3, res=300)
ru<-DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5,label.size = 3)
print(ru)
dev.off()



tiff("UMAP_dataset.tiff", units="in", width=4, height=3, res=300)
ru<-DimPlot(pbmc, reduction = "umap", group.by = 'DataSet', pt.size = 0.5,label.size = 3)
print(ru)
dev.off()

# how manually change names of specific clusters
#scrna <-RenameIdents(scrna, "17" = "plt")
#scrna[["celltype"]] <- Idents(scrna)


table(pbmc@meta.data$celltype,pbmc@meta.data$DataSet)
table(pbmc@meta.data$celltype,pbmc@meta.data$Mutation)
M=table(pbmc@meta.data$celltype,pbmc@meta.data$orig.ident)
write.csv(M, "M.csv")


# different types of plots : VlnPlot(), FeaturePlot(), RidgePlot() and CellScatter()
VlnPlot(pbmc, features = c("LGALS1","LGALS2","LGALS3","LGALS4"),pt.size = 0.2,split.by="DataSet" )
VlnPlot(pbmc, features =  c("nFeature_RNA", "nCount_RNA", "percent.mt"),pt.size = 0)

cellmarkers<- c("CD79A", "CCR7", "CD8A", "GNLY","LYZ", "CD14", "PPBP","PF4","CD34")

tiff("Featureplot.tiff", units="in", width=7.5, height=6, res=300)
fp <-FeaturePlot(pbmc, features = cellmarkers)
print(fp);
dev.off();

tiff("vlnplot.tiff", units="in", res=300, height=6, width=8);
vp <-VlnPlot(pbmc, pt.size=0,features = 'LGALS1',split.by = 'DataSet')
print(vp);
dev.off();

FeaturePlot(pbmc, features = c("LGALS1","LGALS2","LGALS3","LGALS4"))

tiff("Dotplot.tiff", units="in", res=300, height=6, width=6);
dt<-DotPlot(pbmc, features = cellmarkers, dot.scale = 6, cols= "BrBG",group.by= "celltype", split.by = "DataSet") +
  RotatedAxis() + theme(axis.text.x=element_text(angle = 90, hjust = 0,vjust=0.5))
print(dt)
dev.off()
}

# add a list of gene sets as an index for plotting

{#inflammation response
list_a =read.delim("H_inflam.txt")
list_a=list_a[-1,]
pbmc <- AddModuleScore(pbmc,
                       features = list(list_a),
                       name="inflammation")
str(pbmc@meta.data)

tiff("Inflammation_1.tiff", units="in", res=300, height=5, width=5);
p = FeaturePlot(pbmc,
            features = "inflammation1", label = TRUE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
print(p)
dev.off()
tiff("Inflammation_2.tiff", units="in", res=300, height=3, width=6);
p = VlnPlot(pbmc,features = "inflammation1",pt.size=0, split.by ='Mutation')
print(p)
dev.off()

# IFN gamma 
list_b =read.delim("H_IFN_g.txt")

list_b=list_b[-1,]

pbmc <- AddModuleScore(pbmc,
                       features = list(list_b),
                       name="IFNg")
str(pbmc@meta.data)
tiff("IFNg_1.tiff", units="in", res=300, height=5, width=5);
p = FeaturePlot(pbmc,
                features = "IFNg1", label = TRUE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
print(p)
dev.off()
tiff("IFNg_2.tiff", units="in", res=300, height=3, width=6);
p = VlnPlot(pbmc,features = "IFNg1",pt.size=0, split.by ='Mutation')
print(p)
dev.off()

# TNF NKFB
list_c =read.delim("H_TNFA_NFKB.txt")

list_c=list_c[-1,]

pbmc <- AddModuleScore(pbmc,
                       features = list(list_c),
                       name="HNFKB")
str(pbmc@meta.data)
tiff("HNFKB_1.tiff", units="in", res=300, height=5, width=5);
p = FeaturePlot(pbmc,
                features = "HNFKB1", label = TRUE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
print(p)
dev.off()
tiff("HNFKB_2.tiff", units="in", res=300, height=3, width=6);
p = VlnPlot(pbmc,features = "HNFKB1",pt.size=0, split.by ='Mutation')
print(p)
dev.off()



# IL6 JAK2 STAT3
list_d =read.delim("H_IL6_JAK_STAT3.txt")

list_d=list_d[-1,]

pbmc <- AddModuleScore(pbmc,
                       features = list(list_d),
                       name="HJAKSTAT")
str(pbmc@meta.data)
tiff("JAKSTAT_1.tiff", units="in", res=300, height=5, width=5);
p = FeaturePlot(pbmc,
            features = "HJAKSTAT1", label = TRUE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
print(p)
dev.off()
tiff("JAKSTAT_2.tiff", units="in", res=300, height=5, width=5);
p = VlnPlot(pbmc,features = "HJAKSTAT1",pt.size=0, split.by ='orig.ident')
print(p)
dev.off()

# T cell activation gene set
list_t =read.delim("POSITIVE_REGULATION_OF_T_CELL_ACTIVATION.txt")
list_t=list_t[-1,]

pbmc <- AddModuleScore(pbmc,
                       features = list(list_t),
                       name="Tcellmarker")
str(pbmc@meta.data)
tiff("T_activation_1.tiff", units="in", res=300, height=5, width=5);
p = FeaturePlot(pbmc,
                features = "Tcellmarker1", label = TRUE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
print(p)
dev.off()
tiff("T_activation_2.tiff", units="in", res=300, height=5, width=5);
p = VlnPlot(pbmc,features = "Tcellmarker1",pt.size=0)
print(p)
dev.off()
}
rm(list=ls())
pbmc<-readRDS("renamed_clustered_PBMC.rds")

tiff("GALECTINS_PBMC.tiff", units="in", res=300, height=5, width=6);
p1 = VlnPlot(pbmc,features = "LGALS1",pt.size=0, split.by ='DataSet')
p2 = VlnPlot(pbmc,features = "LGALS3",pt.size=0, split.by ='DataSet')
p3 = FeaturePlot(pbmc,features = "LGALS1")
p4 = FeaturePlot(pbmc,features = "LGALS3")
print((p1/p2)|(p3/p4))
dev.off()

tiff("GALECTINS_DOTPLOT.tiff", units="in", res=300, height=4, width=4.5);
LGALS = c("LGALS1","LGALS2","LGALS3","LGALS4","LGALS8","LGALS9","LGALS12")
p=DotPlot(pbmc, features = LGALS, dot.scale = 6, cols= "BrBG", split.by = "DataSet") +
  RotatedAxis() + theme(axis.text.x=element_text(angle = 90, hjust = 0,vjust=0.5))
print(p)
dev.off()

FeatureScatter(object = pbmc, feature1 = 'LGALS1', feature2 = 'inflammation1')

FeatureScatter(object = pbmc, feature1 = 'LGALS1', feature2 = 'IFNg1')

FeatureScatter(object = pbmc, feature1 = 'LGALS1', feature2 = 'HNFKB1')
FeatureScatter(object = pbmc, feature1 = 'LGALS1', feature2 = 'HJAKSTAT1')


Mono<- subset(x = pbmc, idents = "Mono")
saveRDS(Mono,"Mono_PBMC.rds")
# mono
libiray(dplyr)
rm(list=ls())
Mono <- readRDS("Mono_PBMC.rds")
Mono <- NormalizeData(Mono, normalization.method = "LogNormalize",
                     scale.factor = 1e4) 
Mono <- FindVariableFeatures(Mono, selection.method = 'vst',
                            nfeatures = 2000)
Mono <- ScaleData(Mono, features = rownames(x = Mono))
Mono <- RunPCA(Mono, features = VariableFeatures(object = Mono)) 
Mono <- FindNeighbors(Mono, dims = 1:10)
Mono <- FindClusters(Mono, resolution = 0.8 )
Mono <- RunTSNE(object = Mono, dims = 1:10);
Mono <- RunUMAP(object = Mono, dims = 1:10);
library(ggplot2)
library(stringr)
library("fgsea")
library('presto')
Mono.genes <- wilcoxauc(Mono, "DataSet")


tiff("Mono_UMAP_1.tiff", units="in", res=300, height=3, width=3.5);
p <- DimPlot(object = Mono, reduction = "umap",  pt.size=1,group.by="DataSet")
print(p)
dev.off()

tiff("Mono_UMAP_2.tiff", units="in", res=300, height=3, width=3.5);
p <- DimPlot(object = Mono, reduction = "umap", pt.size=1, label = T)
print(p)
dev.off()



str(PLT@meta.data)
tiff("Mono_inf_3.tiff", units="in", res=300, height=3, width=6);
p =VlnPlot(Mono,features="PPBP", pt.size = 0)
print(p)
dev.off()


p1 = FeaturePlot(Mono,
                features = "inflammation1", label = TRUE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

p2 = VlnPlot(Mono,features = "inflammation1",pt.size=0.2, split.by = 'DataSet')

p3 =VlnPlot(Mono,features="PPBP", pt.size = 0)


tiff("Mono_IFNg_1.tiff", units="in", res=300, height=3, width=3.5);
p = FeaturePlot(Mono,
            features = "IFNg1", label = TRUE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
print(p)
dev.off()

VlnPlot(Mono,features = "IFNg1",pt.size=0.2, split.by = 'DataSet')


VlnPlot(Mono,features="PPBP", pt.size = 0)




tiff("PAG_2.tiff", units="in", res=300, height=5, width=5);
p = VlnPlot(Mono,features = "IFNg1",pt.size=0)
print(p)
dev.off()

p=DotPlot(Mono, features = "PPBP", dot.scale = 6, cols= "BrBG",group.by= "RNA_snn_res.0.8", split.by = "DataSet") +
  RotatedAxis() + theme(axis.text.x=element_text(angle = 90, hjust = 0,vjust=0.5))
print(p)





VlnPlot(Mono,features="PPBP", pt.size = 0)
Monomarkers= FindAllMarkers(PLT,only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
table(Idents(PLT))
new.mono.ids <- c("ET", "ET", "ET", "HI", "HI", "HI","B","B",
                     "B","B","B","B","B")
names(new.mono.ids) <- levels(PLT)
PLT <- RenameIdents(PLT, new.mono.ids)

Monomarker012= Monosubclustermarkers %>%
  filter(cluster %in% c(0,1,2))
Monomarker345= Monomarkers %>%
  filter(cluster %in% c(3,4,5))
write.csv(Monomarker345, "Monomarker345.csv")

  #save(file="allcelltype.markers.top10.Rdata")


golist=read.delim("go0019221.txt")

golist=golist[,2]

golisttogo=golist[golist %in% Monomarker012$gene]

jpeg("Monodotplotcytokine-mediated signaling pathway.jpg", units="in", res=300, height=4, width=8);
DotPlot(pbmc, features = golisttogo, dot.scale = 8, cols= "BrBG",group.by= "celltype", split.by = "DataSet") +
  RotatedAxis()
dev.off()


DotPlot(PLT,features= "CXCL1")
write.csv(golist, "mono012markersincytokineresponselist0019221.csv")
"SPI1" %in% Monomarkers$gene
VlnPlot(PLT,"NFKBIA")+
  geom_boxplot()

FeaturePlot(PLT,"NFKBIA")

Hlist=read.delim("IFN gamma gene list.txt")

Hlist=Hlist[-1,]

Hlisttogo=Hlist[Hlist %in% commongenes]
library(ggplot2)
DotPlot(PLT,features=Hlisttogo, split.by = 'DataSet') + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_colour_gradient2(low = "#FF00FF", mid = "#000000", high = "#FFFF00")
str(pbmc@meta.data$celltype)
Idents(pbmc)
VlnPlot(pbmc,'CD14',group.by = 'celltype',split.by = 'DataSet')

genes <- Reduce(intersect, list(A, B, C, D，E))





