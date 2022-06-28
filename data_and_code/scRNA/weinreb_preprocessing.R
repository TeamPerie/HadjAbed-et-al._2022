
rm(list = ls())
library(Seurat)
library(Matrix)
library(clustree)

setwd("/Users/jasoncosgrove/Dropbox (Team_Perie)/Jason/Experiments/Dry_Lab/JCB12_Weinreb/my_analysis/in_vitro")

matrix <- readMM("stateFate_inVitro_normed_counts.mtx")
gene.names <- read.table("stateFate_inVitro_gene_names.txt",header = F)
meta.data <- read.table("stateFate_inVitro_metadata.txt",header = T)
clones <- readMM("stateFate_inVitro_clone_matrix.mtx")
#neu.mono.traj <- read.table("stateFate_inVitro_neutrophil_monocyte_trajectory.txt",header = T)


#create a unique id for each celltype and update the dimnames of hte matrix
meta.data$cell_id <- paste(meta.data$Cell_barcode,meta.data$Library,sep = "_")
dimnames(matrix) <- list(meta.data$cell_id,gene.names$V1)

rownames(meta.data) <- meta.data$cell_id

#process the clone information
clones <- as.matrix(clones)
rownames(clones) <- meta.data$cell_id
colnames(clones) <- 1:ncol(clones)


#convert the celltype to a factor which will make life easier later on. 
meta.data$Cell_type_annotation <- as.factor(meta.data$Cell_type_annotation)

#for each clone we want a probability associated with each fate
clonefates <- matrix(0, ncol = 11, nrow = 5864)
colnames(clonefates) <- levels(meta.data$Cell_type_annotation)


#for each clone, get all of the rows which are non-zero and find the celltypes
#of interest, store the % of each celltype per clone in the clonefates matrix
# also update each cells metadata file to state whcih clone they are associated with. 
meta.data$clone_id <- NA

for(i in 1:ncol(clones)){
  
  #subset the metadata to get only the cells that belong to this clone
  meta.data.subset <- meta.data[meta.data$cell_id %in% names(clones[,i][clones[,i] == 1]),]
  #for each cell update the clonal information
  meta.data[meta.data$cell_id %in% names(clones[,i][clones[,i] == 1]),"clone_id"] <- i
  #save the proportion of each celltype associated with each clone to the clonefates matrix
  clonefates[i,] <- table(meta.data.subset$Cell_type_annotation) / sum(table(meta.data.subset$Cell_type_annotation))
}


#assign a fate to each clone based on which outcome has the highest probability
clone.bias <- colnames(clonefates)[apply(clonefates,1,which.max)]


###now we have done all of the basic processing lets create our seurat object
sobj <- CreateSeuratObject(counts = t(matrix), project = "weinreb", min.cells = 100, min.features = 100,meta.data = meta.data)

#add the spring coordinates as a dimension reduction object
spring <- as.matrix(meta.data[,7:8])
colnames(spring) <- paste0("SPRING_", 1:2)
sobj[["spring"]] <- CreateDimReducObject(embeddings = spring, key = "SPRING_", assay = DefaultAssay(sobj))



output <- matrix("no_clone",ncol = 1, nrow = nrow(sobj@meta.data))
rownames(output) <- rownames(sobj@meta.data)
colnames(output) <- "clone_bias"
for(cell in 1:nrow(sobj@meta.data)){
  #get teh clone id for that cell
  lin <- clone.bias[sobj@meta.data$clone_id[cell]][[1]]
  if(!is.na(lin)){  output[cell,1] <- lin }
}

sobj@meta.data$clone_fate <- output
sobj@meta.data$is.barcoded <- sobj@meta.data$clone_fate != "no_clone"
DimPlot(sobj, reduction = "spring", group.by = "is.barcoded")

cols = c("grey","grey","grey","grey","grey","grey",
         "grey","blue","red","grey","grey","grey")



tiff("/users/jasoncosgrove/Desktop/spring.tiff",units="in", width=6, height=4, res=300)


cols = c("dark blue","#56B4E9","#0072B2",'#F15854',"#2d543d",
         "#B276B2","#CB2027",
         "#617a89","#0b53c1","#60BD68","grey90","#117733")
DimPlot(sobj, reduction = "spring", group.by = "Cell_type_annotation",
        cols = cols,raster = FALSE)
dev.off()


cols = c("white","white","white","white","white","white",
         "white", "blue","red","white","grey80")



sobj.neu.mo <- subset(sobj,(clone_fate ==  "Monocyte" | clone_fate ==  "Neutrophil" | clone_fate ==  "Undifferentiated"))
tiff("/users/jasoncosgrove/Desktop/spring.tiff",units="in", width=6, height=4, res=300)

mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(15)
DimPlot(sobj.neu.mo, reduction = "spring", group.by = "Cell_type_annotation",raster = FALSE,cols = mycolors)
dev.off()

tiff("/users/jasoncosgrove/Desktop/spring2.tiff",units="in", width=6, height=4, res=300)
DimPlot(sobj, reduction = "spring", group.by = "orig.ident",
        cols = c("grey"),raster = FALSE)
dev.off()

cols = c("grey","orange","grey","grey","grey","grey",
         "grey","blue","red","orange","grey","grey")
DimPlot(sobj, reduction = "spring", group.by = "Cell_type_annotation",cols = cols)

clone.fate.simplified <- sobj@meta.data$clone_fate

clone.fate.simplified[clone.fate.simplified == 'Baso'] <- "Myeloid"
clone.fate.simplified[clone.fate.simplified == 'Ccr7_DC'] <- "Myeloid"
clone.fate.simplified[clone.fate.simplified == 'Eos'] <- "Myeloid"
clone.fate.simplified[clone.fate.simplified == 'Mast'] <- "Myeloid"
clone.fate.simplified[clone.fate.simplified == 'Monocyte'] <- "Myeloid"
clone.fate.simplified[clone.fate.simplified == 'Neutrophil'] <- "Myeloid"

sobj@meta.data$clone_fate_simplified <- clone.fate.simplified

cols = c("red","blue","green","orange","light grey","black","dark grey")
DimPlot(sobj, reduction = "spring", group.by = "clone_fate_simplified",cols = cols)


sobj.d2 <- subset(sobj,subset = Time_point == 2)



DimPlot(sobj.d2, reduction = "spring", group.by = "clone_fate_simplified",cols = cols)

load("genesets/geneset_Robjects/metabolic_signatures_DRAG.Rda")

sobj.d2 <- AddModuleScore(sobj.d2, features = list(metabolic.signatures))




lsk <- subset(sobj,Starting_population == "Lin-Kit+Sca1+")
lk <- subset(sobj,Starting_population == "Lin-Kit+Sca1-")

df.lsk <- table(lsk@meta.data$clone_id,
lsk@meta.data$Cell_type_annotation)

df.lk <- table(lk@meta.data$clone_id,
                lk@meta.data$Cell_type_annotation)


write.csv(df.lsk,"/Users/jasoncosgrove/Desktop/weinreb_barcodes_lsk.csv")
write.csv(df.lk,"/Users/jasoncosgrove/Desktop/weinreb_barcodes_lk.csv")


intersect(lsk@meta.data$clone_id,lk@meta.data$clone_id)
