#A script to analyse RNA featurecounts gene count matrix, comparing this to a bedfile of likely L1 locations, finding unique expression and looking for correlations in age and neuron type
#
# featurecounts_GRCh37d5_ERCC92 is a gene count matrix
# overlap2.bed is a bed file containing likely L1 genes as output from BEDtools
# import these files into working directory
# allbrains_celltype_annot_with_ANN_predictions.tsv a tabe seperated file containing classification of the cells 

install.packages("splitstackshape")
library(splitstackshape) #import splitstackshape
library(tidyverse) #import tidyverse
library(ggplot2) #import ggplot2
library(gridExtra) #import gridExtra
require("devtools") #import devtools
install_github('tallulandrews/M3Drop') #install M3drop
library(M3Drop) # import M3 drop
library(corrplot) # import corrplot
library(ltm) # import ltm

# Read data in to R

feature = read.table("Downloads/featurecounts_GRCh37d5_ERCC92.txt", header = TRUE, check.names = FALSE)

# Separate out columns that contain more than one value

chrom = as.data.frame(feature$Chr)
names(chrom) <- c("Chr")
startn = as.data.frame(feature$Start)
names(startn) <- c("Start")
finish = as.data.frame(feature$End)
names(finish) <- c("End")
strand = as.data.frame(feature$Strand)
names(strand) <- c("Strand")

chrom <- cSplit(chrom, "Chr", ";")
startn <- cSplit(startn, "Start", ";")
finish <- cSplit(finish, "End", ";")
strand <- cSplit(strand, "Strand", ";")

#finding the largest region the gene could cover

startnt <- as.data.frame(rowMins(data.matrix(startn), na.rm = TRUE, dim. = dim(startn))) #min start site
names(startnt) <-c("Min")
finish <- as.data.frame(rowMaxs(data.matrix(finish), na.rm = TRUE, dim. = dim(finish)))#max finsih site
names(finish) <-c("Max")

# Create bed file with largest distances for genes

new_df <- cbind.data.frame(chrom$Chr_0001, startnt$Min, finish$Max, strand$Strand_0001, feature$Geneid)
names(new_df) <- c("Chromosome", "Start", "End", "Strand", "Name")

# write out bedfile

write.table(new_df, file = "expressed_genes.bed", quote = FALSE, sep = "\t")


# Read in overlap data from BEDtools run in command line

filterbed<-as.data.frame(read.table("overlap2.bed", header = FALSE, sep = "\t"))
geneid<-as.vector(filterbed$V10) # identify overlapped gene names from bedfile
unigene <- unique(geneid) #ensure genes are unique

feature2 <- feature[,-1] #remove geneid to make rownames
rownames(feature2) <- feature[,1] # make geneid rownames
feature21 <- feature2[6:578] # isolate count data
shortnames<-regmatches(names, gregexpr("\\d{4}_[[:upper:]]\\d{2}", names)) # Get well identification to use as column name
names(feature21)<-shortnames # use shortnames for column names
l1mat <- feature21[unigene,] # isolate L1 crossover genes
filterno<- feature21[!(rownames(feature21) %in% unigene),] # isolate genes not known L1 insertion location

# make a list of spike in identifiers

spike <- regmatches(rownames(feature21), regexpr("ERCC-[0-9]{5}", rownames(feature21)))

# use regex to create list of cells from each brain using well IDs (column names)

brain19 <- regmatches(shortnames, regexpr("2440_[[:upper:]]\\d{2}", shortnames)) #brain 19
brain22 <- regmatches(shortnames, regexpr("244[18]_[[:upper:]]\\d{2}", shortnames)) #brain 22
brain21 <- regmatches(shortnames, regexpr("2444_[[:upper:]]\\d{2}", shortnames)) # brain 21
brain20 <- regmatches(shortnames, regexpr("2445_[[:upper:]]\\d{2}", shortnames)) # brain 20
brain23 <- regmatches(shortnames, regexpr("2449_[[:upper:]]\\d{2}", shortnames)) # brain 23

#create L1 brain matrices

l1brain <- function(brainlist, L1_matrix){
  brainl1<-L1_matrix[,brainlist]
  brainl12 <- brainl1[rowSums(brainl1)>0,]
  return(brainl12)
}

#L1 brain matrices for each brain

brain19l1 <- l1brain(brain19, l1mat)
brain22l1 <- l1brain(brain22, l1mat)
brain21l1 <- l1brain(brain21, l1mat)
brain20l1 <- l1brain(brain20, l1mat)
brain23l1 <- l1brain(brain23, l1mat)

#separate entire gene count matrix for each brain

brain19all<-feature21[,brain19]
brain22all<-feature21[,brain22]
brain21all<-feature21[,brain21]
brain20all<-feature21[,brain20]
brain23all<-feature21[,brain23]

#genes that are not identified as known L1 locations

brain19no<-filterno[,brain19]
brain22no<-filterno[,brain22]
brain21no<-filterno[,brain21]
brain20no<-filterno[,brain20]
brain23no<-filterno[,brain23]

# a function to plot gene distributions in each brain for L1 genes, not L1 genes and the whole gene matrix from the relevant gene count matrices

plot_dist_graph <- function(brain19mat, brain22mat, brain21mat, brain20mat, brain23mat){
  col0a <- as.data.frame(colSums(brain19mat))
  names(col0a)<-c("Total_count")
  p0<-ggplot(col0a, aes(x=Total_count)) + geom_histogram(aes(y=..density..))+geom_density(alpha=.2, fill="#FF6666")+labs(title="Brain 19",x="Total gene count", y = "Density")
  col1a <- as.data.frame(colSums(brain22mat))
  names(col1a)<-c("Total_count")
  p1<-ggplot(col1a, aes(x=Total_count)) + geom_histogram(aes(y=..density..))+geom_density(alpha=.2, fill="#FF6666")+labs(title="Brain 22",x="Total gene count", y = "Density")
  col2a <- as.data.frame(colSums(brain21mat))
  names(col2a)<-c("Total_count")
  p2<-ggplot(col2a, aes(x=Total_count)) + geom_histogram(aes(y=..density..))+geom_density(alpha=.2, fill="#FF6666")+labs(title="Brain 21",x="Total gene count", y = "Density")
  col3a <- as.data.frame(colSums(brain20mat))
  names(col3a)<-c("Total_count")
  p3<-ggplot(col3a, aes(x=Total_count)) + geom_histogram(aes(y=..density..))+geom_density(alpha=.2, fill="#FF6666")+labs(title="Brain 20",x="Total gene count", y = "Density")
  col5a <- as.data.frame(colSums(brain23mat))
  names(col5a)<-c("Total_count")
  p4<-ggplot(col5a, aes(x=Total_count)) + geom_histogram(aes(y=..density..))+geom_density(alpha=.2, fill="#FF6666")+labs(title="Brain 23",x="Total gene count", y = "Density")
  ouput<-list(p0=p0, p1=p1, p2=p2, p3=p3,p4=p4)
  return(output)
  }

#plot histograms for all brains showing gene distribution

allbrain<-plot_dist_graph(brain19all, brain22all, brain21all, brain20all, brain23all)
grid.arrange(allbrain$p0,allbrain$p3,allbrain$p2,allbrain$p1,allbrain$p4, top = "Histogram plots to show gene count distribution for all genes in neurons for each brain")

l1brain<-plot_dist_graph(brain19l1, brain22l1, brain21l1, brain20l1, brain23l1)
grid.arrange(l1brain$p0,l1brain$p3,l1brain$p2,l1brain$p1,l1brain$p4, top = "Histogram plots to show gene count distribution for l1 genes in neurons for each brain")

nobrain<-plot_dist_graph(brain19no, brain22no, brain21no, brain20no, brain23no)
grid.arrange(nobrain$p0,nobrain$p3,nobrain$p2,nobrain$p1,nobrain$p4, top = "Histogram plots to show gene count distribution for all genes in neurons for each brain")

########finding unique expression in each brain nd plotting gene expression in unique genes######

# a function to identify top expression in each brain and elimnate gene with no expression

top_expression <- function(L1brainmat){
  topbrain <- L1brainmat[,colSums(L1brainmat)>=quantile(colSums(L1brainmat), c(.95))]
  topbrain <- topbrain[rowSums(topbrain)>0,]
  topbrain <- topbrain[rowSums(topbrain)>1000,]
  return(topbrain)
}

top19 <- top_expression(brain19l1)
top22 <- top_expression(brain22l1)
top21 <- top_expression(brain21l1)
top20 <- top_expression(brain20l1)
top23 <- top_expression(brain23l1)

#a function to identify the unique genes highly expressed for each brain

unique_genes <- function(target_brainmatL1, brainmatL11, brainmatL12,brainmatL13,brainmatL14,){
  brainuni <- setdiff(rownames(arget_brainmatL1), rownames(brainl11))
  brainuni <- setdiff(brainuni, rownames(brainl12))
  brainuni <- setdiff(brainuni, rownames(brainl13))
  brainuni <- setdiff(brainuni, rownames(brainl14))
  return(brainuni)
}

brain19uni<-unique_genes(brain19l1, brain22l1, brain21l1, brain20l1, brain23l1)
brain22uni<-unique_genes(brain22l1, brain19l1, brain21l1, brain20l1, brain23l1)
brain21uni<-unique_genes(brain21l1, brain22l1, brain19l1, brain20l1, brain23l1)
brain20uni<-unique_genes(brain20l1, brain22l1, brain21l1, brain19l1, brain23l1)
brain23uni<-unique_genes(brain23l1, brain22l1, brain21l1, brain20l1, brain19l1)

# a function to create a matrix of the unique genes for each brain

top_matrix <- function(L1brainmat, unique_list){
  topbrain <- L1brainmat[rowSums(L1brainmat)>0,colSums(L1brainmat)>0]
  topbrain <- topbrain[unique_list,]
  return(topbrain)
}

topbrainmat19 <- top_matrix(brain19l1,brain19uni)
topbrainmat22 <- top_matrix(brain22l1,brain22uni)
topbrainmat21 <- top_matrix(brain21l1,brain21uni)
topbrainmat20 <- top_matrix(brain20l1,brain20uni)
topbrainmat23 <- top_matrix(brain23l1,brain23uni)

# a function to use the top brain matrix and find the total counts for each gene and how many cells express the gene

topgene <- function(top_brain_matl){
  Total<- rowSums(top_brain_matl);
  Numberofcells <- rowSums(top_brain_matl != 0);
  top_brain_matl1<-cbind(top_brain_matl,Total);
  top_brain_matl2<- cbind(top_brain_matl1, Numberofcells)
  return(top_brain_matl2)
}

tot19 <- topgene(topbrainmat19)
tot22 <- topgene(topbrainmat22)
tot21 <- topgene(topbrainmat21)
tot20 <- topgene(topbrainmat20)
tot23 <- topgene(topbrainmat23)

#a function to prepare plots of unique expression

unique_plot<- function(tot){
p <- ggplot(data=tot, aes(x=rownames(tot), y=Total)) +
  geom_bar(stat="identity", width=0.5) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
  return(p)
}

plot19<-unique_plot(tot19)
plot22<-unique_plot(tot22)
plot21<-unique_plot(tot21)
plot20<-unique_plot(tot20)
plot23<-unique_plot(tot23)

#assign labels to plots and plot

plot19 + labs(title="Unique L1 candidate genes brain 19", x="Unique genes", y = "Total count")
plot22 + labs(title="Unique L1 candidate genes brain 22", x="Unique genes", y = "Total count")
plot21 + labs(title="Unique L1 candidate genes brain 21", x="Unique genes", y = "Total count")
plot20 + labs(title="Unique L1 candidate genes brain 20", x="Unique genes", y = "Total count")
plot23 + labs(title="Unique L1 candidate genes brain 23", x="Unique genes", y = "Total count")


################################ M3 drop + Brennecke ##################

#prepare data with normalisation to create expression matrix 
expr_mat <- M3DropConvertData(feature21, is.counts = TRUE)
expr_mat <- as.data.frame(expr_mat)
par(mfrow = c(1,1))

spike <- regmatches(rownames(feature21), regexpr("ERCC-[0-9]{5}", rownames(feature21))) #identify spike ins from original dataframe usinf regex

#get highly variable gene from Brennecke analysis using the expression matrix and spike ins identified earlier

Brennecke_HVG <- BrenneckeGetVariableGenes(expr_mat, spikes = spike,fdr = 0.01, minBiolDisp = 0.5)

HVG_genes <- Brennecke_HVG$Gene # isolate highly variable genes for comparison

#get drop out genes with M3drop using the expression matrix

M3Drop_genes <- M3DropFeatureSelection(
  expr_mat,
  mt_method = "fdr",
  mt_threshold = 0.01
)

M3Drop_genes <- M3Drop_genes$Gene # isolate drop out genes
dropspike <- regmatches(M3Drop_genes, regexpr("ERCC-[0-9]{5}", M3Drop_genes)) # identify spike ins in drop out genes


M3drop_genesnospike <- setdiff(M3Drop_genes, dropspike)

################################################ correlation analysis ################################################

# read in cell type annotation
cell_type <- read.table("Documents/allbrains_celltype_annot_with_ANN_predictions.tsv", header = TRUE, check.names = FALSE)

# identify inhibitory cells
inhib_cell_type <- cell_type[cell_type$Neuronal_type == "Inhibitory",]

#identify ecitatory cells
excib_cell_type <- cell_type[cell_type$Neuronal_type == "Excitatory",]

#create a list of cell names according to classification

inhib <- inhib_cell_type$SampleName 
exci <- excib_cell_type$SampleName

#get cell names for each classification using regex
platinhib <- regmatches(inhib, regexpr("244\\d_[[:upper:]]\\d{2}", inhib))
platexci <- regmatches(exci, regexpr("244\\d_[[:upper:]]\\d{2}", exci))


# define cells into age categories according to the brain they originated from
"<50" <- regmatches(shortnames, regexpr("244[18]_[[:upper:]]\\d{2}", shortnames))
"50-60" <- regmatches(shortnames, regexpr("244[49]_[[:upper:]]\\d{2}", shortnames))
">60" <- regmatches(shortnames, regexpr("244[50]_[[:upper:]]\\d{2}", shortnames))




#create dataframe to run correlations on and log normalise the data
corr_plot_l1 <- function(matrix, inhib, exci, age1, age2, age3){
  l1matcut <- matrix[rowSums(matrix)>1000,]
  l1len <- feature2[rownames(l1matcut),]
  dataframs <- cbind(rowSums(l1matcut), rowSums(l1matcut[,inhib]), rowSums(l1matcut[,exci]),rowSums(l1matcut[,`<50`]), rowSums(l1matcut[,`50-60`]), rowSums(l1matcut[,`>60`]))
  colnames(dataframs) <- c("All_Neurons","Inhibitory","Excitatory", "<50", "50-60", ">60")
  lgdataframs<- log10(dataframs + 1)
  sclgdataframs <- scale(lgdataframs)
  correl <- cor(sclgdataframs)
  return(correl)
}


#plot L1 correlation plot and not L1 correlation plot alongside each other
par(mfrow=c(1,2))
nol1cor <-corr_plot_l1(nol1mat,platinhib, platexci,`<50`,`50-60`,`>60`)
corrplot(nol1cor, type = "upper", title = "Not L1 genes correlation", diag = FALSE,tl.col = "black", cl.lim = c(0, 1))
l1cor <- corr_plot_l1(l1mat,platinhib, platexci,`<50`,`50-60`,`>60`)
corrplot(l1cor, type = "upper", title = "L1 genes correlation", diag = FALSE,tl.col = "black", cl.lim = c(0, 1))



#create dataframe to compare length of gene to liklihood of L1 insertion
notl1names<- !(rownames(feature2) %in% unigene)
lenmat<- as.data.frame(as.numeric(as.character(feature$Length)), row.names = rownames(feature2))
lenmat<- cbind(lenmat, as.integer(as.logical(notl1names)))
colnames(lenmat) <- c("Length", "L1")

attach(lenmat)
biserial.cor(Length, L1)
cor.test(Length, L1) #run correlation test

#plot boxplot showing gene length for L1 genes and non-L1 genes
lenmat$L1<-as.factor(lenmat$L1)
fill <- "#4271AE"
ggplot(lenmat, aes(x=L1, y=Length)) + geom_boxplot(fill=fill) +
  scale_y_continuous(trans='log2')

#plot total unique genes for each brain in a scatter plot
uni_lengths <- c(length(brain0uni),length(brain1uni),length(brain2uni),length(brain3uni),length(brain5uni))
age_brain <- c(65,41,53,72,59)
age_length <- cbind(age_brain, uni_lengths)
age_length<-as.data.frame(age_length)
ggplot(age_length, aes(x=age_brain, y=uni_lengths)) + geom_point() + labs(x="Brain age", y="Number of unique genes")
