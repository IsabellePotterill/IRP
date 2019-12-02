#a script to filter xTEA results and compare number of candidates in boxplots

#set working directory to xTEA_analysis
#this folder contains text files of the ouput of xTEA 'candidate_disc_filtered_cns.txt' for each cell and for the bulk data


library(tidyverse) #import tidyverse
library(gridExtra) #import gridExtra

list.filenames<-list.files(pattern=".txt$") #get file names
list.data<-list() #initiate large list

#a for loop to read each file into the large list
for (i in 1:length(list.filenames))
{
  list.data[[i]]<-read.delim(list.filenames[i], header=FALSE, row.names=NULL)
}

#assign each file in the larg list its correspoding name from above
names(list.data)<-list.filenames

#name columns of each dataframe and filter results
for (i in 1:length(list.filenames))
{
  colnames(list.data[[i]])<-c("Chr","refined-pos","lclip-pos", "rclip-pos", "TSD", "nalclip","narclip", "naldisc", "nardisc", "transducction", "side", "consensus")
  list.data[[i]] <- list.data[[i]][list.data[[i]]$side == "two_side_tprt_both",]
  list.data[[i]] <- list.data[[i]][list.data[[i]]$consensus == "hit_end_of_consensus",]
}


# ensure the left clip read is lower than the right clipped read and if not switch them around to create bed file. Output bedfile.

for (i in 1:length(list.filenames)){
  name <- paste((list.filenames[i]),".bed", sep = "")
  Data <- list.data[[i]]
  Data[which(Data$`rclip-pos` < Data$`lclip-pos`), 3:4] <- Data[which(Data$`rclip-pos` < Data$`lclip-pos`), c(4,3)] 
  df <- data.frame(Data$Chr, Data$`lclip-pos`, Data$`rclip-pos`)
  write.table(df,name,quote = FALSE,append = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

#get bulk data file names
bulk.filenames<-list.files(pattern="bulk.txt$")
bulk.data<-list()

#import bulk data files to separate large list
for (i in 1:length(bulk.filenames))
{
  bulk.data[[i]]<-read.delim(bulk.filenames[i], header=FALSE, row.names=NULL)
}
#name each file in the large list
names(bulk.data)<-bulk.filenames
#name the columns of each dataframe in the large list
for (i in 1:length(bulk.filenames))
{
  colnames(bulk.data[[i]])<-c("Chr","refined-pos","lclip-pos", "rclip-pos", "TSD", "nalclip","narclip", "naldisc", "nardisc", "transducction", "side", "consensus")
}

#create vectors of the number of rows of each dataframe in the original imported data
vec <-as.data.frame(lapply(orig.data,nrow))
vec <-as.data.frame(t(vec))
#create vectors of the number of rows of each dataframe in the cut down data imported data
list_vec <-as.data.frame(lapply(list.data,nrow))
list_vec<-as.data.frame(t(list_vec))
#create vectors of the number of rows of each dataframe in the original bulk data
bulk_vec <-as.data.frame(lapply(bulk.data,nrow))
bulk_vec<-as.data.frame(t(bulk_vec))

#bind individual vectors to cerate a dataframe representing xTEA results before and after cut offs
data_size_vector<-rbind(vec,bulk_vec)
data_size_vector<-cbind(data_size_vector, list_vec$V1)
colnames(data_size_vector)<-c("candidates", "cut_off_candidates")

#get final data from bed files
final.filenames<-list.files(pattern=".txt.bed$")
final.data<-list()
#import final data from bed file crossover
for (i in 1:length(final.filenames))
{
  
  final.data[[i]]<-read.delim(final.filenames[i], header=FALSE, row.names=NULL, blank.lines.skip = FALSE)
}

names(final.data)<-final.filenames


#ensure files are not empty
agevec <- vector()
for (file in list.files(,"*.txt.bed$")){
  if (file.size(file) == 0) {
    a<-(file.size(file))
  }else{
    Table <- read.table(file)
    a<-(as.integer(nrow(Table)))
  }
  agevec <- append(agevec, a)
}
agevec<- append(agevec, rep(NA, 5))

#manipulate data into form for plotting
data_size_vector<-cbind(data_size_vector, agevec)
data_size_vector<-cbind(data_size_vector, brain_age)
brain_number<-c(19,19,19,19,19,19,22,22,22,22,22,21,21,21,21,20,20,20,23,23,23,23,23,23,22,22,22,22,22,23,19,20,21,22)
data_size_vector<-cbind(data_size_vector,brain_number)
data_size_vector<-data_size_vector[,c("brain_number","candidates","cut_off_candidates","agevec")]
age1<-c(65, 65, 65, 65, 65, 65, 41, 41, 41, 41, 41, 53, 53, 53, 53, 72, 72, 72, 59, 59, 59, 59, 59, 59, 41, 41, 41, 41, 41, 59, 65, 72, 53, 41)
data_size_vector<-cbind(data_size_vector, age1)
data_size_vector<-data_size_vector[,c("brain_number","age","candidates","cut_off_candidates","crossover_candidates", "age1")]

#plot 3 boxplots showing the relative number of candidates
data_size_vector$age<-as.factor(data_size_vector$age)
g <- ggplot(data_size_vector, aes(x=age, y=data_size_vector$candidates, fill=age)) + geom_boxplot() + scale_fill_brewer(palette="Dark2") + labs(Title="a)", x="Age", y = "Number of candidates")
h <- ggplot(data_size_vector, aes(x=age, y=data_size_vector$cut_off_candidates, fill = age)) + geom_boxplot() + scale_fill_brewer(palette="Dark2") + labs(Title="b)", x="Age", y = "Number of candidates")
k <- ggplot(data_size_vector, aes(x=age, y=data_size_vector$crossover_candidates, fill = age)) + geom_boxplot() + scale_fill_brewer(palette="Dark2") + labs(Title="c)", x="Age", y = "Number of candidates")
grid.arrange(g,h,k, ncol=3)
