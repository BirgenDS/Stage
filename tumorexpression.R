library("zoo")
library("gplots")
library("tidyr")

source("http://www.bioconductor.org/biocLite.R")
biocLite("limma")
biocLite("edgeR")

setwd("/Users/Birgen/Documents/Howest/Stage/CMGG/RNA-seq_CNV")

source("eSNPKaryotyping/R/CreatVCF.R")
source("eSNPKaryotyping/R/DeletionTable.R")
source("eSNPKaryotyping/R/EditVCF.R")
source("eSNPKaryotyping/R/Edit_dbSNP_Files.R")
source("eSNPKaryotyping/R/Sort_major_minor.R")

# wrapperscript
output_dir <- "/Users/Birgen/Documents/Howest/Stage/CMGG/RNA-seq_CNV/KIRC/TCGA-CZ-5467-01A-01R-1503-07_eSNPout/"
input_dir <- "/Users/Birgen/Documents/Howest/Stage/CMGG/RNA-seq_CNV/"
cancer_type <- "KIRC"
Organism = "Human"
cnv_dir <- "/Users/Birgen/Documents/Howest/Stage/CMGG/RNA-seq_CNV/TCGA_CNVdata"
Ylim = 3

setwd("/Users/Birgen/Documents/Howest/Stage/CMGG/RNA-seq_CNV")

wrapper <- function(output_dir, input_dir, cancer_type, cnv_dir){
  table = read.delim(paste(output_dir,"variantTable.csv",sep=""))
  table$chr=as.numeric(table$chr)
  table=table[order(table$chr,table$position),]
  table=table[table$chr>0,]
  return(table)
}



# MajorMinor calculation
Table = wrapper(output_dir, input_dir, cancer_type, cnv_dir)

MajorMinor <- function(Table,minDP,maxDP,minAF){
  Table[is.na(Table)] = 0
  newTable = Table[Table$DP >= minDP,]
  newTable = newTable[newTable$DP <= maxDP,]
  AF1 = newTable$AD1/newTable$DP
  AF2 = newTable$AD2/newTable$DP
  newTable = data.frame("chr" = newTable$chr, "position" = newTable$position, "AF1" = AF1, "AF2" = AF2)
  frequncyTable = newTable[newTable$AF1 >= minAF,]
  frequncyTable = frequncyTable[frequncyTable$AF2 >= minAF,]
  orderedTable = Sort_major_minor(data=frequncyTable, col1=3, col2=4)
  MajorMinor = orderedTable$AF1/orderedTable$AF2
  orderedTable["MajorMinor"] = MajorMinor
  return(orderedTable)
}



# Plotgenome
table2 = MajorMinor(Table = Table,minDP = 20,maxDP = 1000,minAF = 0.2)
orderedTable = table2

PlotGenome<-function(Output_Directory,orderedTable,Window,Ylim,Organism,CancerType,CNV_Directory){
  newDir = substring(Output_Directory,1)
  setwd(newDir)
  options(bitmapType='cairo')
  
  window=Window
  #moving median
  
  if(Organism == "Human"){
    centromere_pos=c(125,93.3,91,50.4,48.4,61,59.9,45.6,49,40.2,53.7,35.8,17.9,17.6,19,36.6,24,17.2,26.5,27.5,13.2,14.7,60.6,12.5)
    centromere_pos=centromere_pos*1000000
    chr_size=c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,155270560,59373566)
  }
  
  chr_total=0
  for (i in 1:(length(chr_size)-1)){
    chr_total=c(chr_total,sum(chr_size[1:i]))  
  }
  
  # chromosome 23 and 24 not in plot
  genome_size=sum(chr_size[1:22])
  
  orderedTable$position=orderedTable$position+chr_total[orderedTable$chr]
  
  dirs <- list.dirs(path = CNV_Directory, full.names = F, recursive = F)
  selected_dir <- ''
  for (dir in dirs){
    pos <- regexpr("_",dir)[1]
    ss <- substr(dir,0,as.numeric(pos)-1)
    if (ss == CancerType){
      selected_dir = dir # lange versie
    }
  }
  file_manifest <- read.delim(paste(paste(CNV_Directory,selected_dir,sep="/"),"file_manifest.txt",sep="/"))
  sf <- substr(newDir,as.numeric(regexpr("TCGA-",newDir)[1]),as.numeric(regexpr("TCGA-",newDir)[1])+14)
  selected <- subset(file_manifest, grepl("\\.hg19\\.seg\\.txt", File.Name)  &  Sample == sf ,select = File.Name)
  selected[] <- lapply(selected, as.character)
  
  cnv_file = paste(paste(paste(CNV_Directory,selected_dir,sep="/"),"CNV_SNP_Array/BI__Genome_Wide_SNP_6/Level_3",sep="/"),selected[1,1],sep="/")
  cnv = read.delim(cnv_file)
  cnv$Chromosome = as.character(cnv$Chromosome)
  cnv$Chromosome[cnv$Chromosome=='X'] = 23
  cnv$Chromosome[cnv$Chromosome=='Y'] = 24
  cnv$Chromosome = as.numeric(cnv$Chromosome)
  
  cnv_new <- t(apply(cnv,1,function(x) {
    chr_huidig<-as.numeric(x[2])
    new_start<-as.numeric(x[3]) + chr_total[chr_huidig]
    new_stop<-as.numeric(x[4]) + chr_total[chr_huidig]
    return(c(as.numeric(x[2]),as.numeric(new_start),as.numeric(new_stop),as.numeric(x[6]),as.numeric(x[5])))
  }))
  cnv_new <- as.data.frame(cnv_new)
  colnames(cnv_new) <- c('Chromosome', 'Start', 'Stop', 'Segment_Mean','Num_Probes')
  starts = cnv_new[,c('Chromosome', 'Start', 'Segment_Mean','Num_Probes')]
  colnames(starts) <- c('Chromosome', 'Position', 'Segment_Mean','Num_Probes')
  ends = cnv_new[,c('Chromosome', 'Stop', 'Segment_Mean','Num_Probes')]
  colnames(ends) <- c('Chromosome', 'Position', 'Segment_Mean','Num_Probes')
  
  cnv_lg<-as.data.frame(rbind(starts, ends)) # plakt end kolom onder start
  cnv_lg <- cnv_lg[order(cnv_lg[,2]),]
  cnv_lg$Score = 0
  cnv_lg$Score[cnv_lg$Segment_Mean>0.2] = 1
  cnv_lg$Score[cnv_lg$Segment_Mean< (-0.2)] = -1
  cnv_lg <- cnv_lg[cnv_lg$Num_Probes>100,]
  
  orderedTable3 <- orderedTable[orderedTable$chr != 23, ]
  orderedTable <- orderedTable3[orderedTable3$chr != 24, ]
  cnv_lg2 <- cnv_lg[cnv_lg$Chromosome != 23, ] 
  cnv_lg <- cnv_lg2[cnv_lg2$Chromosome != 24, ]
  
  # Table with values from first plot
  orderedTableMedian <- rollmedian(orderedTable$MajorMinor,window, fill = NA)
  orderedTablePostition <- rollmedian(orderedTable$position,window, fill = NA)
  orderedTableRoll = data.frame(orderedTable$chr,orderedTablePostition, orderedTableMedian)
  colnames(orderedTableRoll) <- c('Chromosome','Position', 'Rollmedian')
  orderedTableRoll <- na.omit(orderedTableRoll)
  
  # table with rolling median > 1.75, in stead of MajorMinor
  orderedTableMedian <- rollmedian(orderedTable$MajorMinor,window, fill = NA)
  orderedTable2 = data.frame(orderedTable$chr,orderedTable$position, orderedTableMedian)
  colnames(orderedTable2) <- c('Chromosome','Position', 'Rollmedian')
  orderedTable2 <- orderedTable2[orderedTable2$Rollmedian >= 1.75, ] # aangepast van 2 naar 1.75 om meer pieken over te houden
  orderedTable2 <- data.frame(orderedTable2$Chromosome,orderedTable2$Position, orderedTable2$Rollmedian)
  colnames(orderedTable2) <- c('Chromosome','Position', 'Rollmedian')
  orderedTableMedian <- na.omit(orderedTable2)
  #write.csv(orderedTableMedian, file = 'ordereTableMedian.csv', row.names=FALSE)
  
  chr_total2 <- chr_total[-1];
  
  # loop to select different imballances
  orderedTableMedian2 <- data.frame(Chromosome=orderedTableMedian$Chromosome[1], Start=orderedTableMedian$Position[1], Stop=c(NA), Rollmedian=c(NA))
  rollmedianlist <- c(orderedTableMedian$Rollmedian[1])
  #i = 2
  for (i in 2:nrow(orderedTableMedian)){
    if(orderedTableMedian$Position[i] > chr_total2[orderedTableMedian$Chromosome[i-1]]){
      orderedTableMedian2$Stop[is.na(orderedTableMedian2$Stop)] = orderedTableMedian$Position[i-1]
      orderedTableMedian2$Rollmedian[is.na(orderedTableMedian2$Rollmedian)] = mean(rollmedianlist)
      newerrow = c(orderedTableMedian$Chromosome[i], orderedTableMedian$Position[i],orderedTableMedian2$Stop[NA], orderedTableMedian2$Rollmedian[NA])
      orderedTableMedian2 <- rbind(orderedTableMedian2,newerrow)
      rollmedianlist <- c()
    } else if(orderedTableMedian$Position[i] - orderedTableMedian$Position[i-1] > 10000000) {
      orderedTableMedian2$Stop[is.na(orderedTableMedian2$Stop)] = orderedTableMedian$Position[i-1]
      orderedTableMedian2$Rollmedian[is.na(orderedTableMedian2$Rollmedian)] = mean(rollmedianlist) # 1 vervangen door rollmedian van start op deze lijn
      newrow = c(orderedTableMedian$Chromosome[i], orderedTableMedian$Position[i],orderedTableMedian2$Stop[NA], orderedTableMedian2$Rollmedian[NA])
      orderedTableMedian2 <- rbind(orderedTableMedian2,newrow)
      rollmedianlist <- c()
    }
    rollmedianlist <- append(rollmedianlist,orderedTableMedian$Rollmedian[i])
    #print(rollmedianlist)
    i = i+1
  }
  
  setwd(output_dir)
  orderedTableMedian2$Stop[is.na(orderedTableMedian2$Stop)] = orderedTableMedian$Position[i-1]
  orderedTableMedian2$Rollmedian[is.na(orderedTableMedian2$Rollmedian)] = mean(rollmedianlist)
  orderedTableMedian2 <- orderedTableMedian2[(orderedTableMedian2$Stop - orderedTableMedian2$Start) > 1000, ]
  orderedTableMedian2 <- orderedTableMedian2[orderedTableMedian2$Chromosome != 23, ]
  orderedTableMedian2 <- orderedTableMedian2[orderedTableMedian2$Chromosome != 24, ]
  write.table(orderedTableMedian2[,1:3], file = 'Regions_new.csv', sep = "\t", row.names=FALSE, quote = FALSE)
  return(list(cnv_lg = cnv_lg, orderedTable = orderedTable, orderedTableMedian2 = orderedTableMedian2[,1:3], chr_total2 = chr_total2))
}



Plotgenomecnv <- PlotGenome(Output_Directory = output_dir,orderedTable = table2,Window = 151,Ylim = 3,Organism = "Human", CancerType = cancer_type, CNV_Directory = cnv_dir)
cnv_lg = Plotgenomecnv$cnv_lg
orderedTable = Plotgenomecnv$orderedTable
Regions = Plotgenomecnv$orderedTableMedian2
chr_total2 = Plotgenomecnv$chr_total2

# Expression
DataPrep <- function(Organism, CancerType){
  
  setwd(paste(input_dir,cancer_type, sep = ""))
  
  datt <- read.delim(paste(cancer_type,"_t_geneENS_counts.txt", sep = ""), row.names = 1)
  
  # data normaliseren
  library("limma")
  library("edgeR")
  
  count_table<-datt
  count_table<-as.data.frame(count_table)
  tn<-1:ncol(datt)
  conditions<-paste("tumor",tn)
  colnames(count_table)<-conditions
  count_table<-na.omit(count_table)
  
  targets <- factor(rep("Tumor",length(datt)))
  #design <- model.matrix(~targets)
  dge <- DGEList(counts=count_table,group=targets)
  dge <- calcNormFactors(dge,method='RLE')
  
  counts_norm<-as.data.frame(cpm(dge, normalized.lib.sizes = TRUE))
  
  if(Organism == "Human"){
    centromere_pos=c(125,93.3,91,50.4,48.4,61,59.9,45.6,49,40.2,53.7,35.8,17.9,17.6,19,36.6,24,17.2,26.5,27.5,13.2,14.7,60.6,12.5)
    centromere_pos=centromere_pos*1000000
    chr_size=c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,155270560,59373566)
  }
  
  chr_total=0
  for (i in 1:(length(chr_size)-1)){
    chr_total=c(chr_total,sum(chr_size[1:i]))  
  }
  
  # chromosome 23 and 24 not in plot
  genome_size=sum(chr_size[1:22])
  
  # lijst van alle genen
  setwd(input_dir)
  allgenes <- read.table("allgenes.txt",header=TRUE,sep="\t")
  allgenes <- allgenes[order(allgenes[,1],allgenes[,2]),]
  colnames(allgenes) <- c('Chromosome', 'Start', 'Stop', 'GeneId')
  
  # vergelijk genen in allgenes met datt
  library(compare)
  dattgenes <- setDT(datt, keep.rownames = TRUE)[]
  comparison <- table(dattgenes[ which( dattgenes$rn %ni% allgenes$GeneId) , "rn"])
  
  # allgenes begint van 0, in deze file telt de positie op, dus chromosoom grote bij optellen
  for(x in 1:22){
    start = allgenes$Start[allgenes$Chromosome == x]
    start = as.numeric(start)
    allgenes$Start[allgenes$Chromosome == x] = start + chr_total[x]
    stop = allgenes$Stop[allgenes$Chromosome == x]
    stop = as.numeric(stop)
    allgenes$Stop[allgenes$Chromosome == x] = stop + chr_total[x]
  }
  setwd(output_dir)
  return(list(datt=datt, counts_norm=counts_norm, chr_total=chr_total, genome_size=genome_size, allgenes=allgenes, centromere_pos=centromere_pos, chr_size=chr_size))
}

# plaats van de code in de data frame zoeken
setwd(output_dir)
code <- basename(output_dir)
code <- gsub("-", ".", code)

Datapreperation <- DataPrep(Organism = "Human",CancerType = cancer_type)

datt = Datapreperation$datt
counts_norm = Datapreperation$counts_norm
allgenes = Datapreperation$allgenes

DataAnalysis <- function(Regions, chr_total, genome_size, allgenes, datt, counts_norm, centromere_pos, chr_size){
  
  data <- Regions
  
  results <- data.frame(Chromosome=numeric(), Start=numeric(), Stop=numeric(), GeneId=character(), stringsAsFactors=FALSE)
  #for(chr in unique(data$Chromosome)){
  for(chr in 1:nrow(data)){
    #chr = 12
    chrom = data$Chromosome[chr]
    
    # zoek alle genen in allgenes die in de regio van data liggen per chr
    dt <- (allgenes[allgenes$Chromosome==chrom, ])
    dt$GeneId <- as.character((dt$GeneId))
    
    for(chrr in 1:nrow(dt)){
      #chrr = 300
      if(dt$Start[chrr] >= data$Start[chr] & dt$Start[chrr] <= data$Stop[chr] & dt$Stop[chrr] >= data$Start[chr] & dt$Stop[chrr] <= data$Stop[chr]){
        newrow = c(Chromosome=data$Chromosome[chr], data$Start[chr], data$Stop[chr], dt$GeneId[chrr])
        results[nrow(results)+1, ] <- newrow
      }
      #chrr = chrr+1
    }
  }
  
  # TUMOR
  
  code2 <- substr(code, 1,16)
  k <- grep(code2, colnames(datt))
  listt <- paste("tumor",k)
  col.numb <- which(colnames(counts_norm) %in% listt)
  
  if(length(k) == 1) {
    exprt <- as.data.frame(counts_norm[,col.numb], row.names(counts_norm))
  } else{
    exprt <- counts_norm[,col.numb]
  }
  
  # eerst overal +1 doen en dan de log nemen
  exprt <- exprt + 1
  exprt[,1:ncol(exprt)] <- log(exprt[1:ncol(exprt)], 2)
  
  if(ncol(exprt) > 1){
    exprt$Mean <- rowMeans(exprt[,1:ncol(exprt)],na.rm = TRUE)
    exprt$Mean <- round(exprt$Mean, digits = 6)
    `%ni%` <- Negate(`%in%`)
    exprt <- subset(exprt,select = names(exprt) %ni% listt)
  }
  
  library("data.table")
  exprt <- setDT(exprt, keep.rownames = TRUE)[]
  
  colnames(exprt) <- c('GeneId', 'Expression')
  
  # searching for GeneId and merging the expression to data
  genes <- exprt[exprt$GeneId %in% results$GeneId,]
  genest <- merge(results, genes, all = TRUE, sort = FALSE)
  #genest <- genest[order(genest[,2]),]
  genest <- na.omit(genest)
  
  #genest <- genest[c("Chromosome", "Position","GeneId", "Expression")]
  genest <- genest[c("Chromosome", "Start", "Stop","GeneId", "Expression")]
  
  # get allgenes of this tumor code
  allgenest <- exprt[exprt$GeneId %in% allgenes$GeneId,]
  allgenest <- merge(allgenes, allgenest, sort = FALSE)
  allgenest <- na.omit(allgenest)
  allgenest <- allgenest[c("Chromosome", "Start", "Stop","GeneId", "Expression")]
  
  # Table with rollmean over expression from allgenest
  allgenestroll <- rollmean(allgenest$Expression,window, fill = NA)
  allgenestroll = data.frame(allgenest$Chromosome,allgenest$GeneId, allgenestroll)
  colnames(allgenestroll) <- c('Chromosome','GeneId', 'Expression')
  allgenestroll <- na.omit(allgenestroll)
  
  # # regions of interest from table with rollmean
  # rollgenes <- allgenestroll[allgenestroll$GeneId %in% results$GeneId,]
  # rollgenes <- merge(results, rollgenes, all = TRUE, sort = FALSE)
  # #genest <- genest[order(genest[,2]),]
  # rollgenes <- na.omit(rollgenes)
  # rollgenes <- rollgenes[c("Chromosome", "Start", "Stop","GeneId", "Expression")]
  
  return(list(data=data, tumor=genest, tumormean=rollgenes))
}

Data_analysis <- DataAnalysis(Regions = Regions, chr_total=Datapreperation$chr_total, genome_size=Datapreperation$genome_size, allgenes=Datapreperation$allgenes, datt=Datapreperation$datt, counts_norm=Datapreperation$counts_nor, centromere_pos=Datapreperation$centromere_pos, chr_size=Datapreperation$chr_size)

# ratio tumor min door normaal
tumor = Data_analysis$tumor
tumormean = Data_analysis$tumormean
data = Data_analysis$data
genome_size = Datapreperation$genome_size
chr_total = Datapreperation$chr_total
chr_size = Datapreperation$chr_size
centromere_pos = Datapreperation$centromere_pos

# density plot
d <- density(allgenestroll$Expression) # returns the density data 
plot(d) # plots the results
allgenesmean = mean(allgenestroll$Expression)

plottable <- data.frame(Chromosome=tumor$Chromosome, Position=tumor$Start, GeneId=tumor$GeneId, ExpressionT=tumor$Expression,ExpressionM=allgenesmean, ExpressionMean=tumor$Expression-allgenesmean)
plottable$Position = as.character(plottable$Position)

# als de expressie van een gen zowel in normaal als in tumor 0 is dan verwijderen
plottable<-plottable[!(plottable$ExpressionT==0),]

# dataframe voor nullijn tussen de data met expressie met rollmedian
data2 <- data.frame(Chromosome=numeric(), Start=numeric(), Stop=numeric())
row1 <- c(Chromosome=data$Chromosome[1], Start=0, Stop=data$Start[1])
data2[nrow(data2)+1, ] <- row1
for(i in 1:nrow(data)){
  newrow = c(data$Chromosome[i], data$Stop[i], data$Start[i+1])
  data2 <- rbind(data2,newrow)
}
data2$Stop[is.na(data2$Stop)] = chr_total2[data$Chromosome[i]]
if(tail(data2$Chromosome, n=1) != 22){
  newrow = c(Chromosome=22, Start = tail(data2$Stop, n=1), Stop = chr_total2[22])
  data2 <- rbind(data2,newrow)
}

data2$Expression = 0

# ROLLMEAN
window = 251

ymin = abs(min(rollmean(plottable$ExpressionMean,window)))
ymin = ceiling(ymin)
ymin = -ymin
ymax = 1
ymax2 = max(rollmean(plottable$ExpressionMean,window))
ymax2 = ceiling(ymax2)
if(ymax2 > ymax){ymax = ymax2}else{ymax=ymax}


# PLOT 1 LINE WITH NO MEAN IN SMALL REGIONS
oneline <- data.frame(X=numeric(),Y=numeric())
for(i in 1:nrow(data2)){
  arow <- c(X=data2$Start[i],Y=data2$Expression[i])
  somerow <- c(X=data2$Stop[i],Y=data2$Expression[i])
  oneline[nrow(oneline)+1, ] <- arow
  oneline[nrow(oneline)+1, ] <- somerow
}


# add rollmedian of regions with expression
for(x in data$Start){
  rollmediandata <- (tumor[tumor$Start==x, ])
  rollmediandata2 <- allgenes[(allgenes$GeneId %in% rollmediandata$GeneId),]
  rollmediandata2 <- merge(rollmediandata, rollmediandata2, by = "GeneId", sort = FALSE)
  rollmediandata2 <- merge(rollmediandata2, plottable, by = "GeneId", sort = FALSE)
  rollmediandata <- data.frame(Chromosome = rollmediandata2$Chromosome.x, Start = rollmediandata2$Start.y, GeneId = rollmediandata2$GeneId, ExpressionMean = rollmediandata2$ExpressionMean)
  rollmediandata$Start <- as.numeric(rollmediandata$Start)
  if(nrow(rollmediandata) >= window){
    dtframe <- data.frame(X=rollmedian(rollmediandata$Start,window),Y=rollmean(rollmediandata$ExpressionRatio,window))
    oneline <- as.data.frame(rbind(oneline, dtframe))
  }
}

oneline <- oneline[order(oneline[,1]),]



#png(filename=paste("MeanExpression2-",substr(code, 1,28), ".png", sep=""), width = 1980, height = 1080, units = "px")

par(mar=c(4, 7, 4, 2))
par(mfrow=c(1,1))
par(las=1)
plot(1, type="n", xlim=c(0, genome_size), ylim=c(ymin,ymax),  xlab="", ylab="Expression\n", xaxt = "n", cex.axis=1.5,cex.lab=1.5)
par(new=TRUE)
plot(oneline$X, oneline$Y, type="l", xlim=c(0, genome_size), ylim=c(ymin,ymax), ylab="", xaxt = "n", yaxt = "n",col="554", xlab="",cex.axis=1.5,cex.lab=1.5, lwd=2.5)

for(i in 1 :22){
  if (i>1){abline(v=chr_total[i],col="gray48")}
  abline(v=chr_total[i]+centromere_pos[i],col="gray55",lty=4)
  #mtext(chr_total[i]+centromere_pos[i], side=1,at=chr_total[i]+centromere_pos[i], cex=0.5)
  mtext(as.character(i),side=3, at=chr_total[i]+chr_size[i]/2,cex=1)
}

par(las=1)
plot(rollmean(allgenest$Start,window),rollmean(allgenest$Expression,window),col="554",pch=15,cex=0.4,ylim=c(0,5), ylab="Expression\n",typ="l",xlab="Chromosomal position",xaxt = "n",xlim=c(1,genome_size), cex.axis=1.5,cex.lab=1.5, lwd=2.5)
for(i in 1 :22){
  if (i>1){abline(v=chr_total[i],col="gray48")}
  abline(v=chr_total[i]+centromere_pos[i],col="gray55",lty=4)
  #mtext(chr_total[i]+centromere_pos[i], side=1,at=chr_total[i]+centromere_pos[i], cex=0.5)
  mtext(as.character(i),side=3, at=chr_total[i]+chr_size[i]/2,cex=1)
}

#dev.off()

setwd(output_dir)

#png(filename=paste("Comparisonplot_new2-",substr(code, 1,28), ".png", sep=""), width = 1980, height = 1080, units = "px")

par(mar=c(5, 9, 4, 2))
par(mfrow=c(3,1))
par(las=1)
plot(rollmedian(orderedTable$position,window),rollmedian(orderedTable$MajorMinor,window),col="554",pch=15,cex=0.4,ylim=c(1,Ylim), ylab="Allelic Ratio\n",typ="l",xlab="",xaxt = "n",xlim=c(1,genome_size), cex.axis=2,cex.lab=3,font.axis=2, lwd=2.5)
abline(1.75,0,col="black")
# Draw chromosome guide lines + chromosoom cijfers
for(i in 1 :22){
  if (i>1){abline(v=chr_total[i],col="gray48")}
  abline(v=chr_total[i]+centromere_pos[i],col="gray55",lty=4)
  #mtext(chr_total[i]+centromere_pos[i], side=1,at=chr_total[i]+centromere_pos[i], cex=0.5)
  mtext(as.character(i),side=3, at=chr_total[i]+chr_size[i]/2,cex=1.5)
}
par(las=1)
plot(cnv_lg$Position,cnv_lg$Segment_Mean,col="28",pch=15,cex=0.4,type="l",xlim=c(1,genome_size),ylim=c(-1.5,1.5), ylab="Copy number\n",xlab="",xaxt = "n", cex.axis=2,cex.lab=3, font.axis=2,lwd=2.5)
abline(0,0,col="black")
# Draw chromosome guide lines + chr cijfers
for(i in 1 :22){
  if (i>1){abline(v=chr_total[i],col="gray48")}
  abline(v=chr_total[i]+centromere_pos[i],col="gray55",lty=4)
  mtext(as.character(i),side=3, at=chr_total[i]+chr_size[i]/2,cex=1.5)
}
# par(las=1)
# plot(plottable$Position, plottable$ExpressionRatio,pch=15,cex=0.4, type = "p",col="554", xlab="Chromosomal position", ylab="Expression\n", xaxt = "n",xlim=c(1,genome_size), ylim=c(ymin,ymax), cex.axis=1.5,cex.lab=1.5, lwd=2.5)
# abline(0,0,col="black")
# # Draw chromosome guide lines + Chromosome cijfers
# for(y in 1 :22){
#   if (y>1){abline(v=chr_total[y],col="gray48")}
#   abline(v=chr_total[y]+centromere_pos[y],col="gray55",lty=4)
#   #mtext(chr_total[y]+centromere_pos[y], side=1,at=chr_total[y]+centromere_pos[y], cex=0.5)
#   mtext(as.character(y),side=3, at=chr_total[y]+chr_size[y]/2,cex=1)
# }

par(las=1)
plot(1, type="n", xlim=c(0, genome_size), ylim=c(ymin,ymax),  xlab="Chromosomal position", ylab="Expression\n", xaxt = "n", cex.axis=2,cex.lab=3,font.axis=2)
par(new=TRUE)
plot(oneline$X, oneline$Y, type="l", xlim=c(0, genome_size), ylim=c(ymin,ymax), ylab="", xaxt = "n",col="554", xlab="",cex.axis=2,cex.lab=3, font.axis=2,lwd=2.5)

for(i in 1 :22){
  if (i>1){abline(v=chr_total[i],col="gray48")}
  abline(v=chr_total[i]+centromere_pos[i],col="gray55",lty=4)
  #mtext(chr_total[i]+centromere_pos[i], side=1,at=chr_total[i]+centromere_pos[i], cex=0.5)
  mtext(as.character(i),side=3, at=chr_total[i]+chr_size[i]/2,cex=1.5)
}

#dev.off()




