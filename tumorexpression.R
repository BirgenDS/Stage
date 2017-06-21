# Script voor het maken van plots met enkel tumor data

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
  
  # chromosome 23 and 24 niet in plot
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
  
  # Verwijderd rijen van chr 23 en 24 uit de tabellen omdat deze niet hoeven geplot te worden
  orderedTable3 <- orderedTable[orderedTable$chr != 23, ]
  orderedTable <- orderedTable3[orderedTable3$chr != 24, ]
  cnv_lg2 <- cnv_lg[cnv_lg$Chromosome != 23, ] 
  cnv_lg <- cnv_lg2[cnv_lg2$Chromosome != 24, ]
  
  # rollmedian over orderedtable$MajorMinor zoals in plot onderaan dit script
  # orderedTableMedian is een tabel met de rijen waar deze rollmedian > 1.75
  orderedTableMedian <- rollmedian(orderedTable$MajorMinor,window, fill = NA)
  orderedTable2 = data.frame(orderedTable$chr,orderedTable$position, orderedTableMedian)
  colnames(orderedTable2) <- c('Chromosome','Position', 'Rollmedian')
  orderedTable2 <- orderedTable2[orderedTable2$Rollmedian >= 1.75, ] # aangepast van 2 naar 1.75 om meer pieken over te houden
  orderedTable2 <- data.frame(orderedTable2$Chromosome,orderedTable2$Position, orderedTable2$Rollmedian)
  colnames(orderedTable2) <- c('Chromosome','Position', 'Rollmedian')
  orderedTableMedian <- na.omit(orderedTable2)
  #write.csv(orderedTableMedian, file = 'ordereTableMedian.csv', row.names=FALSE)
  
  # chr_total lijst die begint bij de grote van chr 1
  chr_total2 <- chr_total[-1];
  
  # loop om de allelic imballance > 1,75 te selecteren per chromosoom
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
  # op de laatste rij de NA waarden invullen
  orderedTableMedian2$Stop[is.na(orderedTableMedian2$Stop)] = orderedTableMedian$Position[i-1]
  orderedTableMedian2$Rollmedian[is.na(orderedTableMedian2$Rollmedian)] = mean(rollmedianlist)
  # Enkel rijen overhouden waar het verschil tussen start en stop > 1000
  orderedTableMedian2 <- orderedTableMedian2[(orderedTableMedian2$Stop - orderedTableMedian2$Start) > 1000, ]
  # Wegschrijven van deze tabel naar een csv bestand in de output_dir
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
  
  # gaat naar de map van het kankertype
  setwd(paste(input_dir,cancer_type, sep = ""))
  
  # inladen data
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
  
  # chromosome 23 and 24 niet in plot
  genome_size=sum(chr_size[1:22])
  
  # lijst van alle genen inladen die gehaald is van ensembl
  setwd(input_dir)
  allgenes <- read.table("allgenes.txt",header=TRUE,sep="\t")
  allgenes <- allgenes[order(allgenes[,1],allgenes[,2]),]
  colnames(allgenes) <- c('Chromosome', 'Start', 'Stop', 'GeneId')
  
  # vergelijken van deze genen in allgenes met de genen in datt om te controleren of degene die daar meer in zitten weldegelijk de genen op chr 23 en 24 zijn
  comparison <- table(dattgenes[ which( dattgenes$rn %ni% allgenes$GeneId) , "rn"])
  
  # allgenes begint per chromosome opnieuw van 0, dus chromosoom grote bij optellen
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

setwd(output_dir)
# code is het laatste deel van de output_dir
code <- basename(output_dir)
# streepjes in code vervangen door punten omdat in de titels van datt de code met punten gescheiden zijn
code <- gsub("-", ".", code)

Datapreperation <- DataPrep(Organism = "Human",CancerType = cancer_type)

datt = Datapreperation$datt
counts_norm = Datapreperation$counts_norm
allgenes = Datapreperation$allgenes

DataAnalysis <- function(Regions, chr_total, genome_size, allgenes, datt, counts_norm, centromere_pos, chr_size){
  
  data <- Regions
  
  # dataframe (results) maken met alle genen die in de regio's van allelic imballences liggen
  # start en stop is bij alle genen die in deze regio liggen de start en de stop van de regio
  results <- data.frame(Chromosome=numeric(), Start=numeric(), Stop=numeric(), GeneId=character(), stringsAsFactors=FALSE)
  #for(chr in unique(data$Chromosome)){
  for(chr in 1:nrow(data)){
    #chr = 12
    chrom = data$Chromosome[chr]
    
    # dataframe met alle genen uit allgenes die in de regio van data liggen per chr
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
  
  # code waarmee we zoeken in datt om ons staal er uit te halen
  code2 <- substr(code, 1,16)
  k <- grep(code2, colnames(datt))
  # de nummers van de kolommen met deze code in een lijst steken met het woord tumor ervoor om hier mee te zoeken in de genormaliseerde data, counts_norm
  listt <- paste("tumor",k)
  col.numb <- which(colnames(counts_norm) %in% listt)
  
  # een aparte data frame maken met de kolommen uit counts_norm
  # if voor als het maar 1 kolom is, else voor als het meerdere kolommen zijn
  if(length(k) == 1) {
    exprt <- as.data.frame(counts_norm[,col.numb], row.names(counts_norm))
  } else{
    exprt <- counts_norm[,col.numb]
  }
  
  # alle getallen +1 doen en dan log2 ervan nemen nemen
  exprt <- exprt + 1
  exprt[,1:ncol(exprt)] <- log(exprt[1:ncol(exprt)], 2)
  
  # als er meerdere kolommen matchen met de code en in exprt zijn gestoken, per rij het gemiddelde nemen
  # getallen afronden op 6 cijfers
  # dataframe maken met enkel de kolom met het gemiddelde per rij
  if(ncol(exprt) > 1){
    exprt$Mean <- rowMeans(exprt[,1:ncol(exprt)],na.rm = TRUE)
    exprt$Mean <- round(exprt$Mean, digits = 6)
    `%ni%` <- Negate(`%in%`)
    exprt <- subset(exprt,select = names(exprt) %ni% listt)
  }
  
  # van de rijnamen een kolom maken
  library("data.table")
  exprt <- setDT(exprt, keep.rownames = TRUE)[]
  
  colnames(exprt) <- c('GeneId', 'Expression')
  
  # dataframe matchen met results op basis van gene-id om de expressiewaarden van die genen er bij te voegen
  genes <- exprt[exprt$GeneId %in% results$GeneId,]
  genest <- merge(results, genes, all = TRUE, sort = FALSE)
  #genest <- genest[order(genest[,2]),]
  genest <- na.omit(genest)
  
  genest <- genest[c("Chromosome", "Start", "Stop","GeneId", "Expression")]
  
  # dataframe maken met alle genen van dit staal, door de expressie aan de genen te matchen
  allgenest <- exprt[exprt$GeneId %in% allgenes$GeneId,]
  # chromosoom, start en stop toevoegen
  allgenest <- merge(allgenes, allgenest, sort = FALSE)
  allgenest <- allgenest[c("Chromosome", "Start", "Stop","GeneId", "Expression")]
  
  # dataframe met rollmean over de expressie
  allgenestroll <- rollmean(allgenest$Expression,window, fill = NA)
  allgenestroll = data.frame(allgenest$Chromosome,allgenest$GeneId, allgenestroll)
  colnames(allgenestroll) <- c('Chromosome','GeneId', 'Expression')
  allgenestroll <- na.omit(allgenestroll)
  
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

# density plot van expressie met rollmean
d <- density(allgenestroll$Expression)
plot(d)
allgenesmean = mean(allgenestroll$Expression)

# dataframe met tumor data, allgenesmean en het verschil van de 2 om te kijken of de rollmean in de regios links of rechts van het gemiddelde ligt
plottable <- data.frame(Chromosome=tumor$Chromosome, Position=tumor$Start, GeneId=tumor$GeneId, ExpressionT=tumor$Expression,ExpressionM=allgenesmean, ExpressionMean=tumor$Expression-allgenesmean)
plottable$Position = as.character(plottable$Position)

# als de expressie van een gen 0 is dan verwijderen
plottable<-plottable[!(plottable$ExpressionT==0),]

# dataframe voor nullijn tussen de regios met expressie
# chromosoom 1 begint bij 1 tot de eerste start en op rij 2 van de eerste stop tot de tweede start en zo schuift alles 1 op
data2 <- data.frame(Chromosome=numeric(), Start=numeric(), Stop=numeric())
row1 <- c(Chromosome=data$Chromosome[1], Start=0, Stop=data$Start[1])
data2[nrow(data2)+1, ] <- row1
for(i in 1:nrow(data)){
  newrow = c(data$Chromosome[i], data$Stop[i], data$Start[i+1])
  data2 <- rbind(data2,newrow)
}
data2$Stop[is.na(data2$Stop)] = chr_total2[data$Chromosome[i]]
# eindpunt van de nullijn invoegen als er geen regio is met expressie van chromosoom 22
if(tail(data2$Chromosome, n=1) != 22){
  newrow = c(Chromosome=22, Start = tail(data2$Stop, n=1), Stop = chr_total2[22])
  data2 <- rbind(data2,newrow)
}

data2$Expression = 0


window = 251

# variabele limiten van de y-as afhankelijk van het grootste en kleinste getal in de expressie kolom
ymin = abs(min(rollmean(plottable$ExpressionMean,window)))
ymin = ceiling(ymin)
ymin = -ymin
ymax = 1
ymax2 = max(rollmean(plottable$ExpressionMean,window))
ymax2 = ceiling(ymax2)
if(ymax2 > ymax){ymax = ymax2}else{ymax=ymax}


# dataframe maken met de gegevens voor de plot
# nullijn aan de hand van data2
oneline <- data.frame(X=numeric(),Y=numeric())
for(i in 1:nrow(data2)){
  arow <- c(X=data2$Start[i],Y=data2$Expression[i])
  somerow <- c(X=data2$Stop[i],Y=data2$Expression[i])
  oneline[nrow(oneline)+1, ] <- arow
  oneline[nrow(oneline)+1, ] <- somerow
}


# rollmedian over de regio's met expressie
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


# apart plotten om te kijken of het klopt
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

# plot met drie lijnen onder elkaar, lijn 1 allelic ratio, lijn 2 copy number, lijn 3 expressie
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




