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
#source("eSNPKaryotyping/R/MajorMinorCalc.R")
source("eSNPKaryotyping/R/Sort_major_minor.R")

# wrapperscript
output_dir <- "/Users/Birgen/Documents/Howest/Stage/CMGG/RNA-seq_CNV/KIRC/TCGA-A3-3306-01A-01R-0864-07_eSNPout/"
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
  return(list(cnv_lg = cnv_lg, orderedTable = orderedTable))
}



# Expressionplot
Plotgenomecnv <- PlotGenome(Output_Directory = output_dir,orderedTable = table2,Window = 151,Ylim = 3,Organism = "Human", CancerType = cancer_type, CNV_Directory = cnv_dir)
cnv_lg = Plotgenomecnv$cnv_lg
orderedTable = Plotgenomecnv$orderedTable

#functie 1
DataPrep <- function(Organism, CancerType){
  
  setwd(paste(input_dir,cancer_type, sep = ""))
  
  datn <- read.delim(paste(cancer_type,"_n_geneENS_counts.txt", sep = ""), row.names = 1)
  datt <- read.delim(paste(cancer_type,"_t_geneENS_counts.txt", sep = ""), row.names = 1)
  
  # aantal tumor stalen in de normale data, dus alle stalen met 01B in het midden mogen verwijderd 
  drop.cols <- grep("01B", colnames(datn))
  if(length(drop.cols) != 0){
    datn <- datn[, -grep("01B", colnames(datn))]
  } else{datn <- datn}
  
  # data normaliseren
  library("limma")
  library("edgeR")
  
  count_table<-cbind(datn,datt)
  count_table<-as.data.frame(count_table)
  nn<-1:ncol(datn)
  tn<-1:ncol(datt)
  conditions<-c(paste("normal",nn),paste("tumor",tn))
  colnames(count_table)<-conditions
  count_table<-na.omit(count_table)
  
  targets <- factor(c(rep("Normal",length(datn)),rep("Tumor",length(datt))))
  design <- model.matrix(~targets)
  dge <- DGEList(counts=count_table,group=targets)
  #dge$counts <- dge$counts[rowSums(dge$counts > ncol(datn)) >= 10,]
  dge <- calcNormFactors(dge,method='RLE')
  
  counts_norm<-as.data.frame(cpm(dge, normalized.lib.sizes = TRUE))
  
  # referentiestaal toevoegen aan datn voor als er geen matching normal is
  counts_norm$mean <- rowMeans(counts_norm[,])
  
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
  return(list(datn=datn, datt=datt, counts_norm=counts_norm, chr_total=chr_total, genome_size=genome_size, allgenes=allgenes, centromere_pos=centromere_pos, chr_size=chr_size))
}

#einde functie 1

# plaats van de code in de data frame zoeken
setwd(output_dir)
filename <- (paste(output_dir, "Regions.csv",sep = ""))
code <- basename(output_dir)
code <- gsub("-", ".", code)

# file_list = list.files(pattern="*.csv")
# data_list <- vector("list", "length" = length(file_list))  


# functie 2
Datapreperation <- DataPrep(Organism = "Human",CancerType = cancer_type)

DataAnalysis <- function(chr_total, genome_size, allgenes, datn, datt, counts_norm, centromere_pos, chr_size){
  
  data <- read.csv(filename, header = TRUE, sep = "\t")
  
  
  # genen uit allgenes halen die in de regio's zitten uit bovenstaande dataframe
  results <- data.frame(Chromosome=numeric(), Start=numeric(), Stop=numeric(), GeneId=character(), stringsAsFactors=FALSE)
  for(chr in 1:nrow(data)){
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

  # NORMAL
  # dataframe from counts_norm with only this code
  
  code1 <- substr(code, 1,12)
  j <- grep(code1, colnames(datn))
  listn <- paste("normal",j)
  col.num <- which(colnames(counts_norm) %in% listn)
  
  if(length(j) == 0) {
    exprn <- as.data.frame(counts_norm[,ncol(counts_norm)], row.names(counts_norm))
  } else if(length(j) == 1){
    exprn <- as.data.frame(counts_norm[,col.num], row.names(counts_norm))
  } else{
    exprn <- counts_norm[,col.num]
  }
  
  # eerst overal +1 doen en dan de log nemen
  exprn <- exprn + 1
  exprn[,1:ncol(exprn)] <- log(exprn[1:ncol(exprn)], 2)
  
  # als meer dan 1 match in data, mean nemen van expressiewaarden
  if(ncol(exprn) > 1){
    exprn$Mean <- rowMeans(exprn[,1:ncol(exprn)],na.rm = TRUE)
    exprn$Mean <- round(exprn$Mean, digits = 6)
    `%ni%` <- Negate(`%in%`)
    exprn <- subset(exprn,select = names(exprn) %ni% listn)
  }
  
  library("data.table")
  exprn <- setDT(exprn, keep.rownames = TRUE)[]
  
  colnames(exprn) <- c('GeneId', 'Expression')
  
  # searching for GeneId and merging the expression to genes in results
  genes <- exprn[exprn$GeneId %in% results$GeneId,]
  genesn <- merge(results, genes, all = TRUE, sort = FALSE)
  #genesn <- genesn[order(genesn[,2]),]
  genesn <- na.omit(genesn)
  
  genesn <- genesn[c("Chromosome", "Start", "Stop","GeneId", "Expression")]
  
  # pasting start en stop together in region
  start = genesn[,c('Chromosome', 'Start',  'GeneId','Expression')]
  stop = genesn[,c('Chromosome', 'Stop', 'GeneId', 'Expression')]
  plottable1<-transform(genesn, Region=paste(Start, Stop, sep="-"))
  plottable1<- plottable1[,c('Chromosome', 'Region', 'GeneId','Expression')]
  plottable1$Region <- as.character(plottable1$Region)
  
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
  
  genest <- genest[c("Chromosome", "Start", "Stop","GeneId", "Expression")]
  
  start2 = genest[,c('Chromosome', 'Start',  'GeneId','Expression')]
  stop2 = genest[,c('Chromosome', 'Stop', 'GeneId', 'Expression')]
  plottable2<-transform(genest, Region=paste(Start, Stop, sep="-"))
  plottable2<- plottable2[,c('Chromosome', 'Region', 'GeneId','Expression')]
  plottable2$Region <- as.character(plottable2$Region)
  
  # NORMAL

  #tabel maken met gemiddelde expressiewaarden per regio
  plottablen <- data.frame(Chromosome=numeric(), Region=numeric(), GeneId=character(), Expression=numeric(), stringsAsFactors=FALSE)
  row1 <- c(Chromosome=plottable1$Chromosome[1], Region=plottable1$Region[1], GeneId=plottable1$GeneId[1] ,Expression=c(NA))
  plottablen[nrow(plottablen)+1, ] <- row1
  Expressionlist <- c(plottable1$Expression[1])
  #a = 2
  for (a in 2:nrow(plottable1)){
    if (plottable1$Region[a-1] != plottable1$Region[a]) {
      plottablen$Expression[is.na(plottablen$Expression)] = mean(Expressionlist)
      newrow = c(plottable1$Chromosome[a], as.character(plottable1$Region[a]),as.character(plottable1$GeneId[a]),plottablen$Expression[NA])
      plottablen <- rbind(plottablen,newrow)
      Expressionlist <- c()
    }
    Expressionlist <- append(Expressionlist,plottable1$Expression[a])
    #a = a+1
  }
  # na de laatste loop ook de mean van expressionlist toevoegen aan de dataframe
  plottablen$Expression[is.na(plottablen$Expression)] = mean(Expressionlist)

  plottablen$Chromosome <- as.numeric(plottablen$Chromosome)
  plottablen$Region <- as.character(plottablen$Region)
  plottablen$GeneId <- as.character(plottablen$GeneI)
  plottablen$Expression <- as.numeric(plottablen$Expression)

  # regio splitsen in start en stop en onder elkaar plakken
  plottablen <- separate(data = plottablen, col = Region, into = c("Start", "Stop"), sep = "\\-")

  start = plottablen[,c('Chromosome', 'Start', 'Expression')]
  colnames(start) <- c('Chromosome', 'Position', 'Expression')
  stop = plottablen[,c('Chromosome', 'Stop', 'Expression')]
  colnames(stop) <- c('Chromosome', 'Position', 'Expression')

  plottablen<-as.data.frame(rbind(start, stop)) # plakt stop kolom onder start
  plottablen$Position <- as.numeric(plottablen$Position)
  plottablen <- plottablen[order(plottablen[,2],plottablen[,1]),]

  # alle chormosomen in tabel zetten en ontbrekende stukken van aanwezige chromosomen
  normal <- data.frame(Chromosome=numeric(), Position=numeric(), Expression=numeric(), stringsAsFactors=FALSE)
  for(x in unique(plottablen$Chromosome)){
    df <- (plottablen[plottablen$Chromosome==x, ])

    for(xx in 2:nrow(df)){
      #xx = 2
      if(nrow(df) == 2){
        position = df$Position[xx-1]
        position = as.numeric(position)
        position1 = position - chr_total[x]
        if(position1 != 0){
          newererrow= c(x,chr_total[x],0)
          newestrow= c(x,df$Position[xx-1],0)
          normal[nrow(normal)+1, ] <- newererrow
          normal[nrow(normal)+1, ] <- newestrow
          startrow= c(x,df$Position[xx],0)
          stoprow= c(x,chr_total[x+1],0)
          normal[nrow(normal)+1, ] <- startrow
          normal[nrow(normal)+1, ] <- stoprow
        } else if(position1 == 0) {
          newererrow= c(x,df$Position[xx],0)
          newestrow= c(x,chr_total[x],0)
          normal[nrow(normal)+1, ] <- newererrow
          normal[nrow(normal)+1, ] <- newestrow
        }
      }
      if(nrow(df) > 2){
        if(df$Expression[xx-1] != df$Expression[xx]){
          position = df$Position[xx-1]
          position = as.numeric(position)
          newererrow= c(x,position,0)
          newestrow= c(x,df$Position[xx],0)
          normal[nrow(normal)+1, ] <- newererrow
          normal[nrow(normal)+1, ] <- newestrow
        } else if(df$Expression[xx-1] == df$Expression[xx]){
          position = df$Position[xx-1]
          position = as.numeric(position)
          position1 = position - chr_total[x]
          if(position1 != 0 & position == df$Position[1]){
            newererrow= c(x,chr_total[x],0)
            newestrow= c(x,df$Position[xx-1],0)
            normal[nrow(normal)+1, ] <- newererrow
            normal[nrow(normal)+1, ] <- newestrow
          } else if(position1 != 0 & df$Position[xx] == tail(df$Position, n=1)) {
            newererrow= c(x,df$Position[xx],0)
            newestrow= c(x,chr_total[x+1],0)
            normal[nrow(normal)+1, ] <- newererrow
            normal[nrow(normal)+1, ] <- newestrow
          }
        }
        #xx = xx+1
      }
      #x = x+1
    }
  }

  for(x in 1:22){
    #x = 1
    if(!(x %in% plottablen$Chromosome)){
      newrow = c(x,chr_total[x],0)
      newerrow = c(x,chr_total[x+1],0)
      plottablen <- rbind(plottablen,newrow,newerrow)
    }
  }

  # rank aan toevoegen en sorteren
  plottablen$Position <- as.numeric(plottablen$Position)
  plottablen <- plottablen[order(plottablen[,2],plottablen[,1]),]
  normal$Rank <- c(2,0)
  plottablen$Rank <- 1

  normal2 <- rbind(plottablen, normal)
  normal2 <- normal2[order(normal2[,1],normal2[,2], normal2[,4]),]


  # TUMOR

  #tabel maken met gemiddelde expressiewaarden per regio
  plottablet <- data.frame(Chromosome=numeric(), Region=numeric(), GeneId=character(), Expression=numeric(), stringsAsFactors=FALSE)
  row1 <- c(Chromosome=plottable2$Chromosome[1], Region=plottable2$Region[1], GeneId=plottable2$GeneId[1] ,Expression=c(NA))
  plottablet[nrow(plottablet)+1, ] <- row1
  Expressionlist <- c(plottable2$Expression[1])
  #z = 2
  for (z in 2:nrow(plottable2)){
    if (plottable2$Region[z-1] != plottable2$Region[z]) {
      plottablet$Expression[is.na(plottablet$Expression)] = mean(Expressionlist)
      newrow = c(plottable2$Chromosome[z], as.character(plottable2$Region[z]),as.character(plottable2$GeneId[z]),plottablet$Expression[NA])
      plottablet <- rbind(plottablet,newrow)
      Expressionlist <- c()
    }
    Expressionlist <- append(Expressionlist,plottable2$Expression[z])
    #z = z+1
  }
  # na de laatste loop ook de mean van expressionlist toevoegen aan de dataframe
  plottablet$Expression[is.na(plottablet$Expression)] = mean(Expressionlist)

  plottablet$Chromosome <- as.numeric(plottablet$Chromosome)
  plottablet$Region <- as.character(plottablet$Region)
  plottablet$GeneId <- as.character(plottablet$GeneI)
  plottablet$Expression <- as.numeric(plottablet$Expression)

  plottablet <- separate(data = plottablet, col = Region, into = c("Start", "Stop"), sep = "\\-")

  start = plottablet[,c('Chromosome', 'Start', 'Expression')]
  colnames(start) <- c('Chromosome', 'Position', 'Expression')
  stop = plottablet[,c('Chromosome', 'Stop', 'Expression')]
  colnames(stop) <- c('Chromosome', 'Position', 'Expression')

  plottablet<-as.data.frame(rbind(start, stop)) # plakt stop kolom onder start
  plottablet$Position <- as.numeric(plottablet$Position)
  plottablet <- plottablet[order(plottablet[,2],plottablet[,1]),]

  # alle chormosomen in tabel zetten en ontbrekende stukken van bestaande chromosomen
  tumor <- data.frame(Chromosome=numeric(), Position=numeric(), Expression=numeric(),stringsAsFactors=FALSE)
  for(x in unique(plottablet$Chromosome)){
    df <- (plottablet[plottablet$Chromosome==x, ])

    for(xx in 2:nrow(df)){
      #xx = 2
      if(nrow(df) == 2){
        position = df$Position[xx-1]
        position = as.numeric(position)
        position1 = position - chr_total[x]
        if(position1 != 0){
          newererrow= c(x,chr_total[x],0)
          newestrow= c(x,df$Position[xx-1],0)
          tumor[nrow(tumor)+1, ] <- newererrow
          tumor[nrow(tumor)+1, ] <- newestrow
          startrow= c(x,df$Position[xx],0)
          stoprow= c(x,chr_total[x+1],0)
          tumor[nrow(tumor)+1, ] <- startrow
          tumor[nrow(tumor)+1, ] <- stoprow
        } else if(position1 == 0) {
          newererrow= c(x,df$Position[xx],0)
          newestrow= c(x,chr_total[x],0)
          tumor[nrow(tumor)+1, ] <- newererrow
          tumor[nrow(tumor)+1, ] <- newestrow
        }
      }
      if(nrow(df) > 2){
        if(df$Expression[xx-1] != df$Expression[xx]){
          position = df$Position[xx-1]
          position = as.numeric(position)
          newererrow= c(x,position,0)
          newestrow= c(x,df$Position[xx],0)
          tumor[nrow(tumor)+1, ] <- newererrow
          tumor[nrow(tumor)+1, ] <- newestrow
        } else if(df$Expression[xx-1] == df$Expression[xx]){
          position = df$Position[xx-1]
          position = as.numeric(position)
          position1 = position - chr_total[x]
          if(position1 != 0 & position == df$Position[1]){
            newererrow= c(x,chr_total[x],0)
            newestrow= c(x,df$Position[xx-1],0)
            tumor[nrow(tumor)+1, ] <- newererrow
            tumor[nrow(tumor)+1, ] <- newestrow
          } else if(position1 != 0 & df$Position[xx] == tail(df$Position, n=1)) {
            newererrow= c(x,df$Position[xx],0)
            newestrow= c(x,chr_total[x+1],0)
            tumor[nrow(tumor)+1, ] <- newererrow
            tumor[nrow(tumor)+1, ] <- newestrow
          }
        }
        #xx = xx+1
      }
      #x = x+1
    }
  }

  for(x in 1:22){
    #x = 1
    if(!(x %in% plottablet$Chromosome)){
      newrow = c(x,chr_total[x],0)
      newerrow = c(x,chr_total[x+1],0)
      plottablet <- rbind(plottablet,newrow,newerrow)
    }
  }
  plottablet$Position <- as.numeric(plottablet$Position)
  plottablet <- plottablet[order(plottablet[,2],plottablet[,1]),]
  tumor$Rank <- c(2,0)
  plottablet$Rank <- 1

  tumor2 <- rbind(plottablet, tumor)
  tumor2 <- tumor2[order(tumor2[,1],tumor2[,2], tumor2[,4]),]
  
  return(list(normal2=normal2, tumor2=tumor2))
}
  
# einde functie 2

# PLOT
Data_analysis <- DataAnalysis(chr_total=Datapreperation$chr_total, genome_size=Datapreperation$genome_size, allgenes=Datapreperation$allgenes, datn=Datapreperation$datn, datt=Datapreperation$datt, counts_norm=Datapreperation$counts_nor, centromere_pos=Datapreperation$centromere_pos, chr_size=Datapreperation$chr_size)

# ratio tumor min door normaal
normal2 = Data_analysis$normal2
tumor2 = Data_analysis$tumor2
genome_size = Datapreperation$genome_size
chr_total = Datapreperation$chr_total
chr_size = Datapreperation$chr_size
centromere_pos = Datapreperation$centromere_pos

plottable <- data.frame(Chromosome=normal2$Chromosome, Position=normal2$Position, ExpressionN=normal2$Expression, ExpressionT=tumor2$Expression,ExpressionRatio=tumor2$Expression-normal2$Expression)
plottable$Position = as.character(plottable$Position)

# als de expressie van een gen zowel in normaal als in tumor 0 is dan verwijderen
#plottable<-plottable[!(plottable$ExpressionN==0 & plottable$ExpressionT==0),]

# Expressieplot: dinamische y-limiet
ymin = abs(min(plottable$ExpressionRatio))
ymin = ceiling(ymin)
ymin = -ymin
ymax = 1
ymax2 = max(plottable$ExpressionRatio)
ymax2 = ceiling(ymax2)
if(ymax2 > ymax){ymax = ymax2}else{ymax=ymax}

# plottable <- plottable[plottable$Chromosome != 23, ] 
# plottable <- plottable[plottable$Chromosome != 24, ]

setwd(output_dir)

#png(filename=paste("Comparisonplot3-",substr(code, 1,28), ".png", sep=""), width = 1980, height = 1080, units = "px")

par(mar=c(4, 7, 4, 2))
par(mfrow=c(3,1))
par(las=1)
plot(rollmedian(orderedTable$position,201),rollmedian(orderedTable$MajorMinor,201),col="554",pch=15,cex=0.4,ylim=c(1,Ylim), ylab="Allelic Ratio\n",typ="l",xlab="",xaxt = "n",xlim=c(1,genome_size), cex.axis=1.5,cex.lab=1.5, lwd=2.5)
abline(1.75,0,col="black")
# Draw chromosome guide lines + chromosoom cijfers
for(i in 1 :22){
  if (i>1){abline(v=chr_total[i],col="gray48")}
  abline(v=chr_total[i]+centromere_pos[i],col="gray55",lty=4)
  #mtext(chr_total[i]+centromere_pos[i], side=1,at=chr_total[i]+centromere_pos[i], cex=0.5)
  mtext(as.character(i),side=3, at=chr_total[i]+chr_size[i]/2,cex=1)
}
par(las=1)
plot(cnv_lg$Position,cnv_lg$Segment_Mean,col="28",pch=15,cex=0.4,type="l",xlim=c(1,genome_size),ylim=c(-1.5,1.5), ylab="Segment Mean\n",xlab="",xaxt = "n", cex.axis=1.5,cex.lab=1.5, lwd=2.5)
abline(0,0,col="black")
# Draw chromosome guide lines + chr cijfers
for(i in 1 :22){
  if (i>1){abline(v=chr_total[i],col="gray48")}
  abline(v=chr_total[i]+centromere_pos[i],col="gray55",lty=4)
  mtext(as.character(i),side=3, at=chr_total[i]+chr_size[i]/2,cex=1)
}
par(las=1)
plot(plottable$Position, plottable$ExpressionRatio,pch=15,cex=0.4, type = "l",col="554", xlab="Chromosomal position", ylab="Expression\n", xaxt = "n",xlim=c(1,genome_size), ylim=c(ymin,ymax), cex.axis=1.5,cex.lab=1.5, lwd=2.5)
abline(0,0,col="black")
# Draw chromosome guide lines + Chromosome cijfers
for(y in 1 :22){
  if (y>1){abline(v=chr_total[y],col="gray48")}
  abline(v=chr_total[y]+centromere_pos[y],col="gray55",lty=4)
  #mtext(chr_total[y]+centromere_pos[y], side=1,at=chr_total[y]+centromere_pos[y], cex=0.5)
  mtext(as.character(y),side=3, at=chr_total[y]+chr_size[y]/2,cex=1)
}

#dev.off()

