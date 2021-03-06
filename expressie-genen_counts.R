library("zoo")
library("gplots")
library(tools)
library(gridExtra)
library(tidyr)

source("http://www.bioconductor.org/biocLite.R")
biocLite("limma")
biocLite("statmod")
biocLite("edgeR")
biocLite("DESeq")
library("statmod")
library("DESeq")

Organism = "Human"

par(mar=c(5, 4, 4, 2))
window=151

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

setwd("~/Documents/Howest/Stage/CMGG/RNA-seq_CNV/BRCA")

datn <- read.delim("BRCA_n_geneENS_counts.txt", row.names = 1)
datt <- read.delim("BRCA_t_geneENS_counts.txt", row.names = 1)

# aantal tumor stalen in de normale data, dus alle stalen met 01B in het midden mogen verwijderd 
datn <- datn[, -grep("01B", colnames(datn))]

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

# lijst van alle genen
setwd("~/Documents/Howest/Stage/CMGG/RNA-seq_CNV")
allgenes <- read.table("allgenes.txt",header=TRUE,sep="\t")
allgenes <- allgenes[order(allgenes[,1],allgenes[,2]),]
colnames(allgenes) <- c('Chromosome', 'Start', 'Stop', 'GeneId')

# allgenes begint van 0, in deze file telt de positie op, dus chromosoom grote bij optellen
for(x in 1:22){
  start = allgenes$Start[allgenes$Chromosome == x]
  start = as.numeric(start)
  allgenes$Start[allgenes$Chromosome == x] = start + chr_total[x]
}
for(x in 1:22){
  stop = allgenes$Stop[allgenes$Chromosome == x]
  stop = as.numeric(stop)
  allgenes$Stop[allgenes$Chromosome == x] = stop + chr_total[x]
}

# plaats van de code in de data frame zoeken
setwd("~/Documents/Howest/Stage/CMGG/RNA-seq_CNV/BRCA/Gene-tabellen")

file_list = list.files(pattern="*.csv")
data_list <- vector("list", "length" = length(file_list))  

# nieuwe dataframe maken met kolommen met huidige code
for(i in seq_along(file_list)){
  i = 1
  
  filename = file_list[[i]]
  data <- read.csv(filename, header = TRUE, sep = "\t")
  
  code <- file_path_sans_ext(filename)
  
  # # allgenes begint van 0, in deze file telt de positie op, dus berekenen -> allgenes aangepast
  # data1 <- data.frame(Chromosome=data$Chromosome[1], Start=c(NA), Stop=c(NA))
  # for(chr in 1:nrow(data)){
  #   chrom = data$Chromosome[chr]
  #   number1 = data$Start[chr]
  #   number2 = data$Stop[chr]
  #   data1$Start[is.na(data1$Start)] = number1 - chr_total[chrom]
  #   data1$Stop[is.na(data1$Stop)] = number2 - chr_total[chrom]
  #   newrow = c(data$Chromosome[chr+1], data1$Start[NA], data1$Stop[NA])
  #   data1 <- rbind(data1,newrow)
  # }
  # data1 <- na.omit(data1)
  
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
  # dataframe from datn with only this code
  tryCatch({
    code1 <- substr(code, 6,17)
    j <- grep(code1, colnames(datn))
    listn <- paste("normal",j)
    
    col.num <- which(colnames(counts_norm) %in% listn)
    exprn <- counts_norm[,col.num]
    
    # eerst overal +1 doen en dan de log nemen
    exprn <- exprn + 1
    exprn[,1:ncol(exprn)] <- log(exprn[1:ncol(exprn)], 2)
    
    if(ncol(exprn) > 1){
      exprn$Mean <- rowMeans(exprn[,1:ncol(exprn)],na.rm = TRUE)
      exprn$Mean <- round(exprn$Mean, digits = 6)
      `%ni%` <- Negate(`%in%`)
      exprn <- subset(exprn,select = names(exprn) %ni% listn)
    }
    
    exprn <- setDT(exprn, keep.rownames = TRUE)[]
    
    colnames(exprn) <- c('GeneId', 'Expression')
    
    # searching for GeneId and merging the expression to data
    genes <- exprn[exprn$GeneId %in% results$GeneId,]
    genesn <- merge(results, genes, all = TRUE, sort = FALSE)
    #genesn <- genesn[order(genesn[,2]),]
    genesn <- na.omit(genesn)
    
    genesn <- genesn[c("Chromosome", "Start", "Stop","GeneId", "Expression")]
    
    #setwd("~/Documents/Howest/Stage/CMGG/RNA-seq_CNV/KIRC/Expressietabellen")
    #write.table(genesn, file=paste(code, "_normal.csv", sep=""),sep = "\t", col.names  = TRUE, row.names = FALSE)
    
  }, error=function(e){cat("ERROR :",code,"not in normal data.\n")
  }, finally = {
    start = genesn[,c('Chromosome', 'Start',  'GeneId','Expression')]
    stop = genesn[,c('Chromosome', 'Stop', 'GeneId', 'Expression')]
    plottable1<-transform(genesn, Region=paste(Start, Stop, sep="-"))
    plottable1<- plottable1[,c('Chromosome', 'Region', 'GeneId','Expression')]
    plottable1$Region <- as.character(plottable1$Region)
  })
  
  # TUMOR
  tryCatch({
    code2 <- substr(code, 6,21)
    k <- grep(code2, colnames(datt))
    listt <- paste("tumor",k)
    
    col.numb <- which(colnames(counts_norm) %in% listt)
    exprt <- counts_norm[,col.numb]
    
    # eerst overal +1 doen en dan de log nemen
    exprt <- exprt + 1
    exprt[,1:ncol(exprt)] <- log(exprt[1:ncol(exprt)], 2)
    
    if(ncol(exprt) > 1){
      exprt$Mean <- rowMeans(exprt[,1:ncol(exprt)],na.rm = TRUE)
      exprt$Mean <- round(exprt$Mean, digits = 6)
      `%ni%` <- Negate(`%in%`)
      exprt <- subset(exprt,select = names(exprt) %ni% listt)
    }
    
    exprt <- setDT(exprt, keep.rownames = TRUE)[]
    
    colnames(exprt) <- c('GeneId', 'Expression')
    
    # searching for GeneId and merging the expression to data
    genes <- exprt[exprt$GeneId %in% results$GeneId,]
    genest <- merge(results, genes, all = TRUE, sort = FALSE)
    #genest <- genest[order(genest[,2]),]
    genest <- na.omit(genest)
    
    genest <- genest[c("Chromosome", "Start", "Stop","GeneId", "Expression")]
    
    #setwd("~/Documents/Howest/Stage/CMGG/RNA-seq_CNV/KIRC/Expressietabellen")
    #write.table(genest, file=paste(code, "_tumor.csv", sep=""),sep = "\t", col.names  = TRUE, row.names = FALSE)
    
  }, error=function(e){cat("ERROR : ",code,"not in tumor data.\n")
  }, finally = {
    start2 = genest[,c('Chromosome', 'Start',  'GeneId','Expression')]
    stop2 = genest[,c('Chromosome', 'Stop', 'GeneId', 'Expression')]
    plottable2<-transform(genest, Region=paste(Start, Stop, sep="-"))
    plottable2<- plottable2[,c('Chromosome', 'Region', 'GeneId','Expression')]
    plottable2$Region <- as.character(plottable2$Region)
  })

  setwd("~/Documents/Howest/Stage/CMGG/RNA-seq_CNV/BRCA/Gene-tabellen")
  
  # als tumor en normaal hebben, maak dan de plots
  if(plottable1$GeneId == plottable2$GeneId) {
    
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
    
    plottablen <- separate(data = plottablen, col = Region, into = c("Start", "Stop"), sep = "\\-")
    
    start = plottablen[,c('Chromosome', 'Start', 'Expression')]
    colnames(start) <- c('Chromosome', 'Position', 'Expression')
    stop = plottablen[,c('Chromosome', 'Stop', 'Expression')]
    colnames(stop) <- c('Chromosome', 'Position', 'Expression')
    
    plottablen<-as.data.frame(rbind(start, stop)) # plakt stop kolom onder start
    plottablen$Position <- as.numeric(plottablen$Position)
    plottablen <- plottablen[order(plottablen[,2],plottablen[,1]),]
    
    # alle chormosomen in tabel zetten en ontbrekende stukken van bestaande chromosomen
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
    
    # ratio tumor min door normaal
    plottable <- data.frame(Chromosome=normal2$Chromosome, Position=normal2$Position, ExpressionN=normal2$Expression, ExpressionT=tumor2$Expression,ExpressionRatio=tumor2$Expression-normal2$Expression)
    #plottable <- replace(plottable, is.na(plottable), 0)
    #plottable <- do.call(data.frame,lapply(plottable, function(x) replace(x, is.infinite(x),0)))
    
    plottable <- plottable[plottable$Chromosome != 23, ] 
    plottable <- plottable[plottable$Chromosome != 24, ]
    
    #png(filename=paste("Expressiemeanplot_counts-",code, ".png", sep=""), width = 1980, height = 1080, units = "px")
    par(mfrow=c(3,1))
    plot(plottable$Position, plottable$ExpressionN,pch=15,cex=0.4, type = "l",col="dimgrey", main = "Normal", xlab="Chromosomal position", ylab="Expression", xaxt = "n",xlim=c(1,genome_size), ylim=c(-1,5)) #
    # Draw chromosome guide lines + Chromosome cijfers
    for(y in 1 :22){
      if (y>1){abline(v=chr_total[y],col="gray48")}
      abline(v=chr_total[y]+centromere_pos[y],col="gray55",lty=4)
      #mtext(chr_total[y]+centromere_pos[y], side=1,at=chr_total[y]+centromere_pos[y], cex=0.5)
      mtext(as.character(y),side=3, at=chr_total[y]+chr_size[y]/2,cex=0.8)
    }
    
    plot(plottable$Position, plottable$ExpressionT,pch=15,cex=0.4, type = "l",col="dimgrey", main = "Tumor", xlab="Chromosomal position", ylab="Expression", xaxt = "n",xlim=c(1,genome_size), ylim=c(-1,5)) #
    # Draw chromosome guide lines + Chromosome cijfers
    for(y in 1 :22){
      if (y>1){abline(v=chr_total[y],col="gray48")}
      abline(v=chr_total[y]+centromere_pos[y],col="gray55",lty=4)
      #mtext(chr_total[y]+centromere_pos[y], side=1,at=chr_total[y]+centromere_pos[y], cex=0.5)
      mtext(as.character(y),side=3, at=chr_total[y]+chr_size[y]/2,cex=0.8)
    }
    
    plot(plottable$Position, plottable$ExpressionRatio,pch=15,cex=0.4, type = "l",col="dimgrey", main = "Ratio", xlab="Chromosomal position", ylab="Expression", xaxt = "n",xlim=c(1,genome_size), ylim=c(-1,5)) #
    # Draw chromosome guide lines + Chromosome cijfers
    for(y in 1 :22){
      if (y>1){abline(v=chr_total[y],col="gray48")}
      abline(v=chr_total[y]+centromere_pos[y],col="gray55",lty=4)
      #mtext(chr_total[y]+centromere_pos[y], side=1,at=chr_total[y]+centromere_pos[y], cex=0.5)
      mtext(as.character(y),side=3, at=chr_total[y]+chr_size[y]/2,cex=0.8)
    }
    
    #dev.off()
  }
  
  #i = i+1

}
