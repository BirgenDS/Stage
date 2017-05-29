library("zoo")
library("gplots")
library(tools)
library(ggplot2)
library(gridExtra)
library(tidyr)

Window = 151
Ylim = 3
Organism == "Human"

par(mar=c(5, 4, 4, 2))
window=Window

if(Organism == "Human"){
  centromere_pos=c(125,93.3,91,50.4,48.4,61,59.9,45.6,49,40.2,53.7,35.8,17.9,17.6,19,36.6,24,17.2,26.5,27.5,13.2,14.7,60.6,12.5)
  centromere_pos=centromere_pos*1000000
  chr_size=c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,155270560,59373566)
}

chr_total=0
for (i in 1:(length(chr_size)-1)){
  chr_total=c(chr_total,sum(chr_size[1:i]))  
}
genome_size=sum(chr_size)

setwd("~/Documents/Howest/Stage/CMGG/RNA-seq_CNV/BRCA")
#setwd("~/Documents/Howest/Stage/CMGG/RNA-seq_CNV/THYM")

datn <- read.table("BRCA_n_geneENS_tpm.txt",header=TRUE,sep="\t")
datt <- read.table("BRCA_t_geneENS_tpm.txt",header=TRUE,sep="\t")

#datn <- read.table("THYM_n_geneENS_tpm.txt",header=TRUE,sep="\t")
#datt <- read.table("THYM_t_geneENS_tpm.txt",header=TRUE,sep="\t")

setwd("~/Documents/Howest/Stage/CMGG/RNA-seq_CNV/BRCA/Gene-tabellen")
#setwd("~/Documents/Howest/Stage/CMGG/RNA-seq_CNV/THYM/Gene-tabellen")

file_list = list.files(pattern="*.csv")
data_list <- vector("list", "length" = length(file_list))  

for(i in seq_along(file_list)){
  i = 1
  filename = file_list[[i]]
  data <- read.csv(filename, header = TRUE, sep = "\t")
  
  code <- file_path_sans_ext(filename)

  # NORMAL
  # dataframe from datn with only this code
  
  tryCatch({
    j <- grep(code, colnames(datn))
    
    exprn <- data.frame(datn$X,datn[,j])
    colnames(exprn) <- c('GeneId', 'Expression')
    
    # searching for GeneId and merging the expression to data
    genes <- exprn[exprn$GeneId %in% data$GeneId,]
    genesn <- merge(data, genes, all = TRUE, sort = FALSE)
    genesn <- genesn[order(genesn[,2]),]
    
    genesn <- genesn[c("Chromosoom", "Start", "Stop","GeneId", "Expression")]
    
    setwd("~/Documents/Howest/Stage/CMGG/RNA-seq_CNV/BRCA/Expressietabellen")
    #setwd("~/Documents/Howest/Stage/CMGG/RNA-seq_CNV/THYM/Expressietabellen")
    #write.table(genesn, file=paste(code, "_normal.csv", sep=""),sep = "\t", col.names  = TRUE, row.names = FALSE)
    
  }, error=function(e){cat("ERROR :",code,"not in normal data.\n")
  }, finally = {
    start = genesn[,c('Chromosoom', 'Start',  'GeneId','Expression')]
    stop = genesn[,c('Chromosoom', 'Stop', 'GeneId', 'Expression')]
    plottable1<-transform(genesn, Region=paste(Start, Stop, sep="-"))
    plottable1<- plottable1[,c('Chromosoom', 'Region', 'GeneId','Expression')]
  })

  # TUMOR
  tryCatch({
    k <- grep(code, colnames(datt))
    exprt <- data.frame(datt$X,datt[,k])
    
    # bij meerdere expressions voor 1 staal het gemiddelde nemen
    if(ncol(exprt) > 2){
      exprt$Mean <- rowMeans(exprt[,2:ncol(exprt)],na.rm = TRUE)
      exprt$Mean <- round(exprt$Mean, digits = 6)
      exprt <- data.frame(exprt$datt.X,exprt$Mean)
    }
    
    colnames(exprt) <- c('GeneId', 'Expression')
    
    genes <- exprt[exprt$GeneId %in% data$GeneId,]
    genest <- merge(data, genes, all = TRUE, sort = FALSE)
    genest <- genest[order(genest[,2]),]
    
    genest <- genest[c("Chromosoom", "Start", "Stop","GeneId", "Expression")]
    
    setwd("~/Documents/Howest/Stage/CMGG/RNA-seq_CNV/BRCA/Expressietabellen")
    #setwd("~/Documents/Howest/Stage/CMGG/RNA-seq_CNV/THYM/Expressietabellen")
    #write.table(genest, file=paste(code, "_tumor.csv", sep=""),sep = "\t", col.names  = TRUE, row.names = FALSE)
    
  }, error=function(e){cat("ERROR : ",code,"  not in tumor data.\n")
  }, finally = {
    start2 = genest[,c('Chromosoom', 'Start',  'GeneId','Expression')]
    stop2 = genest[,c('Chromosoom', 'Stop', 'GeneId', 'Expression')]
    plottable2<-transform(genest, Region=paste(Start, Stop, sep="-"))
    plottable2<- plottable2[,c('Chromosoom', 'Region', 'GeneId','Expression')]
  })

  setwd("~/Documents/Howest/Stage/CMGG/RNA-seq_CNV/BRCA/Gene-tabellen")
  #setwd("~/Documents/Howest/Stage/CMGG/RNA-seq_CNV/THYM/Gene-tabellen")

  # als tumor en normaal hebben, maak dan de plots
  if(plottable1$GeneId == plottable2$GeneId) {
    
    #tabel maken met gemiddelde expressiewaarden per regio
    plottablen <- data.frame(Chromosome=plottable1$Chromosoom[1], Region=plottable1$Region[1], GeneId=plottable1$GeneId[1] ,Expression=c(NA))
    Expressionlist <- c(plottable1$Expression[1])
    i = 2
    for (i in 2:nrow(plottable1)){
      if (plottable1$Region[i-1] != plottable1$Region[i]) {
        plottablen$Expression[is.na(plottablen$Expression)] = mean(Expressionlist)
        newrow = c(plottable1$Chromosoom[i], as.character(plottable1$Region[i]),as.character(plottable1$GeneId[i]),plottablen$Expression[NA])
        plottablen <- rbind(plottablen,newrow)
        Expressionlist <- c()
      }
      Expressionlist <- append(Expressionlist,plottable1$Expression[i])
      i = i+1
      # nog toevoegen dat hij bij de laatste loop ook de mean van expressionlist toevoegd aan de dataframe
    }
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
    plottablen <- plottablen[order(plottablen[,1]),]
    
    # alle chormosomen in tabel zetten, als er geen expressie van is, dan is deze 0
    for(x in 1 :24){
      x = 1
      if(!(x %in% plottablen$Chromosome)){
        newerrow = c(plottablen$Chromosome <- x, plottablen$Expression <- 0)
        plottablen2 <- rbind(plottablen,newerrow)
      }
    }
    
    #png(filename=paste("Expressiemeanplot-",code, ".png", sep=""), width = 1980, height = 1080, units = "px")
    #par(mfrow=c(2,1))
    plot(plottablen$Position, plottablen$Expression, type = "l")
    # Draw chromosome guide lines + chromosoom cijfers
    for(y in 1 :24){
      if (y>1){abline(v=chr_total[y],col="gray48")}
      abline(v=chr_total[y]+centromere_pos[y],col="gray55",lty=4)
      #mtext(chr_total[y]+centromere_pos[y], side=1,at=chr_total[y]+centromere_pos[y], cex=0.5)
      mtext(as.character(y),side=3, at=chr_total[y]+chr_size[y]/2,cex=0.8)
    }
    
    #dev.off()
  }
  
  i = i+1
}
