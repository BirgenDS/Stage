library(tools)
library(ggplot2)
library(gridExtra)

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
    
  }, error=function(e){cat("ERROR :",code,"not in normal data.\n")})

  # TUMOR
  tryCatch({
    k <- grep(code, colnames(datt))
    
    exprt <- data.frame(datt$X,datt[,k])
    colnames(exprt) <- c('GeneId', 'Expression')
    
    genes <- exprt[exprt$GeneId %in% data$GeneId,]
    genest <- merge(data, genes, all = TRUE, sort = FALSE)
    genest <- genest[order(genest[,2]),]
    
    genest <- genest[c("Chromosoom", "Start", "Stop","GeneId", "Expression")]
    
    setwd("~/Documents/Howest/Stage/CMGG/RNA-seq_CNV/BRCA/Expressietabellen")
    #setwd("~/Documents/Howest/Stage/CMGG/RNA-seq_CNV/THYM/Expressietabellen")
    #write.table(genest, file=paste(code, "_tumor.csv", sep=""),sep = "\t", col.names  = TRUE, row.names = FALSE)
    
  }, error=function(e){cat("ERROR : ",code,"  not in tumor data.\n")})

  setwd("~/Documents/Howest/Stage/CMGG/RNA-seq_CNV/BRCA/Gene-tabellen")
  #setwd("~/Documents/Howest/Stage/CMGG/RNA-seq_CNV/THYM/Gene-tabellen")
  
  start = genesn[,c('Chromosoom', 'Start',  'GeneId','Expression')]
  stop = genesn[,c('Chromosoom', 'Stop', 'GeneId', 'Expression')]
  plottable1<-transform(genesn, Region=paste(Start, Stop, sep="-"))
  plottable1<- plottable1[,c('Chromosoom', 'Region', 'GeneId','Expression')]
  
  start2 = genest[,c('Chromosoom', 'Start',  'GeneId','Expression')]
  stop2 = genest[,c('Chromosoom', 'Stop', 'GeneId', 'Expression')]
  plottable2<-transform(genest, Region=paste(Start, Stop, sep="-"))
  plottable2<- plottable2[,c('Chromosoom', 'Region', 'GeneId','Expression')]
  
  # plot 
  plot1 <- ggplot(data = plottable1, aes(x=GeneId, y=Expression, cex.axis(0.7), group=1))+
    geom_line()+
    labs(title = "Normal")
  #print(plot1)
  
  plot2 <- ggplot(data = plottable2, aes(x=GeneId, y=Expression, cex.axis(0.7), group=1))+
    geom_line()+
    labs(title = "Tumor")
  #print(plot2)
  
  png(filename=paste("Expressieplot-",code, ".png", sep=""), width = 1980, height = 1080, units = "px")
  grid.arrange(plot1, plot2, nrow=2, ncol=1)
  dev.off()
  
  # nog oplossing zoeken voor als hij alleen maar een tumor tabel kan maken en geen normaal
  
  i = i+1
}


