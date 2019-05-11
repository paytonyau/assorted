#########################################################
### A) Installing and loading required packages
#########################################################

if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}


#########################################################
### B) Reading in data and transform it into matrix format
#########################################################

data <- read.csv("C:/miRNA/Heatmap_CRC.csv")
rnames <- data[,1]                            # assign labels in column 1 to "rnames"
mat_data <- data.matrix(data[,2:ncol(data)])  # transform column 2-5 into a matrix
rownames(mat_data) <- rnames                  # assign row names


#########################################################
### C) Customizing and plotting the heat map
#########################################################

# creates a own color palette from red to green
my_palette <- colorRampPalette(c("green", "yellow", "red"))(n = 299)

# (optional) defines the color breaks manually for a "skewed" color transition
col_breaks = c(seq(0,2,length=100),  # for red
               seq(2.1,9.9,length=100),           # for yellow
               seq(10,99,length=100))             # for green

# creates a 5 x 5 inch image
jpeg("C:/miRNA/Heatmap.jpeg",    # create jpeg for the heat map        
    width = 5*400,        # 5 x 300 pixels
    height = 5*400,
    res = 1000,            # 300 pixels per inch
    pointsize = 2)        # smaller font size

## Example A
heatmap.2(mat_data,
          main = "Correlation", # heat map title
          density.info="none",  # turns off density plot inside color legend
          margins =c(10,20),      # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="row",     # only draw a row dendrogram
          cexRow = 1,         # row label font size
          cexCol = 1,         #col label font size
          lhei = c(0.1,1),      #colour key size
          lwid = c(0.25,1),     #colour key size
          cellnote = mat_data,  # same data set for cell labels
          breaks=col_breaks,    # enable color transition at specified limits
          notecol="black",      # change font color of cell labels to black
          symkey=FALSE,
          srtCol=0,
          trace="none",
          Colv="NA")            # turn off column clustering

## Example B
heatmap.2(mat_data, 
          col=redgreen(75), 
          scale="row", 
          key=TRUE, 
          symkey=FALSE,
          density.info="none",
          cexRow=1,
          cexCol=1,
          margins=c(6,11),  
          trace="none",
          srtCol=90)

dev.off()               # close the jpeg device
