#########################
# Helper functions for 
# MicroExplorer Shiny App
#########################

# +++++++++++++++++++++++
# csv2dat
# +++++++++++++++++++++++
# returns countData, taxaData, sampleData
# from uploaded csv files

csv2dat <- function(countFilePath, taxaFilePath, sampleFilePath) {
  
  countData <- read_csv(countFilePath) %>% 
    rename(otu = 1) %>%
    column_to_rownames(var="otu")
  
  taxaData <- read_csv(taxaFilePath) %>%
    rename(otu = 1) %>%
    column_to_rownames(var="otu")
  
  sampleData <- read_csv(sampleFilePath) %>%
    rename(sname = 1) %>%
    column_to_rownames(var="sname")
  
  dataList <- list("countData" = countData,
                    "taxaData" = taxaData,
                    "sampleData"= sampleData)
  
  return(dataList)
}

# +++++++++++++++++++++++
# biom2dat
# +++++++++++++++++++++++
# returns countData, taxaData, sampleData
# from input BIOM file
biom2dat <- function(biomFilePath) {
}


# +++++++++++++++++++++++
# validateInputs
# +++++++++++++++++++++++
# Validate input files.
# checks if - 
#  i) sample names match in countData and sampleData 
#  AND
# ii) taxa names match in countData and taxaData

validateInputs <- function(countData, taxaData, sampleData) {
  
  snameCheck <- is_empty(setdiff(colnames(countData), rownames(sampleData))) &
                  is_empty(setdiff(rownames(sampleData), colnames(countData)))
  
  otuCheck <- is_empty(setdiff(rownames(countData), rownames(taxaData))) & 
                is_empty(setdiff(rownames(taxaData), rownames(countData)))
  
  if (snameCheck & otuCheck) {
    msg <- "Valid Input files! Please proceed to filtering step."
  } else if (snameCheck == FALSE) {
    msg <- "ERROR: Sample names do not match in count data and sample data. Please re-upload your data."
  } else if (outCheck == FALSE) {
    msg <- "ERROR: OTU IDs do not match in count data and taxa data. Please re-upload your data"
  } else {
    msg <- "ERROR: Unable to validate input files. Please re-upload your data."
  }
  
  dataList = list("msg" = msg,
                 "countData" = countData,
                 "taxaData" = taxaData,
                 "sampleData" = sampleData
                 )

  return(dataList) 
}


# +++++++++++++++++++++++
# filterData
# +++++++++++++++++++++++
# filters the valid data based on user inputs: 
#   i) min seq depth
#  ii) min taxa abundance
# iii) min taxa prevalence

filterData <- function(countData, taxaData, sampleData, seqDepth, minTaxaAbund, minTaxaPrev) {
  
  # remove samples with seq depth below threshold
  countDataMinDepth <- countData %>%
    select(names(which(colSums(countData) >= seqDepth)))
  
  # transform counts to percentage 
  countDataMinDepthPerc <- countDataMinDepth * 100 / colSums(countDataMinDepth)[col(countDataMinDepth)]
  
  # filter low abundance/prevalence taxa
  taxa2keep <- replace(countDataMinDepthPerc, countDataMinDepthPerc < minTaxaAbund, NA) %>%
    rownames_to_column(var="otu") %>%
    filter(rowSums(is.na(.))/ncol(.) * 100.0 < (100.0 - minTaxaPrev)) %>%
    column_to_rownames(var="otu") %>%
    rownames()

  # filtered data
  countDataFiltered <- countDataMinDepth %>%
    rownames_to_column(var="otu") %>%
    filter(otu %in% taxa2keep) %>%
    column_to_rownames(var="otu")
  taxaDataFiltered <- taxaData %>%
    rownames_to_column(var="otu") %>%
    filter(otu %in% taxa2keep) %>%
    column_to_rownames(var="otu")
  sampleDataFiltered <- sampleData %>%
    rownames_to_column(var="sample") %>%
    filter(sample %in% colnames(countDataMinDepth)) %>%
    column_to_rownames(var="sample")
  
  dataList = list("countData" = countDataFiltered,
                  "taxaData" = taxaDataFiltered,
                  "sampleData" = sampleDataFiltered
  )
  
  return(dataList) 
}


# +++++++++++++++++++++++
# sortSamplesByDissimilarity
# +++++++++++++++++++++++
# sorts samples by dissimilarity  
# cluster bray dissimilarity distances using ward.D2 method
# input proportional data: Samples in rows, Taxa in column, fill be relative abundance
# returns sample order
sortSamplesByDissimilarity <- function(propData) {
  bcdist <- vegdist(propData, method="bray")
  hclustBC <- hclust(bcdist, method="ward.D2")
  sorder <- hclustBC$labels[c(hclustBC$order)]
  return(sorder)
}


# +++++++++++++++++++++++
# sortSamplesByDescTaxaAbund
# +++++++++++++++++++++++
# sorts samples by decreasing taxa abundances  
# returns sample order
sortSamplesByDescTaxaAbund <- function(propData) {
  ## identify topTaxa
  maxTaxa <- apply(propData, 1, which.max)
  maxProp <- propData[matrix(c(1:nrow(propData),maxTaxa), ncol=2)]
  maxTaxa <- colnames(propData)[maxTaxa]
  maxTaxaDF <- data.frame(cbind(rownames(propData), maxTaxa))
  #maxTaxa.uniq <- unique(maxTaxaDF$maxTaxa)
  tbl <- data.frame(table(maxTaxaDF$maxTaxa))
  taxaOrder <- tbl[order(-tbl$Freq),]
  taxaOrder <- as.character(taxaOrder$Var1)
  taxaOrder <- c(taxaOrder[!(taxaOrder %in% c("OtherTaxa", "unclassified"))], "OtherTaxa", "unclassified") 
  
  # set default sampleOrder
  sampleOrder <- vector(mode="character", length=0)
  
  for (x in taxaOrder) {
    samples <- maxTaxaDF[maxTaxaDF$maxTaxa == x,]
    if (length(samples$V1) > 1 ){
      tmp <- propData[rownames(propData) %in% samples$V1,]
      tmp.melt <- melt(tmp)
      tmp.melt <- tmp.melt[order(-tmp.melt$value, tmp.melt$Var2, tmp.melt$Var1),]
      tmp.melt$Var1 <- factor(tmp.melt$Var1, levels = unique(tmp.melt$Var1))
      if ( x != "unclassified"){
        sampleOrder <- c(sampleOrder, levels(tmp.melt$Var1))
      } else {
        sampleOrder <- c(sampleOrder, rev(levels(tmp.melt$Var1)))
      }  
    } else {
      sampleOrder <- c(sampleOrder, as.character(samples$V1))
    }
  }
  
  ### return the custom sample order
  return(sampleOrder)
}



# +++++++++++++++++++++++
# myStackedBarPlot
# +++++++++++++++++++++++
# This function generates a stacked bar plot
# of community composition
myStackedBarPlot <- function(countData, taxaData, sampleData, 
                             taxaLevel, taxa2Plot, numTaxa2Plot = NULL, 
                             sortMethod, facetField = "None"){
  
  # get taxa level
  mytaxaData <- taxaData %>%
    rownames_to_column(var="otu") %>%
    select(otu, taxaLevel)
  # create metadata
  metadata <- sampleData %>%
    rownames_to_column(var="Sample")
  
  # create the plot data frame
  plotDF <- countData %>%
    rownames_to_column(var="otu") %>%
    melt(.) %>%
    left_join(mytaxaData, by="otu") %>%
    set_colnames(c("OTU", "Sample", "Count", "Taxa")) %>%
    mutate(Taxa2 = ifelse(grepl("unclassified", Taxa), "unclassified", Taxa)) %>%
    group_by(Sample, Taxa2) %>%
    tally(Count) %>%
    mutate(RelAb = n / sum(n) * 100.0) %>%
    left_join(., metadata, by="Sample")
  
  # propData
  propData <- acast(plotDF, Sample~Taxa2, value.var = "RelAb", fill=0.0)
  
  # get sampleOrder
  if (sortMethod == "Cluster by Dissimilarity"){
    sorder <- sortSamplesByDissimilarity(propData)
  }
  # set sample levels
  plotDF$Sample <- factor(plotDF$Sample, levels = sorder)
  
  # plot
  p <- ggplot(plotDF, aes(Sample, RelAb, fill=Taxa2)) + 
        geom_bar(stat="identity", position="stack") + 
        theme_bw()
  
  # add facet
  if (!(facetField == "None")) {
    p <- p + facet_wrap(~eval(parse(text=facetField)), scales = "free_x")
  }
  
  ggplotly(p)
}


# +++++++++++++++++++++++
# plotSeqDepth
# +++++++++++++++++++++++
plotSeqDepth <- function(countData) {
  seqDepth <- colSums(countData)
  hist(seqDepth, breaks = 20)
}
