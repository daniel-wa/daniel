# Here is a typical way to create a colored Heatmap with the R package: heatmap3

library(svDialogs) # used to select options
library(taxadb)
library(RColorBrewer)
library(makeFlow)
library("seqinr")

# to_num takes character vectors and converts gaps ("-") to zeroes.
# Anything else is a 1.0
to_num <- function(z) {
  r <- rep(0, length(z))
  r[z != "-"]  <-  1
  return(r)
}

# presence_names goes to a named file with fasta names.
# It returns the names as a list and a structure, presence which is the
# sequence alignment converted to presence / absence.
presence_names <- function(fasta_name) {
  a <- seqinr::read.fasta (fasta_name)
  presence <- t(sapply(a, to_num)) # convert to numbers and transpose
  row.names(presence)  <-  as.vector(seqinr::getAnnot(a), mode="character")
  
  nameset <- row.names (presence)
  #   Remove everything that is in front of the species name and behind it
  nameset <- sub(".*\\[", "", nameset) 
  nameset <- sub("\\]", "", nameset)
  nameset <- tolower(nameset)
  
  return (list(nameset=nameset, presence=presence))
}

mymain <- function() {
  Sys.setenv(CONTENTID_REGISTRIES="https://hash-archive.carlboettiger.info")
  
  
  filepathway <- dlg_open(title = "Select a Fasta file.")$res
  if (nchar (filepathway) < 1) {
    stop ('Did not get a valid path for sequence file')
  }
  
  t <- presence_names (filepathway)
  nameset = t$nameset
  presence = t$presence
  rm(filepathway)
  
  #   Filter_name() maps the taxonomic information to the nameset. 
  #   Since NCBI has more unique mappings, we use NCBI
  
  taxnameset <- filter_name(nameset, provider = getOption("taxadb_default_provider","ncbi")) 
  
  
  #   Ask which taxonomic rank he would like to select (and checked for typos)
  taxrank <- ""
  while (taxrank != "kingdom" && taxrank != "phylum" && taxrank != "class" && taxrank != "order" && taxrank != "family" && taxrank != "genus"){
    taxrank <- dlgInput("Which taxonomic rank do you want to select? Please choose:
                      \n'kingdom', 'phylum', 'class', 'order', 'family', 'genus'.")$res
  }
  
  #   Remove all columns from taxanameset except input (species name), Sort (shows duplicates)
  #   and chosen taxrank 
  taxnameset <- taxnameset[ ,c("sort","input",taxrank)]
  
  #   Remove duplicates by using sort as primary key 
  taxnameset <- taxnameset[!duplicated(taxnameset$sort), ]
  
  #   delete "sort" column
  nametaxa <- taxnameset[ ,c("input",taxrank)]
  rm(taxnameset)
  
  #   sort the groups from most seen group to less seen group
  nametaxacopy <- data.frame(nametaxa[complete.cases(nametaxa),])
  names(nametaxacopy)[names(nametaxacopy) == taxrank] <- "taxlevel"
  nametaxacopy <- transform(nametaxacopy, freq= ave(seq(nrow(nametaxacopy)), taxlevel, FUN=length))
  nametaxacopy <- nametaxacopy[order(nametaxacopy$freq, decreasing = TRUE),]
  names(nametaxacopy)[names(nametaxacopy) == "taxlevel"] <- taxrank
  
  #   Look up how many (unique) groups there are in chosen taxrank for selected Fasta file
  groupset <- data.frame(unique(nametaxacopy[1:nrow(nametaxacopy),2]))
  rm(nametaxacopy)
  
  #   How many groups should be displayed from your taxonomic selection
  if (nrow(groupset) > 50) {
    telluser <- paste("How many groups out of", nrow(groupset), 
                      "should be displayed from your taxonomic selection? (max.49 and > 11 not for Colorblind)") 
  } else if (nrow(groupset) > 12) {
    telluser <- paste("How many groups out of", nrow(groupset), 
                      "should be displayed from your taxonomic selection? (> 11 not for Colorblind)") 
  } else {
    telluser <- paste("How many groups out of", nrow(groupset), 
                      "should be displayed from your taxonomic selection?")
  }
  
  groups <- dlgInput(telluser)$res
  groups <- as.numeric(groups)
  
  #   Check for typos
  while (is.na(groups) == T || groups < 1 || groups > 49 || groups > nrow(groupset) ||round(groups) != groups){
    groups <- dlgInput(telluser)$res
    groups <- as.numeric(groups)
  } 
  rm(telluser)
  
  #   Only if the user does not want to see all groups of his selected taxrank, 
  #   the remaining groups will be renamed as "Others" 
  if (groups != nrow(groupset)){
    groups = groups + 1
    for (e in groups:nrow(groupset)){
      nametaxa[nametaxa == groupset[e,1]] <- "Others"
    }
    groups = groups - 1
    rm(e)
  }
  
  #Should the "Others" group be even displayed?
  if (groups != nrow(groupset)){
    other <- ""
    while (other != "yes" && other != "no"){
    other <- dlgInput("Should the group 'Others' be displayed? 'yes' or 'no'?")$res
    }
    if(other == "no"){
    nametaxa[nametaxa == "Others"] <- NA
    rm(other)
    }
  }
  
  #   Transparency of colors is chosen by length of nameset
  if(length(nameset) > 500) {
    transparency <- 0.12
    cexRow <- 0.8
  } else if (length(nameset) > 20 && length(nameset) < 501) {
    transparency <- 0.2
    cexRow <- 1
  } else {
    transparency <- 0.3
    cexRow <- 1.2
  }
  
  #   Here all groups are sorted again, because the user could have changed some groups to "Others"  
  nametaxacopy <- data.frame(nametaxa[complete.cases(nametaxa),])
  names(nametaxacopy)[names(nametaxacopy) == taxrank] <- "taxlevel"
  nametaxacopy <- transform(nametaxacopy, freq= ave(seq(nrow(nametaxacopy)), taxlevel, FUN=length))
  nametaxacopy <- nametaxacopy[order(nametaxacopy$freq, decreasing = TRUE),]
  names(nametaxacopy)[names(nametaxacopy) == "taxlevel"] <- taxrank
  groupset <- data.frame(unique(nametaxacopy[1:nrow(nametaxacopy),2]))
  rm(taxrank, nametaxacopy) 
  
  #   Here some parameters for the legend and colors are saved depending on chosen number for groups
  if (groups > 40 && groups < 51){
    color <- makeFlow:::addAlpha(rainbow(nrow(groupset)), transparency)
    max_name_length <- 19
    cex <- 0.55
    ncol <- 5
    lwd <- 5
  } else if (groups > 30 && groups < 41){
    color <- makeFlow:::addAlpha(rainbow(nrow(groupset)), transparency)
    max_name_length <- 19
    cex <- 0.59
    ncol <- 5
    lwd <- 5
  } else if (groups > 20 && groups < 31){
    color <- makeFlow:::addAlpha(rainbow(nrow(groupset)), transparency)
    max_name_length <- 25
    cex <- 0.66
    ncol <- 4
    lwd <- 5
  } else if (groups > 12 && groups < 21){
    color <- makeFlow:::addAlpha(rainbow(nrow(groupset)), transparency)
    max_name_length <- 25
    cex <- 0.8
    ncol <- 3
    lwd <- 5
  } else if (groups > 7 && groups < 13){
    color <- makeFlow:::addAlpha(brewer.pal(nrow(groupset),"Paired"), transparency)
    max_name_length <- 30
    cex <- 0.8
    ncol <- 2
    lwd <- 7
    
  } else {
    color <- makeFlow:::addAlpha(brewer.pal(nrow(groupset),"Dark2"), transparency)
    max_name_length <- 35
    cex <- 0.8
    ncol <- 1
    lwd <- 9
  }
  
  rm(groups, transparency)
  
  #   Extract all groups for the legend 
  legendgroups <- c(groupset[1:nrow(groupset), ])
  legendgroups <- substr(legendgroups, 1, max_name_length)
  
  #   Create a data.frame with all groups and their unique input names
  for (i in 1:nrow(groupset)){
    inputgroup <- subset(nametaxa, nametaxa[1:nrow(nametaxa),2] == groupset[i,])
    inputgroup <- inputgroup[ ,c("input")]
    names(inputgroup)[names(inputgroup) == "input"] <- groupset[i,]
    inputgroup <- inputgroup[!duplicated(inputgroup), ]
    if (i == 1) {
      allinputgroup <- inputgroup
    } else {
      allinputgroup <- qpcR:::cbind.na(allinputgroup,inputgroup)
    }
  }
  rm(inputgroup, groupset, nametaxa, i)
  
  max_name_length <- 20 # arbitrarily chop off long species names
  
  
  #   In this scheme, wheat means a gap and red means there are residues present     
  col <- c("white", "red")
  
  hh <- hclust(dist(presence))
  
  #   If the User doesn't want to display all groups make the "Others" grey
  if("Others" %in% colnames(allinputgroup) == TRUE){
    number <- grep("Others", colnames(allinputgroup))
    number <- as.numeric(number)
    color[number] <- "#96969633"
    rm(number)
  }
  
  #   create data frame to show coloring in the Heatmap   
  for (j in 1:length(allinputgroup)){
    for (i in 1:nrow(allinputgroup)){
      ID <- allinputgroup[i,j]
      bc_list  <- which(grepl(ID, hh$labels, ignore.case=T))
      colour <- color[j]
      
      ID  <- data.frame(as.vector(sapply(bc_list, rep, ncol(presence))),
                        rep(seq_len(ncol(presence)), length(bc_list)),
                        rep(colour, length(bc_list) * ncol(presence)),
                        rep(1, length(bc_list) * ncol(presence))
      )
      
      if (i == 1 & j ==1){
        colorCell <- ID
      } else {
        colorCell <- rbind(colorCell, ID)
      }
    }
  }
  rm(i, j, allinputgroup, ID, colour)
  
  #   We have to make the names shorter, so put whatever seems appropriate here.
  #   The names start off like "wp_0543 ...[homo sapiens]     
  shorten_names <- function(nameset) {
    nameset <- sub(".*\\[", "", nameset)
    nameset <- sub("\\]", "", nameset)
    return(substr(nameset, 1, max_name_length))
  }
  #   Names are stored in the hclust object, but the heatmap function
  #   takes them from the data array (presence).
  row.names(presence) <- shorten_names(row.names(presence))
  
  heatmap3::heatmap3(presence, Colv = NA, Rowv = as.dendrogram(hh),
                     showRowDendro = T, scale = "none",
                     cexRow = cexRow,
                     legendfun = function() 
                       heatmap3::showLegend(legend = legendgroups, 
                                            lwd = lwd, cex = cex, 
                                            col = color, ncol = ncol),
                     col = col, highlightCell = colorCell)
  
}

mymain()