library(svDialogs) # used to select options
library(taxadb)
library(RColorBrewer)
library(makeFlow)
library(data.table)
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
  
  #   Filter_name() maps the taxonomic information to the nameset. 
  #   Since NCBI has more unique mappings, we use NCBI
  
  taxnameset <- filter_name(nameset, provider = getOption("taxadb_default_provider","ncbi")) 
  
  #   Ask how many taxonomic ranks the User would like to select (and checked for typos)
  taxranknum <- ""
  while (is.na(taxranknum) == T || taxranknum < 1 || taxranknum > 6 ||round(taxranknum) != taxranknum){
    taxranknum <- dlg_input("How many taxranks do you want to show in the Heatmap? (max. 6)")$res
    taxranknum <- as.numeric(taxranknum)
  }
  
  #   Create data frame with as many rows, as taxonomic ranks the User wanted
  taxranklist <- data.frame(taxranklist = c("a","b","c","d","e","f"),"a","a")
  taxranklist <- data.frame(taxranklist[1:taxranknum, ])
  
  #   Ask which taxonomic ranks should be selected
  for (i in 1:taxranknum) {
    while (taxranklist[i,1] != "kingdom" && taxranklist[i,1] != "phylum" && taxranklist[i,1] != "class" && 
           taxranklist[i,1] != "order" && taxranklist[i,1] != "family" && taxranklist[i,1] != "genus"){
      taxranklist[i,1] <- dlgInput("Which taxonomic rank do you want to select? Please choose:
                  \n'kingdom', 'phylum', 'class', 'order', 'family', 'genus'.")$res
    }
  }
  
  #   Remove all columns from taxanameset except input (species name), Sort (shows duplicates)
  #   and chosen taxrank 
  taxnameset <- taxnameset[ ,c("sort","input",taxranklist[1:nrow(taxranklist),1])]
  
  #   Remove duplicates by using sort as primary key 
  taxnameset <- taxnameset[!duplicated(taxnameset$sort), ]
  
  #   delete "sort" column
  nametaxa <- taxnameset[ ,c(taxranklist[1:nrow(taxranklist),1])]
  taxnameset <- taxnameset[ ,c("input",taxranklist[1:nrow(taxranklist),1])]
  nametaxacopy <- nametaxa
  
  #   Look up how many (unique) groups there are in chosen taxrank for selected Fasta file
  for (i in 1:taxranknum){
    # sort the groups from most seen group to less seen group
    setDT(nametaxacopy)[,freq := .N, by = c(taxranklist[i,1])]
    
    # Sort
    sortedgroups <- nametaxacopy[order(freq, decreasing = T),]
    groupset <- data.frame(unique(sortedgroups[1:nrow(sortedgroups),1]))
    nametaxacopy <- nametaxacopy[ ,-1]
    
    #   Exclude "NA" from our List 
    groupset <- data.frame(groupset[complete.cases(groupset),])
    
    #   How many groups should be displayed from your taxonomic selection
    if (nrow(groupset) > 50) {
      telluser <- paste("How many groups out of", nrow(groupset), 
                        "should be displayed from", taxranklist[i,1], "(max.49 and > 11 not for Colorblind)") 
    } else if (nrow(groupset) > 12) {
      telluser <- paste("How many groups out of", nrow(groupset), 
                        "should be displayed from", taxranklist[i,1], "(> 11 not for Colorblind)") 
    } else {
      telluser <- paste("How many groups out of", nrow(groupset), 
                        "should be displayed from", taxranklist[i,1])
    }
    
    groups <- dlgInput(telluser)$res
    groups <- as.numeric(groups)
    
    #   Check for typos
    while (is.na(groups) == T || groups < 1 || groups > 49 || groups > nrow(groupset) ||round(groups) != groups){
      groups <- dlgInput(telluser)$res
      groups <- as.numeric(groups)
    }
    taxranklist[i,2] <- groups
    
    #   Only if the user does not want to see all groups of his selected taxrank, 
    #   the remaining groups will be renamed as "Others" 
    if (groups != nrow(groupset)){
      groups = groups + 1
      taxranklist[i,2] <- groups
      for (e in groups:nrow(groupset)){
        nametaxa[nametaxa == groupset[e,1]] <- "Others"
      }
      rm(e)
    } 
    #  Create colors for every column
  
    ncolor <- taxranklist[i,2]
    ncolor <- as.numeric(ncolor)
    if (ncolor > 12 && ncolor < 51){
      color <- c(rainbow(ncolor))
    } else if (ncolor > 7 && ncolor < 13){
      color <- c(brewer.pal(ncolor,"Paired"))
    } else if (ncolor > 2){
      color <- c(brewer.pal(ncolor,"Dark2"))
    } else if (ncolor == 2){
      color <- c("#1B9E77", "#D95F02")
    } else {
      color <- "#1B9E77"
    }
  #Here you still have to create a Code to create a list for the colors!
    
  
  }
  
  
  
  #Should the "Others" group be even displayed?
  other <- ""
  while (other != "yes" && other != "no"){
    other <- dlgInput("Should the group 'Others' be displayed? 'yes' or 'no'?")$res
  }
  
  if(other == "no"){
    nametaxa[nametaxa == "Others"] <- NA
  }
  
  #Special formatting of nametaxa for the pheatmap. If you are interested in a different 
  #hierarchy “class” should be replaced with the category you are interested in 
  #(see above "kingdom","phylum"...). The object class will be the title in the diagram. 
  #Since identical names are not allowed in row.names(), we use make.unique() to number identical names.
  taxrank <- nametaxa
  annotation_row = data.frame(taxrank)
  row.names(annotation_row) = make.unique(nameset)
  
  #   In this scheme, wheat means a gap and red means there are residues present     
  col <- c("wheat", "red")
  
  hh <- hclust(dist(presence, method="binary"), method = "complete")
  
  rownames(presence) <- paste(nameset)
  rownames(taxnameset) <- colnames(presence)
  
  colspace <- round(ncol(presence) / 20)
  labels_col <- seq (from = 1, to = ncol(presence), by = colspace)
  if (length(nameset) < 21){
    show_rownames <- T
  } else {
    show_rownames <- F
  }
  pheatmap::pheatmap(presence,
                     cluster_cols = F, cluster_rows = hh,
                     color = col,
                     legend = F,
                     annotation_row = annotation_row,
                     show_colnames = T,
                     show_rownames = show_rownames,
                     labels_col = labels_col
                     
  )
}

mymain()