#### Initialize ----

zap()
setwd("/Users/elb4003/Work/PostDocWCM/Projects/MetaboPanCan/MetaboRNA/")

library(readxl)
library(tidyverse)
library(magrittr)
library(purrr)
library(ggrepel)
library(ggpubr)
library(plotly)
library(htmlwidgets)
library(ggpubr)
library(ggforce)
library(grid)
library(moonBook)
library(webr)
library(scales)

library(readxl)
library(pheatmap)
library(colorRamps)
library(RColorBrewer)
library(extrafont)
library(patchwork)
library(ggbreak)
library(xlsx)

# check if directory to save results exists, otherwise create
if(!dir.exists(paste0(Sys.Date(),sep=""))) {
  dir.create(paste0(Sys.Date(),sep=""))
}

#### Load Data ----

# load master mapping file
mapping_file <- read.csv2(file = "Files2Share/MasterMapping_MetImmune_03_16_2022.csv", sep=",") 

load("Files2Share/PreprocessedData.Rdata")

load("Files2Share/Annotations_Metabolon_Metabolites.Rdata")
load("Files2Share/Annotations_Pathways.Rdata")

load("Files2Share/Results_Concordance_MetaAnalysis.Rdata")
load("Files2Share/Results_FisherPathwayEnrichment.Rdata")

mapping_color <- data.frame(OldName= c("BRCA_Terunuma","BRCA_Tang","COAD","CPTAC_GBM","DLBCL",
                                       "HCC","LiCa1","LiCa2","OV","PDAC",
                                       "PRAD","RC12","RC18","RC18_2","RC20"),
                            NewName = c("BRCA1","BRCA2","COAD","GBM","DLBCL",
                                        "HurthleCC","HCC","ICC","OV","PDAC",
                                        "PRAD","ccRCC1","ccRCC2","ccRCC3","ccRCC4")) %>%
  dplyr::arrange(NewName) %>%
  dplyr::mutate(Color = c("#e5b6a5","#c3c1a4","#e5e3c0","#bbe1b9","#88beab",
                        "#97b1ab","#abc5bf","#abe7e1","#c6e6e7","#8fcde1",
                        "#88aee1","#a6b7d3","#d2d5ed","#c8b7de","#e6bbcc"))

#### Set plotting parameters  ----

# x <- "T"
# significance threshold on the FDR adjusted p-values
pcut <- 0.01
# cutoff on the number of datasets a metabolite must be present in
n_data_cut <- 7

# rename variable
res_conc <- conc
remove(conc)

#### Write Significant GMIs to file ----

write.table(x=res_conc %>%
              dplyr::filter(n_dataset>n_data_cut) %>%
              dplyr::filter(padj<pcut),
            file=sprintf("%s/TableS1.txt", Sys.Date()),
            append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)

#### Define Functions ----

# generate piedonut ----
# modified from the webr PieDonut function
PieDonut_edited <- function (data, mapping, start = getOption("PieDonut.start", 
                                                              0), addPieLabel = TRUE, addDonutLabel = TRUE, showRatioDonut = TRUE, 
                             showRatioPie = TRUE, ratioByGroup = TRUE, showRatioThreshold = getOption("PieDonut.showRatioThreshold", 
                                                                                                      0.02), labelposition = getOption("PieDonut.labelposition", 
                                                                                                                                       2), labelpositionThreshold = 0.1, r0 = getOption("PieDonut.r0", 
                                                                                                                                                                                        0.3), r1 = getOption("PieDonut.r1", 1), r2 = getOption("PieDonut.r2", 
                                                                                                                                                                                                                                               1.2), explode = NULL, selected = NULL, explodePos = 0.1, 
                             color = "white", pieAlpha = 0.8, donutAlpha = 1, maxx = NULL, 
                             showPieName = TRUE, showDonutName = FALSE, title = NULL, 
                             pieLabelSize = 4, donutLabelSize = 3, titlesize = 5, explodePie = TRUE, 
                             explodeDonut = FALSE, use.label = TRUE, use.labels = TRUE, 
                             family = getOption("PieDonut.family", ""), palette=NULL) 
{
  (cols = colnames(data))
  if (use.labels) 
    data = moonBook::addLabelDf(data, mapping)
  count <- NULL
  if ("count" %in% names(mapping)) 
    count <- getMapping(mapping, "count")
  count
  pies <- donuts <- NULL
  (pies = getMapping(mapping, "pies"))
  if (is.null(pies)) 
    (pies = getMapping(mapping, "pie"))
  if (is.null(pies)) 
    (pies = getMapping(mapping, "x"))
  (donuts = getMapping(mapping, "donuts"))
  if (is.null(donuts)) 
    (donuts = getMapping(mapping, "donut"))
  if (is.null(donuts)) 
    (donuts = getMapping(mapping, "y"))
  if (!is.null(count)) {
    df <- data %>% group_by(.data[[pies]]) %>% dplyr::summarize(Freq = sum(.data[[count]]))
    df
  }
  else {
    df = data.frame(table(data[[pies]]))
  }
  colnames(df)[1] = pies
  df$end = cumsum(df$Freq)
  df$start = dplyr::lag(df$end)
  df$start[1] = 0
  total = sum(df$Freq)
  df$start1 = df$start * 2 * pi/total
  df$end1 = df$end * 2 * pi/total
  df$start1 = df$start1 + start
  df$end1 = df$end1 + start
  df$focus = 0
  if (explodePie) 
    df$focus[explode] = explodePos
  df$mid = (df$start1 + df$end1)/2
  df$x = ifelse(df$focus == 0, 0, df$focus * sin(df$mid))
  df$y = ifelse(df$focus == 0, 0, df$focus * cos(df$mid))
  df$label = df[[pies]]
  df$ratio = df$Freq/sum(df$Freq)
  if (showRatioPie) {
    df$label = ifelse(df$ratio >= showRatioThreshold, 
                      paste0(df$label, " ", df$Freq), 
                      as.character(df$label))
  }
  df$labelx = (r0 + r1)/2 * sin(df$mid) + df$x
  df$labely = (r0 + r1)/2 * cos(df$mid) + df$y
  if (!is.factor(df[[pies]])) 
    df[[pies]] <- factor(df[[pies]])
  df
  if(!is.null(palette)){
    stopifnot("Number of colors not correct"=length(palette)==nrow(df))
    mainCol = palette
  }else{
    mainCol = gg_color_hue(nrow(df))
  }
  df$radius = r1
  df$radius[df$focus != 0] = df$radius[df$focus != 0] + df$focus[df$focus != 
                                                                   0]
  df$hjust = ifelse((df$mid%%(2 * pi)) > pi, 1, 0)
  df$vjust = ifelse(((df$mid%%(2 * pi)) < (pi/2)) | (df$mid%%(2 * 
                                                                pi) > (pi * 3/2)), 0, 1)
  df$segx = df$radius * sin(df$mid)
  df$segy = df$radius * cos(df$mid)
  df$segxend = (df$radius + 0.05) * sin(df$mid)
  df$segyend = (df$radius + 0.05) * cos(df$mid)
  df
  if (!is.null(donuts)) {
    subColor = webr::makeSubColor(mainCol, no = length(unique(data[[donuts]])))
    subColor
    data
    if (!is.null(count)) {
      df3 <- as.data.frame(data[c(donuts, pies, count)])
      colnames(df3) = c("donut", "pie", "Freq")
      df3
      df3 <- eval(parse(text = "complete(df3,donut,pie)"))
      df3$Freq[is.na(df3$Freq)] = 0
      if (!is.factor(df3[[1]])) 
        df3[[1]] = factor(df3[[1]])
      if (!is.factor(df3[[2]])) 
        df3[[2]] = factor(df3[[2]])
      df3 <- df3 %>% arrange(.data$pie, .data$donut)
      a <- df3 %>% spread(.data$pie, value = .data$Freq)
      a = as.data.frame(a)
      a
      rownames(a) = a[[1]]
      a = a[-1]
      a
      colnames(df3)[1:2] = c(donuts, pies)
    }
    else {
      df3 = data.frame(table(data[[donuts]], data[[pies]]), 
                       stringsAsFactors = FALSE)
      colnames(df3)[1:2] = c(donuts, pies)
      a = table(data[[donuts]], data[[pies]])
      a
    }
    a
    df3
    df3$group = rep(colSums(a), each = nrow(a))
    df3$pie = rep(1:ncol(a), each = nrow(a))
    total = sum(df3$Freq)
    total
    df3$ratio1 = df3$Freq/total
    df3
    if (ratioByGroup) {
      df3$ratio = scales::percent(df3$Freq/df3$group)
    }
    else {
      df3$ratio <- scales::percent(df3$ratio1)
    }
    df3$end = cumsum(df3$Freq)
    df3
    df3$start = dplyr::lag(df3$end)
    df3$start[1] = 0
    df3$start1 = df3$start * 2 * pi/total
    df3$end1 = df3$end * 2 * pi/total
    df3$start1 = df3$start1 + start
    df3$end1 = df3$end1 + start
    df3$mid = (df3$start1 + df3$end1)/2
    df3$focus = 0
    if (!is.null(selected)) {
      df3$focus[selected] = explodePos
    }
    else if (!is.null(explode)) {
      selected = c()
      for (i in 1:length(explode)) {
        start = 1 + nrow(a) * (explode[i] - 1)
        selected = c(selected, start:(start + nrow(a) - 
                                        1))
      }
      selected
      df3$focus[selected] = explodePos
    }
    df3
    df3$x = 0
    df3$y = 0
    df
    if (!is.null(explode)) {
      explode
      for (i in 1:length(explode)) {
        xpos = df$focus[explode[i]] * sin(df$mid[explode[i]])
        ypos = df$focus[explode[i]] * cos(df$mid[explode[i]])
        df3$x[df3$pie == explode[i]] = xpos
        df3$y[df3$pie == explode[i]] = ypos
      }
    }
    df3$no = 1:nrow(df3)
    df3$label = df3[[donuts]]
    if (showRatioDonut) {
      if (max(nchar(levels(df3$label))) <= 2) 
        df3$label = paste0(df3$label, " ", df3$Freq)
      else df3$label = paste0(df3$label, " ", df3$Freq)
    }
    df3$label[df3$ratio1 == 0] = ""
    df3$label[df3$ratio1 < showRatioThreshold] = ""
    df3$hjust = ifelse((df3$mid%%(2 * pi)) > pi, 1, 0)
    df3$vjust = ifelse(((df3$mid%%(2 * pi)) < (pi/2)) | 
                         (df3$mid%%(2 * pi) > (pi * 3/2)), 0, 1)
    df3$no = factor(df3$no)
    df3
    labelposition
    if (labelposition > 0) {
      df3$radius = r2
      if (explodeDonut) 
        df3$radius[df3$focus != 0] = df3$radius[df3$focus != 
                                                  0] + df3$focus[df3$focus != 0]
      df3$segx = df3$radius * sin(df3$mid) + df3$x
      df3$segy = df3$radius * cos(df3$mid) + df3$y
      df3$segxend = (df3$radius + 0.05) * sin(df3$mid) + 
        df3$x
      df3$segyend = (df3$radius + 0.05) * cos(df3$mid) + 
        df3$y
      if (labelposition == 2) 
        df3$radius = (r1 + r2)/2
      df3$labelx = (df3$radius) * sin(df3$mid) + df3$x
      df3$labely = (df3$radius) * cos(df3$mid) + df3$y
    }
    else {
      df3$radius = (r1 + r2)/2
      if (explodeDonut) 
        df3$radius[df3$focus != 0] = df3$radius[df3$focus != 
                                                  0] + df3$focus[df3$focus != 0]
      df3$labelx = df3$radius * sin(df3$mid) + df3$x
      df3$labely = df3$radius * cos(df3$mid) + df3$y
    }
    df3$segx[df3$ratio1 == 0] = 0
    df3$segxend[df3$ratio1 == 0] = 0
    df3$segy[df3$ratio1 == 0] = 0
    df3$segyend[df3$ratio1 == 0] = 0
    if (labelposition == 0) {
      df3$segx[df3$ratio1 < showRatioThreshold] = 0
      df3$segxend[df3$ratio1 < showRatioThreshold] = 0
      df3$segy[df3$ratio1 < showRatioThreshold] = 0
      df3$segyend[df3$ratio1 < showRatioThreshold] = 0
    }
    df3
    del = which(df3$Freq == 0)
    del
    if (length(del) > 0) 
      subColor <- subColor[-del]
    subColor
  }
  p <- ggplot() + theme_no_axes() + coord_fixed()
  if (is.null(maxx)) {
    r3 = r2 + 0.3
  }
  else {
    r3 = maxx
  }
  p1 <- p + geom_arc_bar(aes_string(x0 = "x", y0 = "y", r0 = as.character(r0), 
                                    r = as.character(r1), start = "start1", end = "end1", 
                                    fill = pies), alpha = pieAlpha, color = color, data = df) + 
    transparent() + scale_fill_manual(values = mainCol) + 
    xlim(r3 * c(-1, 1)) + ylim(r3 * c(-1, 1)) + guides(fill = FALSE)
  if ((labelposition == 1) & (is.null(donuts))) {
    p1 <- p1 + geom_segment(aes_string(x = "segx", y = "segy", 
                                       xend = "segxend", yend = "segyend"), data = df) + 
      geom_text(aes_string(x = "segxend", y = "segyend", 
                           label = "label", hjust = "hjust", vjust = "vjust"), 
                size = pieLabelSize, data = df, family = family)
  }
  else if ((labelposition == 2) & (is.null(donuts))) {
    p1 <- p1 + geom_segment(aes_string(x = "segx", y = "segy", 
                                       xend = "segxend", yend = "segyend"), data = df[df$ratio < 
                                                                                        labelpositionThreshold, ]) + 
      geom_text(aes_string(x = "segxend", y = "segyend", label = "label", hjust = "hjust", vjust = "vjust"), size = pieLabelSize, data = df[df$ratio < labelpositionThreshold, ], family = family) + 
      geom_text(aes_string(x = "labelx", y = "labely", label = "label"), size = pieLabelSize, data = df[df$ratio >= labelpositionThreshold, ], family = family)
  }
  else {
    p1 <- p1 + geom_text(aes_string(x = "labelx", y = "labely", 
                                    label = "label"), size = pieLabelSize, data = df, 
                         family = family)
  }
  if (showPieName) 
    p1 <- p1 + annotate("text", x = 0, y = 0, label = pies, 
                        size = titlesize, family = family)
  p1 <- p1 + theme(text=element_text(size=7))
  if (!is.null(donuts)) {
    if (explodeDonut) {
      p3 <- p + geom_arc_bar(aes_string(x0 = "x", y0 = "y", 
                                        r0 = as.character(r1), r = as.character(r2), 
                                        start = "start1", end = "end1", fill = "no", 
                                        explode = "focus"), alpha = donutAlpha, color = color, 
                             data = df3)
    }
    else {
      p3 <- p + geom_arc_bar(aes_string(x0 = "x", y0 = "y", 
                                        r0 = as.character(r1), r = as.character(r2), 
                                        start = "start1", end = "end1", fill = "no"), 
                             alpha = donutAlpha, color = color, data = df3)
    }
    p3 <- p3 + transparent() + scale_fill_manual(values = subColor) + 
      xlim(r3 * c(-1, 1)) + ylim(r3 * c(-1, 1)) + guides(fill = FALSE)
    p3
    if (labelposition == 1) {
      p3 <- p3 + geom_segment(aes_string(x = "segx", y = "segy", 
                                         xend = "segxend", yend = "segyend"), data = df3) + 
        geom_text(aes_string(x = "segxend", y = "segyend", 
                             label = "label", hjust = "hjust", vjust = "vjust"), 
                  size = donutLabelSize, data = df3, family = family)
    }
    else if (labelposition == 0) {
      p3 <- p3 + geom_text(aes_string(x = "labelx", y = "labely", 
                                      label = "label"), size = donutLabelSize, data = df3, 
                           family = family)
    }
    else {
      p3 <- p3 + geom_segment(aes_string(x = "segx", y = "segy", 
                                         xend = "segxend", yend = "segyend"), data = df3[df3$ratio1 < 
                                                                                           labelpositionThreshold, ]) + 
        geom_text(aes_string(x = "segxend", y = "segyend", label = "label", hjust = "hjust", vjust = "vjust"), size = donutLabelSize, data = df3[df3$ratio1 < labelpositionThreshold, ], family = family) + 
        geom_text(aes_string(x = "labelx", y = "labely", 
                             label = "label"), size = donutLabelSize, data = df3[df3$ratio1 >= labelpositionThreshold, ], family = family)
    }
    if (!is.null(title)) 
      p3 <- p3 + annotate("text", x = 0, y = r3, label = title, 
                          size = titlesize, family = family)
    else if (showDonutName) 
      p3 <- p3 + annotate("text", x = (-1) * r3, y = r3, 
                          label = donuts, hjust = 0, size = titlesize, 
                          family = family)
    p3 <- p3 + theme(text=element_text(size=7))
    grid.newpage()
    print(p1, vp = viewport(height = 1, width = 1))
    print(p3, vp = viewport(height = 1, width = 1))
  }
  else {
    p1
  }
}


# compute all intersections based on a named list ----
multiIntersect <- function (lst) {
  
  # assemble list of all items
  items = unique(unlist(lst))
  # assign each sublist  
  assignmat = sapply(lst, function(x){!is.na(match(items,x))});
  # set names
  rownames(assignmat) = items
  colnames(assignmat) = names(lst)
  
  as.data.frame(assignmat)
}

# hue color palette ----
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# function to move files ----
# from https://stackoverflow.com/questions/10266963/moving-files-between-folders
my.file.rename <- function(from, to) {
  todir <- dirname(to)
  if (!isTRUE(file.info(todir)$isdir)) dir.create(todir, recursive=TRUE)
  file.rename(from = from,  to = to)
}

# copying functions that should be in ggpubr but apparently are not exported ----
# check if object is ggtexttable
is_ggtexttable <- function(tab){
  !is.null(attr(tab, "ggtexttableGrob"))
}
# Transform a table grob in ggtexttable like object
as_ggtexttable <- function(tabgrob){
  res <- as_ggplot(tabgrob)
  attr(res, "ggtexttableGrob") <- tabgrob
  res
}
# Return the same class as the input data,
# which can be either ggtextable or a gridExtra::tableGrob
tab_return_same_class_as_input <- function(tabgrob, input){
  if(is_ggtexttable(input)){
    return(as_ggtexttable(tabgrob))
  }
  else if(is_tablegrob(input)){
    return(tabgrob)
  }
  tabgrob
}
# Extract tableGrob from ggtexttable()
get_tablegrob <- function(tab){
  if(is_ggtexttable(tab)){
    tabgrob <- attr(tab, "ggtexttableGrob")
  }
  else if(is_tablegrob(tab)){
    tabgrob <- tab
  }
  else{
    stop("tab should be an object from either ggpubr::ggtexttable() or gridExtra::tableGrob().")
  }
  tabgrob
}
# add title to table
tab_add_title <- function(tab, text,  face = NULL, size = NULL, color = NULL,
                          family = NULL, padding = unit(1.5,"line"),
                          just = "left", hjust = NULL, vjust = NULL){
  library("gtable")
  tabgrob <- get_tablegrob(tab)
  text <- grid::textGrob(
    text, x = 0.02, just = just, hjust = hjust, vjust = vjust,
    gp = grid::gpar(fontsize = size, fontface = face, fontfamily = family, col = color)
  )
  # Add row at the top
  tabgrob <- gtable::gtable_add_rows(
    tabgrob, heights = grid::grobHeight(text) + padding,
    pos = 0
  )
  tabgrob <- gtable::gtable_add_grob(
    tabgrob, list(text),
    t = 1, b = 1, l = 1, r = ncol(tabgrob)
  )
  tab_return_same_class_as_input(tabgrob, input = tab)
}
# add footnote to table
tab_add_footnote <- function(tab, text,  face = NULL, size = NULL, color = NULL,
                             family = NULL, padding = unit(1.5,"line"),
                             just = "right", hjust = NULL, vjust = NULL){
  library("gtable")
  tabgrob <- get_tablegrob(tab)
  text <- grid::textGrob(
    text, x = 0.95, just = just, hjust = hjust, vjust = vjust,
    gp = grid::gpar(fontsize = size, fontface = face, fontfamily = family, col = color)
  )
  # Add row at the bottom
  tabgrob <- gtable::gtable_add_rows(
    tabgrob, heights = grid::grobHeight(text) + padding,
    pos = -1
  )
  tabgrob <- gtable::gtable_add_grob(
    tabgrob, list(text),
    t = nrow(tabgrob), b = nrow(tabgrob), l = 1, r = ncol(tabgrob)
  )
  tab_return_same_class_as_input(tabgrob, input = tab)
}

# moveme function taken from stackoverflow ----
# https://stackoverflow.com/questions/3369959/moving-columns-within-a-data-frame-without-retyping
moveme <- function (invec, movecommand) {
  movecommand <- lapply(strsplit(strsplit(movecommand, ";")[[1]], 
                                 ",|\\s+"), function(x) x[x != ""])
  movelist <- lapply(movecommand, function(x) {
    Where <- x[which(x %in% c("before", "after", "first", 
                              "last")):length(x)]
    ToMove <- setdiff(x, Where)
    list(ToMove, Where)
  })
  myVec <- invec
  for (i in seq_along(movelist)) {
    temp <- setdiff(myVec, movelist[[i]][[1]])
    A <- movelist[[i]][[2]][1]
    if (A %in% c("before", "after")) {
      ba <- movelist[[i]][[2]][2]
      if (A == "before") {
        after <- match(ba, temp) - 1
      }
      else if (A == "after") {
        after <- match(ba, temp)
      }
    }
    else if (A == "first") {
      after <- 0
    }
    else if (A == "last") {
      after <- length(myVec)
    }
    myVec <- append(temp, values = movelist[[i]][[1]], after = after)
  }
  myVec
}

#### Create ausiliary variables ----

# create a mapping dataframe
df <- met_all %>% sapply(ncol) %>% as.data.frame %>% 
  dplyr::rename(n=".") %>% 
  tibble::rownames_to_column("dataset") %>%
  dplyr::arrange(dataset) %>%
  dplyr::mutate(corcol=paste(dataset,"_cor",sep=""))
df %<>%
  dplyr::mutate(corcol=ifelse(grepl(pattern = "_Tumor", x=df$dataset, fixed = T),
                              gsub(pattern = "_Tumor",replacement = "",x = corcol),corcol)) %>%
  # change one entry manually to match result dataframe
  dplyr::mutate(corcol=ifelse(corcol=="RC18_2_Normal_cor","RC18_Normal_2_cor",corcol)) %>%
  dplyr::mutate(tissue=ifelse(grepl(pattern = "_Tumor", x=df$dataset, fixed = T),"Tumor","Normal"))

# merge renaming and coloring variables
df %<>% 
  dplyr::mutate(cohort=sub("_[^_]+$", "", dataset)) %>% 
  dplyr::left_join(mapping_color, by=c("cohort"="OldName"))

#### Figure 1 - Dataset Overview ----

# Sample Size Overview
pdf(sprintf("%s/Figure1A.pdf", Sys.Date()), width=5, height=5, onefile = F)
print(
  sapply(met_all, ncol) %>% unlist %>%
  as.data.frame %>% 
  tibble::rownames_to_column("key") %>%
  dplyr::rename(value=".") %>% 
  dplyr::mutate(tissue=ifelse(grepl(pattern = "Tumor", x = key),"Tumor","Normal")) %>% 
  dplyr::mutate(cohort=ifelse(tissue=="Tumor",
                              str_remove(string=key,pattern = "_Tumor"),
                              str_remove(string=key,pattern = "_Normal"))) %>%
  dplyr::left_join(mapping_color, by=c("cohort"="OldName")) %>% 
  PieDonut_edited(aes(NewName,tissue,count=value),
           r0=0.3,showPieName=FALSE,showRatioThreshold = F,
           explodeDonut = F, selected=c(1:15),
           explodePie = F, explode = c(1:15), explodePos = 0.3,
           labelposition=1, 
           palette = mapping_color %>% dplyr::pull(Color), pieAlpha = 1, donutAlpha = 1),
  vp = viewport(height = 1, width = 1))
dev.off()

# Transcripts overlap across datasets
p_rna <- lapply(mapping_color$OldName %>% {names(.)=.;.}, function(x){
  # get list of metabolites per cohort
  n_met=rna_all %>% lapply(rownames) %>% multiIntersect() %>% 
    dplyr::select(any_of(c(paste0(x,"_Tumor"),paste0(x,"_Normal")))) %>% 
    dplyr::mutate(n=rowSums(.)) %>% 
    dplyr::mutate(count=ifelse(n>0,1,0)) %>% 
    dplyr::filter(count==1) %>% 
    rownames
}) %>% 
  multiIntersect %>% 
  dplyr::mutate(ntot=rowSums(.)) %>% 
  dplyr::select(ntot) %>% 
  table(dnn="ntot") %>% as.data.frame %>% 
  dplyr::mutate(ntot=as.numeric(ntot)) %>%
  dplyr::arrange(-ntot) %>% 
  dplyr::mutate(n_cumul = cumsum(Freq)) %>%
  dplyr::mutate(omics="Transcriptomics") %>% 
  ggplot(aes(x=factor(ntot), y=n_cumul, group=omics)) +
  geom_line(color="grey60") +
  geom_point() +
  xlab("Number of datasets") +
  ylab("Number of common Transcripts") +
  ylim(c(0,45000)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text=element_text(size=7), axis.text=element_text(size=7), axis.title=element_text(size=7),legend.position="none")
# save to file
pdf(sprintf("%s/FigureS4A.pdf", Sys.Date()), width=3.5, height=3.5, onefile = F)
print(p_rna)
dev.off()

# Metabolites overlap across datasets
p_met <- lapply(mapping_color$OldName %>% {names(.)=.;.}, function(x){
  # get list of metabolites per cohort
  n_met=met_all %>% lapply(rownames) %>% multiIntersect() %>% 
    dplyr::select(any_of(c(paste0(x,"_Tumor"),paste0(x,"_Normal")))) %>% 
    dplyr::mutate(n=rowSums(.)) %>% 
    dplyr::mutate(count=ifelse(n>0,1,0)) %>% 
    dplyr::filter(count==1) %>% 
    rownames
}) %>% 
  multiIntersect %>%
  dplyr::mutate(ntot=rowSums(.)) %>%
  dplyr::select(ntot) %>% 
  table(dnn="ntot") %>% as.data.frame %>% 
  dplyr::mutate(ntot=as.numeric(ntot)) %>%
  dplyr::arrange(-ntot) %>% 
  dplyr::mutate(n_cumul = cumsum(Freq)) %>% 
  dplyr::mutate(omics="Metabolomics") %>%
  ggplot(aes(x=factor(ntot), y=n_cumul, group=omics)) +
  geom_line(color="grey60") +
  geom_point() +
  xlab("Number of datasets") +
  ylab("Number of common Metabolites") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text=element_text(size=7), axis.text=element_text(size=7), axis.title=element_text(size=7),legend.position="none")
# save to file
pdf(sprintf("%s/FigureS4B.pdf", Sys.Date()), width=3.5, height=3.5, onefile = F)
print(p_met)
dev.off()

#### Pathway Distance Analysis  ----

# short distance: GEM dist <=2
dt <- res_conc %>%
  dplyr::filter(n_dataset>n_data_cut) %>%
  dplyr::mutate(is_significant=ifelse(padj<pcut,1,0) %>% factor(levels = c(1,0))) %>%
  dplyr::mutate(is_close=ifelse(dist_GEM_min <=2 & !is.na(dist_GEM_min),1,0) %>% factor(levels = c(1,0))) 

contab <- dt %>%
  dplyr::select(is_significant, is_close) %>%
  table
fisher.test(contab, alternative="greater")

m_breaks <- as.numeric(c(seq(0,0.05,0.005), seq(0.1,1,0.1)))
p_proximal_barplot <- contab %>% as.data.frame %>%
  dplyr::group_by(is_significant) %>%
  dplyr::mutate(perc=Freq/sum(Freq)) %>%
  dplyr::mutate(is_close=factor(is_close,levels=c(0,1))) %>%
  ggplot(aes(x=is_significant, y=perc, fill=is_close, label=sprintf("%.2f%%",perc*100))) +
  geom_bar(position="stack", stat="identity") +
  scale_y_continuous(breaks = m_breaks) +
  geom_text(size = 3, position = position_stack(vjust = 1)) +
  scale_fill_manual(values=c("violetred4","salmon")) +
  theme_minimal() +
  theme(text=element_text(size=7), axis.text=element_text(size=7), axis.title=element_text(size=7), strip.text.y = element_text(angle=0)) +
  theme(legend.position="bottom") +
  xlab("Significant") +
  ylab("Percentage") +
  scale_y_break(c(0.025, 0.95))

pdf(sprintf("%s/Figure3A.pdf", Sys.Date()), height = 4, width = 2, onefile=F)
print(p_proximal_barplot)
dev.off()

#### Piecharts ----

# piecharts
p_pie_proximal <- contab %>% as.data.frame %>%
  dplyr::filter(is_close==1) %>%
  dplyr::mutate(label=ifelse(is_significant==1,"significant","non significant")) %>%
  dplyr::mutate(prop = Freq / sum(.$Freq) *100) %>% 
  dplyr::mutate(text=sprintf("%.2f%%",prop)) %>% 
  ggplot(aes(x="", y=Freq, fill=label)) +
  geom_bar(stat="identity", width=1, color="white") +
  geom_text(aes(y = prop*prop, x=1, label = text)) +
  coord_polar("y", start=0) +
  scale_fill_manual(values = c("gray80","firebrick")) +
  theme_void() +
  theme(text=element_text(size=7)) +
  ggtitle("Proximal gene-metabolite pairs")
# save to file
pdf(sprintf("%s/Piechart_proximal.pdf", Sys.Date()), height = 2, width = 3, onefile=F)
print(p_pie_proximal)
dev.off()

# piecharts
p_pie_significant <- contab %>% as.data.frame %>%
  dplyr::filter(is_significant==1) %>%
  dplyr::mutate(label=ifelse(is_close==1,"proximal","non proximal")) %>%
  dplyr::mutate(prop = Freq / sum(.$Freq) *100) %>% 
  dplyr::mutate(text=sprintf("%.2f%%",prop)) %>% 
  ggplot(aes(x="", y=Freq, fill=label)) +
  geom_bar(stat="identity", width=1, color="white") +
  geom_text(aes(y = prop*prop/2, x=1, label = text)) +
  coord_polar("y", start=0) +
  scale_fill_manual(values = c("gray80","firebrick")) +
  theme_void() +
  theme(text=element_text(size=7)) +
  ggtitle("Significant GMIs")
# save to file
pdf(sprintf("%s/Piechart_significant.pdf", Sys.Date()), height = 2, width = 3, onefile=F)
print(p_pie_significant)
dev.off()

#### Hubs missingness ----

# compute number of missing values for all metabolites in each dataset
miss_met <- met_all[grep(pattern = "Tumor",x=names(met_all), value = T)] %>%
  lapply(function(x){
    x %>% apply(1, function(m){
      sum(m==min(m))
    }) %>% t %>% as.data.frame 
  }) %>% {Reduce(f = full_join, x = .)} %>% 
  as.data.frame
rownames(miss_met) <- names(met_all[grep(pattern = "Tumor",x=names(met_all), value = T)])
# compute total number of missing values
tot_miss <- miss_met %>% colSums(na.rm = T) 

# compute number of measurements for each metabolite
tot_n <- met_all[grep(pattern = "Tumor",x=names(met_all), value = T)] %>%
  lapply(function(x){
    x %>% apply(1, function(m){
      length(m)
    }) %>% t %>% as.data.frame 
  }) %>% {Reduce(f = full_join, x = .)} %>% 
  as.data.frame %>%
  colSums(na.rm = T)

# select metabolites that are used for analysis
mm <- res_conc %>%
  dplyr::filter(n_dataset>n_data_cut) %>% 
  dplyr::pull(metabolite) %>% unique

dt_miss <- miss_met %>%
  t %>% as.data.frame %>%
  dplyr::mutate(missing_n=rowSums(., na.rm = T)) %>%
  tibble::rownames_to_column("metabolite") %>%
  dplyr::left_join(tot_n %>% as.data.frame %>% dplyr::rename(tot_n=".") %>% tibble::rownames_to_column("metabolite"),
                   by="metabolite") %>%
  dplyr::mutate(missing_frac=missing_n/tot_n) %>%
  dplyr::mutate(hub=ifelse(metabolite %in% topmet, "yes","no")) %>%
  dplyr::mutate(meta_analysis=ifelse(metabolite %in% mm, "yes","no")) 

# dt_miss %>%
#   dplyr::filter(meta_analysis=="yes") %>%
#   ggplot(aes(x=hub,y=missing_frac)) +
#   geom_boxplot(outlier.shape=NA) +
#   geom_point(position=position_jitter(width = 0.2)) +
#   theme_bw() 

wilcox.test(x=dt_miss %>% dplyr::filter(meta_analysis=="yes") %>% dplyr::filter(hub=="yes") %>% dplyr::pull(missing_frac),
            y=dt_miss %>% dplyr::filter(meta_analysis=="yes") %>% dplyr::filter(hub=="no") %>% dplyr::pull(missing_frac))$p.value

# save missing data to file
write.table(x=dt_miss %>% 
              dplyr::filter(meta_analysis=="yes"),
            file=sprintf("%s/Table_metabolite_missing.txt", Sys.Date()),
            append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)

#### Volcano Plot ----

# volcano plot
p_volc <- res_conc %>%
  dplyr::filter(n_dataset>n_data_cut) %>%
  # dplyr::left_join(dt_ccle %>% dplyr::select(metabolite,gene,color), by=c("metabolite","gene")) %>%
  dplyr::mutate(repl=case_when(padj>=pcut ~ "not_significant",
                               padj<pcut & !(dist_GEM_min %in% c(1,2)) ~ "significant",
                               padj<pcut & dist_GEM_min %in% c(1,2) ~ "significant_proximal")) %>%
  dplyr::mutate(label=ifelse(padj<1E-10 & dist_GEM_min<=2,
                             sprintf("%s - %s", gene, metabolite),NA)) %>% 
  # order points for plotting
  dplyr::mutate(var=case_when(repl=="not_significant"~1,repl=="significant"~2,repl=="significant_proximal"~3)) %>%
  dplyr::arrange(var) %>% 
  ggplot(aes(x=2*concordance-1, y=-log10(padj),
             label=label)) +
  geom_point(aes(color=repl), alpha=0.7, size=0.3) +
  # geom_text_repel(size=2) +
  xlab("2*concordance-1") +
  ylab("-log10(adjusted p-value)") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 36)) +
  scale_color_manual(values=c("significant_proximal"="gray30","significant"="gray30","not_significant"="gray80")) +
  geom_hline(yintercept = -log10(pcut), linetype="dashed", color="gray30") +
  geom_vline(xintercept = 0, linetype="dashed", color="gray30") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.y = element_text(color = "black"), axis.text.x = element_text(color = "black"),
        text=element_text(size=7), axis.text=element_text(size=7), axis.title=element_text(size=7),legend.position="none")

# save as png
ggsave(p_volc,filename = sprintf("%s/Figure2A.png", Sys.Date()), width=2, height=2, units = "in")

# concordance distribution
p_conc_distr <- res_conc %>%
  dplyr::filter(n_dataset>7) %>%
  dplyr::filter(padj<0.01) %>% 
  tidyr::pivot_longer(cols = grep(colnames(res_conc), pattern="_Tumor", value = T), names_to = "dataset", values_to = "conc") %>% 
  dplyr::filter(gene %in% c("IDO1","GDA","CD38") & metabolite %in% c("kynurenine","guanine","nicotinamide ribonucleotide (NMN)")) %>% 
  dplyr::mutate(label=sprintf("%s_%s",gene,metabolite)) %>% 
  dplyr::mutate(dataset=str_remove(dataset, pattern = "_Tumor") %>% factor(levels=mapping_color$OldName)) %>%
  ggplot(aes(x=label, y=2*conc-1)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(color=dataset), position = position_jitter(width = 0.2)) +
  scale_color_manual(values = mapping_color$Color) +
  ylab("Concordance") +
  xlab("") +
  ylim(c(-1,1)) +
  theme_bw() +
  theme(text=element_text(size=7), axis.text=element_text(size=7), 
        axis.title=element_text(size=7), 
        legend.position = "none")

# save as pdf
pdf(sprintf("%s/FigureS1.pdf", Sys.Date()), width=2.5, height=2.5)
print(p_conc_distr)
dev.off()

#### Gene Prioritization ----

# # metabolites with at least one proximal GMI
# proximal_met <- res_conc %>%
#   dplyr::filter(n_dataset>n_data_cut) %>%
#   dplyr::filter(padj<pcut) %>%
#   dplyr::filter(dist_GEM_min <= 2) %>% 
#   dplyr::pull(metabolite) %>% unique

# select metabolites whose strongest GMI is proximal
proximal_met <- res_conc %>%
  # filter significant
  dplyr::filter(n_dataset>n_data_cut) %>%
  dplyr::filter(padj<pcut) %>%
  # group by metabolite
  dplyr::group_by(metabolite) %>%
  # select strongest GMI
  dplyr::filter(padj==min(padj)) %>% 
  dplyr::ungroup() %>% 
  # keep only proximal GMIs
  dplyr::filter(dist_GEM_min<=2) %>%
  # extract metabolites
  dplyr::pull(metabolite)

dt_proximal <- res_conc %>%
  dplyr::filter(n_dataset>n_data_cut) %>%
  dplyr::filter(dist_GEM_min <= 2) %>% 
  dplyr::filter(metabolite %in% proximal_met) %>% 
  dplyr::select(metabolite,gene,padj) %>% 
  group_by(metabolite) %>%
  # dplyr::mutate(n_sign=sum(padj<pcut), d1=-log10(sort(padj)[1]), d2=-log10(sort(padj)[n_sign+1]), diff=d1-d2) %>%
  dplyr::mutate(n_sign=sum(padj<pcut), d1=-log10(sort(padj)[1]), d2=-log10(pcut), diff=d1-d2) %>%
  dplyr::mutate(label=ifelse(padj==min(padj), as.character(gene), NA)) %>% 
  ungroup() %>% 
  # dplyr::filter(diff>3) %>%
  dplyr::arrange(-d1) %>%
  dplyr::mutate(metabolite=factor(metabolite, levels = metabolite %>% unique)) %>%
  dplyr::mutate(color=ifelse(padj<pcut, "significant", "not significant"))

dt_not_proximal <- res_conc %>%
  dplyr::filter(n_dataset>n_data_cut) %>%
  dplyr::filter(padj<pcut) %>%
  dplyr::filter(dist_GEM_min > 2 | is.na(dist_GEM_min)) %>% 
  dplyr::filter(metabolite %in% unique(dt_proximal$metabolite)) %>% 
  dplyr::select(metabolite,gene,padj) %>% 
  group_by(metabolite) %>%
  # dplyr::mutate(n_sign=sum(padj<pcut), d1=-log10(sort(padj)[1]), d2=-log10(sort(padj)[n_sign+1]), diff=d1-d2) %>%
  dplyr::mutate(n_sign=sum(padj<pcut), d1=-log10(sort(padj)[1]), d2=-log10(pcut), diff=d1-d2) %>%
  dplyr::mutate(label=ifelse(padj==min(padj), as.character(gene), NA)) %>% 
  ungroup() %>% 
  # dplyr::filter(!is.na(label)) %>%
  dplyr::arrange(-d1) %>%
  dplyr::mutate(metabolite=factor(metabolite, levels = metabolite %>% unique)) %>%
  dplyr::mutate(color="not proximal")

dt_prio <- dt_proximal %>%
  rbind(dt_not_proximal) %>%
  # dplyr::filter(padj<pcut) %>%
  dplyr::group_by(metabolite) %>%
  dplyr::arrange(padj) %>%
  # only keep top 3 strongest GMIs per metabolite
  # dplyr::filter(row_number() %in% c(1:10)) %>% 
  # compute difference to strongest GMI
  dplyr::mutate(gap=-log10(padj[row_number()==1])+log10(padj[row_number()==2]))

p_proximal <- dt_prio %>%
  ggplot(aes(x=-log10(padj), y=metabolite, label=label)) +
  geom_point(aes(color=color), size=0.1) +
  geom_text_repel(size=2) +
  scale_color_manual(values=c("not_significant"="gray80",
                              "significant"="firebrick",
                              "not proximal"="black")) +
  geom_vline(xintercept = -log10(pcut), linetype="dashed", color="gray60") +
  ylab("") +
  xlab("-log10(adjusted p-value)") +
  coord_flip() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text=element_text(size=7), axis.text=element_text(size=7), axis.title=element_text(size=7), legend.position="none",
        axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(angle = 90, vjust = 1, hjust=1, color="black")) 

# p_proximal <-  dt_prio %>%
#   ggplot(aes(x=-log10(padj), y=metabolite, label=label)) +
#   geom_point(aes(color=color), size=0.1) +
#   geom_text_repel(size=2) +
#   scale_color_manual(values=c("gray80","firebrick","black")) +
#   ylab("") +
#   xlab("-log10(adjusted p-value)") +
#   coord_flip() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"),
#         text=element_text(size=7), axis.text=element_text(size=7), axis.title=element_text(size=7), legend.position="none",
#         axis.text.y = element_text(color = "black"),
#         axis.text.x = element_text(angle = 90, vjust = 1, hjust=1, color="black"))

pdf(sprintf("%s/Figure2B.pdf", Sys.Date()), width=7, height=5)
print(p_proximal)
dev.off()

save(dt_prio,file=sprintf("%s/Figure2BData.Rdata", Sys.Date()))

#### Scatter Plots ----

pp <- data.frame(metabolite= c("kynurenine","tryptophan",
                               "glutathione, oxidized (GSSG)","cystine",
                               "guanine","taurine")) %>%
  dplyr::left_join(dt_proximal, by="metabolite") %>%
  dplyr::filter(!is.na(label)) %>%
  dplyr::select(metabolite, gene)

# add negative controls
pp %<>%
  rbind(data.frame(metabolite="kynurenine",gene="AFMID")) %>%
  rbind(data.frame(metabolite="tryptophan",gene="AFMID"))

plist <- lapply(1:nrow(pp), function(k){
  
  met_name <- pp %>% dplyr::pull(metabolite) %>% .[k]
  gene_name <- pp %>% dplyr::pull(gene) %>% .[k]
  
  dd <- lapply(met_all[df %>% dplyr::filter(tissue=="Tumor") %>% dplyr::pull(dataset)] %>% names %>% {names(.)=.;.}, function(x){
    # R.utils::printf("%s: %s\n",x, met_name %in% rownames(met_all[[x]]))
    if(met_name %in% rownames(met_all[[x]])){
      data.frame(met=met_all[[x]] %>% t %>% as.data.frame %>% dplyr::pull(met_name) %>% scale,#rank %>% {./length(.)},
                 gene= rna_all[[x]] %>% t %>% as.data.frame %>% dplyr::pull(gene_name) %>% scale,#rank %>% {./length(.)},
                 dataset=df %>% dplyr::filter(dataset==x) %>% dplyr::pull(NewName))
    }
  }) %>%
    {do.call(rbind,.)} %>% as.data.frame
  
  dd %>%
    # dplyr::filter(dataset %in% c("ccRCC2","ccRCC4","BRCA2","HurthleCC","OV")) %>%
    ggplot(aes(x=met, y=gene)) +
    geom_point(aes(color=dataset), alpha=0.5, size=0.2) +
    xlab(met_name) +
    ylab(gene_name) +
    geom_smooth(aes(color=dataset), method=MASS::rlm, se=F, size=0.5) +
    # geom_smooth(method="lm", se=F, size=0.5, color="gray40") +
    scale_color_manual(values=mapping_color$Color %>% setNames(mapping_color$NewName)) +
    ggtitle(sprintf("%s - %s", met_name, gene_name)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text=element_text(size=7), axis.text=element_text(size=7), axis.title=element_text(size=7),
          legend.position="none",
          axis.text.y = element_text(color = "black"),
          axis.text.x = element_text(color="black"))
})

# Validation data from Priolo et al. (2018) doi: 10.1073/pnas.1710849115
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6142242/
val_met <- c("glutathione disulfide-nega")
ggt1_ko <- read_excel(path="Files2Share/pnas.1710849115.sd06.xlsx", sheet = "scaled data")

group <- ggt1_ko[1,2:ncol(ggt1_ko)] %>% t %>% as.data.frame %>%
  tibble::rownames_to_column("ID") %>%
  dplyr::rename(group=V1)

dt_val <- ggt1_ko %>%
  dplyr::filter(`...1`!="Label") %>%
  dplyr::rename(metabolite="...1") %>%
  dplyr::filter(metabolite %in% val_met) %>%
  tibble::column_to_rownames("metabolite") %>%
  t %>% as.data.frame %>%
  tibble::rownames_to_column("ID") %>%
  dplyr::left_join(group, by="ID") %>% 
  dplyr::filter(group %in% c("siCtrl","siGGT1")) %>%
  tidyr::pivot_longer(cols=val_met) %>% 
  dplyr::mutate(value=as.numeric(value))

sapply(val_met %>% {names(.)=.;.}, function(m){
  x <- dt_val %>%
    dplyr::filter(group=="siCtrl") %>%
    dplyr::filter(name==m) %>%
    dplyr::pull(value)
  y <- dt_val %>%
    dplyr::filter(group=="siGGT1") %>%
    dplyr::filter(name==m) %>%
    dplyr::pull(value)
  
  wilcox.test(x=x, y=y)$p.value
}) %>% as.data.frame()

plist_GGT1 <- lapply(c("glutathione disulfide-nega") %>% {names(.)=.;.}, function(m){
  m <- "glutathione disulfide-nega"
  x <- dt_val %>%
    dplyr::filter(group=="siCtrl") %>%
    dplyr::filter(name==m) %>%
    dplyr::pull(value)
  y <- dt_val %>%
    dplyr::filter(group=="siGGT1") %>%
    dplyr::filter(name==m) %>%
    dplyr::pull(value)
  
  dt_val %>%
    dplyr::filter(name==m) %>%
    dplyr::mutate(group=factor(group, levels=c("siGGT1","siCtrl"))) %>%
    ggplot(aes(x=group, y=value)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2) +
    ggtitle(sprintf("p-value:%.2e",wilcox.test(x=x, y=y)$p.value)) +
    xlab("") +
    ylab(m) +
    coord_flip() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text=element_text(size=7), axis.text=element_text(size=7), axis.title=element_text(size=7),
          legend.position="none",
          axis.text.y = element_text(color = "black"),
          axis.text.x = element_text(color="black"))
})

# save to file
pdf(sprintf("%s/Figure2CH.pdf", Sys.Date()), width=5, height=7, onefile = F)
print(
ggarrange(plotlist=list(plist[[1]],plist[[3]],plist[[11]],plist[[12]],plist[[5]],plist_GGT1[[1]]), 
          ncol=2, nrow=3,
          labels = LETTERS[3:8],
          common.legend = T)
)
dev.off()
  
#### Top Correlated Metabolites ----

# get ranked order of metabolites and genes for plotting
metlevels <- res_conc %>%
  dplyr::filter(n_dataset>n_data_cut) %>%
  dplyr::filter(padj < pcut) %>%
  dplyr::select(metabolite) %>% 
  table(dnn="metabolite") %>% as.data.frame %>%
  dplyr::arrange(-Freq) %>%
  dplyr::pull(metabolite)

# select top k correlated metabolites and genes
k <- 10
topmet <- metlevels[1:k] %>% as.character()

dt_mostcorr <- res_conc %>%
  dplyr::filter(n_dataset>n_data_cut) %>%
  dplyr::filter(padj < pcut) %>%
  dplyr::select(metabolite) %>%
  group_by(metabolite) %>%
  dplyr::summarize(Freq=n()) %>%
  dplyr::mutate(label=ifelse(metabolite %in% topmet, sprintf("%s (%d)", metabolite, Freq), NA)) %>%
  dplyr::arrange(Freq) %>% 
  dplyr::mutate(metabolite=factor(metabolite, levels = metlevels))

ppw_met <- dt_mostcorr %>%  
  ggplot(aes(x=as.integer(metabolite), y=Freq, group=1)) +
  # geom_line() +
  # geom_point(data=. %>% dplyr::filter(!is.na(label)), size=0.5) +
  geom_point(aes(color=ifelse(!is.na(dt_mostcorr$label),"hub","non_hub")),size=0.1) +
  geom_text_repel(aes(label=label), size=2, vjust=-0.5) +
  xlab("Sorted Metabolites") +
  ylab("Number of correlations") +
  scale_y_continuous(expand = c(0, 0), limits=c(0,1500)) +
  scale_x_continuous(expand = c(0, 0), limits=c(-1,nrow(dt_mostcorr))) +
  scale_color_manual(values=c("hub"="firebrick", "non_hub"="gray40")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.y = element_text(color = "black"), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        text=element_text(size=7), axis.text=element_text(size=7), axis.title=element_text(size=7),
        legend.position="none")


pdf(sprintf("%s/Figure4A.pdf", Sys.Date()), width=3, height=3, onefile=F)
print(ppw_met)
dev.off()

#### Top Correlated Genes ----

genelevels <- res_conc %>%
  dplyr::filter(n_dataset>n_data_cut) %>%
  dplyr::filter(padj < pcut) %>%
  dplyr::select(gene) %>% 
  table(dnn="gene") %>% as.data.frame %>%
  dplyr::arrange(-Freq) %>%
  dplyr::pull(gene)

# select top k correlated metabolites and genes
k <- 10
topgene <- genelevels[1:k] %>% as.character()

dt_mostcorr_gene <- res_conc %>%
  dplyr::filter(n_dataset>n_data_cut) %>%
  dplyr::filter(padj < pcut) %>%
  dplyr::select(gene) %>%
  group_by(gene) %>%
  dplyr::summarize(Freq=n()) %>%
  dplyr::mutate(label=ifelse(gene %in% topgene, sprintf("%s (%d)", gene, Freq), NA)) %>%
  dplyr::arrange(Freq) %>% 
  dplyr::mutate(gene=factor(gene, levels = genelevels))

ppw_gene <- dt_mostcorr_gene %>%  
  ggplot(aes(x=as.integer(gene), y=Freq, group=1)) +
  geom_point(aes(color=ifelse(!is.na(dt_mostcorr_gene$label),"hub","non_hub")),size=0.1) +
  geom_text_repel(aes(label=label), size=2, vjust=-0.5) +
  xlab("Sorted Genes") +
  ylab("Number of correlations") +
  scale_y_continuous(expand = c(0, 0), limits=c(0,25)) +
  scale_x_continuous(expand = c(0, 0), limits=c(-50,nrow(dt_mostcorr_gene))) +
  scale_color_manual(values=c("hub"="firebrick", "non_hub"="gray40")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.y = element_text(color = "black"), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        text=element_text(size=7), axis.text=element_text(size=7), axis.title=element_text(size=7),
        legend.position="none")


pdf(sprintf("%s/Figure4B.pdf", Sys.Date()), width=3, height=3, onefile=F)
print(ppw_gene)
dev.off()

#### Pathway Enrichment Results ----

# remove metabolites with no significant results
metlist <- res_met %>%
  dplyr::filter(padj<0.01) %>%
  dplyr::select(metabolite) %>%
  table(dnn="metabolite") %>% as.data.frame %>%
  dplyr::filter(Freq>=1) 

# remove pathways with no hit
pwlist <- res_met %>%
  dplyr::filter(padj<0.01) %>%
  dplyr::select(pathway) %>%
  table(dnn="pathway") %>% as.data.frame %>%
  dplyr::filter(Freq>=1) 

annot_row <- metlist %>% 
  dplyr::select(-Freq) %>%
  dplyr::left_join(anno_data %>% dplyr::select(metabolite, H_SUPER_PATHWAY), 
                   by="metabolite") %>%
  dplyr::rename(SUPER_PATHWAY=H_SUPER_PATHWAY) %>%
  # add number of connections
  dplyr::left_join(res_conc %>% 
                     dplyr::filter(n_dataset>7) %>%
                     dplyr::filter(padj<0.01) %>%
                     dplyr::select(metabolite) %>% 
                     table(dnn="metabolite") %>% as.data.frame, by="metabolite") %>%
  dplyr::rename(degree=Freq) %>%
  # dplyr::mutate(degreeLog=log10(degree)) %>%
  # dplyr::select(-degree) %>%
  tibble::column_to_rownames("metabolite")

# create df for heatmap
phmat <- res_met %>%
  dplyr::filter(metabolite %in% metlist$metabolite) %>%
  dplyr::filter(pathway %in% pwlist$pathway) %>%
  dplyr::select(pathway, metabolite, score) %>%
  tidyr::spread(key = metabolite, value=score) %>% 
  dplyr::arrange(pathway) %>%
  tibble::column_to_rownames("pathway")

annot_col <- res_met %>%
  dplyr::mutate(pw_class=ifelse(Class=="Metabolism", "Metabolism",
                                ifelse(Group=="Immune system","Immune system", "Other"))) %>%
  dplyr::select(pathway,pw_class) %>%
  dplyr::distinct() %>%
  tibble::column_to_rownames("pathway")

ph <- phmat %>% t %>% as.data.frame %>% 
  pheatmap(cluster_rows = T, cluster_cols = T, clustering_method = "ward.D2",
           cutree_rows = 4, cutree_cols = 7, show_colnames = T,
           annotation_row = annot_row,
           annotation_col = annot_col,
           annotation_colors = list(pw_class=c(Metabolism="thistle",
                                               "Immune system"="olivedrab3",
                                               Other="wheat"),
                                    degree=colorRampPalette(c("#E5E5E5","black"))(10)),
           breaks = c(0,2,5,10,15,20,25),
           color = c("#E5E5E5","#D4A4A4","#D4A4A4","#CB8383","#C36363","#BA4242","#B22222"))

# save to file
pdf(sprintf("%s/Figure4C.pdf", Sys.Date()), width=12, height = 8)
print(ph)
dev.off()

# save pathway annotations to file
res_met %>%
  dplyr::select(Class,Group,ID,Name) %>% 
  distinct %>% 
  write.xlsx2(file=sprintf("%s/TableS3.xlsx", Sys.Date()), sheetName = "KEGGpathways", 
              col.names = TRUE, row.names = F, append = FALSE)

#### TvsN results ----

# load RDS files
ss_gene <- readRDS("../Data/tumor_vs_normal_gene_aggregated_summary.rds")
ss_met <- readRDS("../Data/tumor_vs_normal_metabolite_aggregated_summary.rds")

# assign names to list items
nm <- sapply(ss_gene, function(x){x$pathway %>% unique})
names(ss_gene) <- nm
nm <- sapply(ss_met, function(x){x$pathway %>% unique})
names(ss_met) <- nm

# merge and clean
dt_new <- lapply(ss_gene, function(x){
  x %>% 
    dplyr::select(-name, -value) %>%
    distinct()
}) %>% {do.call(rbind,.)} %>%
  rbind(lapply(ss_met, function(x){
    x %>% 
      dplyr::select(-name, -value) %>%
      distinct()
  }) %>% {do.call(rbind,.)} 
  ) 

# find pathways with no metabolites or gene measured in 2 dataset or more
pw_remove <- dt_new %>%
  # dplyr::filter(omics=="metabolite") %>% 
  dplyr::group_by(pathway, omics) %>% 
  dplyr::summarise(var=sum(numOfItemsMeasured!=0)) %>% 
  dplyr::filter(var<=5) %>%
  dplyr::pull(pathway)

dt_new %>%
  dplyr::select(omics, pathway, dataset, numOfItemsMeasured) %>%
  tidyr::pivot_wider(id_cols = c(omics,pathway), names_from = dataset, values_from = numOfItemsMeasured) %>% 
  dplyr::arrange(pathway) %>% 
  View

# remove those pathways
dt_new %<>%
  dplyr::filter(!(pathway %in% pw_remove))

# save pathway annotations to file
dt_new %>%
  dplyr::arrange(pathway) %>% 
  write.xlsx2(file=sprintf("%s/TableS3.xlsx", Sys.Date()), sheetName = "Tumor_vs_Normal", 
              col.names = TRUE, row.names = F, append = FALSE)

# create score per dataset
dt_new %<>%
  dplyr::mutate(score_DF=ifelse(numOfItemsMeasured!=0,(numOfItemsUp+numOfItemsDown)/numOfItemsMeasured,0))

# compute correlation between gene and metabolite scores per pathways
dt_new_cor <- dt_new %>%
  dplyr::select(-numOfTotalItemsInPathway,-numOfItemsMeasured,
                -numOfItemsUp,-numOfItemsDown) %>%
  tidyr::pivot_wider(names_from = "omics", values_from = "score_DF") %>% 
  dplyr::group_by(pathway) %>%
  dplyr::mutate(cor=cor.test(gene,metabolite, method="spearman")$estimate , 
                p.value=cor.test(gene,metabolite, method="spearman")$p.value,
                ci_lower=DescTools::SpearmanRho(gene,metabolite, conf.level=0.95)[2] %>% as.vector,
                ci_upper=DescTools::SpearmanRho(gene,metabolite, conf.level=0.95)[3] %>% as.vector) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(p.value)

# add p-value adjustment
dt_new_cor %<>%
  dplyr::left_join(dt_new_cor %>%
                     dplyr::select(pathway, cor, p.value, ci_lower, ci_upper) %>%
                     distinct %>%
                     dplyr::mutate(padj=p.adjust(p.value, method = "BH"))) %>%
  mutate_at(vars(pathway), funs(factor(.,levels=unique(.)))) 

#### Plot Results ----

# plot DF score histogram per omics
dt_new %>%
  ggplot(aes(x=score_DF)) +
  geom_histogram(bins=20) +
  theme_bw() +
  facet_wrap(~omics)

# plot all pathways and datasets together
dt_new %>%
  dplyr::select(-numOfTotalItemsInPathway,-numOfItemsMeasured,
                -numOfItemsUp,-numOfItemsDown) %>%
  tidyr::pivot_wider(names_from = "omics", values_from = "score_DF") %>% 
  ggplot(aes(x=metabolite, y=gene)) +
  geom_point(aes(color=dataset), alpha=0.5) +
  geom_smooth(aes(color=dataset), method=MASS::rlm, se=F, size=0.5) +
  scale_color_manual(values=setNames(mapping_color$Color[1:7], 
                                     nm=mapping_color$OldName[1:7])) +
  theme_bw() +
  ggtitle("Differential Score across datasets and pathways")

# plot scores per pathway
dt_new_cor %>%
  ggplot(aes(x=metabolite, y=gene)) +
  geom_point(aes(color=dataset),alpha=0.5) +
  geom_smooth(aes(color=ifelse(padj<0.15,"significant","non significant")),method=MASS::rlm, se=F, size=0.5) +
  scale_color_manual(values=setNames(mapping_color$Color, nm=mapping_color$OldName)) +
  theme_bw() +
  xlab("Metabolomics DF score") +
  ylab("Transcriptomics DF score") +
  facet_wrap(~pathway)

# display correlations per pathway
dt_new_cor %>%  
  dplyr::select(pathway,cor,ci_lower,ci_upper,p.value,padj) %>%
  distinct() %>%
  dplyr::left_join(kegg_hierarchy, by=c("pathway"="Name")) %>%
  dplyr::mutate(sign=ifelse(ci_lower>0 | ci_upper<0,"significant","non significant")) %>%
  dplyr::arrange(-cor) %>%
  mutate_at(vars(pathway), funs(factor(.,levels=unique(.)))) %>%
  ggplot(aes(x=pathway,y=cor)) +
  geom_point(aes(color=Group)) +
  geom_errorbar(aes(ymin=ci_lower, ymax=ci_upper, color=sign)) +
  geom_hline(yintercept = 0, linetype="dashed", color="gray80") +
  scale_color_manual(values=setNames(mapping_color$Color[8:9], nm=mapping_color$OldName[8:9])) +
  theme_bw() +
  coord_flip() +
  ggtitle("Spearman correlation of gene-metabolite DF scores across datasets")

#### Ed's figure ----

genesum = data.frame()
for (ii in 1:length(ss_gene)){
  temp = ss_gene[[ii]][,c(3,5:9)]
  temp = temp[-which(duplicated(temp)),]
  genesum = rbind(genesum,temp)
  
}

genesum$DF = (genesum$numOfItemsUp-genesum$numOfItemsDown)/genesum$numOfItemsMeasured
#genesum = genesum[-which(is.na(genesum$DF)),] # no missing data
genesum$omics = 'Transcript'

# score each pathway by its mean
pwayscore = c()
for (ii in unique(genesum$pathway)){
  pwayscore[ii] = mean(genesum[which(genesum$pathway == ii),'DF'])
}

metsum = data.frame()
for (ii in 1:length(ss_gene)){
  temp = ss_met[[ii]][,c(3,5:9)]
  temp = temp[-which(duplicated(temp)),]
  metsum = rbind(metsum,temp)
}
metsum$DF = (metsum$numOfItemsUp-metsum$numOfItemsDown)/metsum$numOfItemsMeasured
metsum = metsum[-which(is.na(metsum$DF)),]
metsum$omics = 'Metabolite'

allsum = rbind(metsum,genesum)
allsum$studypway = paste(allsum$dataset,allsum$pathway,sep =':')

# calculate the color of the segments
dfdiff = data.frame(matrix(NA,0,4),stringsAsFactors = FALSE)
for (ii in unique(allsum$studypway)){
  tempg = allsum[which(allsum$studypway == ii & allsum$omics == 'Transcript'),]
  tempm = allsum[which(allsum$studypway == ii & allsum$omics == 'Metabolite'),]
  dfdiff[ii,] = c(strsplit(ii,'\\:')[[1]][1],strsplit(ii,'\\:')[[1]][2],tempg[1,'DF'],tempm[1,'DF'])
}
colnames(dfdiff) = c('dataset','pathway','DFgene','DFmet')
dfdiff$DFgene = as.numeric(dfdiff$DFgene)
dfdiff$DFmet = as.numeric(dfdiff$DFmet)
dfdiff$sign = ifelse(sign(dfdiff$DFgene) == sign(dfdiff$DFmet),1,-1)
dfdiff$x = 'Metabolite'
dfdiff$xend = 'Transcript'

pwayscore2 = c()
for (ii in unique(genesum$pathway)){
  pwayscore2[ii] = sum(dfdiff[which(dfdiff$pathway == ii),'sign'],na.rm = TRUE)
}

# order the pathways 
allsum$pathway = factor(allsum$pathway,levels = names(pwayscore2)[order(pwayscore2)])

# order the pathways 
pdf(sprintf("%s/TvsN_Dotplot.pdf", Sys.Date()),height = 8,width = 12)
print(
allsum %>%
  dplyr::filter(pathway %in% (dt_new_cor %>% dplyr::pull(pathway) %>% unique)) %>% 
  ggplot(aes(x=omics, y=pathway,color = DF)) + 
  geom_point(aes(size=log10(numOfItemsMeasured))) + 
  scale_size_continuous(range = c(0.1,3)) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_color_gradient2() +
  theme_minimal() + 
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()) +
  xlab("") +
  ylab("") +
  geom_segment(data = dfdiff[which(dfdiff$sign == 1),],aes(x = x,xend = xend,y = pathway,yend = pathway),color = 'black',linetype = 'solid') +
  geom_segment(data = dfdiff[which(dfdiff$sign == -1),],aes(x = x,xend = xend,y = pathway,yend = pathway),color = 'gray',linetype = 'dashed') +
  facet_wrap(~dataset,nrow = 1) 
)
dev.off()

genesum = data.frame()
for (ii in 1:length(ss_gene)){
  temp = ss_gene[[ii]][,c(3,5:9)]
  if (any(duplicated(temp))){
    temp = temp[-which(duplicated(temp)),]
  }
  genesum = rbind(genesum,temp)
}

genesum$DF = (genesum$numOfItemsUp-genesum$numOfItemsDown)/genesum$numOfItemsMeasured
#genesum = genesum[-which(is.na(genesum$DF)),] # no missing data
genesum$omics = 'Transcript'

# score each pathway by its mean
pwayscore = c()
for (ii in unique(genesum$pathway)){
  pwayscore[ii] = mean(genesum[which(genesum$pathway == ii),'DF'])
}

metsum = data.frame()
for (ii in 1:length(ss_gene)){
  temp = ss_met[[ii]][,c(3,5:9)]
  if (any(duplicated(temp))){
    temp = temp[-which(duplicated(temp)),]
  }
  metsum = rbind(metsum,temp)
}
metsum$DF = (metsum$numOfItemsUp-metsum$numOfItemsDown)/metsum$numOfItemsMeasured
metsum = metsum[-which(is.na(metsum$DF)),]
metsum$omics = 'Metabolite'

allsum = rbind(metsum,genesum)
allsum$studypway = paste(allsum$dataset,allsum$pathway,sep =':')

# calculate the color of the segments
dfdiff = data.frame(matrix(NA,0,4),stringsAsFactors = FALSE)
for (ii in unique(allsum$studypway)){
  tempg = allsum[which(allsum$studypway == ii & allsum$omics == 'Transcript'),]
  tempm = allsum[which(allsum$studypway == ii & allsum$omics == 'Metabolite'),]
  dfdiff[ii,] = c(strsplit(ii,'\\:')[[1]][1],strsplit(ii,'\\:')[[1]][2],tempg[1,'DF'],tempm[1,'DF'])
}
colnames(dfdiff) = c('dataset','pathway','DFgene','DFmet')
dfdiff$DFgene = as.numeric(dfdiff$DFgene)
dfdiff$DFmet = as.numeric(dfdiff$DFmet)
dfdiff$sign = ifelse(sign(dfdiff$DFgene) == sign(dfdiff$DFmet),1,-1)
dfdiff$x = 'Metabolite'
dfdiff$xend = 'Transcript'

pwayscore2 = c()
for (ii in unique(genesum$pathway)){
  pwayscore2[ii] = sum(dfdiff[which(dfdiff$pathway == ii),'sign'],na.rm = TRUE)
}

pwlist <- (dt_new_cor %>% dplyr::pull(pathway) %>% as.character %>% unique)

# filter pathways
pwayscore2 <- pwayscore2[names(pwayscore2) %in% pwlist]
allsum %<>%
  dplyr::filter(pathway %in% pwlist)

# add score to dataframe
allsum %<>%
  dplyr::left_join(pwayscore2 %>% 
                     as.data.frame %>% 
                     dplyr::rename(score=".") %>% 
                     tibble::rownames_to_column("pathway"),
                   by="pathway") %>%
  # order pathways by score
  dplyr::arrange(score) %>%
  dplyr::mutate(pathway=factor(pathway, levels=pathway %>% unique))


pdf("TvsN_Dotplot.pdf",height = 8,width = 12)
print(
  allsum %>% 
    ggplot(aes(x=omics, y=pathway,color = DF)) + 
    geom_point(aes(size=log10(numOfItemsMeasured))) + 
    scale_size_continuous(range = c(0.1,3)) +
    scale_x_discrete(guide = guide_axis(angle = 90)) +
    scale_color_gradient2() +
    theme_minimal() + 
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank()) +
    xlab("") +
    ylab("") +
    geom_segment(data = dfdiff[which(dfdiff$sign == 1),],aes(x = x,xend = xend,y = pathway,yend = pathway),color = 'black',linetype = 'solid') +
    geom_segment(data = dfdiff[which(dfdiff$sign == -1),],aes(x = x,xend = xend,y = pathway,yend = pathway),color = 'gray',linetype = 'dashed') +
    facet_wrap(~dataset,nrow = 1) 
)
dev.off()

save(pwlist, file=sprintf("%s/Pathwaylist.Rdata", Sys.Date()))

#### Pathway Class Enrichment Analysis ----

# plot only significant pathways
dt_new_cor %>%
  dplyr::filter(padj<0.15) %>% 
  ggplot(aes(x=metabolite, y=gene)) +
  geom_point(aes(color=dataset),alpha=0.5) +
  geom_smooth(aes(color=ifelse(padj<0.15,"significant","non significant")),method=MASS::rlm, se=F, size=0.5) +
  scale_color_manual(values=setNames(mapping_color$Color, nm=mapping_color$OldName)) +
  theme_bw() +
  xlab("Metabolomics DF score") +
  ylab("Transcriptomics DF score") +
  facet_wrap(~pathway,ncol = 5)

# visualize significant pathways
dt_new_cor %>% 
  dplyr::filter(padj<0.15) %>% 
  dplyr::left_join(kegg_hierarchy, by=c("pathway"="Name")) %>% 
  dplyr::select(Group, pathway, cor, p.value, padj) %>% 
  distinct() %>% View
dplyr::filter(Group=="Amino acid metabolism") %>% 
  dplyr::pull(pathway) %>% View

# visualize results for the two pathway classes of interest
dt_new_cor %>% 
  dplyr::left_join(kegg_hierarchy, by=c("pathway"="Name")) %>% 
  dplyr::select(Group, pathway, cor, p.value, padj) %>% 
  distinct() %>% 
  dplyr::filter(Group %in% (dt_new_cor %>% 
                              dplyr::filter(padj<0.05) %>% 
                              dplyr::left_join(kegg_hierarchy, by=c("pathway"="Name")) %>%
                              dplyr::pull(Group) %>% unique)) %>% View

# visualize Group distribution across all pathways tested
dt_new %>% 
  dplyr::left_join(kegg_hierarchy, by=c("pathway"="Name")) %>% 
  dplyr::select(pathway, Group) %>% 
  distinct() %>% 
  dplyr::select(Group) %>%
  table

# enrichment dataset
dt_pw_enr <- dt_new_cor %>%
  dplyr::left_join(kegg_hierarchy, by=c("pathway"="Name")) %>%
  dplyr::select(pathway,Group,padj) %>%
  dplyr::distinct() %>%
  dplyr::mutate(is_sign=ifelse(padj<0.05,1,0) %>% factor(levels = c(1,0)),
                is_aa=ifelse(Group=="Amino acid metabolism",1,0) %>% factor(levels = c(1,0)),
                is_carb=ifelse(Group=="Carbohydrate metabolism",1,0) %>% factor(levels = c(1,0)))

# Fisher's test for Amino acid metabolism
dt_pw_enr %>%  
  dplyr::select(is_sign,is_aa) %>% 
  table() %>% fisher.test %>% .$p.value

# Fisher's test for Carbohydrate metabolism
dt_pw_enr %>%  
  dplyr::select(is_sign,is_carb) %>% 
  table() %>% fisher.test %>% .$p.value

#### Paper figures ----

# plot scores for selected pathways
pw <- "Citrate cycle (TCA cycle)"

pdf(sprintf("%s/TvsN_Pathway_%s.pdf", Sys.Date(), pw), width=2.5, height=2.5, onefile = F)
print(
  dt_new_cor %>%
    dplyr::filter(pathway == pw) %>%
    ggplot(aes(x=metabolite, y=gene, label=dataset)) +
    geom_point(aes(color=dataset),alpha=0.5) +
    geom_smooth(aes(color=ifelse(padj<0.05,"significant","non significant")),method=MASS::rlm, se=F, size=0.2) +
    scale_color_manual(values=setNames(mapping_color$Color, nm=mapping_color$OldName)) +
    geom_text_repel(size=2) +
    theme_bw() +
    xlab("Metabolomics DF score") +
    ylab("Transcriptomics DF score") +
    ggtitle(pw) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text=element_text(size=7), axis.text=element_text(size=7), axis.title=element_text(size=7),
          legend.position="none",
          axis.text.y = element_text(color = "black"),
          axis.text.x = element_text(color="black"))
)
dev.off()

pdf(sprintf("%s/TvsN_Pathway_Cor.pdf", Sys.Date()), width=5, height=3, onefile = F)
print(
  dt_new_cor %>%  
    dplyr::select(pathway,cor,ci_lower,ci_upper,p.value,padj) %>%
    distinct() %>%
    dplyr::left_join(kegg_hierarchy, by=c("pathway"="Name")) %>%
    dplyr::mutate(sign=ifelse(ci_lower>0 | ci_upper<0,"significant","non significant")) %>%
    dplyr::arrange(-cor) %>%
    dplyr::mutate(label=ifelse(sign=="significant", pathway, NA)) %>% 
    mutate_at(vars(pathway), funs(factor(.,levels=unique(.)))) %>%
    ggplot(aes(x=pathway,y=cor, label=label)) +
    geom_point(aes(color=Group)) +
    geom_errorbar(aes(ymin=ci_lower, ymax=ci_upper, color=sign)) +
    geom_hline(yintercept = 0, linetype="dashed", color="gray60") +
    scale_color_manual(values=setNames(mapping_color$Color[8:9], nm=mapping_color$OldName[8:9])) +
    theme_bw() +
    # coord_flip() +
    geom_label_repel(box.padding = 0.5,size=2, fill="white") +
    ggtitle("Spearman correlation of gene-metabolite DF scores across datasets") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text=element_text(size=7), axis.text=element_text(size=7), axis.title=element_text(size=7),
          legend.position="none",
          axis.text.x = element_blank(), axis.ticks.x=element_blank(),
          axis.text.y = element_text(color = "black"))
)
dev.off()

# # save results to file
# library(xlsx)
# res <- dt_new_cor %>%  
#   dplyr::select(pathway,cor,ci_lower,ci_upper,p.value,padj) %>%
#   distinct() %>%
#   dplyr::left_join(kegg_hierarchy, by=c("pathway"="Name")) %>%
#   dplyr::arrange(-cor) %>% 
#   dplyr::select(-match_name) %>% as.data.frame
# write.xlsx2(x=res, file=sprintf("%s/TableS17.xlsx", Sys.Date()), sheetName = "TvsN_PathwaySpearman", 
#            col.names = TRUE, row.names = F, append = FALSE)


