markers <- c("CD19","CD3", "CD8","CD4","CD45RA","CD45RO","CD16","CD14","CD56","CD57","CD11b","CD11c","CD27","CD28","CD192_CCR2",
             "CD194_CCR4","CD195_CCR5","CD197_CCR7","CD86","HLA-DR","CD183_CXCR3","CD185_CXCR5","CD123","TGFb","IFNg","IFNB","IL-10",
             "IL-17A","IL-2","IL-1b","TNFa","T-bet","GranzymeB","CD25","CD80","CD33","CD38","CD44","CD127","PD-1","Mip1B_CCL4")

format_data <- function(abund, pheno){
  # function to format input cell abundances and cell phenotypes to a list output
  # containing abundances, phenotypes and sample key
  # abund = "../Data preprocessing/formatted data/220919_cluster_abundances_raw.csv"
  # pheno = "../Data preprocessing/formatted data/220919_cluster_phenotypes_raw.csv"
  abund = data.table::fread(abund)
  abund = setNames(data.frame(t(abund[,-1])), unlist(abund[,1]))
  
  pheno = data.table::fread(pheno)
  colnames(pheno)[1:2] <- c("sample", "Cluster")
  colnames(pheno) <- gsub(" median", "", colnames(pheno))
  
  return(list(abund, pheno))
}


nicer_volcano <- function(main_title, legend_title, name_1, name_2){
  volcano_df = res@results
  row.names(volcano_df) = NULL
  # cluster_size = setNames(tibble::rownames_to_column(data.frame(res@cluster.size)), c("cluster", "cluster_size"))
  # volcano_df = dplyr::left_join(cluster_size, 
  #                               res@results, 
  #                               by="cluster")
  volcano_df$log2FC = log(volcano_df$mean.samples1/volcano_df$mean.samples2, 2)
  volcano_df

  cols <- c("B" = "#0944a8", "I" = "#008000", "T" = "#a11313", "ns" = "#999999")  
  volcano <- ggplot(volcano_df, aes(log2FC, -log(pvalue,10))) +
    geom_point(aes(color = cell, size = Count)) + 
    scale_color_manual(values=cols) +
    xlab(paste("log2FC\n", name_2, " < enriched > ", name_1)) +
    ylab("-log10pval") +
    ggtitle(main_title) +
    geom_text_repel(size=3.5, data=subset(volcano_df, significant=="TRUE"), 
                    aes(y=-log(pvalue,10), label = annotated), 
                    box.padding = 0.2, force=20, point.padding = 0.1, 
                    min.segment.length=0.05, max.overlaps = 20) +
    theme_minimal() +
    theme(axis.text=element_text(size=8),
          axis.title=element_text(size=9, face="bold"), 
          legend.text=element_text(size=8, face="italic"), 
          legend.title=element_text(size=8, face="italic"), 
          plot.title = element_text(size=10, hjust = 0.5)) +
    guides(color = guide_legend(order=1),
           size = guide_legend(order=2, title=legend_title)) +
    geom_hline(yintercept=-log(0.05,10), linetype="dashed", color="#800020") +
    geom_vline(xintercept=log(1.5, 2), linetype="dashed", color="#800020") +
    geom_vline(xintercept=-(log(1.5, 2)), linetype="dashed", color="#800020") + 
    lemon::scale_x_symmetric(mid=0)
  
  show(volcano)
  return(volcano)
}

nasty_function_for_PCPs <- function(Results,
                                    samples        = NULL,
                                    clusters       = NULL,
                                    markers        = NULL,
                                    show.mean      = "both",
                                    show.on_device = TRUE,
                                    sort.markers   = TRUE) {
  
  if (is.null(Results)) {
    stop("Error in phenoViewer: 'Results' parameter can not be NULL")
  } else if (class(Results)[1] != "Results") {
    stop("Error in phenoViewer: 'Results' parameter must be a 'Results' object")
  }
  
  if(length(Results@marker.names) == 0){
    stop("Error in phenoViewer: 'Results' object must contain phenotypes")
  }
  
  if (is.null(samples)) {
    samples     <- Results@sample.names
    data        <- Results@cluster.phenotypes
    cluster.abundances <- Results@cluster.abundances
  } else if (!all(samples %in% Results@sample.names)) {
    stop("Error in phenoViewer: 'samples' parameter must contains only samples names\n Unknown sample names: ",
         paste(setdiff(unique(samples), Results@sample.names), collapse = " "))
  } else {
    data               <- subset(Results@cluster.phenotypes, sample %in% samples, drop = FALSE)
    cluster.abundances <- Results@cluster.abundances[, samples, drop = FALSE]
  }
  
  data <- stats::na.omit(data)
  
  if (is.null(clusters)) {
    stop("Error in phenoViewer: 'clusters' parameter is required")
  } else if (all(clusters %in% Results@cluster.names)) {
    if (typeof(clusters) != "character") {
      stop("Error in phenoViewer: 'clusters' parameter must be a character vector")
    }
    clusters        <- unique(clusters)
    clusters.select <- data[, "cluster"] %in% clusters
    data            <- data[clusters.select,]
    cluster.abundances     <- cluster.abundances[clusters,]
  } else {
    stop("Error in phenoViewer:\nUnknown clusters : ", paste(setdiff(unique(clusters), Results@cluster.names), collapse = " "))
  }
  
  data <- plyr::ddply(data, c("sample"), function(df) {
    apply(df[, 3:ncol(df)], 2, mean, na.rm = TRUE)
  }) 
  
  if (is.null(markers)) {
    markers <- Results@marker.names
  } else if (all(markers %in% Results@marker.names)) {
    markers <- unique(markers)
    data <- data[, c("sample", markers)]
  } else {
    stop("Error in phenoViewer: Unknown markers :", paste(setdiff(unique(markers), Results@marker.names), collapse = " "))
  }
  
  if (show.mean != "none" && show.mean != "both" && show.mean != "only") {
    stop("Error in phenoViewer: 'show.mean' parameter must contain only one of these : 'none', 'both' or 'only'")
  }
  
  if (!is.logical(show.on_device)) { stop("Error in phenoViewer: 'show.on_device' parameter must be a logical") }
  
  data           <- reshape2::melt(data, id = c("sample"), stringsAsFactors = FALSE)
  colnames(data) <- c("samples", "marker", "value")
  
  names.palette  <- unique(Results@cluster.phenotypes[, markers]$sample)
  palette        <- ggcolors(length(names.palette))
  names(palette) <- names.palette
  
  assignments <- Results@assignments
  
  if (!is.null(assignments)) {
    
    order       <- unique(assignments$bc)
    assignments <- assignments[samples, , drop = FALSE]
    data$bc <- assignments[data$samples, "bc"]
    order       <- intersect(order, unique(assignments$bc))
    data$bc <- factor(data$bc, levels = order)
    
    names.palette  <- unique(assignments$bc)
    palette        <- ggcolors(length(names.palette))
    names(palette) <- names.palette
    
  } else if (is.element("bc", colnames(assignments))) {
    warning("Warning in phenoViewer: 'assignments' slot do not contain the column 'bc' in the provided 'Results' object. Consequently, the samples names will be used in remplacement")
  } else {
    warning("Warning in phenoViewer: 'assignments' slot in the provided 'Results' object is absent. Consequently, the samples names will be used in remplacement")
  }
  
  
  if(sort.markers==TRUE){
    clustering.markers  <- Results@clustering.markers[Results@clustering.markers %in% markers]
    ordered.markers     <- c(gtools::mixedsort(clustering.markers),gtools::mixedsort(setdiff(Results@marker.names[Results@marker.names %in% markers], clustering.markers)))
    bold.markers        <- ifelse(is.element(ordered.markers, clustering.markers), "bold", "plain")
    colored.markers     <- ifelse(is.element(ordered.markers, clustering.markers), "blue", "black")
    data$marker         <- factor(data$marker, 
                                  levels = ordered.markers,
                                  #levels = markers,
                                  ordered = TRUE)
  }else{
    clustering.markers  <- Results@clustering.markers
    ordered.markers     <- markers
    bold.markers        <- ifelse(is.element(ordered.markers, clustering.markers), "bold", "plain")
    colored.markers     <- ifelse(is.element(ordered.markers, clustering.markers), "blue", "black")
    data$marker         <- factor(data$marker, levels = ordered.markers, ordered = TRUE)
  }
  
  Results@bounds <- Results@bounds[, colnames(Results@bounds) %in% markers]   
  for (i in seq_len(nrow(data))) {
    data[i, "lower.bound"] <- Results@bounds[1, as.character(data[i, "marker"])]
    data[i, "upper.bound"] <- Results@bounds[2, as.character(data[i, "marker"])]
  }
  
  cells.number <- sum(colSums(cluster.abundances))
  
  title    <- paste("Pheno Viewer - cluster: ", paste0(clusters, collapse = ", "), " (", format(cells.number, big.mark = " "), " cells)", sep = "")
  bounds   <- as.numeric(row.names(Results@bounds))
  subtitle <- paste0("Grey ribbon displays from ", (bounds[1] * 100), "% to ", (bounds[2] * 100), "% percentiles of the range expression")
  
  max.value <- -1
  min.value <- -1
  
  max.value <- max(c(data$value, data$upper.bound), na.rm = TRUE)
  min.value <- min(c(data$value, data$lower.bound), na.rm = TRUE)
  
  max.value <-  max.value * (1 + sign(max.value) * 0.1)
  min.value <-  min.value * (1 - sign(min.value) * 0.1)
  
  
  means <- plyr::ddply(data,
                       c("marker"),
                       function(df){mean(df$value, na.rm = TRUE)})
  colnames(means) <- c("marker", "means")
  
  data_means <- data.frame(marker = 0, means= 0, clusters = 0)
  tmp_clusters<- unique(cluster.phenotypes$Cluster) ###### make sure the cluster.phenotypes file column name is "Cluster" and not "cluster"
  for(i in tmp_clusters){
    tmp_data<- Results@cluster.phenotypes[, c("sample", "cluster", markers)]
    tmp_clusters.select <- tmp_data[, "cluster"] %in% i
    tmp_data <- tmp_data[tmp_clusters.select,]
    
    tmp_data <- plyr::ddply(tmp_data, c("sample"), function(df) {
      apply(df[, 3:ncol(df)], 2, mean, na.rm = TRUE)
    }) 
    
    tmp_data           <- reshape2::melt(tmp_data, id = c("sample"), stringsAsFactors = FALSE)
    colnames(tmp_data) <- c("samples", "marker", "value")
    
    tmp_means <- plyr::ddply(tmp_data,
                             c("marker"),
                             function(df){mean(df$value, na.rm = TRUE)})
    colnames(tmp_means) <- c("marker", "means")
    tmp_means$clusters = i
    
    data_means = rbind(data_means, tmp_means)
  }
  data_means = data_means[-1, ]
  data_means$marker = factor(data_means$marker, levels = ordered.markers)
  #data_means$marker = substr(data_means$marker, 2, 100000) ##
  #data_means = data_means[order(data_means$marker, decreasing = TRUE), ] ##
  plot <- ggplot2::ggplot(data = data_means) +
    ggplot2::ggtitle(bquote(atop(.(title), atop(italic(.(subtitle)), ""))))
  
  plot  <- plot + ggplot2::geom_line(ggplot2::aes_string(x = "marker", y = "means", group = "clusters"),
                                     size = 0.5, #changes size of background lines
                                     alpha = 1,
                                     color = "#CCCCCC")+ 
    ggplot2::scale_y_continuous(limits = c(min.value, max.value), breaks = round(seq(0, max.value, by = 1), 0)) +
    ggplot2::theme_bw()
  
  plot <- plot + ggplot2::geom_line(data  = means,
                                    ggplot2::aes_string(x = "marker", y = "means", group = 1),
                                    #group = 1,
                                    linetype = "solid",
                                    size  = 1,
                                    color = "#21BE75") 
  
  plot <- plot + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, face = bold.markers, color = colored.markers)) +
    ggplot2::theme(legend.text = ggplot2::element_text(size = 6),
                   legend.key  = ggplot2::element_blank(),
                   plot.title  = ggplot2::element_text(hjust=0.5)) +
    ggplot2::xlab("markers") +
    ggplot2::ylab("marker expressions") +
    ggplot2::guides(col = ggplot2::guide_legend(ncol = 1))
  
  
  grid::grid.draw(plot)
  invisible(plot)
}