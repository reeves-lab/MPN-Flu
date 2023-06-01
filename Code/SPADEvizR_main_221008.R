library(glue)
library("readxl")
library(ggplot2)
library(ggrepel)
library(data.table)
library(SPADEVizR)
library(plotflow)
source("./Code/SPADEvizR_comparisons_221008.R")
source("./Code/SPADEvizR_functions_221101.R")


#.................................PARALLEL COORDINATE PLOTS ON UNCORRECTED DATA.................................#
abund_file = "./Data/omiqData_formatted/230324_cluster_abundances.csv"
pheno_file = "./Data/omiqData_formatted/230324_cluster_phenotypes.csv"
formatted_list = format_data(abund=abund_file, pheno=pheno_file)
cluster.abundances = formatted_list[[1]]
cluster.phenotypes = formatted_list[[2]]

results <- importResultsFromTables(cluster.abundances = cluster.abundances, 
                                   cluster.phenotypes = as.data.frame(cluster.phenotypes))

for(i in 1:nrow(cluster.abundances)){
  jpeg(paste("./Graphs/Parallel Coordinate Plots/All/", rownames(cluster.abundances)[i], ".jpeg", sep = ""),
       width=2000,
       height=1500, 
       res = 300)
  nasty_function_for_PCPs(results, clusters=rownames(cluster.abundances)[i])
  dev.off()
}

# Cell specific PCPs
abund_file = "./Data/omiqData_formatted/230324_cluster_abundances.csv"
pheno_file = "./Data/omiqData_formatted/230324_cluster_phenotypes.csv"
formatted_list = format_data(abund=abund_file, pheno=pheno_file)
cluster.abundances = formatted_list[[1]]
cluster.phenotypes = formatted_list[[2]]
cluster.abundances = cluster.abundances[grep('^T', rownames(cluster.abundances)),]
cluster.phenotypes = cluster.phenotypes[grep('^T', cluster.phenotypes$Cluster),]

results <- importResultsFromTables(cluster.abundances = cluster.abundances, 
                                   cluster.phenotypes = as.data.frame(cluster.phenotypes))

for(i in 1:nrow(cluster.abundances)){
  jpeg(paste("./Graphs/Parallel Coordinate Plots/Cell specific/", rownames(cluster.abundances)[i], ".jpeg", sep = ""),
       width=2000,
       height=1500, 
       res = 300)
  nasty_function_for_PCPs(results, clusters=rownames(cluster.abundances)[i])
  dev.off()
}


#.................................DIFFERENTIAL CLUSTER ABUNDANCE ANALYSIS................................#
abund_file = "./Data/omiqData_formatted/230324_cluster_abundances.csv"
pheno_file = "./Data/omiqData_formatted/230324_cluster_phenotypes.csv"

formatted_list = format_data(abund=abund_file, pheno=pheno_file)
cluster.abundances = formatted_list[[1]]
cluster.phenotypes = formatted_list[[2]]

results <- importResultsFromTables(cluster.abundances = cluster.abundances, 
                                   cluster.phenotypes = as.data.frame(cluster.phenotypes))

metadata <- read_excel("./Data/221027_metadata.xlsx", sheet = "Combined")
cols = c("Sample ID", "Diagnosis")
metadata = metadata[ , names(metadata) %in% cols]
metadata["Diagnosis"][is.na(metadata["Diagnosis"])] <- 'Healthy'
metadata$Diagnosis <- gsub('MPN', 'PV_ET', metadata$Diagnosis)
metadata$`Sample ID` <- gsub('$', '_Unstim', metadata$`Sample ID`)

annotations <- read_excel("./Data/Cluster_phenotypes.xlsx")
colnames(annotations)[1] <- 'cluster'


# iterate over comparisons list of lists and get DACs for each, concat them into results_df
directory = "./Graphs/Differential Abundance/"
results_df = data.frame()
for (i in group_effect){
  # i = group_effect[[2]]
  comparison = i[[1]]
  name1 = i[[2]]
  name2 = i[[3]]
  
  samples1 = metadata[metadata[comparison] == name1,]$`Sample ID`
  samples2 = metadata[metadata[comparison] == name2,]$`Sample ID`
  
  # if there are enough observations in samples1 and samples2, then execute DAC
  if(length(samples1) > 1 & length(samples2) > 1 ){
    res <- identifyDAC(results, condition1 = samples1, condition2 = samples2, method.adjust = "fdr", th.fc = 1.5)
    table_res <- res@results[res@results$significant == TRUE, c(1,6,7)]
    
    if (dim(table_res)[1] != 0) {
      table_res$comparison = glue("{name2}_vs_{name1}")
      results_df = rbind(results_df, table_res)
    }
    
    res@results <- dplyr::left_join(res@results, annotations[, c('cluster', 'Cell type', 'Count')], by='cluster')
    res@results$annotated <- paste0(res@results$`Cell type`, ' (', res@results$cluster, ')')
    res@results$cell <- gsub('^(.){1}.*?$', '\\1', res@results$cluster)
    res@results$cell[res@results$significant == FALSE] <- 'ns'

    # source("./Code/SPADEvizR_functions_221101.R")
    tiff(glue('{directory}/{name2}_vs_{name1}.tiff'), width = 6, height = 5.5, units = 'in', res = 300)
      nicer_volcano(legend_title = "Cluster size", 
                    main_title = glue("{name2}_vs_{name1}"), 
                    name_1 = name1, name_2 = name2)
    dev.off()
    
  }else{
    print("Oh no")
  }
}

write.csv(results_df, './Data/DAC_results.csv', row.names = FALSE)
