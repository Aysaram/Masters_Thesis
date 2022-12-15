# Moves this to a new script file in folder function
## Run chemdiv library

#BiocManager::install("fmcsR") # Needs to be run to be able to install chemodiv
# install.packages("chemodiv")
library(chemodiv)
library(grid)
library(gridExtra)
# First run chemoDivCheck, to see that data is correctly inputted. 
library(ggtree)

source("Scripts/Functions/Mol_net_graph.R")

Chemodiv_build <- function(Data_in) {
  To_run_samp_data <- Data_in %>% 
    mutate(Sample = if_else(is.na(`Sample descriptor`), Sample_name, `Sample descriptor`),
           Area.from.MS = as.numeric(Area.from.MS)) %>%
    filter(!is.na(Compound.Name), !is.na(Area.from.MS), !is.na(Smiles),
           !is.na(ST_InChiKey)) %>%  # Describe how many how many is removed.
    group_by(Compound.Name, Sample) %>% mutate(Row_number = row_number()) %>% ungroup() %>% 
    filter(Row_number < 2) %>% 
    dplyr::select(Compound.Name, Area.from.MS, Sample) %>% 
    group_by(Sample) %>% mutate(Total_sum = sum(Area.from.MS, na.rm = T)) %>% ungroup() %>% 
    mutate(Ratio = Area.from.MS/Total_sum) %>% dplyr::select(-Total_sum, -Area.from.MS) %>% 
    arrange(Compound.Name) %>% 
    pivot_wider(names_from = Compound.Name, values_from = Ratio) %>% 
    mutate_all(~replace(., is.na(.), 0))
  
  To_run_comp_data <- Data_in %>% 
    mutate(Sample = if_else(is.na(`Sample descriptor`), Sample_name, `Sample descriptor`),
           Area.from.MS = as.numeric(Area.from.MS)) %>%
    filter(!is.na(Compound.Name), !is.na(Area.from.MS), !is.na(Smiles),!is.na(ST_InChiKey)) %>% 
    group_by(Compound.Name, ST_InChiKey, Smiles) %>% 
    summarise() %>% ungroup() %>% 
    rename(smiles = Smiles, compound = Compound.Name, inchikey = ST_InChiKey) %>% 
    mutate(inchikey = stri_replace_first_fixed(inchikey, "InChIKey=", ""))
  
  Return_list <- list()
  Return_list[[1]] <- To_run_samp_data
  Return_list[[2]] <- To_run_comp_data
  
  return(Return_list)
}
# incomplete needs to stop if error
NPclass <- function(Data_in = Sampdata) {
  
  Sample_data <- Data_in[[1]][2:length(Data_in[[1]])] # Subsets the lenght
  chemodiv::chemoDivCheck(Sample_data, Data_in[[2]])
  
  # If ok run this: If not return error # needs code
  NPC_tibble <- NPCTable(Data_in[[2]]) # Runs the NPclassifier 
  return(NPC_tibble)
}
# Runs NPclassifier
Mol_network <- function(Data_in) {
  Sample_data <- Data_in[[1]][2:length(Data_in[[1]])] # Subsets the lenght
  
  Compound_dissim <- compDis(Data_in[[2]])
  
  Mol_network <- molNet(compDisMat = Compound_dissim$fingerDisMat) # This function is broken, need to select the dissim index from the list
  Return_list <- list()
  Return_list[[1]] <- Sample_data
  Return_list[[2]] <- Compound_dissim
  Return_list[[3]] <- Mol_network
  return(Return_list)
}
# Builds data tibbles correctly
Sampdata <- Chemodiv_build(Data_in = MS_data %>% 
                             filter(!(`Sample descriptor` %in% c("Fingerprint Blend 1", "Fingerprint Blend 2", "Fingerprint Blend 3"))))
# Runs the NPClassifier
NP_classified <- NPclass()
Mol_net <- Mol_network(Data_in = Sampdata)

Remove_guides <- function(x) {
  x + guides(fill = "none", size = "none", width = "none", alpha = "none") +
    ggraph::scale_edge_width(guide = "none")
} # Remove any guides from ggplot object

Mol_net_plot <- Molnet_plot2(sampleData = Mol_net[[1]], 
             networkObject = Mol_net[[3]]$networkObject, 
             groupData = Sampdata[[1]]$Sample, # Group_data
             npcTable = NP_classified, 
             layout = "kk")
Mol_net_plot_large <- Molnet_plot2(sampleData = Mol_net[[1]], 
                             networkObject = Mol_net[[3]]$networkObject, 
                             #groupData = Sampdata[[1]]$Sample, # Group_data
                             npcTable = NP_classified, 
                             layout = "kk")

Molnet_plot_large <- grid.arrange(Mol_net_plot[[4]] %>% Remove_guides(),
                                  Mol_net_plot[[2]] %>% Remove_guides(),
                                  Mol_net_plot[[3]] %>% Remove_guides(),
                                  Mol_net_plot[[1]] %>% Remove_guides(),
                                  ggpubr::get_legend(Mol_net_plot_large, position = NULL),
                                  layout_matrix = rbind(c(1,1,2,2,5,5,5),
                                                        c(1,1,2,2,5,5,5),
                                                        c(3,3,4,4,5,5,5),
                                                        c(3,3,4,4,5,5,5)))

# Calculate key values

Sample_data <- Sampdata[[1]][2:length(Sampdata[[1]])] # Subsets the lenght

HillDiv <- calcDiv(sampleData = Mol_net[[1]], compDisMat = Mol_net[[2]]$fingerDisMat) #calculates alpha diversity and evenness.

calcBetaDiv(sampleData = Mol_net[[1]], compDisMat = Mol_net[[2]]$fingerDisMat) #calculates beta diversity
calcDivProf(sampleData = Mol_net[[1]], compDisMat = Mol_net[[2]]$fingerDisMat) # are used to calculate phytochemical diversity in different ways, using both traditional indices and Hill numbers. 
# lookup what q is
Bray_cut_dissim <- sampDis(sampleData = Mol_net[[1]], compDisMat = Mol_net[[2]])

#ChemoDivPlot
Div_plot <- calcDiv(sampleData = Mol_net[[1]], compDisMat = Mol_net[[2]]) 
DivProf <- calcDivProf(sampleData = Mol_net[[1]], 
                       compDisMat =  Mol_net[[2]]$fingerDisMat) # are used to calculate phytochemical diversity in different ways, using both traditional indices and Hill numbers. 

Chem_Div_plot <- chemoDivPlot(
  compDisMat = Mol_net[[2]]$fingerDisMat,
  divData = Div_plot,
  divProfData = DivProf,
  sampDisMat = Bray_curt_dissim$BrayCurtis,
  groupData = Sampdata[[1]]$Sample
)
gridExtra::grid.arrange(Chem_Div_plot)

