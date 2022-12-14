# Packages to be loaded.
# For chemodiv additional steps might be needed 
# install.packages("webchem")
library("tidyverse")
library("openxlsx")
library("xtable")
library("webchem")
library("stringi")
library("egg")
library("scales")
source("Scripts/Function_read_in_EAD_v1.R")
source("Scripts/Function_read_in_MS_v1.R")
source("Scripts/Function_normalization_v1.R")
source("Scripts/Heatmap_functions_v2.R")

# Reads in experimental data, data in is described. Also normalizes data across replicate with GC_EAD
# Returns key compounds as well as all Data. 
source("Scripts/Vegetables_dipteran/Read_in.R")

# EAD_data these steps needs to be checked.
source("Scripts/Vegetables_dipteran/EAD_data.R")

# Some extra key values 
source("Scripts/Vegetables_dipteran/Key_values.R")

MS_data <- Norm_cheminfo %>% filter(Data == "MS")

# Gives and prints list of compounds from all of the data
# Key_to_extract <- Norm_cheminfo %>%
#   group_by(Compound.Name, CAS.INCHI) %>%
#   summarise() %>% ungroup()
# write.xlsx(Key_to_extract, file = "Output/Vegetables_Diptera/Simple_Compound_Key.xlsx", sheetName = "Compounds_Class", append = FALSE)

source("Scripts/Vegetables_dipteran/Kovats.R")

# Runs chemodiv
source("Scripts/Vegetables_dipteran/Chemodiv_v2.R")
Molnet_plot_large %>% 
  ggsave(plot = ., filename = "Output/Vegetables_diptera/Molnet_plot_large.jpg", units = "mm", height = 240, width = 280, dpi = 350)
Molnet_plot_large %>% 
  ggsave(plot = ., filename = "Output/Vegetables_diptera/Molnet_plot_large.pdf", units = "mm", height = 240, width = 280)

# Plots heatmap 
source("Scripts/Vegetables_dipteran/MS_heatmap.R")
Dendro_heatmap_leg %>% ggsave(plot = ., filename = "Output/Vegetables_diptera/Heatmap_dendro.jpg",
                              height = 350, width = 210, units = "mm", dpi = 450)
Dendro_heatmap_leg %>% ggsave(plot = ., filename = "Output/Vegetables_diptera/Heatmap_dendro.pdf",
                              height = 350, width = 210, units = "mm")
