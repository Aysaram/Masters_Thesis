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
source("Scripts/Functions/Function_read_in_EAD_v1.R")
source("Scripts/Functions/Function_read_in_MS_v1.R")
source("Scripts/Functions/Function_normalization_v1.R")
source("Scripts/Functions/Heatmap_functions_v2.R")

# Reads in experimental data, data in is described. Also normalizes data across replicate with GC_EAD
# Returns key compounds as well as all Data. 
source("Scripts/Read_in.R")

# EAD_data these steps needs to be checked.
source("Scripts/EAD_data.R")
ggsave(plot = Combined_plot, filename = "Output/Fig_4.jpg", 
       units = "mm", height = 250, width = 260, dpi = 300, bg = "white") #simply saves the graph as a jpeg in the given format
ggsave(plot = Combined_plot, filename = "Output/Fig_4.pdf", 
       units = "mm", height = 250, width = 260, dpi = 300, bg = "white") #simply saves the graph as a jpeg in the given format

# Some extra key values 
source("Scripts/Key_values.R")

MS_data <- Norm_cheminfo %>% filter(Data == "MS")

# Gives and prints list of compounds from all of the data
# Key_to_extract <- Norm_cheminfo %>%
#   group_by(Compound.Name, CAS.INCHI) %>%
#   summarise() %>% ungroup()
# write.xlsx(Key_to_extract, file = "Output/Vegetables_Diptera/Simple_Compound_Key.xlsx", sheetName = "Compounds_Class", append = FALSE)

source("Scripts/Kovats.R")

# Runs chemodiv
source("Scripts/Chemodiv_v2.R")
Molnet_plot_large %>% 
  ggsave(plot = ., filename = "Output/Fig_3.jpg", units = "mm", height = 240, width = 280, dpi = 350)
Molnet_plot_large %>% 
  ggsave(plot = ., filename = "Output/Fig_3.pdf", units = "mm", height = 240, width = 280)

# Plots heatmap 
source("Scripts/MS_heatmap.R")
Dendro_heatmap_leg %>% ggsave(plot = ., filename = "Output/Fig_2.jpg",
                              height = 350, width = 210, units = "mm", dpi = 450)
Dendro_heatmap_leg %>% ggsave(plot = ., filename = "Output/Fig_2.pdf",
                              height = 350, width = 210, units = "mm")
