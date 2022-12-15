#Normalizes data 
# New file normalises and filters data 
#Normalises the data from the different sets
Norm_cheminfo <- All_data %>% left_join(Key_compounds_class_data) %>% #left joins the data with the key
  mutate(Compound.Name = if_else(is.na(Checked_name)|Checked_name == "", Compound.Name, Checked_name)) %>% #mutates compound name to the checked name if there is a checked name to overwrite the original name
  mutate(CAS.INCHI = if_else(is.na(Corrected_CAS)|Corrected_CAS == "", CAS.INCHI, Corrected_CAS)) %>% #same as last row but with CAS
  mutate(Functional_class = as.factor(Functional_class)) #makes the functional class to a factor

#Check <- All_data %>% filter(!is.na(Compound.Name)) %>% group_by(Compound.Name,Data, CAS.INCHI) %>% summarise() 

#Filters out the data that comes from EAD
EAD_data <- Norm_cheminfo %>% filter(Data == "EAD") %>% #filters the data to only be data labelelled "EAD"
  mutate(Sample_name = fct_relevel(as.factor(Sample_name),
                                   c("Yellow onion", "Napa", "Radish", "Carrot"#, "Fingerprint Blend 1", "Fingerprint Blend 2", "Fingerprint Blend 3")
                                     ))) #makes the different samples into a string

# Organises EAD data, compound according to number of time detected

EAD_host_data_order <- EAD_data %>% group_by(Compound.Name, Species, Sample_name) %>%  
  summarise(Normalized_response = sum(Normalized_response, na.rm = T)) %>% 
  ungroup() %>% mutate(N = 1) %>% 
  group_by(Compound.Name) %>% 
  summarise(Number = sum(N[Normalized_response > 0]), 
            Sum_strength = sum(Normalized_response, na.rm = T)) %>% ungroup() %>% 
  arrange(-Number, -Sum_strength) %>% mutate(Compound.Name = fct_inorder(Compound.Name))

#Joins the EAD data with the compound/class key
EAD_host_data_functional <- EAD_data %>% group_by(Compound.Name) %>% 
  summarise() %>% ungroup() %>% 
  left_join(Key_compounds_class_data) %>% 
  mutate(Compound.Name = fct_rev(fct_relevel(as.factor(Compound.Name), levels(EAD_host_data_order$Compound.Name))))


# Creates a plot for the different functional classes in the key with colours
Functional_class_plot <- EAD_host_data_functional %>%
  ggplot() +
  geom_tile(aes(x = 1, y = Compound.Name, fill = Functional_class)) +
  theme(legend.position = "right", axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.title = element_blank())
# #ggsave(plot = Functional_class_plot, filename = "Output/Vegetables_Diptera/Functional_classes.jpg",
#        units = "mm", height = 300, width = 300, dpi = 300, bg = "white")

Legend_functional <- ggpubr::get_legend(Functional_class_plot)

Functional_class_plot_no_leg <- Functional_class_plot + guides(fill = "none")
#EAD responses heatmap changing EAD_data to be joined with key messed up the plot somehow
EAD_heatmap <- Chosing_Heatmap(EAD_data, #chooses the data set
                               Heatmap_type_In = "Ave", #takes the data from the different replicates and averages them
                               Order = levels(EAD_host_data_order$Compound.Name), #not sure
                               Xvar = "Species") + # makes it so the X-axis is the species instead of stating the replicate names
  theme(axis.text.x = element_text(face = "italic"),legend.position = "right", axis.ticks.y = element_blank(), #makes the test on the X-axis italic and removes the text on the Y-axis
        axis.text.y = element_blank())  + 
  scale_fill_gradientn("Normalised\nresponse", 
                       colours = c("black", "#a6cee3", "#fb9a99","#e31a1c"), 
                       values = rescale(c(0,1,2,4)))
Legend_heatmap <- ggpubr::get_legend(EAD_heatmap)

EAD_heatmap_no_leg <- EAD_heatmap + guides(fill = "none")

library(grid)
library(gridExtra)
Combined_plot <- grid.arrange(egg::ggarrange(Functional_class_plot_no_leg, EAD_heatmap_no_leg, nrow = 1, widths = c(1,5)),
             Legend_functional, Legend_heatmap,
             layout_matrix = rbind(c(1,1,1,1,1,1,2,2),
                                   c(1,1,1,1,1,1,2,2),
                                   c(1,1,1,1,1,1,3,NA),
                                   c(1,1,1,1,1,1,NA,NA)))
 ggsave(plot = Combined_plot, filename = "Output/Vegetables_Diptera/Heatmap_All_temp.jpg", 
      units = "mm", height = 250, width = 260, dpi = 300, bg = "white") #simply saves the graph as a jpeg in the given format

#Filters out the data that comes from MS
