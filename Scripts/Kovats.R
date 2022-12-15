  MS_data$CAS.INCHI
#Kovats list
Kovats_list <- EAD_data %>% 
  dplyr::select(-Closest.published.kovats, -`Published.kovats(closest)`) %>% 
  left_join(MS_data %>% dplyr::select(CAS.INCHI, `Published.kovats(closest)`)) %>%
  mutate(GC.MS.KOVATS = suppressWarnings(as.numeric(GC.MS.KOVATS)),
         Closest.published.kovats = as.numeric(`Published.kovats(closest)`)) %>%  
  group_by(Compound.Name,CAS.INCHI, Sample_name, Species, GC.MS.KOVATS) %>%
  summarise(Mean_EAD_Kovats = mean(as.numeric(GC.EAD.KOVATS, na.rm = T)),
            Mean_closest = mean(Closest.published.kovats, na.rm = T)) %>% 
  ungroup() %>%
  pivot_wider(names_from = Species, values_from = Mean_EAD_Kovats) %>% 
  mutate(Order = if_else(GC.MS.KOVATS == 0, `Delia Antiqua`, GC.MS.KOVATS)) %>% 
  arrange(Order) %>% 
  mutate(Compound.Name = fct_inorder(Compound.Name)) %>% 
  dplyr::select(-Order)

print(xtable(Kovats_list), include.rownames = F)

#Kovats_list_EAD %>%  write.xlsx("Output/Vegetables_Diptera/Kovats.list.xlsx")

#xtable(Kovats_list_EAD)
# 
# Kovats_list_MS <- MS_data %>% mutate(as.numeric(`Published.kovats(closest)`)) %>% 
#   group_by(Compound.Name, `Sample descriptor`, `GC-MS.KOVATS`, (`Published.kovats(closest)`)) %>%
#   summarise(Mean_MS_Kovats = mean(as.numeric(`GC-MS.KOVATS`))) %>% 
#   ungroup() %>% 
#   pivot_longer(names_to = "Variable", values_to = "Kovats", 4:5)
# 
# 
# Kovats_TEMP <- MS_data %>% dplyr::select(Compound.Name, CAS.INCHI, `Sample descriptor`, `GC-MS.KOVATS`)  %>% group_by(Compound.Name, `Sample descriptor`) %>% 
#   mutate(Number =  n( )) %>% filter(Number < 2) %>% 
#   group_by(Compound.Name) %>% 
#   mutate(Mean = mean(`GC-MS.KOVATS`)) %>% 
#   ungroup() %>% 
#   arrange(Mean)
