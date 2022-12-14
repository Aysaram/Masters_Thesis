
#Checks for a EAD response greater than 0 and puts a label Yes or No depending on the value
Gives_response <- Norm_cheminfo %>% group_by(Species, Sample_name, Compound.Name) %>% 
  summarise(Is_ead = sum(Normalized_response > 0)) %>%  ungroup() %>% 
  mutate(Is_ead1 = if_else(Is_ead > 0, "Yes", "No"))

#Checks the data to see if there is a "Yes" from the EAD data and and joins it
Ms_EAD_response <- Norm_cheminfo %>% ungroup() %>% 
  left_join(Gives_response %>% filter(Is_ead1 == "Yes")) %>%  
  left_join(Key_compounds_class_data) %>%
  mutate(Sample_name = if_else(is.na(`Sample descriptor`), Sample_name, `Sample descriptor`)) %>%  
  mutate(One = 1, Area.from.MS = as.numeric(Area.from.MS))

Ms_response <- Ms_EAD_response %>% filter(Data == "MS") %>% 
  group_by(Sample_name) %>% 
  summarise(Sum_area = sum(Area.from.MS, na.rm = T),
         Number_of_comp = sum(One[Area.from.MS > 0], na.rm = T)) %>% ungroup()
  
EAD_response <- Ms_EAD_response %>% filter(Data != "MS") %>% 
  group_by(Sample_name, Species, Compound.Name, CAS.INCHI) %>% 
  summarise() %>% ungroup() %>% left_join(Ms_EAD_response %>% filter(Data == "MS") %>% 
                                            dplyr::select(-Species)) %>% 
  group_by(Sample_name, Species) %>% 
  summarise(EAD_area = sum(Area.from.MS, na.rm = T),
            Number_of_EAD_comp = sum(One[Area.from.MS > 0], na.rm = T)) %>% ungroup() 

Key_values <- EAD_response %>% left_join(Ms_response) %>% 
  mutate(EAD_area = if_else(Species == "Drosophila Melanogaster" & Sample_name == "Yellow onion", 0, EAD_area)) %>% 
  mutate(Ratio_area_resp = EAD_area/Sum_area,
         Ratio_comp_resp = Number_of_EAD_comp/Number_of_comp) 

#Doesn't work but it supposed to take the area from MS and summarises how big of a percentage in total gives a response, I guess that having more sample than one breaks it
# Key_numbers <- Ms_response %>% mutate(One = 1) %>%
#   mutate(Sample_name = if_else(is.na(`Sample descriptor`), Sample_name, `Sample descriptor`)) %>% 
#   group_by(Sample_name) %>% 
#   mutate(Area.from.MS = as.numeric(Area.from.MS), Sum_area = sum(Area.from.MS, na.rm = T),
#          Number_of_comp = sum(One[Area.from.MS > 0], na.rm = T)) %>% ungroup() %>% 
#   group_by(Sample_name, Species, Is_ead) %>% 
#   summarise(Perc_of_total_gives_resp = sum(Area.from.MS)/mean(Sum_area),
#             How_many = n(), Total = mean(Number_of_comp)) %>% ungroup() %>% 
#   mutate(Ratio = How_many / Total)
# 
# 
# List_HS <- MS_data %>%  filter(!(`Sample descriptor` %in% c("Fingerprint Blend 1", "Fingerprint Blend 2", "Fingerprint Blend 3"))) %>% 
#   group_by(Compound.Name, Functional_class, `Sample descriptor`) %>% summarise()
# 
# Overlap_Compounds_Samples <- List_HS %>% filter(!is.na(Compound.Name)) %>% group_by(Compound.Name) %>%  mutate(Number = n()) 
