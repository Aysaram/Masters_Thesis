# Reads in data 
# Function
rEAD_data <- function(Location, Volatilome_name) {
  Template_EAD_data <- rEADin(Location, Volatilome_name)
  Template_MS_data <- rMSin(Location,  Volatilome_name)
  Template_EAD_data_normalized <- Template_EAD_data %>% NORmalization()
  All_data <- Template_MS_data 
  Template_EAD_data_normalized <- Template_EAD_data_normalized %>% mutate(Data = "EAD") %>%
    mutate_all(as.character)
  # Fix issues with naming on volatilome data sheet mutates all to characters
  Template_MS_data <- Template_MS_data %>%  mutate(Data = "MS") %>% 
    dplyr::rename('CAS.INCHI' = 'CAS/INCHI') %>% 
    mutate_all(as.character)
  # Bind MS data together with the normalized GC-EAD dataset
  Template_EAD_MS_data <- Template_MS_data %>% bind_rows(Template_EAD_data_normalized)
  
  Template_EAD_MS_data %>% return()
}

#Reads in onion data
Onion_data <- rEAD_data("Data/Vegetables_Diptera/Onion_Analysis.xlsx", 
                        Volatilome_name = "Volatilome_220316_Onion_Mix")

#Reads in napa data
Napa_data <- rEAD_data("Data/Vegetables_Diptera/Napa_Analysis.xlsx", 
                       Volatilome_name = "Volatilome_220316_Napa_Mix")

#Reads in radish data
Radish_data <- rEAD_data("Data/Vegetables_Diptera/Radish_Analysis.xlsx", 
                         Volatilome_name = "Volatilome_220328_Radish_Mix")

#Reads in carrot data
Carrot_data <- rEAD_data("Data/Vegetables_Diptera/Carrot_Analysis.xlsx", 
                         Volatilome_name = "Volatilome_220328_Carrot_Mix")

# Droso1_data <- rEAD_data("Data/Vegetables_Diptera/DrosoMix1_Analysis.xlsx",
#                          Volatilome_name = "Volatilome_DrosoMix1")
# 
# Droso2_data <- rEAD_data("Data/Vegetables_Diptera/DrosoMix2_Analysis.xlsx", 
#                          Volatilome_name = "Volatilome_DrosoMix2")
# 
# Droso3_data <- rEAD_data("Data/Vegetables_Diptera/DrosoMix3_Analysis.xlsx", 
#                          Volatilome_name = "Volatilome_DrosoMix3")

#Reads in the key with functional classification 
Key_compounds_class_data <- read.xlsx("Data/Vegetables_Diptera/Key_Compound_names/Key_compounds_class.xlsx") %>% 
  mutate(CAS.INCHI = stri_replace_all_fixed(CAS.INCHI, " ", ""))

#Combines data into one tibble, commented out part messes with Gives_response
All_data <- bind_rows(Onion_data, Carrot_data, Radish_data, Napa_data#, Droso1_data, Droso2_data, Droso3_data)
                      ) %>% mutate(Normalized_response = if_else(is.na(Normalized_response), 0, as.numeric(Normalized_response))) %>% 
  mutate(CAS.INCHI = stri_replace_all_fixed(CAS.INCHI, " ", ""))

remove(Onion_data, Napa_data, Radish_data, Carrot_data #, Droso1_data, Droso2_data, Droso3_data
       )
