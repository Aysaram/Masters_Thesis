library(tidyverse)
library(scales)

Chosing_Heatmap <- function(Normalized_data_in = EAD_host_data, 
                            Heatmap_type_In = "Ave", # Can either be 'Full' or 'Ave
                            Naming = "Compound.Name",
                            Facet = "Sample_name", 
                            Xvar = c("Species", "Sex","Recording_point", "Organ", "Variable_1"), 
                            Order = "Yes",
                            Color_scaling =  "Eco letters") {
  
  Color_scale <- switch(Color_scaling,
                        "Eco letters" =  c("black", "#a6cee3", "#fb9a99","#e31a1c"),
                        "Fire (cont)" = c("black", "#ffffb2", "#fed976", "#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#b10026"),
                        "Water (cont)" = c("#ffffcc", "#c7e9b4", "#7fcdbb", "#41b6c4", "#1d91c0", "#225ea8", "#0c2c84"),
                        "Eart (cont)" = c("#ffffd4", "#fee391", "#fec44f", "#fe9929", "#ec7014", "#cc4c02", "#8c2d04"),
                        "Soft (qual)" = c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69"), 
                        "Sharp (qual)" = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628"), 
                        "Deep (qual)" = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d"))
  
  var <- rlang::sym(Naming)
  Facet3 <- rlang::sym(Facet)
  Xvar <- rlang::syms(Xvar)
  
  Full_data <- Normalized_data_in  %>%   
    group_by(Sample_name, Sex, Species, Recording_point,!!var, Organ, Variable_1, Replicate_number, Normalized_response) %>% 
    summarise() %>% ungroup() %>%  
    spread(key = !!var, value = Normalized_response) %>% 
    mutate_if(is.numeric, funs(ifelse(is.na(.), 0, .)))
  Full_data_without_na <- Full_data %>% 
    gather(key = !!var, value = Normalized_response, 
           8:length(Full_data)) %>% 
    mutate(Xaxis = paste(!!!Xvar))
  
  Averaged_data <- Full_data_without_na %>% 
    mutate(Number_of_reps = max(as.numeric(Replicate_number))) %>% # Calculates averages for number of data points 
    ungroup() %>% 
    group_by(!!!Xvar, !!var, !!Facet3) %>% 
    summarise(Normalized_response = sum(Normalized_response)/Number_of_reps) %>% 
    ungroup() %>% 
    mutate(Xaxis = paste(!!!Xvar),
           Variable = !!var) %>% 
    mutate(Variable  = fct_relevel(Variable, rev(Order),after = Inf))
  
  Full_data_plot <- Full_data_without_na %>%  
    group_by(!!!Xvar, !!var, !!Facet3, Replicate_number) %>% 
    summarise(Normalized_response = mean(Normalized_response)) %>% 
    ungroup() %>% 
    mutate(Xaxis = paste(!!!Xvar, Replicate_number),
           Variable = !!var) %>% 
    mutate(Variable  = fct_relevel(Variable, rev(Order),after = Inf)) %>% 
    mutate(Variable = tolower(Variable)) %>% filter(!is.na(Species))
  # 
  # Order <- Averaged_data %>% filter(!is.na(Normalized_response)) %>% group_by( Variable) %>% 
  #   summarise(Number = n(), Response = sum(Normalized_response)) %>% ungroup() %>% 
  #   arrange(Response) %>% arrange(Number) %>%  
  #   mutate(Variable= fct_inorder(Variable))
  
  # Averaged_data <- switch(Order,
  #        "No" = Averaged_data)
  
  
  # Full_data_without_na$Compound.Name <- Full_data_without_na$Compound.Name  %>% fct_relevel(levels(Normalized_data_in$Compound.Name))
  Heatmap_function <- function(Data_In,  Facet2 = Facet3) {
    
    Data_In %>%
      ggplot(aes(x = Xaxis, y = Variable)) +
      geom_tile(aes(fill = Normalized_response)) + 
      scale_fill_gradientn("Normalized\nresponse", 
                           colours = Color_scale, 
                           values = rescale(c(0,1,2,4))) +
      facet_grid(cols = vars(!!Facet2),  scales = "free_x", space = "free_x") + # What variable to call? Should it be dynamic or fixed?
      theme(strip.placement = "outside",
            axis.text.x = element_text(angle = 90,hjust = 1, vjust = .5),
            axis.title.y = element_blank(), axis.title.x = element_blank())%>%  
      return()
    
  }
  
  Heatplot <- switch(Heatmap_type_In,
                     "Full" = Heatmap_function(Full_data_plot), 
                     "Ave" =  Heatmap_function(Averaged_data,Facet3)) 
  
  Heatplot %>% return()
}

