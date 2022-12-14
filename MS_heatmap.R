# Heatmap to run after Chemodiv
compDisMatClust <- stats::hclust(stats::as.dist(Mol_net[[2]]$fingerDisMat), 
                                 method = "average")
compDisMatClustDend <- stats::as.dendrogram(compDisMatClust)

GGtree_plot <- ggtree::ggtree(compDisMatClustDend) +
  ggtree::geom_tiplab() +
  hexpand(1)

# Plot data
From_plot <- GGtree_plot$data %>% 
  filter(isTip == T) %>% 
  arrange(y)
GGtree_plot2 <- ggtree::ggtree(compDisMatClustDend) + 
  scale_x_reverse()

Ordered <- MS_data  %>% 
  filter(!(`Sample descriptor` %in% c("Fingerprint Blend 1", "Fingerprint Blend 2", "Fingerprint Blend 3"))) %>% 
  mutate(Sample = if_else(is.na(`Sample descriptor`), Sample_name, `Sample descriptor`),
         Area.from.MS = as.numeric(Area.from.MS)) %>%
  filter(!is.na(Compound.Name), !is.na(Area.from.MS), !is.na(Smiles),
         !is.na(ST_InChiKey)) %>%  
  mutate(Compound.Name = as.factor(Compound.Name)) %>% 
  mutate(Compound.Name = fct_relevel(Compound.Name, From_plot$label)) 
Functional <- Ordered %>% left_join(NP_classified, by = c(Compound.Name = "compound"))
Color_scale_pathway <- c("#E7298A", "#66A61E", "#7570B3",  "#A6761D", "#1B9E77", "#D95F02")
Order_pathway <- c("Alkaloids", "Amino acids and Peptides","Carbohydrates","Fatty acids",
                   "Polyketides", "Shikimates and Phenylpropanoids","Terpenoids")
Functions_class_heatmap  <- Functional %>% mutate(pathway = fct_relevel(as.factor(pathway), Order_pathway)) %>% 
  mutate(Compound.Name = fct_relevel(Compound.Name, From_plot$label)) %>%  
  ggplot(aes(1, Compound.Name)) +
  geom_tile(aes(fill = pathway)) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) +
  scale_fill_manual("Pathway", values = Color_scale_pathway) +
  theme(legend.position = "bottom") +
  guides(fill =  guide_legend(title.position="top", title.hjust = 0.5, ncol = 1))  +
  ylab("Compound name")


Heatmap <- Ordered %>%  mutate(Sample = fct_relevel(as.factor(Sample),
                                                    c("Yellow onion", "Carrot", "Napa", "Radish"))) %>% 
  ggplot(aes(Sample, Compound.Name)) +
  geom_tile(aes(fill = as.numeric(Area.from.MS))) +
  theme(legend.position = "bottom") +
  scale_fill_continuous("Area from MS", breaks = c(1e+07, 5e+07)) +
  guides(fill =  guide_colourbar(title.position="top", title.hjust = 0.5)) +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
        axis.ticks.y = element_blank(), axis.title.x = element_blank())

Legends <- gridExtra::grid.arrange(ggpubr::get_legend(Functions_class_heatmap, position = NULL),
                                   ggpubr::get_legend(Heatmap, position = NULL), nrow = 1)


Dendro_heatmap <- egg::ggarrange(Functions_class_heatmap %>% Remove_guides(), 
                                 Heatmap %>% Remove_guides(), GGtree_plot2, ncol = 3, widths = c(.5,3,2))
Dendro_heatmap_leg <- gridExtra::grid.arrange(Dendro_heatmap,Legends, ncol = 1,  layout_matrix = rbind(c(1,1,1,1),
                                                                                                       c(1,1,1,1),
                                                                                                       c(1,1,1,1),
                                                                                                       c(1,1,1,1),
                                                                                                       c(NA,2,2,NA)))
