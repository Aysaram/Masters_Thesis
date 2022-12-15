
molNetPlot <-   function(sampleData,
                       networkObject,
                       groupData = NULL,
                       npcTable = NULL,
                       plotNames = FALSE,
                       layout = "kk")
                       {
                         if (!is.null(groupData) && plotNames) {
                           stop("Names can only be plotted without grouping data.")
                         }
                         if (is.null(npcTable)) {
                           message("It is recommended to include an npcTable for an improved\n            network visualization.")
                         }
                         if (is.data.frame(groupData)) {
                           groupData <- as.vector(groupData[, 1])
                         }
                         if (is.null(groupData) && is.null(npcTable) && !plotNames) {
                           compoundMean <- colMeans(sampleData)
                           p1 <- ggraph::ggraph(graph = networkObject, layout = layout) + 
                             ggraph::geom_edge_link(ggplot2::aes(width = .data$weight), 
                                                    edge_color = "grey40") + ggraph::scale_edge_width(range = c(0.3, 
                                                                                                                3), name = "Molecular similarity") + ggraph::geom_node_point(ggplot2::aes(color = compoundMean), 
                                                                                                                                                                             size = 16) + ggplot2::scale_colour_viridis_c() + 
                             ggplot2::labs(color = "Proportion", width = "Molecular similarity") + 
                             ggplot2::theme(legend.title = ggplot2::element_text(size = 16), 
                                            legend.text = ggplot2::element_text(size = 14), 
                                            panel.background = ggplot2::element_blank(), 
                                            legend.key = ggplot2::element_blank())
                           networkList <- list(p1)
                         }
                         else if (is.null(groupData) && is.null(npcTable) && plotNames) {
                           compoundMean <- colMeans(sampleData)
                           p1 <- ggraph::ggraph(graph = networkObject, layout = layout) + 
                             ggraph::geom_edge_link(ggplot2::aes(width = .data$weight), 
                                                    edge_color = "grey40") + ggraph::scale_edge_width(range = c(0.3, 
                                                                                                                3), name = "Molecular similarity") + ggraph::geom_node_point(ggplot2::aes(color = compoundMean), 
                                                                                                                                                                             size = 16) + ggplot2::scale_colour_viridis_c() + 
                             ggplot2::labs(color = "Proportion", width = "Molecular similarity") + 
                             ggraph::geom_node_label(ggplot2::aes(label = .data$name), 
                                                     nudge_x = 0, nudge_y = 0.2) + ggplot2::theme(legend.title = ggplot2::element_text(size = 16), 
                                                                                                  legend.text = ggplot2::element_text(size = 14), panel.background = ggplot2::element_blank(), 
                                                                                                  legend.key = ggplot2::element_blank())
                           networkList <- list(p1)
                         }
                         else if (is.null(groupData) && !is.null(npcTable) && !plotNames) {
                           npcTable$col <- NA
                           npcTable$col[is.na(npcTable$pathway)] <- "#666666"
                           npcTable$col[npcTable$pathway == "Alkaloids"] <- "#E7298A"
                           npcTable$col[npcTable$pathway == "Amino acids and Peptides"] <- "#66A61E"
                           npcTable$col[npcTable$pathway == "Carbohydrates"] <- "#F0E442"
                           npcTable$col[npcTable$pathway == "Fatty acids"] <- "#7570B3"
                           npcTable$col[npcTable$pathway == "Polyketides"] <- "#A6761D"
                           npcTable$col[npcTable$pathway == "Shikimates and Phenylpropanoids"] <- "#1B9E77"
                           npcTable$col[npcTable$pathway == "Terpenoids"] <- "#D95F02"
                           incPath <- unique(npcTable$pathway)[!is.na(unique(npcTable$pathway))]
                           incPathOrder <- incPath[order(incPath)]
                           legCol <- npcTable$col[match(incPathOrder, npcTable$pathway)]
                           compoundMean <- colMeans(sampleData)
                           p1 <- ggraph::ggraph(graph = networkObject, layout = layout) + 
                             ggraph::geom_edge_link(ggplot2::aes(width = .data$weight), 
                                                    edge_color = "grey40") + ggraph::scale_edge_width(range = c(0.3, 
                                                                                                                3), name = "Molecular similarity") + ggraph::geom_node_point(ggplot2::aes(color = npcTable$pathway, 
                                                                                                                                                                                          fill = npcTable$pathway, size = compoundMean), shape = 21) + 
                             ggraph::geom_node_point(ggplot2::aes(size = compoundMean), 
                                                     shape = 21, color = npcTable$col, fill = npcTable$col) + 
                             ggplot2::scale_fill_manual(values = legCol) + ggplot2::scale_size(range = c(8, 
                                                                                                         16)) + ggplot2::labs(fill = "Pathway", width = "Molecular similarity", 
                                                                                                                              size = "Proportion") + ggplot2::guides(color = "none") + 
                             ggplot2::theme(legend.title = ggplot2::element_text(size = 16), 
                                            legend.text = ggplot2::element_text(size = 14), 
                                            panel.background = ggplot2::element_blank(), 
                                            legend.key = ggplot2::element_blank())
                           networkList <- list(p1)
                         }
                         else if (is.null(groupData) && !is.null(npcTable) && plotNames) {
                           npcTable$col <- NA
                           npcTable$col[is.na(npcTable$pathway)] <- "#666666"
                           npcTable$col[npcTable$pathway == "Alkaloids"] <- "#E7298A"
                           npcTable$col[npcTable$pathway == "Amino acids and Peptides"] <- "#66A61E"
                           npcTable$col[npcTable$pathway == "Carbohydrates"] <- "#F0E442"
                           npcTable$col[npcTable$pathway == "Fatty acids"] <- "#7570B3"
                           npcTable$col[npcTable$pathway == "Polyketides"] <- "#A6761D"
                           npcTable$col[npcTable$pathway == "Shikimates and Phenylpropanoids"] <- "#1B9E77"
                           npcTable$col[npcTable$pathway == "Terpenoids"] <- "#D95F02"
                           incPath <- unique(npcTable$pathway)[!is.na(unique(npcTable$pathway))]
                           incPathOrder <- incPath[order(incPath)]
                           legCol <- npcTable$col[match(incPathOrder, npcTable$pathway)]
                           compoundMean <- colMeans(sampleData)
                           p1 <- ggraph::ggraph(graph = networkObject, layout = layout) + 
                             ggraph::geom_edge_link(ggplot2::aes(width = .data$weight), 
                                                    edge_color = "grey40") + ggraph::scale_edge_width(range = c(0.3, 
                                                                                                                3), name = "Molecular similarity") + ggraph::geom_node_point(ggplot2::aes(color = npcTable$pathway, 
                                                                                                                                                                                          fill = npcTable$pathway, size = compoundMean), shape = 21) + 
                             ggraph::geom_node_point(ggplot2::aes(size = compoundMean), 
                                                     shape = 21, color = npcTable$col, fill = npcTable$col) + 
                             ggplot2::scale_fill_manual(values = legCol) + ggplot2::scale_size(range = c(8, 
                                                                                                         16)) + ggplot2::guides(color = "none") + ggplot2::labs(fill = "Pathway", 
                                                                                                                                                                width = "Molecular similarity", size = "Proportion") + 
                             ggraph::geom_node_label(ggplot2::aes(label = .data$name), 
                                                     nudge_x = 0, nudge_y = 0.2) + ggplot2::theme(legend.title = ggplot2::element_text(size = 16), 
                                                                                                  legend.text = ggplot2::element_text(size = 14), panel.background = ggplot2::element_blank(), 
                                                                                                  legend.key = ggplot2::element_blank())
                           networkList <- list(p1)
                         }
                         else if (!is.null(groupData) && !is.null(npcTable)) {
                           npcTable$col <- NA
                           npcTable$col[is.na(npcTable$pathway)] <- "Unknown"
                           npcTable$col[npcTable$pathway == "Alkaloids"] <- "Alkaloids"
                           npcTable$col[npcTable$pathway == "Amino acids and Peptides"] <- "Amino acids and Peptides"
                           npcTable$col[npcTable$pathway == "Carbohydrates"] <- "Carbohydrates"
                           npcTable$col[npcTable$pathway == "Fatty acids"] <- "Fatty acids"
                           npcTable$col[npcTable$pathway == "Polyketides"] <- "Polyketides"
                           npcTable$col[npcTable$pathway == "Shikimates and Phenylpropanoids"] <- "Shikimates and Phenylpropanoids"
                           npcTable$col[npcTable$pathway == "Terpenoids"] <- "Terpenoids"
                           incPath <- unique(npcTable$pathway)[!is.na(unique(npcTable$pathway))]
                           incPathOrder <- incPath[order(incPath)]
                           legCol <- npcTable$col[match(incPathOrder, npcTable$pathway)]
                           compoundMean <- stats::aggregate(sampleData, by = list(Group = groupData), 
                                                            mean)
                           compoundMeanTrans <- t(compoundMean[, 2:ncol(compoundMean)])
                           colnames(compoundMeanTrans) <- compoundMean$Group
                           compoundMeanTrans <- as.data.frame(compoundMeanTrans)
                           networkList <- list()
                           Breaks <- c("Unknown", "Alkaloids", "Amino acids and Peptides", "Carbohydrates",
                           "Fatty acids", "Polyketides","Shikimates and Phenylpropanoids","Terpenoids")
                           for (j in 1:ncol(compoundMeanTrans)) {
                             networkList[[colnames(compoundMeanTrans)[j]]] <- local({
                               j <- j
                               groupCol <- as.data.frame(npcTable$col)
                               colnames(groupCol) <- "col"
                               for (row in 1:nrow(compoundMeanTrans)) {
                                 if (compoundMeanTrans[row, j] == 0) {
                                   groupCol$col[row] <- "transparent"
                                 }
                               }
                               compSimMat <- matrix(data = NA, nrow = nrow(compoundMeanTrans), 
                                                    ncol = nrow(compoundMeanTrans))
                               colnames(compSimMat) <- rownames(compoundMeanTrans)
                               rownames(compSimMat) <- rownames(compoundMeanTrans)
                               for (i in 1:nrow(compoundMeanTrans)) {
                                 compSimMat[i, ] <- networkObject[i]
                               }
                               linkedComps <- as.data.frame(matrix(data = NA, 
                                                                   nrow = length(compSimMat[compSimMat > 0])/2, 
                                                                   ncol = 2))
                               colnames(linkedComps) <- c("Comp1", "Comp2")
                               l <- 1
                               for (i in 1:nrow(compSimMat)) {
                                 for (k in i:ncol(compSimMat)) {
                                   if (compSimMat[i, k] > 0) {
                                     linkedComps$Comp1[l] <- rownames(compSimMat)[i]
                                     linkedComps$Comp2[l] <- colnames(compSimMat)[k]
                                     l <- l + 1
                                   }
                                 }
                               }
                               linkCol <- rep("grey40", nrow(linkedComps))
                               for (i in 1:nrow(linkedComps)) {
                                 row1 <- match(linkedComps$Comp1[i], rownames(compoundMeanTrans))
                                 row2 <- match(linkedComps$Comp2[i], rownames(compoundMeanTrans))
                                 if (compoundMeanTrans[row1, j] == 0 | compoundMeanTrans[row2, 
                                                                                         j] == 0) {
                                   linkCol[i] <- "transparent"
                                 }
                               }
                               compSimMatTri <- compSimMat[upper.tri(compSimMat)]
                               compSimMatTriPos <- compSimMatTri[compSimMatTri > 
                                                                   0]
                               nodes <- data.frame(name = rownames(compoundMeanTrans))
                               links <- data.frame(from = linkedComps$Comp1, 
                                                   to = linkedComps$Comp2, weight = compSimMatTriPos, 
                                                   linkCol = linkCol)
                               networkObjectMan <- tidygraph::tbl_graph(nodes = nodes, 
                                                                        edges = links)
                               Color_scale <- c('grey80',
                                                '#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628')
                               Breaks
                               Level <- bind_cols(as_tibble(Color_scale), as_tibble(Breaks))
                               colnames(Level) <- c("Color", "Pathway")
                            
                               # Color_scale <-  c(Color_scale[1:(length(unique(npcTable$col)))], "transparent")
                               # c(Color_scale[1:5], "transparent")
                               p1 <- ggraph::ggraph(graph = networkObjectMan, 
                                                    layout = layout) + 
                                 ggraph::geom_edge_link(ggplot2::aes(width = .data$weight,
                                                                     color = .data$linkCol)) +
                                 ggraph::scale_edge_width(range = c(0.3, 2.5), name = "Molecular similarity") +
                                 ggraph::geom_node_point(ggplot2::aes(color = npcTable$pathway, fill = npcTable$pathway,
                                                                      size = compoundMeanTrans[,j]), shape = 21) +
                                 ggraph::scale_edge_color_manual(limits = as.factor(links$linkCol),
                                                                 values = links$linkCol, guide = "none") + 
                                 ggraph::geom_node_point(ggplot2::aes(size = compoundMeanTrans[, j],
                                                                      fill = groupCol$col, color = groupCol$col), 
                                                         shape = 21) +
                                  ggplot2::scale_fill_manual(values = Level$Color, 
                                                             breaks = Level$Pathway,na.value = "transparent"
                                                             ) + 
                                  ggplot2::scale_color_manual(values = c(rep("black", 10)), 
                                 breaks = Level$Pathway,na.value = "transparent") + 
                                 ggplot2::scale_size(range = c(4, 10)) + ggplot2::guides(color = "none") + 
                                 ggplot2::labs(fill = "Pathway", width = "Molecular similarity", 
                                               size = "Proportion") + ggplot2::ggtitle(colnames(compoundMeanTrans)[j]) + 
                                 ggplot2::theme(legend.title = ggplot2::element_text(size = 10), 
                                                legend.text = ggplot2::element_text(size = 8), 
                                                panel.background = ggplot2::element_blank(), 
                                                legend.key = ggplot2::element_blank()) 
                               p1
                             })
                           }
                         }
                         else {
                           compoundMean <- stats::aggregate(sampleData, by = list(Group = groupData), 
                                                            mean)
                           compoundMeanTrans <- t(compoundMean[, 2:ncol(compoundMean)])
                           colnames(compoundMeanTrans) <- compoundMean$Group
                           compoundMeanTrans <- as.data.frame(compoundMeanTrans)
                           networkList <- list()
                           for (j in 1:ncol(compoundMeanTrans)) {
                             networkList[[colnames(compoundMeanTrans)[j]]] <- local({
                               j <- j
                               groupShape <- data.frame(shape = rep("Pos", nrow(compoundMeanTrans)))
                               for (row in 1:nrow(compoundMeanTrans)) {
                                 if (compoundMeanTrans[row, j] == 0) {
                                   groupShape$shape[row] <- "Zero"
                                 }
                               }
                               compSimMat <- matrix(data = NA, nrow = nrow(compoundMeanTrans), 
                                                    ncol = nrow(compoundMeanTrans))
                               colnames(compSimMat) <- rownames(compoundMeanTrans)
                               rownames(compSimMat) <- rownames(compoundMeanTrans)
                               for (i in 1:nrow(compoundMeanTrans)) {
                                 compSimMat[i, ] <- networkObject[i]
                               }
                               linkedComps <- as.data.frame(matrix(data = NA, 
                                                                   nrow = length(compSimMat[compSimMat > 0])/2, 
                                                                   ncol = 2))
                               colnames(linkedComps) <- c("Comp1", "Comp2")
                               l <- 1
                               for (i in 1:nrow(compSimMat)) {
                                 for (k in i:ncol(compSimMat)) {
                                   if (compSimMat[i, k] > 0) {
                                     linkedComps$Comp1[l] <- rownames(compSimMat)[i]
                                     linkedComps$Comp2[l] <- colnames(compSimMat)[k]
                                     l <- l + 1
                                   }
                                 }
                               }
                               linkCol <- rep("transparent", nrow(linkedComps))
                               for (i in 1:nrow(linkedComps)) {
                                 row1 <- match(linkedComps$Comp1[i], rownames(compoundMeanTrans))
                                 row2 <- match(linkedComps$Comp2[i], rownames(compoundMeanTrans))
                                 if (compoundMeanTrans[row1, j] == 0 | compoundMeanTrans[row2, 
                                                                                         j] == 0) {
                                   linkCol[i] <- "grey90"
                                 }
                               }
                               compSimMatTri <- compSimMat[upper.tri(compSimMat)]
                               compSimMatTriPos <- compSimMatTri[compSimMatTri > 
                                                                   0]
                               nodes <- data.frame(name = rownames(compoundMeanTrans))
                               links <- data.frame(from = linkedComps$Comp1, 
                                                   to = linkedComps$Comp2, weight = compSimMatTriPos, 
                                                   linkCol = linkCol)
                               networkObjectMan <- tidygraph::tbl_graph(nodes = nodes, 
                                                                        edges = links)
                               p1 <- ggraph::ggraph(graph = networkObjectMan, 
                                                    layout = layout) + 
                                 ggraph::geom_edge_link(ggplot2::aes(width = .data$weight, 
                                                                     color = .data$linkCol),edge_alpha = .4) + 
                                 ggraph::scale_edge_width(range = c(0.3, 2.5), name = "Molecular similarity", 
                                                            guide = "none" ) + 
                                 ggraph::geom_node_point(ggplot2::aes(color = compoundMeanTrans[, j], 
                                                                      shape = groupShape$shape), 
                                                         guide = "none", size = 10) + 
                                 ggraph::scale_edge_color_manual(limits = as.factor(links$linkCol), 
                                                                 values = links$linkCol, guide = "none") + 
                                 ggplot2::scale_colour_viridis_c() + 
                                 ggplot2::labs(color = "Proportion", 
                                               width = "Molecular similarity") + 
                                 ggplot2::ggtitle(colnames(compoundMeanTrans)[j]) + 
                                 ggplot2::guides(shape = "none") +
                                 ggplot2::theme(legend.title = ggplot2::element_text(size = 10), 
                                                                                  legend.text = ggplot2::element_text(size = 8), 
                                                                                  panel.background = ggplot2::element_blank(), 
                                                                                  legend.key = ggplot2::element_blank())
                           #    print(p1)
                               p1
                             })
                           }
                         }
                         # gridExtra::grid.arrange(grobs = networkList, ncol = ceiling(sqrt(length(networkList))))
                         return(networkList)
}
  

Molnet_plot2 <- function (sampleData, networkObject, groupData = NULL, npcTable = NULL, 
          plotNames = FALSE, layout = "kk") 
{
  if (!is.null(groupData) && plotNames) {
    stop("Names can only be plotted without grouping data.")
  }
  if (is.null(npcTable)) {
    message("It is recommended to include an npcTable for an improved\n            network visualization.")
  }
  if (is.data.frame(groupData)) {
    groupData <- as.vector(groupData[, 1])
  }
  if (is.null(groupData) && is.null(npcTable) && !plotNames) {
    compoundMean <- colMeans(sampleData)
    p1 <- ggraph::ggraph(graph = networkObject, layout = layout) + 
      ggraph::geom_edge_link(ggplot2::aes(width = .data$weight), 
                             edge_color = "grey40") + ggraph::scale_edge_width(range = c(0.3, 
                                                                                         3), name = "Molecular similarity") + ggraph::geom_node_point(ggplot2::aes(color = compoundMean), 
                                                                                                                                                      size = 16) + ggplot2::scale_colour_viridis_c() + 
      ggplot2::labs(color = "Proportion", width = "Molecular similarity") + 
      ggplot2::theme(legend.title = ggplot2::element_text(size = 16), 
                     legend.text = ggplot2::element_text(size = 14), 
                     panel.background = ggplot2::element_blank(), 
                     legend.key = ggplot2::element_blank())
    networkList <- list(p1)
  }
  else if (is.null(groupData) && is.null(npcTable) && plotNames) {
    compoundMean <- colMeans(sampleData)
    p1 <- ggraph::ggraph(graph = networkObject, layout = layout) + 
      ggraph::geom_edge_link(ggplot2::aes(width = .data$weight), 
                             edge_color = "grey40") + ggraph::scale_edge_width(range = c(0.3, 
                                                                                         3), name = "Molecular similarity") + ggraph::geom_node_point(ggplot2::aes(color = compoundMean), 
                                                                                                                                                      size = 16) + ggplot2::scale_colour_viridis_c() + 
      ggplot2::labs(color = "Proportion", width = "Molecular similarity") + 
      ggraph::geom_node_label(ggplot2::aes(label = .data$name), 
                              nudge_x = 0, nudge_y = 0.2) + ggplot2::theme(legend.title = ggplot2::element_text(size = 16), 
                                                                           legend.text = ggplot2::element_text(size = 14), panel.background = ggplot2::element_blank(), 
                                                                           legend.key = ggplot2::element_blank())
    networkList <- list(p1)
  }
  else if (is.null(groupData) && !is.null(npcTable) && !plotNames) {
    npcTable$col <- NA
    npcTable$col[is.na(npcTable$pathway)] <- "#666666"
    npcTable$col[npcTable$pathway == "Alkaloids"] <- "#E7298A"
    npcTable$col[npcTable$pathway == "Amino acids and Peptides"] <- "#66A61E"
    npcTable$col[npcTable$pathway == "Carbohydrates"] <- "#F0E442"
    npcTable$col[npcTable$pathway == "Fatty acids"] <- "#7570B3"
    npcTable$col[npcTable$pathway == "Polyketides"] <- "#A6761D"
    npcTable$col[npcTable$pathway == "Shikimates and \n Phenylpropanoids"] <- "#1B9E77"
    npcTable$col[npcTable$pathway == "Terpenoids"] <- "#D95F02"
    incPath <- unique(npcTable$pathway)[!is.na(unique(npcTable$pathway))]
    incPathOrder <- incPath[order(incPath)]
    legCol <- npcTable$col[match(incPathOrder, npcTable$pathway)]
    compoundMean <- colMeans(sampleData)
    p1 <- ggraph::ggraph(graph = networkObject, layout = layout) + 
      ggraph::geom_edge_link(ggplot2::aes(width = .data$weight), 
                             edge_color = "grey40") + ggraph::scale_edge_width(range = c(0.3, 
                                                                                         3), name = "Molecular similarity") + ggraph::geom_node_point(ggplot2::aes(color = npcTable$pathway, 
                                                                                                                                                                   fill = npcTable$pathway, size = compoundMean), shape = 21) + 
      ggraph::geom_node_point(ggplot2::aes(size = compoundMean), 
                              shape = 21, color = npcTable$col, fill = npcTable$col) + 
      ggplot2::scale_fill_manual(values = legCol) + ggplot2::scale_size(range = c(8, 
                                                                                  16)) + ggplot2::labs(fill = "Pathway", width = "Molecular similarity", 
                                                                                                       size = "Proportion") + ggplot2::guides(color = "none") + 
      ggplot2::theme(legend.title = ggplot2::element_text(size = 16), 
                     legend.text = ggplot2::element_text(size = 14), 
                     panel.background = ggplot2::element_blank(), 
                     legend.key = ggplot2::element_blank())
    networkList <- p1
  }
  else if (is.null(groupData) && !is.null(npcTable) && plotNames) {
    npcTable$col <- NA
    npcTable$col[is.na(npcTable$pathway)] <- "#666666"
    npcTable$col[npcTable$pathway == "Alkaloids"] <- "#E7298A"
    npcTable$col[npcTable$pathway == "Amino acids and Peptides"] <- "#66A61E"
    npcTable$col[npcTable$pathway == "Carbohydrates"] <- "#F0E442"
    npcTable$col[npcTable$pathway == "Fatty acids"] <- "#7570B3"
    npcTable$col[npcTable$pathway == "Polyketides"] <- "#A6761D"
    npcTable$col[npcTable$pathway == "Shikimates and \n Phenylpropanoids"] <- "#1B9E77"
    npcTable$col[npcTable$pathway == "Terpenoids"] <- "#D95F02"
    incPath <- unique(npcTable$pathway)[!is.na(unique(npcTable$pathway))]
    incPathOrder <- incPath[order(incPath)]
    legCol <- npcTable$col[match(incPathOrder, npcTable$pathway)]
    compoundMean <- colMeans(sampleData)
    p1 <- ggraph::ggraph(graph = networkObject, layout = layout) + 
      ggraph::geom_edge_link(ggplot2::aes(width = .data$weight), 
                             edge_color = "grey40") + ggraph::scale_edge_width(range = c(0.3, 
                                                                                         3), name = "Molecular similarity") + ggraph::geom_node_point(ggplot2::aes(color = npcTable$pathway, 
                                                                                                                                                                   fill = npcTable$pathway, size = compoundMean), shape = 21) + 
      ggraph::geom_node_point(ggplot2::aes(size = compoundMean), 
                              shape = 21, color = npcTable$col, fill = npcTable$col) + 
      ggplot2::scale_fill_manual(values = legCol) + ggplot2::scale_size(range = c(8, 
                                                                                  16)) + ggplot2::guides(color = "none") + ggplot2::labs(fill = "Pathway", 
                                                                                                                                         width = "Molecular similarity", size = "Proportion") + 
      ggraph::geom_node_label(ggplot2::aes(label = .data$name), 
                              nudge_x = 0, nudge_y = 0.2) + ggplot2::theme(legend.title = ggplot2::element_text(size = 16), 
                                                                           legend.text = ggplot2::element_text(size = 14), panel.background = ggplot2::element_blank(), 
                                                                           legend.key = ggplot2::element_blank())
    networkList <- p1
    }
  else if (!is.null(groupData) && !is.null(npcTable)) {
    npcTable$col <- NA
    npcTable$col[is.na(npcTable$pathway)] <- "#666666"
    npcTable$col[npcTable$pathway == "Alkaloids"] <- "#E7298A"
    npcTable$col[npcTable$pathway == "Amino acids and Peptides"] <- "#66A61E"
    npcTable$col[npcTable$pathway == "Carbohydrates"] <- "#F0E442"
    npcTable$col[npcTable$pathway == "Fatty acids"] <- "#7570B3"
    npcTable$col[npcTable$pathway == "Polyketides"] <- "#A6761D"
    npcTable$col[npcTable$pathway == "Shikimates and Phenylpropanoids"] <- "#1B9E77"
    npcTable$col[npcTable$pathway == "Terpenoids"] <- "#D95F02"
    incPath <- unique(npcTable$pathway)[!is.na(unique(npcTable$pathway))]
    incPathOrder <- incPath[order(incPath)]
    legCol <- npcTable$col[match(incPathOrder, npcTable$pathway)]
    compoundMean <- stats::aggregate(sampleData, by = list(Group = groupData), 
                                     mean)
    compoundMeanTrans <- t(compoundMean[, 2:ncol(compoundMean)])
    colnames(compoundMeanTrans) <- compoundMean$Group
    compoundMeanTrans <- as.data.frame(compoundMeanTrans)
    networkList <- list()
    for (j in 1:ncol(compoundMeanTrans)) {
      networkList[[colnames(compoundMeanTrans)[j]]] <- local({
        j <- j
        groupCol <- as.data.frame(npcTable$col)
        colnames(groupCol) <- "col"
        for (row in 1:nrow(compoundMeanTrans)) {
          if (compoundMeanTrans[row, j] == 0) {
            groupCol$col[row] <- "#FFFFFF"
          }
        }
        compSimMat <- matrix(data = NA, nrow = nrow(compoundMeanTrans), 
                             ncol = nrow(compoundMeanTrans))
        colnames(compSimMat) <- rownames(compoundMeanTrans)
        rownames(compSimMat) <- rownames(compoundMeanTrans)
        for (i in 1:nrow(compoundMeanTrans)) {
          compSimMat[i, ] <- networkObject[i]
        }
        linkedComps <- as.data.frame(matrix(data = NA, 
                                            nrow = length(compSimMat[compSimMat > 0])/2, 
                                            ncol = 2))
        colnames(linkedComps) <- c("Comp1", "Comp2")
        l <- 1
        for (i in 1:nrow(compSimMat)) {
          for (k in i:ncol(compSimMat)) {
            if (compSimMat[i, k] > 0) {
              linkedComps$Comp1[l] <- rownames(compSimMat)[i]
              linkedComps$Comp2[l] <- colnames(compSimMat)[k]
              l <- l + 1
            }
          }
        }
        linkCol <- rep("grey40", nrow(linkedComps))
        for (i in 1:nrow(linkedComps)) {
          row1 <- match(linkedComps$Comp1[i], rownames(compoundMeanTrans))
          row2 <- match(linkedComps$Comp2[i], rownames(compoundMeanTrans))
          if (compoundMeanTrans[row1, j] == 0 | compoundMeanTrans[row2, 
                                                                  j] == 0) {
            linkCol[i] <- "transparent"
          }
        }
        compSimMatTri <- compSimMat[upper.tri(compSimMat)]
        compSimMatTriPos <- compSimMatTri[compSimMatTri > 
                                            0]
        nodes <- data.frame(name = rownames(compoundMeanTrans))
        links <- data.frame(from = linkedComps$Comp1, 
                            to = linkedComps$Comp2, weight = compSimMatTriPos, 
                            linkCol = linkCol)
        networkObjectMan <- tidygraph::tbl_graph(nodes = nodes, 
                                                 edges = links)
        p1 <- ggraph::ggraph(graph = networkObjectMan, 
                             layout = layout) + 
          ggraph::geom_edge_link(ggplot2::aes(width = .data$weight, 
                                              color = .data$linkCol, alpha = .data$linkCol)) + 
          ggraph::scale_edge_width(range = c(0.1, 2.5), 
                                   name = "Molecular similarity") + 
          ggraph::scale_edge_alpha_manual(values = c(.5,0), guide = "none") +
          ggraph::geom_node_point(ggplot2::aes(color = npcTable$pathway,fill = npcTable$pathway, size = compoundMeanTrans[, j]), shape = 21) + 
          ggraph::scale_edge_color_manual(limits = as.factor(links$linkCol),  values = links$linkCol, guide = "none") + 
          ggraph::geom_node_point(ggplot2::aes(size = compoundMeanTrans[,  j]), shape = 21, color = npcTable$col, fill = groupCol$col) + 
          ggplot2::scale_fill_manual(values = legCol) + 
          ggplot2::scale_size(range = c(4, 10)) + ggplot2::guides(color = "none") + 
          ggplot2::labs(fill = "Pathway", width = "Molecular similarity", 
                        size = "Proportion") + ggplot2::ggtitle(colnames(compoundMeanTrans)[j]) + 
          ggplot2::theme(legend.title = ggplot2::element_text(size = 10), 
                         legend.text = ggplot2::element_text(size = 8), 
                         panel.background = ggplot2::element_blank(), 
                         legend.key = ggplot2::element_blank())
        p1
      })
    }
  }
  else {
    compoundMean <- stats::aggregate(sampleData, by = list(Group = groupData), 
                                     mean)
    compoundMeanTrans <- t(compoundMean[, 2:ncol(compoundMean)])
    colnames(compoundMeanTrans) <- compoundMean$Group
    compoundMeanTrans <- as.data.frame(compoundMeanTrans)
    networkList <- list()
    for (j in 1:ncol(compoundMeanTrans)) {
      networkList[[colnames(compoundMeanTrans)[j]]] <- local({
        j <- j
        groupShape <- data.frame(shape = rep("Pos", nrow(compoundMeanTrans)))
        for (row in 1:nrow(compoundMeanTrans)) {
          if (compoundMeanTrans[row, j] == 0) {
            groupShape$shape[row] <- "Zero"
          }
        }
        compSimMat <- matrix(data = NA, nrow = nrow(compoundMeanTrans), 
                             ncol = nrow(compoundMeanTrans))
        colnames(compSimMat) <- rownames(compoundMeanTrans)
        rownames(compSimMat) <- rownames(compoundMeanTrans)
        for (i in 1:nrow(compoundMeanTrans)) {
          compSimMat[i, ] <- networkObject[i]
        }
        linkedComps <- as.data.frame(matrix(data = NA, 
                                            nrow = length(compSimMat[compSimMat > 0])/2, 
                                            ncol = 2))
        colnames(linkedComps) <- c("Comp1", "Comp2")
        l <- 1
        for (i in 1:nrow(compSimMat)) {
          for (k in i:ncol(compSimMat)) {
            if (compSimMat[i, k] > 0) {
              linkedComps$Comp1[l] <- rownames(compSimMat)[i]
              linkedComps$Comp2[l] <- colnames(compSimMat)[k]
              l <- l + 1
            }
          }
        }
        linkCol <- rep("grey40", nrow(linkedComps))
        for (i in 1:nrow(linkedComps)) {
          row1 <- match(linkedComps$Comp1[i], rownames(compoundMeanTrans))
          row2 <- match(linkedComps$Comp2[i], rownames(compoundMeanTrans))
          if (compoundMeanTrans[row1, j] == 0 | compoundMeanTrans[row2, 
                                                                  j] == 0) {
            linkCol[i] <- "grey90"
          }
        }
        compSimMatTri <- compSimMat[upper.tri(compSimMat)]
        compSimMatTriPos <- compSimMatTri[compSimMatTri > 
                                            0]
        nodes <- data.frame(name = rownames(compoundMeanTrans))
        links <- data.frame(from = linkedComps$Comp1, 
                            to = linkedComps$Comp2, weight = compSimMatTriPos, 
                            linkCol = linkCol)
        networkObjectMan <- tidygraph::tbl_graph(nodes = nodes, 
                                                 edges = links)
        p1 <- ggraph::ggraph(graph = networkObjectMan, 
                             layout = layout) + ggraph::geom_edge_link(ggplot2::aes(width = .data$weight, 
                                                                                    color = .data$linkCol)) + ggraph::scale_edge_width(range = c(0.3, 
                                                                                                                                                 2.5), name = "Molecular similarity") + ggraph::geom_node_point(ggplot2::aes(color = compoundMeanTrans[, 
                                                                                                                                                                                                                                                       j], shape = groupShape$shape), size = 10) + 
          ggraph::scale_edge_color_manual(limits = as.factor(links$linkCol), 
                                          values = links$linkCol, guide = "none") + 
          ggplot2::scale_colour_viridis_c() + ggplot2::labs(color = "Proportion", 
                                                            width = "Molecular similarity") + ggplot2::ggtitle(colnames(compoundMeanTrans)[j]) + 
          ggplot2::guides(shape = "none") + ggplot2::theme(legend.title = ggplot2::element_text(size = 10), 
                                                           legend.text = ggplot2::element_text(size = 8), 
                                                           panel.background = ggplot2::element_blank(), 
                                                           legend.key = ggplot2::element_blank())
        print(p1)
      })
    }
  }
  return(networkList)
}

