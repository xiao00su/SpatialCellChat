#' Compute the network centrality scores allowing identification of 
#' dominant senders, receivers, mediators and influencers in all inferred communication networks
#' NB: This function was previously named as `netAnalysis_signalingRole`.  
#' The previous function `netVisual_signalingRole` is now named as `netAnalysis_signalingRole_network`.
#' 
#' @param object CellChat object; If object = NULL, USER must provide `net`
#' @param net compute the centrality measures on a specific signaling network given by a 2 or 3 dimemsional array net
#' @param slot.name the slot name of object that is used to compute centrality measures of signaling networks
#' @param signaling.name provide signalings name (LR pairs or pathways) to subset the communication net
#' @param do.group set `do.group = TRUE` when computing centrality of each cell group; 
#' set `do.group = FALSE` when computing centrality of each individual cell
#' @param thresh threshold of the p-value for determining significant interaction
#' @param degree.only To speed up computation, run `netAnalysis_computeCentrality` with `degree.only = T` 
#' when user only needs "outdeg_unweighted", "indeg_unweighted","outdeg","indeg","page_rank" to do analysis
#' @importFrom methods slot
#' @importFrom future.apply future_lapply
#'
#' @return slot.name下, do.group=TRUE返回centr, 
#' do.froup = FALSE, 返回centr.cell
#' 
#' @export
#' 
netAnalysis_computeCentralityX <- function(object = NULL, net=NULL, slot.name = "net", 
                                           signaling.name = NULL, do.group = F, 
                                           thresh = 0.05, degree.only=T){
  # 真正用于下游计算的是: 
  # 1. do.group = TRUE, object@net$prob的list形式
  # 2. do.group = FALSE, object@net$tmp$prob.cell
  # 计算centr的核心函数computeCentralityLocal, 其输入数据是list中的一个元素
  
  if (is.null(net)) { #### 输入数据是 object
    if (do.group) {   #### cell cluster
      prob <- methods::slot(object, slot.name)$prob  # object@net$prob
      pval <- methods::slot(object, slot.name)$pval
      pval[prob == 0] <- 1
      prob[pval >= thresh] <- 0 # 过滤
      # 将3Darray转成list
      net <- BiocGenerics::lapply(X = seq_len(dim(prob)[[3]]), 
                                  FUN = function(i){prob[ , ,i,drop=T]} )
      names(net) <- dimnames(prob)[[3]]
      node.names <- dimnames(prob)[[1]]
    } else {  ##### cell
      if (is.null(methods::slot(object, slot.name)$tmp$prob.cell)) {
        if (slot.name == "net") {
          stop( cli.symbol(2), "Please run `computeCommunProb` to compute 
          the communication probability/strength between any interacting individual cells!" )
        } else if (slot.name == "netP") {
          stop(cli.symbol(2), "Please run `computeCommunProbPathway` to compute 
          the communication probability/strength  between any interacting individual cells!" )
        }
      } else {
        net <- methods::slot(object, slot.name)$tmp$prob.cell # a list
        #### node.names <- spatstat.sparse::dimnames.sparse3Darray(methods::slot(object, slot.name)$prob.cell)[[1]]
        node.names <- colnames(object@data.signaling) # cell_id
      }
    } 
  } else {  ##### 输入数据是net
    if (length(dim(net)) == 2) {  #### net是2Darray, ???
      node.names=seq_len(NROW(net))
      net=list("net_temp"=net)
    } else if(length(dim(net)) == 3) { #### net是3Darray, 此步骤将3Darray转为list
      net.names <- dimnames(net)[[3]]
      node.names <- dimnames(net)[[1]]
      net = BiocGenerics::lapply(X = seq_len(dim(net)[[3]]), FUN = function(i){net[ , ,i,drop=T]} )
      names(net) <- net.names
    }
  }# is.null(net)?
  
  # 指定signaling.name
  if (!is.null(signaling.name)) {
    if (!all(signaling.name %in% names(net) )) {
      stop("Please check the input `signaling.name` because some are not the significant signaling!") 
    }
    net <- net[signaling.name] # [] will return a list
  }
  
  signaling.name <- names(net)
  
  N <- dim( net[[1]] )[1]
  nrun <- length(signaling.name)
  
  if(degree.only == T){
    centr.name <- c("outdeg_unweighted", "indeg_unweighted","outdeg","indeg","page_rank")
    centr.all <-  my_future_sapply(
      X = 1:nrun,
      FUN = function(x) {
        net0 <- net[[x]]
        centr.x <- computeCentralityLocal(net0, degree.only = degree.only)
        gc() 
        return(centr.x)
      },
      simplify = TRUE
    )
  } else {
    centr.name <- c("outdeg_unweighted", "indeg_unweighted",
                    "outdeg", "indeg", "hub", "authority", "eigen",
                    "page_rank", "betweenness", "flowbet", "info")
    centr.all <-  my_future_sapply(X = 1:nrun, FUN = function(x) { 
      net0 <- net[[x]]
      centr.x <- computeCentralityLocal(net0,degree.only = degree.only)
      gc()
      return(centr.x)
      },
      simplify = TRUE
    )
  }

  centr.all <- reticulate::array_reshape(centr.all, c(nrow(centr.all)/N,N, nrun), order = "F")

  dimnames(centr.all) <- list(centr.name, node.names,signaling.name)
  cat(cli.symbol(1), "Computing Net Centrality is done.\n")
  
  if (is.null(object)) { 
    return(centr.all)
  } else { 
    if (do.group) {
      methods::slot(object, slot.name)[["centr"]] <- centr.all
    } else { methods::slot(object, slot.name)[["centr.cell"]] <- centr.all }
    
    return(object)
  }
}

#' Compute and visualize the contribution of each ligand-receptor pair in the overall signaling pathways
#'
#' @param object CellChat object
#' @param signaling a signaling pathway name
#' @param signaling.name alternative signaling pathway name to show on the plot
#' @param do.group set `do.group = TRUE` when only showing enriched signaling based on cell group-level communication; 
#' set `do.group = FALSE` when only showing enriched signaling based on individual cell-level communication
#' @param width the width of individual bar
#' @param vertex.receiver a numeric vector giving the index of the cell groups as targets in the first hierarchy plot
#' @param thresh threshold of the p-value for determining significant interaction
#' @param return.data whether return the data.frame consisting of the predicted L-R pairs and their contribution
#' @param x.rotation rotation of x-label
#' @param title the title of the plot
#' @param font.size font size of the text
#' @param font.size.title font size of the title
#' @importFrom dplyr select
#' @importFrom ggplot2 ggplot geom_bar aes coord_flip scale_x_discrete element_text theme ggtitle
#' @importFrom cowplot ggdraw draw_label plot_grid
#'
#' @return
#' @export
#'
#' @examples
netAnalysis_contributionX <- function(object, signaling, signaling.name = NULL, 
                                      do.group = TRUE, width = 0.1, vertex.receiver = NULL, 
                                      thresh = 0.05, return.data = FALSE, 
                                      x.rotation = 0, title = "Contribution of each L-R pair", 
                                      font.size = 10, font.size.title = 10) {
  pairLR <- searchPair(signaling = signaling, pairLR.use = object@LR$LRsig, 
                       key = "pathway_name", matching.exact = T, pair.only = T)
  pair.name.use = select(object@DB$interaction[rownames(pairLR),],"interaction_name_2")
  
  if (is.null(signaling.name)) {signaling.name <- signaling }
  if (do.group) {
    net <- object@net
    pairLR.use.name <- dimnames(net$prob)[[3]]
    pairLR.name <- intersect(rownames(pairLR), pairLR.use.name)
    pairLR <- pairLR[pairLR.name, ]
    prob <- net$prob
    pval <- net$pval
    prob[pval > thresh] <- 0
    cell.level <- FALSE
  } else {
    net <- object@net
    tl <- net$tmp$prob.cell
    
    if (!is.list(tl) || length(tl) == 0L) {
      stop("net$tmp$prob.cell 不存在或为空: 请先运行 computeCommunProb, 或改用 do.group = TRUE")
    }
    
    pairLR.use.name <- names(tl)
    pairLR.name <- intersect(rownames(pairLR), pairLR.use.name)
    pairLR <- pairLR[pairLR.name, ]
    layer.sum <- vapply(tl[pairLR.name], function(m) sum(m@x), numeric(1))
    names(layer.sum) <- pairLR.name
    prob.layers <- NULL
    cell.level <- TRUE
  }

  if (cell.level) { pairLR.name.use <- pairLR.name[layer.sum != 0]
  } else if (length(pairLR.name) > 1) {
    pairLR.name.use <- pairLR.name[apply(prob[,,pairLR.name], 3, sum) != 0]
  } else { pairLR.name.use <- pairLR.name[sum(prob[,,pairLR.name]) != 0] }
  
  
  if (length(pairLR.name.use) == 0) {
    stop(paste0('There is no significant communication of ', signaling.name))
  } else { pairLR <- pairLR[pairLR.name.use,] }
  
  if (cell.level) {
    prob.layers <- tl[pairLR.name.use]
    pSum.all  <- vapply(prob.layers, function(m) sum(m@x), numeric(1))
    pSum.max  <- sum(pSum.all)
    N.cells   <- nrow(tl[[1]])
    n.layer   <- length(prob.layers)
  } else {
    prob <- prob[,,pairLR.name.use]

    if (length(dim(prob)) == 2) {
      prob <- replicate(1, prob, simplify="array")
      dimnames(prob)[3] <- pairLR.name.use
    }
    prob <-(prob-min(prob))/(max(prob)-min(prob))
    n.layer <- dim(prob)[3]
  }
  
  if (is.null(vertex.receiver)) {
    if (cell.level) { pSum <- pSum.all
    } else {
      pSum <- apply(prob, 3, sum)
      pSum.max <- sum(prob)
    }
    
    pSum <- pSum/pSum.max
    pSum[is.na(pSum)] <- 0
    y.lim <- max(pSum)

    pair.name <- if (cell.level) pairLR.name.use else unlist(dimnames(prob)[3])
    pair.name <- factor(pair.name, levels = unique(pair.name))
    
    if (!is.null(pairLR.name.use)) {
      pair.name <- pair.name.use[as.character(pair.name),1]
      pair.name <- factor(pair.name, levels = unique(pair.name))
    }
    
    mat <- pSum
    df1 <- data.frame(name = pair.name, contribution = mat)
    
    if(nrow(df1) < 10) {
      df2 <- data.frame(name = as.character(1:(10-nrow(df1))), contribution = rep(0, 10-nrow(df1)))
      df <- rbind(df1, df2)
    } else { df <- df1 }
    
    df <- df[order(df$contribution, decreasing = TRUE), ]
    # df$name <- factor(df$name, levels = unique(df$name))
    df$name <- factor(df$name,levels=df$name[order(df$contribution, decreasing = TRUE)])
    df1$name <- factor(df1$name,levels=df1$name[order(df1$contribution, decreasing = TRUE)])
    gg <- ggplot(df, aes(x=name, y=contribution)) + geom_bar(stat="identity", width = 0.7) +
      theme_classic() + 
      theme(axis.text.y = element_text(angle = x.rotation, hjust = 1, size=font.size, colour = 'black'), 
            axis.text=element_text(size=font.size),
            axis.title.y = element_text(size= font.size), 
            axis.text.x = element_blank(), 
            axis.ticks = element_blank()) +
      xlab("") + ylab("Relative contribution") + ylim(0,y.lim) + coord_flip() + 
      theme(legend.position="none") +
      scale_x_discrete(limits = rev(levels(df$name)), 
                       labels = c(rep("", max(0, 10-nlevels(df1$name))),rev(levels(df1$name))))
    
    if (!is.null(title)) {
      gg <- gg + ggtitle(title) + 
        theme(plot.title = element_text(hjust = 0.5, size = font.size.title))
    }
    gg
    
  } else {
    pn3 <- if (cell.level) pairLR.name.use else unlist(dimnames(prob)[3])
    pair.name <- factor(pn3, levels = unique(pn3))
    # show all the communications
    if (cell.level) {
      pSum <- pSum.all
    } else {
      pSum <- apply(prob, 3, sum)
      pSum.max <- sum(prob)
    }
    
    pSum <- pSum/pSum.max
    pSum[is.na(pSum)] <- 0
    y.lim <- max(pSum)
    
    df<- data.frame(name = pair.name, contribution = pSum)
    gg <- ggplot(df, aes(x=name, y=contribution)) + geom_bar(stat="identity",width = 0.2) +
      theme_classic() + theme(axis.text=element_text(size=10),axis.text.x = element_text(angle = x.rotation, hjust = 1,size=8),
                              axis.title.y = element_text(size=10)) +
      xlab("") + ylab("Relative contribution") + ylim(0,y.lim)+ ggtitle("All")+ theme(plot.title = element_text(hjust = 0.5))#+
    
    # show the communications in Hierarchy1
    if (cell.level) {
      # 每层在指定细胞(列)子集上的总和, 稀疏直读 O(nnz), 不稠密化
      pSum <- if (n.layer > 1) {
        vapply(prob.layers, function(m) sum(Matrix::colSums(m)[vertex.receiver]), numeric(1))
      } else {
        sum(Matrix::colSums(prob.layers[[1]])[vertex.receiver])
      }
    } else if (dim(prob)[3] > 1) {
      pSum <- apply(prob[,vertex.receiver,], 3, sum)
    } else {
      pSum <- sum(prob[,vertex.receiver,])
    }
    
    pSum <- pSum/pSum.max
    pSum[is.na(pSum)] <- 0
    
    df<- data.frame(name = pair.name, contribution = pSum)
    gg1 <- ggplot(df, aes(x=name, y=contribution)) + geom_bar(stat="identity",width = 0.2) +
      theme_classic() + 
      theme(axis.text = element_text(size=10),
            axis.text.x = element_text(angle = x.rotation, hjust = 1,size=8), 
            axis.title.y = element_text(size=10)) +
      xlab("") + ylab("Relative contribution") + ylim(0,y.lim)+ ggtitle("Hierarchy1") + 
      theme(plot.title = element_text(hjust = 0.5))#+

    
    # show the communications in Hierarchy2
    
    if (cell.level) {
      cols2 <- setdiff(seq_len(N.cells), vertex.receiver)
      pSum <- if (n.layer > 1) {
        vapply(prob.layers, function(m) sum(Matrix::colSums(m)[cols2]), numeric(1))
      } else {
        sum(Matrix::colSums(prob.layers[[1]])[cols2])
      }
    } else if (dim(prob)[3] > 1) {
      pSum <- apply(prob[,setdiff(1:dim(prob)[1],vertex.receiver),], 3, sum)
    } else {
      pSum <- sum(prob[,setdiff(1:dim(prob)[1],vertex.receiver),])
    }
    pSum <- pSum/pSum.max
    pSum[is.na(pSum)] <- 0
    
    df<- data.frame(name = pair.name, contribution = pSum)
    gg2 <- ggplot(df, aes(x=name, y=contribution)) + geom_bar(stat="identity", width=0.9) +
      theme_classic() + 
      theme(axis.text = element_text(size=10), axis.title.y = element_text(size=10),
            axis.text.x = element_text(angle = x.rotation, hjust = 1,size=8)) +
      xlab("") + ylab("Relative contribution") + ylim(0,y.lim) + 
      ggtitle("Hierarchy2")+ theme(plot.title = element_text(hjust = 0.5))#+
    
    #scale_x_discrete(limits = c(0,1))
    title <- cowplot::ggdraw() + 
      cowplot::draw_label(paste0("Contribution of each signaling in ", signaling.name, " pathway"), fontface='bold', size = 10)
    gg.combined <- cowplot::plot_grid(gg, gg1, gg2, nrow = 1)
    gg.combined <- cowplot::plot_grid(title, gg.combined, ncol = 1, rel_heights=c(0.1, 1))
    gg <- gg.combined
    gg
  }
  if (return.data) {
    df <- subset(df, contribution > 0)
    return(list(LR.contribution = df, gg.obj = gg))
  } else {
    return(gg)
  }
}




#' Rank signaling networks based on the information flow or the number of interactions
#'
#' This function can also be used to rank signaling from certain cell groups to other cell groups
#'
#' @param object CellChat object
#' @param slot.name the slot name of object that is used to compute centrality measures of signaling networks
#' @param measure "weight" or "count". "weight": comparing the total interaction weights (strength); "count": comparing the number of interactions;
#' @param mode "single","comparison"
#' @param comparison a numerical vector giving the datasets for comparison; a single value means ranking for only one dataset and two values means ranking comparison for two datasets
#' @param do.group set `do.group = TRUE` when only showing enriched signaling based on cell group-level communication; set `do.group = FALSE` when only showing enriched signaling based on individual cell-level communication
#' @param color.use defining the color for each cell group
#' @param stacked whether plot the stacked bar plot
#' @param sources.use a vector giving the index or the name of source cell groups
#' @param targets.use a vector giving the index or the name of target cell groups.
#' @param signaling a vector giving the signaling pathway to show
#' @param pairLR a vector giving the names of L-R pairs to show (e.g, pairLR = c("IL1A_IL1R1_IL1RAP","IL1B_IL1R1_IL1RAP"))
#' @param signaling.type a char giving the types of signaling from the three categories c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact")
#' @param do.stat whether do a paired Wilcoxon test to determine whether there is significant difference between two datasets. Default = FALSE
#' @param cutoff.pvalue the cutoff of pvalue when doing Wilcoxon test; Default = 0.05
#' @param tol a tolerance when considering the relative contribution being equal between two datasets. contribution.relative between 1-tol and 1+tol will be considered as equal contribution
#' @param thresh threshold of the p-value for determining significant interaction
#'
#' @param do.flip whether flip the x-y axis
#' @param x.angle,y.angle,x.hjust,y.hjust parameters for rotating and spacing axis labels
#' @param axis.gap whetehr making gaps in y-axes
#' @param ylim,segments,tick_width,rel_heights parameters in the function gg.gap when making gaps in y-axes
#' e.g., ylim = c(0, 35), segments = list(c(11, 14),c(16, 28)), tick_width = c(5,2,5), rel_heights = c(0.8,0,0.1,0,0.1)
#' https://tobiasbusch.xyz/an-r-package-for-everything-ep2-gaps
#' @param show.raw whether show the raw information flow. Default = FALSE, showing the scaled information flow to provide compariable data scale; When stacked = TRUE, use raw information flow by default.
#' @param return.data whether return the data.frame consisting of the calculated information flow of each signaling pathway or L-R pair
#' @param x.rotation rotation of x-labels
#' @param title main title of the plot
#' @param bar.w the width of bar plot
#' @param font.size font size
#' @param legend.size the size of legend
#' @param legend.text.size the text size on the legend
#' @param legend.position parameters for configurating the plot
#' @param legend.spacing a two-elements vector respectively specifying legend.key.spacing.x and legend.key.spacing.y for spacing apart legend key-label pairs
#' @import ggplot2
#' @importFrom methods slot
#' @return
#' @export
#'
#' @examples
rankNetX <- function(object, slot.name = "netP", measure = c("weight","count"), 
                     mode = c("comparison", "single"), comparison = c(1,2), 
                     do.group = TRUE, color.use = NULL, stacked = FALSE,  
                     sources.use = NULL, targets.use = NULL,  signaling = NULL, 
                     pairLR = NULL, signaling.type = NULL, do.stat = FALSE, 
                     cutoff.pvalue = 0.05, tol = 0.05, thresh = 0.05, 
                     show.raw = FALSE, return.data = FALSE, 
                     x.rotation = 90, title = NULL, bar.w = 0.75, font.size = 8, 
                     do.flip = TRUE, x.angle = NULL, y.angle = 0, x.hjust = 1,y.hjust = 1,
                     axis.gap = FALSE, ylim = NULL, segments = NULL, 
                     tick_width = NULL, rel_heights = c(0.9,0,0.1), 
                     legend.size = 0.1, legend.text.size = 8, 
                     legend.position = "top", legend.spacing = c(2, -8)) {
  measure <- match.arg(measure)
  mode <- match.arg(mode)
  options(warn = -1)
  object.names <- names(methods::slot(object, slot.name))
  
  if (measure == "weight") { ylabel = "Information flow"
  } else if (measure == "count") { ylabel = "Number of interactions" }
  
  if (mode == "single") {
    object1 <- methods::slot(object, slot.name)  # object@netP
    
    if (do.group) { # cell cluster
      prob = object1$prob
      prob[object1$pval > thresh] <- 0
      if (measure == "count") { prob <- 1*(prob > 0) }
      
      if (!is.null(sources.use)) {
        if (is.character(sources.use)) {
          if (all(sources.use %in% dimnames(prob)[[1]])) {
            sources.use <- match(sources.use, dimnames(prob)[[1]])
          } else {
            stop("The input `sources.use` should be cell group names or a numerical vector!")
          }
        }
        idx.t <- setdiff(1:nrow(prob), sources.use)
        prob[idx.t, , ] <- 0
      }
      
      if (!is.null(targets.use)) {
        if (is.character(targets.use)) {
          if (all(targets.use %in% dimnames(prob)[[1]])) {
            targets.use <- match(targets.use, dimnames(prob)[[2]])
          } else {
            stop("The input `targets.use` should be cell group names or a numerical vector!")
          }
        }
        idx.t <- setdiff(1:nrow(prob), targets.use)
        prob[ ,idx.t, ] <- 0
      }
      
      if (sum(prob) == 0) { stop("No inferred communications for the input!") }

      pSum <- apply(prob, 3, sum)
      pSum.original <- pSum
      
      if (measure == "weight") {
        pSum <- -1/log(pSum)
        pSum[is.na(pSum)] <- 0
        idx1 <- which(is.infinite(pSum) | pSum < 0)
        values.assign <- seq(max(pSum)*1.1, max(pSum)*1.5, length.out = length(idx1))
        position <- sort(pSum.original[idx1], index.return = TRUE)$ix
        pSum[idx1] <- values.assign[match(1:length(idx1), position)]
      } else if (measure == "count") { pSum <- pSum.original }

      pair.name <- names(pSum)
    } else {
      tl <- object1$tmp$prob.cell
      
      if (!is.list(tl) || length(tl) == 0L) { 
        stop("Individual cell-level `net$tmp$prob.cell` is empty; nothing to rank!") }
      
      cell.names <- colnames(object@data.signaling)
      n.cell <- nrow(tl[[1]])
      if (!is.null(cell.names) && length(cell.names) != n.cell) {
        stop("ncol(data.signaling) = ", length(cell.names), " != nrow(prob.cell) = ", n.cell,
             "; the cell order/mapping is inconsistent, refusing to match by name!")
      }

      if (!is.null(sources.use)) {
        if (is.character(sources.use)) {
          if (!is.null(cell.names) && all(sources.use %in% cell.names)) {
            sources.use <- match(sources.use, cell.names)
          } else {
            stop("The input `sources.use` should be cell names or a numerical vector!")
          }
        }
        if (length(sources.use) == n.cell) { sources.use <- NULL }
      }

      if (!is.null(targets.use)) {
        if (is.character(targets.use)) {
          if (!is.null(cell.names) && all(targets.use %in% cell.names)) {
            targets.use <- match(targets.use, cell.names)
          } else {
            stop("The input `targets.use` should be cell names or a numerical vector!")
          }
        }
        if (length(targets.use) == n.cell) { targets.use <- NULL }
      }
      
      layerSum <- function(m) {
        if (is.null(sources.use) && is.null(targets.use)) {
          return(if (measure == "count") sum(m@x > 0) else sum(m@x))
        }
        sub <- if (is.null(sources.use)) {
          m[, targets.use, drop = FALSE]
        } else if (is.null(targets.use)) {
          m[sources.use, , drop = FALSE]
        } else {
          m[sources.use, targets.use, drop = FALSE]
        }
        if (measure == "count") sum(sub@x > 0) else sum(sub@x)
      }
      
      pSum.original <- vapply(tl, layerSum, numeric(1))
      names(pSum.original) <- names(tl)
      
      if (sum(pSum.original) == 0) { stop("No inferred communications for the input!") }
      ## 与原函数一致：cell 级 weight 不做 -1/log 变换，直接用 raw 信息流
      pSum <- pSum.original
      pair.name <- names(pSum)
    }
    
    df<- data.frame(name = pair.name, contribution = pSum.original, contribution.scaled = pSum, group = object.names[comparison[1]])
    idx <- with(df, order(df$contribution))
    df <- df[idx, ]
    df$name <- factor(df$name, levels = as.character(df$name))
    for (i in 1:length(pair.name)) {
      df.t <- df[df$name == pair.name[i], "contribution"]
      if (sum(df.t) == 0) {
        df <- df[-which(df$name == pair.name[i]), ]
      }
    }
    
    if (!is.null(signaling.type)) {
      LR <- subset(object@DB$interaction, annotation %in% signaling.type)
      if (slot.name == "netP") { signaling <- unique(LR$pathway_name)
      } else if (slot.name == "net") { pairLR <- LR$interaction_name }
    }
    
    if ((slot.name == "netP") && (!is.null(signaling))) {
      df <- subset(df, name %in% signaling)
    } else if ((slot.name == "netP") &&(!is.null(pairLR))) {
      stop("You need to set `slot.name == 'net'` if showing specific L-R pairs ")
    }
    
    if ((slot.name == "net") && (!is.null(pairLR))) {
      df <- subset(df, name %in% pairLR)
    } else if ((slot.name == "net") && (!is.null(signaling))) {
      stop("You need to set `slot.name == 'netP'` if showing specific signaling pathways ")
    }
    
    gg <- ggplot(df, aes(x=name, y=contribution.scaled)) + 
      geom_bar(stat="identity",width = bar.w) +
      theme_classic() + 
      theme(axis.text=element_text(size=10), axis.text.x = element_blank(), 
            axis.ticks.x = element_blank(), axis.title.y = element_text(size=10)) +
      xlab("") + ylab(ylabel) + coord_flip()
    
    if (!is.null(title)) { gg <- gg + ggtitle(title)+ theme(plot.title = element_text(hjust = 0.5)) }
    
  } else if (mode == "comparison") {
    prob.list <- list()
    pSum <- list()
    pSum.original <- list()
    pair.name <- list()
    idx <- list()
    pSum.original.all <- c()
    object.names.comparison <- c()
    for (i in 1:length(comparison)) {
      object.list <- methods::slot(object, slot.name)[[comparison[i]]]
      if (do.group) {
        prob <- object.list$prob
        prob[object.list$pval > thresh] <- 0
        if (measure == "count") { prob <- 1*(prob > 0) }
        prob.list[[i]] <- prob
        
        if (!is.null(sources.use)) {
          if (is.character(sources.use)) {
            if (all(sources.use %in% dimnames(prob)[[1]])) {
              sources.use <- match(sources.use, dimnames(prob)[[1]])
            } else {
              stop("The input `sources.use` should be cell group names or a numerical vector!")
            }
          }
          idx.t <- setdiff(1:nrow(prob), sources.use)
          prob[idx.t, , ] <- 0
        }
        
        if (!is.null(targets.use)) {
          if (is.character(targets.use)) {
            if (all(targets.use %in% dimnames(prob)[[1]])) {
              targets.use <- match(targets.use, dimnames(prob)[[2]])
            } else {
              stop("The input `targets.use` should be cell group names or a numerical vector!")
            }
          }
          idx.t <- setdiff(1:nrow(prob), targets.use)
          prob[ ,idx.t, ] <- 0
        }
        
        if (sum(prob) == 0) { stop("No inferred communications for the input!") }
        pSum.original[[i]] <- apply(prob, 3, sum)
      } else {
        ## cell level: list(dgCMatrix) 逐层求和，不物化 3D array
        tl <- object.list$tmp$prob.cell
        if (!is.list(tl) || length(tl) == 0L) {
          stop("Individual cell-level `tmp$prob.cell` is empty for dataset ", comparison[i], "; nothing to rank!")
        }
        
        prob.list[[i]] <- tl
        cell.names <- colnames(object@data.signaling)
        n.cell <- nrow(tl[[1]])
        
        if (!is.null(cell.names) && length(cell.names) != n.cell) {
          stop("dataset ", comparison[i], ": ncol(data.signaling) = ", length(cell.names),
               " != nrow(prob.cell) = ", n.cell,
               "; the cell order/mapping is inconsistent, refusing to match by name!")
        }
        
        if (!is.null(sources.use)) {
          if (is.character(sources.use)) {
            if (!is.null(cell.names) && all(sources.use %in% cell.names)) {
              sources.use <- match(sources.use, cell.names)
            } else {stop("The input `sources.use` should be cell names or a numerical vector!")}
          }
          if (length(sources.use) == n.cell) { sources.use <- NULL }
        }
        
        if (!is.null(targets.use)) {
          if (is.character(targets.use)) {
            if (!is.null(cell.names) && all(targets.use %in% cell.names)) {
              targets.use <- match(targets.use, cell.names)
            } else { stop("The input `targets.use` should be cell names or a numerical vector!")}
          }
          if (length(targets.use) == n.cell) { targets.use <- NULL }
        }
        
        layerSum <- function(m) {
          if (is.null(sources.use) && is.null(targets.use)) {
            return(if (measure == "count") sum(m@x > 0) else sum(m@x))
          }
          sub <- if (is.null(sources.use)) {
            m[, targets.use, drop = FALSE]
          } else if (is.null(targets.use)) {
            m[sources.use, , drop = FALSE]
          } else {
            m[sources.use, targets.use, drop = FALSE]
          }
          if (measure == "count") sum(sub@x > 0) else sum(sub@x)
        }
        
        pSum.original[[i]] <- vapply(tl, layerSum, numeric(1))
        names(pSum.original[[i]]) <- names(tl)
        if (sum(pSum.original[[i]]) == 0) { stop("No inferred communications for the input!") }
      }
      
      if (measure == "weight") {
        if (do.group) {
          pSum[[i]] <- -1/log(pSum.original[[i]])
          pSum[[i]][is.na(pSum[[i]])] <- 0
          idx[[i]] <- which(is.infinite(pSum[[i]]) | pSum[[i]] < 0)
          pSum.original.all <- c(pSum.original.all, pSum.original[[i]][idx[[i]]])
        } else {
          pSum[[i]] <- pSum.original[[i]]
          pSum.original.all <- c(pSum.original.all, pSum.original[[i]])
        }
      } else if (measure == "count") { pSum[[i]] <- pSum.original[[i]] }
      pair.name[[i]] <- names(pSum.original[[i]])
      object.names.comparison <- c(object.names.comparison, object.names[comparison[i]])
    }
    
    if (measure == "weight" & do.group == TRUE) {
      values.assign <- seq(max(unlist(pSum))*1.1, max(unlist(pSum))*1.5, length.out = length(unlist(idx)))
      position <- sort(pSum.original.all, index.return = TRUE)$ix
      for (i in 1:length(comparison)) {
        if (i == 1) {
          pSum[[i]][idx[[i]]] <- values.assign[match(1:length(idx[[i]]), position)]
        } else {
          pSum[[i]][idx[[i]]] <- values.assign[match(length(unlist(idx[1:i-1]))+1:length(unlist(idx[1:i])), position)]
        }
      }
    }
    
    pair.name.all <- as.character(unique(unlist(pair.name)))
    df <- list()
    for (i in 1:length(comparison)) {
      df[[i]] <- data.frame(name = pair.name.all, contribution = 0, contribution.scaled = 0, 
                            group = object.names[comparison[i]], row.names = pair.name.all)
      df[[i]][pair.name[[i]],3] <- pSum[[i]]
      df[[i]][pair.name[[i]],2] <- pSum.original[[i]]
    }
    
    
    contribution.relative <- list()
    for (i in 1:(length(comparison)-1)) {
      contribution.relative[[i]] <- as.numeric(format(df[[length(comparison)-i+1]]$contribution/df[[1]]$contribution, digits=1))
      contribution.relative[[i]][is.na(contribution.relative[[i]])] <- 0
    }
    
    names(contribution.relative) <- paste0("contribution.relative.", 1:length(contribution.relative))
    
    for (i in 1:length(comparison)) {
      for (j in 1:length(contribution.relative)) {
        df[[i]][[names(contribution.relative)[j]]] <- contribution.relative[[j]]
      }
    }
    
    df[[1]]$contribution.data2 <- df[[length(comparison)]]$contribution
    
    if (length(comparison) == 2) {
      idx <- with(df[[1]], order(-contribution.relative.1, contribution, -contribution.data2))
    } else if (length(comparison) == 3) {
      idx <- with(df[[1]], order(-contribution.relative.1, -contribution.relative.2,contribution, -contribution.data2))
    } else if (length(comparison) == 4) {
      idx <- with(df[[1]], order(-contribution.relative.1, -contribution.relative.2, -contribution.relative.3, contribution, -contribution.data2))
    } else {
      idx <- with(df[[1]], order(-contribution.relative.1, -contribution.relative.2, -contribution.relative.3, -contribution.relative.4, contribution, -contribution.data2))
    }
    
    for (i in 1:length(comparison)) {
      df[[i]] <- df[[i]][idx, ]
      df[[i]]$name <- factor(df[[i]]$name, levels = as.character(df[[i]]$name))
    }
    df[[1]]$contribution.data2 <- NULL
    
    df <- do.call(rbind, df)
    df$group <- factor(df$group, levels = object.names.comparison)
    
    if (is.null(color.use)) {
      color.use =  ggPalette(length(comparison))
    }
    
    df$group <- factor(df$group, levels = rev(levels(df$group)))
    color.use <- rev(color.use)
    

    if (do.stat & !do.group) {
      stop("`do.stat` is not applicable to individual cell-level networks: 
           it would densify an N x N matrix per L-R pair (N = 165k -> ~216 GB each)! 
           Please set `do.stat = FALSE`")
    }
    
    if (do.stat & length(comparison) == 2) {
      for (i in 1:length(pair.name.all)) {
        if (nrow(prob.list[[j]]) != nrow(prob.list[[1]])) {
          stop("Statistical test is not applicable to datasets with different cellular compositions! 
               Please set `do.stat = FALSE`")
        }
        prob.values <- matrix(0, nrow = nrow(prob.list[[1]]) * nrow(prob.list[[1]]), ncol = length(comparison))
        for (j in 1:length(comparison)) {
          if (pair.name.all[i] %in% pair.name[[j]]) {
            prob.values[, j] <- as.vector(prob.list[[j]][ , , pair.name.all[i]])
          } else {
            prob.values[, j] <- NA
          }
        }
        prob.values <- prob.values[rowSums(prob.values, na.rm = TRUE) != 0, , drop = FALSE]
        if (nrow(prob.values) >3 & sum(is.na(prob.values)) == 0) {
          pvalues <- wilcox.test(prob.values[ ,1], prob.values[ ,2], paired = TRUE)$p.value
        } else {
          pvalues <- 0
        }
        pvalues[is.na(pvalues)] <- 0
        df$pvalues[df$name == pair.name.all[i]] <- pvalues
      }
    }
    
    
    if (length(comparison) == 2) {
      if (do.stat) {
        colors.text <- ifelse((df$contribution.relative < 1-tol) & (df$pvalues < cutoff.pvalue), 
                              color.use[2], 
                              ifelse((df$contribution.relative > 1+tol) & df$pvalues < cutoff.pvalue, color.use[1], "black"))
      } else {
        colors.text <- ifelse(df$contribution.relative < 1-tol, color.use[2], 
                              ifelse(df$contribution.relative > 1+tol, color.use[1], "black"))
      }
    } else {
      message("The text on the y-axis will not be colored for the number of compared datasets larger than 3!")
      colors.text = NULL
    }
    
    for (i in 1:length(pair.name.all)) {
      df.t <- df[df$name == pair.name.all[i], "contribution"]
      if (sum(df.t) == 0) {
        df <- df[-which(df$name == pair.name.all[i]), ]
      }
    }
    
    if ((slot.name == "netP") && (!is.null(signaling))) {
      df <- subset(df, name %in% signaling)
    } else if ((slot.name == "netP") &&(!is.null(pairLR))) {
      stop("You need to set `slot.name == 'net'` if showing specific L-R pairs ")
    }
    if ((slot.name == "net") && (!is.null(pairLR))) {
      df <- subset(df, name %in% pairLR)
    } else if ((slot.name == "net") && (!is.null(signaling))) {
      stop("You need to set `slot.name == 'netP'` if showing specific signaling pathways ")
    }
    
    if (stacked) {
      gg <- ggplot(df, aes(x=name, y=contribution, fill = group)) + 
        geom_bar(stat="identity",width = bar.w, position ="fill") # +
      # xlab("") + ylab("Relative information flow") #+ theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
      #  scale_y_discrete(breaks=c("0","0.5","1")) +
      if (measure == "weight") {
        gg <- gg + xlab("") + ylab("Relative information flow")
      } else if (measure == "count") {
        gg <- gg + xlab("") + ylab("Relative number of interactions")
      }
      
      gg <- gg + geom_hline(yintercept = 0.5, linetype="dashed", color = "grey50", size=0.5)
    } else {
      if (show.raw) {
        gg <- ggplot(df, aes(x=name, y=contribution, fill = group)) + 
          geom_bar(stat="identity",width = bar.w, position = position_dodge(0.8)) +
          xlab("") + ylab(ylabel) #+ coord_flip()#+ theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
      } else {
        gg <- ggplot(df, aes(x=name, y=contribution.scaled, fill = group)) + 
          geom_bar(stat="identity",width = bar.w, position = position_dodge(0.8)) +
          xlab("") + ylab(ylabel) #+ coord_flip()#+ theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
      }
      
      if (axis.gap) {
        gg <- gg + theme_bw() + theme(panel.grid = element_blank())
        gg.gap::gg.gap(gg,
                       ylim = ylim,
                       segments = segments,
                       tick_width = tick_width,
                       rel_heights = rel_heights)
      }
    }
    gg <- gg +  CellChat_theme_opts() + theme_classic()
    if (do.flip) {
      gg <- gg + coord_flip() + theme(axis.text.y = element_text(colour = colors.text))
      if (is.null(x.angle)) {
        x.angle = 0
      }
      
    } else {
      if (is.null(x.angle)) {
        x.angle = 45
      }
      gg <- gg + scale_x_discrete(limits = rev) + 
        theme(axis.text.x = element_text(colour = rev(colors.text)))
      
    }
    
    gg <- gg + theme(axis.text=element_text(size=font.size), 
                     axis.title = element_text(size=font.size))
    gg <- gg + scale_fill_manual(name = "", values = color.use)
    gg <- gg + guides(fill = guide_legend(reverse = TRUE))
    gg <- gg + theme(legend.position = legend.position,
                     legend.key.spacing.y = unit(legend.spacing[2], 'pt'), 
                     legend.key.spacing.x = unit(legend.spacing[1], 'pt')) +
      theme(legend.title = element_blank(), 
            legend.key.size = unit(legend.size, "inches"), 
            legend.text = element_text(size = legend.text.size, margin = margin(l = 0)))# ,
    gg <- gg + theme(axis.text.x = element_text(angle = x.angle, hjust=x.hjust),
                     axis.text.y = element_text(angle = y.angle, hjust=y.hjust))
    if (!is.null(title)) {
      gg <- gg + ggtitle(title)+ theme(plot.title = element_text(hjust = 0.5))
    }
  }
  
  if (return.data) {
    df$contribution <- abs(df$contribution)
    df$contribution.scaled <- abs(df$contribution.scaled)
    return(list(signaling.contribution = df, gg.obj = gg))
  } else {
    return(gg)
  }
}


#' Identify all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway
#'
#' @param net,LR,DB object@net object@LR object@DB
#' @param signaling a char vector containing signaling pathway names for searching
#' @param enriched.only whether only return the identified enriched signaling genes in the database. Default = TRUE, returning the significantly enriched signaling interactions
#' @param thresh threshold of the p-value for determining significant interaction
#' @param do.group set `do.group = TRUE` when only showing enriched signaling based on cell group-level communication; set `do.group = FALSE` when only showing enriched signaling based on individual cell-level communication
#' @importFrom dplyr select
#'
#' @return a list: list(geneLR, pairLR.name.use)
extractEnrichedLR_internalX <- function(net, LR, DB, signaling, enriched.only = TRUE, 
                                        thresh = 0.05, do.group = TRUE) {
  pairLR <- searchPair(signaling = signaling, pairLR.use = LR$LRsig, 
                       key = "pathway_name", matching.exact = T, pair.only = T)
  pairLR.name.use = dplyr::select(DB$interaction[rownames(pairLR),],"interaction_name")
  
  if (enriched.only) {
    if (do.group) {
      pairLR.use.name <- dimnames(net$prob)[[3]]
      pairLR.name <- intersect(rownames(pairLR), pairLR.use.name)
      pairLR <- pairLR[pairLR.name, ]
      prob <- net$prob
      pval <- net$pval
      prob[pval > thresh] <- 0
      if ("LR.sig" %in% names(net) == FALSE) stop("Please run the `aggregateNet` function!", '\n')
      LR.sig <- net$LR.sig
    } else {
      #### pairLR.use.name <- dimnames(net$prob.cell)[[3]] 
      pairLR.use.name <- names(net$tmp$prob.cell)
      pairLR.name <- intersect(rownames(pairLR), pairLR.use.name)
      pairLR <- pairLR[pairLR.name, ]
      #### prob <- net$prob.cell  
      if ("LR.sig.cell" %in% names(net) == FALSE) stop("Please run the `aggregateNet` function!", '\n')
      LR.sig <- net$LR.sig.cell
    }
    
    pairLR.name.use <- intersect(pairLR.name, LR.sig)
    if (length(pairLR.name.use) == 0) {
      message(paste0('There is no significant communication of ', signaling))
    } else {
      pairLR <- pairLR[pairLR.name.use,]
    }
  }
  
  geneL <- unique(pairLR$ligand)
  geneR <- unique(pairLR$receptor)
  geneL <- extractGeneSubset(geneL, DB$complex, DB$geneInfo)
  geneR <- extractGeneSubset(geneR, DB$complex, DB$geneInfo)
  geneLR <- c(geneL, geneR)
  return(list(geneLR, pairLR.name.use))
}


#' Identify all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway
#'
#' @param object CellChat object
#' @param signaling a char vector containing signaling pathway names for searching
#' @param geneLR.return whether return the related signaling genes of enriched L-R pairs
#' @param enriched.only whether only return the identified enriched signaling genes in the database. 
#' Default = TRUE, returning the significantly enriched signaling interactions
#' @param thresh threshold of the p-value for determining significant interaction
#' @param do.group set `do.group = TRUE` when only showing enriched signaling based on cell group-level communication; 
#' set `do.group = FALSE` when only showing enriched signaling based on individual cell-level communication
#' @param geneInfo a dataframe with gene official symbol (there should be one column named `Symbol`)
#' @param complex_input signaling complex information from CellChatDB
#' @importFrom dplyr select
#'
#' @return The returned value depends on the input argument:
#' When `geneLR.return = FALSE`, it returns a data frame containing the significant interactions (L-R pairs)
#' When `geneLR.return = TRUE`, it returns a list, the first element is a data frame 
#' containing the significant interactions (L-R pairs), 
#' and the second is a vector containing the related signaling genes of enriched L-R pairs, 
#' which can be used for examining the gene expression pattern using the function \code{\link{plotGeneExpression}}
#'
#' @export
#'
extractEnrichedLRX <- function(object, signaling, geneLR.return = FALSE, 
                               enriched.only = TRUE, thresh = 0.05, do.group = TRUE, 
                               geneInfo = NULL, complex_input = NULL) {
  DB <- object@DB
  if (is.null(geneInfo)) {
    geneInfo = DB$geneInfo
  } else {
    DB$geneInfo = geneInfo
  }
  if (is.null(complex_input)) {
    complex_input = DB$complex
  } else {
    DB$complex = complex_input
  }
  pairLR.all <- c()
  geneLR.all <- c()
  net0 <- slot(object, "net")
  for (ii in 1:length(signaling)) {
    signaling.i <- signaling[ii]
    if (object@options$mode == "single") {
      net <- net0
      LR <- object@LR
      res <- extractEnrichedLR_internalX(net, LR, DB, signaling = signaling.i, 
                                         enriched.only = enriched.only, thresh = thresh, do.group = do.group)
    } else {
      geneLR.t <- c()
      pairLR.t <- c()
      for (i in 1:length(net0)) {
        net <- net0[[i]]
        LR <- object@LR[[i]]
        res.t <- extractEnrichedLR_internalX(net, LR, DB, signaling = signaling.i, 
                                             enriched.only = enriched.only, thresh = thresh, do.group = do.group)
        geneLR.t <- BiocGenerics::union(geneLR.t, as.character(res.t[[1]]))
        pairLR.t <- BiocGenerics::union(pairLR.t, as.character(res.t[[2]]))
      }
      res <- list(geneLR.t, pairLR.t)
    }
    geneLR.all <- c(geneLR.all, as.character(res[[1]]))
    pairLR.all <- c(pairLR.all, as.character(res[[2]]))
  }
  pairLR.all <- data.frame(interaction_name = pairLR.all, stringsAsFactors = FALSE)
  
  if (geneLR.return) {
    return(list(pairLR = pairLR.all, geneLR = geneLR.all))
  } else {
    return(pairLR.all)
  }
}



#' Compute the cell-cell communication field (outgoing/incoming)
#'
#' @param object SpatialCellChat object
#' @param slot.name the slot name of object that is used for analysis
#' @param signaling.name alternative signaling pathway name used for analysis. If NULL, all the signaling pathways will be used
#' @param top the cutoff for selecting top signalings
#' @param sparse whether to apply a sparse strategy for selecting top signalings
#'
#' @return A SpatialCellChat object
#' @export
computeCommunFieldX <- function(object, slot.name = "netP", signaling.name = NULL, 
                                top = 0.8, sparse = T ){
  if (is.null(methods::slot(object, slot.name)$tmp$prob.cell)) {
    if (slot.name == "net") {
      stop( cli.symbol(2), "Please run `computeCommunProb` to compute the 
      communication probability/strength between any interacting individual cells!" )
    } else if (slot.name == "netP") {
      stop( cli.symbol(2), "Please run `computeCommunProbPathway` to compute the 
            communication probability/strength between any interacting individual cells!")
    }
  } else {
    net <- methods::slot(object, slot.name)$tmp$prob.cell # a list
    
    if (slot.name == "net") {
      LRsig.use.idx <- object@net$tmp$LRsig.use.idx
      net <- net[LRsig.use.idx] # use LRsigs whose prob.sum > 0
    }
    
    signaling.use <- names(net)
    #### cell.names <- spatstat.sparse::dimnames.sparse3Darray(methods::slot(object, slot.name)$prob.cell)[[1]]
    cell.names <- colnames(object@data.signaling)
  }
  
  nC <- length(cell.names)
  
  if (!is.null(signaling.name)) {
    if (is.numeric(signaling.name)) signaling.name <- signaling.use[signaling.name]
    net <- net[signaling.name] # [] will return a list
  } else { signaling.name <- signaling.use }
  
  data.spatial <- object@images$coordinates
  # temp_coordinates = data.spatial
  # data.spatial[, 1] = temp_coordinates[, 2]
  # data.spatial[, 2] = temp_coordinates[, 1]
  data.spatial <- as.matrix(data.spatial)
  colnames(data.spatial) <- c("x_cent", "y_cent")
  
  nrun <- length(signaling.name)
  
  cf.all.outgoing_ <- my_future_lapply(
    X = 1:nrun,
    FUN = function(x) {
      
      net0 <- net[[x]] # [[]] will return one element(communication matrix)
      
      ### compute cf.outgoing ###
      idx.outgoing.cells <- which(Matrix::rowSums(net0) > 0)
      cf.outgoing <- sapply(
        X = idx.outgoing.cells,
        FUN = function(x) {
          res <- sort(net0[x, , drop = T], decreasing = TRUE, index.return = TRUE)
          if(sparse){
            idx.outgoing <- which(cumsum(res$x) < top * sum(res$x))
          } else {
            if (cumsum(res$x)[[2]] == cumsum(res$x)[[1]]) {
              idx.outgoing <- c(1)
            } else if (cumsum(res$x)[[3]] == cumsum(res$x)[[2]]) {
              idx.outgoing <- c(1, 2)
            } else {
              idx.outgoing <- which(cumsum(res$x) <= top * sum(res$x))
            }
          }
          idx.end <- res$ix[idx.outgoing]
          
          if (length(idx.end) > 0) {
            direction <-
              data.spatial[idx.end, , drop = F] - data.spatial[rep.int(x, times = length(idx.end)), , drop =F]
            prob.outgoing <- res$x[idx.outgoing]
            direction <- purrr::map_dfr(
              .x = 1:NROW(direction),
              .f = function(x) {
                norm.x <- norm(direction[x, , drop = T], type = "2")
                if (norm.x > 0) {
                  return(direction[x, , drop = T] * (prob.outgoing[x] / norm.x))
                } else if (norm.x == 0) {
                  return(direction[x, , drop = T])
                }
              }
            )
            direction <- apply(direction, 2L, FUN = sum)
            return(direction) # dim: 1,2
          } else {
            return(c(0, 0)) # dim: 1,2
          }
        }
      )
      
      idx.i <- rep(idx.outgoing.cells, times = 2)
      idx.j <-
        c(rep.int(1, times = length(idx.outgoing.cells)), rep.int(2, times = length(idx.outgoing.cells)))
      value.x <- c(cf.outgoing[1, , drop = T], cf.outgoing[2, , drop = T])
      
      cf.outgoing.sparse <- Matrix::sparseMatrix(
        i = idx.i,
        j = idx.j,
        x = value.x,
        dims = c(nC, 2),
        dimnames = list(cell.names, c("x_cent", "y_cent")),
        index1 = T
      ) # matrix (dim: NROW(net0) x 2 or nC x 2)
      
      return(cf.outgoing.sparse)
    },
    simplify = F,
    hint.message = "Computing outgoing..."
  )
  
  cf.all.outgoing <-
    my_as_sparse3Darray(cf.all.outgoing_, nonzero = T) # list=>3Darray
  rm(cf.all.outgoing_)
  gc()
  
  
  cf.all.incoming_ <- my_future_lapply(
    X = 1:nrun,
    FUN = function(x) {
      net0 <- net[[x]]
      
      ### compute cf.incoming ###
      idx.incoming.cells <- which(Matrix::colSums(net0) > 0)
      cf.incoming <- sapply(
        X = idx.incoming.cells,
        FUN = function(x) {
          res <- sort(net0[, x, drop = T], decreasing = TRUE, index.return = TRUE)
          if(sparse){
            idx.incoming <- which(cumsum(res$x) < top * sum(res$x))
          } else {
            if (cumsum(res$x)[[2]] == cumsum(res$x)[[1]]) {
              idx.incoming <- c(1)
            } else if (cumsum(res$x)[[3]] == cumsum(res$x)[[2]]) {
              idx.incoming <- c(1, 2)
            } else {
              idx.incoming <- which(cumsum(res$x) <= top * sum(res$x))
            }
          }
          idx.start <- res$ix[idx.incoming]
          
          if (length(idx.start) > 0) {
            direction <-
              # data.spatial[rep.int(x, times = length(idx.start)), , drop = F] - data.spatial[idx.start, , drop = F]
              data.spatial[rep.int(x, times = length(idx.start)), , drop = F] - data.spatial[idx.start, , drop = F]
            prob.incoming <- res$x[idx.incoming]
            direction <- purrr::map_dfr(
              .x = 1:NROW(direction),
              .f = function(x) {
                norm.x <- norm(direction[x, , drop = T], type = "2")
                if (norm.x > 0) {
                  return(direction[x, , drop = T] * (prob.incoming[x] / norm.x))
                } else if (norm.x == 0) {
                  return(direction[x, , drop = T])
                }
              }
            )
            direction <- apply(direction, 2L, FUN = sum)
            return(direction) # dim: 1,2
          } else {
            return(c(0, 0)) # dim: 1,2
          }
        }
      )
      
      idx.i <- rep(idx.incoming.cells, times = 2)
      idx.j <-
        c(rep.int(1, times = length(idx.incoming.cells)), rep.int(2, times = length(idx.incoming.cells)))
      value.x <- c(cf.incoming[1, , drop = T], cf.incoming[2, , drop = T])
      
      cf.incoming.sparse <- Matrix::sparseMatrix(
        i = idx.i,
        j = idx.j,
        x = value.x,
        dims = c(nC, 2),
        dimnames = list(cell.names, c("x_cent", "y_cent")),
        index1 = T
      ) # matrix (dim: NROW(net0) x 2 or nC x 2)
      return(cf.incoming.sparse)
    },
    simplify = F,
    hint.message = "Computing incoming..."
  )
  
  cf.all.incoming <-
    my_as_sparse3Darray(cf.all.incoming_, nonzero = T) # list=>3Darray
  rm(cf.all.incoming_)
  gc()
  
  
  dimnames(cf.all.outgoing) <-
    list(cell.names, colnames(data.spatial), signaling.name)
  dimnames(cf.all.incoming) <- dimnames(cf.all.outgoing)
  
  cf.all <-
    list(outgoing = cf.all.outgoing, incoming = cf.all.incoming)
  methods::slot(object, slot.name)$field <- cf.all
  return(object)
}


#' Relabel cell groups and re-run downstream inference
#'
#' @description
#' `relabelSpatialCellChatX` 是 [relabelSpatialCellChat()] 的迁移版本。
#' 逻辑完全不变，只是把四个下游函数切换为不再依赖 `net$prob.cell` 的 *X 版本。
#'
#' @inheritParams relabelSpatialCellChat
#' @export
relabelSpatialCellChatX <- function(
    object,
    labelSet,
    group.by = NULL,
    nboot = 100,
    min.cells.sr = 10,
    min.percent = 0.1,
    do.permutation = T,
    min.cells = 10,
    thresh = 0.05
){
  if(is.null(group.by)){
    group.lab <- as.character(object@idents)
    group.level <- levels(object@idents)
  } else {
    if (!(group.by %in% colnames(object@meta))) {
      stop("The 'group.by' is not a column name in the `meta` slot, which will be used for cell grouping.")
    } else {
      group.lab <- as.character(object@meta[[group.by]])
      group.level <- unique(group.lab)
    }
  }
  
  if(is.null(names(labelSet))){
    stop("Please provide the names of cell groups which will be dropped
    by setting the `names` of labelSet.\n")
  }
  group.drop <- names(labelSet)
  if(!all(group.drop%in%group.level)) {
    stop("Please check the `labelSet`, make sure the cell groups
    to be dropped are in the `group.by` column.\n")
  }
  
  tmp.df <- data.frame(
    pre = group.lab,
    new = labelSet[group.lab]
  )
  tmp.df <- tmp.df %>% mutate(
    new = if_else(is.na(new), pre, new)
  )
  
  object@meta[["new.ident"]] <- tmp.df[["new"]]
  
  if(!is.null(object@net[["tmp"]]$cell.type.decomposition)){
    cell.type.decomposition <- as.matrix(object@net[["tmp"]]$cell.type.decomposition)
    cell.type.decomposition.left <- cell.type.decomposition[, base::setdiff(group.level, group.drop), drop = F]
    
    cell.type.decomposition.new <- lapply(
      X = unique(labelSet),
      FUN = function(i){
        group.drop.i <- names(labelSet[labelSet == i])
        cell.type.decomposition.i <-
          Matrix::rowSums(cell.type.decomposition[, group.drop.i, drop = F])
        return(cell.type.decomposition.i)
      }
    )
    names(cell.type.decomposition.new) <- unique(labelSet)
    cell.type.decomposition.new <- as.data.frame(cell.type.decomposition.new)
    cell.type.decomposition.new <- cbind(cell.type.decomposition.left, as.matrix(cell.type.decomposition.new))
    object@meta[["new.ident"]] <- factor(object@meta[["new.ident"]], levels = colnames(cell.type.decomposition.new))
    object@idents <- object@meta[["new.ident"]]
    
    object <- computeAvgCommunProb_VisiumX(          # was: computeAvgCommunProb_Visium
      object,
      cell.type.decomposition = cell.type.decomposition.new,
      nboot = nboot,
      do.permutation = do.permutation
    )
  } else {
    object@meta[["new.ident"]] <- factor(object@meta[["new.ident"]], levels = unique(object@meta[["new.ident"]]))
    object@idents <- object@meta[["new.ident"]]
    object <- computeAvgCommunProbX(                 # was: computeAvgCommunProb
      object,
      nboot = nboot,
      min.cells.sr = min.cells.sr,
      min.percent = min.percent,
      do.permutation = do.permutation
    )
  }
  
  object <- filterCommunicationX(                    # was: filterCommunication
    object,
    min.cells = min.cells,
    min.links = NULL,
    min.cells.sr = NULL
  )
  
  object <- computeCommunProbPathwayX(               # was: computeCommunProbPathway
    object,
    thresh = thresh
  )
  
  object <- aggregateNetX(object)                    # was: aggregateNet
  
  return(object)
}


#' Update cell cluster labels and re-run pathway-level aggregation
#'
#' @description
#' `updateClusterLabelsX` 是 [updateClusterLabels()] 的迁移版本。
#' 逻辑完全不变，只把三个下游函数切换为 *X 版本。
#'
#' @inheritParams updateClusterLabels
#' @export
updateClusterLabelsX <- function(object, old.cluster.name = NULL, new.cluster.name = NULL,
                                 new.order = NULL, new.cluster.metaname = "new.labels") {
  if (is.null(old.cluster.name)) {
    old.cluster.name <- levels(object@idents)
  }
  if (new.cluster.metaname %in% colnames(object@meta)) {
    stop("Please define another `new.cluster.metaname` as it exists in `colnames(object@meta)`!")
  }
  if (!is.null(new.cluster.name)) {
    labels.new <- plyr::mapvalues(object@idents, from = old.cluster.name, to = new.cluster.name)
    object@meta[[new.cluster.metaname]] <- labels.new
    object <- setIdent(object, ident.use = new.cluster.metaname, display.warning = FALSE)
  } else {
    new.cluster.metaname <- NULL
    cat("Only reorder cell groups but do not rename cell groups!")
  }
  
  if (!is.null(new.order)) {
    object <- setIdent(object, ident.use = new.cluster.metaname, levels = new.order, display.warning = FALSE)
  }
  message("We now re-run `computeCommunProbPathwayX`, `aggregateNetX`, and `netAnalysis_computeCentralityX`...")
  object <- computeCommunProbPathwayX(object)        # was: computeCommunProbPathway
  ## calculate the aggregated network by counting the number of links or summarizing the communication probability
  object <- aggregateNetX(object)                    # was: aggregateNet
  # network importance analysis
  object <- netAnalysis_computeCentralityX(object, slot.name = "netP")  # was: netAnalysis_computeCentrality
  return(object)
}
