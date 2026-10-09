#' @description 用my_as_sparse3Darray将list形式的 转换成3D稀疏array, 并存入net槽的prob.cell
#' @export
net3Darray <- function(object, use.raw = FALSE ) {
  Prob.cell <- my_as_sparse3Darray(object@net$tmp$prob.cell)  
  if (use.raw) { cell <- colnames(object@data.raw) } else { cell <- colnames(object@data)}  
  dimnames(Prob.cell) <- list(cell, cell, names(object@net$tmp$prob.cell) ) 
  object@net$prob.cell <- Prob.cell
  return(object)
}

#' 把分批次算出来的细胞级通信中间结果（每批一个 net$tmp）拼回主对象
#' @description 
#' @export
combinePathway <- function(chat=NULL, tmps=NULL ) {
  for ( i in seq_along(tmps) ) {
    chat@net$tmp$prob.cell <- c(chat@net$tmp$prob.cell, tmps[[i]]$prob.cell)
    chat@net$tmp$Lavg <- rbind(chat@net$tmp$Lavg, tmps[[i]]$Lavg)
    chat@net$tmp$Ravg <- rbind(chat@net$tmp$Ravg, tmps[[i]]$Ravg)
    chat@net$tmp$ligand <- c(chat@net$tmp$ligand, tmps[[i]]$ligand)
    chat@net$tmp$receptor <- c(chat@net$tmp$receptor, tmps[[i]]$receptor)
  }
  orderX <- match(chat@LR$LRsig$interaction_name, names(chat@net$tmp$prob.cell))
  chat@net$tmp$prob.cell <- chat@net$tmp$prob.cell[orderX]
  chat@net$tmp$Lavg <-  chat@net$tmp$Lavg[orderX, , drop = FALSE]
  chat@net$tmp$Ravg <-  chat@net$tmp$Ravg[orderX, , drop = FALSE]
  chat@net$tmp$ligand <-  chat@net$tmp$ligand[orderX]
  chat@net$tmp$receptor <-  chat@net$tmp$receptor[orderX]
  return(chat)
}


#' computeCommunProbX
#'
#' @description
#' Compute the communication probability/strength between any interacting individual cells
#'
#' @param object CellChat object
#' @param LR.use A subset of ligand-receptor interactions used in inferring communication network
#' @param raw.use Whether use the raw data (i.e., `object@data.signaling`) or the projected data (i.e., `object@data.project`).
#' Set raw.use = FALSE to use the projected data when analyzing single-cell data with shallow sequencing depth because the projected data could help to reduce the dropout effects of signaling genes, in particular for possible zero expression of subunits of ligands/receptors.
#' @param Kh Parameter in Hill function
#' @param n Parameter in Hill function
#' @param distance.use Whether to use distance constraints to compute communication probability. Setting `distance.use = TRUE` indicates that the cell-cell communication probability is inversely proportional to the computed distance.
#' @param tol NULL or Numeric. set a distance tolerance when computing cell-cell distances and contact adjacent matrix. By default `tol = NULL` means distance tolerance equals to spot.size/2 or cell.diameter/2.
#' @param interaction.range The maximum interaction/diffusion length of ligands (Unit: microns). This hard threshold is used to filter out the connections between spatially distant individual cells
#' @param scale.distance A scale or normalization factor for the spatial distances when setting `distance.use = TRUE`. For example, scale.distance equals 1, 0.1, 0.01, 0.001, 0.11, or 0.011. We choose this values such that the minimum value of the scaled distances is in [1,2]. This value is not necessary when setting `distance.use = FALSE`.
#' @param use.AGAN Boolean. Whether to take agonist (AG) and antagonist (AN) into consideration when calculating the intercellular communication probability
#' @param contact.range Numeric. The interaction range (Unit: microns) to restrict the contact-dependent signaling.
#' For spatial transcriptomics in a single-cell resolution, `contact.range` is approximately equal to the estimated cell diameter (i.e., the cell center-to-center distance), which means that contact-dependent and juxtacrine signaling can only happens when the two cells are contact to each other.
#' Typically, `contact.range = 10`, which is a typical human cell size. However, for low-resolution spatial data such as 10X visium, it should be the cell center-to-center distance (i.e., `contact.range = 100` for visium data).  The function \link{computeCellDistance} can compute the center-to-center distance.
#' @param contact.dependent Boolean. Whether determining spatially proximal cell groups based on the `contact.range`. By default `contact.dependent = TRUE` when inferring contact-dependent and juxtacrine signaling (including ECM-Receptor and Cell-Cell Contact signaling classified in CellChatDB$interaction$annotation).
#' If only focusing on `Secreted Signaling`, the `contact.dependent` will be automatically set as FALSE except for `contact.dependent.forced = TRUE`.
#' @param contact.dependent.forced Boolean. Whether forcing to determine spatially proximal cell regions based on the `contact.range` for all cells/spots in the CellChat Object.
#' Users can set `contact.dependent.forced = TRUE` to turn `Secreted Signaling` interactions into a contact manner.
#'
#' @return CellChat object
#' @export
#'
#' @examples
computeCommunProbX <- function(object, LR.use = NULL, raw.use = TRUE, 
                               Kh = 0.5, n = 1, 
                               distance.use = TRUE, tol = NULL, 
                               interaction.range = 250, 
                               use.AGAN = T, scale.distance = 0.01,
                               contact.dependent = TRUE, contact.range = 10,
                               contact.dependent.forced = FALSE){
  
  #  选择 data.signaling 或者 data.project
  if (raw.use) {
    data <- object@data.signaling
    # scale the elements
    data@x <- data@x/max(data@x)
    data.use <- as.matrix(data)
  } else {
    data <- object@data.project
    # scale
    data.use <- data/max(data)
  }
  
  # 提取指定的LR对  
  if (is.null(LR.use)) {
    pairLR.use <- object@LR$LRsig
  } else {
    if (length(unique(LR.use$annotation)) > 1) { 
      LR.use$annotation <- factor(LR.use$annotation, 
                                  levels = c( "Secreted Signaling", "ECM-Receptor", 
                                              "Non-protein Signaling", "Cell-Cell Contact" ) ) 
      LR.use <- LR.use[order(LR.use$annotation), , drop = FALSE]
      LR.use$annotation <- as.character(LR.use$annotation) 
    }
    pairLR.use <- LR.use
  }
  
  complex_input <- object@DB$complex  # 复合物
  cofactor_input <- object@DB$cofactor # cofactor
  
  ptm = Sys.time()
  
  pairLRsig <- pairLR.use
  group <- object@idents
  geneL <- as.character(pairLRsig$ligand) # ligand
  geneR <- as.character(pairLRsig$receptor) # receptor
  nLR <- nrow(pairLRsig) 
  numCluster <- nlevels(group)
  if (numCluster != length(unique(group))) {
    stop(cli.symbol(2),"Please check `unique(object@idents)` and ensure that the factor levels are correct!\n 
         You may need to drop unused levels using 'droplevels' function. e.g.,\n 
         `meta$labels = droplevels(meta$labels, exclude = setdiff(levels(meta$labels),unique(meta$labels)))`")
  }
  
  nC <- ncol(data.use)
  
  # working on spatial transcriptomic data and preferring to infer interactions between individual cell
  cat(cli.symbol(),"Analyzing spatial transcriptomic data and preferring to 
      infer interactions between individual cells...\n")
  
  
  data.spatial <- BiocGenerics::as.data.frame(object@images$coordinates)
  
  ### previous ###
  # spot.size.fullres <- object@images$scale.factors$spot
  # spot.size <- object@images$scale.factors$spot.diameter
  
  ratio <- object@images$spatial.factors[["ratio"]]
  if(is.null(tol)) tol <- object@images$spatial.factors[["tol"]] else NULL
  
  # 计算cell-to-cell distances
  # res is a list object, containing `d.spatial` matrix and `adj.contact` matrix!
  res <- computeCellDistance(coordinates = data.spatial, ratio = ratio,
                             interaction.range = interaction.range,
                             contact.range = contact.range, tol = tol)
  # long-range distance
  d.spatial <- res$d.spatial
  # short-range distance adjacent matrix for contact-dependent and juxtacrine signaling
  adj.contact <- res$adj.contact
  gc()
  
  if (distance.use) {
    cat(paste0(cli.symbol(),"Run CellChat on spatial transcriptomic data 
    using distances as constraints <<< [", Sys.time(), "]\n"))
    d.spatial@x <- d.spatial@x * scale.distance
    d.min <- min(d.spatial@x, na.rm = F) # d.spatial@x has no NA
    if (d.min < 1) {
      cat(cli.symbol(),"The suggested minimum value of scaled distances is in [1,2],
          and the calculated value here is ", d.min,"\n")
      stop(cli.symbol(2),"Please increase the value of `scale.distance` and 
           use a value that is slighly smaller than ", format(1/d.min, digits = 2) ,"\n")
    }
    P.spatial <- createPspatialFrom_dspatial(d.spatial,distance.use = T)
    d.spatial@x <- d.spatial@x / scale.distance
  }  else {
    cat(paste0(cli.symbol(),"Run CellChat on transcriptomic imaging data 
               without distances as constraints <<< [", Sys.time(), "]\n"))
    P.spatial <- createPspatialFrom_dspatial(d.spatial,distance.use = F)
    
  }
  rm(d.spatial);gc()
  
  # set a flag for contact.dependent signaling
  all.contact.dependent <- FALSE
  all.diffusible <- FALSE
  if (contact.dependent.forced == TRUE) {
    # cat(cli.symbol(),"Run with `contact.dependent.forced = T` \n")
    cat(cli.symbol(),"Force to run CellChat in a `contact-dependent` manner for 
        all L-R pairs including secreted signaling.\n")
    P.spatial <- P.spatial * adj.contact
    nLR1 <- nLR
    all.contact.dependent <- TRUE
  } else{ # contact.dependent.forced == F
    # cat(cli.symbol(),"Run with `contact.dependent.forced = F` \n")
    if (contact.dependent == TRUE && length(unique(pairLRsig$annotation))>0 ) {
      if (all(unique(pairLRsig$annotation) == c("Cell-Cell Contact") )) {
        # all interactions in `pairLRsig` are contact-dependent signaling
        cat(cli.symbol(),"All the input L-R pairs are `Cell-Cell Contact` signaling. 
            Run CellChat in a contact-dependent manner. \n")
        P.spatial <- P.spatial * adj.contact
        nLR1 <- nLR
        all.contact.dependent <- TRUE
      } else if (all(unique(pairLRsig$annotation) %in% c("Secreted Signaling", "ECM-Receptor", "Non-protein Signaling"))) {
        # all interactions in `pairLRsig` are not contact-dependent signaling
        cat(cli.symbol(),"Molecules of all the input L-R pairs are diffusible. 
            Run CellChat in a diffusion manner based on the `interaction.range`.\n")
        nLR1 <- nLR
        all.diffusible <- TRUE
      } else {
        # Interactions in `pairLRsig` have both contact-dependent signaling and secreted signaling
        cat(cli.symbol(),"The input L-R pairs have both secreted signaling and contact-dependent signaling. 
            Run CellChat in a contact-dependent manner for `Cell-Cell Contact` signaling, 
            and in a diffusion manner based on the `interaction.range` for other L-R pairs. \n")
        nLR1 <- max(which(pairLRsig$annotation %in% c("Secreted Signaling", "ECM-Receptor", "Non-protein Signaling")))
      }
    } else { # contact.dependent == F or `object@LR$LRsig` does not have `annotation` column, take all interactions as `Secreted Signaling` interactions
      cat(cli.symbol(),"Run CellChat in a diffusion manner based on the `interaction.range` for all L-R pairs. \n")
      cat(cli.symbol(3),"Set`contact.dependent = TRUE` if preferring a contact-dependent manner for `Cell-Cell Contact` signaling. \n")
      nLR1 <- nLR
    }
  }
  
  
  # compute the expression of ligand or receptor
  dataLavg <- computeExpr_LR(geneL, data.use, complex_input)
  dataRavg <- computeExpr_LR(geneR, data.use, complex_input)
  
  # 
  # take account into the effect of co-activation and co-inhibition receptors
  dataRavg.co.A.receptor <- computeExpr_coreceptor(cofactor_input, data.use, pairLRsig, type = "A")
  dataRavg.co.I.receptor <- computeExpr_coreceptor(cofactor_input, data.use, pairLRsig, type = "I")
  dataRavg <- dataRavg * dataRavg.co.A.receptor/dataRavg.co.I.receptor
  rm(dataRavg.co.A.receptor,dataRavg.co.I.receptor);gc();
  
  
  # 激动剂和拮抗剂
  # compute the expression of agonist and antagonist
  # index.agonist <- which(!is.na(pairLRsig$agonist) & pairLRsig$agonist != "")
  # index.antagonist <- which(!is.na(pairLRsig$antagonist) & pairLRsig$antagonist != "")
  
  #  将 稀疏矩阵SparseMat 与 一个稠密向量 DenseVec 进行两次逐元素相乘
  # 内存效率高：只修改非零元素，保持稀疏结构
  # 速度快：直接在 @x 槽位上操作，避免复制整个矩阵
  # chat@data.signaling@i 非0元素的列索引
  myElementwiseProduct_fast <- function(SparseMat, DenseVec) {
    SparseMat@x <- SparseMat@x * DenseVec[SparseMat@i + 1] *   # 非0元素的列索引, +1 从0-base转成1-base
      DenseVec[rep(seq_len(ncol(SparseMat)) - 1, diff(SparseMat@p)) + 1]   # 非0元素的行索引
    SparseMat
  }
  
  # Compute the communication probability/strength between any interacting individual cells for each LR pair
  sp <- summary(P.spatial) 
  # 创建一个空的稀疏矩阵模板，其稀疏结构（非零元素的位置）与 稀疏矩阵 P.spatial 完全相同，但所有值初始化为 0
  # template <- sparseMatrix(i = sp$i, j = sp$j,  x = numeric(length(sp$i)), dims = dim(P.spatial) )
  options(future.stdout = FALSE)
  
  print('准备处理每对LR')
  # 使用 future.apply 进行并行,逐一处理每对LR
  Prob.cell_ <- with_progress({
    
    p <- progressr::progressor(along = seq_len(nLR))
    
    future.apply::future_lapply(
      X = seq_len(nLR), # 遍历所有配体-受体对 (nLR 个)
      future.seed = TRUE,
      
      FUN = function(i) {
        
        # accès direct (rapide)
        x_ <- dataLavg[i, ]
        y_ <- dataRavg[i, ]
        
        # calcul vectorisé
        # dataLR <- x_[sp$i + 1] * y_[sp$j + 1]
        dataLR <- x_[sp$i] * y_[sp$j] # sp是1base的
        dataLR <- dataLR^n / (Kh^n + dataLR^n)
        
        # sparse rapide via template
        # P1_Pspatial <- template
        # P1_Pspatial@x <- dataLR * sp$x
        P1_Pspatial <- P.spatial  # 直接复制，保持完全相同的结构和顺序
        P1_Pspatial@x <- dataLR * sp$x  # 替换值
        
        # contact
        if (i > nLR1) { P1_Pspatial@x <- P1_Pspatial@x * adj.contact@x }
        
        # cas simple
        if (!use.AGAN || all(P1_Pspatial@x == 0, na.rm = TRUE)) { 
          result <- P1_Pspatial 
        } else {
          # 计算激动剂效应
          data.agonist <- computeExpr_agonist(data.use = data.use, pairLRsig, cofactor_input, 
                                              index.agonist = i, Kh = Kh, n = n )
          P_ <- myElementwiseProduct_fast(P1_Pspatial, data.agonist)
          
          # 计算拮抗剂效应
          data.antagonist <- computeExpr_antagonist(data.use = data.use, pairLRsig, cofactor_input, 
                                                    index.antagonist = i, Kh = Kh, n = n )
          result <- myElementwiseProduct_fast(P_, data.antagonist)
        }
        
        result@x[abs(result@x) < 0.0001] <- 0
        result@x[is.na(result@x)] <- 0
        result <- Matrix::drop0(result) # 移除 0 值，增加稀疏度
        
        p(sprintf("i=%d", i)) # 更新进度条
        result # 返回结果
      }
    )
  })
  
  print('All LR pair is done')
  # 每个 LR pair 产出的 P1_Pspatial 是一个 N×N 稀疏矩阵，要全部驻留内存才能合成 3D array。
  # bind the Prob.cell `list` => a `sparse3Darray`
  # then the Prob.cell's shape will be (nC,nC,nLR)
  
  # Prob.cell <- my_as_sparse3Darray(Prob.cell_)
  # cat(cli.symbol(),"The number of cells and L-R pairs in Dim(Prob.cell):",dim(Prob.cell),"\n")
  
  # set `Prob.cell`'s names
  # dimnames(Prob.cell) <- list(colnames(data.use), colnames(data.use), rownames(pairLRsig))
  
  names(Prob.cell_) <- rownames(pairLRsig)
  
  # Tmp <- list(prob.cell = Prob.cell_,Lavg=dataLavg,Ravg=dataRavg) # !important, for parallel iteration
  # net <- list(prob.cell = Prob.cell, tmp = Tmp)
  
  # 关键优化 4：拒绝三份拷贝！只保留 CellChat 必需的最少变量
  # 将 Tmp 中的 prob.cell 指向已经给出的对象或直接赋 NULL，避免重构数据副本
  # 记录受体和配体的名字, 为Lvag和Ravg的rownames, 为combine pathway batch作准备
  Tmp <- list(prob.cell = Prob.cell_, Lavg = dataLavg, Ravg = dataRavg, ligand=geneL, receptor=geneR)   
  net <- list(prob.cell = NULL, tmp = Tmp)  
  # 释放内存垃圾
  rm(Prob.cell_)
  gc()
  
  execution.time = Sys.time() - ptm
  object@options$run.time <- as.numeric(execution.time,
                                        units = "secs")
  object@images[["result.computeCellDistance"]] <- res # 保存细胞距离
  object@options$parameter <- list(
    raw.use = raw.use, ratio = ratio,tol = tol, Kh = Kh, n = n, nLR = nLR, nLR1 = nLR1,
    scale.distance = scale.distance, use.AGAN = use.AGAN, distance.use = distance.use,
    interaction.range = interaction.range, contact.dependent = contact.dependent,
    contact.range = contact.range, contact.dependent.forced = contact.dependent.forced,
    all.contact.dependent = all.contact.dependent, all.diffusible = all.diffusible )
  
  object@net <- net
  cat(paste0(cli.symbol(symbol = "success")," CellChat inference is done. 
             Parameter values are stored in `object@options$parameter` <<< [", Sys.time(), "]", "\n"))
  
  return(object)
}


#' @title filterProbabilityX
#' @description
#' Filter out statistically non-significant communication probability at the level of individual cells after running \link{computeCommunProb}
#'
#' @param object CellChat object
#' @param nboot Numeric. The number of bootstrap samples, 100 by default.
#' @param seed.use Integer. The random seed used when taking a sample from the communication probabilities
#' @param thresh Numeric. The threshold for defining significant individual cell-cell communication at (1-thresh) of a shuffled distribution of each L-R pair.
#'
#' @return CellChat object
#' @export
filterProbabilityX <- function(object, nboot = 100, seed.use = 666L, thresh = 0.05 ){
  quantile.prob <- 1-thresh
  if(quantile.prob==0){
    cat(cli.symbol(1),"Do not filter any CCC probability!")
    return(object)
  } else { # quantile.prob<1
    if (is.null(object@net$tmp$prob.cell)) {
      stop(
        cli.symbol(2),
        "Please run `computeCommunProb` to compute the communication probability/strength 
        between any interacting individual cells! "
      )
    } else {
      cat(paste0(cli.symbol(), "Filter out non-significant communication with a probability quantile being ",
                 quantile.prob, " for each L-R pair... \n"))
      prob.cell_ <- object@net$tmp$prob.cell
      pair.LR.use <- names(prob.cell_)
      #### cell.names <- spatstat.sparse::dimnames.sparse3Darray(object@net$prob.cell)[[1]]
      cell.names <- colnames(object@data.signaling)
    }
    
    d.spatial <- object@images$result.computeCellDistance$d.spatial
    Matrix::diag(d.spatial) <- 1
    adj.contact <- object@images$result.computeCellDistance$adj.contact
    
    nLR <- object@options$parameter$nLR
    nLR1 <- object@options$parameter$nLR1
    nC <- NROW(d.spatial)
    
    set.seed(seed.use)
    if (object@options$parameter[["all.contact.dependent"]] == TRUE) {
      d.spatial <- adj.contact
    }
    
    # dim(permutation) = nboot x nLR
    permutation <- replicate(nLR, base::sample(x = 1:nC, size = nboot, replace = F))
    
    prob.cell_ <- my_future_lapply(X = 1:nLR, FUN = function(i) {
      
      if (i <= nLR1) { d_spatial <- d.spatial
      } else { d_spatial <- adj.contact }
      
      sample.cells <- permutation[ ,i,drop=T]
      Prob.cell.i <- prob.cell_[[i]]
      
      sample.prob.cell.i <- purrr::map(.x = sample.cells,
                                       .f = function(cell.index) {
                                         prob.index <- which(d_spatial[cell.index, ,drop = T] > 0)
                                         nboot.prob.cell.i <- Prob.cell.i[cell.index,prob.index, drop = T] # get a dense vec
                                         return(nboot.prob.cell.i)
                                       }) %>% unlist()
      nboot.quantile <- quantile(sample.prob.cell.i, probs = quantile.prob)
      gc()
      
      if(nboot.quantile == 0){
        # sparse enough, return directly
        return(Prob.cell.i)
      } else if (nboot.quantile > 0) {
        # filter `Prob.cell.i` to make it sparse enough
        Prob.cell.i <- scMatrixTruncation(Prob.cell.i,cutoff = nboot.quantile,remain.cutoff.v = T,repr = "C")
      }
      
      return(Prob.cell.i)
    }, simplify = F, hint.message = "filtering...")
    
    # 用net3Darray,分开做, 减少内存消耗
    # prob.cell <- my_as_sparse3Darray(prob.cell_)
    # dimnames(prob.cell) <- list(cell.names, cell.names, pair.LR.use)
    # object@net$prob.cell <- prob.cell
    
    names(prob.cell_) <- pair.LR.use
    object@net$tmp$prob.cell <- prob.cell_
    
    cat(cli.symbol(1), "Filtering is done.\n")
    return(object)
  } 
}


#' Filter cell-cell communication if there are only few number of cells in certain cell groups or only few interactions
#'
#' @param object CellChat object
#' @param min.cells the minimum number of cells required in each cell group for filtering cell group-level communication
#' @param min.links the minimum number of links/interactions required in the ligand-receptor pair for filtering individual cell-level communication
#' @param min.cells.sr the minimum number of cells required as senders or receivers for filtering individual cell-level communication
#' @return CellChat object with an updated slot net
#' @export
#'
filterCommunicationX <- function(object, min.cells = 10, min.links = 5, min.cells.sr = 5) {
  
  if (!is.null(min.cells)) {
    message("Filter cell-group level communication...",'\n')
    net <- object@net
    cell.excludes <- which(as.numeric(table(object@idents)) < min.cells)
    if (length(cell.excludes) > 0) {
      cat(cli.symbol(),"The cell-cell communication related with the following cell groups 
          are excluded due to the few number of cells: ", levels(object@idents)[cell.excludes],'\n')
      # dim(net$prob) = nCellGroup x nCellGroup x nPairLRsig
      net$prob[cell.excludes,,] <- 0
      net$prob[,cell.excludes,] <- 0
      if (!is.null(net$pval)) {
        net$pval[net$prob == 0] <- 1
      }
      object@net <- net
    }
    rm(net)
    gc()
  }
  
  #### if ("prob.cell" %in% names(object@net)) {
  if (!is.null(object@net$tmp$prob.cell)) {
    net <- object@net
    #### prob.cell <- net$prob.cell
    prob.cell_ <- net$tmp$prob.cell # a list
    
    if (!is.null(min.links) | !is.null(min.cells.sr)) {
      message("Filter individual cell-level communication...",'\n')
      
      prob.sum <- purrr::map_dbl(.x = prob.cell_, 
                                 .f = function(Mat){return(length(Mat@x))} )
      
      ##### dimArr <- dim(prob.cell)
      nCells <- nrow(prob.cell_[[1]])
      nLRs <- length(prob.cell_)
      LRnames <- names(prob.cell_)
      dimArr <- c(nCells, nCells, nLRs)
      
      # define a allzero matrix (CsparseMatrix)
      AllzeroMat <- Matrix::sparseMatrix(
        i = integer(0),
        j = integer(0),
        x = numeric(0),
        repr = "C", # default repr in CellChat
        dims = dimArr[c(1, 2)],
        # dimnames = dns[c(1,2)] # too large, not use!
        dimnames = list(NULL,NULL),
        index1 = T # i and j are interpreted as 1-based indices, following the R convention
      )
      
      # filter communication according to min.links
      gc()
      if (!is.null(min.links)) {
        cat(cli.symbol(),"Filter communication according to min.links...\n")
        idx.signaling.excludes <- which((prob.sum < min.links) & (prob.sum > 0))
        if (length(idx.signaling.excludes) > 0) {
          cat("The cell-cell communication related with #", length(idx.signaling.excludes),
              'L-R pairs are excluded due to the few number of interactions.','\n')
          
          # prob.cell[,,idx.signaling.excludes] <- 0
          pb <- utils::txtProgressBar(min = 0, max = length(idx.signaling.excludes), style = 3, file = stderr(), width = 80);i=0
          for (x in idx.signaling.excludes) {
            prob.cell_[[x]] <- AllzeroMat
            utils::setTxtProgressBar(pb = pb, value = (i=i+1))
          } # forloop
          close(con = pb)
          
        }
      }
      
      # filter communication according to min.cell.sr
      gc()
      if (!is.null(min.cells.sr)) {
        cat(cli.symbol(),"Filter communication according to min.cell.sr... \n")
        ##### pathways0 <- dimnames(prob.cell)[[3]] # L-R pairs' names
        pathways0 <- LRnames
        if (is.null(min.links)) {
          pathways <- pathways0[prob.sum > 0]
        } else {
          pathways <- pathways0[prob.sum >= max(1, min.links)] # Prevent `min.links` from being smaller than 1
        }
        
        if (length(pathways)<1){
          NULL # not filter
        } else {
          pathways.remove <- pbapply::pblapply(
            X = seq_len(length(pathways)),
            FUN = function(x) {
              prob.cell.i <- prob.cell_[[ pathways[[x]] ]] > 0 # `prob.cell.i` is a logical sparse matrix
              if ((sum(Matrix::rowSums(prob.cell.i) > 0) < min.cells.sr) | (sum(Matrix::colSums(prob.cell.i) > 0) < min.cells.sr)) {
                gc()
                return(pathways[[x]])
              }
            }
          )
          pathways.remove <- unlist(pathways.remove) # vec2vec
          pathways.remove.idx <- which(pathways0 %in% pathways.remove)
          
          if (length(pathways.remove) > 0) {
            cat(
              "The cell-cell communication related with #",
              length(pathways.remove),
              'L-R pairs are excluded due to the few number of sending/receiving cells.',
              '\n'
            )
            
            pb <- utils::txtProgressBar(min = 0, max = length(pathways.remove.idx), style = 3, file = stderr(),width = 80);i=0
            for (x in pathways.remove.idx) {
              prob.cell_[[x]] <- AllzeroMat
              utils::setTxtProgressBar(pb = pb, value = (i=i+1))
            } # forloop
            close(con = pb)
          }
        }
      }
      
      # update the obj
      net$tmp$prob.cell <- prob.cell_ # a list
      
      #### dns <- dimnames(prob.cell) # dimnames
      #### net$prob.cell <- my_as_sparse3Darray(prob.cell_)
      #### dimnames(net$prob.cell) <- dns
      
      object@net <- net
      cat(paste0(cli.symbol(1), 'Filtering cell-cell communication is done.<<< [', Sys.time(),']', '\n'))
      
    } # !is.null(min.links) | !is.null(min.cells.sr)
  } else {
    stop( cli.symbol(2), "Please run `computeCommunProb` to compute the 
    communication probability/strengthbetween any interacting individual cells!")
  }
  return(object)
}


#' Compute average communication probabilities of pairwise cell groups for one particular ligand-receptor pair/signaling pathway
#'
#' @param prob a matrix of communication probabilities for pairwise individual cells for one particular ligand-receptor pair/signaling pathway
#' @param dataLR a nCell*2 data matrix of a given pair of ligand-receptor.
#' @param group Character vector. Cell group information used for computing averaged communication probabilities
#' @param min.percent Numeric from 0 to 1. Minimum percentage of expressed ligands or receptors per cell group 
#' to require for computing the group-level signaling. Default is 0.1.
#' @param min.cells.sr Integer greater than 0. Minimum number of cells required as senders or receivers per cell group 
#' for computing the group-level signaling. Default is 5.
#'
#' @return Returns a matrix containing the interaction weights between any two cell groups.
#' @export
computeAvgCommunProb_LR_AvgX <- function (prob, group, dataLR = NULL, min.percent = 0.1, 
                                          min.cells.sr = 5 ){
  if (!is.factor(group)) group <- factor(group)
  gi <- as.integer(group); G <- nlevels(group); N <- nrow(prob)
  Z <- Matrix::sparseMatrix(i = seq_len(N), j = gi, x = 1, dims = c(N, G))
  
  d01 <- 1 * (dataLR > 0)
  m   <- cbind(tapply(d01[, 1], gi, mean), tapply(d01[, 2], gi, mean))

  dpe <- 1 * (format(as.data.frame(m), digits = 1) >= min.percent)
  Pp  <- Matrix::crossprod(matrix(dpe[, 1], nrow = 1), matrix(dpe[, 2], nrow = 1))
  
  if (sum(Pp) == 0) {
    Prob.avg <- Pp
  } else {
    Prob.avg <- Matrix::crossprod(Z, prob) %*% Z
    pb <- prob; pb@x <- rep.int(1, length(pb@x))
    Prob.scale.factor <- Matrix::crossprod(Z, pb) %*% Z
    Prob.avg <- Prob.avg / Prob.scale.factor
    Prob.avg[is.nan(Prob.avg)] <- 0
    
    cs <- rowsum(cbind(Matrix::rowSums(pb), Matrix::colSums(pb)), group, reorder = TRUE)
    cs <- 1 * (cs >= min.cells.sr)
    cells.sr <- Matrix::crossprod(matrix(cs[, 1], nrow = 1), matrix(cs[, 2], nrow = 1))
    
    Prob.avg <- Prob.avg * Pp * cells.sr
  }
  dimnames(Prob.avg) <- list(levels(group), levels(group))
  as.matrix(Prob.avg)
}

#' Compute average communication probabilities of pairwise cell groups for 
#' one particular ligand-receptor pair/signaling pathway
#'
#' @param prob a matrix of communication probabilities for pairwise individual cells 
#' for one particular ligand-receptor pair/signaling pathway
#' @param dataLR a nCell*2 data matrix of a given pair of ligand-receptor.
#' @param group Character vector. Cell group information used for computing averaged communication probabilities
#' @param min.percent Numeric from 0 to 1. Minimum percentage of expressed ligands or receptors 
#' per cell group to require for computing the group-level signaling. Default is 0.1.
#' @param min.cells.sr Integer greater than 0. Minimum number of cells required 
#' as senders or receivers per cell group for computing the group-level signaling. 
#' Default is 5.
#'
#' @return Returns a matrix containing the interaction weights between any two cell groups.
#' @export
computeAvgCommunProb_LR_SumX <- function (prob, group, dataLR = NULL,
                                          min.percent = 0.1, min.cells.sr = 5) {
  if (!is.factor(group)) group <- factor(group)
  gi <- as.integer(group); G <- nlevels(group); N <- nrow(prob)
  
  # cell.type.mat <- model.matrix(~group-1)                      
  Zs <- Matrix::sparseMatrix(i = seq_len(N), j = gi, x = 1, dims = c(N, G))
  
  d01 <- 1 * (dataLR > 0)
  m   <- cbind(tapply(d01[, 1], gi, mean), tapply(d01[, 2], gi, mean))   
  dataLR_percent <- 1 * (format(as.data.frame(m), digits = 1) >= min.percent)
  Prob_percent <- Matrix::crossprod(matrix(dataLR_percent[, 1], nrow = 1),
                                    matrix(dataLR_percent[, 2], nrow = 1))
  
  if (sum(Prob_percent) == 0) {
    Prob.avg <- Prob_percent
  } else {
    # Prob.avg <- Matrix::crossprod(x = cell.type.mat, y = prob %*% cell.type.mat)  
    Prob.avg <- Matrix::crossprod(Zs, prob) %*% Zs 
    pb <- prob; pb@x <- rep.int(1, times = length(pb@x))
    cs <- rowsum(cbind(Matrix::rowSums(pb), Matrix::colSums(pb)), group, reorder = TRUE) 
    cs <- 1 * (cs >= min.cells.sr)
    cells.sr <- Matrix::crossprod(matrix(cs[, 1], nrow = 1), matrix(cs[, 2], nrow = 1))
     
    Prob.avg <- Prob.avg * Prob_percent * cells.sr
  }
  dimnames(Prob.avg) <- list(levels(group), levels(group))
  as.matrix(Prob.avg)
}


#' Compute group-level cell-cell communication
#' 算法原理: 将每LRpair的每个细胞通讯概率(net$tmp$.prob.cell下的一个元素matrix), 
#' 将cell cluster 内所有细胞的通讯概率sum|mean, 每个LRpair生成一个以cell cluster为行列名的matrix
#' @param object SpatialCellChat object with communication probabilities for pairwise individual cells
#' @param group.by cell group information used for computing average communication probabilities
#' @param avg.type methods for integrating communication probabilities per cell group.
#' sum:总通讯量; avg: 每条连接的平均强度
#' @param type methods for computing the average gene expression per cell group.
#' By default = "triMean", defined as a weighted average of the distribution's median and 
#' its two quartiles (https://en.wikipedia.org/wiki/Trimean); 
#' When setting `type = "truncatedMean"`, a value should be assigned to 'trim'. See the function `base::mean`.
#' @param trim the fraction (0 to 0.25) of observations to be trimmed from each end of x before the mean is computed.
#' @param do.permutation whether performing permutation test
#' @param nboot the number of permutations
#' @param seed.use set a random seed. By default, set the seed to 1.
#' @param colocalization.use whether filtering out spatially distant cell groups 
#' based on colocalization analysis between any cell groups
#' @param thresh.colo removal of cell-cell communication with no significant colocalizations (fdr < 0.05)
#' @inheritParams computeAvgCommunProb_LR_Avg
#' @inheritParams computeAvgCommunProb_LR_Sum
#'
#' @return A CellChat object with updated slot 'net':
#' object@net$prob is the inferred group-level communication probability (strength) array, 
#' where the first, second and third dimensions represent 
#' a source group, target group and ligand-receptor pair, respectively.
#' object@net$pval is the corresponding p-values of each interaction
#' @export
#'
computeAvgCommunProbX <- function(object, group.by = NULL, avg.type = c("avg","sum"),
                                 min.percent = 0.1, min.cells.sr = 5, 
                                 do.permutation = T, nboot = 100,
                                 seed.use = 1L, colocalization.use = F, thresh.colo = 0.05) {
  if (is.null(group.by)) {
    group <- object@idents
  } else {
    if (!(group.by %in% colnames(object@meta))) {
      stop("The 'group.by' is not a column name in the `object@meta`, which will be used for cell grouping.")
    } else {group <- object@meta[[group.by]]}    
    if (!is.factor(group)) {group <- factor(group)}
  }
  
  cat(cli.symbol(),"The cell groups used for averaging cell-cell communication are ", cli::col_red(levels(group)), '\n')
  
  numCluster <- nlevels(group)
  if (numCluster != length(unique(group))) {
    stop("Please check `unique(object@idents)` and ensure that the factor levels are correct!
         You may need to drop unused levels using 'droplevels' function. e.g.,
         `meta$labels = droplevels(meta$labels, exclude = setdiff(levels(meta$labels),unique(meta$labels)))`")
  }
  
  #### if ( is.null(object@net$prob.cell) ) {
  if ( is.null(object@net$tmp) ) {
    stop(cli.symbol(2),"Please run `computeCommunProb` to compute 
         the communication probability/strength between any interacting individual cells! ")
  } else {
    #### prob.cell <- object@net$prob.cell
    prob.cell_ <- object@net$tmp$prob.cell # a list
  }
  
  nC <- nrow(prob.cell_[[1]])
  #### LRsig <- dimnames(prob.cell)[[3]]
  LRsig <- names(prob.cell_)
  nLR <- length(LRsig)
  
  interaction_input <- object@DB$interaction
  complex_input <- object@DB$complex
  cofactor_input <- object@DB$cofactor
  
  pairLRsig <- interaction_input[LRsig, , drop = FALSE]
  dataLavg <- object@net$tmp$Lavg
  dataRavg <- object@net$tmp$Ravg
  
  avg.type <- match.arg(avg.type)
  if(avg.type=="avg"){ computeAvgCommunProb_LR <- computeAvgCommunProb_LR_AvgX
  } else if (avg.type=="sum"){ computeAvgCommunProb_LR <- computeAvgCommunProb_LR_SumX }
  
  gc()
  
  if (colocalization.use) {
    data.spatial <- object@images$coordinates
    pval.colo = computeColocalization(coordinates = data.spatial, group = group, nboot = nboot, seed.use = seed.use)
  } else { pval.colo <- matrix(0, nrow = numCluster, ncol = numCluster) }
  
  cat(paste0(cli.symbol(),'Compute group-level cell-cell communication... <<< [', Sys.time(),']'),'\n')
  
  Prob <- array(0, dim = c(numCluster,numCluster,nLR))
  Pval <- array(1, dim = c(numCluster,numCluster,nLR))
  dimnames(Prob) <- list(levels(group), levels(group), rownames(pairLRsig))
  dimnames(Pval) <- dimnames(Prob)
  
  set.seed(seed.use)
  
  # retain dim-3, sum up dim-1 && dim-2, `prob.sum` stores each LR's number of cell-level links/interactions
  prob.sum <- purrr::map_dbl(.x = prob.cell_, .f = function(Mat){return(length(Mat@x))})
  names(prob.sum) <- LRsig
  object@net$tmp$LRsig.CCC.counts <- prob.sum
  LRsig.use.idx <- which(prob.sum > 0)
  object@net$tmp$LRsig.use.idx <- LRsig.use.idx
  gc()
  
  if(length(LRsig.use.idx) < 1){ stop("Each LR pair does not have any cell-level links/interactions.") }
  
  cat(cli.symbol(),"compute the average signaling per cell group...\n")
  
  Prob.avg_ <- my_future_sapply(
    X = seq_len(length(LRsig.use.idx)),
    FUN = function(x) {
      i <- LRsig.use.idx[[x]] # one LR pair index
      # compute the average signaling per cell group
      prob.cell.i <- prob.cell_[[i]]
      dataLR_temp <- cbind(dataLavg[i, ], dataRavg[i, ])
      Prob.avg <- computeAvgCommunProb_LR(prob.cell.i, group = group, dataLR = dataLR_temp, 
                                          min.percent = min.percent, min.cells.sr = min.cells.sr )
      Prob.avg[pval.colo > thresh.colo] <- 0
      gc()
      return(Prob.avg)
    },
    simplify = F # return a list
  )
  
  print('All pathways are done.')
  
  for (x in seq_len( length(LRsig.use.idx) ) ) {
    i <- LRsig.use.idx[[x]]
    Prob[ , , i] <- Prob.avg_[[x]]
  }
  
  # update `prob.sum` & `LRsig.use.idx` to do permutation
  # retain dim-3, sum up dim-1 && dim-2, `prob.sum` stores each LR's number of group-level links/interactions before permutation
  prob.sum <- apply(Prob > 0, 3, sum) # return a named vector
  # each LR's number of group-level links/interactions
  object@net$tmp$LRsig.GGC.counts <- prob.sum
  
  LRsig.use.idx <- which(prob.sum > 0)
  if (do.permutation) {
    cat(paste0(cli.symbol(),'Perform permutation test for group-level communication... <<< [', Sys.time(),']'),'\n')
    permutation <- replicate(nboot, sample.int(nC, size = nC))
    
    # compute the average signaling per cell group after permutation
    Pval_ <- my_future_lapply(
      X = seq_len(length(LRsig.use.idx)), # LRsig.use.idx is a numeric vector
      FUN = function(x){i <- LRsig.use.idx[[x]]
                        prob.cell.i <- prob.cell_[[i]]
                        dataLR_temp <- cbind(dataLavg[i,], dataRavg[i,])
                        Pnull <- as.vector(Prob[ , , i])
                        
                        Pboot <- sapply( X = 1:nboot, FUN = function(nE) {
                          groupboot <- group[permutation[, nE]]
                          Pboot.avg <- computeAvgCommunProb_LR(prob.cell.i, group = groupboot,
                                                               dataLR = dataLR_temp,
                                                               min.percent = min.percent, 
                                                               min.cells.sr = min.cells.sr )  
                          return(as.vector(Pboot.avg)) } )
                        
                        gc()
                        
                        Pboot <- matrix(unlist(Pboot), nrow=length(Pnull), ncol = nboot, byrow = FALSE)
                        nReject <- rowSums(Pboot - Pnull > 0)
                        p = nReject/nboot
                        Pval.i <- matrix(p, nrow = numCluster, ncol = numCluster, byrow = FALSE)
                        return(Pval.i) },
      
      simplify = F, # return a list
      hint.message = "do permutation..." )
    
    for (x in seq_len(length(LRsig.use.idx))) {
      i <- LRsig.use.idx[[x]]       # get correct index
      Pval[ , , i] <- Pval_[[x]]    # update the values
    }
    
    Pval[Prob == 0] <- 1    
  } else { Pval <- NULL }
  
  # Pval[Prob == 0] <- 1
  # dimnames(Prob) <- list(levels(group), levels(group), rownames(pairLRsig))
  # dimnames(Pval) <- dimnames(Prob)
  
  object@net$prob <- Prob
  object@net$pval <- Pval
  
  object@options$parameter$min.percent <- min.percent
  object@options$parameter$min.cells.sr <- min.cells.sr
  object@options$parameter$do.permutation <- do.permutation  
  object@options$parameter$nboot <- nboot
  object@options$parameter$avg.type <- avg.type
  object@options$parameter$seed.use <- seed.use
  object@options$parameter$colocalization.use <- colocalization.use
  
  object@options$parameter$thresh.colo <- thresh.colo
  # object@net$tmp$Lavg <- NULL;object@net$tmp$Ravg <- NULL; # clean the cache
  
  if (colocalization.use) { object@images$colocalization <- pval.colo }
  cat(paste0(cli.symbol(symbol = "success"),'Inference of group-level cell-cell communication is done. 
             Parameter values are stored in `object@options$parameter` <<< [', Sys.time(),']'))
  return(object)
}


#' Compute the communication probability on signaling pathway level by summarizing all related ligands/receptors
#'
#' @param object CellChat object
#' @param net A list from object@net; If net = NULL, net = object@net
#' @param pairLR.use A dataframe giving the ligand-receptor interactions; If pairLR.use = NULL, pairLR.use = object@LR$LRsig
#' @param thresh threshold of the p-value for determining significant interaction
#' @param do.group whether to compute the group-level signaling based on the cell group information in `object@idents`
#' @param do.cell whether to compute the individual-cell signaling at signaling pathway level. 
#' This works when "prob.cell" exists in `object@net`.
#' @param cell.names character vector. 细胞名，需与 `tmp$prob.cell` 各层的行列顺序一致。
#' 如果输入数据是object, 无需提供, 默认提取colnames(object@data.signaling)`；
#' 若输入数据为net, `object` 为 `NULL`，则需要提供。
#' 
#' @return A CellChat object with updated slot 'netP':
#' 含 `pathways`, `prob`, `pathways.cell`, `prob.cell`, `tmp`)；
#' 当 `object = NULL` 时返回 `netP` list。
#' @export
computeCommunProbPathwayX <- function(object = NULL, net = NULL, pairLR.use = NULL, 
                                      thresh = 0.05, do.group = TRUE, do.cell = TRUE, 
                                      cell.names = NULL) {
  if (is.null(net)) {
    if (is.null(object)) { stop(cli.symbol(2), "Please provide either `object` or `net`!") }
    net <- object@net
  }
  
  if (is.null(pairLR.use)) {
    if (is.null(object)) { stop(cli.symbol(2), "Please provide either `object` or `pairLR.use`!") }
    pairLR.use <- object@LR$LRsig
  }
  
  # ---------------- group level ----------------无改动
  if (do.group) {
    if (is.null(net$prob)) {stop("Please run `computeAvgCommunProb` to compute the group-level signaling!")}
    cat(cli.symbol(), "Compute the communication probability between cell groups 
        at signaling pathway level by summarizing all related ligands/receptors...\n")
    prob <- net$prob
    prob[net$pval >= thresh] <- 0
    pairLR.use <- pairLR.use[rownames(pairLR.use) %in% dimnames(prob)[[3]], , drop = FALSE]  ####
    pathways <- unique(pairLR.use$pathway_name)
    
    group <- factor(pairLR.use$pathway_name, levels = pathways) #####
    prob.pathways <- aperm(apply(prob, c(1, 2), by, group, sum), c(2, 3, 1))  ####
    
    pathways.sig <- pathways[apply(prob.pathways, 3, sum) != 0]
     # 保留组间通信非零的通路，并按总通信量降序排列
    prob.pathways.sig <- prob.pathways[, , pathways.sig, drop = FALSE]
    idx <- sort(apply(prob.pathways.sig, 3, sum), decreasing = TRUE, index.return = TRUE)$ix
    pathways.sig <- pathways.sig[idx]
    prob.pathways.sig <- prob.pathways.sig[, , idx, drop = FALSE]
  } else {
    pathways.sig <- NULL
    prob.pathways.sig <- NULL
  }
  
  netP <- list(pathways = pathways.sig, prob = prob.pathways.sig)
  
  # ---------------- individual-cell level ----------------
  if (do.cell) {
    prob.cell_ <- net$tmp$prob.cell   # a named list；替代原 net$prob.cell
    pairLR.use <- pairLR.use[rownames(pairLR.use) %in% names(prob.cell_), , drop = FALSE]
    
    if (!is.null(prob.cell_)) {
      
      pathways <- unique(pairLR.use$pathway_name)
      nC <- nrow(prob.cell_[[1]])
      
      # 细胞名：显式传入优先，其次 object，最后退化为索引
      if (is.null(cell.names) && !is.null(object)) {
        cell.names <- colnames(object@data.signaling)
      }
      
      if (is.null(cell.names)) {
        cell.names <- as.character(seq_len(nrow(prob.cell_[[1]])))
        message("Cell names are not available; using cell indices as names.")
      }
    
      cat(cli.symbol(), "Compute the communication probability between individual cells 
          at signaling pathway level by summarizing all related ligands/receptors...\n")
      
      gc()
      
      # 信号通路下, 所有受体/配体求和
      prob.all <- pbapply::pbsapply( X = pathways,
        FUN = function(one_pathway) {
          one_pathway_LRpair <- rownames(pairLR.use[pairLR.use$pathway_name == one_pathway,,drop=FALSE])
          prob.cell.i <- prob.cell_[one_pathway_LRpair]  # a list
          
          if (length(prob.cell.i) == 0L) {
            return(Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), dims = c(nC, nC))) 
          }
        
          Reduce(`+`, prob.cell.i)   # 同一 pathway 内各 LR 层求和
        },
        simplify = FALSE              # 返回 list
      )
      
      names(prob.all) <- pathways
      
      # 每个 pathway 的总通信概率：原 marginSumsSparse(MARGIN = 3)
      prob.sum <- vapply(prob.all, function(m) sum(m@x), numeric(1))
      names(prob.sum) <- pathways
      
      # 注意：filterCommunication 后可能出现全零层，prob.sum = 0 是合法的
      PathwaySig.use.idx <- which(prob.sum > 0)
      
      if (length(PathwaySig.use.idx) == 0L) {
        cat(cli.symbol(), "No pathway has non-zero individual cell-level communication.\n")
      } else {
        idx <- sort(prob.sum[PathwaySig.use.idx], decreasing = TRUE, index.return = TRUE)$ix
        PathwaySig.sort.idx <- PathwaySig.use.idx[idx]
        pathways.sig.cell <- pathways[PathwaySig.sort.idx]
        
        cat(cli.symbol(), "Subset the pathways with non-zero communication probability and 
            arrange them in a decreasing order based on the total communication probabilities ...\n")
        
        prob.cell.pathways.sig_ <- prob.all[PathwaySig.sort.idx]   # a list
        names(prob.cell.pathways.sig_) <- pathways.sig.cell
        prob.cell.pathways.sig <- my_as_sparse3Darray(prob.cell.pathways.sig_) # 3Darray
        
        dimnames(prob.cell.pathways.sig) <- list(cell.names, cell.names, pathways.sig.cell)
        cat(cli.symbol(), "The number of cells and pathways in Dim(prob.cell.pathways) are :",
            dim(prob.cell.pathways.sig), "\n")
        
        Tmp <- list(prob.cell = prob.cell.pathways.sig_)   # important, for parallel iteration
        netP$pathways.cell <- pathways.sig.cell
        netP$prob.cell <- prob.cell.pathways.sig
        netP$tmp <- Tmp
      }
    }
  }
  
  # group-level: pathways;prob
  # individual cell-level: pathways.cell;prob.cell
  if (is.null(object)) {
    cat(cli.symbol(1), "Computing the communication probability on signaling pathway level is done. \n")
    return(netP)
  } else {
    object@netP <- netP
    cat(cli.symbol(1), "Computing the communication probability on signaling pathway level is done. \n")
    return(object)
  }
}




#' Calculate the aggregated network by counting the number of links or summarizing the communication probability
#'
#' @param object CellChat object
#' @param sources.use,targets.use,signaling,pairLR.use Please check the description in function \code{\link{subsetCommunication}}
#' @param remove.isolate whether removing the isolate cell groups without any interactions when applying \code{\link{subsetCommunication}}
#' @param thresh threshold of the p-value for determining significant interaction
#' @param return.object whether return an updated CellChat object
#' @importFrom  dplyr group_by summarize groups
#' @importFrom stringr str_split
#'
#' @return Return an updated CellChat object:
#' `object@net$count` is a matrix: rows and columns are sources and targets respectively, 
#' and elements are the number of interactions between any two cell groups. 
#' `object@net$weight` is also a matrix containing the interaction weights between any two cell groups
#' `object@net$sum` is deprecated. Use `object@net$weight`
#'
#' @export
#'
aggregateNetX <- function(object, sources.use = NULL, targets.use = NULL, 
                         signaling = NULL, pairLR.use = NULL, remove.isolate = TRUE, 
                         thresh = 0.05, return.object = TRUE) {
  net <- object@net
  
  # 将LRpair的cell cluster层数据进行求和(aggregate)
  if (is.null(sources.use) & is.null(targets.use) & is.null(signaling) & is.null(pairLR.use)) {
    prob <- net$prob
    pval <- net$pval
    pval[prob == 0] <- 1
    prob[pval >= thresh] <- 0
    net$count <- apply(prob > 0, c(1,2), sum)
    net$weight <- apply(prob, c(1,2), sum)
    net$weight[is.na(net$weight)] <- 0
    net$count[is.na(net$count)] <- 0
    net$LR.sig <- dimnames(prob)[[3]][apply(prob, 3, sum) > 0]
  } else {
    df.net <- subsetCommunication(object, slot.name = "net", 
                                  sources.use = sources.use, targets.use = targets.use,
                                  signaling = signaling, pairLR.use = pairLR.use, 
                                  thresh = thresh)
    df.net$source_target <- paste(df.net$source, df.net$target, sep = "_")
    df.net2 <- df.net %>% group_by(source_target) %>% summarize(count = n(), .groups = 'drop')
    df.net3 <- df.net %>% group_by(source_target) %>% summarize(prob = sum(prob), .groups = 'drop')
    df.net2$prob <- df.net3$prob
    a <- stringr::str_split(df.net2$source_target, "_", simplify = T)
    df.net2$source <- as.character(a[, 1])
    df.net2$target <- as.character(a[, 2])
    cells.level <- levels(object@idents)
    
    if (remove.isolate) {
      message("Isolate cell groups without any interactions are removed. To block it, set `remove.isolate = FALSE`")
      df.net2$source <- factor(df.net2$source, levels = cells.level[cells.level %in% unique(df.net2$source)])
      df.net2$target <- factor(df.net2$target, levels = cells.level[cells.level %in% unique(df.net2$target)])
    } else {
      df.net2$source <- factor(df.net2$source, levels = cells.level)
      df.net2$target <- factor(df.net2$target, levels = cells.level)
    }
    
    count <- tapply(df.net2[["count"]], list(df.net2[["source"]], df.net2[["target"]]), sum)
    prob <- tapply(df.net2[["prob"]], list(df.net2[["source"]], df.net2[["target"]]), sum)
    net$count <- count
    net$weight <- prob
    net$weight[is.na(net$weight)] <- 0
    net$count[is.na(net$count)] <- 0
  }
  
  # 将LRpair的cell层数据进行求和(aggregate)
  if ( length(net$tmp$prob.cell) > 0L) {
    tl <- net$tmp$prob.cell                    
    N  <- nrow(tl[[1]])
    
    ii <- unlist(lapply(tl, function(m) m@i), use.names = FALSE) + 1L
    jj <- unlist(lapply(tl, function(m) rep.int(seq_len(ncol(m)), diff(m@p))), use.names = FALSE)
    xx <- unlist(lapply(tl, function(m) m@x), use.names = FALSE)
    
    net$weight.cell <- Matrix::sparseMatrix(i = ii, j = jj, x = xx, dims = c(N, N), index1 = TRUE)
    net$count.cell  <- Matrix::sparseMatrix(i = ii, j = jj, x = rep.int(1, length(ii)), dims = c(N, N), index1 = TRUE)
    net$LR.sig.cell <- names(tl)[vapply(tl, function(m) sum(m@x), numeric(1)) != 0]
    
    rm(ii, jj, xx); invisible(gc(FALSE))
    
    if (!is.null(sources.use) | !is.null(targets.use) | !is.null(signaling) | !is.null(pairLR.use)) {
      message("Subsetting cells or signaling is not applicable to individual cell-based `prob.cell`!", '\n')
    }
  }
  
  if (return.object) {
    object@net <- net
    return(object)
  } else { return(net) }
}




