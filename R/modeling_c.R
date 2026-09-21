#' @description 用my_as_sparse3Darray将list形式的 转换成3D稀疏array, 并存入net槽的prob.cell
#' @export
net3Darray <- function(object, use.raw = FALSE ) {
  Prob.cell <- my_as_sparse3Darray(object@net$tmp$prob.cell)  
  if (use.raw) { cell <- colnames(object@data.raw) } else { cell <- colnames(object@data)}  
  dimnames(Prob.cell) <- list(cell, cell, names(object@net$tmp$prob.cell) ) 
  object@net$prob.cell <- Prob.cell
  return(object)
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
computeCommunProbX <- function (
    object,
    LR.use = NULL,
    raw.use = TRUE,
    Kh = 0.5,
    n = 1,
    distance.use = TRUE,
    tol = NULL, # will be removed in the future
    interaction.range = 250,
    scale.distance = 0.01,
    use.AGAN = T,
    contact.dependent = TRUE,
    contact.range = 10,
    contact.dependent.forced = FALSE) {
  
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
        
        result@x[abs(result@x) < 0.00001] <- 0
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
  Tmp <- list(prob.cell = Prob.cell_, Lavg = dataLavg, Ravg = dataRavg)   
  net <- list(prob.cell = NULL, tmp = Tmp)  
  # 释放内存垃圾
  rm(Prob.cell_)
  gc()
  
  execution.time = Sys.time() - ptm
  object@options$run.time <- as.numeric(execution.time,
                                        units = "secs")
  object@images[["result.computeCellDistance"]] <- res
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
filterProbabilityX <- function (
    object,
    nboot = 100,
    seed.use = 666L,
    thresh = 0.05
){
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
      
      if (i <= nLR1) {
        d_spatial <- d.spatial
      } else {
        d_spatial <- adj.contact
      }
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
    prob.cell <- my_as_sparse3Darray(prob.cell_)
    dimnames(prob.cell) <- list(cell.names, cell.names, pair.LR.use)
    names(prob.cell_) <- pair.LR.use
    object@net$prob.cell <- prob.cell
    object@net$tmp$prob.cell <- prob.cell_
    
    cat(cli.symbol(1), "Filtering is done.\n")
    return(object)
  } # whether to filter out
}

