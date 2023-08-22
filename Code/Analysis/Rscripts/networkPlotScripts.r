
trophiclevels = function(net){
  #Get adjacency matrix 
  mat <- get.adjacency(net, sparse=F)
  
  #Detect basal node
  basal <- rownames(subset(mat, apply(mat, 2, sum)==0) & apply(mat, 1, sum)!= 0)
  #Compute the shortest path to basal node
  paths_prey <- suppressWarnings(shortest.paths(graph = net, v= V(net), to = V(net)[basal], 
                                                mode = "in", weights = NULL, algorithm = "unweighted"))
  
  paths_prey[is.infinite(paths_prey)] <- NA
  shortest_paths <- suppressWarnings(as.matrix(apply(paths_prey, 1, min, na.rm=TRUE)))
  #for species with no prey apart of them
  shortest_paths[is.infinite(shortest_paths)] <- NA
  
  # Shortest TL
  sTL <- 1 + shortest_paths
  
  # Compute average prey trophic level
  # inspired from cheddar package calculation
  W <- t(mat)
  rs <- rowSums(W)
  W <- W/matrix(rs, ncol = ncol(W), nrow = nrow(W))
  W[0 == rs, ] <- 0
  I <- diag(ncol(W))
  tl0<-rowSums(I - W)
  result <- tryCatch(solve(I - W), error = function(e) e)
  if ("error" %in% class(result)) {
    avtl <- rep(NA, ncol(pm))
    names(avtl) <- colnames(pm)
  }
  else {
    avtl <- rowSums(result)
  }
  
  # Short-weighted TL is the average between 
  # Shortest TL and average prey TL
  SWTL <- (sTL + avtl)/2
  
  return(SWTL)
}

plotfw = function(net, col=NULL, lab=NULL, size=NULL,
                   nylevel=7, maxsize=10, labcex=0.01,
                   ynum=6, ylab= "Trophic Level", ...){
  n <- vcount(net)
  if (!is.null(col)){
    V(net)$color <- col
  } else{
    V(net)$color <- rainbow(vcount(net))
  }
  if (!is.null(lab)){
    V(net)$name <- lab
  }
  if (!is.null(size)){
    V(net)$size <- size
  } else {
    V(net)$size <- maxsize/2
  }
  
  tl <- trophiclevels(net)
  dgpred <- tl

  bks <- c(0.9, seq(1.9, max(tl), length.out = nylevel))
  ynod <- cut(tl, breaks = bks, include.lowest = TRUE, 
              labels = 1:(length(bks)-1))

  maxx <- max(table(ynod))
  xnod <- rep(0,n)
  for (i in 1:nylevel){
    l <- sum(ynod==i)
    
    ltr <- (l/maxx)**(1/2)*maxx
    if (l>1) {
      xnod[ynod==i] <- seq(-ltr,ltr,length.out = l)
    } else {
      xnod[ynod==i] <- 0
    }
  }
  
  coo <- cbind(xnod,tl)
  
  #y axis with 1 and continuous axis from 2 to max TL.
  yax <- c(1, seq(2, max(tl), length.out = ynum-1))
  labax <- round(yax, 1)
  #rescale xax between -1 and 1
  laby <- (yax-min(yax))/(max(yax)-min(yax))*2 - 1
  
  plot(net, layout= coo, vertex.label.color="black", 
       vertex.label.cex=labcex, ...)
  axis(2, at = laby, labels = labax)
  mtext(ylab, side = 2, line=2)
  res <- data.frame(coo, "size"= V(net)$size, "color"= V(net)$color)
  names(res) <- c("x", "y", "size", "color")
  row.names(res) <- V(net)$name
  invisible(res)
}