get_groups <- function(paired_t_test_result, alpha = 0.05, rm.subset = FALSE) {
  
  # Get a matrix of the alpha values
  temp <- paired_t_test_result$p.value
  
  # Make a square matrix to populate with the alpha values.
  n <- nrow(temp)
  mat.names <- c(colnames(temp), rownames(temp)[n])
  my.mat <- matrix(data = NA, nrow = n+1, ncol = n+1)
  colnames(my.mat) <- mat.names
  rownames(my.mat) <- mat.names

    # Add diagonal.
  for (i in 1:nrow(my.mat)) {
    my.mat[ i, i] <- 0
  }
  
  # Get vector of p.values
  stat <- na.exclude(as.vector(paired_t_test_result$p.value))

  # Add other cells to square matrix.
  k=1
  for (j in 1:(nrow(my.mat)-1)) {
    for (i in ((j+1):nrow(my.mat))) {
      my.mat[i,j] <-  my.mat[j,i] <- stat[k]
      k=k+1
    }
  }

  # For each column, get list of treatments not significantly different.
  grp <- list()
  trts <- colnames(my.mat)
  for (i in 1:ncol(my.mat)) {
    grp[[i]] <-c(trts[i], names(which(my.mat[ , i] > alpha)))
  }
  
  # Remove groups that are sub-sets of other groups
  k <- 0
  del <- vector()
  for (i in 1:(length(grp)-1)) {
    for ( j in (i+1):length(grp)) {
      if (!rm.subset) {
        if (setequal(grp[[i]], grp[[j]])) {
          k <- k+1
          del[k] <- j
        }
      }
      else {
        if (all(is.element(grp[[i]], grp[[j]]))) {
          k <- k+1
          del[k] <- i
        }
        else if (all(is.element(grp[[j]], grp[[i]]))) {
          k <- k+1
          del[k] <- j
        }
      }

    }
  }

  del <- unique(del)
  del <- del[order(del, decreasing = TRUE)]

  for (i in 1:length(del)) {
    grp[[del[i]]] <- NULL
  }
  
  return(list(groups = grp, p.matrix = my.mat))
  
}