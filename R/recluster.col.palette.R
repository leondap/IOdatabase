
recluster.col.palette<-function (mat, palette,st = T) {
  mat2 <- mat
  mat3 <- mat2
  mat3[, 2] <- (mat2[, 2] - min(mat2[, 2]))/(max(mat2[, 
                                                      2]) - min(mat2[, 2]))
  mat3[, 1] <- (mat2[, 1] - min(mat2[, 1]))/(max(mat2[, 
                                                      1]) - min(mat2[, 1]))
  colour <- array(data = 0, dim = c(dim(mat3)[1], (dim(mat3)[2]) + 
                                      3))
  division=nrow(palette$red)
  for (t in 1:dim(mat3)[1]) {
    colour[t, 1] <- mat3[t, 1]
    colour[t, 2] <- mat3[t, 2]
    location<-c(mat3[t, 1]*division,mat3[t, 2]*division)
    colour[t, 3] <- palette$red[round((mat3[t, 1]*(division-1))+1,0),round((mat3[t, 2]*(division-1))+1,0)]
    colour[t, 4] <- palette$green[round((mat3[t, 1]*(division-1))+1,0),round((mat3[t, 2]*(division-1))+1,0)]
    
    
    colour[t, 5] <- palette$blue[round((mat3[t, 1]*(division-1))+1,0),round((mat3[t, 2]*(division-1))+1,0)]
    
    
    rownames(colour) <- rownames(mat)
    if (!st) {
      colour[, 1] <- mat2[, 1]
      colour[, 2] <- mat2[, 2]
    }
  }
  return(colour)
}
