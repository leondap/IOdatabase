
recluster.col2<-function (mat, st = T) {
    mat2 <- mat
    mat3 <- mat2
  	mat3[, 2] <- (mat2[, 2] - min(mat2[, 2]))/(max(mat2[, 
            2]) - min(mat2[, 2]))
       mat3[, 1] <- (mat2[, 1] - min(mat2[, 1]))/(max(mat2[, 
            1]) - min(mat2[, 1]))
       colour <- array(data = 0, dim = c(dim(mat3)[1], (dim(mat3)[2]) + 
        3))
    for (t in 1:dim(mat3)[1]) {
        colour[t, 1] <- mat3[t, 1]
        colour[t, 2] <- mat3[t, 2]
        colour[t, 5] <- round((1 - max(mat3[t, 2], mat3[t, 1])) * 
            255)
        colour[t, 4] <- round(mat3[t, 1] * 255)
        colour[t, 3] <- round(mat3[t, 2] * 255)
        rownames(colour) <- rownames(mat)
        if (!st) {
            colour[, 1] <- mat2[, 1]
            colour[, 2] <- mat2[, 2]
        }
    }
    return(colour)
}
