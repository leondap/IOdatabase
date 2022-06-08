
palette.col<-function(a,b,c,d,size=10){
res<-NULL
rmatrix<-matrix(NA,size,size)
gmatrix<-rmatrix
bmatrix<-rmatrix

argb<- as.vector(col2rgb(a))
brgb<- as.vector(col2rgb(b))
crgb<- as.vector(col2rgb(c))
drgb<- as.vector(col2rgb(d))
rvect<-c(argb[1],brgb[1],crgb[1],drgb[1])
stepr1<-(rvect[2]-rvect[1])/(size-1)
stepr2<-(rvect[4]-rvect[3])/(size-1)


gvect<-c(argb[2],brgb[2],crgb[2],drgb[2])
stepg1<-(gvect[2]-gvect[1])/(size-1)

stepg2<-(gvect[4]-gvect[3])/(size-1)


bvect<-c(argb[3],brgb[3],crgb[3],drgb[3])
stepb1<-(bvect[2]-bvect[1])/(size-1)

stepb2<-(bvect[4]-bvect[3])/(size-1)



for(c in 1:size){
for(r in 1:size){
#c<-1
#r<-1
x1<-rvect[1]+((c-1)*stepr1)
x2<-rvect[3]+((c-1)*stepr2)
y1<-(r-1)/(size-1)
rmatrix[r,c]<-round((x1*(1-y1))+(x2*y1),0)

x1<-gvect[1]+((c-1)*stepg1)
x2<-gvect[3]+((c-1)*stepg2)
gmatrix[r,c]<-round((x1*(1-y1))+(x2*y1),0)

x1<-bvect[1]+((c-1)*stepb1)
x2<-bvect[3]+((c-1)*stepb2)
bmatrix[r,c]<-round((x1*(1-y1))+(x2*y1))

}
}
res$red<-rmatrix
res$green<-gmatrix
res$blue<-bmatrix

return(res)
}






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
