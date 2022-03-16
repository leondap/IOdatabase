

recluster.col3D<-function (mat, st = FALSE, rot = FALSE) {

rx<-range(mat[,1])
ry<-range(mat[,2])
rz<-range(mat[,3])
x<-abs(range(mat[,1])[1]-range(mat[,1])[2])
y<-abs(range(mat[,2])[1]-range(mat[,2])[2])
z<-abs(range(mat[,3])[1]-range(mat[,3])[2])    
max<-max(x,y,z)

red<-((mat[,1]-rx[1])/(rx[2]-rx[1]))*255
blue<-((mat[,2]-ry[1])/(ry[2]-ry[1]))*255
green<-((mat[,3]-rz[1])/(rz[2]-rz[1]))*255


colour <- matrix(NA,nrow(mat),5)
colour[,1]<-mat[,1]
colour[,2]<-mat[,2]
colour[,3]<-floor(red)
colour[,4]<-floor(green)
colour[,5]<-floor(blue)

return(colour)
}
