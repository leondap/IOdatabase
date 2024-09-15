
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

