
ptest.format=function(x,decimals,scale){
  return(sprintf(paste0("%0.",decimals,"f"),x*scale))
}


ptestci.format=function(x,lb,ub,decimals,scale){
  x=sprintf(paste0("%0.",decimals,"f"),x*scale)
  lb=sprintf(paste0("%0.",decimals,"f"),lb*scale)
  ub=sprintf(paste0("%0.",decimals,"f"),ub*scale)
  return(paste0(x," (",lb,", ",ub, ")"))
}

# replace NA with --
repNA=function(x){
  x=as.character(x)
  x[is.na(x)]="--"
  return(x)
}

