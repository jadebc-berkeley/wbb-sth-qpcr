
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

# make a 2x2 table with the count and percent
# in parentheses
# inputs: 
# tab = 2x2 table of counts
# ptab = 2x2 table of percentages
make2x2 = function(tab, ptab, label){
  a = tab[1,1]
  b = tab[1,2]
  c = tab[2,1]
  d = tab[2,2]
  
  ptab = ptab*100
  ptabf = sprintf("%0.0f", ptab)
  ap = ptabf[1]
  bp = ptabf[3]
  cp = ptabf[2]
  dp = ptabf[4]
  
  out = rbind(
    cbind(
      label, "qPCR +", "qPCR -"
    ),
    cbind(
      "Kato-Katz +",
      paste0(a, " (", ap, "%)"),
      paste0(b, " (", bp, "%)")
    ),
    cbind(
      "Kato-Katz -",
      paste0(c, " (", cp, "%)"),
      paste0(d, " (", dp, "%)")
    )
  )
  
  return(as.data.frame(out))
  
}