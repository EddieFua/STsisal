clean_count = function(st_count){
  st_count= log(st_count+1)
  idx = apply(st_count, 1, function(x)
    any(x != 0))
  st_count = st_count[idx, ]
  #delete genes detected in less than 5%of pixels
  idy = apply(st_count, 1, function(x)
    (sum(x != 0) / ncol(st_count) > 0.05))
  st_count = st_count[idy, ]
  return(st_count)
}