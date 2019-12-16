function percentile_16, array
  n= 0L
  n = n_elements(array)
  sarray = array(sort(array)) ; sort
  n25 = 0L
  n25 = floor(0.16*n) 
  out = sarray(n25)
  return, out
end
