function percentile_84, array
  n= 0L
  n = n_elements(array)
  sarray = array(sort(array)) ; sort
  n75 = 0L
  n75 = floor(0.84*n) 
  out = sarray(n75)
  
  return, out
end
