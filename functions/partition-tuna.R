partition_tuna <- function(i, dat, test_set) {
  if (i == 'all') {
    training <-  dat
    
    test <- dat
    
  } else {
    training <- slice(dat, -test_set[i])
    
    test <- slice(dat, test_set[i])
    
  }
  
  out = list(training = training, test = test)
  
  return(out)
}
