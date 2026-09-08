

FixRepeatSymbols <- function(data, repeatVals) {
  
  # Loop over 
  i <- 1
  while (i < length(data)) {
    
    # Define current level of variable
    level <- data[[i]]
    
    # Define next value
    j <- i + 1
    nextLevel <- data[[j]]
    
    # If the next level is symbol representing the repeated value
    if (nextLevel %in% repeatVals & i <= length(data)) {
      
      # Loop and replace any next value with previous
      while (nextLevel %in% repeatVals & j <= length(data)) {
        
        # Set this value to previous level
        data[[j]] <- level
        
        # Update j
        j <- j + 1
        
        # Check next level
        nextLevel <- data[[j]]
        
      }
      
      # Update i to howver many j is past previous i
      i <- i + (j-i)
      
      # Else
    } else {
      
      # Else iterete by 1
      i <- i + 1
      
    }
    
  }
  
  # Return
  return(data)
  
}
