# Read in data sheet
load("truncated_data.Rda")
# tr[[i]] now refers to the data sheet for the i'th restaurant

# Pickup Model
pickup <- function(data,h){
  # Log+1 transformation (remove info beyond h in this step)
  log_data <- apply(data+1,2,log)[,1:(h+1)]
  # Split into training data and forecasting values
  top <- log_data[1:(nrow(log_data)-h),]
  # Calculate grossing up factors
  # Take the mean of each column and then take the difference of those means
  means <- apply(top,2,mean)
  devs <- rev(diff(rev(means)))
  # Values that are grossed up from when forecasting
  end <- tail(log_data,n=h)
  diag.end <- end[upper.tri(end,diag=FALSE)]
  # Forecast
  # Extract diagonal values to forecast from
  forecast.v <- diag(end[,2:(h+1)])
  # Lower triangle matrix of gross up factors
  dev.mat<-t(matrix(rep(devs,(h)),nrow=h))
  dev.mat <- dev.mat[lower.tri(dev.mat,diag=TRUE)]
  # Add column of zeros for grossing up purposes
  dev.mat.2 <- cbind(dev.mat,rep(0,h))
  # Add diagonal values to development factor matrix
  # Extract diagonal elements as matrix
  # Include 0 column
  beg.elem <- cbind(rep(0,h),diag(diag(end)))
  # Add
  tot <- dev.mat.2 + beg.elem
  # Cumulate gross up factors and forecasting values
  fin <- t(apply(tot,1,function(x){rev(cumsum(rev(x)))}))
  # Extract forecasted column
  fc <- fin[,1]
  # Exponentiate
  fc <- exp(fc)-1
  # Return forecast
  return(fc)
}