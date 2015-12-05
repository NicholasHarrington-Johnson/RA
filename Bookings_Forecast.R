# Read in data sheet
load("truncated_data.Rda")
# tr[[i]] now refers to the data sheet for the i'th restaurant

# Pickup Method
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
  diag.end <- end*upper.tri(end,diag=FALSE)
  # Forecast
  # Extract diagonal values to forecast from
  forecast.v <- diag(end[,2:(h+1)])
  # Lower triangle matrix of gross up factors
  dev.mat<-t(matrix(rep(devs,(h)),nrow=h))
  dev.mat <- dev.mat*lower.tri(dev.mat,diag=TRUE)
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
# Regression models with ARIMA Errors
booking.arima <- function(data,h,phols=0,full.fit=0,knots=NA){
  # Take training set
  new.data <- data[(1:(nrow(data)-h)),]
  # Apply time series properties
  new.data$pubd <- ts(new.data$pubd,end=tsp(data$pubd)[1]+((nrow(data)-h)/365))
  new.data$pubd <- ts(new.data$pubd,start=tsp(data$pubd)[1],frequency=tsp(data$pubd)[3])
  # Perform log transformation
  log.attend <- apply(new.data+1,2,log)[,1]
  # Fix frequency to 7
  log.attend <- ts(log.attend,start=1,frequency=7)
  # Remove zeros from data
  ###
  # Note that this step appears to create bad forecasts
  # BOTH in spline model AND in regular ARIMA type models
  log.attend[log.attend<1] <- NA
  # Create x regressors as required
  # Take Public Holiday Variables if required
  if(phols==1){
    xdums <- new.data[,(ncol(new.data)-2):ncol(new.data)]
    colnames(xdums) <- c("pubd","pubi","pubny")
    # Figure out if there are public holidays in sample
    total.xdum <- colSums(xdums)
    # Remove dummy variables if there are no public holidays in sample
    if (total.xdum[1]==0) {xdums$pubd <- NULL}
    if (total.xdum[2]==0) {xdums$pubi <- NULL}
    if (total.xdum[3]==0) {xdums$pubny <- NULL} 
    # Fit Regression Model with Spline Terms as required
    if (is.na(knots)){
      # Fit Regression Model with Public Holidays and ARIMA Errors
      fit <- auto.arima(log.attend,xreg=xdums)
    } else if (knots==0){
      # Attach unaltered bookings column to forecasts
      spline.terms <- new.data[,(h+1)]
    } else if (knots>0){
      # Knot values set at quantiles
      # Zeros are removed when selecting knots
      knot.pos <- quantile(data[,(h+1)][data[,(h+1)]>0],prob=((1:knots)/(knots+1)))
      # Whole b_t,j is used to fit spline
      full.spline <- ns(data[,(h+1)],knots=knot.pos)
      # Tail of full.spline is saved for forecasting
      spline.terms <- full.spline[1:nrow(new.data),]
    }
    # Fit model to spline terms
    if(!is.null(xdums)&!is.null(spline.terms)){
      spline.xdums <- cbind(xdums,spline.terms)
      fit <- auto.arima(log.attend,xreg=spline.xdums)
    }   
  } else if (phols==0) {
    # Fit Regression Model with Spline Terms as required
    if (is.na(knots)){
      # Fit ARIMA Model
      fit <- auto.arima(log.attend)
    } else if (knots==0){
      # Attach unaltered bookings column to forecasts
      spline.terms <- new.data[,(h+1)]
    } else if (knots>0){
      # Knot values set at quantiles
      # Zeros are removed when selecting knots
      knot.pos <- quantile(data[,(h+1)][data[,(h+1)]>0],prob=((1:knots)/(knots+1)))
      # Whole b_t,j is used to fit spline
      full.spline <- ns(data[,(h+1)],knots=knot.pos)
      # Tail of full.spline is saved for forecasting
      spline.terms <- full.spline[1:nrow(new.data),]
    }
    # Fit model to spline terms
    if(!is.na(knots)){
      fit <- auto.arima(log.attend,xreg=spline.terms)
    }
  }
  # Create forecasting x regressors as required
  if(phols==1){
    xfor <- data[((nrow(data)-h+1):nrow(data)),(ncol(data)-2):ncol(data)]
    colnames(xfor) <- c("pubd","pubi","pubny")
    # Remove x regressors as required
    if (total.xdum[1]==0) {xfor$pubd <- NULL} 
    if (total.xdum[2]==0) {xfor$pubi <- NULL}
    if (total.xdum[3]==0) {xfor$pubny <- NULL}
    # Forecast with Public Holiday Dummies and Spline Terms as required
    if (is.na(knots)){
      fc <- forecast(fit,h,xreg=xfor) 
    }
    else if (knots==0){
      spline.for <- tail(data[,(h+1)],n=h)
    } else if (knots>0){
      spline.for <- tail(full.spline,h)
    }
    if (!is.null(xfor)&!is.null(spline.for)){
      xfor <- cbind(xfor,spline.for)
      fc <- forecast(fit,h,xreg=xfor)
    }
  } else if (phols==0){
    # Forecast without public holidays
    if (is.na(knots)){
      fc <- forecast(fit,h) 
    }
    else if (knots==0){
      spline.for <- tail(data[,(h+1)],n=h)
    } else if (knots>0){
      spline.for <- tail(full.spline,h)
    }
    if (!is.na(knots)){
      xfor <- spline.for
      fc <- forecast(fit,h,xreg=xfor)
    }
  }
  # Exponentiate
  fc$mean <- exp(fc$mean)-1
  fc$lower <- exp(fc$lower)-1
  fc$upper <- exp(fc$upper)-1
  fc$x <- ts(new.data[,1],frequency=365)
  tsp(fc$x)<-tsp(new.data$pubd)
  fc$mean <- ts(fc$mean, end =tsp(data$pubd)[2], frequency=365)
  fc$mean <- ts(fc$mean, start = (tsp(fc$x)[2]+(1/365)), frequency=365)
  tsp(fc$upper) <- tsp(fc$lower) <- tsp(fc$mean) 
  if(full.fit==1){
    obj <- list(fit,fc)
    return(obj)
  }
  return(fc)
}
# Function that can iterate through past bookings
h_step <- function(h,...){
  fc <- booking.arima(h,...)
  obj<-fc$mean[h]
  return(obj)
}
# Function that iterates through past bookings
full.spline.regression <- function(data,h=14,phols=0,knots=0){
  # Define vector to forecast over
  days <- c(1:h)
  # Iterate through vector of forecast horizons
  fc <- sapply(days,function(horizon) {h_step(h=horizon,data=data,phols=phols,knots=knots)})  
  # Apply time series properties to forecast
  fc <- ts(fc,start=tsp(data$pubd)[2]-(h/365),end=tsp(data$pubd)[2],frequency=365)
  return(fc)
}
