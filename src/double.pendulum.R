library(e1071)
library(nnet)
source("utils.R")


#################################################################################################
# Use Euler's method to resolve the movement equations (see the formula in the lab3.pdf).
#################################################################################################
next_thetas <- function(theta1=0, theta2=0, h=0, p.theta1=0, p.theta2=0) {

  mass <- 0.5
  L <<- 0.5
  g <- 9.81
  #first order derivate of theta1 (angle of the first pendulum)
  theta1.n <- 0
  #first order derivate of theta2 (angle of the second pendulum)
  theta2.n <- 0
  #first pendulum moment
  p.theta1.n <- 0
  #second pendulum moment
  p.theta2.n <- 0

  # Commons
  ml2 <- mass * (L**2)
  denom <- 16 - 9 * cos(theta1-theta2)**2

  # Theta1 i+1
  f1 <- ( 6 / ml2 ) * (2*p.theta1 - 3*cos(theta1-theta2) * p.theta2) / denom
  theta1.n <- theta1 + h * f1

  # Theta2 i+1
  f2 <- ( 6 / ml2 ) * (8*p.theta2 - 3*cos(theta1-theta2) * p.theta1) / denom
  theta2.n <- theta2 + h * f2

  # Moment1 i+1
  p.theta1.n <- p.theta1 - h * ml2/2 * ( f1 * f2 * sin(theta1-theta2) + 3 * (g/L) * sin(theta1) )

  # Moment2 i+1
  p.theta2.n <- p.theta2 - h * ml2/2 * ( -f1 * f2 * sin(theta1-theta2) + (g/L) * sin(theta2) )

  result <- list( theta=c(theta1.n, theta2.n), moment=c(p.theta1.n,p.theta2.n) )
  result
}
####################################################################################################
# Use this function to create the dataset for training
####################################################################################################

create_dataset <- function(theta1=0, theta2=0, h=0, p.theta1=0, p.theta2=0)
{
  T1 <- vector()
  T2 <- vector()

  X <- array(dim=c(1000,4))
  Y <- array(dim=c(1000,2))

  for (i in 1:1002) {

    result <- next_thetas(theta1=theta1, theta2=theta2, h=h, p.theta1=p.theta1, p.theta2=p.theta2)

    theta1=result$theta[1]
    theta2=result$theta[2]
    p.theta1=result$moment[1]
    p.theta2=result$moment[2]

    T1 <- c(T1,theta1)
    T2 <- c(T2,theta2)

    if(i>2) {
      #######
      
      X[i-2,] <- c(T1[i-2],T1[i-1],T2[i-2],T2[i-1])
      Y[i-2,] <- c(T1[i],T2[i])
      
      ########
    }
  }

  dataset <- list(X=X,Y=Y,T1=T1,T2=T2)

  dataset
}
################################################################
#
##################################################################
prediction.pendulum.by.ann <- function(ann.model, theta1=0, theta2=0, h=0, p.theta1=0, p.theta2=0)
{
  iter <- 1000
  T1 <- NULL
  T2 <- NULL
  T1.pred <- c(0,0)
  T2.pred <- c(0,0)

  for (i in 1:iter) {

    result <- next_thetas(theta1=theta1, theta2=theta2, h=h, p.theta1=p.theta1, p.theta2=p.theta2)

    theta1=result$theta[1]
    theta2=result$theta[2]
    p.theta1=result$moment[1]
    p.theta2=result$moment[2]

    T1 <- c(T1,theta1)
    T2 <- c(T2,theta2)

    if(i>2) {
      
      #######
      #
      # ADD YOUR CODE HERE
      #
      ########
    }
    # plot.thetas(1:i,T1,T2,TN1,TN2,abs(T1-TN1),abs(T2-TN2), dynamic=TRUE)
  }

  ret <- list(T1=T1,T2=T2,TN1=T1.pred,TN2=T2.pred)

  ret
}

prediction.pendulum.by.svr <- function(svr1.model, svr2.model ,theta1=0, theta2=0, h=0, p.theta1=0, p.theta2=0)
{
  iter <- 1000
  T1 <- NULL
  T2 <- NULL
  T1.pred <- c(0,0)
  T2.pred <- c(0,0)

  for (i in 1:iter) {

    result <- next_thetas(theta1=theta1, theta2=theta2, h=h, p.theta1=p.theta1, p.theta2=p.theta2)

    theta1=result$theta[1]
    theta2=result$theta[2]
    p.theta1=result$moment[1]
    p.theta2=result$moment[2]

    T1 <- c(T1,theta1)
    T2 <- c(T2,theta2)

    if(i>2) {
      #######
      #
      # ADD YOUR CODE HERE
      #
      ########
      
      }
    # plot.thetas(1:i,T1,T2,TN1,TN2,abs(T1-TN1),abs(T2-TN2), dynamic=TRUE)
    
  }

  ret <- list(T1=T1,T2=T2,TN1=T1.pred,TN2=T2.pred)

  ret
}

run_pendulum_experiment <- function(){

  ## Initial Values
  theta1=pi/2
  theta2=pi/2
  h=0.01
  p.theta1=0
  p.theta2=0

  #values of theta1 predicted by ANN
  theta1.ann <- c(0,0)
  #values of theta2 predicted by ANN
  theta2.ann <- c(0,0)

  #values of theta1 predicted by SVR
  theta1.svr <- c(0,0)
  #values of theta2 predicted by SVR
  theta2.svr <- c(0,0)

  ## Create the initial dataset
  ## Please complete the function "create_dataset"
  mydataset <- create_dataset(theta1=theta1, theta2=theta2, h=h, p.theta1=p.theta1, p.theta2=p.theta2)
  T1 <- mydataset$T1 #"real" theta1
  T2 <- mydataset$T2 #"real" theta2
  X <- mydataset$X
  Y <- mydataset$Y

  ## spliting the datasets for training & testing phases

  split_perc <- 0.6

  splits <- split.data(mydataset$X, split_perc)
  X.train <- splits$train
  X.test <- splits$test

  splits <- split.data(mydataset$Y, split_perc)
  Y.train <- splits$train
  Y.test <- splits$test


  ## ANN
  m <- nnet(x=X.train, y=Y.train, size=4, maxit=500, linout=TRUE,abstol = 1.0e-4, reltol = 1.0e-4)
 
  # predict for ANN. Call the method predict.ann

  prediction <- prediction.pendulum.by.ann(ann.model=m,theta1=theta1, theta2=theta2, h=h, p.theta1=p.theta1, p.theta2=p.theta2)

  TN1 <- prediction$TN1 # predicted theta1
  TN2 <- prediction$TN2 #"predicted" theta2

  ## To plot the results we need the variables T1,T2,TN1,TN2

  plot.double.pendulum.init(L, pen.width=0.01,name="Pendulum ANN")
  t_init <- 500
  plot.double.pendulum(T1[t_init], T2[t_init], L, col="blue")
  plot.double.pendulum(TN1[t_init], TN2[t_init], L, col="red")
  for(i in (t_init+1):1000){

    plot.double.pendulum(T1[i-1], T2[i-1], L, col="white")
    init.polygon( H=1.5, W=2, e=0.05, pen.width=0.01)
    plot.double.pendulum(T1[i], T2[i], L, col="blue")

    plot.double.pendulum(TN1[i-1], TN2[i-1], L, col="white")
    init.polygon( H=1.5, W=2, e=0.05, pen.width=0.01)
    plot.double.pendulum(TN1[i], TN2[i], L, col="red")

    Sys.sleep(0.1)
    cat(paste(i, "\n"))
  }

  ## Uncomment if you want to plot thetas
  #plot.thetas(500:1000,T1,T2,TN1,TN2,abs(T1-TN1),abs(T2-TN2))

  ## SVR

  m1 <- svm(x=X.train, y=Y.train[,1])
  m2 <- svm(x=X.train, y=Y.train[,2])

  # predict for SVR

  prediction <- prediction.pendulum.by.svr(svr1.model = m1, svr2.model = m2, theta1=theta1, theta2=theta2, h=h, p.theta1=p.theta1, p.theta2=p.theta2)

  TN1 <- prediction$TN1 # predicted theta1
  TN2 <- prediction$TN2 #"predicted" theta2

  ## To plot the results we need the variables T1,T2,TN1,TN2

  plot.double.pendulum.init(L, pen.width=0.01,name="Pendulum SVR")
  t_init <- 500
  plot.double.pendulum(T1[t_init], T2[t_init], L, col="blue") # "real" values
  plot.double.pendulum(TN1[t_init], TN2[t_init], L, col="red") # predicted values
  
  for(i in (t_init+1):1000){

    plot.double.pendulum(T1[i-1], T2[i-1], L, col="white")
    init.polygon( H=1.5, W=2, e=0.05, pen.width=0.01)
    plot.double.pendulum(T1[i], T2[i], L, col="blue")

    plot.double.pendulum(TN1[i-1], TN2[i-1], L, col="white")
    init.polygon( H=1.5, W=2, e=0.05, pen.width=0.01)
    plot.double.pendulum(TN1[i], TN2[i], L, col="red")

    Sys.sleep(0.1)
    cat(paste(i, "\n"))
  }

  ## Uncomment if you want to plot thetas
  #plot.thetas(500:1000,T1,T2,TN1,TN2,abs(T1-TN1),abs(T2-TN2))
}
