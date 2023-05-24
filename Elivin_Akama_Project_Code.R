##########################################################################################
#Determining the Accuracy of Fast Fourier Transform Technique#
#compound Poisson distribution#
##########################################################################################
########Panjer's recursion algorithm using compound Poisson distribution#### ##############
#simulating individual claim probabilities based on the negative binomial distribution
#Initialize the parameters of the negative binomial distribution.
#Initialize parameter of poisson distribution
rm(list = ls())
n <- 1023
r =5
p = 0.5
lambda<-20

prob_dist<- dnbinom(1:n, r, p)/sum(dnbinom(1:n, r, p))

Panjer.Poisson <- function (prob_dist_X , lambda)
  # prob_dist_X = vector of probabilities for the claims, prob_dist_X[1]= P(X =1),we assume that P(X = 0) = 0 
  # If the sum of p is greater than 1 or any of the values in p are less than 0,
  # then it stops the execution of the function and returns an error message.
{ 
  if (sum(prob_dist_X )>1||any(prob_dist_X <0)) stop("prob_dist_X  parameter not a density")
  #If the product of lambda and the sum of p is greater than 727.
  #if it is, then it stops the execution of the function and returns an error message.
  if (lambda * sum(prob_dist_X ) > 727) stop("Underflow")
  #Initialize the probability generating function of poisson distribution evaluated at t=0. 
  cumul <- f <- exp(-lambda * sum(prob_dist_X ))
  #Initialize the iterations
  r <- length(prob_dist_X )
  s <- 0
  repeat
  { 
    s <- s+1
    m <- min(s, r)
    # compute the next value of the Panjer recursion
    #head(prob_dist_X, m) returns the first m and rev(tail(f, m)) returns the last m elements
    last <- lambda / s * sum(1:m * head(prob_dist_X ,m) * rev(tail(f,m)))
    #Appends the newly calculated value to the f vector and updates cumul
    f <- c(f,last)
    cumul <- cumul + last
    #It loops until the cumul is close to 1
    if (cumul > 0.99999999) break 
  }
  return(f) 
}
panj<-Panjer.Poisson(prob_dist, lambda)
panj
#Calculating the mean of the distribution obtained from the Panjer' recursion.
panj_mean<- sum(panj*(0:(length(panj)-1)))
panj_mean
#Calculating the theoretical mean of the total probabilities distribution.
theoretic_mean <- lambda*sum(prob_dist*(1:(length(prob_dist))))
theoretic_mean
#Plotting the distribution obtained from the Panjer' recursion.
plot(panj, type = "l", xlab = "k", ylab = "P(S=n)")

##########fast Fourier Transform using compound Poisson distribution#######################

#fft is an in-built function that computes total claim probabilities,
f <- Re(fft(exp(lambda*(fft(c(0,prob_dist))-1)), inverse=TRUE))/n
f
#Calculating the mean of the distribution obtained from the FFT
fft_mean<- sum(f*(0:(length(f)-1)))
fft_mean
#Initialize the length of the total probabilities obtained from the 2 distributions,
nn = min(length(panj), length(f))
#calculating the variation norm of probability measures
sum(abs(panj[1:nn]-f[1:nn]))
#Plotting the distributions of total claims obtained from the two techniques
plot(panj, type = "l", xlab = "k", ylab = "P(S=k)", col = "red")
lines(f, col = "blue")
legend("right", legend=c("Panjer's", "FFT"),
       col=c("red", "blue"), lty=1:2, cex=1.5)



##########################################################################################
#Determining the Accuracy of Fast Fourier Transform Technique#
#compound negative binomial distribution#
##########################################################################################
########Panjer's recursion algorithm using compound negative binomial distribution#### ###
#simulating individual claim probabilities based on the negative binomial Distribution
#Initialize the parameters of the negative binomial distribution and length of claim probabilities
rm(list = ls())
n <- 4000
r = 25
p = 0.5
prob_dist<- dnbinom(1:n, r, p)/sum(dnbinom(1:n, r, p))
panjer_neg_binom <- function(prob_dist_X, p, r) {
  #If the sum of p is greater than 1 or any of the values in
  #p are less than 0, then it stops the execution of the function and returns an error message.
  if (sum(prob_dist_X) > 1 || any(prob_dist_X < 0)) stop("p parameter not a probability mass function")
  #Initialize the length of the probability disribution
  z <- length(prob_dist_X)
  #Initialize the first term of the Panjer's recursion 
  cumul <- (p)^r
  f <- c((p)^r)
  #Initialize r,s to enter the loop for recursion
  s <- 0
  repeat {
    s <- s + 1
    m <- min(s, z)
    # compute the next value of the Panjer recursion
    panj <- (1-p)+(((1-p)*(r-1)*(1:m))/s)
    last <-sum(panj * head(prob_dist_X, m) * rev(tail(f, m)))
    #last is appended to the f vector and updated
    f <- c(f, last)
    cumul <- cumul + last
    #The function should run until it is too close to 1
    if (cumul > 0.99999999) break
  }
  return(f)
}
panj<-panjer_neg_binom(prob_dist, p, r)
sum(panj)
#calculating the mean of the distribution obtained from the Panjer' recursion
panj_mean<- sum(panj*(0:(length(panj)-1)))
panj_mean
#Calculating the theoretical mean of the total probabilities distribution
theoretic_mean <- (r*(1-p)/p)*sum(prob_dist*(1:(length(prob_dist))))
theoretic_mean
####FFT technique using compound negative binomial distribution#############
#computing total claim probabilities using fft which is an in-built function
f <- Re(fft((p/(1 - (1 - p) * fft(c(0,prob_dist))))^r, inverse = TRUE))/n
f
#Calculating the mean of the distribution obtained from the FFT
fft_mean<- sum(f*(0:(length(f)-1)))
fft_mean
#calculating the variation norm of probability measures
nn = min(length(panj), length(f))
sum(abs(panj[1:nn]-f[1:nn]))
#Plotting the distributions of total claims obtained from the two techniques
plot(panj, type = "l", xlab = "k", ylab = "P(S=k)", col = 2)
lines(f, col = 1)
legend("right", legend=c("Panjer's", "FFT"),
       col=c(2, 1), lty=1:2, cex=1.5)




##########################################################################################
#Determining the Accuracy of Normal Approximation Technique#
#compound Poisson distribution#
##########################################################################################
#simulating individual claim probabilities and initialize the parameters
rm(list = ls())
n <- 1028
r =5
p = 0.1
lambda<-100
prob_dist<- dnbinom(1:n, r, p)/sum(dnbinom(1:n, r, p))
Panjer.Poisson <- function (prob_dist_X , lambda)
 
{ 
  if (sum(prob_dist_X )>1||any(prob_dist_X <0)) stop("prob_dist_X  parameter not a density")
  if (lambda * sum(prob_dist_X ) > 727) stop("Underflow")
  #Initialize the first term of the Panjer's recursion 
  cumul <- f <- exp(-lambda * sum(prob_dist_X ))
  # initialize the iterations
  r <- length(prob_dist_X )
  s <- 0
  repeat
  { 
    s <- s+1
    m <- min(s, r)
    # compute the next value of the Panjer recursion
    last <- lambda / s * sum(1:m * head(prob_dist_X ,m) * rev(tail(f,m)))
    #last is appended to the f vector and updated
    f <- c(f,last)
    cumul <- cumul + last
    #The function should run until it is too close to 1
    if (cumul > 0.99999999) break 
  }
  return(f) 
}

ff<-Panjer.Poisson(prob_dist, lambda)
ff
nn <- length(ff)
#Calculating the theoretical mean of the total probabilities distribution
theoretic_mean <- lambda*sum(prob_dist*(1:(length(prob_dist))))
theoretic_mean
plot(ff, type = "l")
#calculates the first and second moments of the prob_dist_X distribution.
mean_X = sum(prob_dist*1:n)
mean_X2 = sum(prob_dist*1:n*1:n)
#calculates the mean of the cumulative distribution.(poisson,mean=variance)
mean_cumul_dist = lambda*mean_X 
#calculates the standard deviation of the cumulative distribution.
sd_cumul_dist = sqrt(lambda*mean_X2)
#Plotting the distributions of total claims obtained from Panjer's recursion and normal approximation
plot(ff, type = "l",xlab = "k", ylab = "P(S=k)", col="black")
lines(dnorm(0:(nn-1), mean_cumul_dist , sd_cumul_dist), col = "red")
#initializes an empty vector
normal_approx = c()
#starts a loop that iterates from 1 to the length of the ff vector.
for(m in 1:nn){
  #calculates the corresponding Z-scores and stores the calculated probabilities
  prob_m = pnorm((m-0.5-mean_cumul_dist)/sd_cumul_dist) - pnorm((m-1.5-mean_cumul_dist)/sd_cumul_dist)
  normal_approx = c(normal_approx, prob_m)
}
#calculating the variation norm of probability measures
sum(abs(normal_approx - ff))
legend("left", legend=c("Panjer's", "normal"),
       col=c("red", "black"), lty=1:2, cex=1.2)




##########################################################################################
#Determining the Accuracy of Shifted Gamma Approximation Technique#
#compound Poisson distribution#
##########################################################################################
#simulating individual claim probabilities based on negative binomial distribution and initialize the parameters
rm(list = ls())
n <- 1028
r = 5
p = 0.1
lambda<-70
prob_dist<- dnbinom(1:n, r, p)/sum(dnbinom(1:n, r, p))
Panjer.Poisson <- function (prob_dist_X , lambda)
{ 
  if (sum(prob_dist_X )>1||any(prob_dist_X <0)) stop("prob_dist_X  parameter not a density")
  if (lambda * sum(prob_dist_X ) > 727) stop("Underflow")
  #Initialize the first term of the Panjer's recursion
  cumul <- f <- exp(-lambda * sum(prob_dist_X ))
  #Initialize the iterations
  r <- length(prob_dist_X )
  s <- 0
  repeat
  { 
    s <- s+1
    m <- min(s, r)
    # compute the next value of the Panjer recursion
    last <- lambda / s * sum(1:m * head(prob_dist_X ,m) * rev(tail(f,m)))
    #last is appended to the f vector and updated
    f <- c(f,last)
    cumul <- cumul + last
    #The function should run until it is too close to 1
    if (cumul > 0.99999999) break 
  }
  return(f) 
}

ff<-Panjer.Poisson(prob_dist, lambda)
ff
nn <- length(ff)
plot(ff, type = "l")
#calculates the first, second and third moments of the prob_dist_X distribution.
mean_X = sum(prob_dist*1:n)
mean_X 
mean_X2 = sum(prob_dist*1:n*1:n)
mean_X2 
mean_X3 = sum(prob_dist*1:n*1:n*1:n)
mean_X3
#calculates the parameters of a shifted gamma distribution
alpha = 4*lambda*(((mean_X2)^3)/((mean_X3)^2))
beta = 2*(mean_X2/mean_X3)
x0 = (lambda * mean_X) - (2*lambda*(((mean_X2)^2)/mean_X3))
#Plotting the distributions of total claims obtained from Panjer's recursion and shifted gamma approximation
plot(ff, type = "l",xlab = "k", ylab = "P(S=k)", col=1)
lines(0:(nn-1),dgamma(0:(nn-1)-x0,shape=alpha,scale= 1/beta), col = 2)
#initializes an empty vector 
sgamma_approx = c()
#starts a loop that iterates from 1 to the length of the ff vector.
for(m in 1:nn){
  #calculates the corresponding Z-scores and stores the calculated probabilities
  prob_m =pgamma(m+0.5-x0, shape=alpha,scale= 1/beta)-pgamma(m-0.5-x0, shape=alpha,scale= 1/beta)
  sgamma_approx = c(sgamma_approx, prob_m)
}
lines(sgamma_approx, col = "green")

print(sgamma_approx)
#calculating the variation norm of probability measures
sum(abs(sgamma_approx - ff))
legend("left", legend=c("Panjer's", "shifted gamma"),
       col=c(1, "green"), lty=1:2, cex=0.9)