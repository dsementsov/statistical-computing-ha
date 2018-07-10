########################################################################
#             Home Assignment in Statistical Computing                 #
########################################################################
#
# List of the group members:
# 1. Last name, name, ID, e-mail
# 2. Last name, name, ID, e-mail
# 3. Last name, name, ID, e-mail
# 4. Last name, name, ID, e-mail
# 5. Last name, name, ID, e-mail
#
########################################################################

# Part I

# Data:
x = c(1.9690,-4.0650,-12.7918,1.7495,-0.5349,-0.5384,1.2888,-1.7306,29.7583,-18.9475,
      -17.9632,2.7976,-2.8081,1.9147,9.7825,-1.4867,-0.2228,2.5757,17.1286,0.7200)
n=20

# 1. (a)


# 1. (b)


# 1. (c)


# 2.


# 3.


########################################################################

# Part II

# 1. (a)

set.seed(90)
num_sims = 1000

mean_vector = c()

for (i in 1:num_sims) {
  x = rnorm(100000)
  h = x^2/(exp(x)-1)
  mean_vector[i] = mean(h)
}

mean_mc = mean(mean_vector)
mean_mc

# 1. (b)
set.seed(90)
num_sims = 1000

h_1_mean = c()
h_2_mean = c()
for (i in 1:num_sims) {
  x_1 = rnorm(100000/2)
  x_2 = -x_1
  h_1 = x_1^2/(exp(x_1)-1)
  h_2 = x_2^2/(exp(x_2)-1)
  h_1_mean[i] = mean(h_1)
  h_2_mean[i] = mean(h_2)
}
h1_mc = mean(h_1_mean)
h2_mc = mean(h_2_mean)

mean_antithetic = (h1_mc+h2_mc)/2
mean_antithetic

# 2.


########################################################################

# Part III

# 1.
library(LaplacesDemon)

# True parameters:
m = c(7,10)    #mean
s = c(0.5,0.5) #standard deviation
p = c(0.7,0.3) #probabilities

Nrep=10000     #number of iterations

true=rnormm(Nrep,p=p,mu=m,sigma=s) #true density

# create a matrix of x, where the MH values will be stored
# each column represents a different starting value
x = cbind(rep(0,Nrep),rep(7,Nrep),rep(15,Nrep)) 


# MH algorithm
for (j in 1:3){      # for each of the starting values
  for (i in 2:Nrep){ # for each iteration
    x[i,j] = rnorm(1,mean=x[i-1,j],sd=0.01)  # generate a value from Normal distribution with a mean of a previous value from x matrix
    target_new = dnormm(x[i,j],p=p,mu=m,sigma=s)  # generate a value from a target distribution, conditional on a newly drawn x from a previous step
    target_old = dnormm(x[i-1,j],p=p,mu=m,sigma=s)# generate a value from a target distribution conditional on a previously drawn x
    alpha = min(1,target_new/target_old)  # calculate acceptance ratio
    u=runif(1)  # generate a value from uniform distribution
    if(u<=alpha){ # if the value of the acceptance ratio is less than the value from a uniform distribution, update the value of x
      x[i,j] = x[i,j]
    } else {      # else, keep the old value
      x[i,j] = x[i-1,j]
    }
  }
}

# 2.

y=density(true, from=0, to=17)

# parameters of a plot
layout(matrix(c(1, 2), 1, 2), widths=c(3, 1))
par(mar=c(5, 1, 2, 1), oma=c(0, 4.1, 0, 0))

# sample path plot for each output
plot(x[,1], ylim=c(0,17), type="s", xlab="iterations", ylab="x values", main="sample path")
lines(x[,2],col="red", type="s")
lines(x[,3],col="blue", type="s")
legend("topright", c("x_0=0", "x_0=7", "x_0=15"), lty=c(1,1,1), col= c("black", "red", "blue"), cex=0.7)
plot(y$y, y$x, ylim=c(0,17), main="", ylab="",xlab="",type="l", axes=FALSE)

# 3.

par(mfrow=c(3,1))
hist(x[,1], probability = T, 25, xlim=c(0,17), ylim=c(0,1.2), main="Realizations, starting value x=0")
lines(density(true), col="red")

hist(x[,2], probability = T, 25, xlim=c(0,17), ylim=c(0,1.2), main="Realizations, starting value x=7")
lines(density(true), col="red")

hist(x[,3], probability = T, 25, xlim=c(0,17), ylim=c(0,1.2), main="Realization, starting value x=15")
lines(density(true), col="red")

# 4.

# As one can see, convergence to the true density does not occur (at least for 10000 iterations)
# The problem is that the standard error of the proposal distribution is too small, leading to the very small steps towards the true density
# Hence, to improve the convergence we suggest using N(x_0, 0.2^2) as a proposal distribution.
# We could also set a true standard error and use N(x_0, 0.5^2) as a proposal, but assuming we do not know the true density, let us take the N(x_0, 0.2^2)
for (j in 1:3){
  for (i in 2:Nrep){
    x[i,j] = rnorm(1,mean=x[i-1,j],sd=0.2)
    target_new = dnormm(x[i,j],p=p,mu=m,sigma=s)
    target_old = dnormm(x[i-1,j],p=p,mu=m,sigma=s)
    alpha = min(1,target_new/target_old)
    u=runif(1)
    if(u<=alpha){
      x[i,j] = x[i,j]
    } else {
      x[i,j] = x[i-1,j]
    }
  }
}


y=density(true, from=0, to=17)
layout(matrix(c(1, 2), 1, 2), widths=c(3, 1))
par(mar=c(5, 1, 2, 1), oma=c(0, 4.1, 0, 0))
plot(x[,1], ylim=c(0,17), type="s", xlab="iterations", ylab="x values", main="sample path")
lines(x[,2],col="red", type="s")
lines(x[,3],col="blue", type="s")
legend("topright", c("x_0=0", "x_0=7", "x_0=15"), lty=c(1,1,1), col= c("black", "red", "blue"), cex=0.7)
plot(y$y, y$x, ylim=c(0,17), main="", ylab="",xlab="",type="l", axes=FALSE)


par(mfrow=c(3,1))
hist(x[,1], probability = T, 25, xlim=c(0,17), ylim=c(0,1.2), main="Realizations, starting value x=0")
lines(density(true), col="red")

hist(x[,2], probability = T, 25, xlim=c(0,17), ylim=c(0,1.2), main="Realizations, starting value x=7")
lines(density(true), col="red")

hist(x[,3], probability = T, 25, xlim=c(0,17), ylim=c(0,1.2), main="Realizations, starting value x=15")
lines(density(true), col="red")

