########################################################################
#             Home Assignment in Statistical Computing                 #
########################################################################
#                                                                      #
# List of the group members:                                           #
# 1. Kindop, Igor, 1020566, Stu210866@mail.uni-kiel.de                 #
# 2. Okuneva, Mariia, 1015904, mokuneva@stat-econ.uni-kiel.de          #    
# 3. Sementsov, Denys, 1001084                                         #
# 4. Zaspa, Uliana, 1021140, uliana.zaspa@hotmail.com                  #
#                                                                      #
########################################################################

#### Part I Univariate Optimization ####

# Data:
x = c(1.9690,-4.0650,-12.7918,1.7495,-0.5349,-0.5384,1.2888,-1.7306,29.7583,-18.9475,
      -17.9632,2.7976,-2.8081,1.9147,9.7825,-1.4867,-0.2228,2.5757,17.1286,0.7200)

# Setup
set.seed(90)
require(LaplacesDemon)
library(LaplacesDemon)

#The log-likelihood function
logL = function(theta){-length(x)*log(pi) - sum(log(1+(x-theta)^2))}

#The first derivative of the logLike function
derL = function(theta){sum(2*(x-theta)/(1+(x-theta)^2))}

#The second derivative of the logLike function
sderL = function(theta){sum(-2*(1-(x-theta)^2)/(1+(x-theta)^2)^2)}

# Finding maximum
maximum = optimize(f = logL, interval = c(-1,1), maximum = TRUE)


# 1. (a)#Bisection method

#first consider boundaries where we are going to estimate the parameters of the LL function

epsilon = 1e-04 #the proximity constant

optim.BM = function(x,epsilon)
{

    #pick the left and rigt boundaries
    a_0 = min(x) #left
    b_0 = max(x) #right
    
    #We assume the solution is located between these points
    #According with the provided data, assume the scale is unity
    
    #Starting midpoint
    c = (a_0 + b_0)/2
    
    #use absolute convergence criterion
    dist.BM = abs(a_0 - c)/(abs(a_0)+epsilon)
    counter.BM = 0 #set the counter at 0
    
    #Use a stopping criteria as a distance between the admissible search interval
    while(dist.BM >= epsilon){
        
        if(derL(a_0)*derL(c)<0)
        {
            b_0 = c
        }
        else
        {
            a_0 = c
        }
        
        
        c = a_0 + (b_0-a_0)/2 #starting midpoint  
        dist.BM = abs(a_0 - c)/(abs(a_0)+epsilon)
        counter.BM = counter.BM +1
        
        if(counter.BM >= 100) 
        { 
            warning("Number of iterations exceeds the limit of 100")
            #break the loop if the number of iterations exceeds 100
            break
        }
    }
    
    #final midpoint
    final.BM = a_0 + (b_0-a_0)/2 #final midpoint
    value.BM = logL(final.BM)    #maximized logLikelihood
    
    return(list(opt = final.BM, value = value.BM, count = counter.BM))
}


# 1. (b) #Newton`s method


guess = 0
optim.NM = function(x,epsilon,guess)
{
    #guess a candidate solution
    theta = guess #guess the location parameter
    counter.NM = 0 #reset the counter
    
    #use modified relative convergence criterion
    dist.NM = abs( derL(theta) / sderL(theta) ) / ( abs(theta) + epsilon ) 
    
    while(dist.NM >= epsilon)
    {
        theta = theta - derL(theta)/sderL(theta) #the updating equation
        dist.NM = abs( derL(theta) / sderL(theta) )
        counter.NM= counter.NM + 1
        if(counter.NM>=100)
        {
            warning("Number of iterations exceeds the limit of 100")
            #break the loop if the number of iterations exceeds 100
            break
        }
    }
    final.NM = theta #calculated parameter value
    value.NM = logL(final.NM) #maximized logLikelihood
    return(list(opt = final.NM, value = value.NM, count = counter.NM))
}


# 1. (c) # Scaled fixed point iteration

optim.FPI = function(x,epsilon, guess)
{
    psi = guess #guess the location parameter
    counter.FPI = 0 #reset the counter
    alpha = 1/8
    #use modified relative convergence criterion
    dist.FPI = abs(alpha*derL(psi))/(abs(psi)+epsilon)
    
    while(dist.FPI>=epsilon)
    {
        
        psi = psi + alpha*derL(psi) #the updating equation
        dist.FPI = abs(alpha*derL(psi))/(abs(psi)+epsilon)
        counter.FPI=counter.FPI+1
        
        if(counter.FPI>=100)
        {
            warning("Number of iterations exceeds the limit of 100")
            #break the loop if the number of iterations exceeds 100
            break
        }
    
    }
    final.FPI = psi #calculated parameter value
    value.FPI = logL(final.FPI) #maximized logLikelihood
    
    return(list(opt = final.FPI, value = value.FPI, count = counter.FPI))
}


# 2. #Here we use the proc.time() function

#Bisection method
ptm = proc.time()
BM.max = optim.BM(x,epsilon) 
time.BM = proc.time()-ptm 
BM.max$time = time.BM[3]

#Newton`s method
ptm = proc.time()
NM.max = optim.NM(x,epsilon,guess)
time.NM = proc.time()-ptm
NM.max$time = time.NM[3]

#Scaled. Fixed point method
ptm = proc.time()
FPI.max = optim.FPI(x,epsilon,guess)
time.FPI = proc.time()-ptm
FPI.max$time = time.FPI[3]


# 3.  #here we combine all the results into the table

tab = matrix(0,3,4)
tab[1,] = c(BM.max$opt,BM.max$value,BM.max$count,BM.max$time)
tab[2,] = c(NM.max$opt,NM.max$value,NM.max$count,NM.max$time)
tab[3,] = c(FPI.max$opt,FPI.max$value,FPI.max$count,FPI.max$time)

tab = as.data.frame(tab)
row.names(tab) = c("Bisection","Newton","Fixed Point")
colnames(tab)  = c("Optimum","Max.LL.Value","Iterations", "Time")
print(tab)
maximum$maximum

########################################################################

#### Part II Monte-Carlo ####

# 1. (a) # Computing a standard Monte Carlo estimator


num_sims = 1000

mean_vector = c()

for (i in 1:num_sims) 
{
    # Drawing 100 000 observations from standart normal
    x = rnorm(100000)
    
    h = x^2/(exp(x)-1)
    # mean of a sample
    mean_vector[i] = mean(h)
}

# Monte Carlo estimate
mean_mc = mean(mean_vector)
mean_mc


# 1. (b) Computing an antihetic estimator (based on the original sample split)

num_sims = 1000

h_1_mean = c()
h_2_mean = c()

for (i in 1:num_sims) 
{
    # Sample of original size / 2
    x_1 = rnorm(100000/2)
    
    # Negative transformation
    x_2 = -x_1
    
    h_1 = x_1^2/(exp(x_1)-1)
    h_2 = x_2^2/(exp(x_2)-1)
    
    h_1_mean[i] = mean(h_1)
    h_2_mean[i] = mean(h_2)
}

h1_mc = mean(h_1_mean)
h2_mc = mean(h_2_mean)

# Antihectic estimator
mean_antithetic = (h1_mc+h2_mc)/2
mean_antithetic


# 2. Finding sigma^2 using importance sampling 

par(mfrow = c(1,1))
# X has the density proportional to q(x)
q = function(x) { exp( -( abs(x)^3 ) / 3) }

plot(q, xlim = c(-5,5), ylim = c(0,2))
# shape of a normal distribution with higher maximum -> propose normal distribution

num_sims = 1000
isampling_mean = c()

# running importance sampling mc
for (i in 1:num_sims)
{
    # Sampling from x
    sample = rnorm(100000)
    
    # Draw from q
    f = q(sample)
    
    # Draw for the proposed function
    g = dnorm(sample)
    
    # Weights
    weights =  f / g
    # Standardizing weights
    weights = weights / sum(weights)
    
    # Using Standardized lets us forget about propotionality constant for q
    # Estimating mean of x^2
    isampling_mean[i] = mean( sample^2 * weights )
}

isampling_mean = na.omit(isampling_mean)  # Ommitting all division by 0
sigma_squared = mean(isampling_mean)
format(sigma_squared, scientific = FALSE)

########################################################################

##### Part III Metropolis - Hastings algorythm ####

# Data setup

# True parameters:
m = c(7,10)    #mean
s = c(0.5,0.5) #standard deviation
p = c(0.7,0.3) #probabilities

Nrep = 10000     #number of iterations

true = rnormm( Nrep, p = p, mu = m, sigma = s) #true density

# create a matrix of x, where the MH values will be stored
# each column represents a different starting value
x = cbind(rep(0, Nrep) ,rep(7, Nrep), rep(15, Nrep)) 


# 1. MH algorithm
for (j in 1:3) # for each of the starting values
{      
    for (i in 2:Nrep)  # for each iteration
    {
        x[i,j] = rnorm(1,mean=x[i-1,j],sd=0.01)  # generate a value from Normal distribution with a mean of a previous value from x matrix
        
        target_new = dnormm(x[i,j],p=p,mu=m,sigma=s)   # generate a value from a target distribution, conditional on a newly drawn x from a previous step
        target_old = dnormm(x[i-1,j],p=p,mu=m,sigma=s) # generate a value from a target distribution conditional on a previously drawn x
        
        alpha = min(1,target_new/target_old)  # calculate acceptance ratio
        u=runif(1)  # generate a value from uniform distribution
        
        if(u <= alpha)
        {      # if the value of the acceptance ratio is less than the value from a uniform distribution, update the value of x
            x[i,j] = x[i,j]
        } 
        else 
        {      # else, keep the old value
            x[i,j] = x[i-1,j]
        }
    }
}

# 2. Plotting

y=density(true, from=0, to=17)

# parameters of a plot
layout(matrix(c(1, 2), 1, 2), widths=c(3, 1))
par(mar=c(5, 1, 2, 1), oma=c(0, 4.1, 0, 0))

# sample path plot for each output
plot(x[,1],
     ylim=c(0,17),
     type="s",
     xlab="iterations",
     ylab="x values",
     main="sample path")
lines(x[,2],col="red", type="s")
lines(x[,3],col="blue", type="s")
legend("topright",
       c("x_0=0", "x_0=7", "x_0=15"),
       lty=c(1,1,1),
       col= c("black", "red", "blue"),
       cex=0.7)
plot(y$y, y$x,
     ylim=c(0,17),
     main="",
     ylab="",
     xlab="",
     type="l",
     axes=FALSE)

# 3. Histograms

par(mfrow=c(3,1))
hist(x[,1], probability = T, 25,
     xlim=c(0,17),
     ylim=c(0,1.2),
     main="Realizations, starting value x=0")
lines(density(true), col="red")

hist(x[,2], probability = T, 25,
     xlim=c(0,17),
     ylim=c(0,1.2),
     main="Realizations, starting value x=7")
lines(density(true), col="red")

hist(x[,3], probability = T, 25,
     xlim=c(0,17),
     ylim=c(0,1.2),
     main="Realization, starting value x=15")
lines(density(true), col="red")

# 4. Suggestion

# As one can see, convergence to the true density does not occur (at least for 10000 iterations)
# The problem is that the standard error of the proposal distribution is too small, leading to the very small steps towards the true density
# Hence, to improve the convergence we suggest using N(x_0, 0.2^2) as a proposal distribution.
# We could also set a true standard error and use N(x_0, 0.5^2) as a proposal, but assuming we do not know the true density, let us take the N(x_0, 0.2^2)
for (j in 1:3)
{
    for (i in 2:Nrep)
    {
        x[i,j] = rnorm(1,mean=x[i-1,j],sd=0.2)
        
        target_new = dnormm(x[i,j],p=p,mu=m,sigma=s)
        target_old = dnormm(x[i-1,j],p=p,mu=m,sigma=s)
        
        alpha = min(1,target_new/target_old)
        u = runif(1)
        
        if(u <= alpha)
        {
            x[i,j] = x[i,j]
        } 
        else 
        {
            x[i,j] = x[i-1,j]
        }
    }
}


y = density(true, from=0, to=17)
layout(matrix(c(1, 2), 1, 2), widths=c(3, 1))
par(mar=c(5, 1, 2, 1), oma=c(0, 4.1, 0, 0))
plot(x[,1],
     ylim=c(0,17),
     type="s",
     xlab="iterations",
     ylab="x values",
     main="sample path")
lines(x[,2],col="red", type="s")
lines(x[,3],col="blue", type="s")
legend("topright", c("x_0=0", "x_0=7", "x_0=15"), lty=c(1,1,1), col= c("black", "red", "blue"), cex=0.7)
plot(y$y, y$x,
     ylim=c(0,17),
     main="",
     ylab="",
     xlab="",
     type="l",
     axes=FALSE)


par(mfrow=c(3,1))
hist(x[,1], probability = T, 25,
     xlim=c(0,17),
     ylim=c(0,1.2),
     main="Realizations, starting value x=0")
lines(density(true), col="red")

hist(x[,2], probability = T, 25,
     xlim=c(0,17),
     ylim=c(0,1.2),
     main="Realizations, starting value x=7")
lines(density(true), col="red")

hist(x[,3], probability = T, 25,
     xlim=c(0,17),
     ylim=c(0,1.2),
     main="Realizations, starting value x=15")
lines(density(true), col="red")
