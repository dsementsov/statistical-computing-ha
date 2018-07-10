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
x <- c(1.9690,-4.0650,-12.7918,1.7495,-0.5349,-0.5384,1.2888,-1.7306,29.7583,-18.9475,
      -17.9632,2.7976,-2.8081,1.9147,9.7825,-1.4867,-0.2228,2.5757,17.1286,0.7200)


#The log-likelihood function
logL = function(theta){-length(x)*log(pi) - sum(log(1+(x-theta)^2))}
#The first derivative of the logLike function
derL = function(theta){sum(2*(x-theta)/(1+(x-theta)^2))}
#The second derivative of the logLike function
sderL = function(theta){sum(-2*(1-(x-theta)^2)/(1+(x-theta)^2)^2)}

# 1. (a)#Bisection method

#first consider boundaries where we are going to estimate the parameters of the LL function

epsilon = 1e-04 #the proximity constant

optim.BM = function(x,epsilon){

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
while(dist.BM>=epsilon){
  if(derL(a_0)*derL(c)<0){b_0 = c
  }else{a_0 = c}
  c = a_0 + (b_0-a_0)/2 #starting midpoint  
  dist.BM = abs(a_0 - c)/(abs(a_0)+epsilon)
  counter.BM = counter.BM +1
  if(counter.BM>=100){
    warning("Number of iterations exceeds the limit of 100")
    #break the loop if the number of iterations exceeds 100
    break}
  }
  #final midpoint
final.BM = a_0 + (b_0-a_0)/2 #final midpoint
value.BM = logL(final.BM)    #maximized logLikelihood
return(list(opt = final.BM, value = value.BM, count = counter.BM))}
#


# 1. (b) #Newton`s method


guess = 0
optim.NM = function(x,epsilon,guess){
#guess a candidate solution
theta = guess #guess the location parameter
counter.NM = 0 #reset the counter
#use modified relative convergence criterion
dist.NM = abs(derL(theta)/sderL(theta))/(abs(theta)+epsilon) 

while(dist.NM>=epsilon){
  theta = theta - derL(theta)/sderL(theta)#the updating equation
  dist.NM = abs(derL(theta)/sderL(theta))
  counter.NM=counter.NM+1
  if(counter.NM>=100){
    warning("Number of iterations exceeds the limit of 100")
    #break the loop if the number of iterations exceeds 100
    break}
    }
final.NM = theta #calculated parameter value
value.NM = logL(final.NM) #maximized logLikelihood
return(list(opt = final.NM, value = value.NM, count = counter.NM))}


# 1. (c)#scaled fixed point iteration

optim.FPI = function(x,epsilon, guess){
psi = guess #guess the location parameter
counter.FPI = 0 #reset the counter
alpha = 1/8
#use modified relative convergence criterion
dist.FPI = abs(alpha*derL(psi))/(abs(psi)+epsilon)

while(dist.FPI>=epsilon){
  psi = psi + alpha*derL(psi) #the updating equation
  dist.FPI = abs(alpha*derL(psi))/(abs(psi)+epsilon)
  counter.FPI=counter.FPI+1
  if(counter.FPI>=100){
    warning("Number of iterations exceeds the limit of 100")
    #break the loop if the number of iterations exceeds 100
    break}
  }
final.FPI = psi #calculated parameter value
value.FPI = logL(final.FPI) #maximized logLikelihood

return(list(opt = final.FPI, value = value.FPI, count = counter.FPI))}


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

# 3. #here we combine all the results into the table
tab = matrix(0,3,4)
tab[1,] = c(BM.max$opt,BM.max$value,BM.max$count,BM.max$time)
tab[2,] = c(NM.max$opt,NM.max$value,NM.max$count,NM.max$time)
tab[3,] = c(FPI.max$opt,FPI.max$value,FPI.max$count,FPI.max$time)

tab = as.data.frame(tab)
row.names(tab) = c("Bisection","Newton","Fixed Point")
colnames(tab)  = c("Optimum","Max.LL.Value","Iterations", "Time")
print(tab)


########################################################################

# Part II

# 1. (a)


# 1. (b)


# 2.


########################################################################

# Part III

# 1.



# 2.



# 3.



# 4.


