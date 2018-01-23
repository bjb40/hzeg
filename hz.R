#linear hazard stuff
#simulates linear hazard and transformations
#Bryce Bartlett

alpha = .01
beta = .002

#failure rate = - S(t)dt/S(t)
hazard = function(alpha,beta,t){
  return(alpha + beta*t)
}

#Survivors function = exp(-int_0^t h(u)du)
survival = function(alpha,beta,t){
  return(exp(-1*(alpha*t + (beta*t^2)/2)))
}

#Cumulative failure probability = 1-S(t)
Failure = function(alpha,beta,t){
  return(1-survival(alpha,beta,t))
}

#pdf (point probability) of failure = -S(t)dt
f_density = function(alpha,beta,t){
  return(-1*survival(alpha,beta,t)*(-alpha - beta*t)) 
}

dat = matrix(NA,80,5)
colnames(dat) = c('age','h(t)','S(t)','F(t)','f(t)')

for(t in 1:80){
  dat[t,'age'] = t
  dat[t,'h(t)'] = hazard(alpha,beta,t)
  dat[t,'S(t)'] = survival(alpha,beta,t)
  dat[t,'f(t)'] = f_density(alpha,beta,t)
  dat[t,'F(t)'] = Failure(alpha,beta,t)
}

#windows(height=300,width=300)

par(mfrow=c(3,2),oma=c(0,0,0,0),mar=c(3,5.5,3,3))
plot(dat[,c('age','h(t)')], type='l')
plot(dat[,c('age','S(t)')], type='l')
plot(dat[,c('age','f(t)')], type='l')
plot(dat[,c('age','F(t)')], type='l')

#random draw---from F(t) for time of death (cumulative probability) using inversion sampling

#1 draw u from U~(0,1)

u = runif(200)

#2 inversion sample from g(u), where F(t) = u (solve for t where F(t) = u) (survvial is same)

g = function(u,alpha,beta){
  #quadratic
  root = sqrt((2*alpha)^2 - 4*beta*2*log(1-u))
  pos = ((2*alpha) + root)/(-2*beta)
  neg = ((2*alpha) - root)/(-2*beta)

  return(c(pos,neg))
} 

draw = matrix(c(u,g(u,alpha,beta)),length(u),3)
colnames(draw)=c('u','pos','neg')

sim = round(draw[,'neg'])
hist(sim,freq=F, main='Simulated Age of Death',
      ,xlim=c(0,80))

pyears=sum(sim)
c.sim = c(NULL) #prep annual observations
#(observation for each year of life)

for(i in 1:length(sim)){
  #print(i)
  #print(1:sim[i])
  c.sim = c(c.sim,1:sim[i])
}

#hist(round(draw[,'neg']), xlim=c(0,80), ylim=c(0,1), freq=FALSE)

hist(c.sim,main=paste("Simulated S(x)\nPop=",length(sim)),
     ylab='Person Years')
