
# The script below contains all models examined in De Lisle et al. 'Extinction and temporal distribution
# of macroevolutionary bursts'. Note that included are two parameters not explored explicitly in the main 
# text (drift, white noise), but reported here for completeness. Incorporating drift in a 
# sensible way is not at all straight forward, as our demographic rescue model artificially holds 
# holds populations at low N when maladapted, inflating rates of evolution by drift (artificially 
# exacerbating apparent effects of extinction), and so we ignore drift (set to 0) in all our models. 



##############################################
######## DISPLACED OPTIMUM MODELS ############
##############################################
library(fields)
library(car)
library(matrixStats)
herit <- read.table("HansenSize.prn", header=T, na.strings="NA")
time <- 10^5  ##for all runs##

##DENSITY DEPENDENCE##

#### peak displacement, without extinction (lambda = 10^6)
nsim<-500
#startN <- 1000 #starting population size# silenced so starting N=K
w<- 3 #width of the adaptive landscape, from Estes and Arnold 2007 #
#h<-.4 #heritability, from Estes and Arnold 2007# silenced to randomize heritability from Hansen #
white<- 0 #white noise variance in optimum if any #
wmax<-1.5 #fitness at the optimum #
P<-1 #within population phenotypic variance #
thetav<- 169 #variance in the optimum through time #
drift <- 0 #drift occurs: yes = 1 no = 0 #
rescue <- 50
lambda<- 10^(-6) # poisson variance of displacement
crit<- 0 #criterion for displacement
plot(0, 0, xlab = expression(paste("Log"[10],"Time (generations)")), ylab=expression(paste("Phenotype"," (",sigma, ")")), xlim = c(0, log10(time)), ylim = c(-10, 10), type = "l", main=expression(paste("Displaced Optimum, ", lambda, " = 10"^"-6"))) #create axes#
lastZ6<- matrix(nrow=nsim,)
VDO1<- matrix(nrow = time/10, ncol = nsim)
for(i in 1:(nsim))
{
  h <- sample(herit$H, 1) # sampling from Hansen dataset on heritability
  K <- sample(seq(10561,12259),1)  # carrying capacity drawn from empirically estimated distribution in Reed et al 2003 #
  Z <- 0 # phenotypic mean #
  N <- K  # population size #
  O <- 0 # position of optimum #
  gen <- 1
  data <- matrix(, nrow=time, ncol=3)
  for(j in 1:time)
  {  	 		
    Z <- Z + (-h*(Z-(O+rnorm(1,0, white)+rnorm(1,0,h*drift/N)))*P/(P+w))   # phenotypic evolution towards optimum #
    N <- N*exp(log(wmax)*(1-N/K)-(((Z -(O+rnorm(1,0, white)))^2+P)/(2*(P+w))))  #L&L 1993 Eq 2 population size, see K&B 1997 Eq 7#
    
    rand<-rpois(1, lambda) # occurrence of displacement event with poisson probability #
    if (rand >crit){
      O <- O+ rnorm(1, mean = 0, thetav^0.5) # change in position of optimum conditioned on event occurring and white noise see E&A, Eq. 9 #
    }
    else {
      O <- O  # position of optimum if event does not occur #
    }
    
    gen <- gen + 1  
    data[j, 1] <- gen
    data[j, 2] <- Z
    data[j, 3] <- N
    if(N < 2) N = 50
  }
  lastZ6[i,] <- Z
  dataShort<-data[seq(1, nrow(data), 10), ]
  VDO1[, i] <- cbind(dataShort[,2])
  points(log10(dataShort[,1]), dataShort[, 2], type="l", col = rgb(red = 0, green = 0, blue = 1, alpha = 0.5))
}
varianceDo1<- matrix(nrow= nrow(VDO1), ncol = 4)
varianceDo1[,1]<-seq(1, nrow(data), 10)
varianceDo1[,2]<-rowVars(VDO1, na.rm = TRUE)
varianceDo1[,3]<-varianceDo1[,2] * (nsim-1) / qchisq(0.05/2, (nsim-1), lower.tail = FALSE)
varianceDo1[,4]<-varianceDo1[,2] * (nsim-1) / qchisq(1 - 0.05/2, (nsim-1), lower.tail = FALSE)


####  peak displacement, with extinction (lambda = 10^6)
nsim<-500
#startN <- 1000 #starting population size# silenced so starting N=K
w<- 3 #width of the adaptive landscape, from Estes and Arnold 2007 #
#h<-.4 #heritability, from Estes and Arnold 2007# silenced to randomize heritability from Hansen #
white<- 0 #white noise variance in optimum if any #
wmax<-1.5 #fitness at the optimum #
P<-1 #within population phenotypic variance #
thetav<- 169 #variance in the optimum through time #
drift <- 0 #drift occurs: yes = 1 no = 0 #
rescue <- 50
lambda<- 10^(-6)# poisson variance of displacement
crit<- 0 #criterion for displacement
extinctions <- 0
lastZE6 <- matrix(nrow=nsim,)
VDOE1<- matrix(nrow = time/10, ncol = nsim)
for(i in 1:(nsim))
{
  h <- sample(herit$H, 1) # sampling from Hansen dataset on heritability
  K <- sample(seq(10561,12259),1)  # carrying capacity drawn from empirically estimated distribution in Reed et al 2003 #
  Z <- 0 # phenotypic mean #
  N <- K  # population size #
  O <- 0 # position of optimum #
  gen <- 1
  data <- matrix(, nrow=time, ncol=3)
  for(j in 1:time)
  {  	 		
    Z <- Z + (-h*(Z-(O+rnorm(1,0, white)+rnorm(1,0,h*drift/N)))*P/(P+w))   # phenotypic evolution towards optimum #
    N <- N*exp(log(wmax)*(1-N/K)-(((Z -(O+rnorm(1,0, white)))^2+P)/(2*(P+w))))  #L&L 1993 Eq 2 population size, see K&B 1997 Eq 7#
    
    rand<-rpois(1, lambda) # occurrence of displacement event with poisson probability #
    if (rand >crit){
      O <- O+ rnorm(1, mean = 0, thetav^0.5) # change in position of optimum conditioned on event occurring and white noise see E&A, Eq. 9 #
    }
    else {
      O <- O # position of optimum if event does not occur #
    }
    
    gen <- gen + 1  
    data[j, 1] <- gen
    data[j, 2] <- Z
    data[j, 3] <- N
    if(N < 2) extinctions=extinctions+1 
    if(N < 2) {break}
  }
  lastZE6[i,] <- Z
  dataShort<-data[seq(1, nrow(data), 10), ]
  VDOE1[, i] <- cbind(dataShort[,2])
  points(log10(dataShort[,1]), dataShort[, 2], type="l", col = rgb(red = 1, green = 0, blue = 0))
}
varianceDoe1<- matrix(nrow= nrow(VDOE1), ncol = 4)
varianceDoe1[,1]<-seq(1, nrow(data), 10)
varianceDoe1[,2]<-rowVars(VDOE1, na.rm = TRUE)
varianceDoe1[,3]<-varianceDoe1[,2] * (nsim-1) / qchisq(0.05/2, (nsim-1), lower.tail = FALSE)
varianceDoe1[,4]<-varianceDoe1[,2] * (nsim-1) / qchisq(1 - 0.05/2, (nsim-1), lower.tail = FALSE)
plot(0, 0, xlim = c(0, time), ylim = c(0, 10), type = "l", xlab = expression(paste("Time (generations)")), ylab=expression(paste("Phenotypic Divergence")), main=expression(paste("DO, ", lambda, " = 10^6")))
points(varianceDo1[,1], varianceDo1[, 2], type="l", col = rgb(red = 0, green = 0, blue = 1, alpha = 0.5))
polygon(x = c(varianceDo1[,1], rev(varianceDo1[,1])),
        y = c(varianceDo1[,3], 
              rev(varianceDo1[,4])),
        col =  adjustcolor(rgb(red = 0, green = 0, blue = 1, alpha = 0.1)), border = NA)
points(varianceDoe1[,1], varianceDoe1[, 2], type="l", col = rgb(red = 1, green = 0, blue = 1, alpha = 0.5))
polygon(x = c(varianceDoe1[,1], rev(varianceDoe1[,1])),
        y = c(varianceDoe1[,3], 
              rev(varianceDoe1[,4])),
        col =  adjustcolor(rgb(red = 1, green = 0, blue = 0, alpha = 0.1)), border = NA)





dze6<-density(lastZE6)
dz6<-density(lastZ6)
plot(dze6, xlim = c(-10, 10), main=expression(paste("Phenotypic divergence after ","10"^"5"," generations")), xlab=expression(paste("Phenotype"," (",sigma, ")")))
polygon(dze6, col=rgb(1,0,0,0.5), border= "black")
lines(dz6)
polygon(dz6, col=rgb(0,0,1,0.4), border= "black")
dat <- data.frame(c(lastZE6, lastZ6), lines = rep(c("ext", "noe"), each = 500))
leveneTest(dat[,1],dat[,2] )



##DENSITY INDEPENDENCE##

#### peak displacement, without extinction (lambda = 10^6)
nsim<-500
#startN <- 1000 #starting population size# silenced so starting N=K
w<- 3 #width of the adaptive landscape, from Estes and Arnold 2007 #
#h<-.4 #heritability, from Estes and Arnold 2007# silenced to randomize heritability from Hansen #
white<- 0 #white noise variance in optimum if any #
wmax<-1.5 #fitness at the optimum #
P<-1 #within population phenotypic variance #
thetav<- 169 #variance in the optimum through time #
drift <- 0 #drift occurs: yes = 1 no = 0 #
rescue <- 50
lambda<- 10^(-6) # poisson variance of displacement
crit<- 0 #criterion for displacement
plot(0, 0, xlab = expression(paste("Log"[10],"Time (generations)")), ylab=expression(paste("Phenotype"," (",sigma, ")")), xlim = c(0, log10(time)), ylim = c(-10, 10), type = "l", main=expression(paste("Displaced Optimum, ", lambda, " = 10"^"-6"))) #create axes#
lastZ6<- matrix(nrow=nsim,)
VDO1<- matrix(nrow = time/10, ncol = nsim)
for(i in 1:(nsim))
{
  h <- sample(herit$H, 1) # sampling from Hansen dataset on heritability
  K <- sample(seq(10561,12259),1)  # carrying capacity drawn from empirically estimated distribution in Reed et al 2003 #
  Z <- 0 # phenotypic mean #
  N <- K  # population size #
  O <- 0 # position of optimum #
  gen <- 1
  data <- matrix(, nrow=time, ncol=3)
  for(j in 1:time)
  {  	 		
    Z <- Z + (-h*(Z-(O+rnorm(1,0, white)+rnorm(1,0,h*drift/N)))*P/(P+w))   # phenotypic evolution towards optimum #
    N <- N*wmax*exp(-((Z -(O+rnorm(1,0, white)))^2)/(2*(P+w)))  
    rand<-rpois(1, lambda) # occurrence of displacement event with poisson probability #
    if (rand >crit){
      O <- O+ rnorm(1, mean = 0, thetav^0.5) # change in position of optimum conditioned on event occurring and white noise see E&A, Eq. 9 #
    }
    else {
      O <- O  # position of optimum if event does not occur #
    }
    
    gen <- gen + 1  
    data[j, 1] <- gen
    data[j, 2] <- Z
    data[j, 3] <- N
    if(N > 10000000) N = 10000000
    if(N < 50) N = 50
  }
  lastZ6[i,] <- Z
  dataShort<-data[seq(1, nrow(data), 10), ]
  VDO1[, i] <- cbind(dataShort[,2])
  points(log10(dataShort[,1]), dataShort[, 2], type="l", col = rgb(red = 0, green = 0, blue = 1, alpha = 0.5))
}
varianceDo1<- matrix(nrow= nrow(VDO1), ncol = 4)
varianceDo1[,1]<-seq(1, nrow(data), 10)
varianceDo1[,2]<-rowVars(VDO1, na.rm = TRUE)
varianceDo1[,3]<-varianceDo1[,2] * (nsim-1) / qchisq(0.05/2, (nsim-1), lower.tail = FALSE)
varianceDo1[,4]<-varianceDo1[,2] * (nsim-1) / qchisq(1 - 0.05/2, (nsim-1), lower.tail = FALSE)


####  peak displacement, with extinction (lambda = 10^6)
nsim<-500
#startN <- 1000 #starting population size# silenced so starting N=K
w<- 3 #width of the adaptive landscape, from Estes and Arnold 2007 #
#h<-.4 #heritability, from Estes and Arnold 2007# silenced to randomize heritability from Hansen #
white<- 0 #white noise variance in optimum if any #
wmax<-1.5 #fitness at the optimum #
P<-1 #within population phenotypic variance #
thetav<- 169 #variance in the optimum through time #
drift <- 0 #drift occurs: yes = 1 no = 0 #
rescue <- 50
lambda<- 10^(-6)# poisson variance of displacement
crit<- 0 #criterion for displacement
extinctions <- 0
lastZE6 <- matrix(nrow=nsim,)
VDOE1<- matrix(nrow = time/10, ncol = nsim)
for(i in 1:(nsim))
{
  h <- sample(herit$H, 1) # sampling from Hansen dataset on heritability
  K <- sample(seq(10561,12259),1)  # carrying capacity drawn from empirically estimated distribution in Reed et al 2003 #
  Z <- 0 # phenotypic mean #
  N <- K  # population size #
  O <- 0 # position of optimum #
  gen <- 1
  data <- matrix(, nrow=time, ncol=3)
  for(j in 1:time)
  {  	 		
    Z <- Z + (-h*(Z-(O+rnorm(1,0, white)+rnorm(1,0,h*drift/N)))*P/(P+w))   # phenotypic evolution towards optimum #
    N <- N*wmax*exp(-((Z -(O+rnorm(1,0, white)))^2)/(2*(P+w)))  
    rand<-rpois(1, lambda) # occurrence of displacement event with poisson probability #
    if (rand >crit){
      O <- O+ rnorm(1, mean = 0, thetav^0.5) # change in position of optimum conditioned on event occurring and white noise see E&A, Eq. 9 #
    }
    else {
      O <- O # position of optimum if event does not occur #
    }
    
    gen <- gen + 1  
    data[j, 1] <- gen
    data[j, 2] <- Z
    data[j, 3] <- N
    if(N > 10000000) N = 10000000
    if(N < 50) extinctions=extinctions+1 
    if(N < 50) {break}
  }
  lastZE6[i,] <- Z
  dataShort<-data[seq(1, nrow(data), 10), ]
  VDOE1[, i] <- cbind(dataShort[,2])
  points(log10(dataShort[,1]), dataShort[, 2], type="l", col = rgb(red = 1, green = 0, blue = 0))
}
varianceDoe1<- matrix(nrow= nrow(VDOE1), ncol = 4)
varianceDoe1[,1]<-seq(1, nrow(data), 10)
varianceDoe1[,2]<-rowVars(VDOE1, na.rm = TRUE)
varianceDoe1[,3]<-varianceDoe1[,2] * (nsim-1) / qchisq(0.05/2, (nsim-1), lower.tail = FALSE)
varianceDoe1[,4]<-varianceDoe1[,2] * (nsim-1) / qchisq(1 - 0.05/2, (nsim-1), lower.tail = FALSE)
plot(0, 0, xlim = c(0, time), ylim = c(0, 10), type = "l", xlab = expression(paste("Time (generations)")), ylab=expression(paste("Phenotypic Divergence")), main=expression(paste("DO, ", lambda, " = 10^6")))
points(varianceDo1[,1], varianceDo1[, 2], type="l", col = rgb(red = 0, green = 0, blue = 1, alpha = 0.5))
polygon(x = c(varianceDo1[,1], rev(varianceDo1[,1])),
        y = c(varianceDo1[,3], 
              rev(varianceDo1[,4])),
        col =  adjustcolor(rgb(red = 0, green = 0, blue = 1, alpha = 0.1)), border = NA)
points(varianceDoe1[,1], varianceDoe1[, 2], type="l", col = rgb(red = 1, green = 0, blue = 1, alpha = 0.5))
polygon(x = c(varianceDoe1[,1], rev(varianceDoe1[,1])),
        y = c(varianceDoe1[,3], 
              rev(varianceDoe1[,4])),
        col =  adjustcolor(rgb(red = 1, green = 0, blue = 0, alpha = 0.1)), border = NA)





dze6<-density(lastZE6)
dz6<-density(lastZ6)
plot(dze6, xlim = c(-10, 10), main=expression(paste("Phenotypic divergence after ","10"^"5"," generations")), xlab=expression(paste("Phenotype"," (",sigma, ")")))
polygon(dze6, col=rgb(1,0,0,0.5), border= "black")
lines(dz6)
polygon(dz6, col=rgb(0,0,1,0.4), border= "black")
dat <- data.frame(c(lastZE6, lastZ6), lines = rep(c("ext", "noe"), each = 500))
leveneTest(dat[,1],dat[,2] )



##############################################
#### Single Displacement (FIGURE 4) ##########
##############################################
time <- 100 
nsim<-5000
#startN <- 1000 #starting population size# silenced so starting N=K
w<- 3 #width of the adaptive landscape, from Estes and Arnold 2007 #
#h<-.4 #heritability, from Estes and Arnold 2007# silenced to randomize heritability from Hansen #
white<- 0 #white noise variance in optimum if any #
wmax<-1.5 #fitness at the optimum #
P<-1 #within population phenotypic variance #
thetav<- 9 #variance in the optimum through time #
drift <- 0 #drift occurs: yes = 1 no = 0 #
rescue <- 50
lambda<- 10^(-7)# poisson variance of displacement
crit<- 0 #criterion for displacement
data_newfig <- matrix(nrow=nsim,ncol = 3)
colnames(data_newfig)<-paste(c('survival', 'final_d', 'initial_d'))
for(i in 1:(nsim))
{
  h <- sample(herit$H, 1) # sampling from Hansen dataset on heritability
  K <- sample(seq(10561,12259),1)  # carrying capacity drawn from empirically estimated distribution in Reed et al 2003 #
  Z <- 0 # phenotypic mean #
  N <- K  # population size #
  O <- rnorm(1, mean = 0, thetav^0.5)
  gen <- 1
  data <- matrix(, nrow=time, ncol=3)
  d = abs(Z-O)
  survival <- 1
  for(j in 1:time)
  {  	 		
    Z <- Z + (-h*(Z-(O+rnorm(1,0, white)+rnorm(1,0,h*drift/N)))*P/(P+w))   # phenotypic evolution towards optimum #
    N <- N*exp(log(wmax)*(1-N/K)-(((Z -(O+rnorm(1,0, white)))^2+P)/(2*(P+w))))  #L&L 1993 Eq 2 population size, see K&B 1997 Eq 7#
    gen <- gen + 1  
    data[j, 1] <- gen
    data[j, 2] <- Z
    data[j, 3] <- N
    if(N < 2) survival=0
    if(N < 2) d = abs(Z-O)
    if(N < 2) {break}
  }
  data_newfig[i, 1] <- survival
  data_newfig[i, 2] <- abs(Z-O)
  data_newfig[i, 3] <- abs(O)
}
data_newfig<-data.frame(data_newfig)
fit = glm(survival ~ initial_d, data=data_newfig, family=binomial)
newdat <- data.frame(initial_d=seq(min(data_newfig$initial_d), max(data_newfig$initial_d),len=nsim))
newdat$survival = predict(fit, data=newdat, type="response")
par(mfrow=c(1,2))
Corner_text <- function(text, location="topright"){
  legend(location,legend=text, bty ="n", pch=NA) 
}
plot(survival ~ initial_d, data= data_newfig, col="black", xlab = expression(paste("Magnitude of peak displacement")), ylab = expression(paste("survival probability")))
Corner_text(text='A')
curve(predict(fit, data.frame(initial_d=x), type="response"), add=TRUE, col="red", lwd=2)
extinct<-subset(data_newfig, survival ==0)
survive<-subset(data_newfig, survival ==1)
p1 <- hist(extinct$final_d, main = paste(""), xlab = expression(paste("|",bar(Z),"-", theta,"| at extinction")))                     
plot(p1, col=rgb(1,0,0,1/4), xlim=c(0,10), add=T)
Corner_text(text='B')






##############################################
######## BROWNIAN MOTION MODELS ##############
##############################################

time <- 10^5  ##for all runs##

### a random walk of the optimum, thetaSD =0.1, without extinction
nsim<-500
#startN <- 1000 #starting population size# silenced so starting N=K
w<- 3 #width of the adaptive landscape, from Estes and Arnold 2007 #
#h<-.4 #heritability, from Estes and Arnold 2007# silenced to randomize heritability from Hansen #
white<- 0 #white noise variance in optimum if any #
wmax<-1.5 #fitness at the optimum #
P<-1 #within population phenotypic variance #
thetav<- 0.01 #variance in the optimum through time #
drift <- 0 #drift occurs: yes = 1 no = 0 #
rescue <- 50
lastZBM1 <- matrix(nrow=nsim,)
VBM1<- matrix(nrow = time/10, ncol = nsim)
counter <- 1
for(i in 1:(nsim))
{
  h <- sample(herit$H, 1) # sampling from Hansen dataset on heritability
  K <- sample(seq(10561,12259),1)  # carrying capacity drawn from empirically estimated distribution in Reed et al 2003 #
  Z <- 0 # phenotypic mean #
  N <- K  # population size #
  O <- 0 # position of optimum #
  gen <- 1
  data <- matrix(, nrow=time, ncol=3)
  for(j in 1:time)
  {  	 		
    Z <- Z + (-h*(Z-(O+rnorm(1,0, white)+rnorm(1,0,h*drift/N)))*P/(P+w))   # phenotypic evolution towards optimum #
    N <- N*exp(log(wmax)*(1-N/K)-(((Z -(O+rnorm(1,0, white)))^2+P)/(2*(P+w))))  #L&L 1993 Eq 2 population size, see K&B 1997 Eq 7#
    O <- O + rnorm(1, mean = 0, sd = thetav^0.5) # position of optimum, with displacement  and white noise see E&A, Eq. 9 #
    gen <- gen + 1  
    data[j, 1] <- gen
    data[j, 2] <- Z
    data[j, 3] <- N
    if(N < 2) N = rescue
  }
  lastZBM1[i,] <- Z
  dataShort<-data[seq(1, nrow(data), 10), ]
  VBM1[, i] <- cbind(dataShort[,2])
}
varianceBm1<- matrix(nrow= nrow(VBM1), ncol = 5)
varianceBm1[,1]<-seq(1, nrow(data), 10)
varianceBm1[,2]<-rowVars(VBM1)
varianceBm1[,5]<-rowSums(!is.na(VBM1))
varianceBm1[,3]<-varianceBm1[,2] * (varianceBm1[,5]-1) / qchisq(0.05/2, (varianceBm1[,5]-1), lower.tail = FALSE)
varianceBm1[,4]<-varianceBm1[,2] * (varianceBm1[,5]-1) / qchisq(1 - 0.05/2, (varianceBm1[,5]-1), lower.tail = FALSE)


### a random walk of the optimum, thetaSD =0.1, with extinction
nsim<-500
#startN <- 1000 #starting population size# silenced so starting N=K
w<- 3 #width of the adaptive landscape, from Estes and Arnold 2007 #
#h<-.4 #heritability, from Estes and Arnold 2007# silenced to randomize heritability from Hansen #
white<- 0 #white noise variance in optimum if any #
wmax<-1.5 #fitness at the optimum #
P<-1 #within population phenotypic variance #
thetav<- 0.01 #variance in the optimum through time #
drift <- 0 #drift occurs: yes = 1 no = 0 #
rescue <- 50
#plot(0, 0, xlab = expression(paste("Log"[10],"Time (generations)")), ylab=expression(paste("Phenotype"," (",sigma, ")")), xlim = c(0, log10(time)), ylim = c(-10, 10), type = "l", main=expression(paste("BME, ", sigma[theta], " = 2.5"^"-9"))) #create axes#
extinctions <- 0
lastZBME1 <- matrix(nrow=nsim,)
VBME1<- matrix(nrow = time/10, ncol = nsim)
counter <- 1
for(i in 1:(nsim))
{
  h <- sample(herit$H, 1) # sampling from Hansen dataset on heritability
  K <- sample(seq(10561,12259),1)  # carrying capacity drawn from empirically estimated distribution in Reed et al 2003 #
  Z <- 0 # phenotypic mean #
  N <- K  # population size #
  O <- 0 # position of optimum #
  gen <- 1
  data <- matrix(, nrow=time, ncol=3)
  
  
  for(j in 1:time)
  {  	 		
    Z <- Z + (-h*(Z-(O+rnorm(1,0, white)+rnorm(1,0,h*drift/N)))*P/(P+w))   # phenotypic evolution towards optimum #
    N <- N*exp(log(wmax)*(1-N/K)-(((Z -(O+rnorm(1,0, white)))^2+P)/(2*(P+w))))  #L&L 1993 Eq 2 population size, see K&B 1997 Eq 7#
    O <- O + rnorm(1, mean = 0, sd = thetav^0.5) # position of optimum, with displacement  and white noise see E&A, Eq. 9 #
    gen <- gen + 1  
    data[j, 1] <- gen
    data[j, 2] <- Z
    data[j, 3] <- N
    if(N < 2) extinctions=extinctions+1 
    if(N <2) {break}
  }
  lastZBME1[i,] <- Z
  dataShort<-data[seq(1, nrow(data), 10), ]
  VBME1[, i] <- cbind(dataShort[,2])
}
varianceBme1<- matrix(nrow= nrow(VBME1), ncol = 5)
varianceBme1[,1]<-seq(1, nrow(data), 10)
varianceBme1[,2]<-rowVars(VBME1,na.rm = TRUE)
varianceBme1[,5]<-rowSums(!is.na(VBME1))
varianceBme1[,3]<-varianceBme1[,2] * (varianceBme1[,5]-1) / qchisq(0.05/2, (varianceBme1[,5]-1), lower.tail = FALSE)
varianceBme1[,4]<-varianceBme1[,2] * (varianceBme1[,5]-1) / qchisq(1 - 0.05/2, (varianceBme1[,5]-1), lower.tail = FALSE)

plot(0, 0, xlim = c(0, time), ylim = c(0, 1500), type = "l", xlab = expression(paste("Time (generations)")), ylab=expression(paste("Phenotypic Divergence")), main=expression(paste("BM, ", sigma[theta], " = 0.1")))
points(varianceBm1[,1], varianceBm1[, 2], type="l", col = rgb(red = 0, green = 0, blue = 1, alpha = 0.5))
polygon(x = c(varianceBm1[,1], rev(varianceBm1[,1])),
        y = c(varianceBm1[,3], 
              rev(varianceBm1[,4])),
        col =  adjustcolor(rgb(red = 0, green = 0, blue = 1, alpha = 0.1)), border = NA)
points(varianceBme1[,1], varianceBme1[, 2], type="l", col = rgb(red = 1, green = 0, blue = 1, alpha = 0.5))
polygon(x = c(varianceBme1[,1], rev(varianceBme1[,1])),
        y = c(varianceBme1[,3], 
              rev(varianceBme1[,4])),
        col =  adjustcolor(rgb(red = 1, green = 0, blue = 0, alpha = 0.1)), border = NA)


BMze2<-density(lastZBME1)
BMz2<-density(lastZBM1)
plot(BMze2, xlim = c(-100, 100), main=expression(paste("Phenotypic divergence after ","10"^"5"," generations")), xlab=expression(paste("Phenotype"," (",sigma, ")")))
polygon(BMze2, col=rgb(1,0,0,0.5), border= "black")
lines(BMz2)
polygon(BMz2, col=rgb(0,0,1,0.4), border= "black")
dat <- data.frame(c(lastZBME1, lastZBM1), lines = rep(c("ext", "noe"), each = 500))
leveneTest(dat[,1],dat[,2] )






