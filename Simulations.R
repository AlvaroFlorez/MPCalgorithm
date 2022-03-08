# Run Functions.R before using this file.

#### Section 5: Effect of adding non-informative surrogates
Sigma.4 = matrix(NA,8,8)
diag(Sigma.4) = 1
Sigma.4[1,3] = Sigma.4[2,4]  = 0.95
Sigma.4[1,5]=Sigma.4[1,7]=Sigma.4[2,6]=Sigma.4[2,8]=Sigma.4[3,5]=Sigma.4[3,7]=Sigma.4[4,6]=Sigma.4[4,8]=0
Sigma.4[5,7] = Sigma.4[6,8] =  0.8
Sigma.4 = as.matrix(forceSymmetric(Sigma.4))

## using the GRS algorithm
ICA.BF = ICA.ContCont.MultS(10000, 100, Sigma.4,Show.Progress=F) # time consuming! ~60min
range(ICA.BF$R2_H)

MIN.ICA.R = ICA.BF$Lower.Dig.Corrs.Sigma[ICA.BF$R2_H == min(ICA.BF$R2_H),]
R.min = matrix(NA,8,8)
R.min[lower.tri(R.min)] = unlist(MIN.ICA.R)
diag(R.min) = 1
R.min = forceSymmetric(R.min,uplo=F)

### Figure 1 Relationship between the simulated ICA and four unidentifiable correlations(800 x 700)
par(mfrow=c(2,2))
plot(ICA.BF$Lower.Dig.Corrs.Sigma[,7],ICA.BF$R2_H,col='#999999',ylab='ICA',xlab=expression(rho),main='(a)')
lw1 = loess(ICA.BF$R2_H ~ ICA.BF$Lower.Dig.Corrs.Sigma[,7])
j <- order(ICA.BF$Lower.Dig.Corrs.Sigma[,7])
points(ICA.BF$Lower.Dig.Corrs.Sigma[ICA.BF$R2_H==min(ICA.BF$R2_H),7],ICA.BF$R2_H[ICA.BF$R2_H==min(ICA.BF$R2_H)], 
       col=1,pch=15,cex=1.2)
lines(ICA.BF$Lower.Dig.Corrs.Sigma[j,7],lw1$fitted[j],lty=2,lwd=2)

plot(ICA.BF$Lower.Dig.Corrs.Sigma[,19],ICA.BF$R2_H,col='#999999',ylab='ICA',xlab=expression(rho),main='(b)')
lw1 = loess(ICA.BF$R2_H ~ ICA.BF$Lower.Dig.Corrs.Sigma[,19])
j <- order(ICA.BF$Lower.Dig.Corrs.Sigma[,19])
points(ICA.BF$Lower.Dig.Corrs.Sigma[ICA.BF$R2_H==min(ICA.BF$R2_H),19],ICA.BF$R2_H[ICA.BF$R2_H==min(ICA.BF$R2_H)], 
       col=1,pch=15,cex=1.2)
lines(ICA.BF$Lower.Dig.Corrs.Sigma[j,19],lw1$fitted[j],lty=2,lwd=2)

plot(ICA.BF$Lower.Dig.Corrs.Sigma[,3],ICA.BF$R2_H,col='#999999',ylab='ICA',xlab=expression(rho),main='(c)')
lw1 = loess(ICA.BF$R2_H ~ ICA.BF$Lower.Dig.Corrs.Sigma[,3])
j <- order(ICA.BF$Lower.Dig.Corrs.Sigma[,3])
points(ICA.BF$Lower.Dig.Corrs.Sigma[ICA.BF$R2_H==min(ICA.BF$R2_H),3],ICA.BF$R2_H[ICA.BF$R2_H==min(ICA.BF$R2_H)], 
       col=1,pch=15,cex=1.2)
lines(ICA.BF$Lower.Dig.Corrs.Sigma[j,3],lw1$fitted[j],lty=2,lwd=2)

plot(ICA.BF$Lower.Dig.Corrs.Sigma[,26],ICA.BF$R2_H,col='#999999',ylab='ICA',xlab=expression(rho),main='(d)')
lw1 = loess(ICA.BF$R2_H ~ ICA.BF$Lower.Dig.Corrs.Sigma[,26])
j <- order(ICA.BF$Lower.Dig.Corrs.Sigma[,26])
points(ICA.BF$Lower.Dig.Corrs.Sigma[ICA.BF$R2_H==min(ICA.BF$R2_H),26],ICA.BF$R2_H[ICA.BF$R2_H==min(ICA.BF$R2_H)], 
       col=1,pch=15,cex=1.2)
lines(ICA.BF$Lower.Dig.Corrs.Sigma[j,26],lw1$fitted[j],lty=2,lwd=2)


### Figures supplementary materials
# other correlations useful-useful correlations (800 x 300)
# Figure S1
par(mfrow=c(1,3))
plot(ICA.BF$Lower.Dig.Corrs.Sigma[,1],ICA.BF$R2_H,col='#999999',
     ylab='ICA',xlab=expression(rho),main='(a)',cex.lab=1.6,cex.axis=1.2,cex.main=1.6)
lw1 = loess(ICA.BF$R2_H ~ ICA.BF$Lower.Dig.Corrs.Sigma[,1])
j <- order(ICA.BF$Lower.Dig.Corrs.Sigma[,1])
points(ICA.BF$Lower.Dig.Corrs.Sigma[ICA.BF$R2_H==min(ICA.BF$R2_H),1],ICA.BF$R2_H[ICA.BF$R2_H==min(ICA.BF$R2_H)], 
       col=1,pch=15,cex=1.2)
lines(ICA.BF$Lower.Dig.Corrs.Sigma[j,1],lw1$fitted[j],lty=2,lwd=2)

plot(ICA.BF$Lower.Dig.Corrs.Sigma[,8],ICA.BF$R2_H,col='#999999',
     ylab='ICA',xlab=expression(rho),main='(b)',cex.lab=1.6,cex.axis=1.2,cex.main=1.6)
lw1 = loess(ICA.BF$R2_H ~ ICA.BF$Lower.Dig.Corrs.Sigma[,8])
j <- order(ICA.BF$Lower.Dig.Corrs.Sigma[,8])
points(ICA.BF$Lower.Dig.Corrs.Sigma[ICA.BF$R2_H==min(ICA.BF$R2_H),8],ICA.BF$R2_H[ICA.BF$R2_H==min(ICA.BF$R2_H)], 
       col=1,pch=15,cex=1.2)
lines(ICA.BF$Lower.Dig.Corrs.Sigma[j,8],lw1$fitted[j],lty=2,lwd=2)

plot(ICA.BF$Lower.Dig.Corrs.Sigma[,14],ICA.BF$R2_H,col='#999999',
     ylab='ICA',xlab=expression(rho),main='(c)',cex.lab=1.6,cex.axis=1.2,cex.main=1.6)
lw1 = loess(ICA.BF$R2_H ~ ICA.BF$Lower.Dig.Corrs.Sigma[,14])
j <- order(ICA.BF$Lower.Dig.Corrs.Sigma[,14])
points(ICA.BF$Lower.Dig.Corrs.Sigma[ICA.BF$R2_H==min(ICA.BF$R2_H),14],ICA.BF$R2_H[ICA.BF$R2_H==min(ICA.BF$R2_H)], 
       col=1,pch=15,cex=1.2)
lines(ICA.BF$Lower.Dig.Corrs.Sigma[j,14],lw1$fitted[j],lty=2,lwd=2)

# other correlations with useful-non informative correlations (1000 x 600)
# Figure S2  
par(mfrow=c(2,4))
plot(ICA.BF$Lower.Dig.Corrs.Sigma[,5],ICA.BF$R2_H,col='#999999',ylab='ICA',
     xlab=expression(rho),main='(a)',cex.lab=1.8,cex.axis=1.4,cex.main=1.8)
lw1 = loess(ICA.BF$R2_H ~ ICA.BF$Lower.Dig.Corrs.Sigma[,5])
j <- order(ICA.BF$Lower.Dig.Corrs.Sigma[,5])
points(ICA.BF$Lower.Dig.Corrs.Sigma[ICA.BF$R2_H==min(ICA.BF$R2_H),5],ICA.BF$R2_H[ICA.BF$R2_H==min(ICA.BF$R2_H)], 
       col=1,pch=15,cex=1.2)
lines(ICA.BF$Lower.Dig.Corrs.Sigma[j,5],lw1$fitted[j],lty=2,lwd=2)

plot(ICA.BF$Lower.Dig.Corrs.Sigma[,10],ICA.BF$R2_H,col='#999999',ylab='ICA',
     xlab=expression(rho),main='(b)',cex.lab=1.8,cex.axis=1.4,cex.main=1.8)
lw1 = loess(ICA.BF$R2_H ~ ICA.BF$Lower.Dig.Corrs.Sigma[,10])
j <- order(ICA.BF$Lower.Dig.Corrs.Sigma[,10])
points(ICA.BF$Lower.Dig.Corrs.Sigma[ICA.BF$R2_H==min(ICA.BF$R2_H),10],ICA.BF$R2_H[ICA.BF$R2_H==min(ICA.BF$R2_H)], 
       col=1,pch=15,cex=1.2)
lines(ICA.BF$Lower.Dig.Corrs.Sigma[j,10],lw1$fitted[j],lty=2,lwd=2)

plot(ICA.BF$Lower.Dig.Corrs.Sigma[,12],ICA.BF$R2_H,col='#999999',ylab='ICA',
     xlab=expression(rho),main='(c)',cex.lab=1.8,cex.axis=1.4,cex.main=1.8)
lw1 = loess(ICA.BF$R2_H ~ ICA.BF$Lower.Dig.Corrs.Sigma[,12])
j <- order(ICA.BF$Lower.Dig.Corrs.Sigma[,12])
points(ICA.BF$Lower.Dig.Corrs.Sigma[ICA.BF$R2_H==min(ICA.BF$R2_H),12],ICA.BF$R2_H[ICA.BF$R2_H==min(ICA.BF$R2_H)], 
       col=1,pch=15,cex=1.2)
lines(ICA.BF$Lower.Dig.Corrs.Sigma[j,12],lw1$fitted[j],lty=2,lwd=2)

plot(ICA.BF$Lower.Dig.Corrs.Sigma[,16],ICA.BF$R2_H,col='#999999',ylab='ICA',
     xlab=expression(rho),main='(d)',cex.lab=1.8,cex.axis=1.4,cex.main=1.8)
lw1 = loess(ICA.BF$R2_H ~ ICA.BF$Lower.Dig.Corrs.Sigma[,16])
j <- order(ICA.BF$Lower.Dig.Corrs.Sigma[,16])
points(ICA.BF$Lower.Dig.Corrs.Sigma[ICA.BF$R2_H==min(ICA.BF$R2_H),16],ICA.BF$R2_H[ICA.BF$R2_H==min(ICA.BF$R2_H)], 
       col=1,pch=15,cex=1.2)
lines(ICA.BF$Lower.Dig.Corrs.Sigma[j,16],lw1$fitted[j],lty=2,lwd=2)

plot(ICA.BF$Lower.Dig.Corrs.Sigma[,18],ICA.BF$R2_H,col='#999999',ylab='ICA',
     xlab=expression(rho),main='(e)',cex.lab=1.8,cex.axis=1.4,cex.main=1.8)
lw1 = loess(ICA.BF$R2_H ~ ICA.BF$Lower.Dig.Corrs.Sigma[,18])
j <- order(ICA.BF$Lower.Dig.Corrs.Sigma[,18])
points(ICA.BF$Lower.Dig.Corrs.Sigma[ICA.BF$R2_H==min(ICA.BF$R2_H),18],ICA.BF$R2_H[ICA.BF$R2_H==min(ICA.BF$R2_H)], 
       col=1,pch=15,cex=1.2)
lines(ICA.BF$Lower.Dig.Corrs.Sigma[j,18],lw1$fitted[j],lty=2,lwd=2)

plot(ICA.BF$Lower.Dig.Corrs.Sigma[,19],ICA.BF$R2_H,col='#999999',ylab='ICA',
     xlab=expression(rho),main='(f)',cex.lab=1.8,cex.axis=1.4,cex.main=1.8)
lw1 = loess(ICA.BF$R2_H ~ ICA.BF$Lower.Dig.Corrs.Sigma[,19])
j <- order(ICA.BF$Lower.Dig.Corrs.Sigma[,19])
points(ICA.BF$Lower.Dig.Corrs.Sigma[ICA.BF$R2_H==min(ICA.BF$R2_H),19],ICA.BF$R2_H[ICA.BF$R2_H==min(ICA.BF$R2_H)], 
       col=1,pch=15,cex=1.2)
lines(ICA.BF$Lower.Dig.Corrs.Sigma[j,19],lw1$fitted[j],lty=2,lwd=2)

plot(ICA.BF$Lower.Dig.Corrs.Sigma[,21],ICA.BF$R2_H,col='#999999',ylab='ICA',
     xlab=expression(rho),main='(g)',cex.lab=1.8,cex.axis=1.4,cex.main=1.8)
lw1 = loess(ICA.BF$R2_H ~ ICA.BF$Lower.Dig.Corrs.Sigma[,21])
j <- order(ICA.BF$Lower.Dig.Corrs.Sigma[,21])
points(ICA.BF$Lower.Dig.Corrs.Sigma[ICA.BF$R2_H==min(ICA.BF$R2_H),21],ICA.BF$R2_H[ICA.BF$R2_H==min(ICA.BF$R2_H)], 
       col=1,pch=15,cex=1.2)
lines(ICA.BF$Lower.Dig.Corrs.Sigma[j,21],lw1$fitted[j],lty=2,lwd=2)


# other correlations with noise-noise correlations (800 x 300)
# Figure S3
par(mfrow=c(1,3))

plot(ICA.BF$Lower.Dig.Corrs.Sigma[,23],ICA.BF$R2_H,col='#999999',ylab='ICA',
     xlab=expression(rho),main='(a)',cex.lab=1.6,cex.axis=1.2,cex.main=1.6)
lw1 = loess(ICA.BF$R2_H ~ ICA.BF$Lower.Dig.Corrs.Sigma[,23])
j <- order(ICA.BF$Lower.Dig.Corrs.Sigma[,23])
points(ICA.BF$Lower.Dig.Corrs.Sigma[ICA.BF$R2_H==min(ICA.BF$R2_H),23],ICA.BF$R2_H[ICA.BF$R2_H==min(ICA.BF$R2_H)], 
       col=1,pch=15,cex=1.2)
lines(ICA.BF$Lower.Dig.Corrs.Sigma[j,23],lw1$fitted[j],lty=2,lwd=2)

plot(ICA.BF$Lower.Dig.Corrs.Sigma[,25],ICA.BF$R2_H,col='#999999',ylab='ICA',
     xlab=expression(rho),main='(b)',cex.lab=1.6,cex.axis=1.2,cex.main=1.6)
lw1 = loess(ICA.BF$R2_H ~ ICA.BF$Lower.Dig.Corrs.Sigma[,25])
j <- order(ICA.BF$Lower.Dig.Corrs.Sigma[,25])
points(ICA.BF$Lower.Dig.Corrs.Sigma[ICA.BF$R2_H==min(ICA.BF$R2_H),25],ICA.BF$R2_H[ICA.BF$R2_H==min(ICA.BF$R2_H)], 
       col=1,pch=15,cex=1.2)
lines(ICA.BF$Lower.Dig.Corrs.Sigma[j,25],lw1$fitted[j],lty=2,lwd=2)

plot(ICA.BF$Lower.Dig.Corrs.Sigma[,28],ICA.BF$R2_H,col='#999999',ylab='ICA',
     xlab=expression(rho),main='(c)',cex.lab=1.6,cex.axis=1.2,cex.main=1.6)
lw1 = loess(ICA.BF$R2_H ~ ICA.BF$Lower.Dig.Corrs.Sigma[,28])
j <- order(ICA.BF$Lower.Dig.Corrs.Sigma[,28])
points(ICA.BF$Lower.Dig.Corrs.Sigma[ICA.BF$R2_H==min(ICA.BF$R2_H),28],ICA.BF$R2_H[ICA.BF$R2_H==min(ICA.BF$R2_H)], 
       col=1,pch=15,cex=1.2)
lines(ICA.BF$Lower.Dig.Corrs.Sigma[j,28],lw1$fitted[j],lty=2,lwd=2)


## using the PC algorithm
ICA.PC.oneSurr = ICA.ContCont.MultS.PC(M = 10000, 100, Sigma.4[1:4,1:4]) 
range(ICA.PC.oneSurr$R2_H)
ICA.PC = ICA.ContCont.MultS.PC(M = 10000, 100, Sigma.4) 
range(ICA.PC$R2_H)
## using the MPC algorithm
# determing probabilities
p=3 # number of surrogates
prob = rep(1/p,p) 
prob = choose(p,1:p)/sum(choose(p,1:p)) # probabilities for each combination of r surrogates
ICA.MPC = ICA.ContCont.MultS.MPC(70000,100,Sigma.4,Seed=123,prob = prob,Save.Corr=T,Show.Progress = T)
range(ICA.MPC$R2_H)

### checking correlations
CorrLabel = c(1,3,5,7,8,10,12,14,16,18,19,21,23,25,26,28)
dataCorr = mapply(function(x){
  correlation = c(ICA.BF$Lower.Dig.Corrs.Sigma[,CorrLabel[x]],
                  ICA.PC$Lower.Dig.Corrs.All[,CorrLabel[x]])
  corr = as.factor(x)
  algorithm = rep(c('BF','PC'),c(10000,10000))
  data.frame(algorithm,corr,correlation)
},x=1:length(CorrLabel),SIMPLIFY = F)

dataCorr = do.call(rbind,dataCorr)


CorrLabel2 = c('T[0]~T[1]','T[0]~S[11]','T[0]~S[21]',
               'T[0]~S[31]','T[1]~S[10]','T[1]~S[20]',
               'T[1]~S[30]','S[10]~S[11]','S[10]~S[21]',
               'S[10]~S[31]','S[11]~S[20]','S[11]~S[30]',
               'S[20]~S[21]','S[20]~S[31]','S[21]~S[30]',
               'S[30]~S[31]')

levels(dataCorr$corr) = CorrLabel2
library(ggplot2)

# Figure S4 Density plots of the unidentifiable correlations (RS and PC algorithm)
ggplot(data=dataCorr, aes(x=correlation, colour =algorithm)) +
  scale_colour_grey(start = 0, end = 0.6)  +
  geom_density(alpha = 1, show.legend = FALSE) + 
  facet_wrap(~ corr,ncol = 4,labeller = label_parsed) +
  stat_function(fun = extraDistr::dnsbeta, 
                args = list(shape1 = 4, shape2 = 4,min=-1,max=1), 
                lty = 2, 
                col = 'red'
  )


##### MPC algorithm 
dataCorr.MPC = mapply(function(x){
  correlation = ICA.MPC$Lower.Dig.Corrs.All[,CorrLabel[x]]
  corr = as.factor(x)
  method = 'MPC'
  data.frame(method,corr,correlation)
},x=1:length(CorrLabel),SIMPLIFY = F)
dataCorr.MPC = do.call(rbind,dataCorr.MPC)

levels(dataCorr.MPC$corr) = CorrLabel2

# Figure S5 Density plots of the unidentifiable correlations using the MPC algorithm
ggplot(data=dataCorr.MPC, aes(x=correlation)) +
  geom_density(alpha = 0.7) + facet_wrap(~ corr,labeller = label_parsed)

### checking ICA densities
ICA.values = c(ICA.BF$R2_H,ICA.PC$R2_H,ICA.MPC$R2_H)
ICA.method = c(rep(c('RS','PC'),each=10000),rep('MPC',70000))

ICA.density4 = data.frame(ICA=ICA.values,algorithm=ICA.method)

Range.4 = by(ICA.density4$ICA,ICA.density4$algorithm,range)
MPC.min = unlist(Range.4)[c(1:3)*2-1]
MPC.max = unlist(Range.4)[c(1:3)*2]

MPC.range.4 = data.frame(min=MPC.min,max=MPC.max,algorithm=c('MPC','PC','RS'))

dev.copy2pdf(file="ICADensties4.pdf",out.type="cairo", width=7, height=4)
# Figure 2 Densities of the ICA computed using the RA, PC and MPC algorithm.
ggplot(data=ICA.density4, aes(x=ICA, fill=algorithm)) +
  geom_density(alpha = 0.5) +
  scale_fill_grey(start = 0, end = 0.95) +
  geom_vline(data=MPC.range.4, aes(xintercept=min),
             linetype=c(1,2,3))  + theme_classic() +
  theme(panel.border =element_rect(fill=NA,color = "black"),
        axis.line = element_line(colour = "black"))

dev.off()


##### Simulated case (Section 7)
N = 200
noise.sur = 13 # (max. number of noise-surrogates)

R.noise = matrix(0.5,2*noise.sur,2*noise.sur)
diag(R.noise) = 1
Dn = diag(rep(c(4,1),noise.sur))
Sigma.noise = Dn%*%R.noise%*%Dn

R.test = invvech(c(1,-0.2, 0.9, -0.5, -0.6, 0.5,1, -0.5, 0.7, 0.4, -0.5, 1, -0.7, 
                   -0.5 ,0.7, 1, 0.6, -0.4,1, -0.2, 1))


Dtest = diag(c(5,2,4,1,4,1))
Sigma.test = Dtest%*%R.test%*%Dtest

eigen(Sigma.test)$values
Surrogate::ICA.ContCont.MultS(1,100,Sigma.test)$R2_H
Surrogate::ICA.ContCont.MultS(1,100,Sigma.test[1:4,1:4])$R2_H

Sigma.test.2 = as.matrix(bdiag(Sigma.test,Sigma.noise)) # complete variance matrix
pp.max = ncol(Sigma.test.2)/2 - 1
set.seed(1234)
Data.complete = mvrnorm(N,rep(0,pp.max*2+2),Sigma.test.2) # simulating the dataset
ind = (1:(pp.max+1))*2 
Data.complete[1:(N/2),ind] =NA 
Data.complete[(N/2+1):N,-ind] =NA
Sigma.sample = var(Data.complete,use='pairwise.complete.obs') # incomplete variance matrix

### Simulation of the ICA for different number of surrogates
pp.all = c(1,2,3,5,10,15)
M=10000
Output.MPC = mapply(function(x){
  pp = pp.all[x]
  if(pp > 5){
    comb.p = choose(pp,1:pp)
    prob = choose(pp,1:pp)/sum(comb.p)*c(1,1,rep(0,length(comb.p)-3),1)
    prob = prob/sum(prob)
  }else{
    prob = choose(pp,1:pp)/sum(choose(pp,1:pp))
    prob = prob/sum(prob)
  }
  M2 = round((M/(prob/choose(pp,1:pp)))[1])
  
  ind = (pp+1)*2
  tic()
  test.1 = ICA.ContCont.MultS.MPC(M2,N,Sigma.sample[1:ind,1:ind],prob=prob,Seed = 1234,
                                  Show.Progress = T,nCores = 8)
  exec =  toc(quiet = T)
  comp.time1 = exec$toc - exec$tic
  tic()
  prob = c(rep(0,length(prob)-1),1)
  test.2 = ICA.ContCont.MultS.MPC(M2,N,Sigma.sample[1:ind,1:ind],prob=prob,Seed = 1234,
                                  Show.Progress = T,nCores = 8)
  exec =  toc(quiet = T)
  comp.time2 = exec$toc - exec$tic
  list(R2H.MPC = test.1$R2_H,R2H.PC = test.2$R2_H)
},x=1:length(pp.all))

## RS algorithm
output.RA1 = Surrogate::ICA.ContCont.MultS(10000, 200, Sigma.sample[1:4,1:4])
output.RA2 = Surrogate::ICA.ContCont.MultS(30000, 200, Sigma.sample[1:6,1:6])
output.RA3 = ICA.ContCont.MultS(70000, 200, Sigma.sample[1:8,1:8])

Output.RA = rbind(quantile(output.RA1$R2_H,c(0,0.5,1)),
                  quantile(output.RA2$R2_H,c(0,0.5,1)),
                  quantile(output.RA3$R2_H,c(0,0.5,1)))

## Figure 3 Range of the ICA computed using the RA algoritm (gray dashed line), PC algorithm
## (gray solid line), and MPC algorithm (black solid line)
plot(NULL,xlim=c(0,15),ylim=c(0,1),ylab='ICA',xlab='Number of surrogates',
     xaxt='n')
for(i in 1:6){
  lines(c(pp.all[i],pp.all[i]),range(Output.MPC[[i*2]]),lwd=2,col='lightgrey')
  lines(c(pp.all[i],pp.all[i])+0.3, range(Output.MPC[[i*2-1]]),lwd=2)
  points(pp.all[i]+0.3,median(Output.MPC[[i*2-1]]),pch=19)
  points(pp.all[i],median(Output.MPC[[i*2]]),pch=19,col='grey')
  if(i <=3){
    lines(c(pp.all[i],pp.all[i])-0.3,Output.RA[i,c(1,3)],lwd=2,col='grey',lty=3)
    points(pp.all[i]-0.3,Output.RA[i,2],pch=19,col='grey')
  }
}
axis(1,c(1,2,3,5,10,15))


### data preparation for Figure S6
ICA = c(unlist(Output.MPC[c(3:12)]),output.RA2$R2_H,output.RA3$R2_H)
algorithm= c(rep(c('MPC','PC'),each=30000),
             rep(c('MPC','PC'),each=70000),
             rep(c('MPC','PC'),each=310000),
             rep(c('MPC','PC'),each=560000),
             rep(c('MPC','PC'),each=1210000),
             rep('RS',100000))
p = c(rep(c(2,3,5,10,15),c(30000,70000,310000,560000,1210000)*2),
      rep(c(2,3),c(30000,70000)))

MPC.plot = data.frame(ICA=ICA,algorithm=algorithm,p=as.factor(p))

Range = by(MPC.plot$ICA,MPC.plot[,2:3],range)
MPC.min = unlist(Range)[c(1:12)*2-1]
MPC.max = unlist(Range)[c(1:12)*2]

MPC.range = data.frame(min=MPC.min,max=MPC.max,algorithm=c(rep(c('MPC','PC','RS'),2),
                                                           rep(c('MPC','PC'),3)),
                       p=as.factor(c(rep(c(2,3),each=3),rep(c(5,10,15),each=2))))


dev.copy2pdf(file="plot.pdf",out.type="cairo", width=10, height=5)
# Figure S5 Density plots of the ICA using the RA, PC, and MPC algorithms.
levels(MPC.plot$p) = c('p==~2','p==~3','p==~5','p==~10','p==~15')
levels(MPC.range$p) = c('p==~2','p==~3','p==~5','p==~10','p==~15')
ggplot(data=MPC.plot, aes(x=ICA,fill=algorithm)) +
  geom_density(alpha = 0.5) + facet_wrap(~ p,labeller = label_parsed) +
  geom_vline(data=MPC.range, aes(xintercept=max, color=algorithm),
             linetype="dashed") +
  geom_vline(data=MPC.range, aes(xintercept=min, color=algorithm),
             linetype="dashed") 

dev.off()

### Figure S6 density plot correlations

geom_hline(data = dummy2, aes(yintercept = Z))

### computation time
pp.all = c(rep(5,30),rep(10,30),rep(15,30))
Output2 = pbmapply(function(x){
  pp = pp.all[x]
  ind = (pp+1)*2
  tic()
  test.1 = ICA.ContCont.MultS.PC.alt(1000,N,Sigma.sample[1:ind,1:ind],prob=NULL,Seed = 123)
  exec =  toc(quiet = T)
  comp.time1 = exec$toc - exec$tic
  tic()
  prob = c(rep(0,pp-1),1)
  test.2 = ICA.ContCont.MultS.PC.alt(1000,N,Sigma.sample[1:ind,1:ind],prob=prob,Seed = 123)
  exec =  toc(quiet = T)
  comp.time2 = exec$toc - exec$tic
  c(quantile(test.1$R2_H,c(0,0.5,1)),
    quantile(test.2$R2_H,c(0,0.5,1)),comp.time1,comp.time2)
},x=1:length(pp.all))

## computation for the MPC and PC algorithm
c(mean(Output2[7,pp.all==5]),mean(Output2[7,pp.all==10]),mean(Output2[7,pp.all==15]))
c(mean(Output2[8,pp.all==5]),mean(Output2[8,pp.all==10]),mean(Output2[8,pp.all==15]))

### simulating several datasets (Supplemantary Materials)
### Section C
noise.sur = 2
Sigma.test = invvech(c(1,0,0.95,0,1,0,0.95,1,0,1))
R.noise = invvech(c(1,0,0.8,0,1,0,0.8,1,0,1))
Sigma.test.2 = as.matrix(bdiag(Sigma.test,R.noise)) # complete variance matrix

pp=3
prob = choose(pp,1:pp)/sum(choose(pp,1:pp))
prob = prob/sum(prob)

M2 = round((10000/(prob/choose(pp,1:pp)))[1])
M = 200
Sim.N50 = pbmapply(function(x){
  set.seed(122+x)
  N=50
  Data.complete = mvrnorm(N,rep(0,3*2+2),Sigma.test.2) # simulating the dataset
  ind = (1:(3+1))*2 
  Data.complete[1:(N/2),ind] =NA 
  Data.complete[(N/2+1):N,-ind] =NA
  Sigma.sample = var(Data.complete,use='pairwise.complete.obs') # incomplete variance matrix
  #opb <- pboptions(type="none")
  #on.exit(pboptions(opb))
  test.1 = ICA.ContCont.MultS.MPC(M2,N,Sigma.sample,prob=prob,Seed = 123,Show.Progress = F)
  test.2 = ICA.ContCont.MultS.MPC(M2,N,Sigma.sample,prob=c(0,0,1),Seed = 123,Show.Progress = F)
  c(quantile(test.1$R2_H,c(0,0.5,1)),quantile(test.2$R2_H,c(0,0.5,1)))
},x=1:M)

Sim.N100 = pbmapply(function(x){
  set.seed(122+x)
  N=100
  Data.complete = mvrnorm(N,rep(0,3*2+2),Sigma.test.2) # simulating the dataset
  ind = (1:(3+1))*2 
  Data.complete[1:(N/2),ind] =NA 
  Data.complete[(N/2+1):N,-ind] =NA
  Sigma.sample = var(Data.complete,use='pairwise.complete.obs') # incomplete variance matrix
  #opb <- pboptions(type="none")
  #on.exit(pboptions(opb))
  test.1 = ICA.ContCont.MultS.MPC(M2,N,Sigma.sample,prob=prob,Seed = 123,Show.Progress = F)
  test.2 = ICA.ContCont.MultS.MPC(M2,N,Sigma.sample,prob=c(0,0,1),Seed = 123,Show.Progress = F)
  c(quantile(test.1$R2_H,c(0,0.5,1)),quantile(test.2$R2_H,c(0,0.5,1)))
},x=1:M)

Sim.N200 = pbmapply(function(x){
  set.seed(122+x)
  N=200
  Data.complete = mvrnorm(N,rep(0,3*2+2),Sigma.test.2) # simulating the dataset
  ind = (1:(3+1))*2 
  Data.complete[1:(N/2),ind] =NA 
  Data.complete[(N/2+1):N,-ind] =NA
  Sigma.sample = var(Data.complete,use='pairwise.complete.obs') # incomplete variance matrix
  #opb <- pboptions(type="none")
  #on.exit(pboptions(opb))
  test.1 = ICA.ContCont.MultS.MPC(M2,N,Sigma.sample,prob=prob,Seed = 123,Show.Progress = F)
  test.2 = ICA.ContCont.MultS.MPC(M2,N,Sigma.sample,prob=c(0,0,1),Seed = 123,Show.Progress = F)
  c(quantile(test.1$R2_H,c(0,0.5,1)),quantile(test.2$R2_H,c(0,0.5,1)))
},x=1:M)

Sim.N500 = pbmapply(function(x){
  set.seed(122+x)
  N=500
  Data.complete = mvrnorm(N,rep(0,3*2+2),Sigma.test.2) # simulating the dataset
  ind = (1:(3+1))*2 
  Data.complete[1:(N/2),ind] =NA 
  Data.complete[(N/2+1):N,-ind] =NA
  Sigma.sample = var(Data.complete,use='pairwise.complete.obs') # incomplete variance matrix
  #opb <- pboptions(type="none")
  #on.exit(pboptions(opb))
  test.1 = ICA.ContCont.MultS.MPC(M2,N,Sigma.sample,prob=prob,Seed = 123,Show.Progress = F)
  test.2 = ICA.ContCont.MultS.MPC(M2,N,Sigma.sample,prob=c(0,0,1),Seed = 123,Show.Progress = F)
  c(quantile(test.1$R2_H,c(0,0.5,1)),quantile(test.2$R2_H,c(0,0.5,1)))
},x=1:M)

Sim.N1000 = pbmapply(function(x){
  set.seed(122+x)
  N=1000
  Data.complete = mvrnorm(N,rep(0,3*2+2),Sigma.test.2) # simulating the dataset
  ind = (1:(3+1))*2 
  Data.complete[1:(N/2),ind] =NA 
  Data.complete[(N/2+1):N,-ind] =NA
  Sigma.sample = var(Data.complete,use='pairwise.complete.obs') # incomplete variance matrix
  #opb <- pboptions(type="none")
  #on.exit(pboptions(opb))
  test.1 = ICA.ContCont.MultS.MPC(M2,N,Sigma.sample,prob=prob,Seed = 123,Show.Progress = F)
  test.2 = ICA.ContCont.MultS.MPC(M2,N,Sigma.sample,prob=c(0,0,1),Seed = 123,Show.Progress = F)
  c(quantile(test.1$R2_H,c(0,0.5,1)),quantile(test.2$R2_H,c(0,0.5,1)))
},x=1:M)


### 
load('SimulationRangeSampleSize.RData')

MIN = c(Sim.N50[1,],Sim.N100[1,],Sim.N200[1,],Sim.N500[1,],Sim.N1000[1,],
        Sim.N50[4,],Sim.N100[4,],Sim.N200[4,],Sim.N500[4,],Sim.N1000[4,])
MEDIAN = c(Sim.N50[2,],Sim.N100[2,],Sim.N200[2,],Sim.N500[2,],Sim.N1000[2,],
           Sim.N50[5,],Sim.N100[5,],Sim.N200[5,],Sim.N500[5,],Sim.N1000[5,])
MAX = c(Sim.N50[3,],Sim.N100[3,],Sim.N200[3,],Sim.N500[3,],Sim.N1000[3,],
        Sim.N50[6,],Sim.N100[6,],Sim.N200[6,],Sim.N500[6,],Sim.N1000[6,])
SS = as.factor(rep(rep(c(50,100,200,500,1000),each=200),2))
METHOD = rep(c('mod.PC','PC'),each=5*200)

MIN.DATA = data.frame(SS,METHOD,MIN)
MED.DATA = data.frame(SS,METHOD,MEDIAN)
MAX.DATA = data.frame(SS,METHOD,MAX)

library(ggplot2)
library(ggpubr)

MIN.FIG = ggplot(data = MIN.DATA, aes(x=SS, y=MIN,fill=METHOD,group=interaction(SS,METHOD)))  + geom_boxplot() +
  theme(legend.position= "none") + scale_fill_grey(start=0.4,end=0.85) + 
  labs(y="Minimum ICA", x = "Sample size")

MED.FIG = ggplot(data = MED.DATA, aes(x=SS, y=MEDIAN,fill=METHOD,group=interaction(SS,METHOD)))  + geom_boxplot() +
  theme(legend.position= "none") + scale_fill_grey(start=0.4,end=0.85) + 
  labs(y="Median ICA", x = "Sample size")

MAX.FIG = ggplot(data = MAX.DATA, aes(x=SS, y=MAX,fill=METHOD,group=interaction(SS,METHOD)))  + geom_boxplot() +
  theme(legend.position= "none") + scale_fill_grey(start=0.4,end=0.85) + 
  labs(y="Maximum ICA", x = "Sample size")

#### Figure S6 Boxplots of the minimum ICA
MIN.FIG
#### Figure S7 Boxplots of the median ICA
MED.FIG
#### Figure S8 Boxplots of the maximum ICA
MAX.FIG

### extra simulation with strong correlations - Section D
ar1_cor <- function(n, rho) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                    (1:n - 1))
  rho^exponent
}
# setting
n.sur = 15
Corr0 = ar1_cor(1+n.sur,0.9)
Corr1 = ar1_cor(1+n.sur,0.8)
R = as.matrix(bdiag(Corr0,Corr1))
Dn = diag(c(5,rep(4,15),2,rep(1,15)))
Sigma = Dn%*%R%*%Dn
Sigma[R==0] = NA
d = ncol(Sigma)
IND = ks::vec(matrix(1:d, ncol = 2), byrow = TRUE)
Sigma.test = as.matrix(Sigma[IND,IND])
pp.all = c(1,2,3,5,10,15)
M=20000

Output.MPC.2 = mapply(function(x){
  pp = pp.all[x]
  if(pp > 5){
    comb.p = choose(pp,1:pp)
    prob = choose(pp,1:pp)/sum(comb.p)*c(1,1,rep(0,length(comb.p)-3),1)
    prob = prob/sum(prob)
  }else{
    prob = choose(pp,1:pp)/sum(choose(pp,1:pp))
    prob = prob/sum(prob)
  }
  M2 = round((M/(prob/choose(pp,1:pp)))[1])
  
  ind = (pp+1)*2
  tic()
  test.1 = ICA.ContCont.MultS.PC.alt(M2,N,Sigma.test[1:ind,1:ind],prob=NULL,Seed = 123,
                                     Show.Progress = T,nCores = 8)
  exec =  toc(quiet = T)
  comp.time1 = exec$toc - exec$tic
  tic()
  prob = c(rep(0,length(prob)-1),1)
  test.2 = ICA.ContCont.MultS.PC.alt(M2,N,Sigma.test[1:ind,1:ind],prob=prob,Seed = 123,
                                     Show.Progress = T,nCores = 8)
  exec =  toc(quiet = T)
  comp.time2 = exec$toc - exec$tic
  list(R2H.MPC = test.1$R2_H,R2H.PC = test.2$R2_H)
},x=1:length(pp.all))

## Figure S9
plot(NULL,xlim=c(0,15),ylim=c(0,1),ylab='ICA',xlab='Number of surrogates',
     xaxt='n')
for(i in 1:6){
  lines(c(pp.all[i],pp.all[i]),range(Output.MPC.2[[i*2]]),lwd=2,col='lightgrey')
  lines(c(pp.all[i],pp.all[i])+0.3, range(Output.MPC.2[[i*2-1]]),lwd=2)
  points(pp.all[i]+0.3,median(Output.MPC.2[[i*2-1]]),pch=19)
  points(pp.all[i],median(Output.MPC.2[[i*2]]),pch=19,col='grey')
}

####### Data analysis (Section 7)
Data.micro = read.csv('Microbiome_data.csv')

True  = Data.micro$True
Surrogates = Data.micro[,-c(1:3)]
treat = Data.micro$Treat

p = ncol(Data.micro) - 3
Sigma.control = cov(Data.micro[treat==-1,-c(1:2)])
Sigma.treat = cov(Data.micro[treat==1,-c(1:2)])
Sigma.NA = matrix(NA,p+1,p+1)
Sigma = rbind(cbind(Sigma.control,Sigma.NA),cbind(Sigma.NA,Sigma.treat))
colnames(Sigma) = rownames(Sigma)=c('T0',paste('S',1:p,0,sep=''),'T1',paste('S',1:p,1,sep=''))
d = (p+1)*2
IND =  ks::vec(matrix(1:d,ncol=2),byrow = TRUE) 
Sigma = Sigma[IND,IND]
N = nrow(Data.micro)

### Surrogate selection
#### step 1 (evaluating one surrogate)
ICA.1 = mapply(function(x){
  indicator = c(1:2,1 + (2*x):(2*x+1))
  Sigma.test = Sigma[indicator,indicator]
  ICA.test = ICA.ContCont.MultS.MPC(10000, N, Sigma.test, Seed = 123,prob=NULL, Show.Progress = FALSE)
  ICA.summary = quantile(ICA.test$R2_H,c(0,0.5,1))
  outcome = c(x,ICA.summary,table(ICA.test$surr.eval.r))
  return(outcome)
},x=1:p)
ICA.1 = t(ICA.1)
head(ICA.1[order(ICA.1[,3],decreasing = T),]) #### select sorrogate S_44
### step 2 (evaluating two surrogates) fixing S_44
surr.eval = 1:p
surr.eval = surr.eval[-44]

ICA.2 = mapply(function(x){
  sur = surr.eval[x]
  label.sur = c(44,sur)
  indicator = c(1:2,89:90,1 + (2*sur):(2*sur+1))
  
  Sigma.test = Sigma[indicator,indicator]
  ICA.test = ICA.ContCont.MultS.MPC(30000, N, Sigma.test, Seed = 123,prob=NULL, Show.Progress = FALSE)
  ICA.summary = quantile(ICA.test$R2_H,c(0,0.5,1))
  outcome = c(label.sur,ICA.summary,table(ICA.test$surr.eval.r))
  return(outcome)
},x=1:length(surr.eval))
ICA.2 = t(ICA.2)
head(ICA.2[order(ICA.2[,4],decreasing = F),]) #### select sorrogate S_44

### step 3 (evaluating three surrogates) fixing S_44 S_17
surr.eval = 1:p
surr.eval = surr.eval[-c(17,44)]

pp = 3
prob = choose(pp,1:pp)/sum(choose(pp,1:pp))
(10000/(prob/choose(pp,1:pp)))[1]

ICA.3 = mapply(function(x){
  sur = surr.eval[x]
  label.sur = c(44,17,sur)
  indicator = c(1:2,89:90,35:36,1 + (2*sur):(2*sur+1))
  
  Sigma.test = Sigma[indicator,indicator]
  ICA.test = ICA.ContCont.MultS.MPC(70000, N, Sigma.test, Seed = 123,prob=prob)
  ICA.summary = quantile(ICA.test$R2_H,c(0,0.5,1),na.rm = T)
  ICA.NA = mean(is.na(ICA.test$R2_H))
  outcome = c(label.sur,ICA.summary,table(ICA.test$surr.eval.r))
  return(outcome)
},x=1:length(surr.eval))
ICA.3 = t(ICA.3)
head(ICA.3[order(ICA.3[,5],decreasing = T),]) #### select sorrogate S_44 S_17 S_37

### step 4 (evaluating four surrogates) fixing S_44 S_17 S_37
surr.eval = 1:p
surr.eval = surr.eval[-c(17,37,44)]

pp = 4
prob = choose(pp,1:pp)/sum(choose(pp,1:pp))
M = (10000/(prob/choose(pp,1:pp)))[1]

ICA.4 = mapply(function(x){
  sur = surr.eval[x]
  label.sur = c(44,17,37,sur)
  indicator = c(1:2,89:90,35:36,75:76,1 + (2*sur):(2*sur+1))
  
  Sigma.test = Sigma[indicator,indicator]
  ICA.test = ICA.ContCont.MultS.MPC(M, N, Sigma.test, Seed = 123,prob=prob)
  ICA.summary = quantile(ICA.test$R2_H,c(0,0.5,1),na.rm = T)
  ICA.NA = mean(is.na(ICA.test$R2_H))
  outcome = c(label.sur,ICA.summary,colMeans(ICA.test$surr.eval.r))
  return(outcome)
},x=1:length(surr.eval))
ICA.4 = t(ICA.4)
head(ICA.4[order(ICA.4[,6],decreasing = T),]) #### select sorrogate S_44 S_17 S_37 S_40

### step 5 (evaluating four surrogates) fixing S_44 S_17 S_37 S_40
surr.eval = 1:p
surr.eval = surr.eval[-c(17,37,40,44)]

pp = 5
prob = choose(pp,1:pp)/sum(choose(pp,1:pp))
M = (1000/(prob/choose(pp,1:pp)))[1]
M =50000
M*prob/choose(pp,1:pp)

ICA.5 = mapply(function(x){
  sur = surr.eval[x]
  label.sur = c(44,17,37,40,sur)
  indicator = c(1:2,89:90,35:36,75:76,81:82,1 + (2*sur):(2*sur+1))
  
  Sigma.test = Sigma[indicator,indicator]
  ICA.test = ICA.ContCont.MultS.MPC(M, N, Sigma.test, Seed = 123,prob=prob)
  ICA.summary = quantile(ICA.test$R2_H,c(0,0.5,1),na.rm = T)
  ICA.NA = mean(is.na(ICA.test$R2_H))
  outcome = c(label.sur,ICA.summary,colMeans(ICA.test$surr.eval.r))
  return(outcome)
},x=1:length(surr.eval))
ICA.5 = t(ICA.5)
head(ICA.5[order(ICA.5[,6],decreasing = T),]) #### select sorrogate S_44 S_17 S_37 S_40 S_30

##### Comparing algorithms up to 5 surrogates
## using the GRS algorithm
dim(Sigma)/2
indicator = c(1:2,89:90,35:36,75:76,81:82,61:62)
Sigma.5 = Sigma[indicator,indicator]

##### GRS algorithm
# one surrogate
GRS.1 = ICA.ContCont.MultS(10000,16,Sigma.5[1:4,1:4])
range(GRS.1$R2_H)

# two surrogate
r = 30
M = 30000/r
Seed = 122 + 1:r
cl <- makeCluster(getOption("cl.cores", 3))
clusterExport(cl=cl, varlist=c('M','r','ICA.ContCont.MultS','Sigma.5','Seed'), envir=environment())
clusterEvalQ(cl=cl,library('tictoc'))
GRS.2 = pblapply(X=1:r,function(X) {
  tic()
  ICA = ICA.ContCont.MultS(M,16,Sigma.5[1:6,1:6],Seed=Seed[X])
  time.f = toc()
  time = time.f$toc - time.f$tic
  fit = c(ICA$R2_H,time)
  return(fit)
}, cl=cl)
stopCluster(cl)

GRS2 = do.call('rbind',GRS.2)
quantile(GRS2[,-1001],c(0,0.5,1))
# three surrogate
r = 70
M = 70000/r
Seed = 122 + 1:r
cl <- makeCluster(getOption("cl.cores", 3))
clusterExport(cl=cl, varlist=c('M','ICA.ContCont.MultS','Sigma.5','Seed'), envir=environment())
clusterEvalQ(cl=cl,library('tictoc'))
GRS.3 = pblapply(X=1:70,function(X) {
  tic()
  ICA = ICA.ContCont.MultS(M,16,Sigma.5[1:8,1:8],Seed=Seed[X])
  time.f = toc()
  time = time.f$toc - time.f$tic
  fit = c(ICA$R2_H,time)
  return(fit)
}, cl=cl)
stopCluster(cl)

GRS3 = do.call('rbind',GRS.3)
quantile(GRS3[,-1001],c(0,0.5,1))

### using the PC algorithm

ICA.1.PC = ICA.ContCont.MultS.MPC(10000,100,Sigma.5[1:4,1:4],Show.Progress = T,Save.Corr = F)
range(ICA.1.PC$R2_H,na.rm = T)
ICA.2.PC = ICA.ContCont.MultS.MPC(30000,100,Sigma.5[1:6,1:6],prob=c(0,1),Show.Progress = T,Save.Corr = F)
range(ICA.2.PC$R2_H,na.rm = T)
ICA.3.PC = ICA.ContCont.MultS.MPC(70000,100,Sigma.5[1:8,1:8],prob=c(0,0,1),Show.Progress = T,Save.Corr = F)
range(ICA.3.PC$R2_H,na.rm = T)
ICA.4.PC = ICA.ContCont.MultS.MPC(150000,100,Sigma.5[1:10,1:10],prob=c(0,0,0,1),Show.Progress = T,Save.Corr = F)
range(ICA.4.PC$R2_H,na.rm = T)
ICA.5.PC = ICA.ContCont.MultS.MPC(150000,100,Sigma.5,prob=c(0,0,0,0,1),Show.Progress = T,Save.Corr = F)
range(ICA.5.PC$R2_H,na.rm = T)

