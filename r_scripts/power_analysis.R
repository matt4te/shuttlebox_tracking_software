

### Conduct a power analysis to determine sample size / preference requirements for significance ###


## Set up shop -----------------


# set working directory
dir <- "C:/Users/Matt/Desktop/Current_Lab_Work/ShadowBox"
setwd(dir)


### JUST GO THE "Shuttle_Group_Analysis.R" SCRIPT AND LOAD THE DATA FILE -- 
      # ALSO MAKE SURE TO CREATE "NORM DIFF" VARIABLE THERE

# load required packages
library("MASS")
library("plyr")
library("ggplot2")
library("reshape2")


## Basic routine for power analysis ---------------


# get the parameters of our sample distribution

f1 <- fitdistr(data$normDiff,"normal")

mean1 <- as.numeric(f1$estimate[1])
sd1 <- as.numeric(f1$estimate[2])
n1 <- nrow(data)

# simulate data with these parameters

norm1 <- rnorm(mean=mean1, sd=sd1, n=n1)

# can we determine significance?

t1 <- t.test(norm1, mu=0, alternative = c("greater"))
t1$p.value


# now we repeat this 500 times and extract a p value each time

nsim=500
pval = numeric(nsim)

for (i in 1:nsim) {
  normSim = rnorm(n=n1, mean=mean1, sd=sd1)
  pval[i] = t.test(normSim, mu=0, alternative = c("greater"))$p.value
}

power = sum(pval < .01)/nsim
power



## Investigate species-level preference for major taxa ------------------

# pull out partitus data only
dataP <- subset(data, fish == "Stegastes partitus")

# test for significance
tP <- t.test(dataP$normDiff, mu=0, alternative = c("greater"))
tP$p.value
# 0.005332692

# pull out diancaeus data only
dataD <- subset(data, fish == "Stegastes spp" | fish == "Stegastes diancaeus")

# test for significance
tD <- t.test(dataD$normDiff, mu=0, alternative = c("greater"))
tD$p.value
# 0.1623




## What if we increase the sample size of both? ----------------

# for partititus 

# first, extract distribution paramters
fP <- fitdistr(dataP$normDiff,"normal")
meanP <- as.numeric(fP$estimate[1])
sdP <- as.numeric(fP$estimate[2])
nP <- nrow(dataP)

# cycle through range of N
Nvec = seq(5,50, by =5)
power.N = numeric(length(Nvec)) #pre-allocate this variable to fill later

# 500 times each
nsim=500
pval = numeric(nsim)

# implement loop
for (j in 1:length(Nvec)) {
  N = Nvec[j]
  for (i in 1:nsim) {
    normSimP = rnorm(n=N, mean=meanP, sd=sdP)
    pval[i] = t.test(normSimP, mu=0, alternative = c("greater"))$p.value
  }
  power.N[j] = sum(pval < 0.01)/nsim
}

power.N

## What if we vary the mean value and N? -------------

mVec <- seq(.1,.4,by=.05)

power.Nm <- matrix(nrow=length(Nvec),ncol=length(mVec))

# implement nested loop
for (j in 1:length(Nvec)) {
  Ns = Nvec[j]
  for (k in 1:length(mVec)) {
    meanS = mVec[k]
    for (i in 1:nsim) {
      normSimS <- rnorm(n=Ns, mean=meanS, sd=sdP)
      pval[i] = t.test(normSimS, mu=0, alternative = c("greater"))$p.value
    }
    power.Nm[j,k] = sum(pval < 0.05)/nsim
  }
}
rownames(power.Nm) <- Nvec
colnames(power.Nm) <- mVec

power.Nm


# plot the power contour

meltP <- melt(power.Nm)
names(meltP) <- c("x", "y", "z")


Pp <- ggplot(meltP, aes(x, y, z=z))  +
  theme_classic()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  geom_tile(aes(fill = z)) + stat_contour(bins=5, size=1, color="blacK") +
  theme(plot.title = element_text(lineheight=.8, face="bold")) +
  xlab("Sample Size (N)") + ylab("Mean Count Difference") +
  geom_segment(aes(x=3,xend=nrow(dataP),y=meanP,yend=meanP, color="red")) +
  geom_segment(aes(x=nrow(dataP),xend=nrow(dataP),y=0.075,yend=meanP, color="red")) +
  scale_colour_discrete(name = "S. partitus", labels=c("Real Data")) +
  guides(fill=guide_legend(title="Statistical Power"))

Pp

# what is actual data p value
tP <- t.test(dataP$normDiff, mu=0, alternative = c("greater"))
tP$p.value


# for diancaeus 

# first, extract distribution paramters
fD <- fitdistr(dataD$normDiff,"normal")
meanD <- as.numeric(fD$estimate[1])
sdD <- as.numeric(fD$estimate[2])
nD <- nrow(dataD)

# cycle through range of N
Nvec = seq(5,50, by =5)
power.Nd = numeric(length(Nvec)) #pre-allocate this variable to fill later

# 500 times each
nsim=500
pval = numeric(nsim)

# implement loop
for (j in 1:length(Nvec)) {
  N = Nvec[j]
  for (i in 1:nsim) {
    normSimD = rnorm(n=N, mean=meanD, sd=sdD)
    pval[i] = t.test(normSimD, mu=0, alternative = c("greater"))$p.value
  }
  power.Nd[j] = sum(pval < 0.01)/nsim
}

power.Nd

# lets also vary the mean slightly

mVec <- seq(.1,.4,by=.05)

power.Nmd <- matrix(nrow=length(Nvec),ncol=length(mVec))

# implement nested loop
for (j in 1:length(Nvec)) {
  Nsd = Nvec[j]
  for (k in 1:length(mVec)) {
    meanSd = mVec[k]
    for (i in 1:nsim) {
      normSimSd <- rnorm(n=Nsd, mean=meanSd, sd=sdD)
      pval[i] = t.test(normSimSd, mu=0, alternative = c("greater"))$p.value
    }
    power.Nmd[j,k] = sum(pval < 0.01)/nsim
  }
}
rownames(power.Nmd) <- Nvec
colnames(power.Nmd) <- mVec

power.Nmd


# plot the power contour

meltd <- melt(power.Nm)
names(meltd) <- c("x", "y", "z")


Pd <- ggplot(meltd, aes(x, y, z=z))  +
  theme_classic()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  geom_tile(aes(fill = z)) + stat_contour(bins=5, size=1) +
  theme(plot.title = element_text(lineheight=.8, face="bold")) +
  xlab("Sample Size (N)") + ylab("Mean Count Difference") +
  guides(fill=guide_legend(title="Statistical Power")) +
  geom_segment(aes(x=3,xend=7,y=.133,yend=.133, color="red")) +
  geom_segment(aes(x=7,xend=7,y=0.075,yend=0.133, color="red")) +
  scale_colour_discrete(name = "S. diancaeus", labels=c("Real Data")) 
  
Pd



## What if we pool all of the minority taxa? -----------------

dataA <- subset(data, fish != "Stegastes partitus" & fish != "Stegastes spp" & fish != "Stegastes diancaeus" & fish != "Stegastes variabilis")

# get the parameters of our sample distribution

fA <- fitdistr(dataA$normDiff,"normal")

meanA <- as.numeric(fA$estimate[1])
sdA <- as.numeric(fA$estimate[2])
nA <- nrow(dataA)

# simulate data with these parameters

normA <- rnorm(mean=meanA, sd=sdA, n=nA)

# can we determine significance?

tA <- t.test(normA, mu=0, alternative = c("greater"))
tA$p.value




# Now we produce one of the nested power simulations

# cycle through range of N
Nvec = seq(5,50, by =5)
power.Na = numeric(length(Nvec)) #pre-allocate this variable to fill later

# 500 times each
nsim=500
pval = numeric(nsim)

# implement loop
for (j in 1:length(Nvec)) {
  N = Nvec[j]
  for (i in 1:nsim) {
    normSimA = rnorm(n=N, mean=meanA, sd=sdA)
    pval[i] = t.test(normSimA, mu=0, alternative = c("greater"))$p.value
  }
  power.Na[j] = sum(pval < 0.01)/nsim
}

power.Na

# lets also vary the mean slightly

mVec <- seq(.2,.5,by=.05)

power.Nma <- matrix(nrow=length(Nvec),ncol=length(mVec))

# implement nested loop
for (j in 1:length(Nvec)) {
  Nsa = Nvec[j]
  for (k in 1:length(mVec)) {
    meanSa = mVec[k]
    for (i in 1:nsim) {
      normSimSa <- rnorm(n=Nsa, mean=meanSa, sd=sdA)
      pval[i] = t.test(normSimSa, mu=0, alternative = c("greater"))$p.value
    }
    power.Nma[j,k] = sum(pval < 0.05)/nsim
  }
}
rownames(power.Nma) <- Nvec
colnames(power.Nma) <- mVec

power.Nma


# plot the power contour

melta <- melt(power.Nma)
names(melta) <- c("x", "y", "z")


Pa <- ggplot(melta, aes(x, y, z=z))  +
  theme_classic()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  geom_tile(aes(fill = z)) + stat_contour(bins=5, size=1,color="black") +
  theme(plot.title = element_text(lineheight=.8, face="bold")) +
  xlab("Sample Size (N)") + ylab("Mean Count Difference") +
  guides(fill=guide_legend(title="Statistical Power")) +
  geom_segment(aes(x=3,xend=nrow(dataA),y=meanA,yend=meanA, color="red")) +
  geom_segment(aes(x=nrow(dataA),xend=nrow(dataA),y=0.175,yend=meanA, color="red")) +
  scale_colour_discrete(name = "All other taxa", labels=c("Real Data")) 

Pa

# actual data p -value ?
tA <- t.test(dataA$normDiff, mu=0, alternative = c("greater"))
tA$p.value

## Lets try for entire Stegastes genus--------------

dataSt <- subset(data, fish == "Stegastes partitus" | fish == "Stegastes spp" | fish == "Stegastes diancaeus" | fish == "Stegastes variabilis")


# first, extract distribution paramters
fSt <- fitdistr(dataSt$normDiff,"normal")
meanSt <- as.numeric(fSt$estimate[1])
sdSt <- as.numeric(fSt$estimate[2])
nSt <- nrow(dataSt)

# cycle through range of N
Nvec = seq(5,50, by =5)
power.N = numeric(length(Nvec)) #pre-allocate this variable to fill later

# 500 times each
nsim=500
pval = numeric(nsim)

# implement loop
for (j in 1:length(Nvec)) {
  N = Nvec[j]
  for (i in 1:nsim) {
    normSimSt = rnorm(n=N, mean=meanSt, sd=sdSt)
    pval[i] = t.test(normSimSt, mu=0, alternative = c("greater"))$p.value
  }
  power.N[j] = sum(pval < 0.01)/nsim
}

power.N

## What if we vary the mean value and N? -------------

mVec <- seq(.1,.4,by=.05)

power.NmSt <- matrix(nrow=length(Nvec),ncol=length(mVec))

# implement nested loop
for (j in 1:length(Nvec)) {
  Nst = Nvec[j]
  for (k in 1:length(mVec)) {
    meanS = mVec[k]
    for (i in 1:nsim) {
      normSimSt <- rnorm(n=Nst, mean=meanS, sd=sdSt)
      pval[i] = t.test(normSimSt, mu=0, alternative = c("greater"))$p.value
    }
    power.NmSt[j,k] = sum(pval < 0.05)/nsim
  }
}
rownames(power.NmSt) <- Nvec
colnames(power.NmSt) <- mVec

power.NmSt


# plot the power contour

meltSt <- melt(power.NmSt)
names(meltSt) <- c("x", "y", "z")


PSt <- ggplot(meltSt, aes(x, y, z=z))  +
  theme_classic()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  geom_tile(aes(fill = z)) + stat_contour(bins=5, size=1,color="black") +
  theme(plot.title = element_text(lineheight=.8, face="bold")) +
  xlab("Sample Size (N)") + ylab("Mean Count Difference") +
  geom_segment(aes(x=3,xend=nrow(dataSt),y=meanSt,yend=meanSt, color="red")) +
  geom_segment(aes(x=nrow(dataSt),xend=nrow(dataSt),y=0.075,yend=meanSt, color="red")) +
  scale_colour_discrete(name = "Stegastes spp", labels=c("Real Data")) +
  guides(fill=guide_legend(title="Statistical Power"))

PSt

# whats actual p -value?

tSt <- t.test(dataSt$normDiff, mu=0, alternative = c("greater"))
tSt$p.value