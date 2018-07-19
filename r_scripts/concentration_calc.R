# script to figure out the concentration of DMS in the shuttlebox, based on....

# volume of one side is equal to 175 mL
# flow rate is equal to .5 ml/sec
# concentration at inflow is 340 nM

# see paper for solution to differential equation that yields the equation in the function below:


# concentration, in mM, as a function of time, in seconds
concentration <- function(t) {
  M = (-0.000000595*exp(-0.00286*t) + 0.000000595)/.175
  mM = M * 1000000
  print(mM)
}


# loop through the course of an experiment

# acclimation period
concentration(300)
# 1.958369     so the concentration is already over half of the stock solution after acclimation period

# entire experiment
c <- vector(length=1200)
for (t in 1:1200){ # from the start of the experiment to the end of 20 minutes, includes acclimation time
  c[t] <- concentration(t)
}



# what was the paramaterized form of the equation?

shuttleConc <- function(V,Cin,F,t) {
  
  # V is the volume of the aquarium, in L
  # Cin is the concentration of the stock solution of odor, in mol/L
  # F is the flow rate of the inlet/outlet, in L/s
  # t is time 
  
  k <- -V*Cin # integration constant
  
  M <- (k*exp(-F/V*t) + V*Cin)/V
  mM <- M * 1000000
  print(mM)
  
}


V = 175 / 1000 # 175 mL
Cin = 3.4 / 1000000 #3.4 mM
F = .5 / 1000 # .5 mL/s


c <- vector(length=1200)
for (t in 1:1200){ # from the start of the experiment to the end of 20 minutes, includes acclimation time
  c[t] <- shuttleConc(V,Cin,F,t)
}

plot(c)

# basically we see that at the end of the acclimation period we're already past half the concentration in
#  the stock solution. so we want the stock solution to be just greater than the desired test concentration

Caccl <- c[300]
Cf <- c[1200]

AcclPerc <- (Caccl/Cf)*100
# 59.49224 %



## publication plot------------------

library("ggplot2")

d <- as.data.frame(c)
d$s <- seq(1:1200)
d$c <- d$c * 100

ggplot(d, aes(x=s,y=c)) + geom_line(size=2) +
  geom_vline(x=300,color="black") +
  xlab("Time(s)") + ylab("Concentration(nM)") +
  theme_bw() +
  theme(axis.title.x = element_text(face="bold", size=10)) +
  theme(axis.title.y = element_text(face="bold", size=10,vjust=1))

