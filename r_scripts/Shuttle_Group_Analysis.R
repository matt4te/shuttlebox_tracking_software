        ### Group-Level Analysis of the ShuttleBox Experiments ###

# to do: - add something to check for a side bias
#        - add a binomial test for end side
#        - look at the variation of activity level / speeds on the two sides


# Must be located in a directory containing the 'experiments' folder with numbered experiments inside. 
# Each should folder have an 'expSummary.csv' file produced by 'Shuttle_Track_Analysis.R'



## Setting up shop -------------------

# set the working directory and store path
dir <- "C:/Users/Matt/Desktop/Current_Lab_Work/ShadowBox"
setwd(dir) # if this is on a network drive you must first open the connection via Rstudio 'file' panel

# load the required libraries
library("plyr")
library("dplyr")
library("ggplot2")
library("gridExtra")
library("circular")
library("MASS")

# define function to concatenate all of the individual statistics files

grab_stats <- function(directory) {
  
  data <- "expSummary" # this is here in case the function is modified to have data as a function argument to load other files (e.g., tracks files)
  experiment = sort(as.numeric(dir(directory)))
  allData = data.frame()
  
  cat(paste("Reading", data,"for experiment    "))

  for (e in experiment) {
    # fancy terminal output
      n = nchar(as.character(e))
      backspaces = paste(rep("\b",n+1),collapse="")
      cat(backspaces,e)
      Sys.sleep(0.02)

    dataFile = paste(directory, "/", e, "/", data, ".csv", sep="")
    if (file.exists(dataFile)) {
      dat = read.table(dataFile, header=TRUE, sep=",", as.is=TRUE)
      dat$expNum = e
      allData = rbind.fill(allData, dat)
    }
  }
  return(allData)
}



## Load the data -----------------

# first data set
# stats <- grab_stats("experiments")# first, the preference statistics
stats <- read.csv("dms_stats.csv", stringsAsFactors=FALSE) # if you've saved the stats file already
log <- read.csv("Shuttle_Log.csv", stringsAsFactors=FALSE) # and the data log
colnames(log)[3] <- "expNum" # change the experiment id variable to keep consistent b/w data frames

data1 <- join(stats,log,by="expNum",type="inner")
data1$date <- as.POSIXct(data1$date,format="%m/%d/%Y")

# use only experiments where the fish crossed back and forth at least 3 times
data1 <- subset(data1, explorer == "TRUE" & status == "ok")


# newer data set
stats2 <- read.csv("dms_stats_2.csv", stringsAsFactors=FALSE) 
log2 <- read.csv("Shuttle_Log2.csv", stringsAsFactors=FALSE) 
colnames(log2)[3] <- "expNum" 
data2 <- join(stats2,log2,by="expNum",type="inner")
data2 <- data2[-7,] # exclude based on notes
data2$date <- as.POSIXct(data2$date,format="%m/%d/%Y")


# newest data set
stats3 <- read.csv("dms_stats_3.csv", stringsAsFactors=FALSE) 
log3 <- read.csv("Shuttle_Log3.csv", stringsAsFactors=FALSE) 
colnames(log3)[3] <- "expNum" 
data3 <- join(stats3,log3,by="expNum",type="inner")
data3$date <- as.POSIXct(data3$date,format="%m/%d/%Y")


# then combine the data

data <- rbind(data1,data2,data3)

data <- data[-2,] # exclude based on notes


## Plots for side preference  -------------

# density of count differences between control/treatment (positive skew indicates treatment preference)
dataMean <- mean(data$countDiff)
ggplot(data,aes(x=countDiff)) + geom_density(size=1,fill="black") +
  geom_vline(xintercept=0, color="red", size=1,alpha=0.5) +
  geom_vline(xintercept=dataMean, color="blue",size=1,alpha=.5) +
  xlab("Chamber Count Difference") + ylab("Probability") +
  theme_bw() +
  theme(axis.title.x = element_text(face="bold", size=16,vjust=-.1)) +
  theme(axis.title.y = element_text(face="bold", size=16,vjust=1)) +
  xlim(-900,900)


# and normalized version
data$normDiff <- data$countDiff / (data$treatCount + data$controlCount)
dataNormMean <- mean(data$normDiff)
dataNormMed <- median(data$normDiff)

ggplot(data,aes(x=normDiff)) + geom_density(size=1,fill="gray", colour="black") +
  geom_vline(xintercept=0, color="black", size=1,alpha=0.5,linetype=2) +
  geom_vline(xintercept=dataNormMean, color="black",size=1,alpha=.5) +
  xlab("Normalized Count Difference") + ylab("Relative Probability") +
  theme_bw() +
  theme(axis.title.x = element_text(face="bold", size=10,vjust=-.1)) +
  theme(axis.title.y = element_text(face="bold", size=10,vjust=1)) +
  theme(axis.text.y=element_blank(),axis.ticks=element_blank()) +
  theme(axis.text.x=element_text(face="bold",size=9)) +
  xlim(-1,1)



# bar graph for percentage of time spent in control V treatment 
controlPerc <- mean(data$controlPerc)
sd1 <- sd(data$controlPerc)
se1 <- sd1/sqrt(length(data$controlPerc))
treatPerc <- mean(data$treatPerc)
sd2 <- sd(data$treatPerc)
se2 <- sd2/sqrt(length(data$treatPerc))
barData <- c(controlPerc,treatPerc,sd1,sd2,se1,se2)
barTable <- as.data.frame(matrix(barData,nrow=2,ncol=3))
colnames(barTable) <- c("meanPerc","sd","se")
barTable$Condition <- c("Control","DMS")

ggplot(barTable, aes(x= Condition, y= meanPerc), colour="blue") + geom_bar(stat="identity",aes(fill=Condition)) +
  geom_errorbar(aes(ymin=meanPerc-se, ymax=meanPerc+se)) + 
  scale_fill_grey(start=.8, end=.4) +
  xlab("Experimental Condition") + ylab("Percent Time") +
  theme_bw() +
  theme(axis.title.x = element_text(face="bold", size=10,vjust=-.5)) +
  theme(axis.title.y = element_text(face="bold", size=10,vjust=1)) +
  theme(axis.text.x = element_text(face="bold",size=9)) 






## Statistics for side preference  -----------------------

# QUESTION: are count differences centered at zero? or is there statistical bias for control or treatment?

# testing preference using count differences 
null <- rnorm(n=nrow(data),mean=0) # null nypothesis - normal distribution with same variance as data
ks.test(data$countDiff,"pnorm",0,sd(data$countDiff), alternative="less")
# here the alternative is LESS, since we're asking whether the cdf of the data is BELOW that of the test
# distribution, which is equivalent to it being shifted to the right (i.e., greater values)

# this may not be an appropriate test in the case where the CDFs cross. see link below:
# http://blog.thegrandlocus.com/2014/09/mis-using-the-ks-test-for-p-hacking

# but note that we fix the standard deviation in the test to be the same as our distributions, so they WON'T ever
# cross - see graphically below

# visually compare the experimental distribution to the test distribution (cumulative distribution functions)
ggplot(data, aes(countDiff)) +
  stat_ecdf(geom="smooth") + 
  stat_function(fun = "pnorm",args=list(mean=0, sd=sd(data$countDiff)), geom= "smooth", colour = "red") +
  xlab("Value") + ylab("Cumulative Probability")

#and the normalized version
ggplot(data, aes(normDiff)) +
  stat_ecdf(geom="smooth") + 
  stat_function(fun = "pnorm",args=list(mean=0, sd=sd(data$normDiff)), geom= "smooth", colour = "red") +
  xlab("Value") + ylab("Cumulative Probability")

# test another way
wilcox.test(data$countDiff,mu=0,alt="greater")
# here the alternative hypothesis is GREATER since were asking directly if the values of our distribution tend
# to be greater than the test distribution

# or this way
t.test(data$countDiff, mu=0, alternative = c("greater"))


# and for the noramlized data (THIS IS WHAT SHOULD BE USED) ****************
t.test(data$normDiff, mu=0, alternative = c("greater"))


# or the non-parametric equivalent 
wilcox.test(data$normDiff,mu=0,alt="greater")



# and also, we can do a paired t test on the percentage in each
shapiro.test(data$treatPerc)
qplot(sample=data$treatPerc)
shapiro.test(data$controlPerc)
qplot(sample=data$controlPerc)

t.test(x=data$treatPerc,y=data$controlPerc,paired=TRUE)






## Plot for side activity differences  -------------

# reduce activity levels to account for sub-sampling
data$controlAct <- data$controlAct / 2
data$treatAct <- data$treatAct / 2

# create a data table for plot 1
controlAct <- mean(data$controlAct)
sd12 <- sd(data$controlAct)
se12 <- sd12/sqrt(length(data$controlAct))
treatAct <- mean(data$treatAct)
sd22 <- sd(data$treatAct)
se22 <- sd22/sqrt(length(data$treatAct))
barData2 <- c(controlAct, treatAct, sd12,sd22, se12,se22)
barTable2 <- as.data.frame(matrix(barData2,nrow=2,ncol=3))
colnames(barTable2) <- c("meanAct","sd","se")
barTable2$Condition <- c("Control","DMS")

# plot activity levels
ggplot(barTable2, aes(x= Condition, y= meanAct), colour="blue") + geom_bar(stat="identity",aes(fill=Condition)) +
  geom_errorbar(aes(ymin=meanAct-se, ymax=meanAct+se)) + 
  xlab("Experimental Condition") + ylab("Activity Level (cm/s)") +
  scale_fill_grey(start=.8, end=.4) +
  theme_bw() +
  theme(axis.title.x = element_text(face="bold", size=10,vjust=-.5)) +
  theme(axis.title.y = element_text(face="bold", size=10,vjust=1)) +
  theme(axis.text.x = element_text(face="bold",size=12)) 



## Plots for turning angles -------------
controlCP <- ggplot(data,aes(x=controlTA)) + coord_polar()  + 
  geom_point(aes(y=1)) + xlab("Control Turning Angles")
treatCP <- ggplot(data,aes(x=treatTA)) + coord_polar() + 
  geom_point(aes(y=1))  + xlab("Treatment Turning Angles")

grid.arrange(controlCP,treatCP,ncol=2)


# plot frequency of large V small turns


# normalize turn frequencies for the amount of time spent on that side

data$ncontrolTF45 <-data$controlTF45 / data$controlCount
data$ntreatTF45 <- data$treatTF45 / data$treatCount

data$ncontrolTF90 <- data$controlTF90 / data$controlCount
data$ntreatTF90 <- data$treatTF90 / data$treatCount
  
# first, small turns
controlTFs <- mean(data$controlTF45)
sdTFs <- sd(data$controlTF45)
seTFs <- sdTFs/sqrt(length(data$controlTF45))
treatTFs <- mean(data$treatTF45)
sdTFs2 <- sd(data$treatTF45)
seTFs2 <- sdTFs2/sqrt(length(data$treatTF45))
barDataTFs <- c(controlTFs, treatTFs, sdTFs,sdTFs2, seTFs,seTFs2)
barTableTFs <- as.data.frame(matrix(barDataTFs,nrow=2,ncol=3))
colnames(barTableTFs) <- c("meanTF","sd","se")
barTableTFs$condition <- c("control","treatment")

ggplot(barTableTFs, aes(x= condition, y= meanTF), colour="blue") + geom_bar(stat="identity",aes(fill=condition)) +
  geom_errorbar(aes(ymin=meanTF-se, ymax=meanTF+se)) + 
  xlab("Experimental Condition") + ylab("Turn Frequency") +
  theme_bw() +
  theme(axis.title.x = element_text(face="bold", size=13,vjust=-.5)) +
  theme(axis.title.y = element_text(face="bold", size=13,vjust=1)) +
  theme(axis.text.x = element_text(face="bold",size=12)) +
  scale_fill_brewer(palette = "Set1")

# noramlized version
ncontrolTFs <- mean(data$ncontrolTF45)
nsdTFs <- sd(data$ncontrolTF45)
nseTFs <- nsdTFs/sqrt(length(data$ncontrolTF45))
ntreatTFs <- mean(data$ntreatTF45)
nsdTFs2 <- sd(data$ntreatTF45)
nseTFs2 <- nsdTFs2/sqrt(length(data$ntreatTF45))
nbarDataTFs <- c(ncontrolTFs, ntreatTFs, nsdTFs,nsdTFs2, nseTFs,nseTFs2)
nbarTableTFs <- as.data.frame(matrix(nbarDataTFs,nrow=2,ncol=3))
colnames(nbarTableTFs) <- c("meanTF","sd","se")
nbarTableTFs$condition <- c("Control","DMS")

ggplot(nbarTableTFs, aes(x= condition, y= meanTF), colour="blue") + geom_bar(stat="identity",aes(fill=condition)) +
  geom_errorbar(aes(ymin=meanTF-se, ymax=meanTF+se)) + 
  xlab("Experimental Condition") + ylab("Turn Frequency") +
  theme_bw() +
  theme(axis.title.x = element_text(face="bold", size=13,vjust=-.5)) +
  theme(axis.title.y = element_text(face="bold", size=13,vjust=1)) +
  theme(axis.text.x = element_text(face="bold",size=12)) +
  scale_fill_brewer(palette = "Set1")

# next, large turns
controlTFl <- mean(data$controlTF90)
sdTFl <- sd(data$controlTF90)
seTFl <- sdTFs/sqrt(length(data$controlTF90))
treatTFl <- mean(data$treatTF90)
sdTFl2 <- sd(data$treatTF90)
seTFl2 <- sdTFs2/sqrt(length(data$treatTF90))
barDataTFl <- c(controlTFl, treatTFl, sdTFl,sdTFl2, seTFl,seTFl2)
barTableTFl <- as.data.frame(matrix(barDataTFl,nrow=2,ncol=3))
colnames(barTableTFl) <- c("meanTF","sd","se")
barTableTFl$condition <- c("control","treatment")

ggplot(barTableTFl, aes(x= condition, y= meanTF), colour="blue") + geom_bar(stat="identity",aes(fill=condition)) +
  geom_errorbar(aes(ymin=meanTF-se, ymax=meanTF+se)) + 
  xlab("Experimental Condition") + ylab("Turn Frequency") +
  theme_bw() +
  theme(axis.title.x = element_text(face="bold", size=13,vjust=-.5)) +
  theme(axis.title.y = element_text(face="bold", size=13,vjust=1)) +
  theme(axis.text.x = element_text(face="bold",size=12)) +
  scale_fill_brewer(palette = "Set1")

# noramlized version large
ncontrolTFl <- mean(data$ncontrolTF90)
nsdTFl <- sd(data$ncontrolTF90)
nseTFl <- nsdTFl/sqrt(length(data$ncontrolTF90))
ntreatTFl <- mean(data$ntreatTF90)
nsdTFl2 <- sd(data$ntreatTF90)
nseTFl2 <- nsdTFl2/sqrt(length(data$ntreatTF90))
nbarDataTFl <- c(ncontrolTFl, ntreatTFl, nsdTFl,nsdTFl2, nseTFl,nseTFl2)
nbarTableTFl <- as.data.frame(matrix(nbarDataTFl,nrow=2,ncol=3))
colnames(nbarTableTFl) <- c("meanTF","sd","se")
nbarTableTFl$condition <- c("control","treatment")

ggplot(nbarTableTFl, aes(x= condition, y= meanTF), colour="blue") + geom_bar(stat="identity",aes(fill=condition)) +
  geom_errorbar(aes(ymin=meanTF-se, ymax=meanTF+se)) + 
  xlab("Experimental Condition") + ylab("Turn Frequency") +
  theme_bw() +
  theme(axis.title.x = element_text(face="bold", size=13,vjust=-.5)) +
  theme(axis.title.y = element_text(face="bold", size=13,vjust=1)) +
  theme(axis.text.x = element_text(face="bold",size=12)) +
  scale_fill_brewer(palette = "Set1")


# TOTAL NORMALIZATION - 45 TURNS ***************
data$TcontrolTF45 <-data$controlTF45 / (data$controlCount+data$treatCount)
data$TtreatTF45 <- data$treatTF45 / (data$controlCount+data$treatCount)

TcontrolTFs <- mean(data$TcontrolTF45)
TsdTFs <- sd(data$TcontrolTF45)
TseTFs <- TsdTFs/sqrt(length(data$TcontrolTF45))
TtreatTFs <- mean(data$TtreatTF45)
TsdTFs2 <- sd(data$TtreatTF45)
TseTFs2 <- TsdTFs2/sqrt(length(data$TtreatTF45))
TbarDataTFs <- c(TcontrolTFs, TtreatTFs, TsdTFs,TsdTFs2, TseTFs,TseTFs2)
TbarTableTFs <- as.data.frame(matrix(TbarDataTFs,nrow=2,ncol=3))
colnames(TbarTableTFs) <- c("meanTF","sd","se")
TbarTableTFs$Condition <- c("Control","DMS")

ggplot(TbarTableTFs, aes(x= Condition, y= meanTF), colour="blue") + geom_bar(stat="identity",aes(fill=Condition)) +
  geom_errorbar(aes(ymin=meanTF-se, ymax=meanTF+se)) + 
  scale_fill_grey(start=.8, end=.4) +
  xlab("Experimental Condition") + ylab("Normalized Turn Frequency") +
  theme_bw() +
  theme(axis.title.x = element_text(face="bold", size=10,vjust=-.5)) +
  theme(axis.title.y = element_text(face="bold", size=10,vjust=1)) +
  theme(axis.text.x = element_text(face="bold",size=9)) 

## Statistics for side activity differences --------------


# QUESTION: is there a difference in the activity level of larvae in the treatment chamber?

# test for difference of mean activity level
t.test(x=data$treatAct,y=data$controlAct,paired=TRUE)

# and nonparametric
wilcox.test(x=data$treatAct,y=data$controlAct,paired=TRUE)




## Statistics for turning angles -----------------


# QUESTION: is there a similar distribution of turning angles in control and treatment?

watson.two.test(data$controlTA, data$treatTA)


# QUESTION: is there a difference in the frequency of turns in the treatment chamber?

# test for difference of mean activity level (normalized version)
t.test(x=data$TtreatTF45,y=data$TcontrolTF45,paired=TRUE)

# and nonparametric
wilcox.test(x=data$treatTF45,y=data$controlTF45,paired=TRUE)




## Equally-weighting each taxa -----------------

# first some turn normalization
data$TcontrolTF45 <-data$controlTF45 / (data$controlCount+data$treatCount)
data$TtreatTF45 <- data$treatTF45 / (data$controlCount+data$treatCount)

weighted <- ddply(data, ~fish, summarize, 
              N = length(ID),
              normDiff = mean(normDiff),
              meanCA = mean(controlAct),
              meanTA = mean(treatAct),
              meanCTF = mean(TcontrolTF45),
              meanTTF = mean(TtreatTF45),
              fish = fish)


t.test(weighted$normDiff, mu=0, alternative = c("greater"))

ggplot(weighted,aes(x=normDiff)) + geom_density(size=1,fill="gray", colour="black") +
  geom_vline(xintercept=0, color="black", size=1,alpha=0.5,linetype=2) +
  geom_vline(xintercept=dataNormMean, color="black",size=1,alpha=.5) +
  xlab("Normalized Count Difference") + ylab("Relative Probability") +
  theme_bw() +
  theme(axis.title.x = element_text(face="bold", size=10,vjust=-.1)) +
  theme(axis.title.y = element_text(face="bold", size=10,vjust=1)) +
  theme(axis.text.y=element_blank(),axis.ticks=element_blank()) +
  theme(axis.text.x=element_text(face="bold",size=9)) +
  xlim(-1,1)


# for activity 
# reduce activity levels to account for sub-sampling
weighted$meanCA <- weighted$meanCA / 2
weighted$meanTA <- weighted$meanTA / 2

# create a data table for plot 1
meanCA <- mean(weighted$meanCA)
sd12 <- sd(weighted$meanCA)
se12 <- sd12/sqrt(length(weighted$meanCA))
meanTA <- mean(weighted$meanTA)
sd22 <- sd(weighted$meanTA)
se22 <- sd22/sqrt(length(weighted$meanTA))
barData2 <- c(meanCA, meanTA, sd12,sd22, se12,se22)
barTable2 <- as.data.frame(matrix(barData2,nrow=2,ncol=3))
colnames(barTable2) <- c("meanAct","sd","se")
barTable2$Condition <- c("Control","DMS")


# plot activity levels
ggplot(barTable2, aes(x= Condition, y= meanAct), colour="blue") + geom_bar(stat="identity",aes(fill=Condition)) +
  geom_errorbar(aes(ymin=meanAct-se, ymax=meanAct+se)) + 
  xlab("Experimental Condition") + ylab("Activity Level (cm/s)") +
  scale_fill_grey(start=.8, end=.4) +
  theme_bw() +
  theme(axis.title.x = element_text(face="bold", size=10,vjust=-.5)) +
  theme(axis.title.y = element_text(face="bold", size=10,vjust=1)) +
  theme(axis.text.x = element_text(face="bold",size=12)) 

# test for difference of mean activity level
t.test(x=weighted$meanTA,y=weighted$meanCA,paired=TRUE)
wilcox.test(x=weighted$meanTA,y=weighted$meanCA,paired=TRUE)

# turning frequency
TcontrolTFs <- mean(weighted$meanCTF)
TsdTFs <- sd(weighted$meanCTF)
TseTFs <- TsdTFs/sqrt(length(weighted$meanCTF))
TtreatTFs <- mean(weighted$meanTTF)
TsdTFs2 <- sd(weighted$meanTTF)
TseTFs2 <- TsdTFs2/sqrt(length(weighted$meanTTF))
TbarDataTFs <- c(TcontrolTFs, TtreatTFs, TsdTFs,TsdTFs2, TseTFs,TseTFs2)
TbarTableTFs <- as.data.frame(matrix(TbarDataTFs,nrow=2,ncol=3))
colnames(TbarTableTFs) <- c("meanTF","sd","se")
TbarTableTFs$Condition <- c("Control","DMS")

ggplot(TbarTableTFs, aes(x= Condition, y= meanTF), colour="blue") + geom_bar(stat="identity",aes(fill=Condition)) +
  geom_errorbar(aes(ymin=meanTF-se, ymax=meanTF+se)) + 
  scale_fill_grey(start=.8, end=.4) +
  xlab("Experimental Condition") + ylab("Normalized Turn Frequency") +
  theme_bw() +
  theme(axis.title.x = element_text(face="bold", size=10,vjust=-.5)) +
  theme(axis.title.y = element_text(face="bold", size=10,vjust=1)) +
  theme(axis.text.x = element_text(face="bold",size=9)) 

# and statistics
t.test(x=weighted$meanTTF,y=weighted$meanCTF,paired=TRUE)
wilcox.test(x=weighted$meanTTF,y=weighted$meanCTF,paired=TRUE)


### FOR MAHI DATA --------

# set directory
dir <- "B:/Science/phd_data/mahi_shuttle"
setwd(dir)

# load stats
stats <- grab_stats("experiments")# first, the preference statistics

# read log
log <- read.csv("Mahi_Shuttle_Log.csv", stringsAsFactors=FALSE) # and the data log
colnames(log)[3] <- "expNum" 

# combine log and stats
dataM <- join(stats,log,by="expNum",type="inner")
dataM$date <- as.POSIXct(dataM$date,format="%m/%d/%Y")

# deal with weird count diff thing
dataM[9,]$treatPerc <- 0
dataM[9,]$treatCount <- 0
dataM[dataM$treatPerc==100,]$countDiff <- dataM[dataM$treatPerc==100,]$treatCount
dataM[dataM$controlPerc==100,]$countDiff <- -dataM[dataM$controlPerc==100,]$controlCount


# normalized side difference plot
dataM$normDiff <- dataM$countDiff / (dataM$treatCount + dataM$controlCount)
dataMNormMean <- mean(dataM$normDiff)
dataMNormMed <- median(dataM$normDiff)

ggplot(dataM,aes(x=normDiff)) + geom_density(size=1,fill="gray", colour="black") +
  geom_vline(xintercept=0, color="black", size=1,alpha=0.5,linetype=2) +
  geom_vline(xintercept=dataMNormMean, color="black",size=1,alpha=.5) +
  xlab("Normalized Count Difference") + ylab("Relative Probability") +
  theme_bw() +
  theme(axis.title.x = element_text(face="bold", size=10,vjust=-.1)) +
  theme(axis.title.y = element_text(face="bold", size=10,vjust=1)) +
  theme(axis.text.y=element_blank(),axis.ticks=element_blank()) +
  theme(axis.text.x=element_text(face="bold",size=9)) +
  xlim(-1,1)


# bar graph for percentage of time spent in control V treatment 
controlPerc <- mean(dataM$controlPerc)
sd1 <- sd(dataM$controlPerc)
se1 <- sd1/sqrt(length(dataM$controlPerc))
treatPerc <- mean(dataM$treatPerc)
sd2 <- sd(dataM$treatPerc)
se2 <- sd2/sqrt(length(dataM$treatPerc))
bardataM <- c(controlPerc,treatPerc,sd1,sd2,se1,se2)
barTable <- as.data.frame(matrix(bardataM,nrow=2,ncol=3))
colnames(barTable) <- c("meanPerc","sd","se")
barTable$Condition <- c("Control","DMS")

ggplot(barTable, aes(x= Condition, y= meanPerc), colour="blue") + geom_bar(stat="identity",aes(fill=Condition)) +
  geom_errorbar(aes(ymin=meanPerc-se, ymax=meanPerc+se)) + 
  scale_fill_grey(start=.8, end=.4) +
  xlab("Experimental Condition") + ylab("Percent Time") +
  theme_bw() +
  theme(axis.title.x = element_text(face="bold", size=10,vjust=-.5)) +
  theme(axis.title.y = element_text(face="bold", size=10,vjust=1)) +
  theme(axis.text.x = element_text(face="bold",size=9)) 


# side preference statistics

# shapiro.test(dataM$normDiff) # the data is  NOT normal

wilcox.test(dataM$normDiff,mu=0,alt="greater") # so use non-parametric test


# side activity differences
controlAct <- mean(dataM$controlAct,na.rm=TRUE)
sd12 <- sd(dataM$controlAct,na.rm=TRUE)
se12 <- sd12/sqrt(length(dataM$controlAct))
treatAct <- mean(dataM$treatAct,na.rm=TRUE)
sd22 <- sd(dataM$treatAct,na.rm=TRUE)
se22 <- sd22/sqrt(length(dataM$treatAct))
barData2 <- c(controlAct, treatAct, sd12,sd22, se12,se22)
barTable2 <- as.data.frame(matrix(barData2,nrow=2,ncol=3))
colnames(barTable2) <- c("meanAct","sd","se")
barTable2$Condition <- c("Control","DMS")

# plot activity levels
ggplot(barTable2, aes(x= Condition, y= meanAct), colour="blue") + geom_bar(stat="identity",aes(fill=Condition)) +
  geom_errorbar(aes(ymin=meanAct-se, ymax=meanAct+se)) + 
  xlab("Experimental Condition") + ylab("Activity Level (cm/s)") +
  scale_fill_grey(start=.8, end=.4) +
  theme_bw() +
  theme(axis.title.x = element_text(face="bold", size=10,vjust=-.5)) +
  theme(axis.title.y = element_text(face="bold", size=10,vjust=1)) +
  theme(axis.text.x = element_text(face="bold",size=12)) 

# stats
wilcox.test(x=dataM$controlAct,y=dataM$treatAct,paired=TRUE)


# normalized turn frequency plot

dataM$TcontrolTF45 <-dataM$controlTF45 / (dataM$controlCount+dataM$treatCount)
dataM$TtreatTF45 <- dataM$treatTF45 / (dataM$controlCount+dataM$treatCount)

TcontrolTFs <- mean(dataM$TcontrolTF45,na.rm=TRUE)
TsdTFs <- sd(dataM$TcontrolTF45,na.rm=TRUE)
TseTFs <- TsdTFs/sqrt(length(dataM$TcontrolTF45))
TtreatTFs <- mean(dataM$TtreatTF45,na.rm=TRUE)
TsdTFs2 <- sd(dataM$TtreatTF45,na.rm=TRUE)
TseTFs2 <- TsdTFs2/sqrt(length(dataM$TtreatTF45))
TbarDataTFs <- c(TcontrolTFs, TtreatTFs, TsdTFs,TsdTFs2, TseTFs,TseTFs2)
TbarTableTFs <- as.data.frame(matrix(TbarDataTFs,nrow=2,ncol=3))
colnames(TbarTableTFs) <- c("meanTF","sd","se")
TbarTableTFs$Condition <- c("Control","DMS")

ggplot(TbarTableTFs, aes(x= Condition, y= meanTF), colour="blue") + geom_bar(stat="identity",aes(fill=Condition)) +
  geom_errorbar(aes(ymin=meanTF-se, ymax=meanTF+se)) + 
  scale_fill_grey(start=.8, end=.4) +
  xlab("Experimental Condition") + ylab("Normalized Turn Frequency") +
  theme_bw() +
  theme(axis.title.x = element_text(face="bold", size=10,vjust=-.5)) +
  theme(axis.title.y = element_text(face="bold", size=10,vjust=1)) +
  theme(axis.text.x = element_text(face="bold",size=9)) 



# and nonparametric
wilcox.test(x=dataM$treatTF45,y=dataM$controlTF45,paired=TRUE)



## Printing figures for publication ----------------


setwd("C:/Users/Matt/Desktop/new_dms_revisions/figures")
# file = “ “ part of the functions below.

tiff("FileName.tiff", height = 12, width = 17, units = 'cm', compression = "lzw", res = 300)

## NOW PRODUCE THE PLOT 

dev.off()

