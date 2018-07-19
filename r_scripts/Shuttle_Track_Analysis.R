     ### Track and Analyze the Tracks from Individual ShuttleBox Trials ###

# to do:
#   - after tracking, plot all of the turning angles (y) versus frequency (x) to see what they're mostly doing
#   - try downscaling the angle stuff to see what we get
#   - change this code to mark if in the middle and not in a side


## SETTING UP SHOP -----------

# Step1: Ensure that your working directory is the directoy of the project. 
# It should include:
# 
# (1) an "experiments" folder with pre-processed data sorted by trial #
# 
# Pre-processing is done in as follows:
#   
# ----- One time only
# 1. Download the 32-bit version of virtual dub (windows only software)
# 2. Download the 'DshowInputDriver' plugin, and place in the plugins32 folder of virtual dub
# 
# ------ For each trial
# 3. Open virtual dub, and open video file recorded from gopro
# 4. Navigate to video -> framerate -> and select to decimate the framerate by 30 (go from 30 fps to 1 fps)
# 5. Manually select the start and end time of the 'trial period', excluding acclimation time and any extra time at end of video
# 6. Navigate to edit -> crop video to selection
# 7. Navigate to file -> export to image sequences (choose .bmp format), and save them in the directory with the original video



# store the path to the project directory for use later 
dir <- "B:/Science/phd_data/mahi_shuttle"

# and make it your working directory
setwd(dir) 


# Step2: Edit the Shuttle_Tracking_Macro.ijm manually 
  # This should be located in 'C:\Program Files (x86)\ImageJ\macros' 

  # the first line should include 'open=.....' where you specify the location of the pictures in your "working" folder (see above)
  # this first line should also include 'number=....' of frames you'd like to analyze (based on the fps and length of trial period)
  # the 4th line should the include the same path inserted into line 1 (keep the chosen filename as is)
  # repeat for lines 7, 10, and 18 (keep the chosen filenames as is)

  # NOTE: keep the formatting the same! (i.e., direction and number of slashes, \\)


# Step3:  load required libraries and define some functions
library("plyr")
library("stringr")


# define a function to monitor the status of ImageJ system calls
check_status <- function(status, message="ImageJ error or unexpected termination", ...) {
  if ( status != 0 ) {
    stop(message, call.=FALSE)
  }
  return(invisible(status))
}

# define function to calculate bearing between two points (constant frame of reference)
angleFun <- function(xx,yy){
  ## xx and yy are the differences in x and y coordinates between two points
  c = 180/pi  
  b<-sign(xx)
  b[b==0]<-1  #corrects for the fact that sign(0) == 0
  tempangle = b*(yy<0)*pi+atan(xx/yy)
  tempangle[tempangle<0]<-tempangle[tempangle<0]+2*pi
  return(tempangle*c)
}


# this is the function you will explicitly call to analyze the data
AnalyzeExperiment <- function(exp, control,s) {
  
#   There are three inputs:
#   1. exp - this is the name of the folder (ideally a number) that houses the image sequences
#   2. control - this is either a 1 or a 2, indiciting that the "control" treatment is on the left(1) or right(2) side
        # note: sometimes there will be 2 treatments rather than treatment V control. in this case... bear this in mind when choosing 1 or 2, as outputs will still read "control" and "treatment"
#   3. s - this is the number of sections that the character "\" would split the path to your images into
        # for example, Y:\shuttle_box\shuttle_bk\experiments\4\sequence0000.bmp" would be split into 6 sections
  
  # friendly message
  cat("\n Moving Experiment", str_c(exp), "data to working directory \n\n")
  
  setwd(dir) # make sure the working directory is reset before each round of analysis
  
  # set the path to the experiment folder's video file
  expPath <- str_c(getwd(), "/experiments/", exp)
  
  # and specify the path of the proper working directory for ImageJ
  workPath <- str_c(getwd(), "/working")
  dir.create(workPath,showWarnings=FALSE)
  
  
  # capture the paths of the image sequence in the experiment directory
  sequence = list.files(path=expPath, pattern= "*.bmp", full.names=TRUE)
  # and isolate the name of the image file separated from the full path
  sepa <- t(as.data.frame(str_split(sequence,"/",n=s)))  # YOU WILL NEED TO CHANGE N=6 based on your path chosen
  
  # create a hard link in the working directory of the images in the experiment directory
  for (i in 1:length(sequence)){
    file.link(from=sequence[i], to=str_c(workPath,"/",sepa[i,s],sep=""))
  }
  
  cat("\n Opening images for calibration and tracking \n\n")
  
  # Make a system call to run the ImageJ macro for tracking on the video in the working directory
  command <- 'java -Xmx1000m -jar "C:/Program Files (x86)/ImageJ/ij.jar" -ijpath "C:/Program Files (x86)/ImageJ/" -macro Shuttle_Tracking_Macro.txt'  # define the system call in advance
  status <- system(command, wait=TRUE, invisible=FALSE)  # store the system call as you would feed it to R
  check_status(status, message=str_c("Error running command:\n\n", command, "\n\nAbort"))  # feed the system call to R through the function defined above that gives status information on the command
  
  
  
#   NOW ALL OF THE TRACKING HAPPENS IN IMAGEJ  
  
  
  
  # prompt for indication that the tracking went well and you're ready to proceed with stats
  check <- as.character(readline("Did the tracking go well? Ready to proceed with analysis? y/n:  ")) 
  
  if (check == "y") {
  # set new R working directory to be the project working directory (where ImageJ has saved its output)
    setwd(workPath)
  } else { stop("Analysis Cancelled due to Bad Tracking") }
  
  
  # load the data from the ImageJ tracking
  center_line_endpoints <- read.table("center_coords.txt")
  left_boundary <- read.table("left_coords.txt")
  right_boundary <- read.table("right_coords.txt")
  tracks <- read.csv("tracks.csv", header=TRUE)
  
  
  # Reformat the center line data
  center_line_endpoints[,2] <- -center_line_endpoints[,2] # make the y-values negative, as is the convention in ImageJ
  colnames(center_line_endpoints) <- c("X", "Y")
  cle <- center_line_endpoints # for clarity below
  
#   # use this only when the chamber separation is diagnal and not horizontal (typical)
#   m <- (cle[2,2]-cle[1,2])/(cle[2,1]-cle[1,1]) # determine the slope of the line
#   b <- cle[1,2]-m*cle[1,1] # and y-intercept
#   # where the equation of the chamber bisecting line is y=mx+b, of course
  
  # find the x-value for the chamber separation (for a horizontally oriented shuttleBox)
  midX <- (cle[2,1]+cle[1,1])/2

  # Reformat the left boundary data
  left_boundary[,2] <- -left_boundary[,2] # make the y-values negative, as is the convention in ImageJ
  colnames(left_boundary) <- c("X", "Y")
  lb <- left_boundary

  # find the x-value for the start of the left chamber
  leftX <- (lb[2,1]+lb[1,1])/2

  # Reformat the right boundary data
  right_boundary[,2] <- -right_boundary[,2] # make the y-values negative, as is the convention in ImageJ
  colnames(right_boundary) <- c("X", "Y")
  rb <- right_boundary

  # find the x-value for the start of the left chamber
  rightX <- (rb[2,1]+rb[1,1])/2


  # Reformat the tracks data
  
  trackz = tracks[,3:7] # subset columns of interest
  colnames(trackz) <- c("frame", "x", "y", "distance", "velocity") # rename them
  trackz$y <- -trackz$y # make the y values negative
  
#   # again, this is only for the diag oriented chamber
#   trackz$thresholdY <- m*trackz$x + b # calculate the y value marking the threshold between the two chambers, based on the corresponding x location of the larva
#   trackz$chamber <- 1 # create a variable to store what side of the chamber the larva was in
#   for (i in 1:nrow(trackz)){
#     if (trackz[i,c("y")] < trackz[i,c("thresholdY")]) { # and determine which side by comparing the y value at that time to the thresholdY value calculated above
#       trackz[i,c("chamber")] <- 2
#     }
#   }

  trackz$chamber <- 1 # create a variable to store what side of the chamber the larva was in, default is LEFT
  for (i in 1:nrow(trackz)){
    if (trackz[i,c("x")] > rightX) { # check if any point's x value is greater than the right boundary
      trackz[i,c("chamber")] <- 2 # if so, then its on the RIGHT (value of 2)
    } else if (trackz[i,c("x")] < leftX) { # check if any point's x value is less than the left boundary
      trackz[i,c("chamber")] <- 1 # if so, then its on the LEFT side (value of 1, same as default)
    } else {
      trackz[i,c("chamber")] <- 3 # otherwise (if its in between those two values), its in the middle of the chambers
    }
  }

  trackz$chamber <- as.factor(trackz$chamber)
  
  
  
  ## Some Summary Statistics for Preference --------------
  
  
  # percentage of time spent on each side
  counts <- table(trackz$chamber) # count the number of occurances in each chamber
  percentage1 <- counts[1]/nrow(trackz)*100 # 1 is the left chamber 
  percentage2 <- counts[2]/nrow(trackz)*100 # 2 is the right chamber
  
  # The difference in counts 
  countDiff <- counts[2] - counts[1]
  
  # and  use that to a assign a + or - value to the countDiff (positive for treatment preference, negative for control preference)
  if (control == 2) {
    countDiff = -countDiff # countDiff is positive for a treatment in the bottom, so we reverse this is the control was in the bottom
  }
  
  # and to assign control/treatment labels to the counts/percentages calculated above
  
  if (control == 1) {
    controlCount <- counts[1]
    treatCount <- counts[2]
  } else if (control ==2) {
    controlCount <- counts[2]
    treatCount <- counts[1]
  }
  
  if (control == 1) {
    controlPerc <- percentage1
    treatPerc <- percentage2
  } else if (control ==2) {
    controlPerc <- percentage2
    treatPerc <- percentage1
  }
  
  # and also record the number of counts in the middle of the chambers
  midCount <- counts[3]
  
  ## Statistical Test for Individual Side Preference ----------------------
  
  if (!is.na(treatCount) & !is.na(controlCount)) {
  
  # Exact Binomial test
  n <- nrow(trackz)
  binom <- binom.test(treatCount, n, p=0.5, alternative = "two.sided")
  bp <- binom$p.value
  
  # Chi-squared
  chi <- chisq.test(c(treatCount,controlCount))
  cp <- chi$p.value
  
  # and save a variable indicating the results of stats
  preference <- logical()
  if (max(bp,cp) < 0.05) {
    preference = TRUE
  } else { preference = FALSE }
  
  } else {
    bp <- NA
    cp <- NA
    preference <- TRUE
    if (is.na(controlCount)) {
      controlCount <- 0 
      controlPerc <- 0
    } else if (is.na(treatCount)) {
      treatCount <- 0
      treatPerc <- 0
    }
  }
  
  ## Determine how many "crossings" from side <-> side the larva made -----------------
  
  crosses <- 0
  for (i in 1:(nrow(trackz)-1)) { # for all points but the last one
    c <- as.numeric(trackz[i,c("chamber")]) - as.numeric(trackz[i+1,c("chamber")]) # look to see if the larva changes chambers at the next step
    if (c != 0) { # if it does
      crosses <- crosses +1 # increment the crosses variable by 1
    }
  }
  
  
  # and set a threshold for how many are necessary to constitute a "good trial"
  if (crosses >= 3) { # 3 crosses seems good? 
    explorer = TRUE
  } else { explorer = FALSE }
  
  
  ## Compare activity level on both sides ----------------
  
  activity <- ddply(trackz, ~chamber, summarise, mean_velocity=mean(velocity),total_distance=sum(distance))
  
  # and assign control/treatment labels to the velocities & distances
  if (control == 1) {
    controlAct <- activity[1,2]
    controlD <- activity[1,3]
    treatAct <- activity[2,2]
    treatD <- activity[2,3]
  } else if (control ==2) {
    controlAct <- activity[2,2]
    controlD <- activity[2,3]
    treatAct <- activity[1,2]
    treatD <- activity[1,3]
  }

  if (treatCount > 0 & controlCount > 0) {  
  
  # test the difference in activity levels 
  speedAnova <- anova(lm(velocity ~ chamber, data=trackz))
  speedP <- speedAnova$"Pr(>F)"[1]
  
  } else {speedP <- NA}
 
  
  ## Compare turning frequency and turning angle --------------
   
  # calculate the deltaX and deltaY for each set of consecutive points
  trackz$xx <- 0 # using new columns
  trackz$yy <- 0
  
  for (i in 1:(nrow(trackz)-1)) { # for all points but the last one
    trackz[i,c("xx")] <- trackz[i+1,c("x")] - trackz[i,c("x")]
    trackz[i,c("yy")] <- trackz[i+1,c("y")] - trackz[i,c("y")]
  }
  
  # now calculate the bearing travelled at each step
  trackz$bearing <- 0
  
  for (i in 1:nrow(trackz)) {
    trackz[i,c("bearing")] <- angleFun(trackz[i,c("xx")],trackz[i,c("yy")]) # use the function above
    if (is.na(trackz[i,c("bearing")])) { # if the larva didn't move at all from one frame to another
      trackz[i,c("bearing")] = trackz[i-1,c("bearing")] # use the previous bearing 
    } 
    if (trackz[i,c("bearing")] == 0) { # if the bearing is north
      trackz[i,c("bearing")] <- 360 # let that value be 360 instead of zero (for use in the angle calculation below)
    }
  }
  
  
  # finally, calculate the turning angle at each step based on consecutive bearings (larva's last bearing frame of reference)
  trackz$turnAngle <- 0
  
  for (i in 1:(nrow(trackz)-1)) {
    angle <- trackz[i+1,c("bearing")] - trackz[i,c("bearing")] # subtract each bearing by the one before it
    if (angle <= 0) {
      angle = angle + 360 # and correct for any negative values
    }
    trackz[i,c("turnAngle")] <- angle # this is the larva's instantaneous turning angle
  }
  
  # note: this calculates turning angle between the current time step and the NEXT step step, not the previous
  
  
  # and summarize the turning angles for the entire deployment
  angling <- ddply(trackz, ~chamber, summarise, mean_angle=mean(turnAngle),total_turns_45=sum(turnAngle>45), total_turns_90=sum(turnAngle>90))
  
  # a "turn" here is arbitrarily defined as a rotation of > 45 degrees and a "big turn" is defined at 90 degrees
  
  # and assign control/treatment labels to the average turning angle and turning frequency
  if (control == 1) {
    controlTA <- angling[1,2]
    controlTF45 <- angling[1,3]
    controlTF90 <- angling[1,4]
    treatTA <- angling[2,2]
    treatTF45 <- angling[2,3]
    treatTF90 <- angling[2,4]
  } else if (control ==2) {
    controlTA <- angling[2,2]
    controlTF45 <- angling[2,3]
    controlTF90 <- angling[2,4]
    treatTA <- angling[1,2]
    treatTF45 <- angling[1,3]
    treatTF90 <- angling[1,4]
  }
  
  
  ## Save all of the track information -------------
  expSummary <- data.frame(controlPerc,treatPerc,controlCount,treatCount,countDiff,midCount,bp,cp,preference,controlAct,treatAct,speedP, controlD, treatD, controlTA,controlTF45,controlTF90, treatTA,treatTF45,treatTF90, crosses,explorer)
  write.csv(expSummary, "expSummary.csv", row.names=FALSE)
  
  # now move all of the work BACK to the experiment folder
  file.rename(from=str_c(workPath,"/center_coords.txt"), to=str_c(expPath,"/center_coords.txt"))
  file.rename(from=str_c(workPath,"/left_coords.txt"), to=str_c(expPath,"/left_coords.txt"))
  file.rename(from=str_c(workPath,"/right_coords.txt"), to=str_c(expPath,"/right_coords.txt"))
  file.rename(from=str_c(workPath,"/tracks.csv"), to=str_c(expPath,"/tracks.csv"))
  file.rename(from=str_c(workPath,"/expSummary.csv"), to=str_c(expPath,"/expSummary.csv"))



  # and delete the working directory contents
  unlink(workPath,recursive=TRUE,force=TRUE)
  
  #print results in terminal
  expSummary
  
  # and set the working directory back to the project's main directory
  setwd(dir)
  
cat("\n Tracking and Analysis for Experiment", str_c(exp), "was completed sucessfully \n\n")

}


## Now Just call the function -------------------

# For example:

# > AnalyzeExperiment(1,2,6)  # which would analyze the 1st experiment folder, specifying the control is on the right, and that "\" splits the path to images into 6 sections