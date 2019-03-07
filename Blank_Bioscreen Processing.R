############ Bioscreen Processing ###############

## This script takes in bioscreen c raw data that has been sep with with 3 technical replicates
## arranged from top to bottom and left to right. A key must be provided that aligns
## the columns of the raw files (using the numbering id, V2, V3, V4... Vn) and the 
## conditions, strains ect that you want to be considered 
## This script will find the averages of the tech replicates, blank all the data, refromat the
## data to long (from wide), produce a scatter and smooth plot of the overall data
## and finally determine the doubling time (in min) for each strain/condition and 
## display the doubling time as a boxplot. This script does not yet test for a significant 
## change in doubling time

################################## 
# Fill out the following variables
###################################

file_name = "Example_BioscreenExperimentData.csv"
key_name = "Example_Bioscreen_key.csv"
saved_file_name ="Example_Tidy_output.csv"
saved_raw_name = "raw.dataforcurver.csv"

var = 6 ## total number of coditions/doses/stains
rep = 6 ## biological replicates 

### Set the working directory to your project folder
setwd("~/Desktop/Bioscreen Example")

## Which columns in the raw data file are the blanks?
## Are your blankes at the bottom row and do you have 6 bio reps? If yes, set as TRUE 
## If your blanks are not, then fill in "blanks" vector
Blank_row = TRUE
blanks <- c("includes time (n+1)")

## How many total columns from the data file do you want to read in?
number_of_wells = 120

## What was the time interval of the Bioscreen run?
interval = 15

##################################################
## If you have more than 3 tech replicates, this scripr will not work. 
## Change for loop in line 77 if there are more than 3 tech rep
###### Needed Packages ######

library(tibble)
library(tidyr)
library(cowplot)
library(growthrates)
library(dplyr)
library(reshape2)

############################

## importing data
bioscreen<-read.csv(file_name)

## remove unused rows off the end
bioscreen_edit <-bioscreen[,1:(number_of_wells+1)]

## Set up Blanks 
if(Blank_row == TRUE) {
  blank_positions <- seq(from = 11, to = ((var*20)+1), by = 10)
} else {
  blank_positions = blanks
}

blanks <- subset(bioscreen, select=c(blank_positions))
blanks$average <- c(rowMeans(blanks))
## can be accesed via blanks[average]

tech_data <- bioscreen_edit
## remove blanks and time 
tech_data <- tech_data[ -(blank_positions) ]
tech_data <- tech_data[-c(1)]

# Find the row average of every three columns
tech_avg = data.frame(id= integer(nrow(tech_data)))
tech_avg_c = data.frame(id= integer(nrow(tech_data)))
tech_l = (var * rep)

i=1 
j=1
for(i in 1:tech_l) {
  tech_avg[,i] <- rowMeans(tech_data[j:(j+2)]) ## number off rep - 1
  j = j+ 3 ## increase to # of rep
  i = i + 1 }

## Blank the data, average tec replicates, subtract blank average 
for (i in 1:tech_l) {
    tech_avg_c[,i] <- (tech_avg[,i] - blanks$average)
}

## plotting prep
total_time = (nrow(tech_avg_c)-1)*15
timeline <- seq(0.0,total_time,interval)
tech_avg_c$Time <- timeline 

## http://bconnelly.net/2014/04/analyzing-microbial-growth-with-r/
write.csv(tech_avg_c, saved_file_name) # for making plate key

tech_avg_c <- melt(tech_avg_c, id=c("Time"), variable.name="Well",
                   value.name="OD600")
platemap <- read.csv(key_name) ## defined locations and strains

Long_data <- inner_join(tech_avg_c, platemap, by="Well")

############### Final long format table ###############

write.csv(Long_data, saved_file_name)

################ Inital Plot of the data ################
## You can change the range of the graph using the filer commands
## The ggplot command must be updated for the x, y and color variables 
########################################################

Long_data %>% filter(Time >= 15) -> filtered_data
filtered_data %>% filter(Time <= 1000 ) -> filtered_data

Curve_plotA<- ggplot(filtered_data, aes(x=(Time/60), y=OD600, color=Dose)) +
  geom_point(size=.50) +
  labs(x="Time (Hours)", y="Abs. at OD600 nm") +
  labs(title = "All Growth Data\n by Strain and Stress ") + 
  theme(axis.text.x=element_text(angle=45, vjust =1, hjust = 1, size = 20)) +
  theme(axis.text.y=element_text(size = 20))+
  theme(axis.title.y=element_text(size = 20))+
  theme(axis.title.x=element_text(size = 20))+
  theme(plot.title = element_text(size = 25))
Curve_plotA 

################# Smooth Plot ##############
## You can change the range of the graph using the filer commands
## The ggplot command must be updated for the x, y and color variables 
########################################################

Long_data %>% filter(Time >= 15) %>% filter(Time <= 1000 ) -> filtered_data

Curve_plotB<- ggplot(data=filtered_data, aes(x=(Time/60), y=OD600, color=Dose)) +
  geom_smooth(size=1.5, span = 0.3) +
  labs(x="Time (Hours)", y="Abs. at OD600 nm") +
  labs(title = "Average Growth \n by Stress") + 
  theme(axis.text.x=element_text(angle=45, vjust =1, hjust = 1, size = 20)) +
  theme(axis.text.y=element_text(size = 20))+
  theme(axis.title.y=element_text(size = 20))+
  theme(axis.title.x=element_text(size = 20))+
  theme(plot.title = element_text(size = 22))
Curve_plotB

##### Export scatter plots ######
## Update file name
##################################
final_figure <- plot_grid(Curve_plotB, Curve_plotA, ncol ncol = 2)
ggsave(plot = final_figure , file = "Example_Growth_curves.pdf", device= "pdf", width= 16, height=6, units="in") 


############ Max growth rate #########
## This section will determine the fastest growth rates for each strain and then 
## produce a box and wisker to discribe the variation in doubling time. 
## This works best if you filter out the later timepoints. 
## Manually check the individual plots that are output to ensure that the correct 
## regon is being measured by growthrates. 
## There are two options here, run both and see which fits your data better
#######################################


################################################################################################################
################     All_slpines       ##################
## fits to the whole curve 
##################################################################

Long_data %>% filter(Time <= 400) -> more100
more100 %>% filter(OD600 >= 0) -> more100

all_splines_data <- all_splines(OD600 ~ Time | Dose + Well, data = more100, spar = .5)
plot(all_splines_data)
coef(all_splines_data)
growthrates <-as.data.frame(coef(all_splines_data))

colnames(growthrates) <-c("y0","Max_Growth_Rate") ##name columns

growthrates<-rownames_to_column(growthrates, var = "Dose") ##make names a column

growthrates %>% separate(Dose, c("Dose","Well")) -> DoublingData

### Filter out an strains that didn't grow

DoublingData %>% filter(Dose != "KD10kGy" ) %>% filter(Dose != "KDsham" ) -> DoublingData


### Doubling time
dt_v = vector('numeric')
len = nrow(DoublingData)
for (i in 1:(len)){
  dt_s = (log(2))/(DoublingData[(i),4])
  dt_v = append(dt_v, dt_s) 
}

DoublingData <- add_column(DoublingData, Double_Time = dt_v, .after = 3)
DoublingData %>% arrange(Dose) -> DoublingData

### Filter out crazy data based on individual plots

########### Boxplot of the doubling times ############
## The ouput file name must be updated. 
#####################################################
growth_plot <- ggplot(data=DoublingData, aes(x=(Dose), Double_Time)) +
  geom_boxplot() +
  labs(x="Stress", y="Doubling Time (min)") +
  labs(title = "Effects of Dsr2 levels on Cell Divison \n Post Ionizng Radiation", size = 20) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x=element_text(angle=45, vjust =1, hjust = 1, size = 15))+
  theme(axis.text.y=element_text( size = 15)) +
  ylim(5,200)
growth_plot

ggsave(plot = growth_plot , file = "Example_Doubling_Times.pdf", device= "pdf", width=6, height=6, units="in")
################################################################################################################
################     Simple Linear         ##################
## This test will also produce the lab time, the draw back is that it finds the fastest rate
##################################################################
Long_data %>% filter(Time <= 400) -> more100
more100 %>% filter(OD600 >= 0) -> more100

easylinear_data <- all_easylinear(OD600 ~ Time | Dose + Well, data = more100)
plot(easylinear_data)
coef(easylinear_data)
growthrates <-as.data.frame(coef(easylinear_data))

colnames(growthrates) <-c("y0","Max_Growth_Rate") ##name columns

growthrates<-rownames_to_column(growthrates, var = "Dose") ##make names a column

growthrates %>% separate(Dose, c("Dose","Well")) -> DoublingData

### Filter out an strains that didn't grow

DoublingData %>% filter(Dose != "KD10kGy" ) %>% filter(Dose != "KDsham" ) -> DoublingData


### Doubling time
dt_v = vector('numeric')
len = nrow(DoublingData)
for (i in 1:(len)){
  dt_s = (log(2))/(DoublingData[(i),4])
  dt_v = append(dt_v, dt_s) 
}

DoublingData <- add_column(DoublingData, Double_Time = dt_v, .after = 3)
DoublingData %>% arrange(Dose) -> DoublingData

### Filter out crazy data based on individual plots

########### Boxplot of the doubling times ############
## The ouput file name must be updated. 
#####################################################
growth_plot <- ggplot(data=DoublingData, aes(x=(Dose), Double_Time)) +
  geom_boxplot() +
  labs(x="Stress", y="Doubling Time (min)") +
  labs(title = "Effects of Dsr2 levels on Cell Divison \n Post Ionizng Radiation", size = 20) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x=element_text(angle=45, vjust =1, hjust = 1, size = 15))+
  theme(axis.text.y=element_text( size = 15)) +
  ylim(5,200)
growth_plot
## <Philip Sweet> philipjsweet@utexas.edu

