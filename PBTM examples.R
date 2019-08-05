
install.packages("devtools")
install.packages("roxygen2")
install.packages("usethis")
install.packages("kableExtra")

devtools::install_github("pedrobello/PBTM")

library("devtools")
library("roxygen2")
library("usethis")
library("kableExtra")

library("PBTM")


devtools::document()

#Population-based Threshold Models Calculator ------------------------------------------------------
#Set working folder - Add the folder location of your data files

setwd("~/Dropbox/01-UCDavis/PBTM/PBTM-RPackage/Data")
setwd("~/Dropbox/01-UCDavis/PBTM/PBTM-RPackage")
setwd("~/Dropbox/01-UCDavis/PBTM")

#Load Seeds data - Add the name of the data file (csv file)
mydata <- read.csv("PBTM-Datasets-R.csv", header=T)
mydata <- read.csv("PBTM-Datasets-R-Empty.csv", header=T)
mydata <- read.csv("PBTM-Datasets-Priming.csv", header=T)
mydata <- read.csv("PBTM-Datasets-R-AfricanR.csv", header=T)
mydata <- read.csv("PBTM-Datasets-R-Q2Paper.csv", header=T)
mydata <- read.csv("PBTM-Datasets-R-Chicory.csv", header=T)

TomatoQ2 <- read.csv("PBTM-R-Dt-Q2Paper.csv", header=T)
DatasetDesc <- read.csv("DatasetDescription.csv", header=T)
Models <- read.csv("Models.csv", header=T)

#Save mydata to .RData files
save(mydata, file="data/mydata.RData")

save(DatasetDesc, file="DatasetDesc.RData")
save(TomatoQ2, file="TomatoQ2.RData")
save(Models, file="Models.RData")


#Filter the data to be used on models
#All data inside the table TreatData will be used on the models
TreatData <- subset(mydata, Seed.lot=="26" & Germ.temp=="25")
TreatData <- subset(mydata, Seed.lot=="26" & Germ.temp=="25" & Germ.time.hours<25)
TreatData <- subset(mydata, Seed.lot=="27" & Germ.wp=="0" & Germ.temp<26)
TreatData <- subset(mydata, Treat.ID < 40 )
TreatData <- subset(mydata, Germ.temp > 24)
TreatData <- subset(mydata, Germ.temp > 19)
TreatData <- subset(mydata, Seed.lot > 0 & Treat.priming.temp == 20)
TreatData <- subset(mydata, Seed.lot > 0)
TreatData <- subset(mydata, Treat.ID == 1 & Germ.wp==0)
TreatData <- subset(mydata, Seed.lot == 1 & Germ.temp<24 & Germ.wp==0)
TreatData <- subset(mydata, Seed.lot == 1 & Germ.temp>19 & Germ.wp==0)
TreatData <- subset(mydata, Seed.lot == 1 & Germ.temp<21 & Germ.temp>10)
TreatData <- subset(mydata, Seed.lot == 1 & Germ.temp<24 & Germ.wp>-0.8 & Germ.wp==0)
TreatData <- subset(mydata, Germ.temp==40)
TreatData <- subset(mydata[ which( mydata$Germ.wp == 0 | mydata$Germ.wp < -0.8) , ], Seed.lot == 1 & Germ.temp==12)

TreatData <- subset(mydata, Germ.temp>8 & Germ.temp<24)
TreatData <- subset(mydata, Germ.temp>8 & Germ.temp<24 & Germ.temp!=16 & Germ.temp!=16)
TreatData <- subset(mydata, Germ.temp>7 & Germ.temp<24 & Germ.temp!=16 & Germ.wp==0)

Data(TomatoQ2)

TreatData <- mydata
TreatData <- TomatoQ2

#Calculates the time to 50% germination (T50) and the germination rate (GR50).
#Function separates the treatments and adds data to the Treatments table.
#Functions uses TreatData as data source.
CalcT50nGR50()

#Plot germination rates over temperature
#After CalcT50nGR50()
PlotGR50vsTemp()

#Choose the model that you want to work on
# (1)Hydropriming model; (2)Suboptimal Hydrotime; (2.5)Supra-optimal Hydrotime; (3)Thermaltime;
# (4)Suboptimal Hydrothermal Time; (5)Supra-optimal Hydrothermal Time
DefineModel(1.1)

#Define maximum germination percentage when seed lot does not germinate at optimum temperature and water potential.
#Identify the correspondent treatment on the function. Use this function at your own discretion.
DefineMaxGerm(SeedLot=1,GermWP=0,GermTEMP=20)
# or use standard maximum germination as 100%
MaxGerm <- 1

#Plot raw data to be used on model
PlotRawData()

#Clean Repetitive Percentages for models(keeps only initial germination observation and remove extra points without increase)
#Adds processed data to Table TreatDataClean
CleanGermData()

#Calculate model parameters and plot raw and predicted data
CALCnPLOTModel()

#Plot model predicted and cleaned raw data
PlotCleanModel()

DefineHTo(9)
FixedTb(4.8)
FixedTo(25)

#Plot model predicted vs raw data
PlotActualvsModelPredicted()

PlotActualvsModelPredictedTreatments()

#Plot normalized data
PlotNormalizedTime()
