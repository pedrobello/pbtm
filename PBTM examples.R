
install.packages("devtools")
install.packages("roxygen2")
install.packages("usethis")
install.packages("kableExtra")
install.packages("plotly")

## Or, install from GitHub
devtools::install_github("pedrobello/PBTM")


library("plotly")
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
mydata <- read.csv("PBTM-R-Dt-Priming.csv", header=T)
mydata <- read.csv("PBTM-Datasets-R-AfricanR.csv", header=T)
mydata <- read.csv("PBTM-Datasets-R-Q2Paper.csv", header=T)
mydata <- read.csv("PBTM-Datasets-R-Chicory.csv", header=T)

PrimingDt <- read.csv("PBTM-R-Dt-Priming.csv", header=T)
TomatoQ2Dt <- read.csv("PBTM-R-Dt-Q2Paper.csv", header=T)
DatasetDesc <- read.csv("DatasetDescription.csv", header=T)
Models <- read.csv("Models.csv", header=T)

#Save mydata to .RData files
save(mydata, file="data/mydata.RData")

save(PrimingDt, file="PrimingDt.RData")
save(DatasetDesc, file="DatasetDesc.RData")
save(TomatoQ2Dt, file="TomatoQ2Dt.RData")
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

TreatData <- subset(mydata, Treat.priming.temp == 20 | Treat.priming.temp == 0)

TreatData <- subset(mydata, Germ.temp>8 & Germ.temp<24)
TreatData <- subset(mydata, Germ.temp>8 & Germ.temp<24 & Germ.temp!=16 & Germ.temp!=16)
TreatData <- subset(mydata, Germ.temp>7 & Germ.temp<24 & Germ.temp!=16 & Germ.wp==0)

TreatData <- subset(TomatoQ2, Treat.ID == 1 & Germ.wp==0)

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
DefineModel(6.2)

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

PlotPrimingMatrix()


TreatData <- PrimingDt

CalcT50nGR50()

TreatsPriming <<- Treatments %>% group_by(Treat.priming.wp, Treat.priming.temp,Treat.priming.duration, GR50) %>% tally()
TreatsPriming <<- TreatsPriming[order(TreatsPriming$Treat.priming.temp,TreatsPriming$Treat.priming.wp,TreatsPriming$Treat.priming.duration),]
TreatsPriming <- subset(TreatsPriming, Treat.priming.wp < 0)

PrTemps <- unique(TreatsPriming$Treat.priming.temp)
PrTempAmt <- length(unique(TreatsPriming$Treat.priming.temp))
PrWPs <- unique(TreatsPriming$Treat.priming.wp)
PrWPAmt <- length(unique(TreatsPriming$Treat.priming.wp))
PrDurs <- unique(TreatsPriming$Treat.priming.duration)
PrDurAmt <- length(unique(TreatsPriming$Treat.priming.duration))

TreatsPrTemp1 <- subset(TreatsPriming, Treat.priming.temp == unique(TreatsPriming$Treat.priming.temp)[1])
TreatsPrTemp2 <- subset(TreatsPriming, Treat.priming.temp == unique(TreatsPriming$Treat.priming.temp)[2])
TreatsPrTemp3 <- subset(TreatsPriming, Treat.priming.temp == unique(TreatsPriming$Treat.priming.temp)[3])
TreatsPrTemp4 <- subset(TreatsPriming, Treat.priming.temp == unique(TreatsPriming$Treat.priming.temp)[4])
TreatsPrTemp5 <- subset(TreatsPriming, Treat.priming.temp == unique(TreatsPriming$Treat.priming.temp)[5])

Temp1 <- matrix(TreatsPrTemp1$GR50, nrow = PrDurAmt, ncol = PrWPAmt)
Temp2 <- matrix(TreatsPrTemp2$GR50, nrow = PrDurAmt, ncol = PrWPAmt)
Temp3 <- matrix(TreatsPrTemp3$GR50, nrow = PrDurAmt, ncol = PrWPAmt)
Temp4 <- matrix(TreatsPrTemp4$GR50, nrow = PrDurAmt, ncol = PrWPAmt)
Temp5 <- matrix(TreatsPrTemp5$GR50, nrow = PrDurAmt, ncol = PrWPAmt)


SurfPlotHTP <- plot_ly(showscale = FALSE ) %>%
  add_surface(x = PrWPs , y = PrDurs, z = Temp1, opacity = 1, colorscale = list(c(0,1),c("rgb(19,50,117)","rgb(0,183,255)")), name = paste0(unique(TreatsPriming$Treat.priming.temp)[1], "C")) %>%
  add_surface(x = PrWPs , y = PrDurs, z = Temp2, opacity = 1, colorscale = list(c(0,1),c("rgb(23,46,17)","rgb(45,209,0)")), name = paste0(unique(TreatsPriming$Treat.priming.temp)[2], "C")) %>%
  add_surface(x = PrWPs , y = PrDurs, z = Temp3, opacity = 1, colorscale = list(c(0,1),c("rgb(232,116,49)","rgb(191,204,43)")), name = paste0(unique(TreatsPriming$Treat.priming.temp)[3], "C")) %>%
  add_surface(x = PrWPs , y = PrDurs, z = Temp4, opacity = 1, colorscale = list(c(0,1),c("rgb(23,46,17)","rgb(45,209,0)")), name = paste0(unique(TreatsPriming$Treat.priming.temp)[4], "C")) %>%
  add_surface(x = PrWPs , y = PrDurs, z = Temp5, opacity = 1, colorscale = list(c(0,1),c("rgb(232,116,49)","rgb(191,204,43)")), name = paste0(unique(TreatsPriming$Treat.priming.temp)[5], "C")) %>%

  layout(
    ## title = "Layout options in a 3d scatter plot",
    scene = list(
      surfacecolor = "rgb(244, 244, 248)",
      xaxis = list(title = "WP"),
      yaxis = list(title = "Duration"),
      zaxis = list(title = "GR50"),
      camera = list(eye = list(x = -1.85, y = 1.95, z = 0.75), center = list(x = 0, y = 0, z = 0), up = list(x = 0, y = 0, z = 1)),
      type = "perspective"
    ))




