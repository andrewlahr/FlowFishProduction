

##TOP####
rm(list=ls())
library(paletteer)
library(coda)
library(plyr)
library(rjags)
library(R2jags)
library(here)
library(ggrepel)
library(glmmTMB)
library(MuMIn)
library(broom)
library(ggcorrplot)
library(gridExtra)
library(gridGraphics)
library(grid)
library(runjags)
library(colorspace)
library(tidyverse)
library(dplyr)

Zscore<-function(x,demean.only=F){
  if(demean.only){
    return((x-mean(x,na.rm=T)))
  }else{
    return((x-mean(x,na.rm=T))/sd(x,na.rm=T))
  }
}


#Read in ennvirontmental Data
load(here::here('../LL/Data/ComputeFlowAndTemps071725.Rdata'))

#Set years to filter abundance and temerpature data
YearSet <- seq(1980,2024)

#Set up directory
dir_name <- 'PVA_ProductionModel'
if(!dir.exists(here::here('JAGS_PVA'))){dir.create(here::here('JAGS_PVA'))}
if(!dir.exists(here::here('JAGS_PVA','JAGSmodels'))){dir.create(here::here('JAGS_PVA','JAGSmodels'))}

if(!dir.exists(here::here('JAGS_PVA','ModelFits'))){dir.create(here::here('JAGS_PVA','ModelFits'))}
if(!dir.exists(here::here('JAGS_PVA','ModelOutput'))){dir.create(here::here('JAGS_PVA','ModelOutput'))}
if(!dir.exists(here::here('Data'))){dir.create(here::here('Data'))}

#Filter out sites where there are greater than 10 years missing from ennvirontmental data
WinterMeanFlowTable_Good<- WinterMeanFlowTable[which(rowSums(!is.na(WinterMeanFlowTable))>10),]
SummerMeanFlowTable_Good <- SummerMeanFlowTable[which(rowSums(!is.na(SummerMeanFlowTable))>10),]
SummerMeanTempTable_Good<-SummerMeanTempTable[which(rowSums(!is.na(SummerMeanTempTable))>10),]

#Read in abundance estimates,
PopEstsLast20 <- read.csv('Data/AllRBTrout_AllYears_PopulationEstimates_wDet_071825.csv') %>% 
  filter(Year %in% YearSet)%>%
  dplyr::mutate(Stream.Section=paste(Stream,Section,sep="."),age=as.numeric(as.character(Age)))%>%
  dplyr::mutate(SampleSeason = ifelse(SurveyStartDOY<182,'Spring','Fall'),SampleSeasonCode = ifelse(SurveyStartDOY<182,1,2))%>%
  arrange(Stream.Section,Year,Season, Year,age)%>%dplyr::mutate(age=age+2)%>%filter(!is.na(Season))%>%
  mutate(Season=if_else(Section=="JerryCreek",'Spring',Season))%>% ##Jerry creek was sampled in the fall and spring in multiple years but only select spring samples to avoid double counting
  mutate(LCL=if_else(LCL<0 ,1,LCL)) #if abundance estimates have a negative LCL make it 1

#Seperate out pop estimates based on ages
EST_All <- PopEstsLast20%>%filter(Year %in% YearSet) %>%xtabs(Est ~ Stream.Section + Year + Season + age  ,data=.)
EST_All[EST_All==0] <- NA

EST_All_SE <- PopEstsLast20 %>% filter(Year %in% YearSet) %>% xtabs(SE ~ Stream.Section + Year+ Season + age ,data=.)
EST_All_lcl <- PopEstsLast20 %>% filter(Year %in% YearSet) %>% xtabs(LCL ~ Stream.Section + Year+ Season + age ,data=.)
EST_All_ucl <- PopEstsLast20 %>% filter(Year %in% YearSet) %>% xtabs(UCL ~ Stream.Section + Year+ Season + age ,data=.)

AGE_All <- PopEstsLast20%>%filter(Year %in% YearSet) %>%  xtabs(Est ~ Stream.Section + Year + Season  ,data=.)
AGE_All[AGE_All==0] <- NA

AGE_All_SE <- PopEstsLast20 %>% filter(Year %in% YearSet) %>% xtabs(SE ~ Stream.Section + Year+ Season  ,data=.)
AGE_All_lcl <- PopEstsLast20 %>% filter(Year %in% YearSet) %>% xtabs(LCL ~ Stream.Section + Year+ Season ,data=.)
AGE_All_ucl <- PopEstsLast20 %>% filter(Year %in% YearSet) %>% xtabs(UCL ~ Stream.Section + Year+ Season ,data=.)



#Read in raw survey data
BT_extract<-read.csv("output/RBT_extractIN_071825.csv",header=T,stringsAsFactors = F)%>%
  dplyr::mutate(Stream.Section=paste(Water,Section,sep='.'))%>%
  dplyr::mutate(age=age)

#Set up to create length weight relationship to estimate biomas from abundance and weight data
#Get all weight data and transform to log length and weight so that we can get length x weight relationships for each section
WeightsIN<-BT_extract%>%
  filter(Weight!=0 & Length!=0 & Weight<3000)%>%
  dplyr::mutate(logLength=log(Length),
                logWeight=log(Weight))%>%
  filter(!is.na(logWeight) & !is.na(Year) )%>%dplyr::mutate(Year=as.factor(as.character(Year)))%>%ungroup()

#calculate lw relationship
lengthWeight<- glmmTMB::glmmTMB(formula=logWeight~logLength+(Stream.Section|Year),data=WeightsIN)
summary(lengthWeight)
saveRDS(lengthWeight,"Data/lengthWeight071825.rds")
lengthWeight<-readRDS("Data/lengthWeight071825.rds")
newData<-BT_extract%>%filter(Weight==0 | is.na(Weight) | Weight>3000)%>%
  dplyr::mutate(logLength=log(Length))%>%
  dplyr::mutate(Stream.Section=paste(Water,Section,sep='.'))%>%
  filter(Stream.Section!=".Casca")

#Predict missing weights
newData$Weight<-exp(predict(readRDS("Data/lengthWeight071825.rds"),newdata=newData))

#Get average weight by age for each section in each year
WeightsOUT<-rbind(WeightsIN%>%dplyr::select(-logWeight),newData)%>%
  dplyr::mutate(Stream=Water,
                Age=age)%>%
  filter(Year %in% YearSet) %>%
  dplyr::mutate(Stream.Section=paste(Stream,Section,sep="."),age=as.numeric(as.character(Age)))%>%
  dplyr::mutate(SampleSeason = ifelse(jday<182,'Spring','Fall'),SampleSeasonCode = ifelse(jday<182,1,2))%>%
  dplyr::mutate(age=age)%>%
  group_by(Year,Stream.Section,age,SampleSeason)%>%
  dplyr::rename(Season=SampleSeason)%>%
  dplyr::summarize(avgWeight=mean(Weight,na.rm = T),
                   sdWeight=sd(Weight,na.rm=T),
                   nWeight=n(),
                   seWeight=sdWeight/sqrt(nWeight),
                   avgLength=mean(Length,na.rm = T),
                   sdLength=sd(Length,na.rm=T),
                   nLength=n(),
                   seLength=sdLength/sqrt(nLength))%>%
  ungroup()


#Create tables for average weights
WEIGHT <- WeightsOUT %>%xtabs(avgWeight ~ Stream.Section + Year + Season + age ,data=.)
WEIGHT_SE <- WeightsOUT%>%filter(Year %in% YearSet) %>% xtabs(seWeight ~ Stream.Section + Year + Season + age ,data=.)
WEIGHT_SE[WEIGHT_SE==0|WEIGHT_SE==1] <- 0

WEIGHT_2 <- WeightsOUT %>% filter(age==2)%>%filter(Year %in% YearSet) %>% xtabs(avgWeight ~ Stream.Section + Year + Season ,data=.)
WEIGHT_2[WEIGHT_2==0|WEIGHT_2==1] <- NA
WEIGHT_2_SE <- WeightsOUT%>% filter(age==2)%>%filter(Year %in% YearSet) %>% xtabs(seWeight ~ Stream.Section + Year + Season ,data=.)
WEIGHT_2_SE[WEIGHT_2_SE==0|WEIGHT_2_SE==1] <- 0

WEIGHT_3 <- WeightsOUT %>% filter(age==3)%>%filter(Year %in% YearSet) %>% xtabs(avgWeight ~ Stream.Section + Year + Season ,data=.)
WEIGHT_3[WEIGHT_3==0|WEIGHT_3==1] <- NA
WEIGHT_3_SE <- WeightsOUT%>% filter(age==3)%>%filter(Year %in% YearSet) %>% xtabs(seWeight ~ Stream.Section + Year + Season ,data=.)
WEIGHT_3_SE[WEIGHT_3_SE==0|WEIGHT_3_SE==1] <- 0

WEIGHT_4 <- WeightsOUT %>% filter(age==4)%>%filter(Year %in% YearSet) %>% xtabs(avgWeight ~ Stream.Section + Year + Season ,data=.)
WEIGHT_4[WEIGHT_4==0|WEIGHT_4==1] <- NA
WEIGHT_4_SE <- WeightsOUT%>% filter(age==4)%>%filter(Year %in% YearSet) %>% xtabs(seWeight ~ Stream.Section + Year + Season ,data=.)
WEIGHT_4_SE[WEIGHT_4_SE==0|WEIGHT_4_SE==1] <- 0


LENGTH <- WeightsOUT %>%filter(Year %in% YearSet) %>% xtabs(avgLength ~ Stream.Section + Year + Season + age ,data=.)
LENGTH_SE <- WeightsOUT%>%filter(Year %in% YearSet) %>% xtabs(seLength ~ Stream.Section + Year + Season + age ,data=.)
LENGTH_SE[LENGTH_SE==0|LENGTH_SE==1] <- 0

##Calculate BioMass by multiplying average weight per fish by the number of fish estimated

BWeight<-rbind(WeightsIN%>%dplyr::select(-logWeight),newData)%>%
  dplyr::mutate(Stream=Water,
                Age=age)%>%
  filter(Year %in% YearSet) %>%
  dplyr::mutate(Stream.Section=paste(Stream,Section,sep="."),age=as.numeric(as.character(Age)))%>%
  dplyr::mutate(SampleSeason = ifelse(jday<182,'Spring','Fall'),SampleSeasonCode = ifelse(jday<182,1,2))%>%
  dplyr::mutate(age=age)%>%
  group_by(Year,Stream.Section,age,SampleSeason)%>%
  dplyr::rename(Season=SampleSeason)

BioMass<- left_join(PopEstsLast20,WeightsOUT%>%dplyr::mutate(Year=as.numeric(as.character(Year))),by=c("Stream.Section","Year","Season","age"))%>%
  dplyr::mutate(B_Est=avgWeight*Est,
                B_SE=((seWeight/avgWeight)+(SE/Est))*B_Est)%>%
  arrange(Stream.Section,Year,Season,age)%>%
  dplyr::mutate(Year=as.factor(Year))

#Create Biomass tables
B <- BioMass%>%filter(Year %in% YearSet) %>%  xtabs(B_Est ~ Stream.Section + Year + Season + age  ,data=.)
B_SE <- BioMass %>% filter(Year %in% YearSet) %>% xtabs(B_SE ~ Stream.Section + Year+ Season + age  ,data=.)
B_SE[B_SE==0|B_SE==1] <- 0

B_All <- BioMass%>%filter(Year %in% YearSet) %>%  xtabs(B_Est ~ Stream.Section + Year + Season  ,data=.)
B_All[B_All==0] <- NA
B_All_SE <- BioMass %>% filter(Year %in% YearSet) %>% xtabs(B_SE ~ Stream.Section + Year+ Season  ,data=.)


AllSections<-PopEstsLast20%>%pull(Stream.Section)%>%unique()

##List of Sectionswith bad abundance time series or omitted from analysis to be filtered from dataset
badSections<-c("Ruby.Greenhorn","BigHole.MaidenRock","Ruby.Seyler","Beaverhead.Anderson","ClarkFork.Milltown","ClarkFork.Turah","Shields.Zimmerman","Bitterroot.Hannon","Bitterroot.Hamilton","Bitterroot.Stevensville",
               "ClarkFork.Bearmouth","ClarkFork.BelowSagerLane",'Ruby.SilverSprings',"Ruby.Vigilante" ,"ClarkFork.MorseRanch","ClarkFork.Phosphate","ClarkFork.Phshack","ClarkFork.WilliamsTavenner","Beaverhead.FishAndGame" )

dir_name2<-("JAGS_PVA/ModelOutput")

AllSectionsIN<-AllSections[which(!AllSections%in%badSections)]


modelScript.Global <- "PVA_ProductionModel_Global.txt"
jagsscript <- cat("
model {

    for(i in 1:RecLag){

    ###Start initial abundance and weight values ####
     NALL[i]<-N4[i]+N2[i]+N3[i];

      N4[i]  ~ dgamma(pow(maxN4,2)/pow(SDmaxN4,2),maxN4/pow(SDmaxN4,2));
      N3[i]  ~ dgamma(pow(maxN3,2)/pow(SDmaxN3,2),maxN3/pow(SDmaxN3,2));
      R[i] ~ dgamma(pow(maxNstwo1,2)/pow(SDmaxNstwo1,2),maxNstwo1/pow(SDmaxNstwo1,2));
      N2[i]<-R[i];
      
      BALL[i]<-B4[i]+B3[i]+B2[i];
      B4[i]  ~ dgamma(pow(maxB4,2)/pow(SDmaxB4,2),maxB4/pow(SDmaxB4,2));
      B3[i]  ~ dgamma(pow(maxB3,2)/pow(SDmaxB3,2),maxB3/pow(SDmaxB3,2));
      BR[i] ~ dgamma(pow(maxBstwo1,2)/pow(SDmaxBstwo1,2),maxBstwo1/pow(SDmaxBstwo1,2));
      B2[i]<-BR[i];
      BAdults[i]<-B4[i]+(B3[i]);

      weight4[i]~dnorm(W_4[i],pow(W_4_SE[i],-2));
      weight3[i]~dnorm(W_3[i],pow(W_3_SE[i],-2));
      weight2[i]~dnorm(W_2[i],pow(W_2_SE[i],-2));

    }


    for(t in (RecLag+1):(nYears)){
      NALL[t]<-N4[t]+R[t]+N3[t]
      BALL[t]<-B4[t]+BR[t]+B3[t];
      NAdults[t]<-N4[t]+(N3[t]);
      BAdults[t]<-B4[t]+(B3[t]);

      N4[t]~dpois((Surv[t]*N4[t-1])+(Surv[t]*N3[t-1])); 
      B4[t]<-N4[t]*weight4[t];

      N3[t]~dpois(N2[t-1]*Surv[t]);
      B3[t]<-N3[t]*weight3[t];

      BR[t]<-R[t]*weight2[t];
      N2[t]<-R[t];

      R[t] <-exp(lR[t]);  #Recruits
      lR[t] ~ dnorm(la[t] + log(BAdults[t-RecLag]) - (b*BAdults[t-RecLag]) + (b1r*covarlag1[t-SummerLags[lagDum1]]+((b12r*pow(covarlag1[t-SummerLags[lagDum1]],2))*lagQuad1)) + (b5r*covarlag5[t-WinterLags[lagDum2]]+((b52r*pow(covarlag5[t-WinterLags[lagDum2]],2))*lagQuad2)) - (pow(sigmaR,2)/2),pow(sigmaR,-2)); #Ricker

      Surv[t] <- exp(lSurv[t])/(1+exp(lSurv[t]));
      lSurv[t] ~ dnorm(lSurv0 + (b1S*covarlag1[t-SurvLag])+ ((b12S*pow(covarlag1[t-SurvLag],2))*lagQuad3) + (b5S*covarlag5[t-(SurvLag-winterlag)])+ ((b52S*pow(covarlag5[t-(SurvLag-winterlag)],2))*lagQuad4)  +(b6S*(NALL[t-1]/100)) ,pow(sigmaSurv,-2)); #log survival with covariates


        la[t]<- la0; #time varying alpha parameter for recruitment

        weight4[t]~dnorm(W_4[t],pow(W_4_SE[t],-2));
        weight3[t]~dnorm(W_3[t],pow(W_3_SE[t],-2));
        weight2[t]~dnorm(W_2[t],pow(W_2_SE[t],-2));
    }
    
  ####Priors###
    b1S ~ dnorm(0,pow(3,-2));
    b12S  ~ dt(0, pow(0.125,-2), 1) T(,0);
    
    b5S ~ dnorm(0,pow(3,-2));
    b52S  ~ dt(0, pow(0.125,-2), 1) T(,0);

    b6S ~ dnorm(0,pow(3,-2));

    b1r~ dnorm(0,pow(3,-2));
    b12r ~ dt(0, pow(0.125,-2), 1) T(,0);

    b5r~ dnorm(0,pow(3,-2));
    b52r ~ dt(0, pow(0.125,-2), 1) T(,0);

    Q ~ dunif(0.0001,1);
    lSurv0 ~ dnorm(0,pow(3,-2));
    la0 ~ dunif(-10,10);
    b ~ dunif(0,0.0001);

    sigmaR ~ dt(0, pow(sigmaRmax,-2), 1) T(0.0001,);
    sigmaSurv ~ dt(0, pow(1,-2), 1) T(0.0001,);


    lagDum1~dcat(pi1)
    lagDum2~dcat(pi2)

    lagQuad2~ dbern(0.5);
    lagQuad1~ dbern(0.5);
    lagQuad3~ dbern(0.5);
    lagQuad4~ dbern(0.5);
  
  ####likelihood####
    for(t in 1:(nYears)){
       YALL[t]~ dnorm(NALL[t],pow(YALL_SE[t],-2));
        Y4[t] ~ dnorm(N4[t],pow(Y4_SE[t],-2));
        Y3[t] ~ dnorm(N3[t],pow(Y3_SE[t],-2));
        Y2[t] ~ dnorm(R[t]*Q,pow(Y2_SE[t],-2));
        logLik[t]<- logdensity.norm(Y2[t],R[t]*Q,pow(Y2_SE[t],-2))+logdensity.norm(Y3[t],N3[t],pow(Y3_SE[t],-2))+logdensity.norm(Y4[t],N4[t],pow(Y4_SE[t],-2))

    }
}
",
file = here::here('JAGS_PVA','JAGSmodels',modelScript.Global))


lapply(AllSectionsIN,FUN=function(SS){
  
  ## build model data####
 
  if(SS=="BigHole.JerryCreek"){
    Season<-PopEstsLast20%>%filter(Stream.Section==SS & Season=="Spring")%>%pull(Season)%>%first()
    
  }else{
    Season<- BT_extract%>%
      dplyr::mutate(Stream.Section=paste(Water,Section,sep='.'))%>%
      filter(Stream.Section==SS)%>%pull(selectSeason)%>%first()
  }
  print(SS)
  print(Season)

  Stream<-strsplit(SS,split = '\\.')[[1]][1]
  Section<-strsplit(SS,split = '\\.')[[1]][2]
  
  dir_name2<-("JAGS_PVA/ModelOutput/")
  
  if(!dir.exists(here::here("JAGS_PVA/ModelOutput/",Stream,Section))){dir.create(here::here("JAGS_PVA/ModelOutput/",Stream,Section))}
  if(!dir.exists(here::here("plots_WinterAvg",Stream))){dir.create(here::here("plots_WinterAvg",Stream))}
  if(!dir.exists(here::here("plots_WinterAvg",Stream,Section))){dir.create(here::here("plots_WinterAvg",Stream,Section))}
  
  ##### Build data for model input#####
  
  Cov1<-SummerMeanFlowTable_Good
  Cov5<-WinterMeanFlowTable_Good
  
  if(SS=="RockCreek.HogBack"){
    coVar1 <- Cov1['RockCreek.Hogback',] %>% log() %>% Zscore()
    covarTable_Real1<-Cov1['RockCreek.Hogback',]
  }else{
    if(SS=="Shields.Zimmerman"){
      coVar1 <- Cov1['Shields.ConvictGrade',] %>% log() %>% Zscore()
      covarTable_Real1<-Cov1['Shields.ConvictGrade',]
    }else{
      coVar1 <- Cov1[SS,] %>% log() %>% Zscore()
      covarTable_Real1<-Cov1[SS,]}}
  if(SS=="RockCreek.HogBack"){
    coVar5 <- Cov5['RockCreek.Hogback',] %>% log() %>% Zscore()
    covarTable_Real5<-Cov5['RockCreek.Hogback',]
  }else{if(SS=="Shields.Zimmerman"){
    coVar5 <- Cov5['Shields.ConvictGrade',] %>% log() %>% Zscore()
    covarTable_Real5<-Cov5['Shields.ConvictGrade',]
  }else{
    coVar5 <- Cov5[SS,] %>% log() %>% Zscore()
    covarTable_Real5<-Cov5[SS,]}
  }
  
  #Pop est data
  Y4_In<-EST_All[SS,,Season,"4"]
  Y4_lcl_In<-EST_All_lcl[SS,,Season,"4"]
  Y4_ucl_In<-EST_All_ucl[SS,,Season,"4"]
  Y4_SE_In<-EST_All_SE[SS,,Season,"4"]
  Y4_In[Y4_In==0]<-NA
  
  Y3_In<-EST_All[SS,,Season,"3"]
  Y3_lcl_In<-EST_All_lcl[SS,,Season,"3"]
  Y3_ucl_In<-EST_All_ucl[SS,,Season,"3"]
  Y3_SE_In<-EST_All_SE[SS,,Season,"3"]
  Y3_In[Y3_In==0]<-NA
  
  Y2_In<-EST_All[SS,,Season,"2"]
  Y2_lcl_In<-EST_All_lcl[SS,,Season,"2"]
  Y2_ucl_In<-EST_All_ucl[SS,,Season,"2"]
  Y2_SE_In<-EST_All_SE[SS,,Season,"2"]
  Y2_In[Y2_In==0]<-NA
  
  YALL_In<-AGE_All[SS,,Season]
  YALL_SE_In<-AGE_All_SE[SS,,Season]
  YALL_lcl_In<-AGE_All_lcl[SS,,Season]
  YALL_ucl_In<-AGE_All_ucl[SS,,Season]
  YALL_In[YALL_In==0]<-NA
  
  #Biomass and Weight data
  M4_In<-B[SS,,Season,"4"]
  M4_SE_In<-B_SE[SS,,Season,"4"]
  M4_In[M4_In==0]<-NA
  
  W4_In<-WEIGHT[SS,,Season,"4"]
  W4_SE_In<-WEIGHT_SE[SS,,Season,"4"]
  
  M3_In<-B[SS,,Season,"3"]
  M3_SE_In<-B_SE[SS,,Season,"3"]
  M3_In[M3_In==0]<-NA
  
  W3_In<-WEIGHT[SS,,Season,"3"]
  W3_SE_In<-WEIGHT_SE[SS,,Season,"3"]
  
  M2_In<-B[SS,,Season,"2"]
  M2_SE_In<-B_SE[SS,,Season,"2"]
  M2_In[M2_In==0]<-NA
  
  W2_In<-WEIGHT[SS,,Season,"2"]
  W2_SE_In<-WEIGHT_SE[SS,,Season,"2"]
  
  MALL_In<-B_All[SS,,Season]
  MALL_SE_In<-B_All_SE[SS,,Season]
  MALL_In[MALL_In==0]<-NA
  
  
  delay2<-sapply(coVar1,FUN=function(x){
    sum(is.na(x))
    
  })
  
  covarLagIn1Real <-covarTable_Real1[which(delay2==0)[1]:which(delay2==0)[length(which(delay2==0))]]
  
  IncludedYears <- names(covarLagIn1Real)
  NumericYears<-as.numeric(IncludedYears)
 
  covarLagIn1<-covarLagIn1Real%>%log()%>%Zscore()#}
  
  MALL_In<-MALL_In[which(names(MALL_In) %in% c(IncludedYears,max(as.numeric(names(MALL_In)))))]
  
  delay1<-sapply(MALL_In,FUN=function(x){
    sum(is.na(x))
    
  })
  if(mean(delay1)==0){
    MALL_In<-MALL_In}else{
      
      if(SS=="Ruby.Greenhorn"|SS=="RockCreek.FishAndGame"){
        MALL_In <- MALL_In[((which(delay1<max(delay1))[2])+0):which(delay1<max(delay1))[length(which(delay1<max(delay1)))]]
        
      }else{
        MALL_In <- MALL_In[((which(delay1<max(delay1))[1])+0):which(delay1<max(delay1))[length(which(delay1<max(delay1)))]]}
    }
  #
  if(max(as.numeric(names(MALL_In)))>=2024){
    IncludedYears <- c(min(as.numeric(names(MALL_In))):2023)
  }else{
    IncludedYears <- (as.numeric(names(MALL_In)))}
  
  NumericYears<- as.numeric(IncludedYears)
  if(SS=='Bighorn.Bighorn'){
    NumericYears<-1995:max(NumericYears)
  }else{NumericYears<-NumericYears}
  
  if(SS=='BigHole.Melrose'){
    NumericYears<-1991:max(NumericYears)
  }else{NumericYears<-NumericYears}
  
  if(SS=='RockCreek.FishAndGame'){
    NumericYears<-2004:max(NumericYears)
  }else{NumericYears<-NumericYears}
  if(SS=='RockCreek.HogBack'){
    NumericYears<-2004:max(NumericYears)
  }else{NumericYears<-NumericYears}
  if(SS=='Madison.Varney'){
    NumericYears<-1995:max(NumericYears)
  }else{NumericYears<-NumericYears}
  CovYears<-as.character(c((min(NumericYears)):(max(as.numeric(names(covarTable_Real1))))))
  covarLagIn1Real<-covarTable_Real1[which(names(covarTable_Real1) %in% as.numeric(CovYears))]
  covarLagIn1<-covarTable_Real1[which(names(covarTable_Real1) %in% c(as.numeric(CovYears)))]%>%log()%>%Zscore()
  
  covarLagIn5Real<-covarTable_Real5[which(names(covarTable_Real5) %in% c((min(as.numeric(CovYears))):(max(as.numeric(CovYears)))))]
  covarLagIn5<-covarLagIn5Real[which(names(covarLagIn5Real) %in% c((min(as.numeric(CovYears))):(max(as.numeric(CovYears)))))]%>%log()%>%Zscore()
  
  RecLag <- 3
  
  winterlag <- ifelse(Season=="Fall",0,1)
  SurvLag <- ifelse(Season=="Fall",0,1)
  
  
  dataIN<-data.frame(c(covarLagIn1Real[1:(length(covarLagIn1Real)-1)]),covarLagIn5Real[2:(length(covarLagIn5Real))]) #starts at winter 1992 and summer 1991 because summer 91 is precurser to winter 92 condtions
  colnames(dataIN)<-c("SummerMeanQ" ,"WinterMeanQ")
  cors<-cor(dataIN,use="pairwise.complete.obs")
  
  Y4_In<-Y4_In[which(names(Y4_In) %in% (NumericYears))]
  Y4_In[Y4_In==0]<-NA
  Y3_In<-Y3_In[which(names(Y3_In) %in% (NumericYears))]
  Y3_In[Y3_In==0]<-NA
  Y2_In<-Y2_In[which(names(Y2_In) %in% (NumericYears))]
  Y2_In[Y2_In==0]<-NA
  YALL_In<-YALL_In[which(names(YALL_In) %in% (NumericYears))]
  YALL_In[YALL_In==0]<-NA
  
  M4_In<-M4_In[which(names(M4_In) %in% (NumericYears))]
  M4_In[M4_In==0]<-NA
  M3_In<-M3_In[which(names(M3_In) %in% (NumericYears))]
  M3_In[M3_In==0]<-NA
  M2_In<-M2_In[which(names(M2_In) %in% (NumericYears))]
  M2_In[M2_In==0]<-NA
  
  W4_In<-W4_In[which(names(W4_In) %in% (NumericYears))]
  for(i in 2:length(W4_In)){
    W4_In[i]<-if_else(W4_In[i]==0,W4_In[i-1],W4_In[i])
  }
  W3_In<-W3_In[which(names(W3_In) %in% (NumericYears))]
  for(i in 2:length(W3_In)){
    W3_In[i]<-if_else(W3_In[i]==0,W3_In[i-1],W3_In[i])
  }
  W2_In<-W2_In[which(names(W2_In) %in% (NumericYears))]
  for(i in 2:length(W2_In)){
    W2_In[i]<-if_else(W2_In[i]==0,W2_In[i-1],W2_In[i])
  }
  
  
  Y4_SE_In<-Y4_SE_In[which(names(Y4_SE_In) %in% (NumericYears))]
  Y4_lcl_In<-Y4_lcl_In[which(names(Y4_lcl_In) %in% (NumericYears))]
  Y4_ucl_In<-Y4_ucl_In[which(names(Y4_ucl_In) %in% (NumericYears))]
  
  Y3_SE_In<-Y3_SE_In[which(names(Y3_SE_In) %in% (NumericYears))]
  Y3_lcl_In<-Y3_lcl_In[which(names(Y3_lcl_In) %in% (NumericYears))]
  Y3_ucl_In<-Y3_lcl_In[which(names(Y3_lcl_In) %in% (NumericYears))]
  
  Y2_SE_In<-Y2_SE_In[which(names(Y2_SE_In) %in% (NumericYears))]
  Y2_lcl_In<-Y2_lcl_In[which(names(Y2_lcl_In) %in% (NumericYears))]
  Y2_ucl_In<-Y2_ucl_In[which(names(Y2_ucl_In) %in% (NumericYears))]
  
  YALL_SE_In<-YALL_SE_In[which(names(YALL_SE_In) %in% (NumericYears))]
  YALL_lcl_In<-YALL_lcl_In[which(names(YALL_lcl_In) %in% (NumericYears))]
  YALL_ucl_In<-YALL_ucl_In[which(names(YALL_ucl_In) %in% (NumericYears))]
  
  M4_SE_In<-M4_SE_In[which(names(M4_SE_In) %in% (NumericYears))]
  M3_SE_In<-M3_SE_In[which(names(M3_SE_In) %in% (NumericYears))]
  M2_SE_In<-M2_SE_In[which(names(M2_SE_In) %in% (NumericYears))]
  MALL_SE_In<-MALL_SE_In[which(names(MALL_SE_In) %in% (NumericYears))]
  
  W4_SE_In<-W4_SE_In[which(names(W4_SE_In) %in% (NumericYears))]
  for(i in 2:length(W4_SE_In)){
    W4_SE_In[i]<-if_else(W4_SE_In[i]==0,W4_SE_In[i-1],W4_SE_In[i])
  }
  W3_SE_In<-W3_SE_In[which(names(W3_SE_In) %in% (NumericYears))]
  for(i in 2:length(W3_SE_In)){
    W3_SE_In[i]<-if_else(W3_SE_In[i]==0,W3_SE_In[i-1],W3_SE_In[i])
  }
  W2_SE_In<-W2_SE_In[which(names(W2_SE_In) %in% (NumericYears))]
  for(i in 2:length(W2_SE_In)){
    W2_SE_In[i]<-if_else(W2_SE_In[i]==0,W2_SE_In[i-1],W2_SE_In[i])
  }
  
  Y4_SE_In[Y4_SE_In==0|is.na(Y4_SE_In)] <- 0.0001
  Y4_lcl_In[Y4_lcl_In==0] <- 0.0001
  Y4_ucl_In[Y4_ucl_In==0] <- 0.0001
  
  Y3_SE_In[Y3_SE_In==0|is.na(Y3_SE_In)] <- 0.0001
  Y3_lcl_In[Y3_lcl_In==0] <- 0.0001
  Y3_ucl_In[Y3_ucl_In==0] <- 0.0001
  
  Y2_SE_In[Y2_SE_In==0|is.na(Y2_SE_In)]<- 0.0001
  Y2_lcl_In[Y2_lcl_In==0] <- 0.0001
  Y2_ucl_In[Y2_ucl_In==0] <- 0.0001
  
  YALL_SE_In[YALL_SE_In==0|is.na(YALL_SE_In)] <- 0.0001
  YALL_lcl_In[YALL_lcl_In==0] <- 0.0001
  YALL_ucl_In[YALL_ucl_In==0] <- 0.0001
  
  ###
  M4_SE_In[M4_SE_In==0] <- 0.0001
  M4_SE_In[is.na(M4_SE_In==0)] <- 0.0001
  
  M3_SE_In[M3_SE_In==0] <- 0.0001
  M3_SE_In[is.na(M3_SE_In==0)] <- 0.0001
  
  M2_SE_In[M2_SE_In==0] <- 0.0001
  M2_SE_In[is.na(M2_SE_In==0)] <- 0.0001
  
  MALL_SE_In[MALL_SE_In==0] <- 0.0001
  MALL_SE_In[is.na(MALL_SE_In==0)] <- 0.0001
  
  W4_SE_In[W4_SE_In==0|is.na(W4_SE_In)] <- 0.0001
  W3_SE_In[W3_SE_In==0|is.na(W3_SE_In)] <- 0.0001
  W2_SE_In[W2_SE_In==0|is.na(W2_SE_In)] <- 0.0001

  maxB<-0.0001
  maxRIN<-max(YALL_ucl_In,na.rm=T)
  

  Badults<-M3_In+M4_In
  
  RecVarLag<-c(1,2,3)
  RecWinterLag<-c(0,1,2)
  RecLagCombos<-expand.grid(RecVarLag,RecWinterLag)
 
  jags.params <- c("BR","B3","B4","B2","dum1","dum2","dum3","dum4","prDum",'lagDum','lagDum2','lagDum1',
                   'lagQuad1','lagQuad2','lagDumIN','lagQuad3','lagQuad4',
                   "N4","N3","N2","NALL",'NAdults','Surv','b1r','b12r',"b5r",'b52r',
                   'b1S','b12S','lSurv0','lSurv','b5S','b52S',"b6S",
                   'la0','la','b','R','sigmaR','sigmaSurv','Q','logLik')

  
  #adding years for predicition###
  if(Season=="Fall"|SS=="Ruby.Vigilante"){
    addYears<-abs(max(NumericYears)-2024)
  }else{
    addYears<-abs(max(NumericYears)-2025)
  }
  if(SS=="RockCreek.HogBack"|SS=="BigHole.JerryCreek"){
    addYears<-0
  }
  
  mcmc<-list()
  mcall<-list()
  mc_ll<-list()
  loos<-list()
  modDUMMY<-list()
  modsOUT<-list(list())
  mcallOUT<-list(list())
  mc_llOUT<-list(list())
  loosOUT<-list(list())
  selTabs<-list()
  loos<-list()
  effectSizeL<-list()
  effectSizeQ<-list()
  sims<-list()
  effectSize<-list()
  dummies<-list()
  time1<-Sys.time()
  simsOUT<-list()
  modName<-list()
  modelsLags<-list()
  if(cors[1,2]<0.5){
    global<-"yes"
  }else{
    global<-"no"
  }
  
  if(global=='yes'){
    lagsIN<-RecLagCombos
    modelsIN<-c('PVA_ProductionModel_AGE_NULL.txt',
                rep('PVA_ProductionModel_Global2d.txt',9))
    }else{
      lagsIN<-rbind(RecLagCombos[1:3,],RecLagCombos[c(1,4,7),])
      modelsIN<-c('PVA_ProductionModel_AGE_NULL.txt',
                  rep('PVA_ProductionModel_SUMMERQ2d.txt',3),
                  rep('PVA_ProductionModel_WINTQ2d.txt',3))
    }


  ####RUN MODELS####
  
  
  jags.data <- list(Y4 = c(Y4_In,rep(NA,addYears)),Y3 = c(Y3_In,rep(NA,addYears)),W_4 = c(W4_In,rep(median(W4_In),addYears)), W_3 = c(W3_In,rep(median(W3_In),addYears)),
                    Y2 = c(Y2_In,rep(NA,addYears)),W_2 = c(W2_In,rep(median(W2_In),addYears)),
                    YALL = c(YALL_In,rep(NA,addYears)),
                    Y4_SE = c(Y4_SE_In,rep(0.0001,addYears)),  W_4_SE = c(W4_SE_In,rep(sd(W4_SE_In),addYears)),
                    Y3_SE = c(Y3_SE_In,rep(0.0001,addYears)),  W_3_SE = c(W3_SE_In,rep(sd(W3_SE_In),addYears)),
                    Y2_SE = c(Y2_SE_In,rep(0.0001,addYears)),  W_2_SE = c(W2_SE_In,rep(sd(W2_SE_In),addYears)),
                    YALL_SE = c(YALL_SE_In,rep(0.0001,addYears)),
                    covarlag1 = c(covarLagIn1),covarlag5 = c(covarLagIn5),
                    nYears=length(YALL_In)+addYears,
                    sigmaRmax=5,
                    maxN4 = mean(Y4_In,na.rm=T),
                    maxN3 = mean(Y3_In,na.rm=T),
                    maxNstwo1 = mean(Y2_In,na.rm=T),
                    SDmaxN4 = sd(Y4_In,na.rm=T),
                    SDmaxN3 = sd(Y3_In,na.rm=T),
                    SDmaxNstwo1 = sd(Y2_In,na.rm=T),
                    maxB4 = mean(M4_In,na.rm=T),
                    maxB3 = mean(M3_In,na.rm=T),
                    maxBstwo1 = mean(M2_In,na.rm=T),
                    SDmaxB4 = sd(M4_In,na.rm=T),
                    SDmaxB3 = sd(M3_In,na.rm=T),
                    SDmaxBstwo1 = sd(M2_In,na.rm=T),
                    RecLag = 3,
                    maxB=maxB,
                    lagProb=c(0.5,0.5,0.5),pi=rep(1/3,3),  SummerLags=c(1,2,3),
                    WinterLags=c(0,1,2),
                    winterlag=winterlag,
                    SurvLag=SurvLag,powerSeason=1,
                    maxR=log(maxRIN))
  if(global=="yes"){
    jags.data <- list(Y4 = c(Y4_In,rep(NA,addYears)),Y3 = c(Y3_In,rep(NA,addYears)),W_4 = c(W4_In,rep(median(W4_In),addYears)), W_3 = c(W3_In,rep(median(W3_In),addYears)),
                      Y2 = c(Y2_In,rep(NA,addYears)),W_2 = c(W2_In,rep(median(W2_In),addYears)),
                      YALL = c(YALL_In,rep(NA,addYears)),
                      Y4_SE = c(Y4_SE_In,rep(0.0001,addYears)),  W_4_SE = c(W4_SE_In,rep(sd(W4_SE_In),addYears)),
                      Y3_SE = c(Y3_SE_In,rep(0.0001,addYears)),  W_3_SE = c(W3_SE_In,rep(sd(W3_SE_In),addYears)),
                      Y2_SE = c(Y2_SE_In,rep(0.0001,addYears)),  W_2_SE = c(W2_SE_In,rep(sd(W2_SE_In),addYears)),
                      YALL_SE = c(YALL_SE_In,rep(0.0001,addYears)),
                      covarlag1 = c(covarLagIn1),covarlag5 = c(covarLagIn5),
                      nYears=length(YALL_In)+addYears,
                      sigmaRmax=5,
                      maxN4 = mean(Y4_In,na.rm=T),
                      maxN3 = mean(Y3_In,na.rm=T),
                      maxNstwo1 = mean(Y2_In,na.rm=T),
                      SDmaxN4 = sd(Y4_In,na.rm=T),
                      SDmaxN3 = sd(Y3_In,na.rm=T),
                      SDmaxNstwo1 = sd(Y2_In,na.rm=T),
                      maxB4 = mean(M4_In,na.rm=T),
                      maxB3 = mean(M3_In,na.rm=T),
                      maxBstwo1 = mean(M2_In,na.rm=T),
                      SDmaxB4 = sd(M4_In,na.rm=T),
                      SDmaxB3 = sd(M3_In,na.rm=T),
                      SDmaxBstwo1 = sd(M2_In,na.rm=T),
                      RecLag = 3,
                      maxB=maxB,
                      lagProb=c(0.5,0.5,0.5),pi1=rep(1/3,3),pi2=rep(1/3,3),  SummerLags=c(1,2,3),
                      WinterLags=c(0,1,2),
                      winterlag=winterlag,
                      SurvLag=SurvLag,powerSeason=1,
                      maxR=log(maxRIN))
    
    mod_lm_global<- R2jags::jags.parallel(jags.data,
                                          parameters.to.save = jags.params,
                                          model.file = here::here('JAGS_PVA','JAGSmodels',"PVA_ProductionModel_Global.txt"),
                                          n.chains = 4,
                                          n.burnin = 750000/4,
                                          n.thin = 10,
                                          n.iter = 750000)
    
    lags<-list()
    for( i in 1:nrow(lagsIN)){
      lags[[i]]<-lagsIN[i,]
    }
    bool_df <- data.frame(
      var1 = sample(c(1, 0), 4, replace = TRUE),
      var2 = sample(c(1, 0), 4, replace = TRUE),
      var3 = sample(c(1, 0), 4, replace = TRUE),
      var4 = sample(c(1, 0), 4, replace = TRUE)
    )
    
    # Display the data frame
    print(bool_df)
    
    # Using expand.grid to get all possible combinations
    all_combinations <- expand.grid(
      var1 = c(1, 0),
      var2 = c(1, 0),
      var3 = c(1, 0),
      var4 = c(1, 0)
    )
    
    # Print the all_combinations data frame
    print(all_combinations)    
    
    bool_df <- data.frame(
      SummerRecQuad = sample(c('yes', 'no'), 4, replace = TRUE),
      WinterRecQuad = sample(c('yes', 'no'), 4, replace = TRUE),
      SummerSurvQuad = sample(c('yes', 'no'), 4, replace = TRUE),
      WinterSurvQuad = sample(c('yes', 'no'), 4, replace = TRUE)
    )
    
    # Display the data frame
    print(bool_df)
    
    # Using expand.grid to get all possible combinations
    all_combinations_names <- expand.grid(
      SummerRecQuad = c('yes', 'no'),
      WinterRecQuad = c('yes', 'no'),
      SummerSurvQuad = c('yes', 'no'),
      WinterSurvQuad = c('yes', 'no')
    )
    
    # Print the all_combinations data frame
    print(all_combinations_names)   
    Inclusion<-lapply(lags,function(ll){
      
      #inlcuision probability of the summer lag to reruitment
      lagInclusionSUM0<-length(mod_lm_global$BUGSoutput$sims.list$lagDum1[which(mod_lm_global$BUGSoutput$sims.list$lagDum1==which(SummerLags==ll$Var1) & mod_lm_global$BUGSoutput$sims.list$lagDum2==which(WinterLags==ll$Var2))])/length(mod_lm_global$BUGSoutput$sims.list$lagDum2==which(WinterLags==ll$Var2))

      #inlcuision probability of the winter lag to reruitment
      lagInclusionWINT0<-length(mod_lm_global$BUGSoutput$sims.list$lagDum2[which(mod_lm_global$BUGSoutput$sims.list$lagDum1==which(SummerLags==ll$Var1) & mod_lm_global$BUGSoutput$sims.list$lagDum2==which(WinterLags==ll$Var2))])/length(mod_lm_global$BUGSoutput$sims.list$lagDum1==which(SummerLags==ll$Var1))
      
      #inlcuision probability of the quadratic component of summer flow recruitment relationship
      quadInclusionRecSUM<-mod_lm_global$BUGSoutput$sims.list$lagQuad1[which(mod_lm_global$BUGSoutput$sims.list$lagDum1==which(SummerLags==ll$Var1) & mod_lm_global$BUGSoutput$sims.list$lagDum2==which(WinterLags==ll$Var2))]%>%mean()
      
      #inlcuision probability of the quadratic component of winter flow recruitment relationship
      quadInclusionRecWINT<-mod_lm_global$BUGSoutput$sims.list$lagQuad2[which(mod_lm_global$BUGSoutput$sims.list$lagDum1==which(SummerLags==ll$Var1) & mod_lm_global$BUGSoutput$sims.list$lagDum2==which(WinterLags==ll$Var2))]%>%mean()
      
      #inlcuision probability of the quadratic component of summer flow survival relationship
      quadInclusionSurvSUM<-mod_lm_global$BUGSoutput$sims.list$lagQuad3[which(mod_lm_global$BUGSoutput$sims.list$lagDum1==which(SummerLags==ll$Var1) & mod_lm_global$BUGSoutput$sims.list$lagDum2==which(WinterLags==ll$Var2))]%>%mean()
      
      #inlcuision probability of the quadratic component of winter flow survival relationship
      quadInclusionSurvWINT<-mod_lm_global$BUGSoutput$sims.list$lagQuad4[which(mod_lm_global$BUGSoutput$sims.list$lagDum1==which(SummerLags==ll$Var1) & mod_lm_global$BUGSoutput$sims.list$lagDum2==which(WinterLags==ll$Var2))]%>%mean()
      
      #determine full inclusion for lag comnination and full quadratic effect
      for(i in 1:nrow(all_combinations_names)){
        all_combinations_names$lagInclusion[i]<-length(mod_lm_global$BUGSoutput$sims.list$lagDum1[which(mod_lm_global$BUGSoutput$sims.list$lagDum1==which(SummerLags==ll$Var1) & mod_lm_global$BUGSoutput$sims.list$lagDum2==which(WinterLags==ll$Var2) & mod_lm_global$BUGSoutput$sims.list$lagQuad1==all_combinations[i,1] &  mod_lm_global$BUGSoutput$sims.list$lagQuad2==all_combinations[i,2] &  mod_lm_global$BUGSoutput$sims.list$lagQuad3==all_combinations[i,3] & mod_lm_global$BUGSoutput$sims.list$lagQuad4==all_combinations[i,4]) ])/length(mod_lm_global$BUGSoutput$sims.list$lagDum1[which( mod_lm_global$BUGSoutput$sims.list$lagDum2==which(WinterLags==ll$Var2)) ])
      }
      
      out<-all_combinations_names%>%mutate(SummerLag=ll$Var1,WinterLag=ll$Var2,
                                           AllInclusion=lagInclusionSUM0,quadInclusionRecSUM=quadInclusionRecSUM,
                                           quadInclusionRecWINT=quadInclusionRecWINT,quadInclusionSurvSUM=quadInclusionSurvSUM,
                                           quadInclusionSurvWINT=quadInclusionSurvWINT)
      return(out)
    })%>%bind_rows()

  }else{
    #For non global models
    
    #Summer Model
    mod_lm_summer<- R2jags::jags.parallel(jags.data,
                                          parameters.to.save = jags.params,
                                          model.file = here::here('JAGS_PVA','JAGSmodels',"PVA_ProductionModel_Summer.txt"),
                                          n.chains = 4,
                                          n.burnin = 750000/4,
                                          n.thin = 10,
                                          n.iter = 750000)
    
    lags<-list()
    for( i in 1:3){
      lags[[i]]<-lagsIN[i,]
    }
    bool_df <- data.frame(
      var1 = sample(c(1, 0), 2, replace = TRUE),
      var3 = sample(c(1, 0), 2, replace = TRUE)
    )
    
    # Display the data frame
    print(bool_df)
    
    # Using expand.grid to get all possible combinations
    all_combinations <- expand.grid(
      var1 = c(1, 0),
      var3 = c(1, 0)
    )
    
    # Print the all_combinations data frame
    print(all_combinations)    
    
    bool_df <- data.frame(
      SummerRecQuad = sample(c('yes', 'no'), 2, replace = TRUE),
      SummerSurvQuad = sample(c('yes', 'no'), 2, replace = TRUE)
    )
    
    # Display the data frame
    print(bool_df)
    
    # Using expand.grid to get all possible combinations
    all_combinations_names <- expand.grid(
      SummerRecQuad = c('yes', 'no'),
      SummerSurvQuad = c('yes', 'no')
    )
    
    # Print the all_combinations data frame
    print(all_combinations_names)   
    InclusionSumer<-lapply(lags,function(ll){
      # ll<-lags[[1]]
      
      #Summer recruitment lag inclusion
      #inlcuision probability of the summer lag to reruitment
      lagInclusionSUM0<-length(mod_lm_summer$BUGSoutput$sims.list$lagDum1[which(mod_lm_summer$BUGSoutput$sims.list$lagDum1==which(SummerLags==ll$Var1) & mod_lm_summer$BUGSoutput$sims.list$lagDum2==which(WinterLags==ll$Var2))])/length(mod_lm_summer$BUGSoutput$sims.list$lagDum2==which(WinterLags==ll$Var2))
      
      #inlcuision probability of the quadratic component of summer flow recruitment relationship
      quadInclusionRecSUM<-mod_lm_summer$BUGSoutput$sims.list$lagQuad1[which(mod_lm_summer$BUGSoutput$sims.list$lagDum1==which(SummerLags==ll$Var1) & mod_lm_summer$BUGSoutput$sims.list$lagDum2==which(WinterLags==ll$Var2))]%>%mean()
 
      #inlcuision probability of the quadratic component of summer flow survival relationship
      quadInclusionSurvSUM<-mod_lm_summer$BUGSoutput$sims.list$lagQuad3[which(mod_lm_summer$BUGSoutput$sims.list$lagDum1==which(SummerLags==ll$Var1) & mod_lm_summer$BUGSoutput$sims.list$lagDum2==which(WinterLags==ll$Var2))]%>%mean()

      for(i in 1:nrow(all_combinations_names)){
        #Recruitment lag inclusion for linear flow effects with quadratic effects excluded
        # i<-1
        all_combinations_names$lagInclusion[i]<-length(mod_lm_summer$BUGSoutput$sims.list$lagDum[which(mod_lm_summer$BUGSoutput$sims.list$lagDum==which(SummerLags==ll$Var1)  & mod_lm_summer$BUGSoutput$sims.list$lagQuad1==all_combinations[i,1] &  mod_lm_summer$BUGSoutput$sims.list$lagQuad2==all_combinations[i,2] ) ])/length(mod_lm_summer$BUGSoutput$sims.list$lagDum)
        
      }
    
      out<-all_combinations_names%>%mutate(SummerLag=ll$Var1,WinterLag=NA,AllInclusion=lagInclusionSUM0,quadInclusionRecSUM=quadInclusionRecSUM, quadInclusionRecWINT=NA,quadInclusionSurvSUM=quadInclusionSurvSUM,quadInclusionSurvWINT=NA)
      return(out)
    })%>%bind_rows()

  #Winter model
    
    mod_lm_winter<- R2jags::jags.parallel(jags.data,
                                          parameters.to.save = jags.params,
                                          model.file = here::here('JAGS_PVA','JAGSmodels',"PVA_ProductionModel_Winter.txt"),
                                          n.chains = 4,
                                          n.burnin = 750000/4,
                                          n.thin = 10,
                                          n.iter = 750000)
    
    lags<-list()
    for( i in c(4,5,6)){
      lags[[i]]<-lagsIN[i,]
    }
    bool_df <- data.frame(
      var2 = sample(c(1, 0), 2, replace = TRUE),
      var4 = sample(c(1, 0), 2, replace = TRUE)
    )
    
    # Display the data frame
    print(bool_df)
    
    # Using expand.grid to get all possible combinations
    all_combinations <- expand.grid(
      var2 = c(1, 0),
      var4 = c(1, 0)
    )
    
    # Print the all_combinations data frame
    print(all_combinations)    
    
    bool_df <- data.frame(
      WinterRecQuad = sample(c('yes', 'no'), 2, replace = TRUE),
      WinterSurvQuad = sample(c('yes', 'no'), 2, replace = TRUE)
    )
    
    # Display the data frame
    print(bool_df)
    
    # Using expand.grid to get all possible combinations
    all_combinations_names <- expand.grid(
      WinterRecQuad = c('yes', 'no'),
      WinterSurvQuad = c('yes', 'no')
    )
    
    # Print the all_combinations data frame
    print(all_combinations_names)   
    InclusionWinter<-lapply(lags[4:6],function(ll){
      #inlcuision probability of the winter lag to reruitment
      lagInclusionWINT0<-length(mod_lm_winter$BUGSoutput$sims.list$lagDum2[which(mod_lm_winter$BUGSoutput$sims.list$lagDum1==which(SummerLags==ll$Var1) & mod_lm_winter$BUGSoutput$sims.list$lagDum2==which(WinterLags==ll$Var2))])/length(mod_lm_winter$BUGSoutput$sims.list$lagDum1==which(SummerLags==ll$Var1))

      #inlcuision probability of the quadratic component of winter flow recruitment relationship
      quadInclusionRecWINT<-mod_lm_winter$BUGSoutput$sims.list$lagQuad2[which(mod_lm_winter$BUGSoutput$sims.list$lagDum1==which(SummerLags==ll$Var1) & mod_lm_winter$BUGSoutput$sims.list$lagDum2==which(WinterLags==ll$Var2))]%>%mean()
      
      #inlcuision probability of the quadratic component of winter flow survival relationship
      quadInclusionSurvWINT<-mod_lm_winter$BUGSoutput$sims.list$lagQuad4[which(mod_lm_winter$BUGSoutput$sims.list$lagDum1==which(SummerLags==ll$Var1) & mod_lm_winter$BUGSoutput$sims.list$lagDum2==which(WinterLags==ll$Var2))]%>%mean()
      
      for(i in 1:nrow(all_combinations_names)){
        #Recruitment lag inclusion for linear flow effects with quadratic effects excluded
        all_combinations_names$lagInclusion[i]<-length(mod_lm_winter$BUGSoutput$sims.list$lagDum[which(mod_lm_winter$BUGSoutput$sims.list$lagDum==which(WinterLags==ll$Var2)  & mod_lm_winter$BUGSoutput$sims.list$lagQuad1==all_combinations[i,1] &  mod_lm_winter$BUGSoutput$sims.list$lagQuad2==all_combinations[i,2] ) ])/length(mod_lm_winter$BUGSoutput$sims.list$lagDum)
        
      }
      out<-all_combinations_names%>%mutate(SummerLag=NA,WinterLag=ll$Var2,AllInclusion=lagInclusionWINT0,quadInclusionRecSUM=NA, quadInclusionRecWINT=quadInclusionRecWINT,quadInclusionSurvSUM=NA,quadInclusionSurvWINT=quadInclusionSurvWINT)
      return(out)
    })%>%bind_rows()
    Inclusion<-rbind(InclusionSumer%>%mutate(WinterRecQuad=NA,WinterSurvQuad=NA),InclusionWinter%>%mutate(SummerRecQuad=NA,SummerSurvQuad=NA))%>%
      select(SummerLag,WinterLag,SummerRecQuad,SummerSurvQuad,WinterRecQuad,WinterSurvQuad,lagInclusion,AllInclusion,quadInclusionRecSUM,quadInclusionRecWINT,quadInclusionSurvSUM,quadInclusionSurvWINT)

  
  write.csv(Inclusion,here::here(dir_name2,paste0(SS,"_RecLagInclusionProbQuad.csv")))
  modelsAndLags<-data.frame(model=modelsIN,SummerLag=lagsIN[,1],WinterLag=lagsIN[,2])
  
  
  ####RUN SELECTED MODELS####
  if(cors[1,2]<0.5 ){
    global<-"yes"
  }else{
    global<-"no"
  }
  
  topLags<-read.csv(here::here(dir_name2,paste0(SS,"_RecLagInclusionProbQuad.csv")))%>%
    mutate(Season=if_else(!is.na(WinterLag) ,"Winter","Summer"))
  topLagsOUT<-read.csv(here::here(dir_name2,paste0(SS,"_RecLagInclusionProbQuad.csv")))%>%
    mutate(Season=if_else(!is.na(WinterLag) ,"Winter","Summer"))%>%group_by(Season)%>%slice_max(lagInclusion)
  
  modelsAndLagsIN<-modelsAndLags
  
  if(global=='yes'){
    
    
    for(i in 2:nrow(modelsAndLagsIN)){  
      # i<-3
      print(i)
      # topLags%>%slice_max(quadInclusionSurvSUM)
      seasonInclusion<-(topLags%>%filter(SummerLag==modelsAndLagsIN$SummerLag[i] & WinterLag==modelsAndLagsIN$WinterLag[i])%>%slice_max(lagInclusion)%>%dplyr::select(SummerRecQuad,WinterRecQuad,SummerSurvQuad,WinterSurvQuad))
      check<-c(seasonInclusion$SummerRecQuad,seasonInclusion$WinterRecQuad,seasonInclusion$SummerSurvQuad,seasonInclusion$WinterSurvQuad)
      modelsAndLagsIN$model[i]<-if_else(identical(check,c('yes','no','no','no')),"PVA_ProductionModel_Global2rs.txt",modelsAndLagsIN$model[i])
      modelsAndLagsIN$model[i]<-if_else(identical(check,c('no','yes','no','no')),"PVA_ProductionModel_Global2rw.txt",modelsAndLagsIN$model[i])
      modelsAndLagsIN$model[i]<-if_else(identical(check,c('yes','yes','no','no')),"PVA_ProductionModel_Global2r.txt",modelsAndLagsIN$model[i])
      modelsAndLagsIN$model[i]<-if_else(identical(check,c('no','yes','no','yes')),"PVA_ProductionModel_Global2w.txt",modelsAndLagsIN$model[i])
      modelsAndLagsIN$model[i]<-if_else(identical(check,c('yes','no','yes','no')),"PVA_ProductionModel_Global2s.txt",modelsAndLagsIN$model[i])
      modelsAndLagsIN$model[i]<-if_else(identical(check,c('yes','no','yes','yes')),"PVA_ProductionModel_Global2srSs.txt",modelsAndLagsIN$model[i])
      modelsAndLagsIN$model[i]<-if_else(identical(check,c('no','yes','yes','yes')),"PVA_ProductionModel_Global2wrSs.txt",modelsAndLagsIN$model[i])
      modelsAndLagsIN$model[i]<-if_else(identical(check,c('yes','no','no','yes')),"PVA_ProductionModel_Global2wrSw.txt",modelsAndLagsIN$model[i])
      modelsAndLagsIN$model[i]<-if_else(identical(check,c('no','no','yes','no')),"PVA_ProductionModel_Global2Ss.txt",modelsAndLagsIN$model[i])
      modelsAndLagsIN$model[i]<-if_else(identical(check,c('no','no','no','yes')),"PVA_ProductionModel_Global2Sw.txt",modelsAndLagsIN$model[i])
      modelsAndLagsIN$model[i]<-if_else(identical(check,c('no','no','yes','yes')),"PVA_ProductionModel_Global2Surv.txt",modelsAndLagsIN$model[i])
      modelsAndLagsIN$model[i]<-if_else(identical(check,c('no','no','no','no')),"PVA_ProductionModel_Global.txt",modelsAndLagsIN$model[i])
      modelsAndLagsIN$model[i]<-if_else(identical(check,c('yes','yes','yes','yes')),"PVA_ProductionModel_Global2.txt",modelsAndLagsIN$model[i])
      
    }
    
    
    for(m in 1:nrow(modelsAndLagsIN)){
      print(paste0(modelsAndLagsIN$model[m],"_",modelsAndLagsIN$SummerLag[m],modelsAndLagsIN$WinterLag[m]))
    
      print(modelsAndLagsIN$model[m])
      mod_lm<- R2jags::jags.parallel(jags.data, 
                                     parameters.to.save = jags.params,
                                     model.file = here::here('JAGS_PVA','JAGSmodels',modelsAndLagsIN$model[m]), 
                                     n.chains = 4, 
                                     n.burnin = 80000/4, 
                                     n.thin = 10,
                                     n.iter = 80000)
      modsOUT[[m]]<-mod_lm
      
    }
    
    yearsIN<-as.numeric(names(Y2_In))
    yearsOUT<-c(min(as.numeric(names(Y2_In))):(max(as.numeric(names(Y2_In)))+addYears))
    
    paramsALLmods<-  lapply(1:nrow(modelsAndLagsIN),FUN=function(mm){
      modIN<-modsOUT[[mm]]
      simsIN<-modIN$BUGSoutput$sims.list
      params<-names(simsIN)
      loadParams<-lapply(params,FUN=function(pp){
        print(pp)
        param<-   data.frame(simsIN[which(names(simsIN)==pp)])%>% apply(.,2,FUN=quantile,c(0.1,0.5,0.9))
        paramOUT<-data.frame(name=paste0("Estimated ",pp),N10=param[1,],N50=param[2,],N90=param[3,])
        if(length(param[1,])==length(yearsOUT)){
          paramOUT$Year<-yearsOUT}else{
            if(length(param[1,])==(length(yearsOUT)-RecLag)){
              paramOUT$Year<-c((min(yearsOUT)+3):max(yearsOUT))
            }else{
              if(length(param[1,])==1){
                paramOUT$Year<-NA}
            }
          }
        return(paramOUT)
      })%>%bind_rows()%>%mutate(Stream.Section=SS,model=modelsAndLagsIN$model[mm],SummerLag=modelsAndLagsIN$SummerLag[mm],WinterLag=modelsAndLagsIN$WinterLag[mm])
    })%>%bind_rows()
    
    
    paramsOUT<-paramsALLmods%>% bind_rows()%>%rbind( data.frame(name=c(rep("Observed NALL",length(YALL_In)),
                                                                       rep("Observed N4",length(Y4_In)),
                                                                       rep("Observed N3",length(Y3_In)),
                                                                       rep("Observed N2",length(Y2_In))),
                                                                N50=c(YALL_In,Y4_In,Y3_In,Y2_In),
                                                                N10=c(YALL_lcl_In,Y4_lcl_In,Y3_lcl_In,Y2_lcl_In),
                                                                N90=c(YALL_ucl_In,Y4_ucl_In,Y3_ucl_In,Y2_ucl_In),
                                                                Year=c(names(YALL_In),names(Y4_In),names(Y3_In),names(Y2_In)),
                                                                Stream.Section=SS,model=NA,
                                                                SummerLag=NA, WinterLag=NA
    )%>%
      dplyr::mutate(SummerLag=NA,WinterLag=NA,model=NA))%>%
      dplyr::mutate(Year=if_else(is.na(Year) & grepl("Estimated N",name),"2024",Year))%>%
      mutate(sig=if_else(N10/N90>0,"1","0"))

    write.csv(paramsOUT,here::here(dir_name2,"csvs_quadratic",paste0(Stream,Section,"allParams.csv")))
    topLagsAndEffects<-read.csv(here::here(dir_name2,paste0(SS,"_RecLagInclusionProbQuad.csv")))%>%
      mutate(Season=if_else(!is.na(WinterLag) ,"Winter","Summer"))%>%left_join(paramsOUT%>%filter(grepl('b1',name)|grepl('b5',name)))
    
    topLagsAndEffects%>%group_by(Season)%>%summarize(sig=max(sig,na.rm=T))
    if("1"%in%c(topLagsAndEffects%>%group_by(Season)%>%summarize(sig=max(sig,na.rm=T))%>%pull(sig))){
      topLagsOUT<- read.csv(here::here(dir_name2,paste0(SS,"_RecLagInclusionProbQuad.csv")))%>%
        mutate(Season=if_else(!is.na(WinterLag) ,"Winter","Summer"))%>%left_join(paramsOUT%>%filter(grepl('b1',name)|grepl('b5',name)))%>%
        filter(sig=='1') %>%group_by(Season)%>%slice_max(AllInclusion)%>%slice_max(lagInclusion)
      
    }else{
      topLagsOUT<-read.csv(here::here(dir_name2,paste0(SS,"_RecLagInclusionProbQuad.csv")))%>%
        mutate(Season=if_else(!is.na(WinterLag) ,"Winter","Summer"))%>%group_by(Season)%>%slice_max(AllInclusion)%>%slice_max(lagInclusion)
      
    }
    paramsOUT<-paramsOUT%>%
      mutate(topMod=if_else(grepl('Global',model) & SummerLag==topLagsOUT$SummerLag[which(topLagsOUT$Season=='Winter')][1] & WinterLag==topLagsOUT$WinterLag[which(topLagsOUT$Season=='Winter')][1],'yes','no'))
    
    topMod<-modsOUT[[which(modelsAndLagsIN$model!="PVA_ProductionModel_AGE_NULL.txt" & modelsAndLagsIN$SummerLag==topLagsOUT$SummerLag[1] & modelsAndLagsIN$WinterLag==topLagsOUT$WinterLag[1])]]
    saveRDS(topMod, here::here("JAGS_PVA","ModelFits",paste0(SS,"_IndicatorVarSel_TopMod.rds")))
    write.csv(paramsOUT,here::here(dir_name2,"csvs_quadratic",paste0(Stream,Section,"allParams.csv")))
    
    
  }else{
    modelsAndLagsINSUMMER<-modelsAndLagsIN%>%filter(grepl('SUMMER',model))
    
    for(i in 1:nrow(modelsAndLagsIN%>%filter(grepl('SUMMER',model)))){
      seasonInclusion<-(topLags%>%mutate(SummerLag=if_else(is.na(SummerLag),1,SummerLag),
                                         WinterLag=if_else(is.na(WinterLag),0,WinterLag))%>%
                          filter(Season=='Summer',SummerLag==modelsAndLagsINSUMMER$SummerLag[i] & WinterLag==modelsAndLagsINSUMMER$WinterLag[i])%>%
                          slice_max(lagInclusion)%>%dplyr::select(SummerRecQuad,WinterRecQuad,SummerSurvQuad,WinterSurvQuad,Season))
      check<-c(seasonInclusion$SummerRecQuad,seasonInclusion$WinterRecQuad,seasonInclusion$SummerSurvQuad,seasonInclusion$WinterSurvQuad)
      
      modelsAndLagsINSUMMER$model[i]<-if_else( identical(check,c('yes',NA,'no',NA)),"PVA_ProductionModel_SUMMERQ2r.txt",modelsAndLagsINSUMMER$model[i])
      modelsAndLagsINSUMMER$model[i]<-if_else(identical(check,c('no',NA,'yes',NA)) ,"PVA_ProductionModel_SUMMERQ2S.txt",modelsAndLagsINSUMMER$model[i])
      modelsAndLagsINSUMMER$model[i]<-if_else(identical(check,c('no',NA,'no',NA)),"PVA_ProductionModel_SUMMERQ.txt",modelsAndLagsINSUMMER$model[i])
      modelsAndLagsINSUMMER$model[i]<-if_else(identical(check,c('yes',NA,'yes',NA)),"PVA_ProductionModel_SUMMERQ2.txt",modelsAndLagsINSUMMER$model[i])
      
    }
    
    modelsAndLagsINWINTER<-modelsAndLagsIN%>%filter(grepl('WINT',model))
    
    for(i in 1:nrow(modelsAndLagsIN%>%filter(grepl('WINT',model)))){
      # i<-2
      seasonInclusion<-(topLags%>%mutate(SummerLag=if_else(is.na(SummerLag),1,SummerLag),
                                         WinterLag=if_else(is.na(WinterLag),0,WinterLag))%>%
                          filter(Season=='Winter',SummerLag==modelsAndLagsINWINTER$SummerLag[i] & WinterLag==modelsAndLagsINWINTER$WinterLag[i])%>%
                          slice_max(lagInclusion)%>%dplyr::select(SummerRecQuad,WinterRecQuad,SummerSurvQuad,WinterSurvQuad,Season))
      check<-c(seasonInclusion$SummerRecQuad,seasonInclusion$WinterRecQuad,seasonInclusion$SummerSurvQuad,seasonInclusion$WinterSurvQuad)
      
      modelsAndLagsINWINTER$model[i]<-if_else( identical(check,c(NA,'yes',NA,'no')),"PVA_ProductionModel_WINTQ2r.txt",modelsAndLagsINWINTER$model[i])
      modelsAndLagsINWINTER$model[i]<-if_else(identical(check,c(NA,'no',NA,'yes')) ,"PVA_ProductionModel_WINTQ2S.txt",modelsAndLagsINWINTER$model[i])
      modelsAndLagsINWINTER$model[i]<-if_else(identical(check,c(NA,'no',NA,'no')),"PVA_ProductionModel_WINTQ.txt",modelsAndLagsINWINTER$model[i])
      modelsAndLagsINWINTER$model[i]<-if_else(identical(check,c(NA,'yes',NA,'yes')),"PVA_ProductionModel_WINTQ2.txt",modelsAndLagsINWINTER$model[i])
      
    }
    modelsAndLagsIN<-rbind(modelsAndLagsINWINTER,modelsAndLagsINSUMMER)
    for(m in 1:nrow(modelsAndLagsIN)){
      print(paste0(modelsAndLagsIN$model[m],"_",modelsAndLagsIN$SummerLag[m],modelsAndLagsIN$WinterLag[m]))
     
      mod_lm<- R2jags::jags.parallel(jags.data, 
                                     parameters.to.save = jags.params,
                                     model.file = here::here('JAGS_PVA','JAGSmodels',modelsAndLagsIN$model[m]), 
                                     n.chains = 4, 
                                     n.burnin = 80000/4, 
                                     n.thin = 10,
                                     n.iter = 80000)
      modsOUT[[m]]<-mod_lm
      
    }
    
    
    yearsIN<-as.numeric(names(Y2_In))
    yearsOUT<-c(min(as.numeric(names(Y2_In))):(max(as.numeric(names(Y2_In)))+addYears))
    
    paramsALLmods<-  lapply(1:nrow(modelsAndLagsIN),FUN=function(mm){
      # mm<-1
      modIN<-modsOUT[[mm]]
      simsIN<-modIN$BUGSoutput$sims.list
      params<-names(simsIN)
      loadParams<-lapply(params,FUN=function(pp){
        print(pp)
        param<-   data.frame(simsIN[which(names(simsIN)==pp)])%>% apply(.,2,FUN=quantile,c(0.1,0.5,0.9))
        paramOUT<-data.frame(name=paste0("Estimated ",pp),N10=param[1,],N50=param[2,],N90=param[3,])
        if(length(param[1,])==length(yearsOUT)){
          paramOUT$Year<-yearsOUT}else{
            if(length(param[1,])==(length(yearsOUT)-RecLag)){
              paramOUT$Year<-c((min(yearsOUT)+3):max(yearsOUT))
            }else{
              if(length(param[1,])==1){
                paramOUT$Year<-NA}
            }
          }
        return(paramOUT)
      })%>%bind_rows()%>%mutate(Stream.Section=SS,model=modelsAndLagsIN$model[mm],SummerLag=modelsAndLagsIN$SummerLag[mm],WinterLag=modelsAndLagsIN$WinterLag[mm])
    })%>%bind_rows()
    
    
    paramsOUT<-paramsALLmods%>% bind_rows()%>%rbind( data.frame(name=c(rep("Observed NALL",length(YALL_In)),
                                                                       rep("Observed N4",length(Y4_In)),
                                                                       rep("Observed N3",length(Y3_In)),
                                                                       rep("Observed N2",length(Y2_In))),
                                                                N50=c(YALL_In,Y4_In,Y3_In,Y2_In),
                                                                N10=c(YALL_lcl_In,Y4_lcl_In,Y3_lcl_In,Y2_lcl_In),
                                                                N90=c(YALL_ucl_In,Y4_ucl_In,Y3_ucl_In,Y2_ucl_In),
                                                                Year=c(names(YALL_In),names(Y4_In),names(Y3_In),names(Y2_In)),
                                                                Stream.Section=SS,model=NA,
                                                                SummerLag=NA, WinterLag=NA
    )%>%
      dplyr::mutate(SummerLag=NA,WinterLag=NA,model=NA))%>%
      dplyr::mutate(Year=if_else(is.na(Year) & grepl("Estimated N",name),"2024",Year))%>%
      mutate(sig=if_else(N10/N90>0,"1","0"))
    topLagsAndEffects<-read.csv(here::here(dir_name2,paste0(SS,"_RecLagInclusionProbQuad.csv")))%>%
      mutate(Season=if_else(!is.na(WinterLag) ,"Winter","Summer"))%>%left_join(paramsOUT%>%filter(grepl('b1',name)|grepl('b5',name))%>%
                                                                                 mutate(WinterLag=if_else(grepl("SUMMER",model),NA,WinterLag),
                                                                                        SummerLag=if_else(grepl("WINT",model),NA,SummerLag)))
    
    topLagsAndEffects%>%group_by(Season)%>%summarize(sig=max(sig))
    if("1"%in%c(topLagsAndEffects%>%filter(Season=='Summer')%>%summarize(sig=max(sig))%>%pull(sig))){
      topLagsOUTsummer<- read.csv(here::here(dir_name2,paste0(SS,"_RecLagInclusionProbQuad.csv")))%>%
        mutate(Season=if_else(!is.na(WinterLag) ,"Winter","Summer"))%>%left_join(paramsOUT%>%filter(grepl('b1',name)|grepl('b5',name))%>%
                                                                                   mutate(WinterLag=if_else(grepl("SUMMER",model),NA,WinterLag),
                                                                                          SummerLag=if_else(grepl("WINT",model),NA,SummerLag)))%>%filter(sig=='1' & Season=="Summer") %>%group_by(Season)%>%slice_max(AllInclusion)%>%slice_max(lagInclusion)
      
    }else{
      topLagsOUTsummer<-read.csv(here::here(dir_name2,paste0(SS,"_RecLagInclusionProbQuad.csv")))%>%
        mutate(Season=if_else(!is.na(WinterLag) ,"Winter","Summer"))%>%group_by(Season)%>%filter(Season=="Summer")%>%slice_max(AllInclusion)%>%slice_max(AllInclusion)%>%slice_max(lagInclusion)
      
    }
    
    if("1"%in%c(topLagsAndEffects%>%filter(Season=='Winter')%>%summarize(sig=max(sig))%>%pull(sig))){
      topLagsOUTwinter<- read.csv(here::here(dir_name2,paste0(SS,"_RecLagInclusionProbQuad.csv")))%>%
        mutate(Season=if_else(!is.na(WinterLag) ,"Winter","Summer"))%>%left_join(paramsOUT%>%filter(grepl('b1',name)|grepl('b5',name))%>%
                                                                                   mutate(WinterLag=if_else(grepl("SUMMER",model),NA,WinterLag),
                                                                                          SummerLag=if_else(grepl("WINT",model),NA,SummerLag)))%>%filter(sig=='1' & Season=="Winter") %>%group_by(Season)%>%slice_max(AllInclusion)%>%slice_max(lagInclusion)
      
    }else{
      topLagsOUTwinter<-read.csv(here::here(dir_name2,paste0(SS,"_RecLagInclusionProbQuad.csv")))%>%
        mutate(Season=if_else(!is.na(WinterLag) ,"Winter","Summer"))%>%group_by(Season)%>%filter(Season=="Winter")%>%slice_max(AllInclusion)%>%slice_max(lagInclusion)
      
    }
    
    topLagsOUT<-rbind(topLagsOUTsummer,topLagsOUTwinter)
    
    
    topModSummer<-modsOUT[[which(grepl('SUMMER',modelsAndLagsIN$model) & modelsAndLagsIN$SummerLag==topLagsOUT$SummerLag[which(topLagsOUT$Season=='Summer')][1] )]]
    saveRDS(topModSummer, here::here("JAGS_PVA","ModelFits",paste0(SS,"_IndicatorVarSel_topModSummer.rds")))
    topModSummer<-modsOUT[[which(grepl('WINT',modelsAndLagsIN$model) & modelsAndLagsIN$WinterLag==topLagsOUT$WinterLag[which(topLagsOUT$Season=='Winter')][1])]]
    saveRDS(topModSummer, here::here("JAGS_PVA","ModelFits",paste0(SS,"_IndicatorVarSel_topModWinter.rds")))
    
    paramsOUT<-paramsOUT%>%
      mutate(topMod=if_else(grepl('SUMMER',model) & SummerLag==topLagsOUT$SummerLag[which(topLagsOUT$Season=='Summer')],'yes','no'),
             topMod=if_else(grepl('WINT',model) & WinterLag==topLagsOUT$WinterLag[which(topLagsOUT$Season=='Winter')],'yes',topMod))
    
    write.csv(paramsOUT,here::here(dir_name2,"csvs_quadratic",paste0(Stream,Section,"allParams.csv")))
  }
  
})
