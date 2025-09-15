
rm(list=ls())

library(rjags)
library(R2jags)
library(here)
library(dplyr)
library(tidyverse)
modelScript.name <- "MIX_Hier.txt"
MIXjagsscript <- cat(
  "
  model {
        # Likelihood:
        for(i in 1:N) {
          y[i] ~ dnorm(mu[clust[i],Year[i]], pow(exp(lsigma[clust[i],Year[i]]),-2))
          clust[i] ~ dcat(pi[1:K,Year[i]]);
        }
     
        # Priors:
        for(A in 1:K){
            mu.hyper[A] ~ dnorm(mu.prior[A],pow(exp(lmu.A.sigma[A]),-2));
            lmu.A.sigma[A] ~ dunif(0.0001,5);
            lsigma.hyper[A] ~ dunif(0.0001,5);
            # lsigma.A.sigma[A] ~ dunif(0.0001,5);
          }
     
        for(Y in 1:nYear){
          for(A in 1:K){
            mu[A,Y] ~ dnorm(mu.hyper[A],pow(exp(lmu.A.sigma[A]),-2)) T(mu.prior.min[A],mu.prior.max[A]);
            # lsigma[A,Y] ~ dnorm(lsigma.hyper[A],pow(exp(lsigma.A.sigma[A]),-2));
            lsigma[A,Y] <- lsigma.hyper[A];
          }
          pi[1:K,Y] ~ ddirch(pi_prior);
        }
       
    }
    ",
  
  file = modelScript.name)

##Read in cyphers to make names in data match
SectionCypher <- read.csv(here::here('Data','SectionSynonyms.csv'),header=T,stringsAsFactors = F)
RiverCypher <- read.csv(here::here('Data','RiverSynonyms.csv'),header=T,stringsAsFactor=F)

#Read in data
RBT_extract <- read.csv(here::here('Data','RainbowExtract_07112025.csv'),header=T,stringsAsFactors = F)


###Change names
inWater<-RBT_extract%>%pull(Water)%>%unique()
inSection<-RBT_extract%>%pull(Section)%>%unique()

if(sum(!(RBT_extract$Section %in% SectionCypher$FWPsection))>0 | sum(!(RBT_extract$Water %in% RiverCypher$FWPriver))>0){
  print('Error: Sections/Rivers appear in data that are not in the cyphers')
}else{
  SectionCypher_in <- SectionCypher %>% na.omit()
  RBT_extract <- RBT_extract %>% filter(Water %in% RiverCypher$FWPriver & Section %in% SectionCypher_in$FWPsection & Species == 'Rainbow Trout')
  RBT_extract$Water <- RiverCypher$NewRiver[match(RBT_extract$Water,RiverCypher$FWPriver)]
  RBT_extract$Section <- SectionCypher_in$NewSection[match(RBT_extract$Section,SectionCypher_in$FWPsection)]
  if(T){
    RBT_extract <- RBT_extract %>% filter(paste(Water,Section,sep='.') !='Ruby.Maloney') #Maloney Replaced with Vigilante
  }}

#Determine the sampling season
SelectSeason<-RBT_extract%>%mutate(Section=if_else(Section=="Disease Sampling 2019 Hildreth",'Hildreth Section',Section))%>%
  mutate(Water.Section=paste(Water,Section,sep=","),
         Water.Year=paste(Water,Year,sep="."))%>%
  mutate(jday=yday(mdy(Date)),
         SampleSeason = ifelse(jday<182,'Spring','Fall'),SampleSeasonCode = ifelse(jday<182,1,2))%>%
  group_by(Water.Section,SampleSeason)%>%
  dplyr::summarize(count=n())%>%
  ungroup()%>%
  spread(SampleSeason,count)%>%
  mutate(Fall=if_else(is.na(Fall),0,Fall),
         Spring=if_else(is.na(Spring),0,Spring))%>%
  mutate(selectSeason=if_else( Spring>Fall ,"Spring","Fall"))%>%
  dplyr::select(Water.Section,selectSeason)


###Only selec the season with the most number of years in timeseries and filter out small fish.
RBT_extractIN<-RBT_extract%>%
  mutate(Section=if_else(Section=="Disease Sampling 2019 Hildreth",'Hildreth Section',Section))%>%
  mutate(Water.Section=paste(Water,Section,sep=","),
         Water.Year=paste(Water,Year,sep="."))%>%
  mutate(jday=yday(mdy(Date)),
         SampleSeason = ifelse(jday<182,'Spring','Fall'),SampleSeasonCode = ifelse(jday<182,1,2))%>%
  left_join(SelectSeason)%>%
  filter(Length > 175)%>%
  filter( Length<1000)%>%
  mutate(selectSeason=if_else(Water.Section=="BigHole,JerryCreek",'Spring',selectSeason))%>%
  filter(SampleSeason==selectSeason)%>%
  ungroup()

WatersIn<-unique(RBT_extractIN%>%mutate(Water.Section2=paste0(Water,".",Section))%>%pull(Water.Section))%>%sort()

YearsIn<-unique(RBT_extractIN%>%group_by(Water.Section)%>%select(Year))

##Read in and summarise otolith data
ALLotos<-read.csv(here::here('Data','ALLotos.csv'))%>%
  mutate(Stream.Section=paste0(Stream,",",Section))%>%filter(Age<5 & Age>1)%>%
  filter(Species=="RB")
ALLotosSUMsection<-ALLotos%>%
  mutate(Stream.Section=paste0(Stream,",",Section))%>%filter(Age<5 & Age>1)%>%
  group_by(Stream.Section,Age)%>%
  summarize(mean_mm=mean(mean_mm,na.rm=T),
            min_mm=mean(min_mm,na.rm=T),
            max_mm=mean(max_mm,na.rm=T))%>%
  mutate(min_mm=if_else(is.na(min_mm),0,min_mm),
         max_mm=if_else(is.na(max_mm),0,max_mm))
ALLotosSUMstream<-ALLotos%>%
  mutate(Stream.Section=paste0(Stream,",",Section))%>%filter(Age<5 & Age>1)%>%
  group_by(Stream,Age)%>%
  summarize(mean_mm=mean(mean_mm,na.rm=T),
            min_mm=mean(min_mm,na.rm=T),
            max_mm=mean(max_mm,na.rm=T))%>%
  mutate(min_mm=if_else(is.na(min_mm),0,min_mm),
         max_mm=if_else(is.na(max_mm),0,max_mm))
ALLotosSUMstate<-ALLotos%>%
  mutate(Stream.Section=paste0(Stream,",",Section))%>%filter(Age<5 & Age>1)%>%
  group_by(Age)%>%
  summarize(mean_mm=mean(mean_mm,na.rm=T),
            min_mm=mean(min_mm,na.rm=T),
            max_mm=mean(max_mm,na.rm=T))%>%
  mutate(min_mm=if_else(is.na(min_mm),0,min_mm),
         max_mm=if_else(is.na(max_mm),0,max_mm))

specificOTOsections<-ALLotos%>%pull(Stream.Section)%>%unique()
specificOTOstream<-ALLotos%>%pull(Stream)%>%unique()
specificOTOsections[which(!specificOTOsections%in%WatersIn)]


mix_OutWater<-tibble(data.frame(Clean=WatersIn))
mix_OutWater$data<-list("")
mix_OutWater$model<-list("")
mix_OutWater$season<-list("")

for(ww in 1:length(unique(WatersIn))){
  tALL1<-Sys.time()
  # ww<-5
  lengthsIn1<-RBT_extractIN%>%filter(Water.Section==WatersIn[ww] & !is.na(Length))
  WaterSection<-lengthsIn1%>%pull(Water.Section)%>%unique()
  Water<-lengthsIn1%>%pull(Water)%>%unique()
  print(WaterSection)
  Season<-lengthsIn1%>%pull(SampleSeason)%>%unique()
  yearset <- unique(lengthsIn1$Year)
  
  lengthsIn2 <- lengthsIn1 %>% filter(Year %in% yearset)
  
  y <- lengthsIn2%>%pull(Length)
  Year <- lengthsIn2 %>% pull(Year) %>% as.factor() %>% as.numeric()
  realYear <- lengthsIn2 %>% pull(Year) 
  
  N <- length(y)
  K<-3
  nYear = length(unique(Year))
  
  if(Season=="Spring"){
    jags.data = list(
      y = y,
      N = N,
      Year = Year,
      nYear = nYear,
      K = K,
      mu.prior=c(250,325,450),
      mu.prior.min=c(200,275,375),
      mu.prior.max=c(275,375,525),
      pi_prior = rep(1,K)
    )}else{
      jags.data = list(
        y = y,
        N = N,
        Year = Year,
        nYear = nYear,
        K = K,
        mu.prior=c(275,375,500),
        mu.prior.min=c(225,325,400),
        mu.prior.max=c(325,400,600),
        pi_prior = rep(1,K)
      )}
  
  if(WaterSection%in%specificOTOsections){
    otosIN<-ALLotosSUMsection%>%filter(Stream.Section==WaterSection)
    jags.data = list(
      y = y,
      N = N,
      Year = Year,
      nYear = nYear,
      K = K,
      mu.prior=c(otosIN$mean_mm[1],otosIN$mean_mm[2],otosIN$mean_mm[3]),
      mu.prior.min= if(otosIN$max_mm[1]==0){
        c(otosIN$mean_mm[1]-50,otosIN$mean_mm[2]-50,otosIN$mean_mm[3]-50)}else{
          c(otosIN$min_mm[1],otosIN$min_mm[2],otosIN$min_mm[3])
        },
      mu.prior.max= if(otosIN$max_mm[1]==0){
        c(otosIN$mean_mm[1]+50,otosIN$mean_mm[2]+50,otosIN$mean_mm[3]+50)}else{
          c(otosIN$max_mm[1],otosIN$max_mm[2],otosIN$max_mm[3])}
      ,
      pi_prior = rep(1,K)
    )
  }else{jags.data<-jags.data}
  
  if(!WaterSection%in%specificOTOsections& Water%in%specificOTOstream){
    otosIN<-ALLotosSUMstream%>%filter(Stream==Water)
    jags.data = list(
      y = y,
      N = N,
      Year = Year,
      nYear = nYear,
      K = K,
      mu.prior=c(otosIN$mean_mm[1],otosIN$mean_mm[2],otosIN$mean_mm[3]),
      mu.prior.min= if(otosIN$max_mm[1]==0){
        c(otosIN$mean_mm[1]-50,otosIN$mean_mm[2]-50,otosIN$mean_mm[3]-50)}else{
          c(otosIN$min_mm[1],otosIN$min_mm[2],otosIN$min_mm[3])
        },
      mu.prior.max= if(otosIN$max_mm[1]==0){
        c(otosIN$mean_mm[1]+50,otosIN$mean_mm[2]+50,otosIN$mean_mm[3]+50)}else{
          c(otosIN$max_mm[1],otosIN$max_mm[2],otosIN$max_mm[3])}
      ,
      pi_prior = rep(1,K)
    )
  }else{jags.data<-jags.data}
  
  
  jags.params <- c("mu","pi",'clust','mu.hyper','lsigma.hyper','lsigma','lmu.A.sigma') #,'lsigma.A.sigma'
  
  mod_nmix_melHeir <-  R2jags::jags.parallel(data = jags.data,
                                             parameters.to.save = jags.params,
                                             model.file = modelScript.name,
                                             n.iter = 6000,
                                             n.thin = 10,
                                             n.burnin=6000/4,
                                             n.chains = 3)
  
  
  mix_OutWater$Clean[[ww]]<-WatersIn[ww]
  mix_OutWater$data[[ww]]<-(tibble(lengthsIn2))
  mix_OutWater$season[[ww]]<-Season
  
  mediansY<-matrix(ncol=5,nrow=nYear*K,byrow=T)
  mediansY[,1]<-rep(c(1:K),times=nYear)
  mediansY[,2]<-mod_nmix_melHeir$BUGSoutput$median$mu[1:3,]
  mediansY[,3]<-mod_nmix_melHeir$BUGSoutput$median$pi[1:3,]
  mediansY[,4]<-exp(mod_nmix_melHeir$BUGSoutput$median$lsigma[1:3,])
  mediansY[,5]<-rep(unique(lengthsIn2 %>% pull(Year)),each=K)
  colnames(mediansY)<-c("K","mu","pi","lsigma","Year")
  
  
  
  mediansA<-matrix(ncol=K,nrow=3,byrow=T)
  mediansA[1,]<-mod_nmix_melHeir$BUGSoutput$median$mu.hyper
  mediansA[2,]<-exp(mod_nmix_melHeir$BUGSoutput$median$lsigma.hyper)
  mediansA[3,]<-exp(mod_nmix_melHeir$BUGSoutput$median$lmu.A.sigma)
  rownames(mediansA)<-c("mu.hyper","sigma.hyper","sigma.sigma.hyper")
  models<-list(mediansY,mediansA,mod_nmix_melHeir)
  
  mix_OutWater$model[[ww]]<-models
  tALL2<-Sys.time()
  totalTimeAll<-tALL2-tALL1
  print(totalTimeAll)
}

#Extract age for each fish
BLT.Age.heir <- lapply(mix_OutWater %>%   
                         pull(Clean), FUN = function(SS){
                           # SS <-"Madison,Varney"
                           Model <- mix_OutWater %>% filter(Clean==SS) %>% 
                             pull(model) 
                           
                           DAT <- mix_OutWater %>% filter(Clean==SS) %>% 
                             pull(data) %>% 
                             .[[1]] %>% 
                             mutate(Age = Model[[1]][[3]]$BUGSoutput$sims.list$clust %>% 
                                      apply(2,FUN=median)-1) %>%
                             mutate(Age = ceiling(Age)) %>% 
                             mutate(Clean = SS)
                           return(DAT)
                           
                         }) %>% bind_rows()

AllYearsDist <- BLT.Age.heir %>% 
  group_by(Clean,Year,Age) %>%
  summarise(N=length(Length),MeanLength = mean(Length),SDLength = sd(Length)) %>%
  mutate(PropAge = N/sum(N)) %>% 
  select(Clean, Year, Age, N, MeanLength, SDLength, PropAge)

Sections<- BLT.Age.heir %>%
  pull(Clean) %>%
  unique() %>%
  sort()

#Make density plot of anuual size at age for each section.
BT.Mix.Plot <-lapply(Sections,FUN = function(SS){
  # SS<-"BigHole,Melrose"
  if(count(BLT.Age.heir%>%filter(Clean==SS))>100){
    Sec.Data <-BLT.Age.heir %>% filter(Clean==SS)
    
    UniYears <- Sec.Data %>% 
      pull(Year) %>% 
      unique() %>% 
      sort()
    
    lapply(UniYears,FUN = function(yy){
      # yy=1995
      print(paste(SS,yy,sep = "."))
      YearDat <- Sec.Data %>% filter(Year == yy)
      YearDist <- AllYearsDist %>% filter(Clean == SS,Year == yy)
      YearDist$SDLength[is.na(YearDist$SDLength)] <- mean(AllYearsDist %>% filter(Clean == SS) %>% pull(SDLength),na.rm = T)
      
      Dx <- lapply(1:nrow(YearDist),FUN=function(x){
        Dx <- rnorm(10000,mean=YearDist$MeanLength[x],sd=YearDist$SDLength[x]) %>%
          density()
        Dx$y <- Dx$y*YearDist$PropAge[x]
        return(Dx)
      })
      
      AllSimDens <- lapply(1:nrow(YearDist),FUN=function(x){
        Rx <- rnorm(YearDist$N[x]*100,mean=YearDist$MeanLength[x],sd=YearDist$SDLength[x])
        return(Rx)
      }) %>% unlist() %>% density()
      
      plotYearDist <- YearDat %>% 
        mutate(Age=as.factor(Age+2))%>%
        ggplot(aes(x=Length)) +
        geom_histogram(aes(y=..density..),binwidth=5,colour=1) + 
        geom_line(data=data.frame(XX=AllSimDens$x,YY=AllSimDens$y),aes(x=XX,y=YY),col='darkgreen',linewidth=1) +
        # scale_fill_manual(values=c("white","lightgray","darkgray"))+
        theme_minimal(base_size=6) + 
        labs(title=paste(SS,yy, sep='.'))
      for(i in 1:length(Dx)){
        plotYearDist <- plotYearDist +
          geom_line(data=data.frame(XX=Dx[[i]]$x,YY=Dx[[i]]$y),aes(x=XX,y=YY),col='red',linewidth=1)
      }
      return(plotYearDist)
    })
  }})


pdf("ALL_LFQMix_heir_density.pdf",width=6,height=8)
lapply(1:length(unique(BLT.Age.heir$Clean)),FUN=function(x){
  gridExtra::marrangeGrob(grobs = BT.Mix.Plot[[x]], nrow = 4, ncol =2)
})
dev.off()


