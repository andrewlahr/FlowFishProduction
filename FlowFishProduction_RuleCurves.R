
##Stable Pop Plot####
if(global=="yes"){
  
  mod<- readRDS(here::here("JAGS_PVA","ModelFits",paste0(SS,"_IndicatorVarSel_TopMod.rds")))
  modIN<-mod$BUGSoutput$sims.list
}else{
  mod<- readRDS(here::here("JAGS_PVA","ModelFits",paste0(SS,"_IndicatorVarSel_topModSummer.rds")))
  modIN<-mod$BUGSoutput$sims.list
  
}

Stock<-seq(0, 5000)
Weight4<-median(W4_In,na.rm=T)
Weight3<-median(W3_In,na.rm=T)

old<-8
SSB<-0
SSB<-SSB+((Stock-(Stock*Survival))*Weight3) #age 3s stat dont go on to be 4s

for(i in 1:old){
  SSB<-SSB+(((Stock*(Survival^i))-(Stock*(Survival^(i+1))))*Weight4)
  
}
SSB<-SSB+((Stock*(Survival^(old+1)))*Weight4)

WeightR<-median(W2_In,na.rm=T)


flowScenariosSummer<-read.csv("Data/FlowScenarios.csv")%>%
  mutate(PercentChange=1+(PercentChange/100))%>%
  filter(FlowPeriod=="Summer")

Flows<-quantile(covarLagIn1Real,seq(0.01,1,length=100),na.rm=T)
flow1 <- Flows%>%log()%>%Zscore()
Flows2<-0.000001

trueMeanFlow<-(covarLagIn1Real)%>%log()%>%mean(na.rm=T)%>%exp()

meanIndex<-as.numeric(gsub(".*?([0-9]+).*","\\1",names(Flows[which.min(abs(Flows-trueMeanFlow))])))
for(x in 1:nrow(flowScenariosSummer)){

  flowScenariosSummer$trueMeanFlow[x]<-Flows[which.min(abs(as.numeric(gsub('%','',names(Flows)))-meanIndex))]*flowScenariosSummer$PercentChange[x]
  flowScenariosSummer$ZMeanFlow[x]<-flow1[which.min(abs(as.numeric(gsub('%','',names(flow1)))-meanIndex))]*flowScenariosSummer$PercentChange[x]
  flowScenariosSummer$meanIndex[x]<-names(Flows[which.min(abs(Flows-flowScenariosSummer$trueMeanFlow[x]))])
}

futureFLowIndexUnique<-c(flowScenariosSummer$meanIndex)

Sigma<-modIN$sigmaR
logAlpha<-modIN$la0
Beta<-modIN$b
bRec1<-modIN$b1r
bRec12<-modIN$b12r

lSurvival<-modIN$lSurv0
bSurv1<-modIN$b1S
bSurv12<-modIN$b12S
bSurv2<-modIN$b6S

mortality<-1-Survival

if( length(bRec12)==0){
  bRec12<-rep(0,length(logAlpha))}else{
    bRec12<-bRec12}

if(length(bSurv12)==0){
  bSurv12<-rep(0,length(logAlpha))}else{
    bSurv12<-bSurv12}

library(parallel)
library(foreach)
library(doParallel)

sampleIN<-sample(c(1:length(logAlpha)),size=1000)
trueMeanFlow<-(covarLagIn1Real)%>%log()%>%mean()%>%exp()

meanIndex<-as.numeric(gsub(".*?([0-9]+).*","\\1",names(Flows[which.min(abs(Flows-trueMeanFlow))])))
flow3<-c(flow1[which(names(flow1)%in% c(futureFLowIndexUnique,'5%','25%',paste0(meanIndex,'%'),'75%','90%'))])

cl1<-makeCluster(5) #make cluster
registerDoParallel(cl1) #register cluster
tALL1<-Sys.time()

StockPlot<-0
StockPlot<-StockPlot+(((Stock/Survival)-Stock)) #age 2s that don't go on to be 3s
StockPlot<-StockPlot+(((Stock-(Stock*Survival)))) #age 3s that dont go on to be 4s

for(i in 1:old){
  StockPlot<-StockPlot+(((Stock*(Survival^i))-(Stock*(Survival^(i+1)))))
  
}
StockPlot<-StockPlot+((Stock*(Survival^(old+1)))) 

test<-foreach(i=1:length(flow3),.packages = c('dplyr'))%dopar%{
  MaxSurplusOUT<-list()
  StockAtMaxSurplus<-list()
  CarryingCapacityOUT<-list()
  
  Ad<-list()
  SurplusProduction<-matrix(nrow=length(StockPlot),ncol=length(sampleIN))
  CarryingCapacity<-list()
  count<-1
  for(x in sampleIN){
    Rs <- lapply(1:length(SSB),FUN=function(S){
    
      R <- exp(rnorm(1000,logAlpha[x] + log(SSB[S]) - (Beta[x]*SSB[S]) + (bRec1[x]*flow3[i])  +(bRec12[x]*(flow3[i]^2))-((Sigma[x]^2)/2),Sigma[x]))
      lSurv <- lSurvival[x] + (bSurv1[x]*flow3[i]) +(bSurv12[x]*(flow3[i]^2))  + (bSurv2[x]*((Stock[S] + (R*Survival) + R)/100))
      
      Survi <- (exp(lSurv)/(1+exp(lSurv)))
      
      ALL<-0
      ALL<-(R-(R*Survi))
      
      for(i in 1:old){
        ALL<-ALL+(((R*(Survi^i))-(R*(Survi^(i+1)))))
      }
      
      ALL<-ALL+((R*(Survi^(old+1))))
      
      return(median(ALL,na.rm=T)*(1/(1-Survival)))
    }) %>% unlist()

    Ad[[count]] <- Rs
    stablePopLine<-data.frame(Stock=StockPlot,N=Ad[[count]])
    
    SurplusProduction[,count]<-c(stablePopLine$N-stablePopLine$Stock)
    
    CarryingCapacity[[count]]<-stablePopLine$Stock[which(stablePopLine$Stock>0)][which.min(abs(stablePopLine$N[which(stablePopLine$Stock>0)]-stablePopLine$Stock[which(stablePopLine$Stock>0)]))]
    
    count<-count+1
  }
  

  MaxSurplusOUT<-SurplusProduction[which.max(SurplusProduction%>%apply(.,1,FUN=median)),] 
  StockAtMaxSurplus<-StockPlot[which.max(SurplusProduction%>%apply(.,1,FUN=median))]
  CarryingCapacityOUT<-CarryingCapacity%>%unlist()
  
  test2<-Ad%>%bind_cols()%>%apply(1,FUN=quantile,c(0.25,0.5,0.75),na.rm=T)
  
  testOUT<-lapply(1:length(SSB),FUN=function(S){
    return(test2[,S])
  })%>%bind_rows()
  
  colnames(testOUT)<-c('N10','N50','N90')
  
  return(list(MaxSurplusOUT,StockAtMaxSurplus,CarryingCapacityOUT,testOUT))
  
}

tALL2<-Sys.time()
totalTimeAll<-tALL2-tALL1
print(totalTimeAll)
stopCluster(cl1) #turn off clusters

realflows<-c(Flows[which(names(Flows)%in% c(futureFLowIndexUnique,'5%','25%',paste0(meanIndex,'%'),'75%','90%'))])#%>%unique()

KandMaxSur<-lapply(1:length(flow3),function(ff){
  meanFF<-which(names(flow3)==paste0(meanIndex,'%'))
  MaxSurplusProduction<-c(test[[ff]][[1]]%>%quantile(c(0.25,0.5,0.75),na.rm=T))
  StockAtMaxSurplus<-test[[ff]][[2]]
  CarryingCapacity<- c(test[[ff]][[3]]%>%quantile(c(0.25,0.5,0.75),na.rm=T))
  deltaK<- c(((test[[ff]][[3]]-test[[meanFF]][[3]])/test[[meanFF]][[3]])%>%quantile(c(0.25,0.5,0.75),na.rm=T))
  
  deltaMS<-c(((test[[ff]][[1]]-test[[meanFF]][[1]])/test[[meanFF]][[1]])%>%quantile(c(0.25,0.5,0.75),na.rm=T))
  flow<-flow3[ff]
  
  return(data.frame(MaxSurplusProduction=MaxSurplusProduction,StockAtMaxSurplus=StockAtMaxSurplus,
                    CarryingCapacity=CarryingCapacity,realflow=realflows[ff],deltaK=deltaK,deltaMS=deltaMS,flow=flow,quantile=c('25','50','75'),percentile=names(flow3[ff])))
})%>%bind_rows()%>%mutate(Stream.Section=SS)

rownames(KandMaxSur)<-NULL

write.csv(KandMaxSur,here::here(dir_name2,'csvs_quadratic',paste0(SS,"_KandMaxSur.csv")))

stablePopOUT<-lapply(1:length(flow3),function(ff){
  
  return( data.frame(test[[ff]][[4]])%>%
            mutate(realflow=realflows[ff],flow=flow3[ff],percentile=names(flow3[ff])))
  
})%>%bind_rows()%>%mutate(Stream.Section=SS,StockPlot=rep(StockPlot,length(flow3)))%>%left_join(  data.frame(
  K50=KandMaxSur$CarryingCapacity[which(KandMaxSur$quantile=='50')],
  K25=KandMaxSur$CarryingCapacity[which(KandMaxSur$quantile=='25')],
  K75=KandMaxSur$CarryingCapacity[which(KandMaxSur$quantile=='75')],
  MS50=KandMaxSur$MaxSurplusProduction[which(KandMaxSur$quantile=='50')],
  MS25=KandMaxSur$MaxSurplusProduction[which(KandMaxSur$quantile=='25')],
  MS75=KandMaxSur$MaxSurplusProduction[which(KandMaxSur$quantile=='75')],
  percentile=names(flow3)))

write.csv(stablePopOUT,here::here(dir_name2,'csvs_quadratic',paste0(SS,"_stablePopOUTBig2.csv")))

futureKforJoin<- data.frame(ForecastPeriod=rep(flowScenariosSummer$ForecastPeriod,1),FlowPeriod=rep(flowScenariosSummer$FlowPeriod,1),
                            HighLow=rep(flowScenariosSummer$HighLow,1),Index=rep(50,6)
)%>%arrange(ForecastPeriod,HighLow)
futureKforJoin$percentile<-NA
futureKforJoin$percentile[which(futureKforJoin$Index==50)]<-c(flowScenariosSummer$meanIndex)


write.csv(futureKforJoin%>%left_join(stablePopOUT),here::here(dir_name2,'csvs_quadratic',paste0(SS,"_stablePopOUTFutureKBig.csv")))

flowsNeeded<-c('5%','25%',paste0(meanIndex,"%"),'75%','90%')
Stream<-strsplit(SS,split = '\\.')[[1]][1]
Section<-strsplit(SS,split = '\\.')[[1]][2]
datIN<-read.csv(here::here(dir_name2,'csvs_quadratic',paste0(SS,"_stablePopOUTBig2.csv")))%>%
  dplyr::select(N50,N10,N90,StockPlot,MS50,K50,percentile,realflow,flow)%>%unique()%>%
  filter(percentile%in%flowsNeeded)%>%
  mutate(Stream.Section=SS)
flowKey<-data.frame(flowName= c("Severe drought - 1/20 yr",'Moderate drought - 1/4 yr','Mean','75th Percentile','90th Percentile'),
                    flowsNeeded=c('5%','25%',paste0(meanIndex,"%"),'75%','90%'))

for(i in 1:nrow(datIN)){
  perIN<- datIN$percentile[i]
  perOUT<-NA
  for(n in 1:nrow(flowKey)){
    perOUT<-if_else(perIN==flowKey$flowsNeeded[n],flowKey$flowName[n],perOUT)
  }
  datIN$flowName[i]<-perOUT
}

N50s<-datIN%>%filter(percentile=='90%')%>%pull(N50)


meanK<-datIN%>%select(K50,MS50,flowName,percentile)%>%unique()

meanPop<-round(read.csv(here::here(dir_name2,"csvs_quadratic",paste0(Stream,Section,"WeightedForeCast.csv")))%>%
                 # filter(topMod=='yes')%>%
                 filter(name=="Estimated NALL")%>%pull(N50w)%>%log()%>%mean%>%exp())

newKs<-datIN%>%
  group_by(flowName)%>%
  reframe(N=c(0,gcplyr::smooth_data(x=StockPlot,y=N50,sm_method='smooth.spline',warn_ungrouped=F)),
          StockPlot=c(0,StockPlot))%>%
  filter(StockPlot>0 & N>0)%>%
  mutate(difference=StockPlot-N)%>%
  group_by(flowName)%>%
  slice_min(abs(difference))%>%
  mutate(K50new=N)%>%select(K50new,flowName)

xlim<-newKs%>%pull(K50new)%>%max()
ylim<-max(datIN$N50[which(datIN$StockPlot<=xlim)])


plot1<-datIN%>%
  group_by(flowName)%>%
  reframe(N=c(0,gcplyr::smooth_data(x=StockPlot,y=N50,sm_method='smooth.spline',warn_ungrouped=F)),
          StockPlot=c(0,StockPlot))%>%
  left_join(newKs)%>%
  mutate(flowName=factor(flowName,levels=rev(c("Severe drought - 1/20 yr",'Moderate drought - 1/4 yr','Mean','75th Percentile','90th Percentile'))))%>%
  ggplot(aes(x=StockPlot,y=N,color=flowName))+
  geom_abline(intercept=0,slope=1,lty=2,lwd=1)+
  geom_linerange(mapping=aes(x=K50new,ymin=0,ymax=K50new,color=flowName),lwd=0.75,show.legend = FALSE)+
  geom_line(,lwd=1)+
  geom_text(aes(x=K50new, y=0, label = signif(K50new,2)),vjust = 1,size=3,angle = -45,hjust = 0,show.legend = FALSE) +
  coord_cartesian(clip = "off") +
  scale_color_manual(name="",values=c("dodgerblue","lightblue","grey","darkorange","red"),
                     labels=c( paste0("90th Percentile \n(",round(Flows[[90]],0)," cfs)"),
                               paste0("75th Percentile \n(",round(Flows[[75]],0)," cfs)"),
                               paste0("Mean \n(",round(Flows[[meanIndex]],0)," cfs)"),
                               paste0("Moderate drought - 1/4 yr \n(",round(Flows[[25]],0)," cfs)"),
                               paste0("Extreme drought - 1/20 yr \n(",round(Flows[[5]],0)," cfs)")))+
  scale_fill_manual(name="",values=c("dodgerblue","lightblue","grey","darkorange","red"),
                    labels=c(paste0("90th Percentile \n(",round(Flows[[90]],0)," cfs)"),
                             paste0("75th Percentile \n(",round(Flows[[75]],0)," cfs)"),
                             paste0("Mean \n(",round(Flows[[meanIndex]],0)," cfs)"),
                             paste0("Moderate drought - 1/4 yr \n(",round(Flows[[25]],0)," cfs)"),
                             paste0("Extreme drought - 1/20 yr \n(",round(Flows[[5]],0)," cfs)")))+
  theme_classic()+
  scale_y_continuous(expand = c(0,NA),limits=c(0,1*ylim))+
  scale_x_continuous(expand = c(0,0),limits=c(0,1.1*xlim))+
  labs(x="N[t]",y="N[t+4]",title=paste0(SS,", RBT","\n Mean Pop Est = ",meanPop))+
  theme(legend.position = "right",legend.text=element_text(size=7),plot.title = element_text(size =10))
