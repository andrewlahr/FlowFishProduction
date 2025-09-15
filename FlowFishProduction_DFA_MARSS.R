library(MuMIn)
library(MARSS)
library(tidyverse)

load(here::here('../LL/data/ComputeFlowAndTemps071725.Rdata'))
dir_name2<-("JAGS_PVA/ModelOutput/")
# tempDat<-SummerMeanTempTable
library(tidyverse)

sectionsIN<-read.csv("AllSectionsSpecies.csv")

species<-c('LL','RBT')


#read in abundance data from model output
densOUT<-lapply(species,function(sp){
  # sp='LL'
  print(sp)

  if(sp=='LL'){
    AllSectionsIN<-sectionsIN%>%filter(species=='Brown Trout')%>%pull(Stream.Section)}else{
      AllSectionsIN<-sectionsIN%>%filter(species=='Rainbow Trout')%>%pull(Stream.Section)
    }
  densities<-lapply(AllSectionsIN,function(SS){
    # SS<-"BigHole.JerryCreek"
    print(SS)
    # sp<-'RBT'
    
    Stream<-strsplit(SS,split = '\\.')[[1]][1]
    Section<-strsplit(SS,split = '\\.')[[1]][2]
    Density<-read.csv(here::here('..',sp,dir_name2,"csvs_quadratic",paste0(Stream,Section,"allParams.csv")))%>%
      filter(topMod=='yes'& name=="Estimated NALL")%>%
      dplyr::select(N50,Year,Stream.Section)%>%mutate(Species=sp)
    
    
    
    return(Density)
  })%>%bind_rows()
  
  
  return(densities)
})%>%bind_rows()

#Create table of brown trout abundance

EstimatedBT <- densOUT %>% filter(Species =='LL')

BT <- xtabs(N50 ~ Stream.Section + Year,EstimatedBT)
BT[BT==0] <- NA

BTz <- BT %>% Zscore()%>% as.matrix()
  
###PCA

BTz

dfa1_BT_e <- MARSS(BTz,model=list(m=1,R='diagonal and equal'),form='dfa')
dfa1_BT_ue <- MARSS(BTz,model=list(m=1,R='diagonal and unequal'),form='dfa')

AICc(dfa1_BT_e)
AICc(dfa1_BT_ue)

TopMod1 <- dfa1_BT_ue

data.frame(Site=row.names(BTz),Loading = coef(TopMod1)$Z) %>% arrange(Loading)

FitVec <- coef(TopMod1,type='matrix')$Z %*% TopMod1$states[1,] %>% as.vector()

summary(lm(as.vector(BTz) ~ FitVec))

#Rainbow trout

EstimatedRB <- densOUT %>% filter(Species =='RBT')

RB <- xtabs(N50 ~ Stream.Section + Year,EstimatedRB)
RB[RB==0] <- NA

RBz <- RB %>% apply(1,FUN=scale) %>% t() %>% as.matrix()

###PCA

RBz

dfa1_RB_e <- MARSS(RBz,model=list(m=1,R='diagonal and equal'),form='dfa')
dfa1_RB_ue <- MARSS(RBz,model=list(m=1,R='diagonal and unequal'),form='dfa')

AICc(dfa1_RB_e)
AICc(dfa1_RB_ue)

TopMod1RB <- dfa1_RB_e

###
SummerFlow <- densOUT%>%filter(Year<=2024)%>%left_join(
  
  read.csv('allN_LL_RBT.csv',header=T,stringsAsFactors = F)%>%
    dplyr::select(Year,Stream.Section,Species,SummerFlow,WinterFlow)%>%
    mutate(Species=if_else(Species=="Brown Trout","LL","RBT"))
  
)%>%filter(!duplicated(.))
SumFlow <- xtabs(SummerFlow ~ Stream.Section + Year,data=EstimatedRB)
SumFlow[SumFlow==0] <- NA
SumFlowZ <- SumFlow %>%
  apply(1,FUN=scale) %>%
  t()
TotalSumFlow <- SumFlowZ %>% colMeans(na.rm=T)
SumFlowRoll <- zoo::rollmean(TotalSumFlow,k=3,align='right',na.pad=T)
plot(SumFlowRoll,type='l')


data.frame(trend=TopMod1$states[1,],se_Low=TopMod1$states[1,]-TopMod1$states.se[1,],
           se_High=TopMod1$states[1,]+TopMod1$states.se[1,],
           Species='Brown Trout',Year=1980:2024)%>%
  rbind(data.frame(trend=TopMod1RB$states[1,],se_Low=TopMod1RB$states[1,]-TopMod1RB$states.se[1,],
                   se_High=TopMod1RB$states[1,]+TopMod1RB$states.se[1,],
                   Species='Rainbow Trout',Year=1980:2024),
        data.frame(trend=SumFlowRoll*4,
                   se_Low=SumFlowRoll*4-seFlowRoll*4,
                   se_High=SumFlowRoll*4+seFlowRoll*4,
                   Species='Summer Flow',Year=1980:2024))%>%
  mutate(Species=as.factor(Species))%>%
  ggplot(aes(x=Year,y=trend,color=Species))+
  geom_line()+
  geom_ribbon(aes(x=Year,ymin=se_Low,ymax=se_High,fill=Species),alpha=0.2,color=NA,show.legend=F)+
  theme_classic()+
  scale_color_manual(values=c('chocolate2','yellowgreen','royalblue'),name='')+
  scale_fill_manual(values=c('chocolate2','yellowgreen','royalblue'))+
  ylab('Trend')


BTloadings <- coef(TopMod1, type = "matrix")$Z
BTloadingsOUT<-data.frame(Stream.Section=rownames(BTloadings),BrownTrout=BTloadings[,1])
RBloadings<- coef(TopMod1RB, type = "matrix")$Z
RBloadingsOUT<-data.frame(Stream.Section=rownames(RBloadings),RainbowTrout=RBloadings[,1])
rownames(RBloadingsOUT)<-NULL
rownames(BTloadingsOUT)<-NULL

lambsOUT<-lapply(species,function(sp){
  print(sp)
  # sp='LL'
  print(sp)
  

  LL_kSections<-sectionsIN%>%filter(species=='Brown Trout')%>%pull(Stream.Section)
  
  RBT_kSections<- AllSectionsIN<-sectionsIN%>%filter(species=='Rainbow Trout')%>%pull(Stream.Section)
    
  densities<-lapply(AllSectionsIN,function(SS){
    print(SS)
    if(sp=="RBT"){
      if( SS%in%RBT_kSections){
        Stream<-strsplit(SS,split = '\\.')[[1]][1]
        Section<-strsplit(SS,split = '\\.')[[1]][2]
        Density<-read.csv(here::here('..',sp,dir_name2,"csvs_quadratic",paste0(Stream,Section,"allParams.csv")))%>%
          filter(topMod=='yes'& name=="Estimated NALL")%>%
          dplyr::select(N50,Year,Stream.Section)%>%mutate(Species=sp)
        Density2<- Density%>%filter(Year<=2023)%>%pull(N50)
        plot(Density2)
        mod<-MARSS::MARSS(log(Density2),model=list(R='zero',U='unequal',Q='unconstrained'))
        lnLam<-coef(mod)$U
        lams<-data.frame(lnLambda=lnLam,
                         Stream.Section=SS,
                         species=sp)
        rownames(lams)<-NULL
      }else{
        
        lams<-data.frame(lnLambda=NA,
                         Stream.Section=SS,
                         species=sp)
        rownames(lams)<-NULL
      }
    }else{
      if( SS%in%LL_kSections){
        Stream<-strsplit(SS,split = '\\.')[[1]][1]
        Section<-strsplit(SS,split = '\\.')[[1]][2]
        Density<-read.csv(here::here('..',sp,dir_name2,"csvs_quadratic",paste0(Stream,Section,"allParams.csv")))%>%
          filter(topMod=='yes'& name=="Estimated NALL")%>%
          dplyr::select(N50,Year,Stream.Section)%>%mutate(Species=sp)
        Density2<- Density%>%filter(Year<=2023)%>%pull(N50)
        plot(Density2)
        mod<-MARSS::MARSS(log(Density2),model=list(R='zero',U='unequal',Q='unconstrained'))
        lnLam<-coef(mod)$U
        lams<-data.frame(lnLambda=lnLam,
                         Stream.Section=SS,
                         species=sp)
        rownames(lams)<-NULL
      }else{
        
        lams<-data.frame(lnLambda=NA,
                         Stream.Section=SS,
                         species=sp)
        rownames(lams)<-NULL
      }
    }
    
    
    return(lams)
  })%>%bind_rows()
  
  
  return(densities)
})%>%bind_rows()

streamsections<-lambsOUT%>%
  pull(Stream.Section)%>%unique()

BTlambdas<-lambsOUT%>%filter(species=='LL')%>%dplyr::select(Stream.Section,lnLambda)%>%dplyr::rename(LLlnLambda=lnLambda)
RBlambdas<-lambsOUT%>%filter(species=='RBT')%>%dplyr::select(Stream.Section,lnLambda)%>%dplyr::rename(RBTlnLambda=lnLambda)


write.csv(left_join(BTloadingsOUT,BTlambdas)%>%left_join(left_join(RBloadingsOUT,RBlambdas)),'lambdasAndLoadings.csv')
