rm(list=ls()) 

library(spatstat)
library(gstat)
library(RColorBrewer)
library(readr)
library(dplyr)
library(tidyr)
library(corrplot)


data<-read.csv("Coaseian Tuna Data/WCPFC_Yearbook2014/YB_WCP_CA.csv")
which_catch = which(grepl("_mt", colnames(data)))

ll_data<-data%>%
  gather('bad_name','catch',which_catch)%>%
  ungroup()%>%
  mutate(species = gsub('_.*$','',bad_name),
         set_type = gsub('.*_','',bad_name)) %>%
  select(-bad_name)%>%
  ungroup()%>%
  group_by(yy,flag,gear,species)%>%
  summarise(sp_catch=sum(catch,na.rm=TRUE))%>%
  ungroup()%>%
  select(yy,gear,flag,species,sp_catch)%>%
  filter(gear=="L")%>%
  filter(species=="bet")%>%
  filter(yy>'1999')%>%
  arrange(flag,yy)%>%
  select(yy,flag,sp_catch)

ps_data<-read.csv("Coaseian Tuna Data/WCPFC_Yearbook2014/YB_WCP_CA.csv")%>%
  gather('bad_name','catch',which_catch)%>%
  ungroup()%>%
  mutate(species = gsub('_.*$','',bad_name),
         set_type = gsub('.*_','',bad_name)) %>%
  select(-bad_name)%>%
  ungroup()%>%
  group_by(yy,flag,gear,species)%>%
  summarise(sp_catch=sum(catch,na.rm=TRUE))%>%
  ungroup()%>%
  select(yy,gear,flag,species,sp_catch)%>%
  filter(gear=="S")%>%
  filter(species=="bet")%>%
  filter(yy>'1999')%>%
  arrange(flag,yy)%>%
  select(yy,flag,sp_catch)



ps_cor<-spread(ps_data,flag,sp_catch)
#ps_cor[is.na(ps_cor)]<-0

ll_cor<-spread(ll_data,flag,sp_catch)
#ll_cor[is.na(ll_cor)]<-0
ll_cor<-as.matrix(ll_cor[,-1])
ps_cor<-as.matrix(ps_cor[,-1])

correlations<-cor(ll_cor,ps_cor,method="pearson",use="pairwise.complete.obs")
correlations<-correlations[,-1]
write.csv(correlations,"Catch Correlations/BET PS_LL correlation coefficients.csv")

png(height=1600, width=1200, file="Catch Correlations/BET LL_PS correlation matrix.png")

corrplot(correlations,na.label=" ", tl.col = "black",method='color',title="BET total catch 2000-2014 (LL / PS correlations)",
         mar=c(0,0,2,0))
 
#vertical axis is longline and horizonal axis is PS
dev.off()


