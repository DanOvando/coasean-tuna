library(raster)
library(rgdal)


months <- c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct', 'Nov', 'Dec')

ps_fun <- function(m){
  
  ps_dat = read.csv(paste('Coaseian Tuna Data/PS_5Y/PS_5Y/PS_5Y_',m,'.csv', sep = ''), stringsAsFactors = F)
  
  ps_dat$month = m
  
  return(ps_dat)
}


ps_5y = bind_rows(lapply(months, ps_fun))
ps_5y<-ps_5y[ps_5y$LON5_D>-126,]
ps_5y<-ps_5y[ps_5y$LON5_D<175,]
ps_5y<-ps_5y[ps_5y$LAT5_D>-45,]
Year=2008
dat<-ps_5y[ps_5y$YY==Year,]
pts<-as.data.frame(dat[,6:5])
coords<-SpatialPoints(pts)
coords<-remove.duplicates(coords)
length(coords)

sstfiles<-list.files(paste(Year,'/cropped',sep=""),pattern=".tif",full.names = T)
# Store all sst in lists
sstlist<-lapply(sstfiles,raster)

sst<-stack(sstlist)

slope<-7.17185E-4
intercept<--2
sst<-(slope*sst)+intercept
sst[sst>45]<-NA

Year<-2008

sst_data<-data.frame("LON5_D"=numeric(),"LAT5_D"=numeric(),"YY"=numeric(),"MM"=numeric(),"sst"=double())
month<-c(1,10,11,12,2,3,4,5,6,7,8,9)
for (i in 1:length(month)){
  months=month[i]
  tmp<-as.data.frame(coords)
  l<-nrow(tmp)
  YY<-as.vector(rep(Year,l))
  MM<-as.vector(rep(months,l))
  temp<-as.vector(extract(sst[[i]],coords,fun=mean,buffer=250000,na.rm=TRUE))
  temp<-as.numeric(temp)
  tmp<-cbind(tmp,YY,MM,temp)
  sst_data<-rbind(sst_data,tmp)
  print(month[i])
}
write.csv(sst_data,paste(Year,"sst_data.csv",sep=" "))


Year<-2009
dat<-ps_5y[ps_5y$YY==Year,]
pts<-as.data.frame(dat[,6:5])
coords<-SpatialPoints(pts)
coords<-remove.duplicates(coords)
length(coords)

sstfiles<-list.files(paste(Year,'/cropped',sep=""),pattern=".tif",full.names = T)
# Store all sst in lists
sstlist<-lapply(sstfiles,raster)
names(sstlist)
sst<-stack(sstlist)



names(sst)
slope<-7.17185E-4
intercept<--2
sst<-(slope*sst)+intercept
sst[sst>45]<-NA
#plot(sst)
sst_data<-data.frame("LON5_D"=numeric(),"LAT5_D"=numeric(),"YY"=numeric(),"MM"=numeric(),"sst"=double())

for (i in 1:length(month)){
  months=month[i]
  tmp<-as.data.frame(coords)
  l<-nrow(tmp)
  YY<-as.vector(rep(Year,l))
  MM<-as.vector(rep(months,l))
  temp<-as.vector(extract(sst[[i]],coords,fun=mean,buffer=250000,na.rm=TRUE))
  temp<-as.numeric(temp)
  tmp<-cbind(tmp,YY,MM,temp)
  sst_data<-rbind(sst_data,tmp)
  print(month[i])
}
write.csv(sst_data,paste(Year,"sst_data.csv",sep=" "))



Year<-2010
dat<-ps_5y[ps_5y$YY==Year,]
pts<-as.data.frame(dat[,6:5])
coords<-SpatialPoints(pts)
coords<-remove.duplicates(coords)
length(coords)
sstfiles<-list.files(paste(Year,'/cropped',sep=""),pattern=".tif",full.names = T)
# Store all sst in lists
sstlist<-lapply(sstfiles,raster)

sst<-stack(sstlist)

slope<-7.17185E-4
intercept<--2
sst<-(slope*sst)+intercept
sst[sst>45]<-NA
#plot(sst)
sst_data<-data.frame("LON5_D"=numeric(),"LAT5_D"=numeric(),"YY"=numeric(),"MM"=numeric(),"sst"=double())

for (i in 1:length(month)){
  months=month[i]
  tmp<-as.data.frame(coords)
  l<-nrow(tmp)
  YY<-as.vector(rep(Year,l))
  MM<-as.vector(rep(months,l))
  temp<-as.vector(extract(sst[[i]],coords,fun=mean,buffer=250,na.rm=TRUE))
  temp<-as.numeric(temp)
  tmp<-cbind(tmp,YY,MM,temp)
  sst_data<-rbind(sst_data,tmp)
  print(month[i])
}
write.csv(sst_data,paste(Year,"sst_data.csv",sep=" "))

Year<-2011
dat<-ps_5y[ps_5y$YY==Year,]
pts<-as.data.frame(dat[,6:5])
coords<-SpatialPoints(pts)
coords<-remove.duplicates(coords)
length(coords)

sstfiles<-list.files(paste(Year,'/cropped',sep=""),pattern=".tif",full.names = T)
# Store all sst in lists
sstlist<-lapply(sstfiles,raster)
names(sstlist)
sst<-stack(sstlist)

names(sst)
slope<-7.17185E-4
intercept<--2
sst<-(slope*sst)+intercept
sst[sst>45]<-NA
#plot(sst)
sst_data<-data.frame("LON5_D"=numeric(),"LAT5_D"=numeric(),"YY"=numeric(),"MM"=numeric(),"sst"=double())

for (i in 1:length(month)){
  months=month[i]
  tmp<-as.data.frame(coords)
  l<-nrow(tmp)
  YY<-as.vector(rep(Year,l))
  MM<-as.vector(rep(months,l))
  temp<-as.vector(extract(sst[[i]],coords,fun=mean,buffer=250000,na.rm=TRUE))
  temp<-as.numeric(temp)
  tmp<-cbind(tmp,YY,MM,temp)
  sst_data<-rbind(sst_data,tmp)
  print(month[i])
}
write.csv(sst_data,paste(Year,"sst_data.csv",sep=" "))

Year<-2012
dat<-ps_5y[ps_5y$YY==Year,]
pts<-as.data.frame(dat[,6:5])
coords<-SpatialPoints(pts)
coords<-remove.duplicates(coords)
length(coords)

dat<-ps_5y[ps_5y$YY==Year,]
pts<-as.data.frame(dat[,6:5])
coords<-SpatialPoints(pts)
coords<-remove.duplicates(coords)
length(coords)
sstfiles<-list.files(paste(Year,'/cropped',sep=""),pattern=".tif",full.names = T)
# Store all sst in lists
sstlist<-lapply(sstfiles,raster)
names(sstlist)
sst<-stack(sstlist)

names(sst)
slope<-7.17185E-4
intercept<--2
sst<-(slope*sst)+intercept
sst[sst>45]<-NA
#plot(sst)
sst_data<-data.frame("LON5_D"=numeric(),"LAT5_D"=numeric(),"YY"=numeric(),"MM"=numeric(),"sst"=double())

for (i in 1:length(month)){
  months=month[i]
  tmp<-as.data.frame(coords)
  l<-nrow(tmp)
  YY<-as.vector(rep(Year,l))
  MM<-as.vector(rep(months,l))
  temp<-as.vector(extract(sst[[i]],coords,fun=mean,buffer=250000,na.rm=TRUE))
  temp<-as.numeric(temp)
  tmp<-cbind(tmp,YY,MM,temp)
  sst_data<-rbind(sst_data,tmp)
  print(month[i])
}
write.csv(sst_data,paste(Year,"sst_data.csv",sep=" "))

Year<-2013
dat<-ps_5y[ps_5y$YY==Year,]
pts<-as.data.frame(dat[,6:5])
coords<-SpatialPoints(pts)
coords<-remove.duplicates(coords)
length(coords)
sstfiles<-list.files(paste(Year,'/cropped',sep=""),pattern=".tif",full.names = T)
# Store all sst in lists
sstlist<-lapply(sstfiles,raster)
names(sstlist)
sst<-stack(sstlist)

names(sst)
slope<-7.17185E-4
intercept<--2
sst<-(slope*sst)+intercept
sst[sst>45]<-NA
#plot(sst)
sst_data<-data.frame("LON5_D"=numeric(),"LAT5_D"=numeric(),"YY"=numeric(),"MM"=numeric(),"sst"=double())

for (i in 1:length(month)){
  months=month[i]
  tmp<-as.data.frame(coords)
  l<-nrow(tmp)
  YY<-as.vector(rep(Year,l))
  MM<-as.vector(rep(months,l))
  temp<-as.vector(extract(sst[[i]],coords,fun=mean,buffer=250000,na.rm=TRUE))
  temp<-as.numeric(temp)
  tmp<-cbind(tmp,YY,MM,temp)
  sst_data<-rbind(sst_data,tmp)
  print(month[i])
}
write.csv(sst_data,paste(Year,"sst_data.csv",sep=" "))


data1<-read.csv("2008 sst_data.csv")
data2<-read.csv("2009 sst_data.csv")
data3<-read.csv("2010 sst_data.csv")
data4<-read.csv("2011 sst_data.csv")
data5<-read.csv("2012 sst_data.csv")
data6<-read.csv("2013 sst_data.csv")

final_file<-rbind(data1,data2,data3,data4,data5,data6)
write.csv(final_file,"wcpfc_sst.csv")
unique(data1$LAT5_D)
nrow(final_file)
804*6
nrow(ps_5y)
?read.csv
#coords<-remove.duplicates(coords)


sstspatial<-SpatialPointsDataFrame(coords,sst_data)
plot(sstspatial)
proj4string(sstspatial) <- CRS("+proj=longlat +datum=WGS84")


library(tmap)
tmap_mode("plot")

data("World")
years<-c(2008:2013)
t=2009

#Creates a panel of spatial PS total catch and bycatch rate by year
sp_plot<- tm_shape(World)+
  tm_fill(col="gray91") +
  tm_borders() +
  tm_shape(sstspatial,is.master=TRUE)+
  #tm_facets("YY")+
  #tm_grid(x=(pts$LON5_D-2.5),y=(pts$LAT5_D-2.5),projection="longlat",col="black",labels.size = .5,
  #       alpha=0.08,lwd=1)+
  tm_bubbles(size='temp', col = 'MM')
,
             border.col = "black", border.alpha = .5,
             style="fixed",scale=1.2,
             breaks=c(-Inf, seq(0, 0.20, by=0.02), Inf),
             palette = c(blu="royalblue", reddish="tomato"),alpha=0.7, contrast=1,legend.size.show=TRUE,showNA=FALSE,
             title.size = paste(t,'Total Tuna PS Catch (mt)'),legend.size.is.portrait=TRUE,
             title.col = '% Bigeye',