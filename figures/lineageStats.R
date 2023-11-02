#!/usr/bin/env Rscript
library(lubridate)
library(ggplot2)
library(dplyr)
library(data.table)
library(versionsort)
library("ggsci")

md=fread("data/meta_shrink3_oct.tsv",header=T,quote = "\"",sep="\t", stringsAsFactors=F)
md$dates=as.Date(md$`Collection date`)

md2023 = md[md$dates>as.Date("2023-01-01") & md$dates<as.Date("2023-10-07"),]
md2023$days=as.Date(lubridate::floor_date(md2023$dates, "day"))
md2023$weeks=as.Date(lubridate::floor_date(md2023$dates, "week"))
md2023$months=as.Date(lubridate::floor_date(md2023$dates, "month"))

totalweeks=md2023 %>%
  group_by(weeks) %>%
  summarise(totalofweek=n())
totalweekspango=md2023 %>%
  group_by(weeks,`short_pango`) %>%
  summarise(pangoofweek=n())

totalweekspangomerge = merge(totalweeks,totalweekspango,by.x="weeks",by.y="weeks")

cutoff=0.01
# we take the pango lineages that are at least once (one month) > cutoff 
cutofflineages=unique(totalweekspangomerge[totalweekspangomerge$pangoofweek/totalweekspangomerge$totalofweek>cutoff,"short_pango"])

md2023$short_pango_clean="XXX.Others"
md2023$full_pango_clean="XXX.Others"
md2023$short_pango_clean[md2023$short_pango%in%cutofflineages]= md2023$short_pango[md2023$short_pango%in%cutofflineages]
md2023$full_pango_clean[md2023$short_pango%in%cutofflineages]= md2023$full_pango[md2023$short_pango%in%cutofflineages]

pangofactors=unique(md2023[,c("short_pango_clean","full_pango_clean")])

# We sort the pangolin alias in the way of version sort according to full pango name (unaliased)
md2023$short_pango_clean=factor(md2023$short_pango_clean, levels = pangofactors$short_pango_clean[rev(ver_order(pangofactors$full_pango_clean))]) 

palette1 <- c("#7d7d7d","yellow","#FF8000",colorRampPalette(c("lightblue", "darkblue"))(13),"#009400","#006400","pink",colorRampPalette(c("red", "darkred"))(4))

# We build the dataframe with lineage labels
# associated to their week with the most abundance (x position)
# and the cumulated abundance (y position of the label)
totalweekspangomerge$ratio = totalweekspangomerge$pangoofweek/totalweekspangomerge$totalofweek
ratiogrouped=totalweekspangomerge%>%group_by(short_pango)%>%filter(ratio == max(ratio))
ratiogroupedok=ratiogrouped[ratiogrouped$short_pango%in%unique(md2023$short_pango_clean),]
lpos=data.frame(matrix(nrow = 0, ncol = 3)) 
for (l1 in rev(levels(md2023$short_pango_clean))){
  short_pango_clean=l1
  cumu=0
  week=ratiogroupedok$weeks[ratiogroupedok$short_pango==l1]
  if(length(week)>0){
    week=week[1]
  }else{
    week=min(ratiogroupedok$weeks)
  }
  for (l2 in rev(levels(md2023$short_pango_clean))){
    r=totalweekspangomerge$ratio[totalweekspangomerge$weeks==week & totalweekspangomerge$short_pango==l2]
    if(length(r)>0){
      r=r[1]
    }else{
      r=0
    }
    if(l1==l2){
      cumu=cumu+r/2
      lpos=rbind(lpos,data.frame(short_pango_clean,week,cumu))
      break
    }else{
      cumu=cumu+r
    }
  }
}

svg(file="ba.2.86_oct.world_dates.svg",width=20,height=10)
ggplot(md2023, aes(x = dates,fill=short_pango_clean,group=short_pango_clean,after_stat(count)))+
  geom_density(adjust=2, position="fill", kernel="cosine",linewidth=0) +
  #geom_text(aes(label=short_pango_clean),position = position_fill(vjust = 0.5))+
  ggtitle("Prop of SARS-CoV-2 lineages per week on GISAID") + 
  xlab("date")+
  scale_fill_manual(values=palette1)+
  theme_bw()+
  theme(text = element_text(size = 15,family = "sans"))+
  geom_text(data=lpos,aes(label=short_pango_clean,x=week,y=cumu),color="black",vjust=0.5,hjust=0.5)
dev.off()


###############################################
## Now the same with Only European sequences ##
###############################################

md2023 = md[md$dates>as.Date("2023-01-01") & md$dates<as.Date("2023-10-07") & grepl("Europe",md$Location),]
md2023$days=as.Date(lubridate::floor_date(md2023$dates, "day"))
md2023$weeks=as.Date(lubridate::floor_date(md2023$dates, "week"))
md2023$months=as.Date(lubridate::floor_date(md2023$dates, "month"))

totalweeks=md2023 %>%
  group_by(weeks) %>%
  summarise(totalofweek=n())
totalweekspango=md2023 %>%
  group_by(weeks,`short_pango`) %>%
  summarise(pangoofweek=n())

totalweekspangomerge = merge(totalweeks,totalweekspango,by.x="weeks",by.y="weeks")

cutoff=0.009
# we take the pango lineages that are at least once > cutoff
cutofflineages=unique(totalweekspangomerge[totalweekspangomerge$pangoofweek/totalweekspangomerge$totalofweek>cutoff,"short_pango"])

md2023$short_pango_clean="XXX.Others"
md2023$full_pango_clean="XXX.Others"
md2023$short_pango_clean[md2023$short_pango%in%cutofflineages]= md2023$short_pango[md2023$short_pango%in%cutofflineages]
md2023$full_pango_clean[md2023$short_pango%in%cutofflineages]= md2023$full_pango[md2023$short_pango%in%cutofflineages]

pangofactors=unique(md2023[,c("short_pango_clean","full_pango_clean")])

# We sort the pangolin alias in the way of version sort according to full pango name (unaliased)
md2023$short_pango_clean=factor(md2023$short_pango_clean, levels = pangofactors$short_pango_clean[rev(ver_order(pangofactors$full_pango_clean))]) 

palette1 <- c("#7d7d7d","yellow","#FF8000",colorRampPalette(c("lightblue", "darkblue"))(13),"#009400","#006400","pink",colorRampPalette(c("red", "darkred"))(4))

# We build the dataframe with lineage labels
# associated to their week with the most abundance (x position)
# and the cumulated abundance (y position of the label)
totalweekspangomerge$ratio = totalweekspangomerge$pangoofweek/totalweekspangomerge$totalofweek
ratiogrouped=totalweekspangomerge%>%group_by(short_pango)%>%filter(ratio == max(ratio))
ratiogroupedok=ratiogrouped[ratiogrouped$short_pango%in%unique(md2023$short_pango_clean),]
lpos=data.frame(matrix(nrow = 0, ncol = 3)) 
for (l1 in rev(levels(md2023$short_pango_clean))){
  short_pango_clean=l1
  cumu=0
  week=ratiogroupedok$weeks[ratiogroupedok$short_pango==l1]
  if(length(week)>0){
    week=week[1]
  }else{
    week=min(ratiogroupedok$weeks)
  }
  for (l2 in rev(levels(md2023$short_pango_clean))){
    r=totalweekspangomerge$ratio[totalweekspangomerge$weeks==week & totalweekspangomerge$short_pango==l2]
    if(length(r)>0){
      r=r[1]
    }else{
      r=0
    }
    if(l1==l2){
      cumu=cumu+r/2
      lpos=rbind(lpos,data.frame(short_pango_clean,week,cumu))
      break
    }else{
      cumu=cumu+r
    }
  }
}

svg(file="ba.2.86_oct.europe_dates.svg",width=20,height=10)
ggplot(md2023, aes(x = dates,fill=short_pango_clean,group=short_pango_clean,after_stat(count)))+
  geom_density(adjust=2, position="fill", kernel="cosine",linewidth=0) +
  #geom_text(aes(label=short_pango_clean),position = position_fill(vjust = 0.5))+
  ggtitle("Prop of SARS-CoV-2 lineages in Europe per week on GISAID") + 
  xlab("date")+
  scale_fill_manual(values=palette1)+
  theme_bw()+
  theme(text = element_text(size = 15,family = "sans"))+
  geom_text(data=lpos,aes(label=short_pango_clean,x=week,y=cumu),color="black",vjust=0.5,hjust=0.5)
dev.off()


#### Statistics for French sequences from Sept 2023
md=fread("data/meta_shrink3_sept_fr.tsv",header=T,quote = "\"",sep="\t", stringsAsFactors=F)
md$dates=as.Date(md$`Collection date`)
mdFranceSept = md[md$dates>=as.Date("2023-09-01") & md$dates<as.Date("2023-10-01") & grepl("France",md$Location),]
