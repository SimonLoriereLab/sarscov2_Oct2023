#!/usr/bin/env Rscript
library(lubridate)
library(ggplot2)
library(dplyr)
library(data.table)
library(versionsort)
library("ggsci")


setwd("/pasteur/zeus/projets/p01/CNRVIR_bioinfoSEQ/users/flemoine/2023_09_BA.2.86/")
md=fread("data/meta_shrink3.tsv",header=T,quote = "\"",sep="\t", stringsAsFactors=F)
md$dates=as.Date(md$`Collection date`)
md2023 = md[md$dates>as.Date("2023-01-01") & md$dates<as.Date("2023-09-10"),]
md2023$days=as.Date(lubridate::floor_date(md2023$dates, "day"))
md2023$weeks=as.Date(lubridate::floor_date(md2023$dates, "week"))
md2023$months=as.Date(lubridate::floor_date(md2023$dates, "month"))


totalweeks=md2023 %>%
  group_by(weeks) %>%
  summarise(totalofweek=n())
#totalweekspango=md2023 %>%
#  group_by(weeks,`Pango lineage`) %>%
#  summarise(pangoofweek=n())
totalweekspango=md2023 %>%
  group_by(weeks,`short_pango`) %>%
  summarise(pangoofweek=n())

totalweekspangomerge = merge(totalweeks,totalweekspango,by.x="weeks",by.y="weeks")

cutoff=0.02
# we take the pango lineages that are at least once > cutoff
cutofflineages=unique(totalweekspangomerge[totalweekspangomerge$pangoofweek/totalweekspangomerge$totalofweek>cutoff,"short_pango"])

#md2023=merge(md2023,totalweeks, by.x="weeks", by.y="weeks")
#md2023=merge(md2023,totalweekspango, by.x=c("weeks","Pango lineage"), by.y=c("weeks","Pango lineage"))
#md2023=merge(md2023,totalweekspango, by.x=c("weeks","short_pango"), by.y=c("weeks","short_pango"))

#md2023[md2023$pangoofweek/md2023$totalofweek<cutoff & md2023$`Pango lineage`!="BA.2.86","Pango lineage"]="Others"
md2023$short_pango_clean="XXX.Others"
md2023$full_pango_clean="XXX.Others"
md2023$short_pango_clean[md2023$short_pango%in%cutofflineages]= md2023$short_pango[md2023$short_pango%in%cutofflineages]
md2023$full_pango_clean[md2023$short_pango%in%cutofflineages]= md2023$full_pango[md2023$short_pango%in%cutofflineages]

pangofactors=unique(md2023[,c("short_pango_clean","full_pango_clean")])

# On trie l'alias façon versions selon le full pango (unaliased)
md2023$short_pango_clean=factor(md2023$short_pango_clean, levels = pangofactors$short_pango_clean[rev(ver_order(pangofactors$full_pango_clean))]) 

#mdfrance=md2023[grep('France', md2023$`Virus name`,'France'),]


palette1 <- c("#7d7d7d","yellow","#FF8000",colorRampPalette(c("lightblue", "darkblue"))(13),"#009400","#006400","pink",colorRampPalette(c("red", "darkred"))(4))

# We build the dataframe with lineage labels
# associated to their week with the most abundance (x position)
# and the cumulated abundance (y position of the label)
totalweekspangomerge$ratio = totalweekspangomerge$pangoofweek/totalweekspangomerge$totalofweek
ratiogrouped=totalweekspangomerge%>%group_by(short_pango)%>%filter(ratio == max(ratio))
ratiogroupedok=ratiogrouped[ratiogrouped$short_pango%in%unique(md2023$short_pango_clean),]
#ratiogrouped[grep("XBB.1.18",ratiogrouped$short_pango),]
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

svg(file="ba.2.86.world_dates.svg",width=20,height=10)
#pdf(file="pango_europe.pdf",width=18,height=10)
#ggplot(md2023, aes(x = weeks,fill=`Pango lineage`,group=`Pango lineage`,after_stat(count)))+
#  geom_density(adjust=2, position="fill", kernel="cosine") +
#  ggtitle("Prop of SARS-CoV-2 sequences per month on Emergen, by Clade") + 
#  xlab("date")
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



# Uniquement les séquences Européennes

md2023 = md[md$dates>as.Date("2023-01-01") & md$dates<as.Date("2023-09-10") & grepl("Europe",md$Location),]
md2023$days=as.Date(lubridate::floor_date(md2023$dates, "day"))
md2023$weeks=as.Date(lubridate::floor_date(md2023$dates, "week"))
md2023$months=as.Date(lubridate::floor_date(md2023$dates, "month"))


totalweeks=md2023 %>%
  group_by(weeks) %>%
  summarise(totalofweek=n())
#totalweekspango=md2023 %>%
#  group_by(weeks,`Pango lineage`) %>%
#  summarise(pangoofweek=n())
totalweekspango=md2023 %>%
  group_by(weeks,`short_pango`) %>%
  summarise(pangoofweek=n())

totalweekspangomerge = merge(totalweeks,totalweekspango,by.x="weeks",by.y="weeks")

cutoff=0.009
# we take the pango lineages that are at least once > cutoff
cutofflineages=unique(totalweekspangomerge[totalweekspangomerge$pangoofweek/totalweekspangomerge$totalofweek>cutoff,"short_pango"])

#md2023=merge(md2023,totalweeks, by.x="weeks", by.y="weeks")
#md2023=merge(md2023,totalweekspango, by.x=c("weeks","Pango lineage"), by.y=c("weeks","Pango lineage"))
#md2023=merge(md2023,totalweekspango, by.x=c("weeks","short_pango"), by.y=c("weeks","short_pango"))

#md2023[md2023$pangoofweek/md2023$totalofweek<cutoff & md2023$`Pango lineage`!="BA.2.86","Pango lineage"]="Others"
md2023$short_pango_clean="XXX.Others"
md2023$full_pango_clean="XXX.Others"
md2023$short_pango_clean[md2023$short_pango%in%cutofflineages]= md2023$short_pango[md2023$short_pango%in%cutofflineages]
md2023$full_pango_clean[md2023$short_pango%in%cutofflineages]= md2023$full_pango[md2023$short_pango%in%cutofflineages]

pangofactors=unique(md2023[,c("short_pango_clean","full_pango_clean")])

# On trie l'alias façon versions selon le full pango (unaliased)
md2023$short_pango_clean=factor(md2023$short_pango_clean, levels = pangofactors$short_pango_clean[rev(ver_order(pangofactors$full_pango_clean))]) 

#mdfrance=md2023[grep('France', md2023$`Virus name`,'France'),]


palette1 <- c("#7d7d7d","yellow","#FF8000",colorRampPalette(c("lightblue", "darkblue"))(13),"#009400","#006400","pink",colorRampPalette(c("red", "darkred"))(4))

# We build the dataframe with lineage labels
# associated to their week with the most abundance (x position)
# and the cumulated abundance (y position of the label)
totalweekspangomerge$ratio = totalweekspangomerge$pangoofweek/totalweekspangomerge$totalofweek
ratiogrouped=totalweekspangomerge%>%group_by(short_pango)%>%filter(ratio == max(ratio))
ratiogroupedok=ratiogrouped[ratiogrouped$short_pango%in%unique(md2023$short_pango_clean),]
#ratiogrouped[grep("XBB.1.18",ratiogrouped$short_pango),]
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

svg(file="ba.2.86.europe_dates.svg",width=20,height=10)
#pdf(file="pango_europe.pdf",width=18,height=10)
#ggplot(md2023, aes(x = weeks,fill=`Pango lineage`,group=`Pango lineage`,after_stat(count)))+
#  geom_density(adjust=2, position="fill", kernel="cosine") +
#  ggtitle("Prop of SARS-CoV-2 sequences per month on Emergen, by Clade") + 
#  xlab("date")
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


