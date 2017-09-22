

#----------------------------------
# cambodia-heatmap.R
#
# heatmap of individual level
# antibody response
#----------------------------------

#----------------------------------
# input files:
#   cambodia_serology_public.rds
#
# output files:
#   cambodia-Ab-heatmap.pdf
#----------------------------------

#----------------------------------
# preamble
#----------------------------------

rm(list=ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)

#----------------------------------
# load the dataset
#----------------------------------
d <- readRDS("~/dropbox/cambodia/data/final/cambodia_serology_public.rds")

#----------------------------------
# subset data to antibodies
#----------------------------------
mbavars <- c("ttmb","sag2a","nie","t24","wb123","bm14","bm33","pvmsp19","pfmsp19")
d <- subset(d,select=c("region","psuid",mbavars)) 

#----------------------------------
# count and recode values <=0 to 1
#----------------------------------
apply(d[mbavars],2,function(x) table(x<=0) )
for(vv in mbavars) {
  d[vv][d[vv]<=0]<-1
}
#----------------------------------
# convert values to log10 scale
#----------------------------------
for(vv in mbavars) {
  d[vv] <- log10(d[vv])
}

#----------------------------------
# calculate mean antibody response
# (used for sorting indivs in heatmap)
#----------------------------------
d$abmean <- rowMeans(d[mbavars])
d <- group_by(d,region) %>%
  arrange(region,abmean) %>%
  mutate(indivrank = 1:n())

#----------------------------------
# gather data into long format
# for plotting, and sort by 
#----------------------------------
dlong <- gather(d,antigen,mfi,-region,-psuid,-abmean,-indivrank)  %>%
  arrange(region,indivrank)

#----------------------------------
# format factors for plotting
#----------------------------------
levels(dlong$region) <- c("Phnom Penh","Southeast","Southwest","West","North")


Abnames <- c("Tetanus toxiod",
             "Toxoplasma gondii SAG2A",
             "Strongyloides stercoralis NIE",
             "Taenia solium T24H",
             "Lymphatic filariasis Wb123",
             "Lymphatic filariasis Bm14",
             "Lymphatic filariasis Bm33",
             "Plasmodium vivax MSP-1(19)",
             "Plasmodium falciparum MSP1(19)")

dlong$ab <- factor(dlong$antigen,levels=rev(mbavars),labels=rev(Abnames))

dlong$abgroup <- factor(NA,levels=c("VPDs","NTDs","Malaria"))
dlong$abgroup[dlong$antigen %in% c("ttmb")] <-"VPDs"
dlong$abgroup[dlong$antigen %in% c("sag2a","nie","t24","wb123","bm14","bm33")] <-"NTDs"
dlong$abgroup[dlong$antigen %in% c("pvmsp19","pfmsp19")] <-"Malaria"

#----------------------------------
# make a heatmap of all antigens
# inspiration: http://www.roymfrancis.com/a-guide-to-elegant-tiled-heatmaps-in-r/
#----------------------------------

#define a colour for fonts
textcol <- "grey20"

p <- ggplot(dlong,aes(x=indivrank,y=ab,fill=mfi)) +
  #facet over village
  facet_grid(abgroup~region,scales='free',space='free',switch="y")+
  
  #add border white colour of line thickness 0.25
  # geom_tile(colour="white",size=0.25)+
  geom_tile() +
  #remove y axis labels, 
  labs(x="Individuals are sorted by mean antibody response within region",y="",title="")+
  #remove extra space
  scale_y_discrete(expand=c(0,0),position="right")+
  scale_x_continuous(expand=c(0,0))+
  # scale_x_discrete(expand=c(0,0),
  #                  breaks=1:9,labels=1:9)+
  #change the scale_fill_manual
  scale_fill_distiller(palette="BuGn",na.value="grey90",
                       direction=0,
                       guide=guide_colorbar(title="log10\nMFI-bg",face='bold'))+
  # scale_fill_manual(values=rev(brewer.pal(7,"YlGnBu")),na.value="grey90",guide=guide_colorbar(title="log10\nMFI-bg",face='bold'))+
  #one unit on x-axis is equal to one unit on y-axis.
  #equal aspect ratio x and y axis
  coord_equal()+
  #set base size for all font elements
  theme_grey(base_size=10)+
  #theme options
  theme(
    
    legend.title=element_text(color=textcol,size=8),
    legend.position = "left",
    #reduce/remove legend margin
    legend.margin = margin(grid::unit(0.1,"cm")),
    #change legend text properties
    legend.text=element_text(colour=textcol,size=7,face="bold"),
    #change legend key height
    legend.key.height=grid::unit(0.8,"cm"),
    #set a slim legend
    legend.key.width=grid::unit(0.4,"cm"),
    
    #set remove x axis ticks
    axis.title.x=element_text(hjust=0),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    #set y axis text colour and adjust vertical justification
    axis.text.y=element_text(size=10,vjust = 0.2,colour=textcol),
    #change axis ticks thickness
    axis.ticks=element_line(size=0.4),
    
    #change title font, size, colour and justification
    plot.title=element_text(colour=textcol,hjust=0,size=12,face="bold"),
    #adjust facet labels
    strip.text.x = element_text(size=10),
    strip.text.y = element_text(size=8),
    #remove plot background
    plot.background=element_blank(),
    #remove plot border
    panel.border=element_blank(),
    #remove plot margins
    # plot.margin=margin(grid::unit(1,"cm"))
  )

#p

ggsave(filename="~/dropbox/cambodia/results/figs/cambodia-Ab-heatmap.pdf",plot = p,device='pdf',width=13,height=4)


