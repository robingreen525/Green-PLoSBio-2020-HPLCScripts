rm(list = ls())
#setup
{
  #go to data directory
  require(ggplot2)
  #require(AUC)
  require(flux)
  require(stringr)
  hplcData<-NULL
  Spectra_Short<-data.frame()
  
  mytheme_biggerfonts=list( theme_bw(),theme(axis.text.x = element_text(face="bold",size=20),axis.text.y = element_text(face="bold",size=20),
                                             strip.text.x=element_text(face="bold",size=15),strip.text.y=element_text(face="bold",size=15)),
                            theme(text=element_text(size=20,face="bold"),legend.position='none'))
  
  mytheme_biggerfonts_withLegend=list( theme_bw(),theme(axis.text.x = element_text(face="bold",size=20),axis.text.y = element_text(face="bold",size=20),
                                                        strip.text.x=element_text(face="bold",size=15),strip.text.y=element_text(face="bold",size=15)),
                                       theme(text=element_text(size=20,face="bold"),legend.position='right'))
  
  blanktheme<- theme(axis.line=element_blank(),axis.text.x=element_blank(),
                     axis.text.y=element_blank(),axis.ticks=element_blank(),
                     axis.title.x=element_blank(),
                     axis.title.y=element_blank(),legend.position="none",
                     panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank(),plot.background=element_blank())#SD+Lys, 
  
}
source('HPLC_Functions.R')

files<-list.files(path='20171030_exportdataq/supernatants/')

for(file in files)
{
  print(file)
  timerep<-as.numeric(substr(strsplit(file,split = '_')[[1]][3],start=1,stop=1))
  sampleclass<-as.character(((strsplit(file,split = '_')[[1]][1])))



    nam<-strsplit(file,split='.txt')[[1]][1]
 
   {
 
    start=6.25
    end=8


    vari='did not work'
    tryCatch({
      #read in data

      r<-read.table(paste('20171030_exportdataq/supernatants/',file,sep='/'),header=T,sep='\t',skip=26)
      colnames(r)<-c('RT','Intensity')
      excitation<-rep(400,times=nrow(r))
      emission<-rep(465,times=(nrow(r)))
      r<-cbind(r,excitation,emission)

      sample<-rep(nam,times=nrow(r))
      sampleclass<-rep(sampleclass,times=nrow(r))
      timerep<-rep(timerep,times=nrow(r))
      minutesIncubation<-rep(20,times=nrow(r))
      r<-cbind(r,sample,sampleclass,minutesIncubation,timerep)
      g<-ggplot(data=r,aes(x=RT,y=Intensity,col=sample))+geom_point()
      print(g)
      r<-AreaIntegrate_MoveDownFromMaxForLocalMin(r,start=start,end=end)
      Spectra_Short<-rbind(Spectra_Short,r)
    },error = function(e) print(vari))
  }



}


if(!(dir.exists('Figures/'))){dir.create('Figures/')}
#
 g<-ggplot(data=Spectra_Short,aes(x=RT,y=CorrectedIntense,col=sampleclass))+
   geom_point()+facet_grid(sampleclass ~ timerep)
print(g)
 ggsave(g,filename = 'Figures/HPLC_Supernatants.png')

mean_areas<-data.frame()

for(samp in unique(Spectra_Short$sample))
{
  print(samp)
  t<-Spectra_Short[which(Spectra_Short$sample==samp),]
  t$Intensity<-t$CorrectedIntense
  #Spectra_baseline_subtract<-rbind(Spectra_baseline_subtract,t)
  a<-auc.mc(x=t$RT,y=t$Intensity)
  avg<-mean(a)
  sd<-sd(a)


  df<-data.frame(samp,avg,sd,t$timerep[1],t$sampleclass[1])
  colnames(df)<-c('Sample','Average_Area','SD_Area','timerep','sampleclass')
  mean_areas<-rbind(mean_areas,df)
}

write.csv(mean_areas,file = 'supernatants.csv')

# # # plot_dat<-mean_areas[which(mean_areas$Sample !='0.03uMGSHinSD_timerep3_20min_3' ),]
# # # plot_dat<-plot_dat[which(plot_dat$uMGSH>=0.03),]
# # # #plot_dat<-plot_dat[which(plot_dat$timerep>1),]
# # # g<-ggplot(data=plot_dat,aes(x=(as.numeric(as.character(plot_dat$uMGSH))),y=(Average_Area)))+geom_jitter(width=0.1,aes(shape=as.factor(timerep),col=as.factor(timerep)))+
# # #   geom_smooth(method='lm',formula=y~x)+scale_y_log10()+scale_x_log10()+annotation_logticks(sides='lb')
# # # print(g)
# # # # dir.create('./Figures')
# # # ggsave(g,filename = './Figures/StdCurve.pdf')
# # # 
# # # 
# # # g<-ggplot(data=plot_dat,aes(x=(as.numeric(as.character(plot_dat$uMGSH))),y=(Average_Area)))+geom_jitter(size=4,width=0.1,aes(shape=as.factor(timerep),col=as.factor(timerep)))+
# # #   geom_smooth(method='lm',formula=y~x)+scale_y_log10()+scale_x_log10()+annotation_logticks(sides='lb')+
# # #   theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
# # # print(g)
# # # ggsave(g,filename = './Figures/StdCurve_all.pdf')