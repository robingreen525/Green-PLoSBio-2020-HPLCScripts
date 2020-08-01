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
#source('HPLC_Functions.R')
require(ggplot2)
r<-read.table('20170915_exportdata/0.1uMGSHinSD_BioRep3_20min_3.txt',header=T,sep='\t',skip=26)

colnames(r)<-c('RT','Intensity')
excitation<-rep(400,times=nrow(r))
emission<-rep(465,times=(nrow(r)))
r<-cbind(r,excitation,emission)

sample<-rep('0.1uMGSHinSD_BioRep3_20min_3',times=nrow(r))
uM.GSH<-rep(0.3,times=nrow(r))
timerep<-rep(1,times=nrow(r))
biorep<-rep(1,time=nrow(r))
minutesIncubation<-rep(20,times=nrow(r))
r<-cbind(r,sample,uM.GSH,minutesIncubation,timerep,biorep)


AreaIntegrate<-function(sample,start=3.3,end=3.55)
{
  print(head(sample))
  sample<-sample[which(sample$RT > start),]
  sample<-sample[which(sample$RT < end),]
  
  #save images for future reference/troubleshooting.
  
  if(!(dir.exists('HPLC_Traces/')))
  {
    dir.create('HPLC_Traces/')
  }
  
  nam<-as.character(r$sample[1]) # get name for folder, all images saved in this folder
  
  g<-ggplot(data=sample,aes(x=RT,y=Intensity,col=sample))+geom_point()
  print(g)
  
  if(!(dir.exists(paste('HPLC_Traces',nam,sep='/'))))
  {
    dir.create((paste('HPLC_Traces',nam,sep='/')))
    
  }
  ggsave(g,filename = paste((paste('HPLC_Traces',nam,sep='/')),'1_initial_trace.png',sep='/'))
  
  
  
  #define max point
  max_Intense<-max(sample$Intensity)
  peak_RT<-sample[which(sample$Intensity == max_Intense),]$RT[1]
  print(peak_RT)

  g<-ggplot(data=sample,aes(x=RT,y=Intensity,col=sample))+geom_point()+
    geom_vline(xintercept=peak_RT)
  print(g)
  ggsave(g,filename = paste((paste('HPLC_Traces',nam,sep='/')),'2_peak_max_find.png',sep='/'))
  #first half min
  t<-sample[which(sample$RT<peak_RT) ,]
  print(head(t))

  g<-ggplot(data=t,aes(x=RT,y=Intensity,col=sample))+geom_point()+
    geom_vline(xintercept=peak_RT)
  print(g)

  
  # algo for finding local min
  {
  y_vals<-t$Intensity
  length_half_trace<-length(y_vals)
  working_range<-c(1:length_half_trace)
  
  #flip because moving right to level
  
  #working_range<-rev(working_range)
 
  #get m=n/2
  n=max(working_range)
  m=floor(n/2)
  
  A_m<-y_vals[m]
  A_m_minus_1<-y_vals[m-1]
  A_m_plus_1<-y_vals[m+1]
  if(A_m_minus_1 < A_m)
  {
   while(A_m_minus_1 < A_m)
   {
     m=m-1
     A_m<-y_vals[m]
     A_m_minus_1<-y_vals[m-1]
     
   }
    
    min_first_half_intense<-y_vals[m]
    min_first_half_RT<-t[which(t$Intensity == min_first_half_intense),]$RT
    
    
  }
  
  else if(A_m_minus_1 > A_m )
  {
    while(A_m_minus_1 > A_m && A_m >=A_m_plus_1)
    {
      m=m+1
      A_m<-y_vals[m]
      A_m_minus_1<-y_vals[m-1]
      
    }
    min_first_half_intense<-y_vals[m]
    min_first_half_RT<-t[which(t$Intensity == min_first_half_intense),]$RT
    
  }
  else if(A_m_minus_1 > A_m && A_m<A_m_plus_1 )
  {
    min_first_half_intense<-y_vals[m]
    min_first_half_RT<-t[which(t$Intensity == min_first_half_intense),]$RT
    
  }
  
  
  }
  
  
  
  # min_first_half_intense<-min(t$Intensity)
  # min_first_half_RT<-t[which(t$Intensity == min_first_half_intense),]$RT
  # 
  g<-ggplot(data=t,aes(x=RT,y=Intensity,col=sample))+geom_point()+
    geom_vline(xintercept=peak_RT)+
    geom_point(aes(x=min_first_half_RT,y=min_first_half_intense),shape=25,size=5)
  print(g)
  ggsave(g,filename = paste((paste('HPLC_Traces',nam,sep='/')),'3_firsthalfpeak_min.png',sep='/'))

  #
  #second half min
  t<-sample[which(sample$RT>peak_RT) ,]
  print(head(t))

  g<-ggplot(data=t,aes(x=RT,y=Intensity,col=sample))+geom_point()+
    geom_vline(xintercept=peak_RT)
  print(g)
  
  # algo for finding local min
  {
    y_vals<-t$Intensity
    length_half_trace<-length(y_vals)
    working_range<-c(1:length_half_trace)
    
    #flip because moving right to level
    
    #working_range<-rev(working_range)
    
    #get m=n/2
    n=max(working_range)
    m=floor(n/2)
    
    A_m<-y_vals[m]
    A_m_minus_1<-y_vals[m-1]
    A_m_plus_1<-y_vals[m+1]
    if(A_m_minus_1 < A_m)
    {
      while(A_m_minus_1 < A_m)
      {
        m=m-1
        A_m<-y_vals[m]
        A_m_minus_1<-y_vals[m-1]
        
      }
      
      min_second_half_intense<-y_vals[m]
      min_second_half_RT<-t[which(t$Intensity == min_second_half_intense),]$RT
      
      
    }
    
    else if(A_m_minus_1 > A_m&& A_m >=A_m_plus_1 )
    {
      while(A_m_minus_1 > A_m)
      {
        m=m+1
        A_m<-y_vals[m]
        A_m_minus_1<-y_vals[m-1]
        
      }
      min_second_half_intense<-y_vals[m]
      min_second_half_RT<-t[which(t$Intensity == min_second_half_intense),]$RT
      
    }
    else if(A_m_minus_1 > A_m && A_m<A_m_plus_1 )
    {
      min_second_half_intense<-y_vals[m]
      min_second_half_RT<-t[which(t$Intensity == min_second_half_intense),]$RT
      
    }
    
    
  }

  #
  g<-ggplot(data=t,aes(x=RT,y=Intensity,col=sample))+geom_point()+
    geom_vline(xintercept=peak_RT)+
    geom_point(aes(x=min_second_half_RT,y= min_second_half_intense),shape=25,size=5)
  print(g)
  ggsave(g,filename = paste((paste('HPLC_Traces',nam,sep='/')),'4_secondhalfpeak_min.png',sep='/'))

  #
  #see whole plot
  g<-ggplot(data=sample,aes(x=RT,y=Intensity,col=sample))+geom_point()+
    geom_vline(xintercept=peak_RT)+
    geom_point(aes(x=min_second_half_RT,y= min_second_half_intense),shape=25,size=5)+
    geom_point(aes(x=min_first_half_RT,y=min_first_half_intense),shape=25,size=5)+
    geom_segment(x=min_first_half_RT,y=min_first_half_intense,xend=min_second_half_RT,yend=min_second_half_intense,col='black')

  print(g)
  ggsave(g,filename = paste((paste('HPLC_Traces',nam,sep='/')),'5_WholeCrossSection.png',sep='/'))

  #
  slope=((min_second_half_intense-min_first_half_intense)/(min_second_half_RT-min_first_half_RT))
  print(slope)

  # #now create the baseline y=slope*x-b
  #
  sample<-sample[which(sample$RT < min_second_half_RT),]
  sample<-sample[which(sample$RT > min_first_half_RT),]
  g<-ggplot(data=sample,aes(x=RT,y=Intensity,col=sample))+geom_point()+
    geom_vline(xintercept=peak_RT)+
    geom_point(aes(x=min_second_half_RT,y= min_second_half_intense),shape=25,size=5)+
    geom_point(aes(x=min_first_half_RT,y=min_first_half_intense),shape=25,size=5)+
    geom_segment(x=min_first_half_RT,y=min_first_half_intense,xend=min_second_half_RT,yend=min_second_half_intense,col='black')

  print(g)
  ggsave(g,filename = paste((paste('HPLC_Traces',nam,sep='/')),'6_WholeCrossSection_withSlopedBaseline.png',sep='/'))


  baseline_y=data.frame()

  for(i in 1:nrow(sample))
  {
    rt=sample[i,]$RT
    first_rt<-sample[1,]$RT
    baseline=slope*(rt-first_rt)+(min_first_half_intense)
    df<-data.frame(baseline)
    colnames(df)<-c('baseline_y')
    baseline_y<-rbind(baseline_y,df)
  }

  sample<-cbind(sample,baseline_y)

  g<-ggplot(data=sample,aes(x=RT,y=Intensity,col=sample))+geom_point()+
    geom_vline(xintercept=peak_RT)+
    geom_point(aes(x=min_second_half_RT,y= min_second_half_intense),shape=25,size=5)+
    geom_point(aes(x=min_first_half_RT,y=min_first_half_intense),shape=25,size=5)+
    geom_segment(x=min_first_half_RT,y=min_first_half_intense,xend=min_second_half_RT,yend=min_second_half_intense,col='black')+
    geom_point(data=sample,aes(x=RT,y=baseline_y,col='blue'))

  print(g)
  ggsave(g,filename = paste((paste('HPLC_Traces',nam,sep='/')),'6_WholeCrossSection_withSlopedBaseline.png',sep='/'))

  #
  # #now correct the peak with baseline
  #
  sample$CorrectedIntense<-sample$Intensity-sample$baseline_y

  g<-ggplot(data=sample,aes(x=RT,y=CorrectedIntense,col=sample))+geom_point()
  print(g)
  ggsave(g,filename = paste((paste('HPLC_Traces',nam,sep='/')),'7_final_corrected.png',sep='/'))
  if(!(dir.exists('./HPLC_Traces/corrected_all/')))
  {
    dir.create('./HPLC_Traces/corrected_all/')
    
  }
  
  
  ggsave(g,filename = paste((paste('HPLC_Traces/corrected_all',nam,sep='/')),'png',sep='.'))
  
  return(sample)
}

AreaIntegrate(r,start=6.75,end=7.4)