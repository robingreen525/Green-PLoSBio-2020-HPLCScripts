#functions used
{

  AreaIntegrate_MoveDownFromMaxForLocalMin<-function(sample,start,end)
  {
    print(head(sample))
    sample<-sample[which(sample$RT > start),]
    sample<-sample[which(sample$RT < end),]
    
    #save images for future reference/troubleshooting.
    print('1')
    nam<-(sample$file[1])
    if(!(dir.exists('HPLC_Traces/')))
    {
      dir.create('HPLC_Traces/')
    }
    
        g<-ggplot(data=sample,aes(x=RT,y=Intensity,col=file))+geom_point()
    print(g)
    
    if(!(dir.exists(paste('HPLC_Traces',nam,sep='/'))))
    {
      dir.create((paste('HPLC_Traces',nam,sep='/')))
      
    }
    ggsave(g,filename = paste((paste('HPLC_Traces',nam,sep='/')),'1_initial_trace.png',sep='/'))
    
    print('2')
    
    
    #define max point
    max_Intense<-max(sample$Intensity)
    peak_RT<-sample[which(sample$Intensity == max_Intense),]$RT[1]
    print(peak_RT)
    
    g<-ggplot(data=sample,aes(x=RT,y=Intensity,col=file))+geom_point()+
      geom_vline(xintercept=peak_RT)
    print(g)
    ggsave(g,filename = paste((paste('HPLC_Traces',nam,sep='/')),'2_peak_max_find.png',sep='/'))
    #first half min
    t<-sample[which(sample$RT<peak_RT) ,]
    print(head(t))
    
    g<-ggplot(data=t,aes(x=RT,y=Intensity,col=file))+geom_point()+
      geom_vline(xintercept=peak_RT)
    print(g)
    
    # algo for finding local min on left side
    {
      y_vals<-t$Intensity
      print(y_vals)
      m<-which((y_vals)==max(y_vals))
      #if(len(as.data.frame(m)>1)){m<-m[1]}
      print(m)
      if(length(m)>1){m<-m[length(m)]}
      
      
      n<-m-1
      A_m<-y_vals[m]
      A_n<-y_vals[n]
      
      if(A_n <=A_m)
      {
        while(A_n<=A_m && n>0)
        {
          m=m-1
          n=m-1
          
          A_m<-y_vals[m]
          A_n<-y_vals[n]
          
        }
        min_first_half_intense<-y_vals[m]
        min_first_half_RT<-t[which(t$Intensity == min_first_half_intense),]$RT[1]
        
        
        
      }
    }
    # min_first_half_intense<-min(t$Intensity)
    # min_first_half_RT<-t[which(t$Intensity == min_first_half_intense),]$RT
    #
    g<-ggplot(data=t,aes(x=RT,y=Intensity,col=file))+geom_point()+
      geom_vline(xintercept=peak_RT)+
      geom_point(aes(x=min_first_half_RT,y=min_first_half_intense),shape=25,size=5)
    print(g)
    ggsave(g,filename = paste((paste('HPLC_Traces',nam,sep='/')),'3_firsthalfpeak_min.png',sep='/'))
    
    #second half min
    t<-sample[which(sample$RT>peak_RT) ,]
    print(head(t))
    
    g<-ggplot(data=t,aes(x=RT,y=Intensity,col=file))+geom_point()+
      geom_vline(xintercept=peak_RT)
    print(g)
    
    # algo for finding local min on right side
    {
      y_vals<-t$Intensity
      print(y_vals)
      m<-which((y_vals)==max(y_vals))
      # if(len(as.data.frame(m)>1)){m<-m[1]}
      print(m)
      print(length(m))
      if(length(m)>1){m<-m[length(m)]}
      
      n<-m+1
      A_m<-y_vals[m]
      A_n<-y_vals[n]
      
      if(A_n <=A_m)
      {
        while(A_n<=A_m && n<length(y_vals))
        {
          m=m+1
          n=m+1
          
          A_m<-y_vals[m]
          A_n<-y_vals[n]
          
        }
        min_second_half_intense<-y_vals[m]
        min_second_half_RT<-t[which(t$Intensity == min_second_half_intense),]$RT[1]
        
        
        
      }
    }
    
    g<-ggplot(data=t,aes(x=RT,y=Intensity,col=file))+geom_point()+
      geom_vline(xintercept=peak_RT)+
      geom_point(aes(x=min_second_half_RT,y= min_second_half_intense),shape=25,size=5)
    print(g)
    ggsave(g,filename = paste((paste('HPLC_Traces',nam,sep='/')),'4_secondhalfpeak_min.png',sep='/'))
    
    #
    #see whole plot
    g<-ggplot(data=sample,aes(x=RT,y=Intensity,col=file))+geom_point()+
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
    g<-ggplot(data=sample,aes(x=RT,y=Intensity,col=file))+geom_point()+
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
    
    g<-ggplot(data=sample,aes(x=RT,y=Intensity,col=file))+geom_point()+
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
    
    g<-ggplot(data=sample,aes(x=RT,y=CorrectedIntense,col=file))+geom_point()
    print(g)
    ggsave(g,filename = paste((paste('HPLC_Traces',nam,sep='/')),'7_final_corrected.png',sep='/'))
    if(!(dir.exists('./HPLC_Traces/corrected_all/')))
    {
      dir.create('./HPLC_Traces/corrected_all/')
      
    }
    
    
    ggsave(g,filename = paste((paste('HPLC_Traces/corrected_all',nam,sep='/')),'png',sep='.'))
    
    return(sample)
  }
  
 }