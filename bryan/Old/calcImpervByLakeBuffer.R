v='calcImpervByLakeBuffer'

#########load packages
libs<-c("sp","rgeos","rgdal","raster","maptools") #list of packages to load
  installLoad<-function(pck)#user defined function
    {
      if(!pck%in%installed.packages()){install.packages(pck,repos="http://rweb.quant.ku.edu/cran/")}
      require(pck, character.only = TRUE)
    }
lapply(libs,function(x) installLoad(x))  #Load/Install require packages
#########
#get the list of COMIDs and their corresponding NLA_IDs
  ID<- read.csv(url('https://raw.github.com/jhollist/cyanoLakes/master/nlalm/nlaid_comid_lut.csv'))

#get the NLCD grid data
  Imperv<-raster('L:/Public/Milstead_Lakes/NLCD2006/nlcd_2006_impervious_2011_edition_2014_03_31.img')
  #image(Imperv)
  #extent(Imperv)

##############################
####function to calculate the impervious cover area and percent in a buffer.
calcImperv1<- function(Buf,Grid){ 
  gc() #release unused memory 
  Mask<-mask(Grid, Buf)  #extract impervious cover pixels in the buffer (all others changed to NA)
  a<-table(getValues(Mask),useNA='ifany')  #the NA are pixels outside of bbox of buffer.  Values of "127" are the real NA
  a<-na.exclude(data.frame(Percent=as.numeric(names(a)),a)[,-2]) #replace values stored as factors with numbers
  a$ImpPix<-a$Percent*a$Freq/100 #processing step-convert frequency distribution to number of impervious pixels
  MissingPix<-ifelse(max(a$Percent)==127,a[a$Percent==127,'Freq'],0)  #number of missing pixels
  TotalPix<-sum(a$Freq)  #total number of pixels in buffer
  PercentNA<-round(MissingPix/TotalPix,3)*100  #percent NA cells in buffer
  BufAreaKm2<-round(sum(a$Freq)*30*30/1000000,3)  #total area of buffers (based on pixel number)-includes NA (value=127)
  a<-subset(a,a$Percent<=100)  #remove pixels with missing values (=127)
  BufAreaKm2Adj<-round(sum(a$Freq)*30*30/1000000,3)  #this doesn't include Pixels with missing data
  ImpervAreaKm2<-round(sum(a$ImpPix)*30*30/1000000,3) #calculate the area of impervious cover
  PercentImperv<-round(sum(ImpervAreaKm2)/sum(BufAreaKm2Adj),4)*100 #calculate the percent of impervious cover 
  data.frame(PercentImperv,ImpervAreaKm2,BufAreaKm2,BufAreaKm2Adj,PercentNA)
}
####function to calculate the impervious cover area and percent in four buffers.
####################Function to calculate the impervious surface around the spatialPointsDataFrame 'Lake'
                    #for 4 buffer distances. 
  #Index is the line number for the lake in the ID df
  #PlotYN enter 'Y' to generate figure or 'N' for text only
calcImperv<-function(Index,PlotYN='Y'){ 
  Start<-Sys.time()  #Record Start Time
  Comid<-ID$COMID[Index]      #Comid of chosen lake (from the df 'ID')
  NLAID<-ID$nlaSITE_ID[Index] #NLA id of chosen lake
  File<-paste('lakemorpho',Comid,sep='')  #identify the df that corresponds to the comid
  #get the spatial data object for the chosen lake
    load(url(paste('https://raw.github.com/jhollist/cyanoLakes/master/nlalm/',File,'.RData',sep=''))) 
    Lake<-get(File)$lake  #write lake polygon to new file
  #define lake buffers (in Meters) to be used.  
    Buffer<-c(NA,300,1500,3000) #first buffer will be equal to max in lake distance caluclated below.
    Buffer[1]<-round(max(get(File)$lakeDistance@data@values,na.rm=T),0)
  #remove the spatial object to clear memory
    rm(list=File)  
    gc() #release unused memory 
  #close connections-r keeps the connection to gitHub open and will eventually be unhappy.
    closeAllConnections() 
  #Associate lake holes (islands) with the correct.  
    slot(Lake, "polygons") <- lapply(slot(Lake, "polygons"), checkPolygonsHoles)
  #reproject lake spdf to match the NLCD grid
    #ESRI USA_Contiguous_Albers_Equal_Area_Conic_USGS_version
    AlbersContigUSGS<-CRS('+proj=aea +x_0=0 +y_0=0 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 
                             +lat_0=23 +units=m +datum=NAD83') 
    Lake<-spTransform(Lake,AlbersContigUSGS)   
  #Imperv is the impervious cover for all lower 48 states.  Reduce to bbox(Lake) + max(buffer) + 1000m
    a<-max(Buffer)+1000 #maximum buffer distance + 1000m
    B<-bbox(Lake) #bounding box for the Lake 
    Extent<-c(B[1,1]-a,B[1,2]+a,B[2,1]-a,B[2,2]+a)
      #Crop grid data for the Largest Lake Buffer + 1000m
        Grid<-crop(Imperv,Extent)
  #Get the buffers (without the lake)
    Buf1<-gDifference(gBuffer(Lake,width=Buffer[1]),Lake)
    Buf2<-gDifference(gBuffer(Lake,width=Buffer[2]),Lake)
    Buf3<-gDifference(gBuffer(Lake,width=Buffer[3]),Lake)
    Buf4<-gDifference(gBuffer(Lake,width=Buffer[4]),Lake)
  #plot grid and buffers    ###optional PlotYN='Y' to plot; PlotYN='N' to return text only.
    if(PlotYN=='Y'){
      image(Grid,main=paste('NLA ID = ',NLAID))
      plot(Buf1,add=T)                  
      plot(Buf2,add=T)
      plot(Buf3,add=T)
      plot(Buf4,add=T)
      title(sub=v)
    }
  #calculate imperv area and percent for each buffer
    a<-data.frame(matrix(nrow=4,ncol=8)) 
    a[1,]<-c(Comid,NLAID,Buffer[1],calcImperv1(Buf1,Grid)) 
    a[2,]<-c(Comid,NLAID,Buffer[2],calcImperv1(Buf2,Grid)) 
    a[3,]<-c(Comid,NLAID,Buffer[3],calcImperv1(Buf3,Grid)) 
    a[4,]<-c(Comid,NLAID,Buffer[4],calcImperv1(Buf4,Grid)) 
  #rename the output columns
    names(a)<-c('Comid','NLA_ID','BufWidthM','PercentImperv','ImpervAreaKm2','BufferAreaKm2',
                  'BufferAreaKm2Adj','PercentNA')
  #print elapsed time
    print(paste(Comid,' Time Elapsed = ',round(Sys.time()-Start),' seconds')) #how long did the process take?
    print('') #extra line
  #return results
    return(a)
}    

calcImperv(1,'N')
calcImperv(2,'Y')



##loop to run all NLA lakes in the df 'ID'
  #create df to store the output
    lakeImperv<-data.frame(matrix(NA,nrow=nrow(ID)*4,ncol=8))  
      names(lakeImperv)<-c('Comid','NLA_ID','BufWidthM','PercentImperv','ImpervAreaKm2','BufferAreaKm2',
                              'BufferAreaKm2Adj','PercentNA')
  #constants for loop and for saving
    B<-1  #row number to start processing lakes; this will be 1 to start.
    N<-nrow(ID)    #last row to process; 
    #N<-3    #last row to process; 
    S<-100                 #save the work every "S" lakes
    Counter<-0             #counter; should be zero to start
  #start loop
    for(i in c(B:N)){
      #calc the imperv for each lake ID[i,]
        tryCatch(lakeImperv[(4*i-3):(i*4),]<-calcImperv(i,'Y'),error=function(e) NA) 
          lakeImperv[(4*i-3):(i*4),1:2]<-ID[i,c(2,4)] #re-add the IDs to find problem lakes
          Counter<-Counter+1
            if (Counter==S) save(lakeImperv, #save file every S iterations
                file='C:/Bryan/PortableApps/R/scripts/cyanoLakes/bryan/lakeImperv.rda') 
            if(Counter==S) Counter<-0
            if (i==N) save(lakeImperv, #save file at end of loop
                file='C:/Bryan/PortableApps/R/scripts/cyanoLakes/bryan/lakeImperv.rda') 
    }#End Loop

lakeImperv[1:16,]

#############old below

#add field to distinguish between fixed width and variable (radius) width buffers.  
#some radius buffers the same as the fixed width.
lakeImperv$BufType<-c('MaxDist','fixed','fixed','fixed')
#table(lakeImperv$BufType)
#table(lakeImperv[lakeImperv$BufWidthM==300,'BufType'])

#save the data
#save(lakeImperv,file="L:/Public/Milstead_Lakes/RData/NLCD2006ImpervLakes_20130212.rda")  #save the last records.
#load("L:/Public/Milstead_Lakes/RData/NLCD2006ImpervLakes_20130212.rda")
#write.table(lakeImperv, file='c:/temp/lakeImperv.csv',row.names=F,sep=',')

# Write data to tblMNKA in fj.mdb
require(RODBC)   #Package RODBC must be installed
con <- odbcConnectAccess("C:/Bryan/EPA/Data/WaterbodyDatabase/WaterbodyDatabase.mdb")
FJ <- sqlSave(con,dat=lakeImperv,tablename='LakeNLCD2006ImperviousCover',append=F,rownames=F)
close(con)


##Data Definitions  n=28,121 (Lake Champlain excluded due to size and complexity)
#WB_ID:  Lake identification number
#BufWidthM: (m) width of lake buffer for calculations
#PercentImperv:  (%) Percent impervious cover in buffer = ImpervAreaKm2/BufferAreaKm2Adj/100
#ImpervAreaKm2:  (Km2) total impervious cover in the buffer
#BufferAreaKm2:  (Km2) Area of entire buffer-number of grid cells*30*30/1000000 
#BufferAreaKm2Adj (Km2) Area of Non-NA buffer-number of grid cells with data*30*30/1000000 
#PercentNA:  (%) (number of NA grid cells)/(total number of grid cells)*100
#BufType: "Fixed" = standard buffer width; "Radius" buffer width = radius of a circle with Area = lake area.



