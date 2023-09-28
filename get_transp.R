{
if (!require("rstudioapi")) install.packages("rstudioapi")
library("rstudioapi") 
setwd(dirname(getActiveDocumentContext()$path))
.libPaths('../lib')
library(rtransparent)
library(oddpub)
library(metareadr, lib.loc='./lib')
library(beepr)
library(plyr)
library(dplyr)
library(tibble)
library(parallel)
library(doParallel)
library(ggplot2)
library(reshape2)
library(ddpcr)
library(stringr, lib.loc='/Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/library')

#install.packages('SparkR')

}


downloads= function(loc){
  pmcidfilename=paste0('./',loc,".csv")
  pmcidlist<-read.delim(pmcidfilename, header = TRUE, sep=',')
  pmcidlist=pmcidlist$PMCID
  pmcnumber<-list()
  for (i in pmcidlist){
    go=str_replace(i,'PMC','')
    pmcnumber=c(pmcnumber,go)
  }  
  
  
  filenames=paste0('./publications', '/PMC',as.character(pmcnumber),'.xml')
  
  mapply(metareadr::mt_read_pmcoa,pmcid=pmcnumber,file_name=filenames)
  
}
downloads('pmc')
mclapply('pmc',downloads)


checkdiff= function(loc,loc1){
  filelist <- list.files(paste0('./',loc,'/'), pattern='*.xml', all.files=FALSE, full.names=FALSE)
  pmcidfilename=paste0("./",loc1,".csv")
  pmcidlist<-read.delim(pmcidfilename, header = TRUE, sep=',')
  pmcidlist=pmcidlist$PMCID
  pmcnumber<-list()
  for (i in pmcidlist){
    go=str_replace(i,'PMC','')
    pmcnumber=c(pmcnumber,go)}
  downloaded=str_remove(filelist,'PMC')
  downloaded=str_remove(downloaded,'.xml')
  return(setdiff(pmcnumber, downloaded))
  
}


downloadspmc=function(pmcnumber){
  
  filenames=paste0('./publications/PMC',as.character(pmcnumber),'.xml')

  mapply(metareadr::mt_read_pmcoa,pmcid=pmcnumber,file_name=filenames)
}

beep()

pmcnumber=checkdiff('publications','pmc')

downloadspmc(pmcnumber)


library(plyr)

dohooptyhoop=function(loc){
filepath='./publications/'
filelist <- as.list(list.files(filepath, pattern='*.xml', all.files=FALSE, full.names=FALSE))
                        
filelist=paste0(filepath, filelist)
#filelist=tail(filelist,10)                    
cores <- detectCores()
registerDoParallel(cores=cores)
                        
return(foreach::foreach(x = filelist,.combine='rbind.fill') %dopar%{
## Use the same library paths as the master R session
#.libPaths(libs[1])
rtransparent::rt_data_code_pmc(x)
})}
dohooptyhooprest=function(loc){
  filepath='./publications/'
  filelist <- as.list(list.files(filepath, pattern='*.xml', all.files=FALSE, full.names=FALSE))
  
  filelist=paste0(filepath, filelist)
  #filelist=tail(filelist,10)                    
  cores <- detectCores()
  registerDoParallel(cores=cores)
  
  return(foreach::foreach(x = filelist,.combine='rbind.fill') %dopar%{
    ## Use the same library paths as the master R session
    #.libPaths(libs[1])
    rtransparent::rt_all_pmc(x)
  })}
rbind.fill()

library(beepr)

code_df=dohooptyhoop() ; beep()
beep()
write.csv(code_df,"./codesharing.csv", row.names = FALSE); beep()

rest_df=dohooptyhooprest() 
write.csv(rest_df,"./resttransp.csv", row.names = FALSE); beep()
