# Key data and analyses of the sample
library(ggplot2)
library(dplyr)
library(DescTools)
library(scales)
library(tidyr)
library(tidylog) # Tidylog provides feedback about dplyr and tidyr operations.
library(stringr)


# LOAD DATA
code_df=read.delim('./codesharing.csv', sep=',')
rest_df=read.delim('./resttransp.csv', sep=',')
combined=merge(x = code_df, y = rest_df, by = "pmid", all = TRUE)
df=read.csv('./pmc_final.csv')

df%>%filter(Publication.Year>2012)%>%nrow()/nrow(df)
median(df$Publication.Year)
quantile(df$Publication.Year,probs=c(0.25,0.75))

table(df$Journal.Book)

sort(table(df$Journal.Book),decreasing=TRUE)[1:5]


comb=merge(x = combined, y = df, by = "pmid", all = TRUE)
comb[is.na(comb)] <- FALSE  

# loading the transparency data
df_val=readxl::read_excel('./Extraction_RKB.xlsx')
exclude=df_val%>%filter(Exclude=='y')%>%select(DOI)
comb=comb%>%filter(!DOI%in%exclude$DOI)

df_val=df_val%>%filter(is.na(Exclude)) 

# Define indicator variables
indicator_vars=c('is_open_code','is_open_data','is_register_pred','is_coi_pred','is_fund_pred')


# RESULTS
## OVERALL TRANSPARENCY
ovr_transp=comb%>% 
  select(!!!indicator_vars)%>%
  dplyr::summarise_all(~sum(.x[.x==TRUE]))
ovr_transp/nrow(comb) 


### Manual validation
# validation section checking congruence of validation to algorithm results
# gets you FP numerator and denominator
# gets you FN numerator and denominator
df_val=(df_val%>%mutate(across(indicator_vars[1:3],~as.factor(.x))))
df_val=(df_val%>%mutate(across(c('code_val','data_val','register_val'),~as.factor(.x))))

validationfun=function(indicator, df_val){
  if (indicator=='code' |indicator== 'data'){indicatorz=paste0('is_open_',indicator)}else{indicatorz=paste0('is_',indicator,'_pred')}
  valcol=paste0(indicator,'_val')
  print(valcol)
  print(indicatorz)
  k=df_val %>%
    group_by(get(indicatorz), get(valcol)) %>%
    count() %>%
    group_by(`get(indicatorz)`)%>%
    dplyr::mutate(total=sum(n))%>%
    mutate(p=n/total)#%>%
    #rename(algorithm=1,validation=2)
  
  
  return(k)
}



# adjusted indicators FP and FN rates
validationfun('code', df_val)
validationfun('data', df_val)
validationfun('register', df_val)

# U(i) × TP + (1 − U(i)) × FN
# insert the obtained values (get indicatorz is the original
# classification get valcol is the corrected)
# corrected code
0.03209877*0.857+(1-0.03209877)*0.0217

# corrected data
0.2123457*0.974+(1-0.2123457)*0.085

# corrected registration
0.008230453*1+(1-0.008230453)*0
  

table(df_val$`Details about COI (actual disclosures)`)
table(df_val$`Details about funding (actual disclosures)`)

# MAIN COMPARISONS
## GRAPHING OVER TIME

# counts the TRUE values per group of year and transparency indicator
total=comb%>%group_by(Publication.Year)%>%count()

indicator_by_year <- 
  comb %>% 
  group_by(Publication.Year)%>%
  # total per year
  dplyr::mutate(total=n())%>%
  select(Publication.Year,total, !!!indicator_vars) %>% 
  # elongate your dataframe
  gather("key", "value", -c(Publication.Year,total)) %>% 
  group_by(Publication.Year,key)%>%
  # gives you total and n for the transparent articles
  dplyr::reframe(n=sum(value[value==TRUE]), total=first(total)) %>%
  #filters the transparent ones, gets proportion of them tot total and Wilson CI
  filter(total>20)%>%
  group_by(Publication.Year, key)%>%
  dplyr::summarise(p=n/total, n=n,total=first(total))%>%
  mutate(lwr.ci=BinomCI(n, total, conf.level=0.95, method='wilson')[,2], 
         upr.ci=BinomCI(n, total, conf.level=0.95, method='wilson')[,3])
  
# excluding reviews TEST
indicators_by_year <- 
  comb %>% 
  filter(is_review.x==FALSE)%>%
  group_by(Publication.Year)%>%
  # total per year
  mutate(total=n())%>%
  select(Publication.Year,total, !!!indicator_vars) %>% 
  # elongate your dataframe
  gather("key", "value", -c(Publication.Year,total)) %>% 
  group_by(Publication.Year,key)%>%
  # gives you total and n for the transparent articles
  dplyr::reframe(n=sum(value[value==TRUE]), total=first(total)) %>%
  #filters the transparent ones, gets proportion of them tot total and Wilson CI
  filter(total>20)%>%
  group_by(Publication.Year, key)%>%
  dplyr::summarise(p=n/total, n=n,total=first(total))%>%
  mutate(lwr.ci=BinomCI(n, total, conf.level=0.95, method='wilson')[,2], 
         upr.ci=BinomCI(n, total, conf.level=0.95, method='wilson')[,3])

# filter for the boundaries you've set as reasonable

tl_comb=comb%>%filter(Publication.Year>2012)%>%filter(Publication.Year<2023)

# Get p-values for fit over time
for(i in indicator_vars){
  print(i)
  p_val=summary(glm(tl_comb[[i]] ~ tl_comb$Publication.Year, family = 'binomial'))$coefficients[,4][2]
  print(summary(glm(tl_comb[[i]] ~ tl_comb$Publication.Year, family = 'binomial'))$coefficients[,4][2])
  print(p_val<0.005)
}

for(i in indicator_vars){
  print(i)
  p_val=summary(glm(comb[[i]] ~ comb$Publication.Year, family = 'binomial'))$coefficients[,4][2]
  print(summary(glm(comb[[i]] ~ comb$Publication.Year, family = 'binomial'))$coefficients[,4][2])
  print(p_val<0.005)
}



indicators_by_year <- 
  df_val %>% 
  #filter(is_review.x==FALSE)%>%
  group_by(`Publication year`)%>%
  # total per year
  mutate(total=n())%>%
  select(`Publication year`,total, initiator) %>% 
  # elongate your dataframe
  gather("key", "value", -c(`Publication year`,total)) %>% 
  group_by(`Publication year`,key)%>%
  # gives you total and n for the transparent articles
  dplyr::reframe(n=count(value=='University'), total=first(total)) %>%
  #filters the transparent ones, gets proportion of them tot total and Wilson CI
  filter(total>20)%>%
  group_by(`Publication year`, key)%>%
  dplyr::summarise(p=n/total, n=n,total=first(total))%>%
  mutate(lwr.ci=BinomCI(n, total, conf.level=0.95, method='wilson')[,2], 
         upr.ci=BinomCI(n, total, conf.level=0.95, method='wilson')[,3])

for(i in indicator_vars){
  print(i)
  p_val=summary(glm(comb[[i]] ~ comb$Publication.Year, family = 'binomial'))$coefficients[,4][2]
  print(summary(glm(comb[[i]] ~ comb$Publication.Year, family = 'binomial'))$coefficients[,4][2])
  print(p_val<0.005)
}

# Gives us p<0.005 for code sharing + conflict of interest. 
# hence the stars added in legend

# Plot the transparency over time
p2 <- 
  indicator_by_year %>% 
  ggplot(aes(x = Publication.Year, y = p, color = key, fill=key)) +
  geom_point(size = 2) +
  geom_smooth(method='glm', method.args=list(family="binomial"), se=FALSE)+  
  scale_color_manual(labels = c("COI statement*", "Funding statement",'Code sharing*','Data sharing','Registration'),values=c('red','blue','black','forestgreen','orange')) +
  scale_fill_manual(labels = c("COI statement*", "Funding statement",'Code sharing*','Data sharing','Registration'),values=c('red','blue','black','forestgreen','orange'))+
  geom_ribbon(aes(ymin=lwr.ci, ymax=upr.ci), alpha=0.2, linetype=0)+
  coord_cartesian()+
  scale_x_continuous(breaks = seq(2013, 2022, 1), limits = c(2013, 2022)) +
  scale_y_continuous(limits = c(-0.1, 1), labels = scales::percent) +
  #scale_color_discrete(drop = FALSE) +
  guides(fill='none') +
  theme_bw()+
  theme(
    panel.grid.minor = element_blank(), 
    legend.position = c(0.83, 0.61),
    legend.title=element_text(size=9,face='bold')
  )+
  labs(y = "Proportion of articles (%)\n", x = "\nYear", color='Transparency indicator\n(p<0.005 for trend*)')

# create figure
tiff('./figures/transparency_time.tiff', units="in", width=5, height=4, res=300, compression = 'lzw')

p2
dev.off()

init=df_val
init[init['initiator']!='University',]['initiator']='FALSE'
init[init['initiator']!='FALSE',]['initiator']='TRUE'
init=init%>%mutate(initiator=as.factor(initiator))
init['Publication.Year']=init['Publication year']
summary(glm(initiator~Publication.Year, data=init,family='binomial'))

fund=df_val%>%filter(!is.na(funder))
fund[fund['funder']!='Public institution',]['funder']='FALSE'
fund[fund['funder']!='FALSE',]['funder']='TRUE'
fund=fund%>%mutate(funder=as.factor(funder))
fund['Publication.Year']=fund['Publication year']
summary(glm(funder~Publication.Year, data=fund,family='binomial'))

code_aov <- aov(code_val ~ `Sample size` + articlelength +`Number of tables`+`Number of figures`+`Number of appendices/supplementary files`, data = df_val)%>%summary()
data_aov <- aov(data_val ~ `Sample size` + articlelength+`Number of tables`+`Number of figures`+`Number of appendices/supplementary files`, data = df_val)%>%summary()
reg_aov <- aov(register_val ~ `Sample size` + articlelength+`Number of tables`+`Number of figures`+`Number of appendices/supplementary files`, data = df_val)%>%summary()
coi_aov <- aov(`Conflict of interest statement` ~ `Sample size` + articlelength+`Number of tables`+`Number of figures`+`Number of appendices/supplementary files`, data = df_val)%>%summary()
fund_aov <- aov(`Funding disclosure` ~ `Sample size` + articlelength+`Number of tables`+`Number of figures`+`Number of appendices/supplementary files`, data = df_val)%>%summary()

code_aov
data_aov
reg_aov
coi_aov
fund_aov


## MANUALLY EXTRACTED CHARACTERISTICS

# Manual validation

# Correlates of transparency in the manually assessed sample
## possibly split into meta-data, sample characteristics, study design...?

## Characteristics of Extraction

# All categorical variables overall transparency
cat_corr=c('design','setting','interventional','funder', 'datatype','initiator','population','pvalabs', 'umeasure')
vars_manual=c('code_val','data_val','register_val','Conflict of interest statement','Funding disclosure')

for(i in cat_corr){
  x=df_val%>%
    group_by(get(i))%>%
    summarise(across(vars_manual,~paste0(sum(.x[.x==TRUE]),
                                       '/',n(),' (',
                                       round(sum(.x[.x==TRUE])/n(),2)*100,
                                       ')')))
  if (!exists("cat_corr_transp")) {cat_corr_transp <- x}
  else {cat_corr_transp <- rbind(cat_corr_transp, x)}
  #x=rbind(x,x)
  
}

l=(nrow(df_val))
# gets you total
for(i in cat_corr){
  x=df_val%>%
    #select(get(i),!!!vars_manual)%>%
    group_by(get(i))%>%
    #mutate_at(test=sum(code_val[code_val=='TRUE']))
    tally()%>%
    mutate(Total=paste0(n,'/',l, '(',round(n/as.numeric(l),2)*100,')'))
  
  
  if (!exists("cat_corr_tot")) {cat_corr_tot <- x}
  else {cat_corr_tot <- rbind(cat_corr_tot, x)}
  
}

cat_corr=cbind(cat_corr_tot, cat_corr_transp[,2:6])

write.csv(cat_corr, './tables/categorical_correlates_transparency.csv')

# FOR NUMERICAL CORRELATES
# get the median(IQR) for the ones that are transparent according to our indicators
num_corr_transp=data.frame()
num_vars=c('Sample size','articlelength',
           'Number of tables','Number of figures',
           'Number of appendices/supplementary files')
df_vals=df_val%>%mutate(across(num_vars,~as.numeric(.x)))
num_val=df_vals%>%select(!!!vars_manual,!!!num_vars)

num_corr_transp=data.frame()
for(i in vars_manual){
  x=num_val%>%
    filter(get(i)==TRUE)%>%
    summarise(across(num_vars,~paste0(median(.x,na.rm=TRUE),
                                      ' (',paste(
                                        quantile(.x,probs=c(0.25,0.75), na.rm=TRUE),
                                        collapse='-'),')')))
  if (!exists("num_corr_transp")) {num_corr_transp <- x}
  else {num_corr_transp <- rbind(num_corr_transp, x)}
  #x=rbind(x,x)
  
  
}

num_corr_transp=t(num_corr_transp)
colnames(num_corr_transp)=vars_manual

overall_num=num_val%>%
  summarise(across(num_vars,~paste0(median(.x,na.rm=TRUE),
                                    ' (',paste(
                                      quantile(.x,probs=c(0.25,0.75), na.rm=TRUE),
                                      collapse='-'),')')))


overall_num=t(overall_num)  
colnames(overall_num)='Overall'
num_corr_transp=cbind(overall_num,num_corr_transp)

write.csv(num_corr_transp,'./tables/numerical_correlates_transparency.csv')

# to see all the very transparent articles there are none with all five but a few with these four
comb%>%filter(is_open_code==TRUE & 
                  is_open_data==TRUE & 
                  #is_register_pred==TRUE & 
                  is_coi_pred==TRUE &
                  is_fund_pred==TRUE)%>%
  write.csv('./tables/Very transparent articles.csv')


comb%>%filter(is_register_pred==TRUE)%>%write.csv('./tables/registered articles.csv')

