
########################################################
#The following script takes the grad_dropout.csv file and
#cleans it to make the data frame mich_grad_dat which looks
#like this:

#   county_code      county_name num_graduating cohort_count graduation_rt percent_male percent_female
# 1           1  michigan,alcona             56           61          91.8         52.0           48.0
# 2           2   michigan,alger             76           85          89.4         55.7           44.3
# 3           3 michigan,allegan           1204         1344          88.8         51.9           48.1
# 4           4  michigan,alpena            259          314          80.9         47.8           52.2
# 5           5  michigan,antrim            227          264          87.9         49.2           50.8
# 6           6  michigan,arenac            125          149          83.9         52.2           47.7
#   percent_econ_disadvan percent_hispanic percent_white percent_black percent_asian
# 1                 100.0              8.2          93.4           0.0           8.2
# 2                  34.1              0.0          84.7           0.0           0.0
# 3                  32.7              7.1          89.9           3.3           1.9
# 4                  41.7              0.0          97.5           0.0           1.6
# 5                  43.6              3.8          94.7           1.9           1.9
# 6                  55.7              6.7          94.6           0.0           0.0

#Columns indicate:
#county_code: unique identifier Michigan uses for each county
#county_name: county_name associated with each county_code
#num_graduating: 
      #Completed high school with a regular diploma in four years or less, or
      #- Early/Middle College participant, completed high school with a regular diploma AND an associate
      #degree, other advanced certificate or up to 60 transferable college credits in five years or less.
      #- The student must be reported with a graduate exit date on or prior to August 31 of his or her cohort year
      #in order to be considered an on-track graduate
#cohort_count: Number of students in the cohort year (2015):
      #Students are placed into a cohort year when they are first identified as ninth graders. Students who transfer
      #into the public education system after ninth grade are placed into the appropriate cohort based on the grade in
      #which the initial Michigan district places them. 
#graduation_rt: % percent of students in the cohort who are graduating
#percent_x: % of the cohort that is x
#percent_econ_disadvan: 
      #A student will be identified in this subgroup if s/he was 1)
      #reported as eligible for supplemental nutrition in any certified collection in the MSDS in the current
      #school year by any reporting entity, 2) directly certified or 3) reported as homeless or migrant. 

#Note that the graduation_rt will not always equal 100*num_graduating/cohort_count due to 
#discrepancies in Michigan's public databases. A row is included for each of Michigan's 
#83 counties. 

#Also note that due to these discrepancies, the mean graduation rate here is inflated to 84-85% 
#rather than the true 79%. This is because the aggregated county-level cohort_count underestimates
#the true statelevel cohort_count that Mighigan records. This difference can't quite be reconciled 
#with data they provide. Regardless, these graduation rates provide a good picture of the spatial
#correlation.

#Data can be found at:
#https://www.mischooldata.org/DistrictSchoolProfiles/StudentInformation/GraduationDropoutRate2.aspx

#More information on Michigan's 4-year graduation rates can be found at:
#https://www.michigan.gov/documents/cepi/Understanding_Michigans_Cohort_Grad-Drop_Rates_2015_496943_7.pdf

#This script produces a .csv file called data_cleaned which is just a .csv version of the
#data frame mich_grad_dat.
########################################################

setwd("~/spatial_project")
#load libraries
library(RColorBrewer)
library(classInt)
library(RgoogleMaps)
suppressMessages(library(fields))
library(akima)
library(maps)
library(mapproj)
library(gstat)
library(geoR)
library(spBayes)
library(fields)
library(MBA)
library(geoRglm)
library(SpatialEpi)
library(CARBayes)
library(stargazer)
library(choroplethr)
library(leaps)
library(maptools)
library(spdep)
library(sqldf)
savepdf <- function(file, width=16, height=10)
{
  fname <- paste("plots/",file,".pdf",sep="")
  pdf(fname, width=width/2.54, height=height/2.54,
      pointsize=10)
  par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.1,1.1))
}
sql<-function(str){
  sqldf()
  ret<-sqldf(str,drv="SQLite")
  sqldf()
  return(ret)
}



grad_dat<-read.csv('data//grad_dropout.csv')

#Eliminate statewide aggregated instances
grad_dat<-subset(grad_dat,CountyCode!=0)

#Only use 2014-2015 student cohort
grad_dat<-subset(grad_dat,CohortYear==2015)

#Change >95 to 97.5 and <5 to 2.5 for analytical purposes
grad_dat$GraduationRate<-as.character(grad_dat$GraduationRate)
grad_dat$DropoutRate<-as.character(grad_dat$DropoutRate)

grad_dat[grad_dat$GraduationRate=='>95%',c('GraduationRate')]<-'97.5'
grad_dat[grad_dat$GraduationRate=='<5%',c('GraduationRate')]<-'2.5'
grad_dat[grad_dat$DropoutRate=='>95%',c('DropoutRate')]<-'97.5'
grad_dat[grad_dat$DropoutRate=='<5%',c('DropoutRate')]<-'2.5'

grad_dat$CohortCount<-as.character(grad_dat$CohortCount)
grad_dat[grad_dat$CohortCount=='<10',c('CohortCount')]<-'5'

grad_dat$OnTrackGraduatedCount<-as.character(grad_dat$OnTrackGraduatedCount)
grad_dat[grad_dat$OnTrackGraduatedCount=='<10',c('OnTrackGraduatedCount')]<-'5'

grad_dat$DropoutCount<-as.character(grad_dat$DropoutCount)
grad_dat[grad_dat$DropoutCount=='<10',c('DropoutCount')]<-'5'

grad_dat$OffTrackContinuingCount<-as.character(grad_dat$OffTrackContinuingCount)
grad_dat[grad_dat$OffTrackContinuingCount=='<10',c('OffTrackContinuingCount')]<-'5'

grad_dat$OtherCompleterCount<-as.character(grad_dat$OtherCompleterCount)
grad_dat[grad_dat$OtherCompleterCount=='<10',c('OtherCompleterCount')]<-'5'

grad_dat<-subset(grad_dat,BuildingName!=0 & Gender=='All Students')

bla<-subset(grad_dat,Subgroup=='All Students') 



temp<-sql("
select
  gd.*
  ,total_students_county
  ,total_students_building
from grad_dat gd inner join
(
  select 
    CountyCode
    ,sum(CohortCount) as total_students_county
  from grad_dat
  where Subgroup='All Students'
  group by 1
) as t
on gd.CountyCode=t.CountyCode
  inner join
(
  select
    BuildingCode
    ,sum(CohortCount) as total_students_building
  from grad_dat
  where Subgroup='All Students'
  group by 1
) as t2
on gd.BuildingCode=t2.BuildingCode

")
  
temp$prop_county<-temp$total_students_building/temp$total_students_county  

building_props<-sql("
select
  BuildingCode
  ,BuildingName
  ,round(sum(case when Subgroup='Male' then CohortCount else 0 end),3)/total_students_building as prop_male
  ,round(sum(case when Subgroup='Female' then CohortCount else 0 end),3)/total_students_building as prop_female
  ,round(sum(case when Subgroup='Economically Disadvantaged' then CohortCount else 0 end),3)/total_students_building as prop_econ_disadvan
  ,round(sum(case when Subgroup='Hispanic' then CohortCount else 0 end),3)/total_students_building as prop_hispanic
  ,round(sum(case when Subgroup='White, not of Hispanic origin' then CohortCount else 0 end),3)/total_students_building as prop_white
  ,round(sum(case when Subgroup='Black, not of Hispanic origin' then CohortCount else 0 end),3)/total_students_building as prop_black
  ,round(sum(case when Subgroup='Asian' then CohortCount else 0 end),3)/total_students_building as prop_asian
from temp
group by 1,2
")

temp2<-sql("
select 
  t.CountyCode
  ,t.CountyName
  ,t.BuildingCode
  ,t.BuildingName
  ,t.prop_county
  ,round(prop_male,2)/(prop_male + prop_female) as prop_male
  ,round(prop_female,2)/(prop_male + prop_female) as prop_female
  ,prop_econ_disadvan
  ,prop_hispanic
  ,prop_white
  ,prop_black
  ,prop_asian
  ,sum(OnTrackGraduatedCount) as num_graduating
  ,sum(CohortCount) as cohort_count
  ,max(GraduationRate) as graduation_rt
from temp t 
  inner join building_props bp on t.BuildingCode=bp.BuildingCode
where Subgroup='All Students' 
group by 1,2,3,4
")



temp3<-sql("
select
  CountyCode as county_code
  ,CountyName as county_name
  ,sum(num_graduating) as num_graduating
  ,sum(cohort_count) as cohort_count
  ,round(sum(prop_county*graduation_rt),1) as graduation_rt
  ,round(100*sum(prop_county*prop_male),1) as percent_male
  ,round(100*sum(prop_county*prop_female),1) as percent_female
  ,round(100*sum(prop_county*prop_econ_disadvan),1) as percent_econ_disadvan
  ,round(100*sum(prop_county*prop_hispanic),1) as percent_hispanic
  ,round(100*sum(prop_county*prop_white),1) as percent_white
  ,round(100*sum(prop_county*prop_black),1) as percent_black
  ,round(100*sum(prop_county*prop_asian),1) as percent_asian
from temp2
group by 1,2
")

#Here, we'll impute Keweenaw County with means
newrow<-c(42,'Keweenaw',
  round(mean(as.numeric(temp3$num_graduating)),1),
  round(mean(as.numeric(temp3$cohort_count)),1),
  round(mean(as.numeric(temp3$graduation_rt)),1),
  round(mean(temp3$percent_male),1),
  round(mean(temp3$percent_female),1),
  round(mean(temp3$percent_econ_disadvan),1),
  round(mean(temp3$percent_hispanic),1),
  round(mean(temp3$percent_white),1),
  round(mean(temp3$percent_black),1),
  round(mean(temp3$percent_asian),1))

temp3<-rbind(temp3,newrow)


temp3$num_graduating<-as.numeric(temp3$num_graduating)
temp3$cohort_count<-as.numeric(temp3$cohort_count)
temp3$graduation_rt<-as.numeric(temp3$graduation_rt)
temp3$percent_male<-as.numeric(temp3$percent_male)
temp3$percent_female<-as.numeric(temp3$percent_female)
temp3$percent_econ_disadvan<-as.numeric(temp3$percent_econ_disadvan)
temp3$percent_hispanic<-as.numeric(temp3$percent_hispanic)
temp3$percent_white<-as.numeric(temp3$percent_white)
temp3$percent_black<-as.numeric(temp3$percent_black)
temp3$percent_asian<-as.numeric(temp3$percent_asian)
temp3$county_code<-as.integer(temp3$county_code)



mich_grad_dat<-sql("select * from temp3 order by county_code asc")

mich_grad_dat$county_name<-paste('michigan,',tolower(mich_grad_dat$county_name),sep='')


write.csv(mich_grad_dat,file='data//data_cleaned.csv',row.names=FALSE)
