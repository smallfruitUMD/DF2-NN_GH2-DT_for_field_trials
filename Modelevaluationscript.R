#This script imports csv files containing leaf wetness duration (LWD) and temperature (Temp) values at 5 minute increments and runs the data through two models, DF2-NN and GH2-DT


#load packages
if (!require(pacman)) install.packages("pacman")
pacman::p_load(tidyverse,
               readr,
               knitr,
               zoo,
               cowplot,
               here,
               kableExtra,
               update = FALSE)




#What is the current phenological stage? pick: pre-berry touch, berry touch (BT to 80% color change), veraison (80% color change to 3-4 weeks to harv), or pre-harvest (3-4 weeks out until harvest) 
growth.stages <- c("pre-berry.touch", "berry.touch", "veraison", "pre-harvest")
stage.risks <- c(0.043,0.133,0.821,0.961)
gstagetable <- data.frame(growth.stages, stage.risks)

#PUT CURRENT STAGE HERE!!!
current.phenological.stage <- "pre-harvest"

current.stage.filtered <- gstagetable %>%
  filter(growth.stages == current.phenological.stage)

current.pheno.risk <- current.stage.filtered[1,2] # this is the current risk to ripe rot based on the phenological stage of the clusters



Wetness.cutoff <- 460 # the recommended cutoff from the LWD sensor of what is determined as wet. 



############## FIELD 1 ############################################
###Load field 1 data
fresh.download.WRV <- read_csv(here("example_datasetrscript.csv"))

####Is data formatted correctly?
#The data should be recorded every 5 minutes and it should be checked first before continuing (TRUE means we're okay). Also note that the first row is always junk.
ifelse(fresh.download.WRV[3,1] - fresh.download.WRV[2,1] == 5, "Yes", FALSE)

# Averaged the LWD data because multiple LWD sensors were deployed in this field
counts.averaged <- fresh.download.WRV %>%
  group_by(DateTime) %>%
  dplyr::mutate(avg_count = mean(c(P6_Counts,P4_Counts))) %>% 
  ungroup()




### Turning leaf wetness values from every 5 minutes into LWD
fresh.sensors6 <- counts.averaged %>%
  dplyr::mutate(is_wet = ifelse(avg_count>=Wetness.cutoff,TRUE,FALSE)) %>% 
  dplyr::mutate(Wet_hr = ifelse(is_wet == TRUE,5/60,0)) %>% 
  dplyr::mutate(is_wet_previous48rows = ifelse(is_wet == TRUE, TRUE, ifelse(rollsumr(Wet_hr, k = 48, fill = NA)>=0.1,TRUE,FALSE))) %>% 
  dplyr::mutate(Group_TFchange = 0) %>% # a new column created for a loop below
  dplyr::mutate(Wet_TFchange = 0) %>% # for loop
  dplyr::mutate(Wet_TFchangelength=0)

fresh.sensors6 <- fresh.sensors6[2:nrow(fresh.sensors6),]

n <- 1
fresh.sensors6[48, "Group_TFchange"] <- 1
for(i in 49:nrow(fresh.sensors6)){
  
  if(fresh.sensors6[i-1, "is_wet_previous48rows"] != fresh.sensors6[i, "is_wet_previous48rows"]){
    n=n+1 #if it changes from false to true or true to false as this iterates, then add 1 to the column
  }
  fresh.sensors6[i, "Group_TFchange"] <- n
}

n <- 1
fresh.sensors6[1, "Wet_TFchange"] <- 1
for(i in 2:nrow(fresh.sensors6)){
  
  if(fresh.sensors6[i-1, "is_wet"] != fresh.sensors6[i, "is_wet"]){
    n=n+1 #if it changes from false to true or true to false as this iterates, then add 1 to the column
  }
  fresh.sensors6[i, "Wet_TFchange"] <- n
}

for(i in 1:nrow(fresh.sensors6)){
  fresh.sensors6[i,"Wet_TFchangelength"] <- sum(fresh.sensors6$Wet_TFchange == fresh.sensors6$Wet_TFchange[i])} 

df.with.infriskGHmodel6 <- fresh.sensors6 %>% 
  group_by(Group_TFchange) %>%  
  mutate(cumsum_Wethr = cumsum(Wet_hr)) %>% #calculate cumulative sum of Wetness hours in each group
  mutate(Wetmax24 = ifelse(cumsum_Wethr>24,24,cumsum_Wethr)) %>% #model was only created with max of 24 hours
  dplyr::mutate(drygaplength = ifelse(is_wet == TRUE, 0, Wet_TFchangelength)) %>%
  dplyr::mutate(LWD = ifelse(drygaplength > 48, Wet_hr, Wetmax24))
  





 ###########GH2-DT model function
  riskfun <- function(LWD, P0_Temp) {         
    infsevrisk <- if(LWD >= 12) { if(P0_Temp >= 22) {if(P0_Temp < 32) {if(P0_Temp >= 27) {if(LWD >= 24) {0.951891614010398} else {0.869729502372144}} else {if(LWD < 18) {0.792281699216377} else {if(LWD < 24) {0.694639118800699} else {0.592852158400932} } } } else {  if(LWD >= 24) {0.787267174337893} else {0.577759247720404} } } else {if(LWD >= 24) { if(P0_Temp < 17) {0.555210636351001}   else {0.205955012972572}  } else {if(P0_Temp >= 17) {if(LWD >= 18) {0.488539106004448} else {0.202824820290162} } else { if(LWD < 18) {0.343484160949503} else {0.057769875235217} } } } } else { if(LWD >= 6) { if(P0_Temp >= 22) {if(P0_Temp < 27) {0.488123722737249} else { if(P0_Temp >= 32) {0.48742666720062} else {0.235331111734057} } } else {0.201449841845058} } else{ if(P0_Temp < 17) {0.232500695949199} else { if(P0_Temp >= 32) {0.227140980894363} else {0.0259175632404413} } } }
    return(infsevrisk)
  }

#Creating a new df with blank column. Then loop the function over the column            
df.with.infriskGH2DT.6 <- df.with.infriskGHmodel6  %>%
  mutate(infsevrisk = 0)

for(i in 1:nrow(df.with.infriskGH2DT.6)){
  WD <- df.with.infriskGH2DT.6$LWD[i]
  Temp <- df.with.infriskGH2DT.6$P0_Temp[i]
  n = riskfun(WD,Temp)
  df.with.infriskGH2DT.6[i,"infsevrisk"] <- n
}


############ Calculating of final ripe rot risk determined from the GH2-DT model, which is the GH2-DT risk multiplied by the risk due to the current phenological stage
df.with.infriskGHmodel6 <- df.with.infriskGH2DT.6 %>%
  dplyr::mutate(risk.temp.lim = ifelse(P0_Temp>10,ifelse(P0_Temp>35,0,infsevrisk),0)) %>% #temp limits
  dplyr::mutate(risk.limited6 = ifelse(risk.temp.lim > 0, ifelse(risk.temp.lim>1,1,risk.temp.lim),0)) %>% #risk limits 
  mutate(ripe.rot.risk = risk.limited6*current.pheno.risk)




#### GH2-DT Summary Table
maxinfriskbydayGHjoined <- df.with.infriskGHmodel6 %>%
  mutate(Month.day = substr(DateTime,6,10)) %>% #just kept the month&day for less redundant values (need to be careful if data is multiple years worth)
  group_by(Month.day) %>% #each day is grouped because it will be summarised over each day
  summarise(GH2DT = max(ripe.rot.risk), Wetnesss_duration_hr2 = max(cumsum_Wethr), Mean_P0_Temperature_C = mean(P0_Temp)) %>% #desired maximum infection risk of the day, the real world wetness duration and average P0_Temperature#only keep the top 10 infection days
  arrange(Month.day) %>% #sort chronologically rather than by the highest infection risk
  kbl(caption = "**Table 1.** Days with the highest infection risk (GREENHOUSE 2 DT model) during the queried period", align = "c") %>%
  kable_classic(full_width = F, html_font = "Cambria")




#simple chart with all data

xmin1=min(df.with.infriskGHmodel6$DateTime) #limits for color shading
xmax1 = max(df.with.infriskGHmodel6$DateTime) # limits for color shading

# base plot
infection.chart.base <- ggplot(NULL, aes(NULL)) +
  labs(x = "Date", y = "Probability of infection") + #label axes
  theme_classic() +
  scale_y_continuous(name="Disease risk", breaks = seq(0,1,0.1), expand = c(0,0)) + #set a custom scale to more easily see the data with tick marks every 0.1
  geom_hline(yintercept=0.65, color = "orange", size = 0.5)+
  geom_hline(yintercept=0.65, color = "red", size = .5)+
  annotate("rect", xmin = xmin1, xmax = xmax1, ymin = 0, ymax = 0.65, alpha = 0.1, fill = "green") +
  annotate("rect", xmin = xmin1, xmax = xmax1, ymin = 0.65, ymax = 0.65, alpha = 0.1, fill = "orange") +
  annotate("rect", xmin = xmin1, xmax = xmax1, ymin = 0.65, ymax = 1, alpha = 0.1, fill = "red") +
  scale_x_datetime(date_breaks = "1 day", expand = c(0,0), date_labels = "%b %d")

#####GH2DT chart
all.data.figureGH.sensors <- infection.chart.base +
  geom_line(data = df.with.infriskGHmodel6, aes(DateTime,ripe.rot.risk), size = 1, alpha = 0.8) +
  ggtitle("GH2DT")

all.data.figureGH.sensors






########### DF2-NN model function


df.with.infriskDF2NN <- fresh.sensors6 %>% 
  group_by(Group_TFchange) %>%  
  mutate(cumsum_Wethr = cumsum(Wet_hr)) %>% #calculate cumulative sum of Wetness hours in each group
  mutate(Wetmax24 = ifelse(cumsum_Wethr>24,24,cumsum_Wethr)) %>% #model was only created with max of 24 hours
  dplyr::mutate(drygaplength = ifelse(is_wet == TRUE, 0, Wet_TFchangelength)) %>%
  dplyr::mutate(LWD = ifelse(drygaplength > 48, Wet_hr, Wetmax24))


#actual model
risk_DF2NN_fun <- function(LWD, P0_Temp) {
  H1_3_1 <- tanh(-70.17431845 + 1.8578573232*P0_Temp + 1.8506739621*LWD)
  
  H1_2_1 <- tanh(18.8373120639994 + -0.00388711424490551*P0_Temp + -0.739676546562124*LWD)
  
  H1_1_1 <-  tanh(-0.686972746892936 + 0.736364917645139*P0_Temp + -1.33325598926744*LWD)
  
  infsevrisk <- (1/(1 + exp(4.12545439783833 + 0.436894051174693 * H1_1_1 + -3.64951612080918 * H1_2_1 + -0.608296811922651 * H1_3_1)))
  return(infsevrisk)
}
  
#First create a new df with blank column. Then loop the function over the column            
df.with.infriskDF2NNb <- df.with.infriskDF2NN  %>%
  mutate(infsevrisk = 0)

for(i in 1:nrow(df.with.infriskDF2NNb)){
  WD <- df.with.infriskDF2NNb$LWD[i]
  Temp <- df.with.infriskDF2NNb$P0_Temp[i]
  n = risk_DF2NN_fun(WD,Temp)
  df.with.infriskDF2NNb[i,"infsevrisk"] <- n
}

############ Calculating of final ripe rot risk determined from the DF2-NN model, which is the DF2-NN risk multiplied by the risk due to the current phenological stage
df.with.infriskDF2NN <- df.with.infriskDF2NNb %>%
  dplyr::mutate(risk.temp.lim = ifelse(P0_Temp>10,ifelse(P0_Temp>35,0,infsevrisk),0)) %>% #temp limits
  dplyr::mutate(risk.limited6 = ifelse(risk.temp.lim > 0, ifelse(risk.temp.lim>1,1,risk.temp.lim),0))  %>% #risk limits 
  mutate(ripe.rot.risk = risk.limited6*current.pheno.risk)


##DF2NN Table
maxinfriskbydayDF2NN <- df.with.infriskDF2NN %>%
  mutate(Month.day = substr(DateTime,6,10)) %>% #just kept the month&day for less redundant values (need to be careful if data is multiple years worth)
  group_by(Month.day) %>% #each day is grouped because it will be summarised over each day
  summarise(DF2NN = max(ripe.rot.risk), Wetnesss_duration_hr2 = max(cumsum_Wethr), Mean_P0_Temperature_C = mean(P0_Temp)) %>% #desired maximum infection risk of the day, the real world wetness duration and average P0_Temperature#only keep the top 10 infection days
  arrange(Month.day) %>% #sort chronologically rather than by the highest infection risk
  kbl(caption = "**Table 2.** Days with the highest infection risk (DETACHED FRUIT 2 NN model) during the queried period", align = "c") %>%
  kable_classic(full_width = F, html_font = "Cambria")





##DF2NN chart
# base plot
infection.chart.base.DF2NN <- ggplot(NULL, aes(NULL)) +
  labs(x = "Date", y = "Probability of infection") + #label axes
  theme_classic() +
  scale_y_continuous(name="Disease risk", breaks = seq(0,1,0.1), expand = c(0,0)) + #set a custom scale to more easily see the data with tick marks every 0.1
  geom_hline(yintercept=0.45, color = "orange", size = 0.5)+
  geom_hline(yintercept=0.45, color = "red", size = .5)+
  annotate("rect", xmin = xmin1, xmax = xmax1, ymin = 0, ymax = 0.45, alpha = 0.1, fill = "green") +
  annotate("rect", xmin = xmin1, xmax = xmax1, ymin = 0.45, ymax = 0.45, alpha = 0.1, fill = "orange") +
  annotate("rect", xmin = xmin1, xmax = xmax1, ymin = 0.45, ymax = 1, alpha = 0.1, fill = "red") +
  scale_x_datetime(date_breaks = "1 day", expand = c(0,0), date_labels = "%b %d")

all.data.DF2NNchart <- infection.chart.base.DF2NN + geom_line(data = df.with.infriskDF2NN, aes(DateTime,ripe.rot.risk), size = 1, alpha = 0.8)+
  ggtitle("DF2NN WRV")

all.data.DF2NNchart


df.with.infriskDF2NNwGH <- df.with.infriskDF2NN
df.with.infriskDF2NNwGH$GH2DTrisk <- df.with.infriskGHmodel6$ripe.rot.risk

##Combined Table
maxinfriskbydaycomb <- df.with.infriskDF2NNwGH %>%
  mutate(Month.day = substr(DateTime,6,10)) %>% #just kept the month&day for less redundant values (need to be careful if data is multiple years worth)
  group_by(Month.day) %>% #each day is grouped because it will be summarised over each day
  summarise(GH2DT = max(GH2DTrisk), DF2NN = max(ripe.rot.risk), LWD_hr = max(cumsum_Wethr), Mean_Temperature_C = mean(P0_Temp)) %>% #desired maximum infection risk of the day, the real world wetness duration and average P0_Temperature#only keep the top 10 infection days
  arrange(Month.day) %>% #sort chronologically rather than by the highest infection risk
  mutate(isGH2DThi = ifelse(GH2DT>=0.65,"YES","no")) %>%
  mutate(isDF2NNhi = ifelse(DF2NN>=0.45,"YES","no")) %>%
  kbl(caption = "Table. Daily highest infection risk value, total hours of leaf wetness duration, mean temperature,and 'YES' if risk crossed the threshold for either model", align = "c", digits = c(3,3,3,1,1,3,3)) %>%
  kable_classic(full_width = F, html_font = "Cambria")


maxinfriskbydaycomb








