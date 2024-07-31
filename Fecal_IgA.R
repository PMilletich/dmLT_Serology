#' Patricia Milletich 
#' July, 2024

library(ggplot2)
library(ggpubr)
library(DescTools)
library(EnvStats)

OG_data = read.csv("EDU_1.csv")
#Fold Change as.Date(as.character(survey$tx_start), format="%Y/%m/%d")
OG_data$STARTDT_1 = as.Date(OG_data$STARTDT, format="%m/%d/%Y")
OG_data$COL_DATE_1 = as.Date(OG_data$COL_DATE, format="%m/%d/%Y")
OG_data$Date_Change = OG_data$COL_DATE_1 - OG_data$STARTDT_1
OG_data$Date_Change_2 = as.numeric(OG_data$Date_Change)

#Participant Excluded because ineligible at baseline: "08ECI088"
#Participant Excluded because no post-baseline specimens: "08ECI040"
#13202 --> 13119
OG_data = subset(OG_data, ! OG_data$PATID %in% c("08ECI088", "08ECI040"))
OG_data$TRTTRUE = ifelse(OG_data$PATID %in% c("08ECI131","08ECI133"), #4 --> 5B
                         "Two Doses of Intradermal dmLT 2.0 mcg", 
                         ifelse(OG_data$PATID %in% c("08ECI108", "08ECI136", "08ECI158", "08ECI179"), #4 --> 5A
                                "One Dose of Intradermal dmLT 2.0 mcg", 
                                ifelse(OG_data$PATID %in% c("08ECI219"), #5C --> 5B
                                       "Two Doses of Intradermal dmLT 2.0 mcg", OG_data$TRTTRUE)))

#Second Vaccination not received; exclude Day 23 – Day 223 in PP Population:
No_Second = subset(OG_data, OG_data$PATID %in% c("08ECI004","08ECI041","08ECI077","08ECI075","08ECI105"))
No_Second = subset(No_Second, No_Second$Date_Change_2 >= 22)
OG_data = subset(OG_data, ! OG_data$SN %in% No_Second$SN)


#Third Vaccination not received; exclude Day 44 – Day 223 in PP Population:
No_Third = subset(OG_data, OG_data$PATID %in% c("08ECI051","08ECI098","08ECI076"))
No_Third = subset(No_Third, No_Third$Date_Change_2 >= 42)
OG_data = subset(OG_data, ! OG_data$SN %in% No_Third$SN)
remove(No_Third); remove(No_Second)

#Rename the Treatment groups 
OG_data$TRTTRUE = gsub("One Dose of Intradermal dmLT 2.0 mcg", "dmLT 2.0 mcg: One Dose", OG_data$TRTTRUE)
OG_data$TRTTRUE = gsub("Two Doses of Intradermal dmLT 2.0 mcg", "dmLT 2.0 mcg: Two Doses", OG_data$TRTTRUE)
OG_data$TRTTRUE = gsub("Three Doses of Intradermal dmLT 2.0 mcg", "dmLT 2.0 mcg", OG_data$TRTTRUE)
OG_data$TRTTRUE = gsub("Intradermal ", "", OG_data$TRTTRUE)

OG_data$Group = ifelse(grepl("Placebo", OG_data$TRTTRUE), "Placebo", OG_data$TRTTRUE)

#No baseline data 
OG_data = subset(OG_data, OG_data$PATID != "08ECI129")

Antibody_subset = subset(OG_data, OG_data$TESTTYPA %in% c("DATI", "TATI")) 
Antibody_subset = subset(Antibody_subset, Antibody_subset$ST == "FRZ")


gm = function(x){
  gm1 = exp(mean(log(x)))
  cil = exp(gm1-(1.96*(sd(log(x), na.rm = T)/sqrt(length(x)))))
  ciupp = exp(gm1+(1.96*(sd(log(x), na.rm = T)/sqrt(length(x)))))
  vec = c(round(gm1,2), round(cil,2), round(ciupp,2))
  return (vec)
}

Final_data = data.frame()
current_Group = "dmLT 0.1 mcg"
for (current_Group in unique(Antibody_subset$Group)) {
  print(current_Group)
  group_subset = subset(Antibody_subset, Antibody_subset$Group == current_Group)
  
  #Pool to the correct days 
  group_subset$Study_Day = ifelse(group_subset$Date_Change_2 %in% 0:1, 1, 
                                  ifelse(group_subset$Date_Change_2 %in% 7:9, 8, 
                                         ifelse(group_subset$Date_Change_2 %in% 21:23, 22, 
                                                ifelse(group_subset$Date_Change_2 %in% 28:30, 29, 
                                                       ifelse(group_subset$Date_Change_2 %in% 42:44, 43, 
                                                              ifelse(group_subset$Date_Change_2 %in% 49:51, 50, 
                                                                     ifelse(group_subset$Date_Change_2 %in% 63:65, 64, 
                                                                            ifelse(group_subset$Date_Change_2 %in% 69:73, 71, 
                                                                                   ifelse(group_subset$Date_Change_2 < 0, -1, "Bad Day")))))))))
  
  
  group_subset = subset(group_subset, group_subset$Study_Day != "Bad Day")

  Calculated = data.frame()
  for (current_pat in unique(group_subset$PATID)) {
    current_pat_subset = subset(group_subset, group_subset$PATID == current_pat)

    for (current_time in unique(current_pat_subset$Study_Day)) {
      current_pat_subset_time = subset(current_pat_subset, current_pat_subset$Study_Day == current_time)
      if (current_time != 1) {
        current_pat_subset$ERESULTA = ifelse(current_pat_subset$ERESULTA == 0.0015, 
                                                  0, current_pat_subset$ERESULTA)
      }
      Total = current_pat_subset_time[current_pat_subset_time$TESTTYPA == "TATI", "ERESULTA"]
      current_pat_subset[current_pat_subset$Study_Day == current_time,"Specific.Total"] = 
        current_pat_subset[current_pat_subset$Study_Day == current_time,"ERESULTA"]/Total
    }
    
    current_pat_subset_Specific = subset(current_pat_subset, current_pat_subset$TESTTYPA %in% c("DATI", "DGTI"))
    
    Baseline = subset(current_pat_subset_Specific, current_pat_subset_Specific$Study_Day == 1)
    if (nrow(Baseline) == 0) {
      print(paste("No Day 1:", current_pat))
      Baseline = subset(current_pat_subset_Specific, current_pat_subset_Specific$Study_Day == -1)
    } else if (nrow(Baseline) > 1) {
      print(paste("Too many baselines:", current_pat))
      Baseline = Baseline[1,]
    }
    Baseline = Baseline$Specific.Total

    current_pat_subset_Specific$Fold.Rise = current_pat_subset_Specific$Specific.Total/Baseline
    Calculated = rbind(Calculated, current_pat_subset_Specific)
  }  
  
  Calculated = subset(Calculated, Calculated$Study_Day > 1)
  Calculated_1 = subset(Calculated, Calculated$ERESULTA == 0.0015)
  
  line_data = data.frame()
  for (current_time in unique(Calculated$Study_Day)) {
    time_subset = subset(Calculated, Calculated$Study_Day == current_time)
    if (length(time_subset$Fold.Rise[time_subset$Fold.Rise>0]) == 1) {
      Fold.Rise_GM = time_subset$Fold.Rise[time_subset$Fold.Rise>0]
      Fold.Rise_Lower = time_subset$Fold.Rise[time_subset$Fold.Rise>0]
      Fold.Rise_Upper = time_subset$Fold.Rise[time_subset$Fold.Rise>0]
    } else {
      Fold.Rise_GM = geoMean(time_subset$Fold.Rise[time_subset$Fold.Rise>0])
      Fold.Rise_SD = geoSD(time_subset$Fold.Rise[time_subset$Fold.Rise>0], 
                           na.rm = T, sqrt.unbiased = TRUE)
      MarginofError = (Fold.Rise_SD-1)/sqrt(length(time_subset$Fold.Rise[time_subset$Fold.Rise>0]))
      Fold.Rise_Lower = Fold.Rise_GM - MarginofError
      Fold.Rise_Upper = Fold.Rise_GM + MarginofError
    }

    Seroconverted = nrow(subset(time_subset, time_subset$Fold.Rise >= 4))/nrow(time_subset)
    Seroconverted = round(Seroconverted, 3)*100
    
    time_row = data.frame("Study_Day" = current_time,
                          "Seroconversion" = Seroconverted,
                          "Fold.Rise_GM" = Fold.Rise_GM,
                          "Fold.Rise_Lower" = Fold.Rise_Lower,
                          "Fold.Rise_Upper" = Fold.Rise_Upper,
                          "Group" = paste(current_Group, " (N=", length(unique(Calculated$PATID)), ")", sep = ""))
    line_data = rbind(line_data, time_row)
    
  }

  line_data$Study_Day = as.numeric(line_data$Study_Day)
  
  Final_data = rbind(Final_data, line_data)
}
#Final_data = data.frame()
Final_data = unique(Final_data)
#Scale the FR by the maximum of the treatments to ignore Placebo error bars 
Scale_value = max(c(Final_data$Fold.Rise_GM, 
                    Final_data$Fold.Rise_Lower,
                    Final_data$Fold.Rise_Upper), na.rm = T)
Final_data$Mean_FR_scaled <- (Final_data$Fold.Rise_GM/Scale_value) * 100
Final_data$Upper_FR_scaled <- (Final_data$Fold.Rise_Upper/Scale_value) * 100
Final_data$Lower_FR_scaled <- (Final_data$Fold.Rise_Lower/Scale_value) * 100

Scale_value = max(c(Final_data$Fold.Rise_GM, 
                    Final_data$Fold.Rise_Lower,
                    Final_data$Fold.Rise_Upper), na.rm = T)
Final_data$Mean_FR_scaled <- (Final_data$Fold.Rise_GM/Scale_value) * 100
Final_data$Upper_FR_scaled <- (Final_data$Fold.Rise_Upper/Scale_value) * 100
Final_data$Lower_FR_scaled <- (Final_data$Fold.Rise_Lower/Scale_value) * 100

Final_data = subset(Final_data, is.na(Final_data$Fold.Rise_GM) == F)

Final_data$Lower_FR_scaled = ifelse(Final_data$Lower_FR_scaled < 0, 0, Final_data$Lower_FR_scaled)
  
Fold.Rise = ggplot() +   
  theme_bw() + 
  facet_wrap(~Group, ncol = 1) + 
  geom_vline(xintercept = 1, color = "lightcoral", linewidth = 1) + 
  geom_vline(data=subset(Final_data, grepl("One Dose", Final_data$Group) == F),
             aes(xintercept=22), colour="lightcoral", linewidth = 1) +
  geom_vline(data=subset(Final_data, grepl("One Dose", Final_data$Group) == F &
                           grepl("Two Dose", Final_data$Group) == F),
             aes(xintercept=43), colour="lightcoral", linewidth = 1) +
  
  geom_bar(data = Final_data, aes(x = Study_Day, y  = Seroconversion), stat = "identity", 
           alpha= 0.5, fill = "grey", color = "grey24", width = 4) + 
  
  geom_hline(yintercept = 0, color = "black")+  
  scale_y_continuous(
    position = "right",
    name = "Percent Seroconverted",
    sec.axis = sec_axis(~./100 * Scale_value, 
                        name= expression(paste("Geometric Mean Titer (95% CI)")))) + 
  geom_line(data = Final_data, aes(x = Study_Day, y = Mean_FR_scaled)) + 
  geom_point(data = Final_data, aes(x = Study_Day, y = Mean_FR_scaled)) +  
  geom_linerange(data = Final_data, aes(x = Study_Day, ymin = Lower_FR_scaled,
                                       ymax = Upper_FR_scaled)) +
  ggtitle("Fecal IgA [dmLT/Total]") + 
  theme(strip.text = element_text(size = 6),
        axis.text.y = element_text(size= 7.5)); Fold.Rise


jpeg("./Figures/Fecal_IgA.jpeg", res = 400, height = 3000, width = 2000)
Fold.Rise
dev.off()
