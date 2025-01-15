# Author: Trish Milletich, PhD

library(ggplot2)
library(ggpubr)
library(DescTools)

#Preset example Values 
current_Test = "DGTI"
current_specimen = "SER"
current_Group = "Cohort 1:dmLT 0.1 µg"

OG_data = read.csv("DataAvailability_Biomarkers.csv")
OG_data = subset(OG_data, OG_data$PATID != "08ECI129")

All_data_groups = data.frame()
Antibody_subset = subset(OG_data, OG_data$TESTTYPA %in% c("DATI", "TATI")) 
Antibody_subset = subset(Antibody_subset, Antibody_subset$ST == "FRZ")

ID_Count = Antibody_subset[,c("Group", "PATID")]
#ID_Count = ID_Count[ID_Count$Group == "Cohort 1:dmLT 0.1 µg",]
ID_Count = unique(ID_Count)
table(ID_Count$Group)

Final_data = data.frame()
current_Group = "Cohort 4+5C: dmLT 2.0 mcg"
for (current_Group in unique(Antibody_subset$Group)) {
  #print(current_Group)
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
  
  group_subset_bad = subset(group_subset, group_subset$Study_Day == "Bad Day")
  
  group_subset = subset(group_subset, group_subset$Study_Day != "Bad Day")
  group_subset$Study_Day= as.numeric(group_subset$Study_Day)
  
  Calculated = data.frame()
  no_register = c()
  
  for (current_pat in unique(group_subset$PATID)) {
    current_pat_subset = subset(group_subset, group_subset$PATID == current_pat)
    current_time =  22
    for (current_time in unique(current_pat_subset$Study_Day)) {
      current_pat_subset_time = subset(current_pat_subset, current_pat_subset$Study_Day == as.numeric(current_time))
      
      if (nrow(current_pat_subset_time) > 2 & length(unique(current_pat_subset_time$Date_Change_2)) > 1) {
        current_pat_subset_time = subset(current_pat_subset_time, current_pat_subset_time$Date_Change_2 == 
                                           max(current_pat_subset_time$Date_Change_2))
      } else if (nrow(current_pat_subset_time) > 2) {
        current_pat_subset_time = subset(current_pat_subset_time, current_pat_subset_time$SN == current_pat_subset_time$SN[1])
      }
      
      current_pat_subset_time = subset(current_pat_subset, current_pat_subset$Study_Day == current_time)
      if (current_time != 1) {
        current_pat_subset_time$ERESULTA = ifelse(current_pat_subset_time$ERESULTA == 0.0015, 
                                             NA, current_pat_subset_time$ERESULTA)
      }
      DATI = current_pat_subset_time[current_pat_subset_time$TESTTYPA == "DATI", "ERESULTA"]
      
      Total = current_pat_subset_time[current_pat_subset_time$TESTTYPA == "TATI", "ERESULTA"]
      current_pat_subset[current_pat_subset$Study_Day == current_time,"Specific.Total"] =  DATI/Total
      
    }
    current_pat_subset_Specific = subset(current_pat_subset, current_pat_subset$TESTTYPA %in% c("DATI", "DGTI"))
    
    Baseline = subset(current_pat_subset_Specific, current_pat_subset_Specific$Study_Day == 1)
    if (nrow(Baseline) == 0) {
      Baseline = subset(current_pat_subset_Specific, current_pat_subset_Specific$Study_Day == -1)
    } else if (nrow(Baseline) > 1) {
      Baseline = Baseline[1,]
    }
    Baseline = Baseline$Specific.Total
    
    current_pat_subset_Specific$Fold.Rise = current_pat_subset_Specific$Specific.Total/Baseline
    Calculated = rbind(Calculated, current_pat_subset_Specific)
    
  }  
  
  
  Calculated = subset(Calculated, Calculated$Study_Day > 0)
  Calculated = Calculated[is.na(Calculated$Specific.Total)== F,]
  Calculated = subset(Calculated, Calculated$Study_Day > 1)

  if (current_Group == "Cohort 1:dmLT 0.1 mcg") {
    N = 12
    Update = "None"
  } else if (current_Group == "Cohort 2:dmLT 0.3 mcg") {
    N = 10
    Update = "None"
  } else if (current_Group == "Cohort 3:dmLT 1.0 mcg") {
    N = 12
    Update = "None"
  } else if (current_Group == "Cohort 4+5C: dmLT 2.0 mcg") {
    N = 18
    Update= data.frame("Study_Day" = c(8, 22, 29, 43, 50, 64, 71),
                       "OG_Sample_Num" = c(11, 13, 13, 12, 10, 12, 11))
  } else if (current_Group == "Cohort 5A:One Dose dmLT 2.0 mcg") {
    N = 15
    Update= data.frame("Study_Day" = c(8, 22, 29),
                       "OG_Sample_Num" = c(9,8,8))
  } else if (current_Group == "Cohort 5B:Two Doses dmLT 2.0 mcg") {
    N = 15
    Update= "None"
  } else if (current_Group == "Placebo") {
    N = 16
    Update= data.frame("Study_Day" = c(8, 22, 29, 43, 50, 64, 71),
                       "OG_Sample_Num" = c(8, 7, 7, 7, 6, 7, 5))
  } else {
    print("Unknown Cohort") 
    break
  }
    print(paste(current_Group, N, length(unique(Calculated$PATID))))

  
  
  
  line_data = data.frame()
  for (current_time in unique(Calculated$Study_Day)) {
    time_subset = subset(Calculated, Calculated$Study_Day == current_time)
    if (typeof(Update) != "list") {
      time_subset_n = nrow(time_subset)
    } else {
      time_subset_n = subset(Update, Update$Study_Day == current_time)$OG_Sample_Num
    }
    
    
    if (length(time_subset$Fold.Rise[time_subset$Fold.Rise>0]) == 1) {
      Fold.Rise_GM = time_subset$Fold.Rise[time_subset$Fold.Rise>0]
      Fold.Rise_Lower = time_subset$Fold.Rise[time_subset$Fold.Rise>0]
      Fold.Rise_Upper = time_subset$Fold.Rise[time_subset$Fold.Rise>0]
      
    } else {
      Geom_mean_FR = suppressWarnings(Gmean(time_subset$Fold.Rise[time_subset$Fold.Rise>0],
                                            method = c("classic"), conf.level = 0.95,
                                            sides = c("two.sided"), na.rm = T))
      
      Fold.Rise_GM = Geom_mean_FR[1]
      Fold.Rise_Lower = Geom_mean_FR[2]
      Fold.Rise_Upper = Geom_mean_FR[3]
    }
    
    Seroconverted = nrow(subset(time_subset, time_subset$Fold.Rise >= 4))/
      time_subset_n
    Seroconverted = round(Seroconverted, 3)*100
    
    time_row = data.frame("Study_Day" = current_time,
                          "Fold.Rise_GM" = Fold.Rise_GM,
                          "Fold.Rise_Lower" = Fold.Rise_Lower,
                          "Fold.Rise_Upper" = Fold.Rise_Upper,
                          "Seroconversion" = Seroconverted,
                          "Fraction" = paste(nrow(subset(time_subset, time_subset$Fold.Rise >= 4)),
                                             time_subset_n, sep = '_'),
                          "Group" = paste(current_Group, " (N=", N, ")", sep = ""))
    line_data = rbind(line_data, time_row)
    
  }
  
  line_data$Study_Day = as.numeric(line_data$Study_Day)
  
  Final_data = rbind(Final_data, line_data)
}
#Final_data = data.frame()
Final_data = unique(Final_data)
#Scale the FR by the maximum of the treatments to ignore Placebo error bars 
Scale_value = max(c(Final_data$Fold.Rise_GM))*5
#Scale_value = Scale_value[2]
Final_data$Mean_FR_scaled <- (Final_data$Fold.Rise_GM/Scale_value) * 100
Final_data$Upper_FR_scaled <- (Final_data$Fold.Rise_Upper/Scale_value) * 100
Final_data$Lower_FR_scaled <- (Final_data$Fold.Rise_Lower/Scale_value) * 100

Final_data = subset(Final_data, is.na(Final_data$Fold.Rise_GM) == F)

#Final_data$Lower_FR_scaled = ifelse(Final_data$Lower_FR_scaled < 0, 0, Final_data$Lower_FR_scaled)

Final_data$Upper_FR_scaled = ifelse(Final_data$Upper_FR_scaled > 100, 100, Final_data$Upper_FR_scaled)
Final_data$Fraction = ifelse(Final_data$Study_Day == "1", "", Final_data$Fraction)
Final_data$Fraction = gsub("_", "/", Final_data$Fraction)

Final_data_1 = Final_data 
Final_data_1$Group = gsub("mcg", "µg", Final_data_1$Group)

Fold.Rise = ggplot() +   
  theme_bw() +
  facet_wrap(~Group, ncol = 1) + 
  geom_vline(xintercept = 1, color = "lightcoral", linewidth = 1) + 
  geom_vline(data=subset(Final_data_1, grepl("One Dose", Final_data_1$Group) == F),
             aes(xintercept=22), colour="lightcoral", linewidth = 1) +
  geom_vline(data=subset(Final_data_1, grepl("One Dose", Final_data_1$Group) == F &
                           grepl("Two Dose", Final_data_1$Group) == F),
             aes(xintercept=43), colour="lightcoral", linewidth = 1) +
  
  geom_bar(data = Final_data_1, aes(x = Study_Day, y  = Seroconversion), stat = "identity", 
           alpha= 0.5, fill = "grey", color = "grey24", width = 4) + 
  
  geom_hline(yintercept = 0, color = "black")+  
  scale_y_continuous(
    position = "right",
    name = "Percent Responders",
    sec.axis = sec_axis(~./100 * Scale_value, 
                        name= expression(paste("Geometric Mean FR (95% CI)")))) + 
  geom_text(data = Final_data_1, aes(x = Study_Day, y = Seroconversion+25, label = Fraction),
            size = 2) + 
  geom_line(data = Final_data_1, aes(x = Study_Day, y = Mean_FR_scaled)) + 
  geom_linerange(data = Final_data_1, aes(x = Study_Day, ymin = Lower_FR_scaled,
                                        ymax = Upper_FR_scaled), color = "grey40") +
  geom_point(data = Final_data_1, aes(x = Study_Day, y = Mean_FR_scaled)) +  
  ggtitle("Fecal IgA [dmLT/Total]") + xlab("Study Day") +
  theme(strip.text = element_text(size = 6),
        axis.text.y = element_text(size= 7.5),          
        strip.text.x = element_text(margin = margin(0.5,0,0.5,0, "mm"))); Fold.Rise


jpeg("./Figures/Final_Figures/Supplemental_Figure_2.jpeg", res = 400, height = 2000, width = 2000)
Fold.Rise
dev.off()

X = subset(Final_data, Final_data$Group == "Cohort 3:dmLT 1.0 µg (N=11)")

