# Author: Trish Milletich, PhD

library(ggeffects)  # install the package first if you haven't already, then load it
library(ggplot2)
library(ggpubr)
library(geepack)
#library(lmerTest)

OG_data = read.csv("DataAvailability_Biomarkers.csv")
OG_data = subset(OG_data, OG_data$ST != "FRZ")
Test_subset = subset(OG_data, OG_data$Date_Change_2 > 1)

#Exploratory 
table(sort(OG_data$Date_Change_2))
table(sort(OG_data$Group))
Dose_1 = subset(Test_subset, Test_subset$Date_Change_2 %in% c(8:23))
Dose_2 = subset(Test_subset, Test_subset$Date_Change_2 %in% c(29:45))
Dose_3 = subset(Test_subset, Test_subset$Date_Change_2 %in% c(50:232))

#Example values 
current_patient = unique(Test_subset$PATID)[1]
current_patient = "08ECI003"

Max_data = data.frame()
i = 1
for (current_patient in unique(Test_subset$PATID)) {
  patient_subset = subset(Test_subset, Test_subset$PATID == current_patient)
  if (patient_subset$Group[1] == "Cohort 5A:One Dose dmLT 2.0 mcg") {
    Time_list = list(c("After Dose 1", 8,232))
  } else if (patient_subset$Group[1] == "Cohort 5B:Two Doses dmLT 2.0 mcg") {
    Time_list = list(c("After Dose 1", 8,23),
                     c("After Dose 2", 29,232))
  } else {
    Time_list = list(c("After Dose 1", 8, 23),
                     c("After Dose 2", 29,45),
                     c("After Dose 3", 50,232))
    }
  
  
  if (patient_subset$Group[1] !=  "Placebo") {
    
    
  #After Dose 2
  for (current_time_list in Time_list) {
    
    Post_dose = subset(patient_subset, patient_subset$Date_Change_2 %in% 
                         c(as.numeric(current_time_list[2]):as.numeric(current_time_list[3])))
    
    current_time = current_time_list[1]

    if (nrow(Post_dose) > 0) {
      #ALS Present on only certain timepoints 
      if (nrow(Post_dose[Post_dose$ST == "PBS" & 
                         Post_dose$TESTTYPA == "DATI",]) > 0) {
        ALS_IgG = Post_dose[Post_dose$ST == "PBS" & 
                              Post_dose$TESTTYPA == "DGTI","ERESULTA"]
        ALS_IgA = Post_dose[Post_dose$ST == "PBS" & 
                              Post_dose$TESTTYPA == "DATI","ERESULTA"]
      } else {
        ALS_IgA = NA; ALS_IgG = NA
      }
      
      
      SER_IgA = max(Post_dose[Post_dose$ST == "SER"& 
                                Post_dose$TESTTYPA == "DATI","ERESULTA"], na.rm = T)
      SER_IgG = max(Post_dose[Post_dose$ST == "SER"& 
                                Post_dose$TESTTYPA == "DGTI","ERESULTA"], na.rm = T)
      
      
      LTTN = max(Post_dose[Post_dose$TESTTYPA == "LTTN","ERESULTA"], na.rm = T)
      
      data = data.frame("TimePoint" = current_time, 
                        "ALS_IgA" = log10(ALS_IgA), 
                        "SER_IgA" = log10(SER_IgA),
                        "ALS_IgG" = log10(ALS_IgG), 
                        "SER_IgG" = log10(SER_IgG), 
                        "LTTN" = log10(LTTN),
                        "Cohort" = patient_subset[1,"Group"], 
                        "PATID" = current_patient,
                        "PATID_Num" = i)
      Max_data = rbind(Max_data, data)
    }
  }
  
  }
  i = i + 1 
} # ID loop

Max_data$Cohort = gsub("mcg", "µg", Max_data$Cohort)
current_pair= c("ALS_IgG", "SER_IgG")
Total_data = data.frame()
for (current_pair in list(c("ALS_IgA", "SER_IgA"), c("LTTN", "SER_IgA"),c("ALS_IgA", "LTTN"),
                          c("ALS_IgG", "SER_IgG"),c("LTTN", "SER_IgG"),c("ALS_IgG", "LTTN"))) {
  
  print(current_pair)
  #Create new data with set columns for each pair 
  new_data = Max_data
  new_data$Pair1 = new_data[,current_pair[1]] 
  new_data$Pair2 = new_data[,current_pair[2]] 
  
  new_data = subset(new_data, is.na(new_data$Pair1) == F & 
                      is.na(new_data$Pair2) == F)
  
  print(length(unique(new_data$PATID)))
  print(nrow(new_data))
  
  #https://www.rdocumentation.org/packages/geepack/versions/1.3.11.1/topics/geeglm
  mixed.lmer <- geeglm(Pair1 ~ Pair2, id=PATID_Num,
                       data = new_data, family = gaussian(), corstr = "exchangeable")
  
  #Summarize model
  sum_model = summary(mixed.lmer)
  
  #append summary data 
  print_df = sum_model$coefficients
  print_df$P = ifelse(print_df$`Pr(>|W|)` < 0.001, "<0.001", round(print_df$`Pr(>|W|)`, 3))
  print_df[,c("Estimate","Std.err","Wald", "Pr(>|W|)")] = round(print_df[,c("Estimate","Std.err","Wald", "Pr(>|W|)")], 2)
  print_df$Pair = paste(current_pair,collapse = ":")
  Total_data = rbind(Total_data, print_df)
  
  #Extract summary values from model 
  Coeff = sum_model$coefficients[,"Estimate"][2]
  StdError = sum_model$coefficients[,"Std.err"][2]
  P = sum_model$coefficient[,"Pr(>|W|)"][2]
  P = ifelse(P < 0.001, "<0.001", paste("=", round(P, 3), sep = ""))
  
  # Extract the prediction data frame
  pred.mm <- ggemmeans(mixed.lmer, terms = c("Pair2"))  
  
  #Rename from original variables 
  title_2 = ifelse(current_pair[2] == "LTTN", "LT Neutralizing", 
                   gsub("SER", "Serum", current_pair[2]))
  title_1 = ifelse(current_pair[1] == "LTTN", "LT Neutralizing", 
                   gsub("SER", "Serum", current_pair[1]))
  
  pred.mm.df = data.frame(pred.mm)
  
  # Plot the predictions
  current_plot = ggplot() +
    geom_line(data = pred.mm, aes(x = x, y = predicted)) +
    geom_ribbon(data = pred.mm.df,
                aes(x = x, ymin = conf.low, ymax = conf.high),
                fill = "grey", alpha = 0.5) +
    
    geom_point(data = new_data,                    
               aes(x = Pair2, y = Pair1, shape = TimePoint), alpha = 0.75) +
    labs(x = paste("log10(", gsub("_", " ", title_2), ")", sep = ""),
         y = paste("log10(", gsub("_", " ", title_1), ")", sep = ""),
         title = paste(gsub("_", " ", title_1), gsub("_", " ", title_2), sep = " vs. "))  +
    annotate("text",
             label = paste("β: ", format(round(Coeff, 2), nsmall=2),
                           "\nSE: ", round(StdError,3),
                           "\np",P, sep = ""),
             x = min(new_data$Pair2, na.rm = T)+ 0.35,
             y = max(new_data$Pair1, na.rm = T) - (max(new_data$Pair1, na.rm = T)/6) , size =4) +
    theme_bw() + theme(legend.text=element_text(size=12)); current_plot
  
  assign(paste(current_pair[1], current_pair[2], "OG", sep = "_"),  current_plot)
}


jpeg("./Figure_6.jpeg", 
     res = 500, height = 3000, width = 5000)
ggarrange(ALS_IgA_SER_IgA_OG, LTTN_SER_IgA_OG, ALS_IgA_LTTN_OG, 
          ALS_IgG_SER_IgG_OG, LTTN_SER_IgG_OG, ALS_IgG_LTTN_OG,
          ncol = 3, nrow = 2, 
          common.legend = T)
dev.off()
