#' Patricia Milletich 
#' July, 2024

library(ggplot2)
library(ggpubr)
library(DescTools)
library(EnvStats)
#Geometric CI: https://stats.stackexchange.com/questions/285173/how-to-calculate-confidence-interval-for-a-geometric-mean


OG_data = read.csv("EDU_1.csv")
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

OG_data_Count = OG_data[,c("PATID", "Group")]
OG_data_Count = unique(OG_data_Count)
table(OG_data_Count$Group)


#No baseline data 
OG_data = subset(OG_data, OG_data$Date_Change_2 > 0)

OG_data_22 = subset(OG_data, OG_data$Date_Change_2 == 21)
OG_data = subset(OG_data, OG_data$PATID %in% OG_data_22$PATID)


Final_data = data.frame()
current_Group = "dmLT 2.0 mcg"

for (current_Group in unique(OG_data$Group)) { 
  print(current_Group)
  current_subset = subset(OG_data, OG_data$Group == current_Group)
  current_subset = subset(current_subset, current_subset$TESTTYPA %in% 
                            c("DATI", "DGTI", "TATI", "TGTI", "LTTN"))
  
  current_antibody = "DATI"
  for (current_antibody in c("DATI", "DGTI")) {
    Antibody = subset(current_subset, current_subset$TESTTYPA %in% c(current_antibody,  "LTTN"))
    
    i = 1
    while (i < 3) {
      Current_list = c("SER", "PBS", "LTTN")
      Pair_1 =Current_list[i]
      
      j = i + 1
      while (j < 4) {
        Pair_2 = Current_list[j]
        if (Pair_1 == "SER") {
          Data_1 = subset(Antibody, Antibody$ST == "SER" & Antibody$TESTTYPA != "LTTN")
          Title_1 = ifelse(current_antibody == "DATI", "IgA_Serum", "IgG_Serum")
        } else if (Pair_1 == "PBS") {
          Data_1 = subset(Antibody, Antibody$ST == "PBS")
          Title_1 = ifelse(current_antibody == "DATI", "IgA_ALS", "IgG_ALS")
        } else if (Pair_1 == "LTTN") {
          Data_1 = subset(Antibody, Antibody$TESTTYPA == "LTTN")
          Title_1 = "LTTN"
        }
        
        if (Pair_2 == "PBS") {
          Data_2 = subset(Antibody, Antibody$ST == "PBS")
          Title_2 = ifelse(current_antibody == "DATI", "IgA_ALS", "IgG_ALS")
        } else if (Pair_2 == "LTTN") {
          Data_2 = subset(Antibody, Antibody$TESTTYPA == "LTTN")
          Title_2 = "LTTN"
        }
        
        Data_1 = Data_1[,c("PATID", "Date_Change_2", "ERESULTA")]
        colnames(Data_1) = c("PATID", "Date_Change_2", Title_1)
        
        
        Data_2 = Data_2[,c("PATID", "Date_Change_2", "ERESULTA")]
        colnames(Data_2) = c("PATID", "Date_Change_2", Title_2)
        
        Data_Combined = merge(Data_1, Data_2)
        
        SPEARMAN = cor.test(x = log10(Data_Combined[,Title_1]), 
                            y = log10(Data_Combined[,Title_2]),
                            method = "spearman", exact = F)
        cor_rs = SPEARMAN$estimate
        cor_p = SPEARMAN$p.value
        cor_n = length(unique(Data_Combined$PATID))
        options(scipen = 1); options(digits = 2)
        current_row = data.frame("Antibody" = current_antibody, 
                                 "Test" = paste(Title_1, "vs.", Title_2), 
                                 "n" = cor_n, 
                                 "Spearman" = round(cor_rs,3),
                                 "p" = format(cor_p, scientific = FALSE), 
                                 "Group" = current_Group)
        Final_data = rbind(Final_data, current_row)
        
        Data_Combined$X1 = Data_Combined[,Title_1]
        Data_Combined$X2 = Data_Combined[,Title_2]
        if (current_Group == "dmLT 2.0 mcg") {
          
          current_plot = ggplot(Data_Combined, aes(x = log10(X2), y = log10(X1))) + 
            geom_point() + 
            theme_bw() + 
            xlab(paste("log10(", gsub("_", " ", Title_2), ")", sep = "")) + 
            ylab(paste("log10(", gsub("_", " ", Title_1), ")", sep = "")) + 
            geom_smooth(method='lm', formula= y~x) + 
            labs(title = paste(gsub("_", " ", Title_1),  "vs", gsub("_", " ", Title_2)),
                 subtitle = paste("R2 =",signif(cor_rs, 2),
                                  "| P =",signif(cor_p, 2))); current_plot
          
          assign(paste(Title_1, Title_2, "plot", sep = "_"), current_plot)
          print(paste(Title_1, Title_2, "plot", sep = "_"))
        }
        j = j + 1
      }
      i = i + 1
    }
  }
  
  ###############################
  #Stool
  ###############################
  current_test = "SER"
  for (current_test in c("SER", "PBS", "LTTN")) {
    
    #Combined Fecal Specific and Total 
    Fecal_Subset_Specific = subset(current_subset, current_subset$TESTTYPA %in% c("DATI") & 
                                     current_subset$ST == "FRZ")
    Fecal_Subset_Total = subset(current_subset, current_subset$TESTTYPA %in% c("TATI") & 
                                  current_subset$ST == "FRZ")
    Fecal_Subset_Specific = Fecal_Subset_Specific[,c("PATID", "Date_Change_2", "ERESULTA")]
    colnames(Fecal_Subset_Specific) = c("PATID", "Date_Change_2", "Specific")
    Fecal_Subset_Total = Fecal_Subset_Total[,c("PATID", "Date_Change_2", "ERESULTA")]
    colnames(Fecal_Subset_Total) = c("PATID", "Date_Change_2", "Total")
    
    Fecal_Subset_Combined = merge(Fecal_Subset_Total, Fecal_Subset_Specific)
    #Ratio between Specific and Total 
    Fecal_Subset_Combined$Specific.Total = Fecal_Subset_Combined$Specific/Fecal_Subset_Combined$Total
    
    #Serum Subset
    
    if (current_test == "SER") {
      Data_3 = subset(current_subset, current_subset$TESTTYPA == "DATI" &
                        current_subset$ST == "SER")
      Title_3 = "IgA_Serum"
    } else if (current_test == "PBS") {
      Data_3 = subset(current_subset, current_subset$TESTTYPA == "DATI" &
                        current_subset$ST == "PBS")
      Title_3 = "IgA_ALS"
    } else if (current_test == "LTTN") {
      Data_3 = subset(current_subset, current_subset$TESTTYPA == "LTTN") 
      Title_3 = "LTTN"
    }
    
    Data_3 = Data_3[,c("PATID", "Date_Change_2", "ERESULTA")]
    colnames(Data_3) = c("PATID", "Date_Change_2", Title_3)
    
    Fecal_Other_Combined = merge(Fecal_Subset_Combined, Data_3)
    
    ######
    #Serum vs Fecal Specific/Total 
    SPEARMAN = cor.test(log10(Fecal_Other_Combined$Specific.Total), 
                        log10(Fecal_Other_Combined[,Title_3]),
                        method = "spearman", exact = F)
    SP_value = SPEARMAN$estimate
    SP_p = SPEARMAN$p.value
    SP_n = length(unique(Fecal_Other_Combined$PATID))
    
    current_row = data.frame("Antibody" = "DATI", 
                             "Test" = paste(gsub("_", " ", Title_3), "vs. Fecal Specific/Total"), 
                             "n" = SP_n, 
                             "Spearman" = round(SP_value,3),
                             "p" = format(SP_p, scientific = FALSE), 
                             "Group" = current_Group)
    Final_data = rbind(Final_data, current_row)
    
    Fecal_Other_Combined$Data3 = Fecal_Other_Combined[,Title_3]
    if (current_Group == "dmLT 2.0 mcg") {
      
      current_plot = ggplot(Fecal_Other_Combined, aes(x = log10(Specific.Total), y = log10(Data3))) + 
        geom_point() + 
        theme_bw() + 
        xlab(paste("log10(IgA Fecal Specific/Total)", sep = "")) + 
        ylab(paste("log10(", gsub("_", " ", Title_3), ")", sep = "")) + 
        geom_smooth(method='lm', formula= y~x) +
        labs(title = paste(gsub("_", " ", Title_3), "vs Fecal dmLT Specific/Total"),
             subtitle = paste("R2 =",signif(SP_value, 2),
                              "| P =",signif(SP_p, 2))); current_plot
      
      assign(paste(Title_3, "_Fecal_Specific.Total_plot", sep = ""), current_plot)
      print(paste(Title_3, "_Fecal_Specific.Total_plot", sep = ""))
    }
    ###############
    #Serum vs Fecal Specific
    SPEARMAN = cor.test(log10(Fecal_Other_Combined$Specific), 
                        log10(Fecal_Other_Combined[,Title_3]),
                        method = "spearman", exact = F)
    SP_value = SPEARMAN$estimate
    SP_p = SPEARMAN$p.value
    SP_n = length(unique(Fecal_Other_Combined$PATID))
    
    current_row = data.frame("Antibody" = "DATI", 
                             "Test" = paste(gsub("_", " ", Title_3), "vs. Fecal Specific"), 
                             "n" = SP_n, 
                             "Spearman" = round(SP_value,3),
                             "p" = format(SP_p, scientific = FALSE), 
                             "Group" = current_Group)
    Final_data = rbind(Final_data, current_row)
    
    Fecal_Other_Combined$Data3 = Fecal_Other_Combined[,Title_3]
    
    if (current_Group == "dmLT 2.0 mcg") {
      current_plot = ggplot(Fecal_Other_Combined, aes(x = log10(Specific), y = log10(Data3))) + 
        geom_point() + 
        theme_bw() + 
        xlab(paste("log10(IgA Fecal Specific)", sep = "")) + 
        ylab(paste("log10(", gsub("_", " ", Title_3), ")", sep = "")) + 
        geom_smooth(method='lm', formula= y~x) +
        labs(title = paste(gsub("_", " ", Title_3), "vs Fecal dmLT Specific"),
             subtitle = paste("R2 =",signif(SP_value, 2),
                              "| P =",signif(SP_p, 2))); current_plot
      
      assign(paste(Title_3, "_Fecal_Specific_plot", sep = ""), current_plot)
      print(paste(Title_3, "_Fecal_Specific_plot", sep = ""))
    }
  }
}

Combined_correlations = ggarrange(IgA_Serum_IgA_ALS_plot, IgA_Serum_LTTN_plot, IgA_ALS_LTTN_plot, 
                                  IgG_Serum_IgG_ALS_plot, IgG_Serum_LTTN_plot, IgG_ALS_LTTN_plot, 
                                  IgA_Serum_Fecal_Specific.Total_plot, IgA_ALS_Fecal_Specific.Total_plot, 
                                  LTTN_Fecal_Specific.Total_plot, 
                                  IgA_Serum_Fecal_Specific_plot, IgA_ALS_Fecal_Specific_plot,
                                  LTTN_Fecal_Specific_plot, 
                                  ncol = 3, nrow = 4)


jpeg("./Figures/Final_Figures/Figure_7.jpeg", 
     res= 425, height = 5000, width = 6000)
Combined_correlations
dev.off()


write.csv(Final_data, "Correlation.csv", row.names = F)
