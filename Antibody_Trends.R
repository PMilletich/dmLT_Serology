#' Patricia Milletich 
#' July, 2024

library(ggplot2)
library(ggpubr)
library(DescTools)

#Preset example Values 
current_Test = "DATI"
current_specimen = "SER"
current_Group = "dmLT 0.1 mcg"
Antibody_list = c("DATI", "DGTI", "LTTN") 

OG_data = read.csv("EDU_1.csv")
#Fold Change as.Date(as.character(survey$tx_start), format="%Y/%m/%d")
OG_data$STARTDT_1 = as.Date(OG_data$STARTDT, format="%m/%d/%Y")
OG_data$COL_DATE_1 = as.Date(OG_data$COL_DATE, format="%m/%d/%Y")
OG_data$Date_Change = OG_data$COL_DATE_1 - OG_data$STARTDT_1
OG_data$Date_Change_2 = as.numeric(OG_data$Date_Change)+1

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

# Other Exclusions for Specific Time Points:
# 08ECI017 - Day 29 visit was -3 days out of window; exclude Day 29
# 08ECI151 - Third Vaccination 6 days out of window; exclude Day 43 through Day 223
# 08ECI214 - Day 202 visit -12 days out of window; exclude Day 202
# 08ECI236 - Day 202 visit 16 days out of window; exclude Day 202
# 08ECI238 - Day 202 visit 16 days out of window; exclude Day 202

#Rename the Treatment groups 
table(OG_data$TRTTRUE)
OG_data$TRTTRUE = gsub("One Dose of Intradermal dmLT 2.0 mcg", "dmLT 2.0 mcg: One Dose", OG_data$TRTTRUE)
OG_data$TRTTRUE = gsub("Two Doses of Intradermal dmLT 2.0 mcg", "dmLT 2.0 mcg: Two Doses", OG_data$TRTTRUE)
OG_data$TRTTRUE = gsub("Three Doses of Intradermal dmLT 2.0 mcg", "dmLT 2.0 mcg", OG_data$TRTTRUE)
OG_data$TRTTRUE = gsub("Intradermal ", "", OG_data$TRTTRUE)

OG_data$Group = ifelse(grepl("Placebo", OG_data$TRTTRUE), "Placebo", OG_data$TRTTRUE)

for (current_Test in Antibody_list) {
  print(current_Test)
  Subset_Current = subset(OG_data, OG_data$TESTTYPA == current_Test)
  Subset_Current = subset(Subset_Current, Subset_Current$ST != "FRZ")
  
  if (current_Test == "DATI") {
    current_test_title = "dmLT-specific IgA ELISA"
  } else if (current_Test == "DGTI") {
    current_test_title = "dmLT-specific IgG ELISA"
  } else if (current_Test == "LTTN") {
    current_test_title = "LT Toxin Neutralization"
  } else {
    current_test_title = "Whoops"
  }
  
  # Loop through SER or PBS as appropriates
  for (current_specimen in unique(Subset_Current$ST)) {
    print(current_specimen)
    
    if (current_specimen == "PBS") {
      current_specimen_title = "ALS"
      FC_Threshold = 2
    } else if (current_specimen == "SER") {
      current_specimen_title = "Serum"
      FC_Threshold = 4
    } else {
      current_specimen_title = "Whoops"
    }
    
    Specimen_Current = subset(Subset_Current, Subset_Current$ST == current_specimen)
    
    graph_data = data.frame()
    Specimen_Current_3 = data.frame()
    for (current_Group in unique(Specimen_Current$Group)) {
      group_Subset = subset(Specimen_Current, Specimen_Current$Group == current_Group)
      Group_N = paste(current_Group, " (N=", length(unique(group_Subset$PATID)), ")", sep = "")
      
      if (current_Group %in% c("dmLT 0.1 mcg","dmLT 0.3 mcg","dmLT 1.0 mcg")) {
        date_list = c(1, 8, 22, 29, 43, 50, 64, 71, 99)
      } else if (current_Group == "dmLT 2.0 mcg") {
        date_list = c(1, 8, 22, 29, 43, 50, 64, 71, 99,223)
      } else if (current_Group == "dmLT 2.0 mcg: One Dose") {
        date_list = c(1, 8, 22, 29, 57, 181)
      } else if (current_Group == "dmLT 2.0 mcg: Two Doses"){
        date_list = c( 1, 8, 22, 29, 43, 50, 78, 202)
      } else if (current_Group == "Placebo") {
        date_list = c(1, 8, 22, 29, 43, 50, 57, 64, 71, 78, 99, 181, 202, 223) 
      } else {
        print("Bad Dates")
      }
      
      for (current_date in date_list) {
        if (current_date %in% c(1, 8, 22, 29, 43, 50, 64)) {
          date_threshold = (current_date-1):(current_date+1)
        } else if (current_date %in% c(57, 78, 99)) {
          date_threshold = (current_date-4):(current_date+4)
        } else if (current_date %in% c(71)) {
          date_threshold = (current_date-2):(current_date+2)
        } else if (current_date %in% c(181, 202, 223)) {
          date_threshold = (current_date-7):(current_date+7)
        }
        
        date_subset = subset(group_Subset, group_Subset$Date_Change_2 %in% date_threshold)
        if(nrow(date_subset) > 0) {
          date_subset$Titer = date_subset$ERESULTA
          date_subset = date_subset[,c("PATID", "Titer","ST","Group")]
          
          baseline = subset(group_Subset, group_Subset$PATID %in% date_subset$PATID & 
                              group_Subset$Date_Change_2 == 1)
          baseline$Baseline = baseline$ERESULTA
          baseline = baseline[,c("PATID", "Baseline","ST","Group")]
          
          date_combined = merge(date_subset, baseline)
          date_combined$Study_day = current_date
          date_combined$Fold.Rise = date_combined$Titer/date_combined$Baseline
          date_combined$Seroconversion = ifelse(date_combined$Fold.Rise >= FC_Threshold, "Seroconversion", NA)
          
          Specimen_Current_3 = rbind(Specimen_Current_3, date_combined)
          
          #Graph Data
          #Serconversion prevalence
          Seroconverted_Rate = round(nrow(subset(date_combined, is.na(date_combined$Seroconversion) == F))/nrow(date_combined),3)*100
          
          #Geometric Means
          Geom_mean = suppressWarnings(Gmean(date_combined$Titer[date_combined$Titer>0],
                                             method = c("classic"), conf.level = 0.95,
                                             sides = c("two.sided"), na.rm = T))
          Geom_mean_FR = suppressWarnings(Gmean(date_combined$Fold.Rise[date_combined$Fold.Rise>0],
                                                method = c("classic"), conf.level = 0.95,
                                                sides = c("two.sided"), na.rm = T))
          
          
          data_mean_row = data.frame("Day" = current_date, 
                                     "Group" = Group_N,
                                     "Seroconverted"  = Seroconverted_Rate,
                                     "Mean_Titer" = Geom_mean[1],
                                     "Lower_Titer" = Geom_mean[2],
                                     "Upper_Titer" = Geom_mean[3], 
                                     "Mean_FR" = Geom_mean_FR[1],
                                     "Lower_FR" = Geom_mean_FR[2],
                                     "Upper_FR" = Geom_mean_FR[3])
          graph_data = rbind(graph_data, data_mean_row)
        }
      } #Date Loop
    } #Group loop
    
    
    graph_data_Cases = subset(graph_data, grepl("Placebo", graph_data$Group) == F)
    set.seed(1)
    
    Max_Titer = max(c(graph_data_Cases$Mean_Titer, graph_data_Cases$Lower_Titer, 
                      graph_data_Cases$Upper_Titer), na.rm= T)
    upper_titer <- log10(Max_Titer) / 100
    
    breakfun <- function(x) {
      round(10^scales::extended_breaks()(log10(x)))
    }
    
    #graph_data$Log_LowerTiter = ifelse(graph_data$Log_LowerTiter <0, 0, graph_data$Log_LowerTiter)
    
    if (Max_Titer > 1000) {
      GeomMean_plot = ggplot() +  
        facet_wrap(~Group, ncol = 1) + 
        theme_bw() + 
        geom_vline(xintercept = 1, color = "lightcoral", linewidth = 1) + 
        geom_vline(data=subset(graph_data, grepl("One Dose", graph_data$Group) == F),
                   aes(xintercept=22), colour="lightcoral", linewidth = 1) +
        geom_vline(data=subset(graph_data, grepl("One Dose", graph_data$Group) == F &
                                 grepl("Two Dose", graph_data$Group) == F),
                   aes(xintercept=43), colour="lightcoral", linewidth = 1) + 
        geom_bar(data = graph_data, aes(x = Day, y  = Seroconverted), stat = "identity", 
                 alpha= 0.25, fill = "grey", color = "grey44", width = 4)+
        geom_hline(yintercept = 1, color = "black", linewidth = 1)+ 
        xlab("Study Day") + 
        scale_y_continuous(
          position = "right",
          name = "Percent Seroconverted",
          sec.axis = sec_axis(~10^ (. * upper_titer), name= expression(paste("Geometric Mean Titer (95% CI)")),
                              breaks = c(1, 10, 100, 1000, 10000, 100000),
                              labels = function(x) format(x, scientific = FALSE))) + 
        geom_line(data = graph_data, aes(x = Day, y=log10(Mean_Titer) / upper_titer)) + 
        geom_point(data = graph_data, aes(x = Day, y = log10(Mean_Titer) / upper_titer))   + 
        geom_linerange(data = graph_data, aes(x = Day, ymin = log10(Lower_Titer) / upper_titer, 
                                             ymax = log10(Upper_Titer) / upper_titer)) + 
        theme(strip.text = element_text(size = 6),axis.text.y = element_text(size= 7.5)) + 
        coord_cartesian(ylim = c(0, 100)); GeomMean_plot
    } else if (Max_Titer > 150) {
      GeomMean_plot = ggplot() +  
        facet_wrap(~Group, ncol = 1) + 
        theme_bw() + 
        geom_vline(xintercept = 1, color = "lightcoral", linewidth = 1) + 
        geom_vline(data=subset(graph_data, grepl("One Dose", graph_data$Group) == F),
                   aes(xintercept=22), colour="lightcoral", linewidth = 1) +
        geom_vline(data=subset(graph_data, grepl("One Dose", graph_data$Group) == F &
                                 grepl("Two Dose", graph_data$Group) == F),
                   aes(xintercept=43), colour="lightcoral", linewidth = 1) + 
        geom_bar(data = graph_data, aes(x = Day, y  = Seroconverted), stat = "identity", 
                 alpha= 0.25, fill = "grey", color = "grey44", width = 4)+
        geom_hline(yintercept = 1, color = "black", linewidth = 1)+ 
        xlab("Study Day") + 
        scale_y_continuous(
          position = "right",
          name = "Percent Seroconverted",
          sec.axis = sec_axis(~10^ (. * upper_titer), 
                              name= expression(paste("Geometric Mean Titer (95% CI)")),
                              breaks = c(1, 10, 100, 1000),
                              labels = function(x) format(x, scientific = FALSE))) +
        geom_line(data = graph_data, aes(x = Day, y=log10(Mean_Titer) / upper_titer)) + 
        geom_point(data = graph_data, aes(x = Day, y = log10(Mean_Titer) / upper_titer))   + 
        geom_linerange(data = graph_data, aes(x = Day, ymin = log10(Lower_Titer) / upper_titer, 
                                             ymax = log10(Upper_Titer) / upper_titer)) + 
        theme(strip.text = element_text(size = 6),axis.text.y = element_text(size= 7.5)) + 
        coord_cartesian(ylim = c(0, 100)); GeomMean_plot
    } else  {
      GeomMean_plot = ggplot() +  
        facet_wrap(~Group, ncol = 1) + 
        theme_bw() + 
        geom_vline(xintercept = 1, color = "lightcoral", linewidth = 1) + 
        geom_vline(data=subset(graph_data, grepl("One Dose", graph_data$Group) == F),
                   aes(xintercept=22), colour="lightcoral", linewidth = 1) +
        geom_vline(data=subset(graph_data, grepl("One Dose", graph_data$Group) == F &
                                 grepl("Two Dose", graph_data$Group) == F),
                   aes(xintercept=43), colour="lightcoral", linewidth = 1) + 
        geom_bar(data = graph_data, aes(x = Day, y  = Seroconverted), stat = "identity", 
                 alpha= 0.25, fill = "grey", color = "grey44", width = 4)+
        geom_hline(yintercept = 1, color = "black", linewidth = 1)+ 
        xlab("Study Day") + 
        scale_y_continuous(
          position = "right",
          name = "Percent Seroconverted",
          sec.axis = sec_axis(~./100 * Max_Titer, 
                              name= expression(paste("Geometric Mean Titer (95% CI)")))) + 
        geom_line(data = graph_data, aes(x = Day, y=(Mean_Titer / Max_Titer)*100)) + 
        geom_point(data = graph_data, aes(x = Day, y = (Mean_Titer / Max_Titer)*100))   + 
        geom_linerange(data = graph_data, aes(x = Day, ymin = (Lower_Titer / Max_Titer)*100, 
                                             ymax = (Upper_Titer / Max_Titer)*100)) + 
        theme(strip.text = element_text(size = 6),axis.text.y = element_text(size= 7.5)) ; GeomMean_plot
    }
    
    ######################
    # Fold Rise 
    ######################
    Max_FR = max(c(graph_data_Cases$Mean_FR, graph_data_Cases$Lower_FR, 
                   graph_data_Cases$Upper_FR), na.rm= T)
    upper_FR <- log10(Max_FR) / 100
    
    
    if (Max_FR > 1000) {
      FR_plot = ggplot() +  
        facet_wrap(~Group, ncol = 1) + 
        theme_bw() + 
        geom_vline(xintercept = 1, color = "lightcoral", linewidth = 1) + 
        geom_vline(data=subset(graph_data, grepl("One Dose", graph_data$Group) == F),
                   aes(xintercept=22), colour="lightcoral", linewidth = 1) +
        geom_vline(data=subset(graph_data, grepl("One Dose", graph_data$Group) == F &
                                 grepl("Two Dose", graph_data$Group) == F),
                   aes(xintercept=43), colour="lightcoral", linewidth = 1) + 
        geom_bar(data = graph_data, aes(x = Day, y  = Seroconverted), stat = "identity", 
                 alpha= 0.25, fill = "grey", color = "grey44", width = 4)+
        geom_hline(yintercept = 1, color = "black", linewidth = 1)+ 
        xlab("Study Day") + 
        scale_y_continuous(
          position = "right",
          name = "Percent Seroconverted",
          sec.axis = sec_axis(~10^ (. * upper_FR), name= expression(paste("Geometric Mean FR (95% CI)")),
                              breaks = c(1, 10, 100, 1000, 10000),
                              labels = function(x) format(x, scientific = FALSE))) + 
        geom_line(data = graph_data, aes(x = Day, y=log10(Mean_FR) / upper_FR)) + 
        geom_point(data = graph_data, aes(x = Day, y = log10(Mean_FR) / upper_FR))   + 
        geom_linerange(data = graph_data, aes(x = Day, ymin = log10(Lower_FR) / upper_FR, 
                                             ymax = log10(Upper_FR) / upper_FR)) + 
        theme(strip.text = element_text(size = 6),axis.text.y = element_text(size= 7.5)) + 
        coord_cartesian(ylim = c(0, 100)); FR_plot
    } else if (Max_FR > 150) {
      FR_plot = ggplot() +  
        facet_wrap(~Group, ncol = 1) + 
        theme_bw() + 
        geom_vline(xintercept = 1, color = "lightcoral", linewidth = 1) + 
        geom_vline(data=subset(graph_data, grepl("One Dose", graph_data$Group) == F),
                   aes(xintercept=22), colour="lightcoral", linewidth = 1) +
        geom_vline(data=subset(graph_data, grepl("One Dose", graph_data$Group) == F &
                                 grepl("Two Dose", graph_data$Group) == F),
                   aes(xintercept=43), colour="lightcoral", linewidth = 1) + 
        geom_bar(data = graph_data, aes(x = Day, y  = Seroconverted), stat = "identity", 
                 alpha= 0.25, fill = "grey", color = "grey44", width = 4)+
        geom_hline(yintercept = 1, color = "black", linewidth = 1)+ 
        xlab("Study Day") + 
        scale_y_continuous(
          position = "right",
          name = "Percent Seroconverted",
          sec.axis = sec_axis(~10^ (. * upper_FR), name= expression(paste("Geometric Mean FR (95% CI)")),
                              breaks = c(1, 10, 100, 1000),
                              labels = function(x) format(x, scientific = FALSE))) +
        geom_line(data = graph_data, aes(x = Day, y=log10(Mean_FR) / upper_FR)) + 
        geom_point(data = graph_data, aes(x = Day, y = log10(Mean_FR) / upper_FR))   + 
        geom_linerange(data = graph_data, aes(x = Day, ymin = log10(Lower_FR) / upper_FR, 
                                             ymax = log10(Upper_FR) / upper_FR)) + 
        theme(strip.text = element_text(size = 6),axis.text.y = element_text(size= 7.5)) + 
        coord_cartesian(ylim = c(0, 100)); FR_plot
    } else {
      FR_plot = ggplot() +  
        facet_wrap(~Group, ncol = 1) + 
        theme_bw() + 
        geom_vline(xintercept = 1, color = "lightcoral", linewidth = 1) + 
        geom_vline(data=subset(graph_data, grepl("One Dose", graph_data$Group) == F),
                   aes(xintercept=22), colour="lightcoral", linewidth = 1) +
        geom_vline(data=subset(graph_data, grepl("One Dose", graph_data$Group) == F &
                                 grepl("Two Dose", graph_data$Group) == F),
                   aes(xintercept=43), colour="lightcoral", linewidth = 1) + 
        geom_bar(data = graph_data, aes(x = Day, y  = Seroconverted), stat = "identity", 
                 alpha= 0.25, fill = "grey", color = "grey44", width = 4)+
        geom_hline(yintercept = 1, color = "black", linewidth = 1)+ 
        xlab("Study Day") + 
        scale_y_continuous(
          position = "right",
          name = "Percent Seroconverted",
          sec.axis = sec_axis(~./100 * Max_FR, 
                              name= expression(paste("Geometric Mean FR (95% CI)")))) + 
        geom_line(data = graph_data, aes(x = Day, y=(Mean_FR / Max_FR)*100)) + 
        geom_point(data = graph_data, aes(x = Day, y = (Mean_FR / Max_FR)*100))   + 
        geom_linerange(data = graph_data, aes(x = Day, ymin = (Lower_FR / Max_FR)*100, 
                                             ymax = (Upper_FR / Max_FR)*100)) + 
        theme(strip.text = element_text(size = 6),axis.text.y = element_text(size= 7.5)) ; FR_plot
    }
    
    
    FR_plot = annotate_figure(FR_plot, top = paste(current_test_title, " - ", current_specimen_title, sep = ""))
    GeomMean_plot = annotate_figure(GeomMean_plot, top = paste(current_test_title, " - ", current_specimen_title, sep = ""))
    
    
    jpeg(paste("./Figures/", current_Test, "_", current_specimen, "_Titer_3.jpeg", sep = ""), 
         res = 500, height = 3000, width = 2000) 
    plot(GeomMean_plot)
    dev.off()
    
    jpeg(paste("./Figures/", current_Test, "_", current_specimen, "_FoldRise_3.jpeg", sep = ""), 
         res = 500, height = 3000, width = 2000) 
    plot(FR_plot)
    dev.off()
  }
}

Subset_Current = subset(OG_data, OG_data$TESTTYPA == "DGTI")

