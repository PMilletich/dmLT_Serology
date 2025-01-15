# Author: Trish Milletich, PhD

library(ggplot2)
library(ggpubr)
library(DescTools)

#Preset example Values 
Antibody_list = c("DATI", "DGTI", "LTTN") 

OG_data = read.csv("DataAvailability_Biomarkers.csv")
OG_data$Group = gsub("mcg", "µg", OG_data$Group)

unique(OG_data$Group)
All_data_groups = data.frame()

current_visNo = "12"
current_Test = "DATI"
current_specimen = "PBS"
for (current_Test in Antibody_list) { # c("DATI")){ # 
  print(current_Test)
  Subset_Current = subset(OG_data, OG_data$TESTTYPA == current_Test)
  
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
    Test_Total_FR = data.frame()
    
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
    current_Group = "Placebo"
    for (current_Group in unique(Specimen_Current$Group)) {
      group_Subset = subset(Specimen_Current, Specimen_Current$Group == current_Group)
      Group_N = paste(current_Group, " (N=", length(unique(group_Subset$PATID)), ")", sep = "")
      
      if (current_Group %in% c("Cohort 1:dmLT 0.1 µg",
                               "Cohort 2:dmLT 0.3 µg",
                               "Cohort 3:dmLT 1.0 µg")) {
        date_list = c(1, 8, 22, 29, 43, 50, 64, 71, 99)
      } else if (current_Group == "Cohort 4+5C: dmLT 2.0 µg") {
        date_list = c(1, 8, 22, 29, 43, 50, 64, 71, 99,223)
      } else if (current_Group == "Cohort 5A:One Dose dmLT 2.0 µg") {
        date_list = c(1, 8, 22, 29, 57, 181)
      } else if (current_Group == "Cohort 5B:Two Doses dmLT 2.0 µg"){
        date_list = c( 1, 8, 22, 29, 43, 50, 78, 202)
      } else if (current_Group == "Placebo") {
        date_list = c(1, 8, 22, 29, 43, 50, 57, 64, 71, 78, 99, 181, 202, 223) 
      } else {
        print("Bad Dates")
      }
      
      
      if (current_specimen == "PBS") {
        date_list = subset(date_list, date_list %in% c(1, 8, 29, 50))
      }
      
      table(group_Subset$Date_Change_2)
      
      for (current_date in date_list) {
        if (current_date %in% c(1, 8, 22, 29, 43, 50, 64)) {
          date_threshold = (current_date-(1+3)):(current_date+(1+3))
        } else if (current_date %in% c(57, 78, 99)) {
          date_threshold = (current_date-(4+3)):(current_date+(4+3))
        } else if (current_date %in% c(71)) {
          date_threshold = (current_date-(2+3)):(current_date+(2+3))
        } else if (current_date %in% c(181, 202, 223)) {
          date_threshold = (current_date-(7+3)):(current_date+(7+3))
        } else {
          print("FUCK")
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
          #print(paste(current_Group, current_date,
          #            nrow(subset(date_combined, is.na(date_combined$Seroconversion) == F)), 
          #            nrow(date_combined), Seroconverted_Rate))
          
          #Geometric Means
          Geom_mean = suppressWarnings(Gmean(date_combined$Titer[date_combined$Titer>0],
                                             method = c("classic"), conf.level = 0.95,
                                             sides = c("two.sided"), na.rm = T))
          Geom_mean_FR = suppressWarnings(Gmean(date_combined$Fold.Rise[date_combined$Fold.Rise>0],
                                                method = c("classic"), conf.level = 0.95,
                                                sides = c("two.sided"), na.rm = T))
          
          
          data_mean_row = data.frame("Day" = current_date, 
                                     "Group" = Group_N,
                                     "n" = nrow(date_combined),
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
          name = "Percent Responders",
          sec.axis = sec_axis(~10^ (. * upper_titer), name= expression(paste("Geometric Mean Titer (95% CI)")),
                              breaks = c(1, 10, 100, 1000, 10000, 100000),
                              labels = function(x) format(x, scientific = FALSE))) + 
        geom_line(data = graph_data, aes(x = Day, y=log10(Mean_Titer) / upper_titer)) + 
        geom_point(data = graph_data, aes(x = Day, y = log10(Mean_Titer) / upper_titer))   + 
        geom_linerange(data = graph_data, aes(x = Day, ymin = log10(Lower_Titer) / upper_titer, 
                                             ymax = log10(Upper_Titer) / upper_titer)) + 
        theme(strip.text = element_text(size = 6),
              strip.text.x = element_text(margin = margin(0.5,0,0.5,0, "mm")),
              axis.text.y = element_text(size= 7.5))+
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
          name = "Percent Responders",
          sec.axis = sec_axis(~10^ (. * upper_titer), 
                              name= expression(paste("Geometric Mean Titer (95% CI)")),
                              breaks = c(1, 10, 100, 1000),
                              labels = function(x) format(x, scientific = FALSE))) +
        geom_line(data = graph_data, aes(x = Day, y=log10(Mean_Titer) / upper_titer)) + 
        geom_point(data = graph_data, aes(x = Day, y = log10(Mean_Titer) / upper_titer))   + 
        geom_linerange(data = graph_data, aes(x = Day, ymin = log10(Lower_Titer) / upper_titer, 
                                             ymax = log10(Upper_Titer) / upper_titer)) + 
        theme(strip.text = element_text(size = 6), 
              strip.text.x = element_text(margin = margin(0.5,0,0.5,0, "mm")),
              axis.text.y = element_text(size= 7.5))+
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
          name = "Percent Responders",
          sec.axis = sec_axis(~./100 * Max_Titer, 
                              name= expression(paste("Geometric Mean Titer (95% CI)")))) + 
        geom_line(data = graph_data, aes(x = Day, y=(Mean_Titer / Max_Titer)*100)) + 
        geom_point(data = graph_data, aes(x = Day, y = (Mean_Titer / Max_Titer)*100))   + 
        geom_linerange(data = graph_data, aes(x = Day, ymin = (Lower_Titer / Max_Titer)*100, 
                                             ymax = (Upper_Titer / Max_Titer)*100)) + 
        theme(strip.text = element_text(size = 6),                
              strip.text.x = element_text(margin = margin(0.5,0,0.5,0, "mm")),                
              axis.text.y = element_text(size= 7.5)) ; GeomMean_plot
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
          name = "Percent Responders",
          sec.axis = sec_axis(~10^ (. * upper_FR), name= expression(paste("Geometric Mean FR (95% CI)")),
                              breaks = c(1, 10, 100, 1000, 10000),
                              labels = function(x) format(x, scientific = FALSE))) + 
        geom_line(data = graph_data, aes(x = Day, y=log10(Mean_FR) / upper_FR)) + 
        geom_point(data = graph_data, aes(x = Day, y = log10(Mean_FR) / upper_FR))   + 
        geom_linerange(data = graph_data, aes(x = Day, ymin = log10(Lower_FR) / upper_FR, 
                                             ymax = log10(Upper_FR) / upper_FR)) + 
        theme(strip.text = element_text(size = 6),                
              strip.text.x = element_text(margin = margin(0.5,0,0.5,0, "mm")),                
              axis.text.y = element_text(size= 7.5))+
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
          name = "Percent Responders",
          sec.axis = sec_axis(~10^ (. * upper_FR), name= expression(paste("Geometric Mean FR (95% CI)")),
                              breaks = c(1, 10, 100, 1000),
                              labels = function(x) format(x, scientific = FALSE))) +
        geom_line(data = graph_data, aes(x = Day, y=log10(Mean_FR) / upper_FR)) + 
        geom_point(data = graph_data, aes(x = Day, y = log10(Mean_FR) / upper_FR))   + 
        geom_linerange(data = graph_data, aes(x = Day, ymin = log10(Lower_FR) / upper_FR, 
                                             ymax = log10(Upper_FR) / upper_FR)) + 
        theme(strip.text = element_text(size = 6),                
              strip.text.x = element_text(margin = margin(0.5,0,0.5,0, "mm")),                
              axis.text.y = element_text(size= 7.5)) +
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
          name = "Percent Responders",
          sec.axis = sec_axis(~./100 * Max_FR, 
                              name= expression(paste("Geometric Mean FR (95% CI)")))) + 
        geom_line(data = graph_data, aes(x = Day, y=(Mean_FR / Max_FR)*100)) + 
        geom_point(data = graph_data, aes(x = Day, y = (Mean_FR / Max_FR)*100))   + 
        geom_linerange(data = graph_data, aes(x = Day, ymin = (Lower_FR / Max_FR)*100, 
                                             ymax = (Upper_FR / Max_FR)*100)) + 
        theme(strip.text = element_text(size = 6), 
              strip.text.x = element_text(margin = margin(0.5,0,0.5,0, "mm")), 
              axis.text.y = element_text(size= 7.5)) ; FR_plot
    }
    
    
    FR_plot = annotate_figure(FR_plot, top = paste(current_test_title, " - ", current_specimen_title, sep = ""))
    GeomMean_plot = annotate_figure(GeomMean_plot, top = paste(current_test_title, " - ", current_specimen_title, sep = ""))
    
    
    jpeg(paste("./Figures/", current_Test, "_", current_specimen, "_Titer_10.jpeg", sep = ""), 
         res = 500, height = 3000, width = 2000) 
    plot(GeomMean_plot)
    dev.off()
    
    jpeg(paste("./Figures/", current_Test, "_", current_specimen, "_FoldRise_10.jpeg", sep = ""), 
         res = 500, height = 3000, width = 2000) 
    plot(FR_plot)
    dev.off()
    graph_data$Test = paste(current_Test, current_specimen)
    All_data_groups = rbind(All_data_groups, graph_data)
    
  }
}
