library(ggplot2)
library(ggridges)
library(viridis)
library(plyr)
library(dplyr)

predictions <- read.csv("DataForAnalysis/predictions.csv")

effector_data <- read.csv("DataForAnalysis/AverageEffectorRepertoires.csv")

actual_repertoires <- read.csv("DataForAnalysis/repertoires.csv")

ani98_phylogroups <- read.csv("DataForAnalysis/ANI_98_phylogroups.csv")
colnames(ani98_phylogroups) <- c("group","phylogroup")

predicted_groups <- unique(predictions$group)

Isolate_T3E_repertoires <- actual_repertoires[ , grepl( "source|Hop|Avr" , names( actual_repertoires ) ) ]

ID_and_Repertoires <- merge(predictions,Isolate_T3E_repertoires,by="source", all=T)

primer_amplifications <- count(ID_and_Repertoires, "primer")

#plot primer classification model performance
ggplot(data=ID_and_Repertoires, aes(x=resolution, y=primer, fill=primer, scale=0.95)) +
  geom_density_ridges(stat = "binline", bins = 5, binwidth = 1,center=0, scale = 0.65, draw_baseline = TRUE, alpha = 1, color = NA) +
  theme_ridges(center_axis_labels = TRUE) + 
  theme(legend.position = "none",panel.grid.major.x = element_blank()) +
  geom_text(
    stat = "bin",
    aes(y = group + 0.7*stat(count/max(count)),
        label = ifelse(stat(count) > 60, stat(count), "")),
    vjust = -0.3, size = 2.5, color = "grey30", binwidth = 1
  ) + 
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(expand = expand_scale(mult = c(0.01, .12)))

ID_and_Repertoires$accuracy <- c(0)
T3E <- colnames(effector_data)
T3E <- T3E[grep( "Hop|Avr" , T3E )]
T3Epredictions <- data.frame('source'='','ANI'='','effector'='', 'present'='', 'correct'='')
T3Epredictions <-T3Epredictions[-1,]

for(i in 1:nrow(ID_and_Repertoires)) {
  unknown_isolate <- ID_and_Repertoires[i,]
  reference_repertoire <- effector_data[effector_data$name==unknown_isolate$group,]
  total_correct <- c()
  for(x in T3E[[1]]){
    if (!x %in% colnames(ID_and_Repertoires)&reference_repertoire[[x]]<0.5){
      present<- 0
      correct <- 1
      }
    else if (!x %in% colnames(ID_and_Repertoires)&reference_repertoire[[x]]>=0.5){
      present <- 0
      correct <- 0
      }
    else if(unknown_isolate[[x]]==0&reference_repertoire[[x]]<0.5){
      present <- 0
      correct <- 1
    }
    else if(unknown_isolate[[x]]==0&reference_repertoire[[x]]>=0.5){
      present <- 0
      correct <- 0
    }
    else if(unknown_isolate[[x]]>0&reference_repertoire[[x]]>=0.5){
      present <- 1
      correct <- 1
    }
    else if(unknown_isolate[[x]]>0&reference_repertoire[[x]]<0.5){
      present <- 1
      correct <- 0
    }
    total_correct <- append(total_correct, correct)
    T3Epredictions[nrow(T3Epredictions) + 1,] = list(unknown_isolate$source,unknown_isolate$resolution,x,present,correct)
  }
  unknown_isolate[["accuracy"]] <- sum(total_correct)/length(total_correct)
  print(i)
  print(sum(total_correct)/length(total_correct))
  ID_and_Repertoires[i,] <- unknown_isolate
  
  gc()
}

mean_prediction_accuracy <- aggregate(ID_and_Repertoires$accuracy, list(ID_and_Repertoires$primer), FUN = function(x) c(median = median(x), range = max(x)-min(x) )) 


ani98_phylogroups<- within(ani98_phylogroups,
                           {phylogroup <- as.character(phylogroup);
                           phylogroup <- ave(phylogroup, group, FUN = toString)})

ID_and_Repertoires_viz_data <- data.frame(ID_and_Repertoires)
ID_and_Repertoires_viz_data$group <- ifelse(nchar(ID_and_Repertoires_viz_data$group) > 19, substring(ID_and_Repertoires_viz_data$group, 1, 19), ID_and_Repertoires_viz_data$group)
ID_and_Repertoires_viz_data <- merge(ID_and_Repertoires_viz_data, ani98_phylogroups, by="group", all.x=FALSE, all.y=FALSE)
ID_and_Repertoires_viz_data <- ID_and_Repertoires_viz_data[!duplicated(ID_and_Repertoires_viz_data[,c('source','primer')]),]
#plot primer T3E prediction ability
ggplot(ID_and_Repertoires, aes(x = primer, y = accuracy)) + 
  geom_violin()


ggplot(data=ID_and_Repertoires, aes(x=accuracy, y=primer, fill=primer, scale=0.95)) +
  stat_density_ridges(from=0.6,to=1.0,quantile_lines = TRUE, quantiles = 2, scale = 0.95, draw_baseline = FALSE, alpha = 1,
                      vline_size = 0.5, vline_color = "black",
                      point_size = 0.4, point_alpha = 1,color = NA,bandwidth = 0.005,center=0) +
  theme_ridges() + 
  theme(legend.position = "none",panel.grid.major.x = element_blank())

# plotting repertoire prediction accuracy by phylogroup
stat_box_data <- function(y) {
  return( 
    data.frame(
      y = 103,  
      label = paste('n =', length(y))
    )
  )
}


ggplot(data=ID_and_Repertoires_viz_data[ID_and_Repertoires_viz_data$resolution>=98 & ID_and_Repertoires_viz_data$primer=="gyrB_Hwang",],aes(x=reorder(phylogroup, -accuracy), y=accuracy*100)) + 
  geom_jitter(width=0.2) + geom_boxplot(alpha=0.2,outlier.alpha=0) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.background = element_blank(),
        panel.grid.major = element_line(color = "grey90"),
        axis.title.x = element_text(vjust=-1.1)) +
  stat_summary(
    fun.data = stat_box_data, 
    color = "grey40",
    geom = "text", 
    hjust = 0.5,
    vjust = 0.5
  ) + xlab("phylogroup") + ylab("percentage of T3E repertoire predicted correctly")

# plotting prediction accuracy by effector protein
T3Epredictions$correct <- as.numeric(T3Epredictions$correct)
T3Epredictions$present <- as.numeric(T3Epredictions$present)
T3Epredictions$Present_wrong <- as.numeric(T3Epredictions$present == 1 & T3Epredictions$correct == 0)
T3Epredictions$Present_right <- as.numeric(T3Epredictions$present == 1 & T3Epredictions$correct == 1)
T3Epredictions$Absent_wrong <- as.numeric(T3Epredictions$present == 0 & T3Epredictions$correct == 0)
T3Epredictions$Absent_right <- as.numeric(T3Epredictions$present == 0 & T3Epredictions$correct == 1)
T3E_prediciton_agg <- aggregate(cbind(Present_right, Present_wrong, Absent_right, Absent_wrong) ~ effector + ANI, data = T3Epredictions,  FUN=mean) 
melted_T3E_agg <- melt(T3E_prediciton_agg[T3E_prediciton_agg$ANI == 98,], id.vars = "effector", measure.vars = c("Present_right", "Present_wrong", "Absent_right", "Absent_wrong"))


t <- melted_T3E_agg[melted_T3E_agg$variable == "Present_right" | melted_T3E_agg$variable == "Absent_right",]
t_agg <- aggregate(value ~ effector, data=t, sum)
t_agg <- t_agg[order(t_agg$value),]
lvls <- t_agg$effector

ggplot(melted_T3E_agg, aes(factor(effector, levels = lvls), y=value, fill=factor(variable, levels=c("Present_right", "Absent_right","Present_wrong", "Absent_wrong")))) +
  geom_bar(position = "fill", stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.background = element_blank(),
        panel.grid.major = element_line(color = "grey90"),
        axis.title.x = element_text(vjust=-1.1))+
  scale_fill_manual(values=c('dodgerblue4', 'lightblue2', 'red4','red')) + guides(fill=guide_legend(title="T3E Predictions"))



ggplot(data=T3E_prediciton_agg[T3E_prediciton_agg$ANI == 98,]) + 
  geom_bar(aes(x=reorder(effector, -Present_right), y=Present_right),stat="identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.background = element_blank(),
        panel.grid.major = element_line(color = "grey90"),
        axis.title.x = element_text(vjust=-1.1)) +
  xlab("Type 3 Effector protein subfamily") + ylab("presence/absense correctly predicted")

  ggplot(data=T3E_prediciton_agg[T3E_prediciton_agg$ANI == 98,]) + 
    geom_point(aes(x=present,y=correct)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.background = element_blank(),
        panel.grid.major = element_line(color = "grey90"),
        axis.title.x = element_text(vjust=-1.1)) 
  