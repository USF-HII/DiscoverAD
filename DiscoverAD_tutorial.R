
library(tidyverse)
library(Boruta)
library(Amelia)
library(cluster)
library(Rtsne)
library(ggrepel)
library(factoextra)
library(reshape2)
library(gridExtra)


#Read in input file
dat<-read.table("DiscoverAD_sample_data.tsv", sep="\t", stringsAsFactors = F, header = T)
head(dat, n = 3)


#Read in variable type file
var_type<-read.table("variable_type.tsv", sep="\t", stringsAsFactors = F, header = T)
head(var_type, n = 5)


#Label patients c-peptide levels as high (0) (>=1.0) or low (<1.0) (c.peptide.ng.ml)
dat$c_peptide<-ifelse((dat$c.peptide.ng.ml<1),1,NA)
dat$c_peptide[(dat$c.peptide.ng.ml>=1)]<-0

#Label as Ab positive if:
# GAD65Ab >0.05 (GAD65Ab.index)
dat$Ab_pos<-ifelse(dat$GAD65Ab.index>.05, 1,0)


#Label patient with weight class (3=obese, 2=overweight, 1=normal)
#Waist circumference >102cm (men) or >88cm (women) (WAIST)
dat$weight_class<-ifelse(((dat$gender=="1" & dat$WAIST>102) | (dat$gender=="2" & dat$WAIST>88)),3,
                         ifelse(((dat$gender=="1" & dat$WAIST>94) | (dat$gender=="2" & dat$WAIST>80)),2,1))

#Set weight class to ordinal variable
dat$weight_class<-ordered(dat$weight_class)

#Label patients with hypertriglyceridemia
#Triglyceride levls >1.7 mmol/L (150.569)
dat$hypertriglyceridemia <- ifelse(dat$trig>150.569,1,0)

#Label patients with low HDL
# HDL <1.03 mmol/L (39.8301) (men) or <1.29 mmol/L (49.8843) (women)
dat$lowHDL<-ifelse(((dat$gender=="1" & dat$hdlc<39.8301) | (dat$gender=="2" & dat$hdlc<49.8843)),1,0)

#Label patients with hypertension (high blood pressure)
# Blood pressure >130/85 mmHg (mcorrsys, mcorrdia)
# Self-reported hypertension (highbp)
dat$hypertension <- ifelse((dat$MCORRSYS>=130 | dat$MCORRDIA>=85 | (dat$HIGHBP==1 & !is.na(dat$HIGHBP))),1,0)

#Label patients with high fasting plasma glucose (>=5.6) or if diabetic
dat$high_fbg<-ifelse(dat$MFBG>=100.8,1,0)
dat$high_fbg<-1


#Label patients with MetS if:
#Waist circumfrence >102cm (men) or >88cm(women) (WAIST)
#And two of the following criteria:
#Hypertriglyceridemia (>1.7 mmol/L (150.569mg/dL))(trig)
#Low HDL ((<1.03 mmol/L (39.8301mg/dL) men, <1.29 mmol/L (49.8843mg/dL) women)) (HDL, meds)
#Hypertension (Systolic bp >= 130mmHg or Diastolic >=85mmHg) (CORRSYS1,CORRDIA1)
#fasting blood glucose >= 5.6 mmol/L (100.8 mg/dL)

for (i in 1:nrow(dat)){
  temp_sum<-sum(dat[i,"hypertriglyceridemia"],dat[i,"lowHDL"],dat[i,"hypertension"],dat[i,"high_fbg"], na.rm = TRUE)
  dat[i,"MetS_sum"]<-temp_sum
  if (dat[i,"weight_class"]==3 & !is.na(dat[i,"weight_class"]) & temp_sum>=2){
    dat[i,"MetS"]<-1
  }
  else {
    dat[i,"MetS"]<-0
  }
}

#Set MetS_sum to ordinal variable
dat$MetS_sum<-ordered(dat$MetS_sum)


#Label as young diagnosis if:
# age of patient at diagnosis (age) was under 21, assign 1 to 'young_diag'
dat$young_diag<-ifelse(dat$age<=21,1,0)

#Label anemic patients
# If low hemoglobin (hgb):
#   Male: <13.5g/dl
#   Female: <12g/dl
# If low hematocrit (hct):
#   Male: 38.8%
#   Female: <34.9%
dat$anemia <- ifelse((((dat$gender==1 & dat$HGB<13.5) | (dat$gender==2 & dat$HGB<12)) & ((dat$gender==1 & dat$HCT<38.8)|(dat$gender==2 & dat$HCT<34.9))),1,0)


#Set a default value for diabetes type
dat$type<-0

#Label patient as Type 1 if:
# patient is autoantibody positive (Ab_pos) AND
# patient is currently taking insulin (insulin_use) OR patient has poor beta cell function
# OR
# patient is autoantibody positive (Ab_pos) AND patient has a c-peptide score <1
dat$type[dat$Ab_pos==1 & (dat$insulin_use==1 | dat$HOMA_beta<=50 | dat$disease_duration<6)]<-1


#Label patient with Type 2 if:
# patient is autoantibody negative (AB_pos) AND
# patient age at diagnosis >= 21 (age) AND
# waist circumference >94cm for men or >80cm for women

dat$type[dat$Ab_pos==0 & dat$age>21 & ((dat$gender=="1" & dat$WAIST>94) | (dat$gender=="2" & dat$WAIST>80))]<-2


dat$type[dat$type==0]<-3

#Add exclusion filter flag
dat$exclusion<-ifelse(dat$type==3, 1, 0)


dat$inclusion <- 0

#Include A+/B+ KPD participants
dat$type[dat$Ab_pos==1 & dat$HOMA_beta>=50 & dat$KETOACID==1]<-3
dat$inclusion<-ifelse(dat$Ab_pos==1 & dat$HOMA_beta>=50 & dat$KETOACID==1, 1, dat$inclusion)

#Include A-/B+ KPD participants
dat$type[dat$Ab_pos==0 & dat$HOMA_beta>=50 & dat$KETOACID==1]<-3
dat$inclusion<-ifelse(dat$Ab_pos==0 & dat$HOMA_beta>=50 & dat$KETOACID==1, 1, dat$inclusion)

#Include A-/B- KPD participants
dat$type[dat$Ab_pos==0 & dat$HOMA_beta<50 & dat$KETOACID==1]<-3
dat$inclusion<-ifelse(dat$Ab_pos==0 & dat$HOMA_beta<50 & dat$KETOACID==1, 1, dat$inclusion)


#Create a binary atypical identification flag
dat$atypical <- ifelse(dat$type==3,1,0)

#Reduce data to only relevant variables for clustering
dat_fs<-dat[, colnames(dat) %in% c(var_type$Variable, "atypical")]

#Reduce data to only potential atypical participants and potential variables for clustering to relevant ones
dat_cluster<-dat[which(dat$atypical==1), colnames(dat) %in% var_type$Variable]

#Set gender to factor
dat_cluster$gender<-as.factor(dat_cluster$gender)


set.seed(123)

#Impute missing values
dat_fs_imputed <- amelia(dat_fs, m=1, parallel = "multicore")


boruta_selection<- Boruta(atypical~., data = dat_fs_imputed$imputations[[1]], doTrace = 2, maxRuns=500)


#Plot boruta selection
gscale<-c("#4D4D4D", "#4D4D4D","#CCCCCC")
plot(boruta_selection, xlab = "", xaxt = "n", col=gscale[as.numeric(boruta_selection$finalDecision)], family="G")
lz<-lapply(1:ncol(boruta_selection$ImpHistory),function(i)
  boruta_selection$ImpHistory[is.finite(boruta_selection$ImpHistory[,i]),i])
names(lz) <- unlist(lapply(colnames(boruta_selection$ImpHistory), str_remove_all, "`"))
Labels <- sort(sapply(lz,median))
axis(side = 1,las=2,labels = names(Labels),
     at = 1:ncol(boruta_selection$ImpHistory), cex.axis = 0.7)
legend("topleft", inset=0.02, legend=c("Confirmed predictor", "Rejected predictor"), fill=gscale[2:3], cex=0.8)
title("Boruta Feature Selection")


#Decide on tentative variables
boruta_selection<-TentativeRoughFix(boruta_selection)


boruta_variables<-unlist(lapply(getSelectedAttributes(boruta_selection), str_remove_all, "`"))

#Filter dat_cluster to only boruta selected variables
dat_cluster<-dat_cluster[,colnames(dat_cluster) %in% boruta_variables]


#Generate vector containing the variable types
ab<-c()
sb<-c()

for (k in 1:ncol(dat_cluster)){
  var1<-colnames(dat_cluster[k])
  binary_type<-var_type[which(var_type$Variable==var1), "Type"]
  if (!is.na(binary_type) & binary_type=="SB"){
    sb<-c(sb, k)
  }
  else if (!is.na(binary_type) & binary_type=="AB"){
    ab<-c(ab, k)
  }
}

type1<-c()

if (!is.null(sb) & !is.null(ab)){
  type1<-list(asymm=ab, symm=sb)
}
if (!is.null(sb) & is.null(ab)){
  type1<-list(symm=sb)
}
if (is.null(sb) & !is.null(ab)){
  type1<-list(asymm=ab)
}

#Create dissimiliarity matrix
gower.daisy.mat<-daisy(dat_cluster, metric="gower", stand= TRUE, type=type1)


#Function to generate silhouette width plots
silwidth<-function(daisy.mat, plot_name){
  sil_width <- c(NA)

  for(i in 2:10){

    pam_fit <- pam(daisy.mat,
                   diss = TRUE,
                   k = i)

    sil_width[i] <- pam_fit$silinfo$avg.width

  }

  # Plot sihouette width (higher is better)
  plot(1:10, sil_width,
       xlab = "Number of clusters",
       ylab = "Silhouette Width",
       main = "Silhouette Width by Cluster for CCHC AD Cohort")
  lines(1:10, sil_width)

  #Return best number of clusters
  best_k<-match(max(sil_width, na.rm = TRUE), sil_width)
  return(best_k)
}

best_k<-silwidth(gower.daisy.mat,'Atypical_sil_width')


print(paste0("Silhouette Width: ", best_k, " clusters"))


#Run PAM clustering using gower matrix and silhouette-derived best k number of clusters
Atypical_pam_fit <- pam(gower.daisy.mat, diss = TRUE, k = best_k)

#Plot silhouette widths
s<-fviz_silhouette(Atypical_pam_fit, label=TRUE) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_grey() +
  scale_color_grey() +
  xlab("Patient #") +
  geom_hline(yintercept=Atypical_pam_fit$silinfo["avg.width"]$avg.width, linetype="dashed") +
  theme(title = element_text(family = 'Arial'))

plot(s)
```

#Generate visual representation of cluster for Atypical
tsne_obj <- Rtsne(gower.daisy.mat, is_distance = TRUE, perplexity=6)

tsne_data <- tsne_obj$Y %>%
  data.frame() %>%
  setNames(c("X", "Y")) %>%
  mutate(cluster = factor(Atypical_pam_fit$clustering),
         name = dat_cluster$BDVISIT,
         pid = names(Atypical_pam_fit$clustering))

tsne<-ggplot(aes(x = X, y = Y), data = tsne_data) + 
  geom_point(aes(shape=cluster), size=4) +
  scale_shape_manual(values=c(15,2,19,5)) +
  labs(fill="Cluster") +
  xlab("Tsne Dimension 1") +
  ylab("Tsne Dimension 2") +
  ggtitle("AD PAM Clusters") +
  theme(plot.title = element_text(hjust = 0.6)) + 
  theme(title = element_text(family = 'Arial')) +
  geom_text_repel(aes(label=pid), point.padding = .25)
plot(tsne)


Atypical_pam_results <- dat_cluster %>%
  mutate(cluster = Atypical_pam_fit$clustering) %>%
  group_by(cluster) %>%
  do(the_summary = summary(.))

print(Atypical_pam_results$the_summary)

#Visualize clusters with radar plots
           
theme(plot.title = element_text(hjust = 0.5))

## Data prep

#Filter to only cluster variables and gender
dat_cluster_radar <- dat_cluster

dat_cluster_radar <- merge(dat_cluster_radar, dat[,c("Mask", "consensus_gender")], by.x="row.names", by.y="Mask")
rownames(dat_cluster_radar) <- dat_cluster_radar$Row.names

#Convert consensus_gender to Gender (Male)
dat_cluster_radar$Gender <- ifelse(dat_cluster_radar$consensus_gender==1,1,0)

#Drop unnecessary variables
dat_cluster_radar <- dplyr::select(dat_cluster_radar, -c("Row.names","consensus_gender"))

#Get variable types
var_types <- tbl_summary(dplyr::select(dat_cluster_radar,-"Cluster"))
var_types <- var_types$table_body
var_types <- data.frame(distinct(var_types[, c("variable", "var_type")]))

#Output table
radar_df <- as.data.frame(matrix(nrow=length(unique(dat_cluster_radar$Cluster)), ncol=nrow(var_types)))
colnames(radar_df) <- var_types$variable

#Get stats for dichotomous variables
for (var_name in var_types[which(var_types$var_type=="dichotomous"),"variable"]) {
  for (cluster_id in unique(dat_cluster_radar$Cluster)){
    dat_cluster_radar_var <- as.data.frame(dat_cluster_radar[which(dat_cluster_radar$Cluster==cluster_id), var_name])
    mean_avg <- mean(dat_cluster_radar_var[,1], na.rm = T)
    radar_df[cluster_id, var_name] <- mean_avg
  }
}

#Get stats for continuous variables
continuous_vars <- dplyr::select(dat_cluster_radar, c(var_types[which(var_types$var_type=="continuous"),"variable"], "Cluster"))

continuous_vars_norm <- as.data.frame(matrix(nrow=nrow(continuous_vars), ncol=ncol(continuous_vars)))
colnames(continuous_vars_norm) <- colnames(continuous_vars)
continuous_vars_norm$Cluster <- continuous_vars$Cluster
rownames(continuous_vars_norm) <- rownames(dat_cluster_radar)

for (j in 1:(ncol(continuous_vars)-1)){
  x <- continuous_vars[,j]
  normalized <- (x-min(x, na.rm = T))/(max(x, na.rm = T)-min(x, na.rm = T))
  continuous_vars_norm[,j] <- normalized
}


for (var_name in var_types[which(var_types$var_type=="continuous"),"variable"]) {
  for (cluster_id in unique(dat_cluster_radar$Cluster)){
    dat_cluster_radar_var <- as.data.frame(continuous_vars_norm[which(continuous_vars_norm$Cluster==cluster_id), var_name])
    median_avg <- median(dat_cluster_radar_var[,1], na.rm = T)
    sd <- sd(dat_cluster_radar_var[,1], na.rm = T)
    radar_df[cluster_id, var_name] <- median_avg
  }
}

radar_df_t <- as.data.frame(t(radar_df))
colnames(radar_df_t)<-paste("Cluster", seq(1,4), sep="_")
radar_df_t$variable <- rownames(radar_df_t)

coord_radar <- function (theta = "x", start = 0, direction = 1) {
  theta <- match.arg(theta, c("x", "y"))
  r <- if (theta == "x") "y" else "x"
  ggproto("CordRadar", CoordPolar, theta = theta, r = r, start = start, 
          direction = sign(direction),
          clip="off",
          is_linear = function(coord) TRUE)
}

radar_long <- reshape2::melt(radar_df_t)
colnames(radar_long) <- c("Variable", "Cluster", "Value")

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

rainbow4 <- gg_color_hue(4)
ybreaks <- seq(0,1,.25)

plot_list <- list()

count<-1

for (cluster_id in unique(radar_long$Cluster)){
  radar_long_tmp <- radar_long[which(radar_long$Cluster==cluster_id),]
  radar_long_tmp <- radar_long_tmp[order(radar_long_tmp$Variable),]
  radar_long_tmp[nrow(radar_long_tmp)+1,] <- radar_long_tmp[1,]
  
  p <- ggplot(data=radar_long_tmp, aes(x=Variable, y=Value)) +
    geom_polygon(aes(group=Cluster, color="Cluster"), color=rainbow4[count], fill = rainbow4[count], alpha=0.5, size = 1, show.legend = FALSE) +
    geom_line(aes(group=Cluster), color=rainbow4[count], size=.75) +
    ylim(-0.25, 1) +
    geom_point(size=1.5, color=rainbow4[count]) +
    coord_radar(start=0) +
    theme_minimal() +
    ylab("") +
    xlab("") +
    geom_text(data = data.frame(x = 0, y = ybreaks, label = ybreaks),
              aes(x = x, y = y, label = label),
              size =4) +
    ggtitle(str_replace(cluster_id, "_"," ")) +
    theme(legend.position = "none",
          axis.ticks = element_blank(),
          axis.text.y = element_blank(),
          panel.grid.major.x = element_line(colour="white"),
          plot.title = element_text(hjust = 0.5, face="bold", size=24),
          axis.text.x = element_text(size=10)
    )
    
  plot_list[[count]]<-p
  
  count <- count + 1
}

plot_list <- lapply(plot_list, ggplotGrob)
x <- arrangeGrob(grobs=plot_list, nrow=2, ncol=2)
grid.arrange(x)
