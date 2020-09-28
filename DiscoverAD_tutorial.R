
library(tidyverse)
library(Boruta)
library(Amelia)
library(cluster)
library(Rtsne)
library(ggrepel)


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


#Include A-/B+ and A-/B- KPD participants
dat$type[dat$Ab_pos==0 & dat$KETOACID==1]<-3
dat$inclusion<-ifelse(dat$Ab_pos==0 & dat$KETOACID==1, 1, 0)

#Include A+/B+ KPD participants
dat$type[dat$Ab_pos==1 & dat$HOMA_beta>50 & dat$KETOACID==1]<-3
dat$inclusion<-ifelse(dat$Ab_pos==1 & dat$HOMA_beta>50 & dat$KETOACID==1, 1, dat$inclusion)


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
plot(boruta_selection, xlab = "", xaxt = "n")
lz<-lapply(1:ncol(boruta_selection$ImpHistory),function(i)
  boruta_selection$ImpHistory[is.finite(boruta_selection$ImpHistory[,i]),i])
names(lz) <- unlist(lapply(colnames(boruta_selection$ImpHistory), str_remove_all, "`"))
Labels <- sort(sapply(lz,median))
axis(side = 1,las=2,labels = names(Labels),
     at = 1:ncol(boruta_selection$ImpHistory), cex.axis = 0.7)


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
plot(Atypical_pam_fit, main=NULL)


#Generate visual representation of cluster for Atypical
tsne_obj <- Rtsne(gower.daisy.mat, is_distance = TRUE, perplexity=6)

tsne_data <- tsne_obj$Y %>%
  data.frame() %>%
  setNames(c("X", "Y")) %>%
  mutate(cluster = factor(Atypical_pam_fit$clustering),
         name = dat_cluster$BDVISIT,
         pid = names(Atypical_pam_fit$clustering))

print(ggplot(aes(x = X, y = Y), data = tsne_data) +
        geom_point(aes(color=cluster), size=2) +
        theme(text = element_text(size=18)) +
        labs(color="Cluster"))


Atypical_pam_results <- dat_cluster %>%
  mutate(cluster = Atypical_pam_fit$clustering) %>%
  group_by(cluster) %>%
  do(the_summary = summary(.))

print(Atypical_pam_results$the_summary)

