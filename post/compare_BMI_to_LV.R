library(ggplot2)
library(psych); library(pastecs)
library(boot); library(car); library(QuantPsyc)
library(psych); library(pastecs)
library(lme4)
library(nlme)
library(effects)
library(plotly)
library(dvmisc)
packageVersion('plotly')
source("/data/pt_life/lampe/NPM/FSLrandomise/1825_no_weight/test_regression/diagnostic_fcns.r")

data<-read.csv('/data/pt_life/LIFE_bl/publications/2017_beyer_metabolic_pls_ohbm_poster/PLS_on_TBM_obesity_markers/tables/748_allLabvals_MMST_T1_qual_checked_log_transforms_all_regressed_with_LV.csv')


##model 0:

res0<- lm(data$lv_brain ~ data$lv_ob)
summary(res0)


res<- lm(data$lv_brain ~ data$X.RES_BMI)
summary(res)

res2 <- lm(data$lv_brain ~ data$X.RES_BMI + data$lv_ob)

anova(res,res2)

##
ggplot(data, aes(x=data$lv_ob, y=data$lv_brain)) + geom_point(size=4, position=position_jitter(width=0.020), alpha=0.75) +    # Use hollow circles
  geom_smooth(method=lm, colour="#66CC99", size=2, se=FALSE) + labs(x="", y="") +
  theme_bw(base_size=20) + theme( axis.text = element_text(colour = "black"))

##coloring according to cognitive performance.
cognition=read.csv('/data/pt_life/LIFE_bl/publications/2017_beyer_metabolic_pls_ohbm_poster/PLS_on_TBM_obesity_markers/tables/748_allLabvals_MMST_T1_qual_checked_log_transforms_all_regressed_with_cognition.csv')

data$sum_exec=cognition$Z_executive
data$sum_memory=cognition$Z_memory
data$sum_processing=cognition$Z_processingspeed

data$group=quant_groups(data$lv_ob, groups = 2)
levels(data$group)=c("low","high")

##"3D" like plot of latent variables and cognition (with discrete labels)
ggplot(data, aes(x=data$lv_brain, y=data$sum_exec, color = data$group)) + 
  geom_point(size=4, position=position_jitter(width=0.020), alpha=0.75) +    # Use hollow circles
  geom_smooth(method=lm, colour="#66CC99", size=2, se=FALSE) + labs(x="LV brain", y="Z[exec function]",colour="LV obesity") +
  theme_bw(base_size=20) + theme( axis.text = element_text(colour = "black"))
  
#continuous colors (without adipo outlier)
jpeg("/home/raid1/fbeyer/Documents/Results/Metabolic_VBM/Brain_exec_function_obesity_withoutoutlier.jpg", quality=100)
ggplot(data[data$lv_ob>-5,], aes(x=data[data$lv_ob>-5,]$lv_brain, y=data[data$lv_ob>-5,]$sum_exec,color = data[data$lv_ob>-5,]$lv_ob)) + 
  geom_point(size=2, position=position_jitter(width=0.020), alpha=0.75) + 
  scale_color_gradientn(na.value = "grey", colours = c("blue", "yellow", "darkgreen")) +
  #scale_color_gradient(low="red", high="yellow") +
  #scale_color_gradient2(midpoint=0, low="blue", mid="white",high="red", space ="Lab" )+
  labs(x="GMV score", y="Z[exec function]",colour="metabolic score") +
  theme_bw(base_size=20) + theme( axis.text = element_text(colour = "black")) +
  geom_smooth(method=lm, colour="black", size=1, se=FALSE) 
dev.off()

#continuous colors (with adipo outlier)
jpeg("/home/raid1/fbeyer/Documents/Results/Metabolic_VBM/Brain_exec_function_obesity_withoutlier.jpg", quality=100)
ggplot(data, aes(x=data$lv_brain, y=data$sum_exec,color = data$lv_ob)) + 
  geom_point(size=2, position=position_jitter(width=0.020), alpha=0.75) + 
  scale_color_gradientn(na.value = "grey", colours = c("blue", "yellow", "darkgreen")) +
  #scale_color_gradient(low="red", high="yellow") +
  #scale_color_gradient2(midpoint=0, low="blue", mid="white",high="red", space ="Lab" )+
  labs(x="GMV score", y="Z[exec function]",colour="metabolic score") +
  theme_bw(base_size=20) + theme( axis.text = element_text(colour = "black")) +
  geom_smooth(method=lm, colour="black", size=1, se=FALSE) 
dev.off()
##memory score
ggplot(data, aes(x=data$lv_brain, y=data$sum_memory, color = data$group)) + 
  geom_point(size=4, position=position_jitter(width=0.020), alpha=0.75) +    # Use hollow circles
  geom_smooth(method=lm, colour="#66CC99", size=2, se=FALSE) + labs(x="GMV score", y="Z[exec function]",colour="metabolic score") +
  theme_bw(base_size=20) + theme( axis.text = element_text(colour = "black"))


ggplot(data, aes(x=data$lv_brain, y=data$sum_processing, color = data$group)) + 
  geom_point(size=2, position=position_jitter(width=0.020), alpha=0.75) +    # Use hollow circles
  geom_smooth(method=lm, colour="#66CC99", size=2, se=FALSE) + labs(x="LV brain", y="Z[exec function]",colour="LV obesity") +
  theme_bw(base_size=20) + theme( axis.text = element_text(colour = "black"))

###obesity measures andcognition
res<- lm(data$sum_exec ~ data$lv_ob)
summary(res)
lm.beta(res) #To obtain the standardized beta estimates
confint(res) #To get confidence intervals for the parameters in the model

plot(x=fitted(res),y=residuals(res))
plot(residuals(res))
hist(residuals(res),breaks=50,probability = T)
x=seq(from=min(residuals(res)),
      to=max(residuals(res)), length.out=100)
lines(x=x,y=dnorm(x, mean=0, sd=sd(residuals(res))))
qqnorm(residuals(res))
qqline(residuals(res))

##model stability-> looks ok
max(abs(dffits(res)))
head(dfbeta(res))
xx=cbind(coefficients(res),coefficients(res)+t(apply(dfbeta(res),MARGIN=2, FUN=range)))
colnames(xx)=c("orig","min","max")
round(xx,5)
max(as.vector(influence(res)$hat))

res<- lm(data$sum_memory ~ data$lv_ob)
summary(res)
res<- lm(data$sum_processing ~ data$lv_ob)
summary(res)
lm.beta(res)

###brain LV
res<- lm(data$sum_exec ~ data$lv_brain)
summary(res)
lm.beta(res) #To obtain the standardized beta estimates
confint(res)
res<- lm(data$sum_memory ~ data$lv_brain)
summary(res)
res<- lm(data$sum_processing ~ data$lv_brain)
summary(res)
lm.beta(res) 


ggplot(data, aes(x=data$lv_ob, y=data$sum_exec)) + 
  geom_point(size=4, position=position_jitter(width=0.020), alpha=0.75) +    # Use hollow circles
  geom_smooth(method=lm, colour="#66CC99", size=2, se=FALSE) + labs(x="", y="") +
  scale_color_gradient2(midpoint=0, low="blue", mid="white",high="red", space ="Lab" )+
  theme_bw(base_size=20) + theme( axis.text = element_text(colour = "black"))

ggplot(data, aes(x=data$lv_brain, y=data$sum_exec)) + 
  geom_point(size=4, position=position_jitter(width=0.020), alpha=0.75) +    # Use hollow circles
  geom_smooth(method=lm, colour="#66CC99", size=2, se=FALSE) + labs(x="", y="") +
  scale_color_gradient2(midpoint=0, low="blue", mid="white",high="red", space ="Lab" )+
  theme_bw(base_size=20) + theme( axis.text = element_text(colour = "black"))

hist(data$sum_exec)
# p <- plot_ly(data, x = ~lv_ob, y = ~lv_brain, z = ~sumscore, color = ~am, colors = c('#BF382A', '#0C4B8E')) %>%
#   add_markers() %>%
#   layout(scene = list(xaxis = list(title = 'Weight'),
#                       yaxis = list(title = 'Gross horsepower'),
#                       zaxis = list(title = '1/4 mile time')))
# p
# install.packages("scatterplot3d")
# library("scatterplot3d") 
# 
# 
# scatterplot3d(data$lv_brain,data$lv_ob,data$sumscore)

###perform median split
data$lv_brain_cat[data$lv_brain<median(data$lv_brain)]="low"
data$lv_brain_cat[data$lv_brain>=median(data$lv_brain)]="high"

data$lv_brain_ob_cat[data$lv_brain<median(data$lv_brain)&data$lv_ob<median(data$lv_ob)]="low"
data$lv_brain_ob_cat[data$lv_brain>=median(data$lv_brain)&data$lv_ob>=median(data$lv_ob)]="high"

res<- lm(data$sumscore ~ data$lv_brain_ob_cat)
summary(res)

res<- lm(data$ZCERAD_TOTAL_WL_dup_corrected ~ data$lv_brain_ob_cat)
summary(res)

res<- lm(data$ZCERAD_dup_corrected ~ data$lv_brain_ob_cat)
summary(res)

res<- lm(data$ZTMTB._over_TMTA ~ data$lv_brain_ob_cat)
summary(res)



##
gm=read.table('/afs/cbs.mpg.de/share/gr_agingandobesity/life_shared/LIFE_comprehensive_list/Freesurfer_Segmentation_Parcellation_results/FS_results_subcor_LIFE.txt',
            sep="\t", header=T)
colnames(gm)[1]="SIC"

library(haven)
sics=read_sav("/data/pt_life/LIFE_bl/publications/2017_beyer_metabolic_pls_ohbm_poster/PLS_on_TBM_obesity_markers/tables/748_allLabvals_MMST_T1_qual_checked_order_of_mergedfile.sav")


sics$lv_ob=data$lv_ob
sics$lv_brain=data$lv_brain
sics$BMI_res=data$X.RES_BMI

sics_gm <- merge(sics, gm,by="SIC", all.x = T)
sics_gm$TotalGrayVol_adj=sics_gm$TotalGrayVol/sics_gm$EstimatedTotalIntraCranialVol.x

##try predictions with total gray matter volume, ajusted for TIV

##
res0<- lm(sics_gm$lv_brain ~ sics_gm$Age_all + sics_gm$sex_bin + sics_gm$lv_ob)
summary(res0)
res<- lm(sics_gm$lv_brain ~ sics_gm$Age_all + sics_gm$sex_bin + sics_gm$BMI_BMI.1)
summary(res)


res00<- lm(sics_gm$TotalGrayVol_adj ~ sics_gm$Age_all + sics_gm$sex_bin)
summary(res00)

res0<- lm(sics_gm$TotalGrayVol_adj ~ sics_gm$Age_all + sics_gm$sex_bin + sics_gm$lv_ob)
summary(res0)

res<- lm(sics_gm$TotalGrayVol_adj ~ sics_gm$Age_all + sics_gm$sex_bin + sics_gm$BMI_BMI.1)
summary(res)
lm.beta(res)

res2<-lm(sics_gm$TotalGrayVol_adj ~ sics_gm$Age_all + sics_gm$sex_bin + sics_gm$BMI_BMI.1 + sics_gm$lv_ob)
summary(res2)
lm.beta(res2)

res3<- lm(sics_gm$TotalGrayVol_adj ~ sics_gm$lv_ob)
summary(res3)

res4<- lm(sics_gm$TotalGrayVol_adj ~ sics_gm$BMI_res)
summary(res4)


anova(res,res2)
