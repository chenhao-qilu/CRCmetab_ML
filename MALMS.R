rm(list = ls())
###dataset prepare work###

library(org.Hs.eg.db)
g_symbol <- mget(colnames(merge), 
                 org.Hs.egSYMBOL, 
                 ifnotfound=NA)
g_symbol<-as.character(g_symbol)
library(GSEABase)
gene<-getGmt("Reactome_Metabolism.gmt")
gene<-sapply(gene,function(x){x@geneIds})
a<-NULL
for( i in c(1:84)){
  b<-gene[[i]]
  a<-union(a,b)
}
kegggene<-unique(a)


####cox####
merge<-merge()[,c("OS","OS.time",metageneintersect)]
Coxoutput <- NULL 
for(i in 3:ncol(merge)){
  g <- colnames(merge)[i]
  cox <- coxph(Surv(OS.time,OS) ~ merge[,i], data = merge) # 单变量cox模型
  coxSummary = summary(cox)
  
  Coxoutput <- rbind.data.frame(Coxoutput,
                                data.frame(gene = g,
                                           HR = as.numeric(coxSummary$coefficients[,"exp(coef)"])[1],
                                           z = as.numeric(coxSummary$coefficients[,"z"])[1],
                                           pvalue = as.numeric(coxSummary$coefficients[,"Pr(>|z|)"])[1],
                                           lower = as.numeric(coxSummary$conf.int[,3][1]),
                                           upper = as.numeric(coxSummary$conf.int[,4][1]),
                                           stringsAsFactors = F),
                                stringsAsFactors = F)
}
save(merge,file = "0321combpfs.rda")
Coxoutput <- Coxoutput[which(Coxoutput$pvalue < 0.05),]
prevar<-intersect(Coxoutput2$gene,metageneintersect)


####machinelearning framework####
library(tidyverse)
library(survival)
library(randomForestSRC)
library(glmnet)
library(plsRcox)
library(superpc)
library(gbm)
library(CoxBoost)
library(survivalsvm)
library(dplyr)
library(tibble)

prevar<-intersect(Coxoutput$gene,colnames(d$GSE29621))
prevar<-intersect(prevar,colnames(d$GSE103479))

trans<-lapply(trans,function(x)scale(x))

mm <- list(merge,gse1,gse2,gse9,gse143,gse161)
names(mm)<-c("meta-cohort","GSE17536","GSE29621","GSE92921","GSE143985","GSE161158")


val_data_list <- lapply(val_data_list,function(x){x[,c('OS.time','OS',prevar)]})
val_data_list<-lapply(val_data_list,function(x){x%>%drop_na()})
val_dd_list <- lapply(val_data_list,function(x){x[,c('OS.time','OS',prevar)]})
val_dd_list<-lapply(val_dd_list,function(x){x%>%drop_na()})
save(val_dd_list,file = "0603shujuji.rda")
result <- data.frame()
est_data <- val_dd_list$`meta-cohort`
pre_var <- colnames(est_data)[-c(1:2)]
est_dd <- est_data[,c('OS.time','OS',pre_var)]


rf_nodesize <- 5
seed <- 123

####RSF ####

a<-randomForestSRC::tune.rfsrc(Surv(OS.time,OS) ~ ., est_dd,
                               ntreeTry = 1000,doBest = T,
                               nodesizeTry = c(1:9, seq(10, 1000, by = 10)),)
a$optimal  
if (library("interp", logical.return = TRUE)) {
  
  ## nice little wrapper for plotting results
  plot.tune <- function(o, linear = TRUE) {
    x <- a$results[,1]
    y <- a$results[,2]
    z <- a$results[,3]
    so <- interp(x=x, y=y, z=z, linear = linear)
    idx <- which.min(z)
    x0 <- x[idx]
    y0 <- y[idx]
    filled.contour(x = so$x,
                   y = so$y,
                   z = so$z,
                   xlim = range(so$x, finite = TRUE) + c(-2, 2),
                   ylim = range(so$y, finite = TRUE) + c(-2, 2),
                   color.palette =
                     colorRampPalette(c("yellow", "red")),
                   xlab = "nodesize",
                   ylab = "mtry",
                   main = "error rate for nodesize and mtry",
                   key.title = title(main = "OOB error", cex.main = 1),
                   plot.axes = {axis(1);axis(2);points(x0,y0,pch="x",cex=1,font=2);
                     points(x,y,pch=16,cex=.25)})
  }
  
  ## plot the surface
  plot.tune(a)
  
}
fit <- rfsrc(Surv(OS.time,OS)~.,data = est_dd,
             nodesize = as.numeric(a$optimal[1]),mtry=as.numeric(a$optimal[2]),ntree=1000,##该值建议多调整  
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
save(fit,file ="0603rsf.rda")
rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=predict(fit,newdata = x)$predicted)})
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- 'RSF'

result <- rbind(result,cc)




####Enet ####

library(glmnet)
x1 <- as.matrix(est_dd[,pre_var])
x2 <- as.matrix(Surv(est_dd$OS.time,est_dd$OS))

for (alpha in seq(0,1,0.1)) {
  set.seed(seed)
  fit = cv.glmnet(x1, x2,family = "cox",alpha=alpha)
  save(fit, file=paste0('0322Enet','[α=',alpha,']',".rda"))
  
  rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type='link',newx=as.matrix(x[,-c(1,2)]),s=fit$lambda.min)))})
  
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('Enet','[α=',alpha,']')
  cc$Model <- paste0('Ridge')
  result <- rbind(result,cc)
}
save(result,file="0322result.rda")

####Enet+step ####
library(glmnet)
x1 <- as.matrix(est_dd[,pre_var])
x2 <- as.matrix(Surv(est_dd$OS.time,est_dd$OS))

for (alpha in c(0,1)) {
  set.seed(seed)
  fit = cv.glmnet(x1, x2,family = "cox",alpha=alpha)
  coef=coef(fit, s = fit$lambda.min)
  index=which(coef != 0)
  rid<-row.names(coef)[index]
  est_dd2 <- est_data[,c('OS.time','OS',rid)]
  val_dd_list2 <- lapply(val_dd_list,function(x){x[,c('OS.time','OS',rid)]})
  
  fit <- step(coxph(Surv(OS.time,OS)~.,est_dd2),direction = "both")
  save(fit,file=paste0('0322Enet','[α=',alpha,']',"+step",".rda"))
  rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=predict(fit,type = 'risk',newdata = x))})
  
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('Enet','[α=',alpha,'] + StepCox')
  cc$Model <- paste0('Lasso + StepCox')
  result <- rbind(result,cc)
}



####StepCox####

fit <- step(coxph(Surv(OS.time,OS)~.,est_dd),direction = "both")
save(fit,file="step0311.rda")  
rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=predict(fit,type = 'risk',newdata = x))})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('StepCox')
result <- rbind(result,cc)
}



####CoxBoost####
library(CoxBoost)
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd[,'OS.time'],est_dd[,'OS'],as.matrix(est_dd[,-c(1,2)]),
                            trace=TRUE,start.penalty=500,parallel = T)
cv.res <- cv.CoxBoost(est_dd[,'OS.time'],est_dd[,'OS'],as.matrix(est_dd[,-c(1,2)]),
                      maxstepno=500,type="verweij",penalty=pen$penalty)
fit <- CoxBoost(est_dd[,'OS.time'],est_dd[,'OS'],as.matrix(est_dd[,-c(1,2)]),
                stepno=cv.res$optimal.step,penalty=pen$penalty)
save(fit,file = "coxboost0310.rda")
rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,newdata=x[,-c(1,2)], newtime=x[,1], newstatus=x[,2], type="lp")))})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost')
result <- rbind(result,cc)

####plsRcox####


set.seed(seed)
cv.plsRcox.res=cv.plsRcox(list(x=est_dd[,pre_var],time=est_dd$OS.time,status=est_dd$OS),nt=10,verbose = FALSE)
fit <- plsRcox(est_dd[,pre_var],time=est_dd$OS.time,event=est_dd$OS,nt=as.numeric(cv.plsRcox.res[5]))
save(fit,file="plsrcox0311.rda")
rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type="lp",newdata=x[,-c(1,2)])))})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('plsRcox')
result <- rbind(result,cc)


####superpc####

data <- list(x=t(est_dd[,-c(1,2)]),y=est_dd$OS.time,censoring.status=est_dd$OS,featurenames=colnames(est_dd)[-c(1,2)],nfold=566)
set.seed(seed)
fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) #default
save(fit,file="superpc0311.rda")
cv.fit <- superpc.cv(fit,data,n.threshold = 20,#default 
                     n.fold = 10,
                     n.components=3,
                     min.features=5,
                     max.features=nrow(data$x),
                     compute.fullcv= TRUE,
                     compute.preval=TRUE)
save(cv.fit,file="superpc cv0311.rda")
rs <- lapply(val_dd_list,function(w){
  test <- list(x=t(w[,-c(1,2)]),y=w$OS.time,censoring.status=w$OS,featurenames=colnames(w)[-c(1,2)])
  ff <- superpc.predict(fit,data,test,threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[,1:2],RS=rr)
  return(rr2)
})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('SuperPC')
result <- rbind(result,cc)


####GBM####

set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~.,data = est_dd,distribution = 'coxph',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 6)
# find index for number trees with minimum CV error
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~.,data = est_dd,distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 8)
save(fit,file="gbm0311.rda")
rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,x,n.trees = best,type = 'link')))})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('GBM')
result <- rbind(result,cc)



####survivalsvm####

fit = survivalsvm(Surv(OS.time,OS)~., data= est_dd, gamma.mu = 1)
save(fit,file="survivalsvm0311.rda")
rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit, x)$predicted))})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('survival-SVM')
result <- rbind(result,cc)


result2 <- result
result2$Model <- gsub('α','a',result2$Model)
result2$Model <- gsub('Enet[a=0]','Ridge',result2$Model)
result2$Cindex<-ifelse(result2$Cindex >=0.5,result2$Cindex,1-result2$Cindex)
result$Cindex<-ifelse(result$Cindex >=0.5,result$Cindex,1-result$Cindex)
export(result,"0518cindexintergration.csv")
library(tidyverse)
result$Model <- str_replace(result$Model, "Enet[α=0]", "Ridge") 
result<-import("cindexintergration.csv")


####cindexbar####
library(ggplot2)
library(ggsci)
library(tidyr)
library(ggbreak)
dd4 <- cc%>%
  dplyr::filter(ID!="meta-cohort")%>%
  dplyr::group_by(Model)%>%
  dplyr::summarise(mean(Cindex))


colnames(dd4)[3]<-"Cindex"
library(ggeasy)
install.packages("ggeasy")
library(ggplot2)
p2<-dd4%>%
  ggplot(aes(Cindex,reorder(Model,Cindex)))+
  geom_bar(width = 0.7,stat = 'summary',fun='mean',fill='orange2')+
  theme_classic()+
  labs(y=NULL)+
  geom_text(aes(label=sprintf("%0.3f", Cindex)), nudge_x = -0.07, color="black", size = 3)+
  easy_remove_y_axis()
p2
pdf("0518cindex.pdf",width =2.5 ,height = 14)
p2
dev.off()


dev.off()
colnames(dd4)[2]<-"Cindex"


####cindexheatmap####
dd1 <- pivot_wider(result,names_from = 'ID',values_from = 'Cindex')%>%as.data.frame()
dd2<-dd1
dd2<-dd2%>%column_to_rownames("Model")

dd4<-arrange(dd4,desc(Cindex))
dd2<-dd2[dd4$Model,]%>%select(-`meta-cohort`)
library(pheatmap)
annotation_col = data.frame(dataset=as.character(c("GSE17536","GSE29621","GSE92921","GSE143985","GSE161158")))
row.names(annotation_col) = colnames(dd2)
p1<-pheatmap(dd1,
             
             annotation_col = annotation_col,
             show_colnames =F,
             # # color = colorRampPalette(c("navy", "white", "red"))(50),
             # annotation_colors = ann_colors[1],
             cluster_cols =F,
             cluster_rows = F,
             display_numbers = T,
             fontsize = 10,
             fontsize_row=10,
             fontsize_col=3,
             # breaks=bk
)
p1
pdf("heatmap.pdf",width = 7,height = 14)
p1
dev.off()

####survival analysis#####
library(survminer)
library(survival)
rid<-fit[["var.names"]]
rs1 <- lapply(rs,function(x){cbind(x,
                                   risk=surv_categorize(surv_cutpoint(x, 
                                                                      time="OS.time",
                                                                      event="OS",
                                                                      variables="RS"))[,3],risk2=ifelse(x$RS>median(x$RS),"high","low"))})
fitsurv <- lapply(rs1,function(x)survfit(Surv(OS.time, OS=="1")~ risk2, x))
ggsurvplot(survfit(Surv(OS_MONTHS, OS_STATUS=="1")~ agegroup,tcgaclinicalem),
           data = tcgaclinicalem,
           pval = T,
           conf.int = F, # 置信区间
           risk.table = TRUE, # 风险表
           risk.table.col = "strata",
           
           # legend.labs = c("high", "low"),
           legend.title = "Risk Group",
           size = 1,
           # xlim = c(0,80), # x轴长度，一般为0-10年
           # break.time.by = 20, # x轴步长为20个月
           surv.median.line = "hv", # 限制垂直和水平的中位生存
           ylab = "Survival probability (%)", # 修改y轴标签
           xlab = "Time (Months)",
           # palette = c("#ED4C67","#D980FA"),
           ggtheme = theme_classic())

p1<-ggsurvplot(fitsurv[[10]],
               data = rs1[[10]],
               pval = T,
               conf.int = TRUE, # 置信区间
               risk.table = TRUE, # 风险表
               risk.table.col = "strata",
               
               legend.labs = c("high", "low"),
               legend.title = "Risk Group",
               size = 1,
               # xlim = c(0,80), # x轴长度，一般为0-10年
               # break.time.by = 20, # x轴步长为20个月
               surv.median.line = "hv", # 限制垂直和水平的中位生存
               ylab = "Survival probability (%)", # 修改y轴标签
               xlab = "Time (Months)",
               palette = c("#ed0000", "#00468b"),
               ggtheme = theme_classic()) #

tcgaclinical$groupcom<-paste0(tcgaclinical$risk2,tcgaclinical$group9)
p1
pdf("0330surv161158.pdf",width=4.5,height=4.5)
p1
dev.off()
dev.new()

####gbmgraph####
library(rio)
gbm<-import("gbm.txt")
fix(gbm)
p9<-gbm%>%
  ggplot(aes(`Relative Importance`))+
  geom_bar(width = 0.7,fill='sky')+
  theme_classic()+
  labs(y=NULL)+
  geom_text(aes(label=sprintf("%0.3f", Cindex)), nudge_x = -0.07, color="black", size = 3)+
  easy_remove_y_axis()
p9
pdf("gbm.pdf",width=8,height=15)
a<-summary(fit2)
pdf("gbm.pdf",width=4,height=6)
a<-summary(fit2)
dev.off()

####roc####
library(timeROC)
library(survival)
library(tidyverse)
ROC <- lapply(rs1,function(x)timeROC(T=x$OS.time,   #结局时间
                                     delta=x$OS,   #结局指标
                                     marker=x$RS,   #预测变量
                                     cause=1,   #阳性结局指标数值
                                     weighting="marginal",   #计算方法，默认为marginal
                                     times=c(12, 36, 60),   #时间点，选取1年，3年和5年的生存率
                                     iid=TRUE))
library(grid)
library(reshape2)
library(ggthemes)

rocnumber<-rocnumber%>%rownames_to_column("time")
str(rocnumber)

rocnumber$dataset<-as.factor(rocnumber$dataset)
rocnumber$time<-as.factor(rocnumber$time)
rocnumber<-melt(rocnumber,id.vars="time",variable.name="dataset",value.name="AUC")
p11<-ggplot(rocnumber,aes(dataset,AUC,fill=time))+
  geom_bar(stat="identity",position="dodge")+
  theme_bw()+
  scale_fill_wsj("rgby", "")+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL))+
  ggtitle("AUC")+
  theme(axis.title = element_blank(),legend.position='none')+
  facet_grid(.~time)+
  coord_flip()+
  geom_text(aes(label=sprintf("%0.3f",AUC )), nudge_y = -0.15, color="white", size = 5)
p11
ggsave("aucall.pdf",width=7,height=4)


ROCmerge<-ROC$merge
rocfinal<-data.frame()
rocframe <- data.frame(ROCmerge$TP,ROCmerge$FP)
rocmerge1 <- data.frame(tp=rocframe$t.60,fp=rocframe$t.60.1)%>%
  rownames_to_column('time')
rocmerge1$time <- paste0('5-Year = 0.720')

rocfinal <- rbind(rocfinal,rocmerge1)

# rocfinal$time <- ifelse(rocfinal$time=='t=1',"1-Year = 0.854",
# ifelse(rocfinal$time=='t=3','3-Year = 0.910','5-Year = 0.898'))

ggplot(rocfinal,aes(fp,tp,color=time))+
  geom_line(size=1)+
  labs(x='1-Specificity',y='Sensitivity',color=NULL)+
  theme_bw(base_rect_size = 1.5)+
  geom_abline(slope = 1,color='grey70')+
  theme(panel.grid =element_blank(),
        axis.text = element_text(size=11),
        axis.title = element_text(size=13),
        legend.text = element_text(size=12),
        legend.position = c(0.995,0.012),
        legend.justification = c(1,0))+
  scale_color_nejm()+
  scale_x_continuous(expand = c(0.01,0.01))+
  scale_y_continuous(expand = c(0.01,0.01))
ggsave(filename = 'merge-timeROC.pdf',width = 4.3,height = 4)


rocnumber<-sapply(ROC,function(x)x[["AUC"]])
rocnumber<-rocnumber%>%as.data.frame()
fix(rocnumber)


####cindex all dataset####
std <-sapply(rs1,function(x){cbind(as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1]),
                                   as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[2]))})


a<-summary(coxph(Surv(OS.time,OS)~RS,rs[[1]]))

CC$se<-std[2,]
CC<-cc[c(1,2,3,7,9,10),]
library(ggplot2)
library(ggsci)
p10<-CC%>%
  ggplot(aes(x=Cindex, y=ID)) + 
  geom_bar(position=position_dodge(), stat="identity",fill="#5ab4ac") +
  geom_errorbar(aes(xmin=Cindex-se, xmax=Cindex+se),
                width=.2, position = position_dodge(0.8),                   # Width of the error bars
                size=0.6)+
  theme_bw(base_rect_size = 1.5)+
  theme(axis.title = element_blank(),
        axis.text = element_text(size=12),
        legend.position = 'none',
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5,size=14),
        strip.text = element_text(size=12),
        axis.ticks.x = element_blank())+
  scale_fill_npg()+
  scale_x_continuous(breaks=seq(0,0.96,0.1),expand = c(0,0.01),limits = c(0,0.96))
p10
pdf("cindexdetail.pdf",width=6,height=4)
p10
dev.off()


####compare-factor####
install.packages("compareC")
library(compareC)
tt <- my
dd <- data.frame()
for (i in colnames(my)[5:13]) {
  fit <- summary(coxph(Surv(OS.time,OS)~get(i),tt))
  CC2 <- fit$concordance[1]%>%as.numeric()
  se <- fit$concordance[2]%>%as.numeric()
  p <- compareC(tt$OS.time,tt$OS,tt$RS,tt[,i])$pval
  dd <- rbind(dd,data.frame(ID=i,C=CC,SE=se,P=p))
}
survmerge<-import("~/R/CRC2/clinical_molecular_public_all.txt")
survmerge<-survmerge%>%column_to_rownames("sample")
surv9<-surv9%>%column_to_rownames("accession")
surv9final<-surv9[rownames(gse9),]
survmerge<-survmerge[rownames(merge),]
survmergefinal<-cbind(merge[,1:2],survmerge[,c(2:7,10,11,12,13)],RS=rs1[["meta-cohort"]][["RS"]])



surv2final<-surv2final%>%select(`t stage`,`m stage (0`,gender,`preoperative chemo`,`adjuvant chemo`,`histologic grade`,`n stage`,`ajcc staging`)
for (i in 1:ncol(surv161final)){
  surv161final[,i]<-as.numeric(as.factor(surv161final[,i]))
}


fix(surv1final)
surv1final$grade<-substr(surv1final$grade,1,1)
identical(rownames(rs1[["GSE143985"]]),rownames(surv143final))
surv161final<-cbind(rs1[["GSE161158"]][,1:2],surv161final,RS=rs1[["GSE161158"]][["RS"]])

surv1final=surv1final %>% mutate(gender=case_when(gender == "male" ~ 1,
                                                  gender == "female" ~ 0))

surv143final$msi_status<-as.numeric(as.factor(surv143final$msi_status))
fix(survmergefinal)
str(survmergefinal)

dd <- data.frame()
clincom<-list(survmergefinal,surv1final,surv2final,surv9final,surv143final,surv161final)
dd<-data.frame()
clin<-lapply(clincom,function(x) for (i in colnames(x)[3:ncol(x)]) {
  a<-x%>%drop_na(i)
  fit <- summary(coxph(Surv(OS.time,OS)~ get(i),a))
  CC2 <- fit$concordance[1]%>%as.numeric()
  se <- fit$concordance[2]%>%as.numeric()
  p <- compareC(a$OS.time,a$OS,a$RS,a[,i])$pval
  dd <- rbind(dd,data.frame(ID=i,C=CC2,SE=se,P=p))
}
return(dd))

dd <- data.frame()
for (i in colnames(surv1final)[3:ncol(surv1final)]) {
  a<-surv1final%>%drop_na(i)
  fit <- summary(coxph(Surv(OS.time,OS)~get(i),a))
  CC2 <- fit$concordance[1]%>%as.numeric()
  se <- fit$concordance[2]%>%as.numeric()
  p <- compareC(a$OS.time,a$OS,a$RS,a[,i])$pval
  dd <- rbind(dd,data.frame(ID=i,C=CC2,SE=se,P=p))
}

dd$ll<-paste0(ifelse(dd$P< 0.001,"***",ifelse(dd$P< 0.01,"**",ifelse(dd$P <0.05,"*",""))))

library(RColorBrewer)
ggplot(dd,aes(ID,C,fill=ID))+
  geom_bar(stat='identity',position=position_dodge(0.8),width=0.6)+
  geom_errorbar(aes(ymax=C+SE,ymin=C-SE),
                width=0.1,position = position_dodge(0.8),size=0.6)+
  theme_bw()+
  # brewer.pal(11,"Set3")+
  ggtitle('C-index Comparison')+
  theme(axis.title = element_blank(),
        axis.text = element_text(size=12),
        legend.position = 'none',
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5,size=14),
        strip.text = element_text(size=12),
        axis.ticks.x = element_blank())+
  scale_fill_igv()+
  scale_y_continuous(expand = c(0,0.01),limits = c(0,0.95))+
  geom_text(aes(y=0.90,label=ifelse(P<0.001,"***",ifelse(P<0.01,"**",ifelse(P<0.05,"*","")))),size=5)

ggsave("clinical1.pdf",width=3.816,height=4)
dev.off()
dev.new()

####compare-othermodel####
library(stats4)
library(BiocGenerics)
library(parallel)

library("AnnotationDbi")

library("org.Hs.eg.db")
o<-list(o1,o2,o3,o4,o5)
o<-lapply(o,function(x){cbind(x[1:2],symbol=mapIds(org.Hs.eg.db,keys=x$gene,column="SYMBOL",
                                                   keytype="ENSEMBL",multiVals="first"))})
library(compareC)
###4.2 计算riskScore
FinalGeneExp = sjj$GSE12945[,lassoGene]
myFun = function(x){crossprod(as.numeric(x),o5$coef)}
riskScore = apply(FinalGeneExp,1,myFun)
tcgaclinical1$o5<-riskScore
uni <- data.frame()
for (i in colnames(tcgaclinical1)[c(21,25:29)]){
  scox <- summary(coxph(Surv(OS.time,OS)~get(i),tcgaclinical1))
  p <- compareC(tcgatpm$OS.time,tcgatpm$OS,tcgaclinical1$RS,tcgaclinical1[,i])$pval
  uni <- rbind(uni,data.frame(A=i,
                              HR=scox$conf.int[,1],
                              HR.95L=scox$conf.int[,3],
                              HR.95R=scox$conf.int[,4],
                              pvalue=scox$coefficients[,5],
                              cindex=scox$concordance[1]%>%as.numeric(),
                              cse=scox$concordance[2]%>%as.numeric(),
                              CP=p))
}



trans<-lapply(trans,function(x)scale(x,center=T,scale = T))
trans<-lapply(trans,function(x)x %>% t() %>% as.data.frame())
library(compareC)

FinalGeneExp = lapply(trans1,function(x)x[,lassoGene])
myFun = function(x){crossprod(as.numeric(x),o5$coef)}
riskScore5 = sapply(FinalGeneExp,function(x)(apply(x,1,myFun)))
tcgaclinical1$o5<-riskScore
trans1<-list(trans$GSE17536,trans$GSE29621,trans$GSE92921,trans$GSE143985,trans$GSE161158)
colnames(tcgaclinical1)[25:29]<-c("Fu","Huang","Dekervel","Li","Busuioc")
identical(rownames(trans1[[1]]),rownames(rs$GSE17536))
uni <- data.frame()
for (i in colnames(o29621)[3:8]){
  scox <- summary(coxph(Surv(OS.time,OS)~get(i),o29621))
  p <- compareC(o29621$OS.time,o29621$OS,o29621$RS,o29621[,i])$pval
  uni <- rbind(uni,data.frame(A=i,
                              HR=scox$conf.int[,1],
                              HR.95L=scox$conf.int[,3],
                              HR.95R=scox$conf.int[,4],
                              pvalue=scox$coefficients[,5],
                              cindex=scox$concordance[1]%>%as.numeric(),
                              cse=scox$concordance[2]%>%as.numeric(),
                              CP=p))
}
unique(uni$ID)
dd <- uni
dd$ll <- ifelse(dd$CP<0.001,'***',ifelse(dd$CP<0.01,'**',ifelse(dd$CP<0.05,'*','')))
rownames(dd) <- NULL

ggplot(dd,aes(cindex,reorder(A,cindex)))+
  geom_errorbarh(aes(xmax=cindex+1.5*cse,xmin=cindex-1.5*cse),color="black",height=0,size=0.7)+
  geom_point(size=4,shape=21,fill=pal_nejm()(8)[6])+
  ylab(NULL)+xlab(NULL)+
  labs(title ="TCGA-CRC")+
  geom_vline(xintercept = 0.6,linetype='dashed',size=0.5,color='grey50')+
  theme_bw(base_rect_size = 1)+
  theme(panel.grid =element_blank(),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12),
        axis.title = element_text(size=13),
        plot.title = element_text(hjust = 0.5,size=15),
        legend.position = 'none',
        strip.text = element_text(size=14))+ 
  geom_text(aes(x=0.89,y=A,label=ll),color='black',size=3,vjust=0.76)+
  scale_x_continuous(breaks = c(0.5,0.7,0.9),limits = c(0.4,0.94))
ggsave("0517comparesigTCGA.pdf",width=2.5,height = 4)




library(ggpubr)
library(tidyverse)
library(ConsensusClusterPlus)

tcgatpm<-import("tcga325tpm.txt")
tcgatpm<-tcgatpm %>% column_to_rownames("gene")
tcgatpm<-log2(tcgatpm+1)
celltypeuse<-xCell.data$spill$K
rs<-xCellAnalysis(tcgatpm,parallel.sz=10) #计算


####mutation####
tcga<-merge[which(substr(rownames(merge),1,4)=="TCGA"),]
mergers<-rs1$`meta-cohort`
tcgars<-b[which(substr(b$sample,1,4)=="TCGA"),]
devtools::install_github("RobinHankin/Brobdingnag")

BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)
library(clusterProfiler)
(.packages())
query <- GDCquery(
  project = "TCGA-COAD",
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  access = "open"
)

GDCdownload(query)

GDCprepare(query, save = T,save.filename = "TCGA-COAD_SNP.Rdata")
library(maftools)
library(tidyverse)
library(rio)
library(NMF)
maf <- read.maf(maffilt,clinicalData = tcgaclinical,isTCGA = T)

library(circlize)

col_fun2 <- colorRamp2(
  c(-1.63, 0, 2.62), 
  c("#a99af7", "white", "#eff79a")
)
fabcolors = list(RS = col_fun2,risk=c("#ed0000", "#00468b"))

mafsele<-subsetMaf(maf,tsb=maf@clinical.data$Tumor_Sample_Barcode)

cg_tmp=data[data$Hugo_Symbol %in% genes,]
laml = read.maf(maf = cg_tmp,clinicalData = tcgaclinical,isTCGA = T)
laml@clinical.data$RS<-as.numeric(laml@clinical.data$RS)
lamlsele<-subsetMaf(laml,tsb=laml@clinical.data$Tumor_Sample_Barcode)
lamlsele@clinical.data$RS<-as.numeric(lamlsele@clinical.data$RS)
oncoplot(maf =lamlsele, fontSize = 0.8 ,
         showTumorSampleBarcodes = F ,
         clinicalFeatures = "risk",
         sortByAnnotation = TRUE,
         removeNonMutated = F,
         # genes=gene,keepGeneOrder=TRUE
         top=34)
table(tcgaclinical$risk2)

str(mafsele@clinical.data$RS)
mafsele@clinical.data$RS<-as.numeric(mafsele@clinical.data$RS)
p15
pdf("oncoplot.pdf",width=50,height=5)
p15
dev.off()

a <- mafsele@data %>% 
  .[,c("Hugo_Symbol","Variant_Classification","Tumor_Sample_Barcode")] %>% 
  as.data.frame() %>% 
  mutate(Tumor_Sample_Barcode = substring(.$Tumor_Sample_Barcode,1,12))
gene <- as.character(unique(a$Hugo_Symbol))
sample <- as.character(unique(a$Tumor_Sample_Barcode))

mat <- as.data.frame(matrix("",length(gene),length(sample),
                            dimnames = list(gene,sample)))
mat_0_1 <- as.data.frame(matrix(0,length(gene),length(sample),
                                dimnames = list(gene,sample)))




for (i in 1:nrow(a)){
  mat[as.character(a[i,1]),as.character(a[i,3])] <- as.character(a[i,2])
}
for (i in 1:nrow(a)){
  mat_0_1[as.character(a[i,1]),as.character(a[i,3])] <- 1
}

mat_0_1<-mat_0_1%>%t()%>%as.data.frame()

mat_0_1<-mat_0_1[mafsele@clinical.data$Tumor_Sample_Barcode,]
identical(rownames(mat_0_1),mafsele@clinical.data$Tumor_Sample_Barcode)
mat_0_1$risk2<-mafsele@clinical.data$groupcom

library(tidyverse)
laml.titv = titv(maf = mafsele, plot = FALSE, useSyn = TRUE)
plotTiTv(res = laml.titv)
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
BiocManager::install("GenomeInfoDb")
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Hsapiens.UCSC.hg19)
library(NMF)
detach("package:BSgenome.Hsapiens.UCSC.hg38")

####mutational signature####
luad.tnm <- trinucleotideMatrix(maf=mafsele,  
          
                                ref_genome="BSgenome.Hsapiens.UCSC.hg38")

apobec_enrich <- plotApobecDiff(tnm=luad.tnm, maf=mafsele)
maf.sign = estimateSignatures(mat = luad.tnm, nTry = 6)
plotCophenetic(res = maf.sign)
dev.off()
maf.sig = extractSignatures(mat = luad.tnm, n =3 )
plotSignatures(nmfRes = maf.sig, title_size = 1.3)
laml.og30.cosm = compareSignatures(nmfRes = maf.sig, sig_db = "legacy")

library(pheatmap)
pheatmap::pheatmap(mat=laml.og30.cosm$cosine_similarities, cluster_rows=FALSE, main="cosine similarity against validated signatures")
dev.off()
com.sig<-maf.sig[[2]]
com.sig<-maf.sig[["nmfObj"]]@fit@H
com.sig<-com.sig%>%t()%>%as.data.frame()
maf2 = tmb(maf = mafsele,
           captureSize = 50,
           logScale = TRUE)   
maf2$Tumor_Sample_Barcode<-as.character(maf2$Tumor_Sample_Barcode)
maf2<-inner_join(maf2,tcgaclinical,by=c("Tumor_Sample_Barcode"="Tumor_Sample_Barcode"))



####drugsensitivity####

library(gbm)
rsact <- sapply(act,function(x){RS=as.numeric(predict(fit,x,n.trees = 8341,type = 'link'))})
surv198$rs<-rsact[[1]]
surv620$rs<-rsact[[4]]

colnames(surv198)[7]<-"response1"
library(ggpubr)
library(oncoPredict)
library(data.table)
library(gtools)
library(reshape2)
library(ggpubr)
th=theme(axis.text.x = element_text(angle = 45,vjust = 0.5))

GDSC2_Expr = readRDS('GDSC2_Expr (RMA Normalized and Log Transformed).rds')
GDSC2_Res = readRDS("GDSC2_Res.rds")
GDSC2_Res <- exp(GDSC2_Res) 
library(rio)
tcgatpm<-import("tcga325tpm.txt")
tcgatpm<-tcgatpm%>%column_to_rownames("gene")
testExpr<-as.matrix(gse3)
colnames(testExpr)=paste0('test',colnames(testExpr))
dim(testExpr)
calcPhenotype(trainingExprData = GDSC2_Expr,
              trainingPtype = GDSC2_Res,
              testExprData = testExpr,
              batchCorrect = 'eb',  #   "eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10, 
              printOutput = TRUE, 
              removeLowVaringGenesFrom = 'rawData' )




####clinicalmodel####
library(rms)
library(survminer)
library(survival)
require(survival)
colnames(clinical)[17]<-"OS.time"
colnames(clinical)[18]<-"OS"
RS<-as.numeric(rs$`meta-cohort`$RS)


f <- cph(Surv(OS.time, OS) ~age+gender+msi+pt+RS, x=T, y=T, surv=T, data=ddist1, time.inc=12)

surv <- Survival(f)

nom2 <- nomogram(f, fun=list(function(x) surv(12, x), function(x) surv(36, x), function(x) surv(60, x)), 
                 lp=F, funlabel=c("1-year survival", "3-year survival", "5-year survival"), 
                 maxscale=100, 
                 fun.at=c(0.99, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3,0.2,0.1))
pdf("0515nomo.pdf",width=8,height=5)
plot(nom2)
dev.off()

library(nomogramFormula)
install.packages("nomogramFormula")
results<-formula_rd(nomogram =nom2 )
ddist1<-ddist1 %>% drop_na()

res<-TotalPoints.rms(rd = ddist1,fit = f,nom = nom2)
res<-TotalPoints.rms(rd = clin395,fit = f,nom = nom2)



####calibrate curve####
f1<-cph(Surv(OS.time, OS) ~age+gender+msi+pt+RS,x=T, y=T, surv=T, time.inc=12,data=ddist1)
sum.surv<-summary(f1)
c_index<-sum.surv$concordance
c_inde
cal <- calibrate(f1, cmethod="KM", method="boot", u=12, m=100, B=1000)
f2 <- cph(Surv(OS.time, OS) ~age+gender+msi+pt+RS,x=T, y=T, surv=T, data=ddist1, time.inc=36)
f3 <- cph(Surv(OS.time, OS) ~age+gender+msi+pt+RS, x=T, y=T, surv=T, data=ddist1, time.inc=60)
cal2<- calibrate(f2, cmethod="KM", method="boot", u=36, m=100, B=1000)
cal3<- calibrate(f3, cmethod="KM", method="boot", u=60, m=100, B=1000)

plot(cal,lwd=2,lty=1,errbar.col=c(rgb(0,118,192,maxColorValue=255)),
     xlim=c(0,1),ylim=c(0,1),xlab="Nomogram-Predicted Probability of DFS",
     ylab="Actual DFS(proportion)", col=c(rgb(192,98,83,maxColorValue=255)))
plot(cal3,lwd=2,lty=1,errbar.col=c(rgb(0,118,192,maxColorValue=255)),
     xlim=c(0,1),ylim=c(0,1),
     col=c(rgb(192,98,83,maxColorValue=255)),add=T)
plot(cal2,lwd=2,lty=1,errbar.col=c(rgb(0,118,192,maxColorValue=255)),
     xlim=c(0.2,1),ylim=c(0,1),add=T)
lines(cal[,c("mean.predicted","KM")],type="b",lwd=2,col=c("#6bb952"), pch=16)
lines(cal2[,c("mean.predicted","KM")],type="b",lwd=2,col=c("#ec748b"), pch=16)
lines(cal3[,c("mean.predicted","KM")],type="b",lwd=2,col=c("#c4a751"), pch=16)
abline(0,1,lty=3,lwd=2,col=c(rgb(0,118,192,maxColorValue=255)))
legend("bottomright",
       c("1-year", "3-year", "5-year"),
       col=c("#6bb952", "#ec748b", "#c4a751"),
       lty=1, lwd=2)
####DCA####
library(ggDCA)
library(rms)
library(ezcox)
library(forestmodel)
#2.从github安装ggDCA
devtools::install_github('yikeshu0611/ggDCA')

library(devtools)
# 
# plot1 <- stdca(rt2,outcome="OS_STATUS",ttoutcome="OS_MONTHS",timepoint = 12,
#                predictors="1-Year OS",smooth= T)
plot1<-dca(f)
f0 <- cph(Surv(OS.time, OS) ~RS , x=T, y=T,  data=clin395)
f1<-cph(Surv(OS.time, OS) ~msi+pt+gender+ age+RS , x=T, y=T,  data=clin395)
f2<-cph(Surv(OS.time, OS) ~msi+pt+gender+ age, x=T, y=T,  data=clin395)
f0 <- cph(Surv(OS.time, OS) ~RS , x=T, y=T,  data=ddist1)
f1<-cph(Surv(OS.time, OS) ~msi+pt+gender+ age+RS , x=T, y=T,  data=ddist1)
f2<-cph(Surv(OS.time, OS) ~msi+pt+gender+ age, x=T, y=T,  data=ddist1)
plot1<-dca(f0,f1,f2)
plot2<-dca(f0)  
ggplot(plot1,lty=1)

c$score<-f0$linear.predictors
c$scorewhole<-f1$linear.predictors


survtcga<-survtcga%>%column_to_rownames("sample")
survtcga1<-survtcga[rownames(dat),]
dat<-cbind(dat,survtcga1$AGE,survtcga1$SEX,survtcga1$AJCC_PATHOLOGIC_TUMOR_STAGE)
fix(dat)
dat1$Stage <- str_replace(dat1$Stage, "TAGE", "tage") 
dat1<-dat
dat1$STAGE[which(dat1$STAGE=="")]=NA
dat1<-dat1%>%drop_na()
dat1$AGE<-as.numeric(dat1$AGE) 
fix(dat1)
coxphmodel1 <- coxph(Surv(OS_MONTHS, OS_STATUS) ~ Age+Sex+Stage+`Risk Score`, dat1)
forest_model(coxphmodel1,
             
             format_options = forest_model_format_options(colour = "#00A896",  # 修改颜色       
                                                          shape = 15,            # 改变形状
                                                          text_size = 5,       # 字体文本的大小
                                                          point_size = 3,        # 森林图中“方框”的大小
                                                          banded = TRUE),
             factor_separate_line = TRUE) 