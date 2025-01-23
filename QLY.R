#Threshold identification

X <- c("Relative_Coverage_Toxic_Plants", "Relative_Quantity_Toxic_Plants", "Relative_AGB_Toxic_Plants", "Pi_Toxic_Plants", "Relative_High_Toxic_Plants")
Y <- c("pH_0_30", "BD_0_30", "WHC_0_30", "TC_Soil_0_30", "TN_Soil_0_30", "TP_Soil_0_30",
       "Simpson", "Shannon_wiener", "Pielou", "Richness", "AGB", "Cover_Plant", "Number", "TC_Plant", "TN_Plant", "TP_Plant",
       "BGB_0_30", "TC_Root_0_30", "TN_Root_0_30", "TP_Root_0_30")

# 读取数据
dat = read.csv('', header = T, sep = ',',check.names = FALSE)

# compute_pairs函数
compute_pairs =  function(X, Y, df, log.y){
  OUT = list()
  for (i in Y){
    for (j in X){
      print(c(i,j,date()))
      dd = data.frame(y = if(log.y == TRUE) {log10(df[,i]*100+1)} else {df[,i]}, x = df[,j], sampleID = df$Order_ID)
      dd = na.omit(dd)
      res = compute(dd)
      OUT$outliers = rbind(data.frame(Yvar = i, Xvar = j, res$outliers), OUT$outliers)
      OUT$gam_test = rbind(data.frame(Yvar = i, Xvar = j, res$gam_test), OUT$gam_test)
      
      OUT$aics = rbind(data.frame(Yvar = i, Xvar = j, res$aics), OUT$aics)
    }
  }
  return(OUT)
}
# compute函数
compute=function(dd){
  OUT=list()
  OUT$outliers=dd
  #check gam and NL
  result=lapply(1:1000, function(bt,dat){
    ddi=sample(1:nrow(dat),size = nrow(dat),replace = TRUE)
    ddi=dat[ddi,]
    mdl0=glm(data=ddi,formula=y~x,family = "gaussian")
    mdl1=gam(data=ddi,formula=y~s(x),family = "gaussian")
    mdl2=glm(data = ddi,formula = y~poly(x,degree = 2),family = "gaussian")
    aics=c(aic.linear=mdl0$aic,aic.gam=mdl1$aic,aic.Quad=mdl2$aic) 
    aics=aics[order(aics)] 
    daics=aics-min(aics) 
    complexity=c("aic.linear","aic.Quad","aic.gam") 
    daics[daics<2]=0
    bf=complexity[min(match(names(daics[daics==0]),complexity))] 
    bf=strsplit(bf,".",fixed = TRUE)[[1]][2] 
    
    out2 = cbind(data.frame(t(aics)), bf)
    
    return(out2)
    
  },dat=dd)
  
  OUT$gam_test = table(sapply(result, '[[','bf'))
  OUT$aics = data.frame(aic.gam = sapply(result, '[[', 'aic.gam'),
                        aic.linear = sapply(result, '[[', 'aic.linear'),
                        aic.Quad = sapply(result, '[[', 'aic.Quad'))
  return(OUT)
  
}

result <- compute_pairs(X, Y, dat, log.y = TRUE) 

gam_test = result$gam_test
gam_test$Var1 = factor(gam_test$Var1, levels = c('gam','Quad', 'lineal'))
gam_test$Xvar = gsub('AridityIndex','Wat bal', gam_test$Xvar)
head(gam_test)
gam_test_plot = ggplot(subset(gam_test,Xvar != 'Coverage'), aes(Xvar, Freq/1000, fill = Var1)) +
  geom_bar(stat = 'identity') +
  geom_hline(yintercept = 0.5, lwd = 1.5, color = 'black') +
  facet_wrap(.~Yvar)+
  labs(x = NULL, y = 'Frequency', fill = 'fitting model') +
  scale_fill_manual(values = c("blue", "red", "green"),  
                    labels = c("gam", "quad", "linear")) +  
  theme_classic()+
  scale_y_continuous(labels = label_percent()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,size = 14, face = "bold"),
        legend.direction = 'horizontal',
        legend.position = 'top',
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        strip.text = element_text(size = 16, face = "bold"))


Threshold_pairs = function(X, Y, df , log.y){
  OUT = list()
  for (i in Y){
    for (j in X){
      print(c(i,j,date()))
      dd = data.frame(y = if(log.y == TRUE) {log10(df[,i]*100+1)} else {df[,i]}, x = df[,j], sampleID = df$Order_ID)
      dd = na.omit(dd) 
      
      
      mdlseg=chngptm(formula.1 = y~1, formula.2 = ~x, dd, save.boot = TRUE, type="segmented", family = "gaussian") 
      mdlstep=chngptm(formula.1 = y~1,formula.2 = ~x, dd, save.boot = TRUE, type="step", family = "gaussian")
      mdlsteg=chngptm(formula.1 = y~1,formula.2 = ~x, dd, save.boot = TRUE, type="stegmented", family = "gaussian")
      
      
      Threshold_seg= mdlseg$chngpt
      Threshold_step= mdlstep$chngpt
      Threshold_steg=mdlsteg$chngpt
      
      
      AIC_seg=mdlseg$best.fit$aic
      AIC_step=mdlstep$best.fit$aic
      AIC_steg=mdlsteg$best.fit$aic
      
      
      results <- data.frame(
        Model = c("Segmented", "Step", "Stegmented"),
        Threshold = c(Threshold_seg, Threshold_step, Threshold_steg),
        AIC = c(AIC_seg, AIC_step, AIC_steg)
      )
      
      OUT$outliers = rbind(data.frame(Yvar = i, Xvar = j, results), OUT$outliers) ####### 它将一个新的数据框添加到 `outliers` 数据框的顶部，以扩展 `outliers` 数据框的行数
      
    }
  }
  return(OUT) 
}


result <-Threshold_pairs(X, Y, dat, log.y = TRUE)  

write.csv(result$outliers)

#Draw an inflection point graph
TP_Plant_data <- dat %>%
  select(Relative_High_Toxic_Plants, TP_Plant) %>%
  mutate(TP_Plant_log = log10(TP_Plant *100 + 1))

TP_Plant_data$group=1

TP_Plant_data$group[TP_Plant_data$Relative_High_Toxic_Plants > 0.32] <- 2

TP_Plant_data$group <- factor(TP_Plant_data$group)

TP_Plant_plot有标签 <- ggplot(TP_Plant_data, aes(x = Relative_High_Toxic_Plants, y = TP_Plant_log)) +
  geom_point(col = "grey", pch = 21, size = 3, fill = "grey") +
  geom_smooth(data = subset(TP_Plant_data, group == "1"), aes(group = group), method = "lm",
              formula = y ~ x, size = 2.5, color = 'blue', se = FALSE) +
  geom_smooth(data = subset(TP_Plant_data, group == "2"), aes(group = group), method = "lm",
              formula = y ~ x, size = 2.5, color = 'red', se = FALSE) +
  geom_smooth(method = "gam", formula = y ~ s(x), color = "black", linetype = "dashed", size = 2.5, se = FALSE) +
  geom_vline(xintercept = 0.32, lty = 2, color = 'black', lwd = 2) +
  xlab("Relative_High_Toxic_Plants") + ylab("") +  
  scale_y_continuous(limits = c(0, NA), labels = number_format(accuracy = 0.01)) +
  theme_classic() +
  theme(axis.text.x = element_text(family = "serif", size = 40, face = "bold"),
        axis.text.y = element_text(family = "serif", size = 40, face = "bold"),
        axis.title.x = element_text(family = "serif", size = 40, face = "bold"),
        axis.title.y = element_blank(),  
        axis.line = element_line(size = 2),
        axis.ticks.length = unit(0.5, "cm"),
        axis.ticks = element_line(size = 2))

#Stability test
data <- read.csv("qly9.25.csv")
data <- data[, c("BGB_0_30", "Relative_Coverage_Toxic_Plants")]
data$BGB_0_30 <- log10(data$BGB_0_30*100+1)

bootthres = function(data, indices, formula = y ~ x, thres = thres) {
  d = data[indices, ]  # 允许自助法选择样本
  mdl = lm(data = d, formula = BGB_0_30 ~ Relative_Coverage_Toxic_Plants)
  slp = coef(mdl)[2]  # 提取斜率
  intcp = coef(mdl)[1]  # 提取截距
  thy = predict.lm(mdl, newdata = data.frame(Relative_Coverage_Toxic_Plants = thres))  # 在阈值处预测
  return(c(slp, intcp, thy))
}

funcdiff = function(data, thres, bootthres = bootthres) {
  dfs = list(before = data[data$Relative_Coverage_Toxic_Plants <= thres, ],
             after = data[data$Relative_Coverage_Toxic_Plants > thres, ])
  
  mdls = lapply(dfs, function(x) return(boot(data = x, statistic = bootthres, R = 200, thres = thres)))
  
  lboot = lapply(mdls, function(x) {
    res = as.data.frame(x$t)
    colnames(res) = c("slope", "intcp", "value")
    return(res)
  })
  
  res = data.frame(thres = rep(thres, 400),
                   position = c(rep("before", 200), rep("after", 200)))
  res$slope = c(lboot$before$slope, lboot$after$slope)
  res$intcp = c(lboot$before$intcp, lboot$after$intcp)
  res$value = c(lboot$before$value, lboot$after$value)
  return(res)
}

resMbio_BGB = funcdiff(data, 0.10, bootthres = bootthres)
write.csv(resMbio_BGB, "BGB.csv", quote = F)

#Random Forest
pH_0_30 <- dat %>%
  select(pH_0_30,
         Relative_Coverage_Toxic_Plants, Relative_Quantity_Toxic_Plants, 
         Relative_AGB_Toxic_Plants,
         Pi_Toxic_Plants, Relative_High_Toxic_Plants) %>%
  mutate(across(c(pH_0_30), 
                ~ log10(. * 100 + 1)))


set.seed(123)

pH_0_30_PE6_rf <- randomForest(pH_0_30 ~ ., data = pH_0_30, importance = TRUE, proximity = TRUE)

pH_0_30_PE6_rf

set.seed(123)
pH_0_30_PE6_perm <- rf.significance(pH_0_30_PE6_rf, pH_0_30[,-1], nperm = 99, ntree = 501)

pH_0_30_PE6_perm

set.seed(123)
pH_0_30_PE6_rfP <- rfPermute(pH_0_30 ~ ., data = pH_0_30, ntree = 500, na.action = na.omit, nrep = 100, num.cores = 1)

pH_0_30_PE6_rfP

pH_0_30_PE6_dat <- importance(pH_0_30_PE6_rfP, sort.by = NULL, decreasing = TRUE)

pH_0_30_PE6_dat
write.csv(pH_0_30_PE6_dat, file = "pH_0_30_PE6_dat.csv", row.names = FALSE)

custom_colors <- c("Sig" = "#0072B2",    
                   "In_sig" = "#D55E00")  


pH_0_30_p <- pH_0_30_PE6_dat %>% 
  as_tibble(rownames = "names") %>% 
  data.frame() %>% 
  mutate(label = if_else(X.IncMSE.pval < 0.001, "***", 
                         if_else(X.IncMSE.pval < 0.01, "**", 
                                 if_else(X.IncMSE.pval < 0.05, "*", "ns"))),
         X.IncMSE = as.numeric(X.IncMSE)) %>% 
  arrange(X.IncMSE) %>% 
  mutate(group = if_else(label == "ns", "In_sig", "Sig"),
         names = forcats::fct_inorder(names)) %>%
  ggplot(aes(x = names, y = X.IncMSE)) + 
  geom_bar(aes(fill = group), stat = "identity") + 
  geom_text(data = . %>% filter(label != "ns"),  
            aes(y = X.IncMSE + 1, label = label), 
            size = 30) +  
  labs(x = "", y = "%IncMSE") + 
  coord_flip() + 
  scale_fill_manual(values = custom_colors) +  
  theme(
    axis.text.x = element_text(size = 60, face = "bold"), 
    axis.text.y = element_text(size = 20, face = "bold"), 
    axis.title.x = element_text(size = 40, face = "bold"),  
    axis.title.y = element_text(size = 40, face = "bold"),  
    panel.background = element_blank(),           
    axis.line = element_line(color = "black", size = 1),     
    axis.ticks = element_line(color = "black", size = 1),  
    axis.ticks.length = unit(0.5, "cm"),        
    plot.title = element_text(size = 20, face = "bold"),  
    legend.position = "none"  
  )

#pearson correlation coefficient
dat = read.csv("qly.csv", header = T, sep = ',', check.names = FALSE)

M2 = cor(dat)
M2

write.csv(M2, "correlation_matrix.csv")

jpeg(filename = "correlation_heatmap.jpeg", width = 1000, height = 1000, res = 100)

corrplot.mixed(M2, 
               lower = 'shade',  
               upper = 'pie',    
               order = 'hclust',
               tl.pos = 'lt',    
               tl.col = "black")  

dev.off()

#Hierarchical partitioning and multiple linear regression

library(ggplot2)
library(MuMIn)
library(performance)
library(rdacca.hp) 
library(patchwork) 
library(RColorBrewer)

datatotal = read.csv("qly.csv", header = T, sep = ',', check.names = FALSE)

#提取目标列
data <- datatotal[, c("pH",
                      "BD",  
                      "SR",
                      "Precipitation",
                      "Temperature",
                      "SWC",
                      "EMF2")]
print(data)

data[,c(1,2,3,4,6)]<-log(data[,c(1,2,3,4,6)]+1)
data[,5]<-log(data[,5]-min(data[,5])+1)
data[,7]<-log(data[,7]-min(data[,7])+1)

data$pH<-(data$pH-mean(data$pH))/sd(data$pH)
data$BD<-(data$BD-mean(data$BD))/sd(data$BD)
data$Precipitation<-(data$Precipitation-mean(data$Precipitation))/sd(data$Precipitation)
data$SR<-(data$SR-mean(data$SR))/sd(data$SR)
data$Temperature<-(data$Temperature-mean(data$Temperature))/sd(data$Temperature)
data$EMF2<-(data$EMF2-mean(data$EMF2))/sd(data$EMF2)
data$SWC<-(data$SWC-mean(data$SWC))/sd(data$SWC)

fit<-lm(EMF2 ~ pH + BD + Precipitation + SR + Temperature  +SWC, data=data)

summary(fit)

r2(fit)

stat.lm <- summary(fit)$coefficients[-1, ]
stat.lm <- data.frame(stat.lm, check.names = FALSE)
stat.lm$env <- rownames(stat.lm)
stat.lm$sig <- ifelse(stat.lm$'Pr(>|t|)'>0.05, '', 
                      ifelse(stat.lm$'Pr(>|t|)'>0.01, '*', 
                             ifelse(stat.lm$'Pr(>|t|)'>0.001, '**', '***')))
stat.lm$label <- paste(stat.lm$env, stat.lm$sig)

stat.lm <- stat.lm[order(stat.lm$Estimate), ]
stat.lm$env <- factor(stat.lm$env, levels = stat.lm$env)

stat.lm

p1 <- ggplot(stat.lm, aes(x = env, y = Estimate, color = env)) +  
  geom_point(size = 10) +  
  geom_errorbar(aes(ymin = Estimate - `Std. Error`, ymax = Estimate + `Std. Error`), width = 0.3, size = 1.5) +  #调整误差帮
  coord_flip() +  
  geom_hline(yintercept = 0, linetype = 2, size = 3) +  
  labs(x = '', y = 'Parameter estimate') +  
  theme(panel.grid = element_blank(), 
        panel.background = element_blank(), 
        axis.line.x = element_line(size = 3),  
        axis.ticks.x = element_line(size = 3),  
        axis.ticks.length = unit(0.4, "cm"),  
        axis.ticks.y = element_blank(), 
        legend.key = element_blank(),
        text = element_text(size = 45, color = "black", face = "bold", family = "Times New Roman"),  
        axis.text.x = element_text(size = 45, color = "black", family = "Times New Roman"),  
        axis.text.y = element_text(size = 45, color = "black", family = "Times New Roman"),  
        legend.text = element_text(size = 45, color = "black", family = "Times New Roman")) +  
  scale_x_discrete(breaks = stat.lm$env, labels = stat.lm$label, position = 'top') +
  scale_y_continuous( expand = c(0, 0),limits = c(-0.8, 0.8), breaks = seq(-0.8, 0.8, by = 0.2))+  
  scale_color_manual(values = c("SWC" = "#F0A73A", 
                                "Temperature" = "#80A6E2", 
                                "SR" = "#F46F43", 
                                "BD" = "red", 
                                "pH" = "#6F6F6F", 
                                "Precipitation" = "#8DD3C7"))  

print(p1)

iv <- list(
  Precipitation = data[c('Precipitation')],
  pH = data[c('pH')],
  BD= data[c('BD')],
  SR = data[c('SR')],
  Temperature = data[c('Temperature')], 
  SWC = data[c('SWC')]
)


hp <- rdacca.hp(data['EMF2'], iv, method = 'RDA', type = 'adjR2')
hp  

hp$Total_explained_variation

hp.ie <- data.frame(hp$Hier.part, check.names = FALSE)
hp.ie$env <- rownames(hp.ie)
hp.ie$env <- factor(hp.ie$env, levels = rev(hp.ie$env))
hp.ie$exp_var <- ''

colnames(hp.ie)[colnames(hp.ie) == "I.perc(%)"] <- "I.perc"

hp.ie

hp.ie$Individual[hp.ie$Individual < 0] <- 0

head(hp.ie)

total_individual <- sum(hp.ie$Individual, na.rm = TRUE)  

hp.ie$I.perc <- (hp.ie$Individual / total_individual) * 100

head(hp.ie)

write.csv(hp.ie, file = "hp_ie低.csv", row.names = FALSE)

p2 <- ggplot(hp.ie, aes(x = exp_var, y = I.perc, fill = env)) + 
  geom_col(width = 0.6) + 
  scale_fill_manual(values = c(
    "SWC" = "#F0A73A", 
    "Temperature" = "#80A6E2", 
    "SR" = "#F46F43", 
    "BD" = "red", 
    "pH" = "#6F6F6F", 
    "Precipitation" = "#8DD3C7"
  )) +  
  theme(panel.grid = element_blank(), 
        panel.background = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.ticks.y = element_line(color = 'black', size = 3), 
        axis.line.y = element_line(color = 'black', size = 3), 
        axis.text.y = element_text(size = 45, color = 'black', face = 'bold'),  
        axis.title.y = element_text(size = 45, color = 'black', face = 'bold'),  
        axis.title.x = element_text(size = 45, face = "bold", color = 'black'), 
        legend.position = 'none',
        axis.ticks.length.y = unit(0.4, "cm")) +  
  labs(x = expression(R^2 == 0.46),  
       y = 'Relative Effect of Estimates (%)', 
       fill = '') +
  scale_x_discrete(position = 'top') +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) 

print(p2)

library(patchwork)

p2 <- p2 + theme(legend.position = 'none')
p2 + p1 + plot_layout(ncol = 2, widths = c(1, 2.5))

library(ggplot2)
library(patchwork)

combined_plot <- p2 + p1 + plot_layout(ncol = 2, widths = c(1, 2.5))
combined_plot
ggsave("1.jpeg", plot = combined_plot, 
       width = 25, height = 16, 
       units = "in", dpi = 600) 



