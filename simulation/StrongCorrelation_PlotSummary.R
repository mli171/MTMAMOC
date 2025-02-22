datafile.dir = "~/Desktop/MTM/MTM_CPT/MTM_CPT/power/simres/strpos/"

nSim = 1000

##### parameters
allDelta = c(-0.3, 0.3, 0.5)
# allprop  = c(0.1, 0.3, 0.5, 0.7, 0.9)
allprop = c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
# allnY    = c(3, 5, 10)
allnY    = 3
# number of DeltaT = 0
allm = c(1, 2, 3, 4)
eg = expand.grid(Delta=allDelta, prop=allprop, nY=allnY, m=allm, stringsAsFactors=FALSE)

allm = c(5, 6)
eg = rbind(eg,
           expand.grid(Delta=allDelta, prop=allprop, nY=allnY, m=allm, 
                       stringsAsFactors=FALSE))
eg$ID = 1:NROW(eg)


res = data.frame(matrix(NA, nrow=NROW(eg)*3, ncol=2))
colnames(res) = c("Power", "Methods")

for(iscen in 1:NROW(eg)){
  
  condition = eg[iscen,]
  load(paste0(datafile.dir, "multicp_configID", condition$ID, "_strpos_model1.RData"))
  
  if(NROW(FinalRes) > 1000){FinalRes = FinalRes[sample(1:NROW(FinalRes), 1000),]}
  
  DeltaT = condition$Delta
  propT  = condition$prop
  nY     = condition$nY
  N      = 365*nY
  tauT   = round(propT*N)
  
  # p1 = which(FinalRes$CpValY.cusum > 0.925)
  # # [1]  7  8 13 14 32 33 37
  # p2 = which(FinalRes$CpValZ.cusum > 0.925)
  # # [1]  3  7 14 32 33 37
  # p3 = which(FinalRes$fit > 19.836)
  # [1] 14 19 37
  
  # FinalRes$CpLocY.cusum[p1]
  # # [1] 513 529 580 549 929 538 520
  # FinalRes$CpLocY.cusum[p2]
  # # [1] 274 513 549 929 538 520
  # FinalRes$tauhat[p3]
  # # [1] 536 657 534
  
  
  res[iscen, ] = c(mean(1*(FinalRes$CpValY.cusum > 1.357)), "CUSUMY")
  res[NROW(eg)+iscen, ] = c(mean(1*(FinalRes$CpValZ.cusum > 1.357)), "CUSUMZ")
  res[NROW(eg)*2+iscen, ] = c(mean(1*(FinalRes$fit.mtm > 17.918), na.rm = TRUE), "MTM")
  
}

res = cbind(do.call("rbind", rep(list(eg), 3)), res)

res$m = factor(res$m, levels = c(1, 6, 2, 5, 3, 4),
               labels = c('Delta[1] *"=" *Delta[2] *"=" *Delta[3] *"=" *Delta[4] *"=" *Delta', 
                          'Delta[3] *"=0; " *Delta[1] *"=" *Delta[2] *"=" *Delta[4] *"=" *Delta',
                          'Delta[1] *"=0; " *Delta[2] *"=" *Delta[3] *"=" *Delta[4] *"=" *Delta', 
                          'Delta[1] *"=" *Delta[3] *"=0; " *Delta[2] *"=" *Delta[4] *"=" *Delta',
                          'Delta[1] *"=" *Delta[2] *"=0; " *Delta[3] *"=" * Delta[4] *"=" *Delta',
                          'Delta[1] *"=" *Delta[2] *"=" *Delta[3] *"=0; " *Delta[4] *"=" *Delta'))
res$Delta = factor(res$Delta, labels = c('Delta*" = -0.3"', 'Delta*" = 0.3"', 'Delta*" = 0.5"'))
res$Methods = factor(res$Methods, labels = c('M[I]', 'tilde(M)[epsilon]', 'Lambda[max]'))
res$Power = as.numeric(res$Power)

library(ggplot2)
library(scales)
p = ggplot(data=res, aes(x=prop, y=Power, color=Methods))
p = p + geom_line(aes(color=Methods), linewidth=1)
p = p + geom_point(aes(shape=Methods, color=Methods), size = 5)
p = p + geom_line(y = 0.05, linetype=2, col = "grey75", linewidth=1)
p = p + ylim(0,1) + xlim(0.2,0.8)
p = p + scale_shape_manual(values=c(2,6,0,17), labels = parse_format())
p = p + scale_color_manual(values=c('red', 'blue', 'orange'), labels = parse_format())
p = p + xlab(expression(tau/n)) + ylab("Power")
p = p + facet_grid(Delta ~ m, labeller = label_parsed)
p = p + theme(axis.text.x = element_text(size=20,color="black"),
              axis.text.y = element_text(size=20,color='black'),
              axis.title.x = element_text(size=30,vjust=1),
              axis.title.y = element_text(size=25,vjust=1),
              strip.text.x = element_text(size=20,,vjust=1),
              strip.text.y = element_text(size=20,,vjust=1),
              strip.background = element_blank())
p = p + theme(legend.key.size = unit(0.3, "in"),
              legend.key.width = unit(0.8, "in"),
              legend.text = element_text(size = 17),
              legend.title= element_blank(),
              panel.spacing = unit(1, "lines"),
              panel.border = element_rect(colour = "black", fill=NA),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              legend.position = "inside", 
              legend.position.inside =  c(0.06, 0.93),
              legend.key = element_rect(colour = NA, fill = NA),
              legend.direction='vertical',
              legend.box='vertical',
              legend.background = element_rect(fill='transparent'),
              legend.box.background = element_blank())
# p = p + guides(color = guide_legend(ncol = 2),
#                shape = guide_legend(ncol = 2),
#                linetype = guide_legend(ncol = 2))
p

ggsave(paste0("~/Desktop/MTM/MTM_CPT/MTM_CPT/power/simres/StrPosComparison_model1.eps"),
       plot = p, width = 18, height = 9, unit = 'in', dpi = 300)





















datafile.dir = "~/Desktop/MTM/MTM_CPT/MTM_CPT/power/simres/strpos/"

nSim = 1000

##### parameters
allDelta = c(-0.3, 0.3, 0.5)
# allprop  = c(0.1, 0.3, 0.5, 0.7, 0.9)
allprop = c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
# allnY    = c(3, 5, 10)
allnY    = 3
# number of DeltaT = 0
allm = c(1, 2, 3, 4)
eg = expand.grid(Delta=allDelta, prop=allprop, nY=allnY, m=allm, stringsAsFactors=FALSE)

allm = c(5, 6)
eg = rbind(eg,
           expand.grid(Delta=allDelta, prop=allprop, nY=allnY, m=allm, 
                       stringsAsFactors=FALSE))
eg$ID = 1:NROW(eg)


res = data.frame(matrix(NA, nrow=NROW(eg)*3, ncol=2))
colnames(res) = c("Power", "Methods")

for(iscen in 1:NROW(eg)){
  
  condition = eg[iscen,]
  load(paste0(datafile.dir, "multicp_configID", condition$ID, "_strpos_model2.RData"))
  
  if(NROW(FinalRes) > 1000){FinalRes = FinalRes[sample(1:NROW(FinalRes), 1000),]}
  
  DeltaT = condition$Delta
  propT  = condition$prop
  nY     = condition$nY
  N      = 365*nY
  tauT   = round(propT*N)
  
  # p1 = which(FinalRes$CpValY.cusum > 0.925)
  # # [1]  7  8 13 14 32 33 37
  # p2 = which(FinalRes$CpValZ.cusum > 0.925)
  # # [1]  3  7 14 32 33 37
  # p3 = which(FinalRes$fit > 19.836)
  # [1] 14 19 37
  
  # FinalRes$CpLocY.cusum[p1]
  # # [1] 513 529 580 549 929 538 520
  # FinalRes$CpLocY.cusum[p2]
  # # [1] 274 513 549 929 538 520
  # FinalRes$tauhat[p3]
  # # [1] 536 657 534
  
  
  res[iscen, ] = c(mean(1*(FinalRes$CpValY.cusum > 0.923)), "CUSUMY")
  res[NROW(eg)+iscen, ] = c(mean(1*(FinalRes$CpValZ.cusum > 0.923)), "CUSUMZ")
  res[NROW(eg)*2+iscen, ] = c(mean(1*(FinalRes$fit.mtm > 19.836), na.rm = TRUE), "MTM")
  
}

res = cbind(do.call("rbind", rep(list(eg), 3)), res)

res$m = factor(res$m, levels = c(1, 6, 2, 5, 3, 4),
               labels = c('Delta[1] *"=" *Delta[2] *"=" *Delta[3] *"=" *Delta[4] *"=" *Delta', 
                          'Delta[3] *"=0; " *Delta[1] *"=" *Delta[2] *"=" *Delta[4] *"=" *Delta',
                          'Delta[1] *"=0; " *Delta[2] *"=" *Delta[3] *"=" *Delta[4] *"=" *Delta', 
                          'Delta[1] *"=" *Delta[3] *"=0; " *Delta[2] *"=" *Delta[4] *"=" *Delta',
                          'Delta[1] *"=" *Delta[2] *"=0; " *Delta[3] *"=" * Delta[4] *"=" *Delta',
                          'Delta[1] *"=" *Delta[2] *"=" *Delta[3] *"=0; " *Delta[4] *"=" *Delta'))
res$Delta = factor(res$Delta, labels = c('Delta*" = -0.3"', 'Delta*" = 0.3"', 'Delta*" = 0.5"'))
res$Methods = factor(res$Methods, labels = c('M[I]', 'tilde(M)[epsilon]', 'Lambda[max]'))
res$Power = as.numeric(res$Power)

library(ggplot2)
library(scales)
p = ggplot(data=res, aes(x=prop, y=Power, color=Methods))
p = p + geom_line(aes(color=Methods), linewidth=1)
p = p + geom_point(aes(shape=Methods, color=Methods), size = 5)
p = p + geom_line(y = 0.05, linetype=2, col = "grey75", linewidth=1)
p = p + ylim(0,1) + xlim(0.2,0.8)
p = p + scale_shape_manual(values=c(2,6,0,17), labels = parse_format())
p = p + scale_color_manual(values=c('red', 'blue', 'orange'), labels = parse_format())
p = p + xlab(expression(tau/n)) + ylab("Power")
p = p + facet_grid(Delta ~ m, labeller = label_parsed)
p = p + theme(axis.text.x = element_text(size=20,color="black"),
              axis.text.y = element_text(size=20,color='black'),
              axis.title.x = element_text(size=30,vjust=1),
              axis.title.y = element_text(size=25,vjust=1),
              strip.text.x = element_text(size=20,,vjust=1),
              strip.text.y = element_text(size=20,,vjust=1),
              strip.background = element_blank())
p = p + theme(legend.key.size = unit(0.3, "in"),
              legend.key.width = unit(0.8, "in"),
              legend.text = element_text(size = 17),
              legend.title= element_blank(),
              panel.spacing = unit(1, "lines"),
              panel.border = element_rect(colour = "black", fill=NA),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              legend.position = "inside", 
              legend.position.inside =  c(0.06, 0.93),
              legend.key = element_rect(colour = NA, fill = NA),
              legend.direction='vertical',
              legend.box='vertical',
              legend.background = element_rect(fill='transparent'),
              legend.box.background = element_blank())
# p = p + guides(color = guide_legend(ncol = 2),
#                shape = guide_legend(ncol = 2),
#                linetype = guide_legend(ncol = 2))
p

ggsave(paste0("~/Desktop/MTM/MTM_CPT/MTM_CPT/power/simres/StrPosComparison_model2.eps"), 
       plot = p, width = 18, height = 9, unit = 'in', dpi = 300)