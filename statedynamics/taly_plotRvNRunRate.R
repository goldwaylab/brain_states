args <- commandArgs(TRUE)
name_root <- ''
numClusters <- 6
basedir <- 'C:/Users/marko/ketamine/brain_states/brain_states-master_1/'
receptor = 'dopamine_D2'
resdir <- switch(
  receptor,
  'Raw' = 'C:/Users/marko/ketamine/Results_No_Physio/Raw/results/',
  'dopamine_D2' = 'C:/Users/marko/ketamine/Results_No_Physio/dopamine_D2/results/',
  'opioid_mu_1' = 'C:/Users/marko/ketamine/Results_No_Physio/opioid_mu_1/results/',
  'opioid_kappa_1' = 'C:/Users/marko/ketamine/Results_No_Physio/opioid_kappa_1/results/',
  'NMDA' = 'C:/Users/marko/ketamine/Results_No_Physio/NMDA/results/',
)


library(ggplot2)
library(R.matlab)
library(RColorBrewer)

masterdir <- resdir

source(paste(basedir,'code/plottingfxns/GeomSplitViolin.R',sep=''))
source(paste(basedir,'code/plottingfxns/plottingfxns.R',sep=''))

restRunRate <- readMat(paste(masterdir,'analyses/transitionprobabilities/KetDwellTime_k',numClusters,name_root,'.mat',sep = ''))$RunRate
nbackRunRate <- readMat(paste(masterdir,'analyses/transitionprobabilities/No_KetDwellTime_k',numClusters,name_root,'.mat',sep = ''))$RunRate
clusterNames <- readMat(paste(masterdir,name_root,'/clusterAssignments/k',numClusters,name_root,'.mat',sep=''))
clusterNames <- unlist(clusterNames$clusterAssignments[[1]][[5]])
clusterColors <- getClusterColors(numClusters)
RNcolors <- c('#005C9F','#FF8400')  

grps <- rbind(matrix('Ket',nrow = nrow(restRunRate),ncol = ncol(restRunRate)),matrix('No Ket',nrow = nrow(nbackRunRate),ncol = ncol(nbackRunRate)))
states <- sapply(1:numClusters, function(K) rep(as.character(K),nrow(restRunRate)))
df.plt <- data.frame(states=as.vector(rbind(states,states)),grps = as.vector(grps),dt=as.vector(rbind(restRunRate,nbackRunRate)))
p <- ggplot(df.plt) + geom_split_violin(aes(x = states, y = dt,fill = grps)) + theme_classic() +
  scale_fill_manual(limits = c('Ket','No Ket'), values = RNcolors) + 
  ylab("Appearance Rate (runs/min)") + xlab("") +theme(text = element_text(size = 8)) +
  theme(legend.title = element_blank()) + 
  scale_x_discrete(limits = 1:numClusters, breaks=1:numClusters, labels = list(clusterNames)) +
  theme(legend.key.size = unit(0.5,'line'))

diffs.full <- lapply(1:numClusters, function(i) t.test(restRunRate[,i],nbackRunRate[,i],paired=TRUE))
print(diffs.full)
diffs <- sapply(diffs.full, function(x) x$p.value)
diffs <- p.adjust(diffs,method = "bonf")
print(diffs)
for(K in 1:numClusters){
  if(diffs[K] < 10^-15){
    p <- p + annotate("text", x = K, y = 1.1*max(rbind(restRunRate,nbackRunRate)),label = "**",color = 'red')
  } else if(diffs[K] < 10^-4){
  	p <- p + annotate("text", x = K, y = 1.1*max(rbind(restRunRate,nbackRunRate)),label = "*",color = 'red')
  }
}

if(numClusters == 5 | numClusters == 6){
  p <- p + theme(axis.text.x = element_text(size=8,colour = clusterColors))
}

save(df.plt,clusterNames,clusterColors,p,diffs.full,file=paste0(masterdir,'analyses/transitionprobabilities/Fig3c__RunRate.RData'))
ggsave(plot = p, filename = paste(masterdir,'analyses/transitionprobabilities/RvNRunRateTime_k',numClusters,'.pdf',sep =""),
  height = 2,width = (numClusters-2) + 0.25, units = "in")