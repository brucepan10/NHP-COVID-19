##########################################################  Figure 5H right ##################################################################### 
library(ggplot2)
library(viridis)
library(data.table)

# input Human kidney matrix
exp=as.data.frame(fread("HumanKidney_matrix.csv",sep=","))
exp_t=as.data.frame(t(exp))
exp_t[1:3,1:3]

# input cluster, umap coordinates information of Human kidney
scanpy <- read.csv(file="Scanpy_cell_info.csv", header=TRUE, sep=",")
scanpy_umap=read.table("Scanpy_X_umap.csv",header = T,sep = ",")
scanpy$X_umap1 <- scanpy_umap$X_umap1
scanpy$X_umap2 <- scanpy_umap$X_umap2
colnames(scanpy)[1]="ID"

all(scanpy$ID==rownames(exp_t))
sub_exp_t<-exp_t[,c("ACE2","IL6R")]
sub_exp_t$X_umap1<-scanpy$X_umap1
sub_exp_t$X_umap2<-scanpy$X_umap2

sum<-sub_exp_t
sum$exp_ACE2=sum$ACE2>0
sum$exp_IL6R=sum$IL6R>0
sum$Both=sum$ACE2>0 & sum$IL6R>0

sum_both=subset(sum,sum$Both=="TRUE")

pdf("Human_kidney_ACE2_IL6R.pdf",5.2,5)
ggplot()+
  geom_point(data=sum,aes(x=X_umap1, y =X_umap2),alpha=1, size=0.0001,colour="#0D0887FF")+
  theme_gray() +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size=25),panel.background = element_blank(),axis.line = element_line(colour = "black"))+
  geom_point(data=sum_both,aes(x=X_umap1, y =X_umap2),alpha=1, size=1,colour="#F0F921FF")+
  ggtitle("HumanKidney:ACE2/IL6R+")+
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.text = element_blank()) + 
  theme(axis.ticks = element_blank()) 
dev.off()
