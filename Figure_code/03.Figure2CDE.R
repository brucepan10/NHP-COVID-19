library(data.table)
library(ggplot2)
library(viridis)
########################################################## 1. Figure 2C #####################################################################
exp=as.data.frame(fread("ACE2_TMPRSS2_normalize_log_mat.csv"))
exp=exp[-1,]
rownames(exp)=exp[,1]
exp=exp[,-1]
exp_t=as.data.frame(t(exp))
info=read.table("organ9_cellinfo.txt",header = T,sep = "\t")
sum=cbind(exp_t,info)
sum$Both=sum$ACE2>0 & sum$TMPRSS2>0
sum_both=subset(sum,sum$Both=="TRUE")
pdf("Figure2C")
ggplot()+
  geom_point(data=sum,aes(x=X_umap1, y =X_umap2),alpha=1, size=0.0001,colour="#0D0887FF")+
  theme_gray() +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size=25),panel.background = element_blank(),axis.line = element_line(colour = "black"))+
  geom_point(data=sum_both,aes(x=X_umap1, y =X_umap2),alpha=1, size=1,colour="#F0F921FF")+
  ggtitle("ACE2/TMPRSS2+")+
  theme_bw() +
  theme(panel.grid =element_blank()) +
  theme(axis.text = element_blank()) + 
  theme(axis.ticks = element_blank())
dev.off()
########################################################## 2. Figure 2D #####################################################################
sub_both=subset(sum,sum$ACE2>0 & sum$TMPRSS2>0)
table_both=as.data.frame(table(sub_both$celltype))
colnames(table_both)=c("celltype","num_both")
table_all_both=merge(table_all,table_both,by="celltype")
table_all_both$Ratio_both=table_all_both$num_both/table_all_both$Total
table_all_both_ss=subset(table_all_both,table_all_both$Ratio_both>0)
table_all_both_ss_order=table_all_both_ss[order(table_all_both_ss$Ratio_both),]
table_all_both_ss_order$celltype=factor(table_all_both_ss_order$celltype,levels = table_all_both_ss_order$celltype)
table_all_both_ss_order_plot=subset(table_all_both_ss_order,table_all_both_ss_order$Ratio_both>0)
pdf("Figure2D")
ggplot(table_all_both_ss_order_plot,aes(celltype,Ratio_both,fill=celltype))+
  geom_bar(stat="identity")+coord_flip()+
  scale_fill_manual(values = viridis(23,direction = -1,option = "C"))+
  theme_bw()+
  theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=0.5,colour="black"))
dev.off()
########################################################## 3. Figure 2E #####################################################################
sum$celltype=factor(sum$celltype,levels = rev(as.character(table_all$celltype)))
colnames(table_all)=c("celltype","Total")
sum$celltype=factor(sum$celltype,levels = as.character(table_all$celltype))
table_all$mean_ACE2=tapply(sum$ACE2,sum$celltype,mean)
table_all$mean_TMPRSS2=tapply(sum$TMPRSS2,sum$celltype,mean)
sub_ACE2=subset(sum,sum$ACE2>0)
sub_TMPRSS2=subset(sum,sum$TMPRSS2>0)
table_ACE2=as.data.frame(table(sub_ACE2$celltype))
table_TMPRSS2=as.data.frame(table(sub_TMPRSS2$celltype))
table_all$num_ACE2=table_ACE2$Freq
table_all$num_TMPRSS2=table_TMPRSS2$Freq
table_all$Ratio_ACE2=table_all$num_ACE2/table_all$Total
table_all$Ratio_TMPRSS2=table_all$num_TMPRSS2/table_all$Total
table_all_sub=table_all
table_all=table_all_sub
table_all_1=table_all[,c(1,3,7)]
table_all_2=table_all[,c(1,4,8)]
table_all_1$Gene="ACE2"
table_all_2$Gene="TMPRSS2"
colnames(table_all_1)=c("celltype","Mean","Ratio","Gene")
colnames(table_all_2)=c("celltype","Mean","Ratio","Gene")
table_plot=rbind(table_all_1,table_all_2)
table_plot$Gene=factor(table_plot$Gene,levels = c("TMPRSS2","ACE2"))
###ACE2 bubble plot 
table_all_1_order=table_all_1[order(table_all_1$Ratio,decreasing = T),]
table_plot$celltype=factor(table_all_1$celltype,levels = table_all_1_order$celltype)
table_all_1_order$celltype=factor(table_all_1_order$celltype,levels = table_all_1_order$celltype)
pdf("Figure2E-ACE2.pdf")
ggplot(table_all_1_order,aes_(x = ~celltype, y = ~Gene, size = ~Ratio))+
  geom_point() +
  aes_string(color="Mean") +
  theme_bw()+
  theme(axis.text.x = element_text(size=8,angle=90, hjust=1, vjust=1))+
  scale_color_gradientn(colours = c("grey",viridis(50,direction = -1,option = "C")))
dev.off()
###TMPRSS2 bubble plot 
table_all_2$celltype=factor(table_all_2$celltype,levels = table_all_1_order$celltype)
pdf("Figure2E-TMPRSS2.pdf")
ggplot(table_all_2,aes_(x = ~celltype, y = ~Gene, size = ~Ratio))+
  geom_point() +
  aes_string(color="Mean") +
  theme_bw()+
  theme(axis.text.x = element_text(size=8,angle=90, hjust=1, vjust=1))+
  scale_color_gradientn(colours = c("grey",viridis(50,direction = -1,option = "C")))
dev.off()
