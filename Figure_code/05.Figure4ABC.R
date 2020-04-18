library(psych)
library(ggrepel)
library(RColorBrewer)
library(viridis)
library(ggplot2)

###calculate cor
nor_t=read.table("rlogTransformation.txt",header = T,sep = "\t")
nor_ACE2=nor_t[14039]
cor_fdr=corr.test(nor_ACE2,nor_t,use = "complete", method = "pearson", adjust = "BH")
r=t(cor_fdr$r)
p=t(cor_fdr$p)
res=cbind(r,p)
colnames(res)=c("cor","padj")
res=as.data.frame(res)
res$gene=rownames(res)
########################################################## 1. Figure 4A #####################################################################
res_Sig_cor_positive=subset(res,res$padj<0.001 & res$cor>0.6)
res_Sig_cor_negative=subset(res,res$padj<0.001 & res$cor< -0.6)
res_other=subset(res,res$padj>0.001 | abs(res$cor) <0.6)
res_Sig_cor_positive$Tag="Sig_cor_positive"
res_Sig_cor_negative$Tag="Sig_cor_negative"
res_other$Tag="Non-Sig"
res_plot=rbind(res_Sig_cor_positive,res_Sig_cor_negative,res_other)
res_plot1=res_plot[-match("ACE2",res_plot$gene),]
la=subset(res_plot1, res_plot$cor > 0.73 & res_plot$padj<0.01)
gene_ANPEP=res_plot1[match("ANPEP",res_plot1$gene),]
la1=rbind(la,gene_ANPEP)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
c=getPalette(9)
pdf("Figure 4A.pdf")
ggplot(res_plot1,aes(cor,-log10(padj)))+ 
  geom_point(aes(color = Tag),size = 1)+
  geom_text_repel(data=la1,aes(x=cor,y=-log10(padj),label=gene),size=3)+
  scale_color_manual(values = c("lightgrey", c[2], c[1])) +
  geom_hline(yintercept = 3,linetype = 2) +
  geom_vline(xintercept = c(-0.6, 0.6), linetype = 2)+
  theme(legend.title=element_blank())+theme(legend.background = element_blank())+theme(legend.key = element_blank())+
  theme(legend.text = element_text(size=20))+theme_bw()+theme(axis.text.x = element_text(size = 20))+theme(axis.text.y = element_text(size = 20))
dev.off()
########################################################## 2. Figure 4B #####################################################################
go=read.table("GO input.txt",sep = "\t")
colnames(go)=c("ID","Name","log10P")
go$log10P=-go$log10P
go_order=go[order(go$log10P,decreasing = T),]
go_order$Name=factor(go_order$Name,levels = rev(go_order$Name))
pdf("Figure 4B.pdf")
ggplot(go_order,aes(Name,log10P,fill=log10P))+
  geom_bar(stat="identity")+coord_flip()+
  theme_bw()+
  scale_fill_viridis(option = "B",direction = -1)+
  theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=0.5,colour="black"))
dev.off()
########################################################## 3. Figure 4C #####################################################################
###TMEM27
p_TMEM27=nor_t[,match(c("ACE2","TMEM27"),colnames(nor_t))]
p_TMEM27$celltype=factor(rownames(p_TMEM27))
pdf("Figure 4C-TMEM27.pdf")
ggplot(p_TMEM27,mapping=aes(x=ACE2,y=TMEM27)) + 
  geom_point() +
  geom_smooth(method = lm,colour=c[2],fill=c[2])+
  theme(legend.background = element_blank())+
  theme(legend.key = element_blank())+
  theme(legend.text = element_text(size=10))+
  theme_bw()+
  ggtitle(paste("cor",as.character(round(res_plot1[match("TMEM27",rownames(res_plot1)),][1],5)),
                "padj",as.character(round(res_plot1[match("TMEM27",rownames(res_plot1)),][2],5)),
                sep = " "))
dev.off()
###IDO2
p_IDO2=nor_t[,match(c("ACE2","IDO2"),colnames(nor_t))]
p_IDO2$celltype=factor(rownames(p_IDO2))
pdf("Figure 4C-IDO2.pdf")
ggplot(p_IDO2,mapping=aes(x=ACE2,y=IDO2)) + 
  geom_point()+
  geom_smooth(method = lm,colour=c[2],fill=c[2])+
  theme(legend.background = element_blank())+
  theme(legend.key = element_blank())+
  theme(legend.text = element_text(size=10))+
  theme_bw()+
  ggtitle(paste("cor",as.character(round(res_plot1[match("IDO2",rownames(res_plot1)),][1],5)),
                "padj",as.character(round(res_plot1[match("IDO2",rownames(res_plot1)),][2],5)),
                sep = " "))
dev.off()
###DNAJC12
p_DNAJC12=nor_t[,match(c("ACE2","DNAJC12"),colnames(nor_t))]
p_DNAJC12$celltype=factor(rownames(p_DNAJC12))
pdf("Figure 4C-DNAJC12.pdf")
ggplot(p_DNAJC12,mapping=aes(x=ACE2,y=DNAJC12)) + 
  geom_point() +
  geom_smooth(method = lm,colour=c[2],fill=c[2])+
  theme(legend.background = element_blank())+
  theme(legend.key = element_blank())+
  theme(legend.text = element_text(size=10))+
  theme_bw()+
  ggtitle(paste("cor",as.character(round(res_plot1[match("DNAJC12",rownames(res_plot1)),][1],5)),
                "padj",as.character(round(res_plot1[match("DNAJC12",rownames(res_plot1)),][2],5)),
                sep = " "))
dev.off()
###ANPEP
p_ANPEP=nor_t[,match(c("ACE2","ANPEP"),colnames(nor_t))]
p_ANPEP$celltype=factor(rownames(p_ANPEP))
pdf("Figure 4C-ANPEP.pdf")
ggplot(p_ANPEP,mapping=aes(x=ACE2,y=ANPEP)) + 
  geom_point() +
  geom_smooth(method = lm,colour=c[2],fill=c[2])+
  theme(legend.background = element_blank())+
  theme(legend.key = element_blank())+
  theme(legend.text = element_text(size=10))+
  theme_bw()+
  ggtitle(paste("cor",as.character(round(res_plot1[match("ANPEP",rownames(res_plot1)),][1],5)),
                "padj",as.character(round(res_plot1[match("ANPEP",rownames(res_plot1)),][2],5)),
                sep = " "))
dev.off()
