########################################################## Figure 5E #####################################################################
library(Seurat)
library("ggplot2")
load("kidney.integrated.Find.RData") # input kidney Rdata
head(kidney.integrated.Find@meta.data)

ACE2_plot<-data.frame(ACE2=as.vector(kidney.integrated.Find@assays$RNA@data["ACE2",]),cluster=as.vector(kidney.integrated.Find@meta.data$cluster))
back<-ACE2_plot

pos<-ACE2_plot$cluster%in%c(0:1,3:20)
ACE2_plot<-ACE2_plot[pos,]
ACE2_plot$cluster<-factor(ACE2_plot$cluster)

### annotation celltype
ACE2_plot$celltype="NNN"
ACE2_plot$celltype[ACE2_plot$cluster==0]="Proximal tubule S1"
ACE2_plot$celltype[ACE2_plot$cluster==4]="Proximal tubule S3"
ACE2_plot$celltype[ACE2_plot$cluster==1]="Thick ascending limb cell"
ACE2_plot$celltype[ACE2_plot$cluster==3]="Thick ascending limb cell"
ACE2_plot$celltype[ACE2_plot$cluster==7]="Thick ascending limb cell"
ACE2_plot$celltype[ACE2_plot$cluster==18]="Thick ascending limb cell"
ACE2_plot$celltype[ACE2_plot$cluster==12]="Thin limb cell"
ACE2_plot$celltype[ACE2_plot$cluster==6]="Principal cell"
ACE2_plot$celltype[ACE2_plot$cluster==5]="Principal cell"
ACE2_plot$celltype[ACE2_plot$cluster==8]="Principal cell"
ACE2_plot$celltype[ACE2_plot$cluster==14]="Principal cell"
ACE2_plot$celltype[ACE2_plot$cluster==9]="Distal convoluted tubule cell"
ACE2_plot$celltype[ACE2_plot$cluster==19]="Distal convoluted tubule cell"
ACE2_plot$celltype[ACE2_plot$cluster==10]="Connecting tubule cell"
ACE2_plot$celltype[ACE2_plot$cluster==11]="Endothelial cell"
ACE2_plot$celltype[ACE2_plot$cluster==13]="Intercalated cell"
ACE2_plot$celltype[ACE2_plot$cluster==15]="Intercalated cell"
ACE2_plot$celltype[ACE2_plot$cluster==17]="Intercalated cell"
ACE2_plot$celltype[ACE2_plot$cluster==16]="Podocyte"
ACE2_plot$celltype[ACE2_plot$cluster==20]="Stromal cell"

## plot proporation bar graph
table_all=as.data.frame(table(ACE2_plot$celltype))
sub_ACE2=subset(ACE2_plot,ACE2_plot$ACE2>0)
table_ACE2=as.data.frame(table(sub_ACE2$celltype))
table_all$num_ACE2=table_ACE2$Freq
table_all$Ratio_ACE2=table_all$num_ACE2/table_all$Freq
pdf("figure5E_monkey_kidney_ACE2_ratio.pdf",11,7)
ggplot(table_all,aes(Var1,Ratio_ACE2,fill=Var1))+	
	geom_bar(stat="identity")+
	theme_classic()+
	theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.y=element_text(size=15))+
	theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "inches"))+
	guides(fill=FALSE)
dev.off()




