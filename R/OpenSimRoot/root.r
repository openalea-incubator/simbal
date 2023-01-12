library(ggplot2)
library(forcats)
library(ggpubr)


####Graphs####

read.csv("r21.csv",header = TRUE, sep = ";", dec = ",")-> root

root$J <- fct_relevel(root$J, c("J1", "J2", "J3","J4","J5","J6", "J7", "J8", "J9", "J10" ))
levels(root$J)

as.factor(root$Modalities)-> root$Modalities
as.factor(root$Sim)-> root$Sim
as.factor(root$J)-> root$J

a<-ggplot(root, aes(x=J, y=Nuptake ,fill=Sim))+
geom_col(width = 0.7, position='dodge') +
theme_bw()+
theme(axis.line = element_line(color = "black", linetype = "solid"),
axis.text.y=element_text(face="bold", size=14),
axis.text.x= element_text(face = "bold", size = 14),
legend.text = element_text(size=14, face="plain"),
legend.title = element_text(size = 14, face="plain"),
axis.title=element_text(size=16,face="plain"),
legend.position="bottom")+
scale_fill_brewer(palette="PuRd") +
labs(title =" ",
subtitle = " ",
caption = "RStudio v4.2.2",
fill="Maize phenotype",
y="+ N uptake ",
x="Jours")
a


read.csv("r38.csv",header = TRUE, sep = ";", dec = ",")-> root

root$J <- fct_relevel(root$J, c("J1", "J2", "J3","J4","J5","J6", "J7", "J8", "J9", "J10" ))
levels(root$J)

as.factor(root$Modalities)-> root$Modalities
as.factor(root$Sim)-> root$Sim
as.factor(root$J)-> root$J

b<-ggplot(root, aes(x=J, y=Nuptake ,fill=Sim))+
geom_col(width = 0.7, position='dodge') +
theme_bw()+
theme(axis.line = element_line(color = "black", linetype = "solid"),
axis.text.y=element_text(face="bold", size=14),
axis.text.x= element_text(face = "bold", size = 14),
legend.text = element_text(size=14, face="plain"),
legend.title = element_text(size = 14, face="plain"),
axis.title=element_text(size=16,face="plain"),
legend.position="bottom")+
scale_fill_brewer(palette="PuRd") +
labs(title =" ",
subtitle = " ",
caption = "RStudio v4.2.2",
fill="Maize phenotype",
y="+ N uptake ",
x="Jours")
b



read.csv("r65.csv",header = TRUE, sep = ";", dec = ",")-> root

root$J <- fct_relevel(root$J, c("J1", "J2", "J3","J4","J5","J6", "J7", "J8", "J9", "J10" ))
levels(root$J)

as.factor(root$Modalities)-> root$Modalities
as.factor(root$Sim)-> root$Sim
as.factor(root$J)-> root$J

c<-ggplot(root, aes(x=J, y=Nuptake ,fill=Sim))+
geom_col(width = 0.7, position='dodge') +
theme_bw()+
theme(axis.line = element_line(color = "black", linetype = "solid"),
axis.text.y=element_text(face="bold", size=14),
axis.text.x= element_text(face = "bold", size = 14),
legend.text = element_text(size=14, face="plain"),
legend.title = element_text(size = 14, face="plain"),
axis.title=element_text(size=16,face="plain"),
legend.position="bottom")+
scale_fill_brewer(palette="PuRd") +
labs(title =" ",
subtitle = " ",
caption = "RStudio v4.2.2",
fill="Maize phenotype",
y="+ N uptake ",
x="Jours")
c


figure <- ggarrange(a, b, c,
                    labels = c("(A)", "(B)", "(C)"),
                    ncol = 2, nrow = 2)
figure


ggsave("stats.pdf",
  plot = figure,
  width = 11, height = 8)

#boxplot

library(ggplot2)
library(forcats)
library(ggpubr)

read.csv("c.csv",header = TRUE, sep = ";", dec = ",")-> root

as.factor(root$Conductivity)-> root$Conductivity
as.factor(root$Sim)-> root$Sim
as.factor(root$J)-> root$J
as.numeric(root$Nuptake)-> root$Nuptake


p<-ggplot(root, aes(x=Conductivity, y=Nuptake, colour=Sim, fill=Sim))+
        geom_point(position=position_jitterdodge(dodge.width=0.7), size=1.2) + 
        geom_boxplot(alpha=0.5, position = position_dodge(width=0.8), fatten=NULL)+ 
  scale_fill_viridis_d() +
  scale_color_viridis_d()+ 
        stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), width=0.65,size = 1.5, linetype = "solid",position = position_dodge(width=0.8))+ 
        ylab("Nuptake")+ 
        labs( x = "Conductivity", y = "Nuptake",
              title ="",
              subtitle = "",
              caption = "",
              color=" ",
              fill=" ") +
  theme_bw()+
        theme(axis.text.x= element_text(face = "plain", size = 14),
              axis.text.y= element_text(face = "plain", size = 14),
              legend.position = "right",
              legend.text = element_text(size=14, face="plain"),
              legend.title = element_text(size = 14, face="plain"),
              axis.title=element_text(size=16,face="plain")) 
p  

ggsave("ppp.pdf",
  plot = p,
  width = 11, height = 8)


####Stats####
#Anova 2 facteurs
library(rstatix)
library(Rmisc)
library(car)

read.csv("r.csv",header = TRUE, sep = ";", dec = ",")-> r

as.factor(r$Modalities)-> r$Modalities
as.factor(r$Sim)-> r$Sim
as.factor(r$J)-> r$J
as.numeric(r$Nuptake)-> r$Nuptake

mydata<-aggregate(r$Nuptake,
                  by=list(J=r$J,
                          Sim=r$Sim),
                  FUN=mean)
d<-as.matrix((mydata))


f<-friedman_test(mydata, x ~ Sim | J)
f#test de Friedman p-value=206.10^-11, significatif
#normalité


pwc<- r %>%
  wilcox_test(Nuptake~Sim)
pwc

write.csv(pwc,"pwc.csv",row.names = TRUE )


#tous est non significatif, pas de diff�rence entre les Sim
#Le test de Wilcoxon entre les traitements ne r�v�le aucunes diff�rences statistiques.

#graph
read.csv("oc357.csv",header = TRUE, sep = ";", dec = ",")-> root

root$J <- fct_relevel(root$J, c("J1", "J2", "J3","J4","J5","J6", "J7", "J8", "J9", "J10","J11","J12","J13","J14","J15" ))
levels(root$J)

as.factor(root$Conductivity)-> root$Conductivity
as.factor(root$Sim)-> root$Sim
as.factor(root$J)-> root$J

a<-ggplot(root, aes(x=J, y=Nuptake ,fill=Sim))+
geom_col(width = 0.7, position='dodge') +
theme_bw()+
theme(axis.line = element_line(color = "black", linetype = "solid"),
axis.text.y=element_text(face="bold", size=14),
axis.text.x= element_text(face = "bold", size = 14),
legend.text = element_text(size=14, face="plain"),
legend.title = element_text(size = 14, face="plain"),
axis.title=element_text(size=16,face="plain"),
legend.position="bottom")+
scale_fill_brewer(palette="PuRd") +
labs(title =" ",
subtitle = " ",
caption = "RStudio v4.2.2",
fill="Maize phenotype",
y="+ N uptake ",
x="Jours")
a

ggsave("a.pdf",
  plot = a,
  width = 11, height = 8)



####Stats2####
#Anova 2 facteurs
library(rstatix)
library(Rmisc)
library(car)

read.csv("octave.csv",header = TRUE, sep = ";", dec = ",")-> r

as.factor(r$Conductivity)-> r$Conductivity
as.factor(r$Sim)-> r$Sim
as.factor(r$J)-> r$J
as.numeric(r$Nuptake)-> r$Nuptake

mydata<-aggregate(r$Nuptake,
                  by=list(J=r$J,
                          Sim=r$Sim),
                  FUN=mean)
d<-as.matrix((mydata))


f<-friedman_test(mydata, x ~ Sim | J)
f#test de Friedman p-value=206.10^-11, significatif
#normalité


pwc<- r %>%
  wilcox_test(Nuptake~Sim)
pwc

write.csv(pwc,"pwc.csv",row.names = TRUE )


#tous est non significatif, pas de diff�rence entre les Sim
#Le test de Wilcoxon entre les traitements ne r�v�le aucunes diff�rences statistiques.



#boxplot

library(ggplot2)
library(forcats)
library(ggpubr)

read.csv("j15.csv",header = TRUE, sep = ";", dec = ",")-> root

as.factor(root$Conductivity)-> root$Conductivity
as.factor(root$Sim)-> root$Sim
as.factor(root$J)-> root$J
as.numeric(root$Nuptake)-> root$Nuptake


p<-ggplot(root, aes(x=Conductivity, y=Nuptake, colour=Sim, fill=Sim))+
        geom_point(position=position_jitterdodge(dodge.width=0.7), size=1.2) + 
        geom_boxplot(alpha=0.5, position = position_dodge(width=0.8), fatten=NULL)+ 
  scale_fill_viridis_d() +
  scale_color_viridis_d()+ 
        stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), width=0.65,size = 1.5, linetype = "solid",position = position_dodge(width=0.8))+ 
        ylab("Nuptake")+ 
        labs( x = "Conductivity", y = "Nuptake",
              title ="",
              subtitle = "",
              caption = "",
              color=" ",
              fill=" ") +
  theme_bw()+
        theme(axis.text.x= element_text(face = "plain", size = 14),
              axis.text.y= element_text(face = "plain", size = 14),
              legend.position = "right",
              legend.text = element_text(size=14, face="plain"),
              legend.title = element_text(size = 14, face="plain"),
              axis.title=element_text(size=16,face="plain")) 
p  

ggsave("pp.pdf",
  plot = p,
  width = 11, height = 8)
