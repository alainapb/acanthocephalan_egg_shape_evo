library(ape)
library(plyr)
library(tidyverse)
library(ouch)
library(dplyr)
library(phytools)
library(geiger)
library(caper)
library(picante)
library(gridExtra)
library(OUwie)
library(patchwork)
library(pmc)
library(car)
library(corHMM)


####ALL the species ####
all<-read.csv("~/Documents/UNL/6Spring_2020/acanth_eggs_MS/ALL_egg_data.csv")
all$Binomial_name<-gsub(" ", "_", all$Binomial_name)
#subset(all, all$Egg_length_avg_um != " " )->all
names(all)
subset(all, all$Species != "")->all

ER_class<-lm(Elongation_ratio~ Class - 1, data=all)
summary(ER_class)

anova(ER_all<-lm(Elongation_ratio ~ Class * Order * Family -1, data = all))
summary(ER_all)


Arch<-subset(all, all$Class == "Archiacanthocephala")
Eo<-subset(all, all$Class == "Eoacanthocephala")
Palae<-subset(all, all$Class == "Palaeacanthocephala")

#Archiacanthocephala
anova(ER_Arch_nest<-lm(Elongation_ratio ~ Order * Family-1, data =Arch)) #only order is significant
anova(ER_Arch_order<-lm(Elongation_ratio ~ Order - 1, data = Arch)) #all significant
anova(ER_Arch_family<-lm(Elongation_ratio ~ Family - 1, data = Arch)) # all significant

#Eoacanthocephala
anova(ER_Eo_nest<-lm(Elongation_ratio ~ Order * Family-1, data = Eo)) #only order is significant
anova(ER_Eo_order<-lm(Elongation_ratio ~ Order - 1, data = Eo)) #all significant
anova(ER_Eo_family<-lm(Elongation_ratio ~ Family - 1, data = Eo)) # all significant

#Palaeacanthocephala
anova(ER_Palae_nest<-lm(Elongation_ratio ~ Order * Family-1, data = Palae)) 
anova(ER_Palae_order<-lm(Elongation_ratio ~ Order - 1, data = Palae)) #all significant
anova(ER_Palae_family<-lm(Elongation_ratio ~ Family - 1, data = Palae)) # all significant



####Figures for manuscript####

#Figure2
b1<-ggplot(filter(all, Class != "Polyacanthocephala"), aes(x=Elongation_ratio))+
  geom_histogram(binwidth = .50, color="black", fill="grey")+
  labs(x="Egg Shape (ER)", y="Frequency")+
  xlim(1, 8)+
  ylim(0,60)+
  theme_classic()+
  facet_wrap(~Class, 3,1)

b2<-ggplot(filter(all, Class != "Polyacanthocephala"), aes(x=Egg_length_avg_um))+
  geom_histogram(binwidth = 20, color="black", fill="grey")+
  xlim(10, 210)+
  labs(x="Egg Length (um)", y="")+
  theme_classic()+
  theme(axis.text.y = element_blank())+
  facet_wrap(~Class, 3,1)

#grid.arrange(b1, b2, ncol=2)
#using patchwork
b1 + b2

##cv=stdev/mean

filter(all, !is.na(Elongation_ratio))%>%
  group_by(Class) %>%
  transmute(shape_mean = mean(Elongation_ratio), 
            shape_sd = sd(Elongation_ratio),
            shape_cv = shape_sd/shape_mean) %>%
  summarise(shape_cv = unique(shape_cv))


filter(all, !is.na(Egg_length_avg_um))%>%
  group_by(Class) %>%
  transmute(length_mean = mean(Egg_length_avg_um), 
            length_sd = sd(Egg_length_avg_um),
            length_cv = length_sd/length_mean) %>%
  summarise(length_cv = unique(length_cv))

#Figure 3
ggplot(filter(all, Class != "Polyacanthocephala"), aes(x=Elongation_ratio, y=Egg_length_avg_um, pch = Class)) +
  geom_point(size = 3, show.legend = F)+
  theme_classic()+
  labs(x="Egg Shape (ER)", y="Egg Size (length, um)")+
  #geom_smooth(method="lm", alpha=0, color="grey")+
  facet_wrap(~Class, 3,1)+
  theme(panel.border = element_rect(fill=NA, color="black", size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size=10),
        text = element_text(size=10))

cor.test(all$Elongation_ratio, all$Egg_length_avg_um)

Arch<-subset(all, all$Class == "Archiacanthocephala")
cor.test(Arch$Elongation_ratio, Arch$Egg_length_avg_um)

Eo<-subset(all, all$Class == "Eoacanthocephala")
cor.test(Eo$Elongation_ratio, Eo$Egg_length_avg_um) 

Palae<-subset(all, all$Class == "Palaeacanthocephala")
cor.test(Palae$Elongation_ratio, Palae$Egg_length_avg_um) 


####Adding in the tree####
tree<-read.nexus("~//Documents/UNL/6Spring_2020/acanth_eggs_MS/acanth_tree/alignments_8_4_20/acanth.trees")
plot(tree, cex=0.3)

tip<-c("Polymorphus_diploinflatus", "Polymorphus_contortus")
plot(drop.tip(tree, tip))

tree1<-drop.tip(tree, tip)

#plot(tree1, cex=0.6)
#nodelabels()

##what species are in the data but not in the tree
setdiff(all$Binomial_name, tree1$tip.label)

##fix mispelled names
all$Binomial_name<-gsub("Andracantha_phalacrococoracis" ,"Andracantha_phalacrocoracis", all$Binomial_name)
all$Binomial_name<-gsub("Echinorhynchus_salmonis_","Echinorhynchus_salmonis", all$Binomial_name)

##Fix mispellings in the tree
tree1$tip.label<-gsub("Pomphorhynchus_bulbocoli", "Pomphorhynchus_bulbocolli", tree1$tip.label)
tree1$tip.label<-gsub("Pseudocorynosoma_anatarium","Pseudocorynosoma_anatrium", tree1$tip.label)
tree1$tip.label<-gsub("Neoechinorhynchus_saginata", "Neoechinorhynchus_saginatus", tree1$tip.label)
tree1$tip.label<-gsub("Pallisentis_naqpurensis", "Pallisentis_nagpurensis", tree1$tip.label)
tree1$tip.label<-gsub("Koronacantha_pectinara", "Koronacantha_pectinaria", tree1$tip.label)

##Now remove species that are in the data but not in the tree
setdiff(all$Binomial_name, tree1$tip.label)->bye
all[-which(all$Binomial_name %in% bye),]->tree_data

##Now make the tree_data in the same order as the tree1 tip labels
tree_data <- tree_data[match(tree1$tip.label, tree_data$Binomial_name),]

tree_data$Order<-as.character(tree_data$Order)


####Check the tree at the Order and Family level####
tree_data$Family<-gsub("Neoechonorhynchidae","Neoechinorhynchidae", tree_data$Family)
family<-tree_data$Family
names(family)<-tree_data$Binomial_name

plot(tree1, cex=0.3)
cols<-setNames(palette()[1:16], sort(unique(family)))
tiplabels(pie=to.matrix(family, sort(unique(family))), piecol=cols, cex=0.2)
add.simmap.legend(colors=cols, prompt=FALSE, x=0.9*par()$usr[1],
                  y=20+max(nodeHeights(tree1)), fsize=0.3)



order<-tree_data$Order
names(order)<-tree_data$Binomial_name

plot(tree1, cex=0.3)
cols<-setNames(palette()[1:8], sort(unique(order)))
tiplabels(pie=to.matrix(order, sort(unique(order))), piecol=cols, cex=0.2)
add.simmap.legend(colors=cols, prompt=FALSE, x=.9*par()$usr[1],
                  y=9+max(nodeHeights(tree1)), fsize=0.3)




#### Prep data for test for phylogenetic signal ####


##subset tree to only include species that we have egg data for :(
phy_data<-subset(tree_data, tree_data$Egg_length_avg_um != " ")
setdiff(tree1$tip.label, phy_data$Binomial_name)->nodata

plot(drop.tip(tree1, nodata))

tree2<-drop.tip(tree1, nodata)

##Make the phy_data in the same order as the tree2 tip labels
phy_data <- phy_data[match(tree2$tip.label, phy_data$Binomial_name),]


#named vector of average ER
shape<-phy_data$Elongation_ratio
names(shape)<-phy_data$Binomial_name

contMap(tree2, shape, outline= TRUE, fsize = .8)

### Test of phylogenetic signal ###

#egg shape
fitContinuous(tree2, shape, model="lambda") ## lambda= 0.840597 ; aicc= 178.6164 ##Best Model
fitContinuous(tree2, shape, model="white") #aicc= 220.7986

#make a table
phy_signal_shape<-c(178.616,220.7986)
names(phy_signal_shape)<-c("lambda", "star")

aic_phy_s <- tibble(aicc=phy_signal_shape, model=names(phy_signal_shape))

shape_signal<-aic_phy_s %>% 
  mutate(deltaaicc=aicc-min(aicc), 
         pow=exp(-.5*deltaaicc),
         sumpow=sum(pow),
         waicc=pow/sumpow) %>%
  dplyr::select(c(2,1,3,6)) %>%
  arrange(desc(waicc))

###figure illustrating phy signal ###

#shape_w_l <- pmc(tree2, shape, "white", "lambda", nboot = 2000, mc.cores = 2)
#saveRDS(shape_w_l, file="~/Documents/UNL/6Spring_2020/acanth_eggs_MS/acanth_eggs_2020/shape_w_l")

shape_w_l<-readRDS("~/Documents/UNL/6Spring_2020/acanth_eggs_MS/acanth_eggs_2020/shape_w_l")


#plot results

dists3 <- data.frame(whitenoise = shape_w_l$null, lambda = shape_w_l$test)
dists3 %>% 
  gather(dist, value) %>%
  ggplot(aes(value, fill = dist)) + 
  geom_density(alpha = 0.5) + 
  geom_vline(xintercept = shape_w_l$lr)+
  scale_x_continuous(limits=c(0,80))+
  scale_y_continuous(limits=c(0,1))+
  labs(fill = " ", x= " ", y= "Density")+
  theme_bw()

### Evolutionary Models ###

fitContinuous(tree2, shape, model="BM") #aicc= 199.00467
fitContinuous(tree2, shape, model="OU") #alpha = 2.7182 ; aicc= 197.21635
fitContinuous(tree2, shape, model= "EB") # alpha = 0; aicc = 201.227

shape_evo<-c(199.00467,197.21635, 201.227)
names(shape_evo)<-c("BM", "OU", "EB")

aic_evo_s <- tibble(aicc=shape_evo, model=names(shape_evo))

shape_evo<-aic_evo_s %>% 
  mutate(deltaaicc=aicc-min(aicc), 
         pow=exp(-.5*deltaaicc),
         sumpow=sum(pow),
         waicc=pow/sumpow) %>%
  dplyr::select(c(2,1,3,6)) %>%
  arrange(desc(waicc))




####Power Analyses for Evolution Models ####

##Egg shape

#shape_BM_OU <- pmc(tree2, shape, "BM", "OU", nboot = 2000, mc.cores = 2)
#saveRDS(shape_BM_OU, file="~/Documents/UNL/6Spring_2020/acanth_eggs_MS/acanth_eggs_2020/shape_BM_OU")

shape_BM_OU<-readRDS("~/Documents/UNL/6Spring_2020/acanth_eggs_MS/acanth_eggs_2020/shape_BM_OU")

#plot results

dists <- data.frame(BM = shape_BM_OU$null, OU = shape_BM_OU$test)
dists %>% 
  gather(dist, value) %>%
  ggplot(aes(value, fill = dist)) + 
  geom_density(alpha = 0.4) + 
  geom_vline(xintercept = shape_BM_OU$lr)+
  labs(fill = " ", x = " ", y = "Density")+
  theme_bw()




#### RUNNING PGLS ####
## Combine all data (egg length, egg shape, membrane #, polar prolongation, int habitat, def habitat) into a single data.frame
#write.csv(dat, file="~/Documents/UNL/6Spring_2020/acanth_eggs_MS/acanth_eggs_2020/combined.csv")
combined<-read.csv("~/Documents/UNL/6Spring_2020/acanth_eggs_MS/acanth_eggs_2020/combined.csv")

dat<-data.frame(taxa=tree2$tip.label,
                shape=phy_data$Elongation_ratio,
                int_hab=as.factor(phy_data$Int_habitat_terr.),
                def_hab=as.factor(phy_data$Def_habitat_terr.),
                comb=combined$combined)

## constructs the data.frame of the taxa names and data in the order of the tree tips
cdat<-comparative.data(data=dat, phy=tree2, names.col = "taxa")
print(cdat)


## model testing of egg shape (0 = aquatic , 1 = terrestrial)

#intermediate host habitat only!

#lambda from fitContinous
anova(pgls(shape~ int_hab -1, cdat, lambda = .84)) # aquatic has more elongated eggs
summary(pgls(shape~ int_hab -1, cdat, lambda = .84)) ##used for S4 table


#definitive host habitat only (0 = aquatic , 1 = semiaquatic , 2 = terrestial)

#lambda from fitContinous
anova(pgls(shape~ def_hab -1, cdat, lambda = .84)) # aquatic has more elongated eggs
summary(pgls(shape~ def_hab -1, cdat, lambda = .84)) # used for S4 TABLE


##Both intermediate and definitve host habitats 

lam_shape<-pgls(shape~ int_hab + def_hab -1, cdat, lambda = .84) #intermediate host habitat p<0.05
anova(lam_shape) # how much of the variance is explained by each predictor 
summary(lam_shape) # whether the coefficents are significantly different from zero


com_shape<-pgls(shape ~ comb -1, cdat, lambda = .84)
anova(com_shape)
summary(com_shape)


#### OUWIE ####

#make Class dataframe for OUwie (names, regimes, trait)
class_data<-data.frame(tree2$tip.label)
class_data$regime<-phy_data$Class
class_data$trait<-shape

# estimate the regimes for internal nodes uses "corHMM",
# uses a single rate model of evolution and marginal reconstruction of ancestral states
p<-rayDISC(tree2, class_data[,c(1,2)], model="ER", node.states="marginal")

# test regime models with Class regimes (1=Arch, 2=Eo, 3=Palae)
OUwie(p$phy, class_data, model="BM1") #random   == aicc = 204.48, optima = 3.61
OUwie(p$phy, class_data, model="OU1") #1 regime == aicc = 190.77, optima = 3.64
OUwie(p$phy, class_data, model="OUM") #3 regime == aicc = 194.83, optima1 = 1.58, optima2 = 2.13, optima3 = 4.25
OUwie(p$phy, class_data, model="OUMV")#3 reg +  == aicc = 187.36, optima1 = 1.75, optima2 = 2.22, optima3 = 4.21


###Def host table ###
class_evo<-c(204.48,190.77, 194.83, 187.36)
names(class_evo)<-c("BM", "OU", "OUM", "OUMV")

aic_class <- tibble(aicc=class_evo, model=names(class_evo))

class_evo_table<-aic_class %>% 
  mutate(deltaaicc=aicc-min(aicc), 
         pow=exp(-.5*deltaaicc),
         sumpow=sum(pow),
         waicc=pow/sumpow) %>%
  dplyr::select(c(2,1,3,6)) %>%
  arrange(desc(waicc))


#make definitive host dataframe for OUwie (names, regimes, trait)
def_data<-data.frame(tree2$tip.label)
def_data$regime<-dat$def_hab
def_data$trait<-shape

# estimate the regimes for internal nodes uses "corHMM",
# uses a single rate model of evolution and marginal reconstruction of ancestral states
p1<-rayDISC(tree2, def_data[,c(1,2)], model="ER", node.states="marginal")

# test regime models with definitive host habitat regimes (0 = aquatic , 1 = semiaquatic , 2 = terrestial)
OUwie(p1$phy, def_data, model="BM1") #random   == aicc = 204.49, optima = 3.61
OUwie(p1$phy, def_data, model="OU1") #1 regime == aicc = 190.77, optima = 3.65
OUwie(p1$phy, def_data, model="OUM") #3 regime == aicc = 190.28, optima0 = 4.16, optima1 = 2.77, optima2 = 1.22
OUwie(p1$phy, def_data, model="OUMV")#3 reg +  == aicc = 180.61, optima0 = 4.20, optima1 = 3.26, optima2 = 2.02


###Def host table ###
def_evo<-c(204.49,190.77, 190.28, 180.61)
names(def_evo)<-c("BM", "OU1", "OU3", "OU3+V")

aic_def <- tibble(aicc=def_evo, model=names(def_evo))

def_evo_table<-aic_def %>% 
  mutate(deltaaicc=aicc-min(aicc), 
         pow=exp(-.5*deltaaicc),
         sumpow=sum(pow),
         waicc=pow/sumpow) %>%
  dplyr::select(c(2,1,3,6)) %>%
  arrange(desc(waicc))


#### Testing a regime of intermediate host habitat ####
int_data<-data.frame(tree2$tip.label)
int_data$regime<-dat$int_hab
int_data$trait<-shape

# estimate the regimes for internal nodes uses "corHMM"
p2<-rayDISC(tree2, int_data[,c(1,2)], model="ER", node.states="marginal")

# run the regime (0 = aquatic , 1 = terrestrial)
OUwie(p2$phy, int_data, model="BM1") #random   == aicc = 204.49, optima = 3.61
OUwie(p2$phy, int_data, model="OU1") #1 regime == aicc = 190.77, optima = 3.64
OUwie(p2$phy, int_data, model="OUM") #2 regime == aicc = 190.16, optima0 = 3.99, optima1 = 1.21
OUwie(p2$phy, int_data, model="OUMV")#2 reg +  == aicc = 177.43, optima0 = 3.99, optima1 = 1.68


###Int host table ###
int_evo<-c(204.49,190.77, 190.16, 177.43)
names(int_evo)<-c("BM", "OU1", "OU2", "OU2+V")

aic_int <- tibble(aicc=int_evo, model=names(int_evo))
int_evo_table<-aic_int %>% 
  mutate(deltaaicc=aicc-min(aicc), 
         pow=exp(-.5*deltaaicc),
         sumpow=sum(pow),
         waicc=pow/sumpow) %>%
  dplyr::select(c(2,1,3,6)) %>%
  arrange(desc(waicc))


#### Testing a regime of definitive and intermediate host habitat combined ####
com_data<-data.frame(tree2$tip.label)
com_data$regime<-dat$comb
com_data$trait<-shape

# estimate the regimes for internal nodes uses "corHMM"
p3<-rayDISC(tree2, com_data[,c(1,2)], model="ER", node.states="marginal")

# run the regime
OUwie(p3$phy, com_data, model="BM1") #random   == aicc = 204.49, optima = 3.61
OUwie(p3$phy, com_data, model="OU1") #1 regime == aicc = 190.77, optima = 3.64
OUwie(p3$phy, com_data, model="OUM") #4 regime == aicc = 198.33, 
OUwie(p3$phy, com_data, model="OUMV")#4 reg +  == aicc = 201.17, 



###Int + Def host table ###
comb_evo<-c(204.49,190.77, 198.33, 201.17)
names(comb_evo)<-c("BM", "OU1", "OU4", "OU4+V")

aic_comb <- tibble(aicc=comb_evo, model=names(comb_evo))
comb_evo_table<-aic_comb %>% 
  mutate(deltaaicc=aicc-min(aicc), 
         pow=exp(-.5*deltaaicc),
         sumpow=sum(pow),
         waicc=pow/sumpow) %>%
  dplyr::select(c(2,1,3,6)) %>%
  arrange(desc(waicc))



##wAICc table

#load data
#waicc_figure <- read_csv("waicc_figure.csv")
#colnames(waicc_figure)<-c("Regime", "Model", "Weight")

#ggplot(waicc_figure)+
#  geom_bar(mapping=aes(x= Regime, y =Weight, color= Model))
