rm(list=ls())
library(plyr)
library(ggplot2)
library(ggpubr) #data visualisation 
library(hrbrthemes)#COLORS GEOM TILE
library(data.table)#Dataset into datatable
library(viridis)
library(rpart)
library(dplyr)
library(tidyr)
library(stringr)

#Load epibiont quantitative dataset----
epi<-read.csv("epi_quant.csv")
epi$Date<-as.Date(epi$Date, "%d/%m/%Y")
epi <- subset(epi,Seaweed.species == "Alaria"|Seaweed.species == "Saccharina")
#epi<-epi[epi_segment$Date < "2022-07-30",]

epi_segment<-subset(epi,Seaweed.species == "Saccharina")
epi_segment<-subset(epi_segment, Segment == "tip"|Segment == "middle"|Segment == "base")
epi_segment<-epi_segment[epi_segment$Date < "2022-07-30",]

#epi<-epi %>% drop_na(blade_area)# this keeps only the info on the 
#epi<-epi %>% drop_na(blade_area)

#phyto_mean<-phyto1 %>% group_by(Sample.no) %>%
  #summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))
#epi<-epi %>% filter(sample.no <57)

#Epibionts by part blade
ggplot(epi_segment, aes(x=factor(Date), y=percent_Hydrozoans, fill=Segment)) +
  geom_boxplot()+
  scale_fill_manual(values=c( "lightgreen","green" ,"darkgreen"))+
  scale_x_discrete(labels=c("7-5-22","22-5-22","29-6-22", "15-7-22","30-7-22" )) +
  theme_test()+
  theme(legend.position = "top", legend.text=element_text(size=14),legend.title=element_text(size=14),
        axis.text.x = element_text(angle = 0, size=12),axis.text.y = element_text(size=12),
        axis.title=element_text(size=12),
        strip.text = element_text(size = 12))+
  labs(x="Date", y="Obelia spp. (Hydrozoa) % blade", fill="Saccharina segment") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    axis.line.x      = element_line(color = "black", size = 0.3),
    axis.line.y      = element_line(color = "black", size = 0.3),
    axis.text.x      = element_text(angle = , hjust = 0.5)
  )

#stats
model1<-lm(percent_Hydrozoans~Segment+factor(Date), data=epi_segment)
anova(model1)

# Blade growth 
ggplot(epi, aes(x=factor(Date), y=blade_area, fill=Seaweed.species)) +
  geom_boxplot()+
  facet_grid(Seaweed.species~., scales="free_y")+
  ylab(expression("Blade Area"~ cm^2))+
  xlab(expression("Date"))+
  scale_fill_manual(values=c( "darkorange", "darkgreen"))+
  scale_x_discrete(labels=c("14-1-22","30-1-22","8-2-22","25-2-22","8-3-22", "21-3-22", "10-4-22","7-5-22","22-5-22","29-6-22", "15-7-22","30-7-22" )) +
  theme_test()+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust=1,size=12),axis.text.y = element_text(size=12),
        axis.title=element_text(size=14),
        strip.text = element_text(size = 17, face = "italic"))

# hydrozoans per blade 
ggplot(epi, aes(x=factor(Date), y=percent_blade_Hydrozoans, fill=Seaweed.species)) +
  geom_boxplot()+
  facet_grid(Seaweed.species~.)+
  ylab(expression("Blade Area"~ cm^2))+
  xlab(expression("Date"))+
  scale_fill_manual(values=c( "darkorange", "darkgreen"))+
  scale_x_discrete(labels=c("14-1-22","30-1-22","8-2-22","25-2-22","8-3-22", "21-3-22", "10-4-22","7-5-22","22-5-22","29-6-22", "15-7-22","30-7-22" )) +
  theme_test()+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust=1,size=12),axis.text.y = element_text(size=12),
        axis.title=element_text(size=14),
        strip.text = element_text(size = 17, face = "italic"))


# Bryozoans per blade 
ggplot(epi, aes(x=factor(Date), y=percent_blade_bryozoans, fill=Seaweed.species)) +
  geom_boxplot()+
  facet_grid(Seaweed.species~.)+
  ylab(expression("Blade Area"~ cm^2))+
  xlab(expression("Date"))+
  scale_fill_manual(values=c( "darkorange", "darkgreen"))+
  scale_x_discrete(labels=c("14-1-22","30-1-22","8-2-22","25-2-22","8-3-22", "21-3-22", "10-4-22","7-5-22","22-5-22","29-6-22", "15-7-22","30-7-22" )) +
  theme_test()+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust=1,size=12),axis.text.y = element_text(size=12),
        axis.title=element_text(size=14),
        strip.text = element_text(size = 17, face = "italic"))

# Tufts per blade 
ggplot(epi, aes(x=factor(Date), y=Tufts_cm2, fill=Seaweed.species)) +
  geom_boxplot()+
  facet_grid(Seaweed.species~.)+
  ylab(expression("Tufts"~ cm^-2))+
  xlab(expression("Date"))+
  scale_fill_manual(values=c( "darkorange", "darkgreen"))+
  scale_x_discrete(labels=c("14-1-22","30-1-22","8-2-22","25-2-22","8-3-22", "21-3-22", "10-4-22","7-5-22","22-5-22","29-6-22", "15-7-22","30-7-22" )) +
  theme_test()+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust=1,size=12),axis.text.y = element_text(size=12),
        axis.title=element_text(size=14),
        strip.text = element_text(size = 17, face = "italic"))




#Load and process eDNA data ---- 
#before you start, remove the column "sequence" from the far right end of the excel dataset and bring it in the front, right after the taxonomic info
otu<-read.csv("Final_edna.csv")
taxonomy<- otu[c(1,7:14)] # subset for taxonomy
taxonomy<-taxonomy %>% mutate(across(1:9, as.factor)) #R recognizes these columns as character and not factors so it will
taxonomy<-as.data.frame(t(taxonomy)) # turns this into dataframe
colnames(taxonomy) <- taxonomy[1,] # assigns otus as column names
samples<- otu[c(1,19:147)] # subset for samples
samples<-as.data.frame(t(samples)) #transpose samples
samples <- tibble::rownames_to_column(samples, "id")
colnames(samples) <- samples[1,]
samples <- samples[-1, ]
#extract metadata information on date, site, exposure state
samples$Site <- str_extract(samples$id, "[^_]+") #extract the site in new column by taking only the part of the sample code that is before the first _
samples$Date <- str_extract(samples$id, "[^_]*$") #extract the year in new column by taking only the part of the sample code that is before the first _
samples$Date <- str_c(samples$Date, "2022", sep="-")
samples<-samples %>% relocate(Date, .after = id)
samples<-samples %>% relocate(Site, .after = id)
samples$Date <- chartr(".", "-", samples$Date)
#samples$Date<-as.Date(samples$Date, "%d-%m-%Y) #turn date variable into a date format
samples$Date<-as.Date(samples$Date, "%d-%m-%Y")
samples$Site <- gsub("SCA", "Scalpay", samples$Site) #fix site names
samples$Site <- gsub("PaB", "Pabay", samples$Site) #fix site names
samples$Site <- gsub("PAB", "Pabay", samples$Site) #fix site names
samples<-samples %>% filter(!str_detect(Site, 'blank')) #delete all rows that are blank samples


#metadata<- otu[ c(1:3) ] 
#species<- otu[ -c(1:3) ]
#colnames(species) <- taxonomy[1,] #assign as column names the tax level  name from the tax info file 
#species<-t(apply(species,1, function(x) tapply(x,colnames(species),sum))) #sum taxa with the same name
#species<-as.data.frame(species)
#names(species)[1] <- "unnasigned"
#colnames(species)  <- paste("OTU",sep="_",colnames(species)) # rename all columns to have the suffix "OTU" for molecular data
#colnames(species) <- gsub("-", "_", colnames(species))
#colnames(species) <- gsub(" ", "_", colnames(species))
#OTU_samples <- cbind(metadata, species)
#OTU_samples<-OTU_samples[colMeans(OTU_samples[3:9580] == 0) <= 0.9] #This will delete the OTUs where more than 90% of the entries are 0.#This means that a the remaining species will be present in at least 30 samples  have at least 30 samples across the two sites that will have


#Load Phyto dataset----
phyto<-read.csv("phyto.csv")
phyto$Date<-as.Date(phyto$Date, "%d/%m/%Y")
phyto$Sample.no<-as.factor(phyto$Sample.no)
phyto1 <- subset(phyto,Site == "Pabay")
phyto_sites<-subset(phyto,Sample.no == "1"| Sample.no == "7")
phyto_mean<-phyto1 %>% group_by(Sample.no) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))

#Load Zoo dataset----
zoo<-read.csv("zoo.csv")
zoo$Date<-as.Date(zoo$Date, "%d/%m/%Y")
zoo$Sample.no<-as.factor(zoo$Sample.no)
zoo1 <- subset(zoo,Site == "Pabay")
zoo_sites<-subset(zoo,Sample.no == "1"| Sample.no == "7")
zoo_mean<-zoo1 %>% group_by(Sample.no) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))



#Combined datasets----
plankton<- merge(phyto_mean,zoo_mean,  by=c("Date"))
plankton_sites<- merge(phyto_sites,zoo_sites,  by=c("SampleCode"))

plankton<-rbind.fill(phyto, zoo)

#Melt for heatmaps the species abundance datasets to allow a column for zooplankton and one for abundance----
#first get rid of any columns that are not needed from the initial dataset
phyto2<-subset((phyto1), select = c(Sample.no, Pseudonitzschia.spp, CHAETOCEROS,RHIZOSOLENIA,CERATIUM, DICTYOCHA, THALASSIOSIRA, Skeletonema.spp, 
                                  Ceratoneis.closterium, Paralia.sulcata, Prorocentrum.micans, PROTOPERIDINIUM, GUINARDIA, 
                                  HETEROCAPSA, ALEXANDRIUM))
zoo2<-subset((zoo1), select = c(Sample.no,Copepods, Appendicularia, Gastropoda, Sabellariidae, Ophiura.larvae, Total_MALACOSTRACA,
                                        Total_HYDROZOA,Total_POLYCHAETAE, total_BALANOID))

# then turn dataset into a table
dt_phyto<-setDT(phyto2)
dt_zoo<-setDT(zoo2)

#then melt the data table
m_phyto<-melt(dt_phyto,id.vars=c("Sample.no"))
m_zoo<-melt(dt_zoo,id.vars=c("Sample.no"))


#renaming the two melted columns
colnames(m_phyto)[2] <- "Phytoplankton"
colnames(m_phyto)[3] <- "Abundance"
m_phyto$Abundance<-as.numeric(m_phyto$Abundance)

colnames(m_zoo)[2] <- "Zooplankton"
colnames(m_zoo)[3] <- "Abundance"
m_zoo$Abundance<-as.numeric(m_zoo$Abundance)

#Heatmap of species abundances----
ggplot(m_phyto, aes(Sample.no, Phytoplankton, fill= log(Abundance+1))) + 
  geom_tile()+
  scale_fill_viridis() +
  guides(fill = guide_colourbar(label = TRUE,
                                ticks = FALSE))+
  theme_classic()+
  ylab("Phytoplankton cells/L")+
  xlab("Sampling Day 15 June 2021 - 10 April 2022")+
  theme(text = element_text(size=14))

ggplot(m_zoo, aes(Sample.no, Zooplankton, fill= log(Abundance+1))) + 
  geom_tile()+
  scale_fill_viridis() +
  guides(fill = guide_colourbar(label = TRUE,
                                ticks = FALSE))+
  theme_classic()+
  ylab(expression("Zooplankton ind.m"^{"-3"}))+
  xlab("Sampling Day 15 June 2021 - 10 April 2022")+
  theme(text = element_text(size=14))

#phyto vs zoo

ggplot() + 
  geom_area(data = plankton, aes(x = Date, y = Total_Phytoplankton), fill = "darkgreen") +
  geom_area(data = plankton, aes(x = Date, y = Total_Zooplankton), fill = "orange") +
  # facet_wrap(~Site.x,ncol = 1,scales="free_y")+
  theme_classic()+
  ylab("Plankton")+
  theme(legend.position="none",text = element_text(size=14))



#Melt for stackplots the species abundance datasets to allow a column for zooplankton and one for abundance----
#first get rid of any columns that are not needed from the initial dataset
phyto3<-subset((phyto_mean), select = c(Date, Skeletonema.spp, Pseudonitzschia.spp,Thalassiosira.rotula,
                                        Chaetoceros.spp,Thalassiosira.spp, Thalassiosira.small,
                                        Ditylum.brightwellii, Guinardia.striata,
                                        Ceratoneis.closterium,Karlodinium.sp,
                                        Chaetoceros.curvisetus,  Thalassiosira.eccentrica,
                                        Paralia.sulcata,Rhizosolenia.styliformis,
                                        Cyclotella.spp, Licmophora,
                                        Rhizosolenia.setigera,Dictyocha.speculum, Alexandrium.spp))
zoo3<-subset((zoo_mean), select = c(Date,Copepod.nauplii,  Copepods,  Bivalvia,
                                    Appendicularia, Gastropoda, Sabellariidae,
                                    Bryozoa.larvae, Polychaetae.larvae,  Fish.egg,
                                    Evadne, Decapoda.larvae, Sea.urchin.larvae,
                                    Ophiura.larvae,  Asteroida.larvae,Ascidian.larvae,
                                    Podon,Obelia,Lizzia.blondina,
                                    Anemone.larvae, Nemertini.larvae, Amphipod))

# then turn dataset into a table
dt_phyto<-setDT(phyto3)
dt_zoo<-setDT(zoo3)

#then melt the data table
m_phyto1<-melt(dt_phyto,id.vars=c("Date"))
m_zoo1<-melt(dt_zoo,id.vars=c("Date"))


#renaming the two melted columns
colnames(m_phyto1)[2] <- "Phytoplankton"
colnames(m_phyto1)[3] <- "Abundance"
m_phyto1$Abundance<-as.numeric(m_phyto1$Abundance)

colnames(m_zoo1)[2] <- "Zooplankton"
colnames(m_zoo1)[3] <- "Abundance"
m_zoo1$Abundance<-as.numeric(m_zoo1$Abundance)

#stacked barplot of multiple species abundances per dates----
ggplot(m_phyto1, aes(fill=Phytoplankton, y=Abundance, x=Date)) + 
  geom_bar(position="fill", stat="identity")+
  theme(legend.position="right",text = element_text(size=14))

ggplot(m_zoo1, aes(fill=Zooplankton, y=Abundance, x=Date)) + 
  geom_bar(position="fill", stat="identity")+
  theme(legend.position="right",text = element_text(size=14))
#phyto vs zoo----
ggplot(data = plankton, aes(x=Date,y=Total_Phytoplankton)) +
  geom_area(fill="#A3A500") +
  geom_line(aes(y=Total_Zooplankton*8),size=1.2, color="red")+ 
  ylab("Phytoplankton cells/L")+
  scale_y_continuous(sec.axis = sec_axis( trans=~./8, name=expression("Zooplankton ind.m"^{"-3"})))+ 
  theme_classic()+
  theme(legend.position="none",text = element_text(size=14))

a<-ggplot(data = phyto1, aes(x = Sample.no, y =Total_Phytoplankton)) + 
  geom_boxplot() +
    theme_classic()+
  xlab("")+
  ylab("Phytoplankton cells/L")+
  theme(legend.position="none",text = element_text(size=14))

b<-ggplot(data = zoo1, aes(x = Sample.no, y =Total_Zooplankton)) + 
  geom_boxplot() +
  theme_classic()+
  ylab(expression("Zooplankton ind.m"^{"-3"}))+
  theme(legend.position="none",text = element_text(size=14))


vjust=c(1,1)
ggarrange(a+ rremove("x.text"),b, heights = c(4, 4),
          labels = c("A", "B"), hjust= 0, vjust = vjust,
          ncol = 1, nrow = 2,align = "v")
  


#scalpay vs pabay
ggplot(data=plankton_sites, aes(x=Site.x, y=Total_Zooplankton)) +
  geom_boxplot()+
    facet_wrap(~Date.x,ncol = 1, scales = "free")+
  ylab("Phytoplankotn cells/L")+
  theme_classic()+
  theme(legend.position = "none",text = element_text(size=14))+scale_fill_brewer(palette="Dark2") 
