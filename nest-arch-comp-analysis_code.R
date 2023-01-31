### Analysis for "Foraging behavior affects nest architecture in a cross-species comparison of ant nests"
## Authors: Sean O'Fallon, Kim Drager, Art Zhao, Andy Suarez, and Noa Pinter-Wollman
## Email Sean O'Fallon at seanofallon@g.ucla.edu with questions on this script and data files

#################################################### START UP ###########################################################33
# clear workspace, load packages
rm(list = ls())
library(igraph)
library(car)
library(lme4)
library(lsmeans)
library(multcomp)
library(dplyr)
library(nlme)
library(ape)
library(caper)
library(phytools)
library(phylolm)
library(MASS)
library(scales)
library(goeveg)
library(gridExtra)
library(knitr)
library(purrr)
library(FSA)
library(ggplot2)
library(grid)
library(lattice)
library(gtable)

####################### Load in data files, set working directory, create/edit objects needed:
################# General
# load nest summaries to reference for species names
nest_summaries<-read.csv("C:/Users/seano/Desktop/Projects/ComparativeAnalysisNests/data-analysis/nest-summaries_upload.csv") 

## load and prep foraging strategies datasheet to be paired with nests
foraging_strats<-read.csv("C:/Users/seano/Desktop/Projects/ComparativeAnalysisNests/data-analysis/foraging-by-species_upload.csv")
# remove excess cols
foraging_strats<-foraging_strats[,c(1,3)]
# fix spelling
foraging_strats[foraging_strats == "Mass Recruitment by Pheromone Trail"] <- "Mass Recruitment"
foraging_strats[foraging_strats == "Mass recruitment by pheromone trail"] <- "Mass Recruitment"
foraging_strats[foraging_strats == "Solitary foraging"] <- "Solitary Foraging"
# Group like categories
foraging_strats[foraging_strats == "Mass recruitment by foraging trail"] <- "Mass Recruitment"
foraging_strats[foraging_strats == "Tandem Running"] <- "Group Recruitment"
foraging_strats[foraging_strats == "Stable Trunk Trail"] <- "Stable Trail"
foraging_strats[foraging_strats == "Long-term Trail Network"] <- "Stable Trail"


### Nest networks for Density and ID entrance chambers
# set working directory to folder w/ edgelist files to create nest networks
setwd("C:/Users/seano/Desktop/Projects/ComparativeAnalysisNests/data-analysis/All Networks")
# make a list of the names of all the files in that folder, now input_files is a list of csv's
input_files<-list.files() 
# saves path to be called in loop (must go to folder w/ network csvs - path must end in "/" to append filename to)
data_folder<-("C:/Users/seano/Desktop/Projects/ComparativeAnalysisNests/data-analysis/All Networks/") 

################# Entrance Chamber Widths, # Chambers
# make empty dataframe to fill with entrance chamber widths
e_chamber_data<-data.frame() 
#loads chamber datasheet to be referenced in loop - change to path to Chamber_Widths_perChamber.csv
chamber_data<-read.csv("C:/Users/seano/Desktop/Projects/ComparativeAnalysisNests/data-analysis/Chamber_Widths_perChamber.csv") 
# clean chamber_data so it only contains chambers
chamber_data<-chamber_data[chamber_data$Structure.Type=="Chamber",]
# load data file w/ widths of entrance chambers for which we don't have a network
networkless_e_chambers<-read.csv("C:/Users/seano/Desktop/Projects/ComparativeAnalysisNests/data-analysis/e-chamber-widths_no-net.csv", header=FALSE)
# load data file w/ ant morphometry measurements - read.csv needs to be called on path to morpho_upload
morph_data<-read.csv("C:/Users/seano/Desktop/Projects/ComparativeAnalysisNests/data-analysis/morpho_upload.csv") #loads morph datasheet to be references in loop


## Extract e. chamber widths from edgelists in All_Networks and widths in Chamber_Widths_perChamber
# loop takes each file in input_files and processes it - representing each nest as network, extracting chamber data from chamber sheet
for (i in 1:length(input_files)){ #For each nest (saved as a csv file) in the folder...
  # get the nest as a network
  path<-paste(data_folder,input_files[i],sep="") #select the filepath for ith nest
  nest<-read.csv(path, header=FALSE) #read this file and convert into a data frame
  edgelist<-as.matrix(nest[c(3,4)]) #in Ant Nest Excel Sheet, the 2-column matrix is in colmns 3 & 4, convert data to a matrix - 2 columns
  graph<-graph_from_edgelist(edgelist,directed=FALSE) #make an igraph object (network) of nest
  # plot(graph) # this will plot the network - uncomment for visualized networks to all pop up
  # based on that network and the chamber data table, pick entrance chamber(s) and record them
  if(nrow(chamber_data[chamber_data$NestID==nest[1,1],])>0){ # checks that ith nest appears in chamber data sheet - skips nest if not
    # Find entrance nodes
    e_nodes<-edgelist[grepl("E", edgelist)]
    ends<-e_nodes[grepl("END", e_nodes)]
    e_nodes<-e_nodes[!e_nodes %in% ends] 
    # return nodes attached to E nodes
    e_chambers<-V(graph)[.nei(e_nodes)] #choose the node(s) connected to entrance as the entrance chambers
    jxns<-e_chambers$name[grepl("J", e_chambers$name)] #IDs nodes connected to entrances that are junctions
    e_chambers<-V(graph)[.nei(c(e_nodes, jxns))] #chambers connected to entrances and chambers connected to junctions directly connected to entrances ID'd as entrance chambers
    # Write entrance chamber widths to e_chamber_data object
    for(j in 1:length(e_chambers)) { #loop for each node ID'd as entrance chamber (few have more than 1)
      if(nrow(chamber_data[chamber_data$Structure.ID==e_chambers$name[j],])>0){ # check that ith entrance chamber is in chamber data sheet - skips nest if not
        e_chamber_data[nrow(e_chamber_data) + 1, 1]<-nest[1,1] #write name of nest in output dataframe (col 1)
        # write species name, name of entrance chamber to output
        if(nrow(nest_summaries[nest_summaries$Nest.ID==nest[1,1],])>0){
          e_chamber_data[nrow(e_chamber_data),2]<-nest_summaries[nest_summaries$Nest.ID==nest[1,1],'Species']
          e_chamber_data[nrow(e_chamber_data),3]<-e_chambers$name[j] #write name of entrance chamber in output dataframe (col 2)
          # write width of entrance chamber to output dataframe (col 4)
          if(sum(chamber_data$NestID==nest$V1[1] & chamber_data$Structure.ID==e_chambers$name[j])>0){
            e_chamber_data[nrow(e_chamber_data), 4]<-chamber_data[chamber_data$NestID==nest$V1[1] &
                                                                    chamber_data$Structure.ID %in% e_chambers$name[j], "Maximum.Width"]
          }
          else{
            e_chamber_data[nrow(e_chamber_data), 4]<-NA
          }
        }
      }
    }
  }
}
# remove na rows
e_chamber_data<-na.omit(e_chamber_data)

# Attach widths that come from non-networked nests
e_chamber_data<-rbind(e_chamber_data, networkless_e_chambers)

# add names to columns
names(e_chamber_data) <- c("Nest.ID","Species","E_chamber_ID","E_chamber_width")

### Scale e. chambers by size of ant
## clean morph_data, make col for max head width
morph_data<-morph_data[,-c(2,3,7)] # get rid of extra columns
morph_data[,'Max.Width']<-NA # make new column for max head size for below loop
morph_data[is.na(morph_data)]<-0 # change NAs to 0 for below loop (note that now NA values should be removed by removing 0s, not NAs)

# make column with max head size (should check head width at eyes vs max if different, report higher)
for (i in 1:nrow(morph_data)) {
  if (morph_data$Max.head.width..if.different.[i]>morph_data$Head.width.at.eyes[i]) {
    morph_data$Max.Width[i]<-morph_data$Max.head.width..if.different.[i]
  }
  else{
    morph_data$Max.Width[i]<-morph_data$Head.width.at.eyes[i]
  }
}

## make object mean_heads for mean head size / species
# first make object w/ just all head widths
head_data<-morph_data[morph_data$Max.Width>0,c('Species.Name','Max.Width')]
# calculate average head sizes, save in mean_heads
mean_heads<-aggregate(head_data$Max.Width, list(head_data$Species.Name), FUN=mean)
# add column names
names(mean_heads)<-c("Species","Mean_Head")

## make object mean_webers for mean webers length / species
# first make object w/ just all webers lengths
webers_data<-morph_data[morph_data$Weber.s.length>0,c('Species.Name','Weber.s.length')]
# calculate average webers lengths, save to mean_webers
mean_webers<-aggregate(webers_data$Weber.s.length, list(webers_data$Species.Name), FUN=mean)
# add column names
names(mean_webers)<-c("Species","Mean_Webers")

# combine mean morpho measures with e. chamber df
e_chamber_data<-merge(mean_heads, e_chamber_data, by="Species")
e_chamber_data<-merge(mean_webers, e_chamber_data, by="Species")

# add empty columns to put in scaled e. chamber widths
e_chamber_data[,'Width_by_Head']<-NA
e_chamber_data[,'Width_by_Webers']<-NA

# loop to divide e. chamber width by head size and webers length for each e. chamber
for(i in 1:nrow(e_chamber_data)){
  e_chamber_data[i,'Width_by_Head']<-e_chamber_data[i,"E_chamber_width"]/e_chamber_data[i,"Mean_Head"]
  e_chamber_data[i,'Width_by_Webers']<-e_chamber_data[i,"E_chamber_width"]/e_chamber_data[i,"Mean_Webers"]
}
# create ec_w_reps w/ individual widths for all species w/ > 2 nests (2 or less are not used in analyses)
ec_w_reps<-e_chamber_data %>%
  group_by(Species) %>%
  filter(n() > 2)
ec_w_reps<-as.data.frame(ec_w_reps)
## For phylogenetic signal calculations, need single value / species
# Make object w/ mean scaled e. chamber values for each species
widths_by_species<-merge(aggregate(ec_w_reps$Width_by_Webers, list(ec_w_reps$Species),mean),
                         aggregate(ec_w_reps$Width_by_Head, list(ec_w_reps$Species),mean),
                         by="Group.1")
# add "genus" col to be paired w/ tree
widths_by_species$genus<-sapply(strsplit(as.character(widths_by_species$Group.1), "\\ "), "[[", 1) 
# add col for number of nests included
for (i in 1:nrow(widths_by_species)) {
  widths_by_species[i,"num.nest"]<-nrow(filter(ec_w_reps, Species==widths_by_species[i,1]))
}
# add cols for coefficient of variance for scaled e. chamber widths / species
widths_by_species<-merge(
  merge(widths_by_species,
        aggregate(ec_w_reps$Width_by_Webers, list(ec_w_reps$Species),cv),by="Group.1"),
  aggregate(ec_w_reps$Width_by_Head, list(ec_w_reps$Species),cv),by="Group.1")
# give each column its name
names(widths_by_species)[1]<-"Species"
names(widths_by_species)[2]<-"webers.width"
names(widths_by_species)[3]<-"width.by.head"
names(widths_by_species)[6]<-"ww.cv"
names(widths_by_species)[7]<-"hw.cv"


################# Density
density_data<-data.frame() # empty df to fill w/ density measurements

# loop takes each file in input_files and processes it - representing each nest as network, extracting density
for (i in 1:length(input_files)){ #For each nest (saved as a csv file) in the folder...
  
  # get the nest as a network
  path<-paste(data_folder,input_files[i],sep="") #select the filepath for ith nest
  nest<-read.csv(path, header=FALSE) #read this file and convert into a data frame
  edgelist<-as.matrix(nest[c(3,4)]) #in Ant Nest Excel Sheet, the 2-column matrix is in colmns 3 & 4, convert data to a matrix - 2 columns
  graph<-graph_from_edgelist(edgelist,directed=FALSE) #make an igraph object (network) of nest
  density_data[nrow(density_data) + 1, 1]<-nest[1,1] #write name of nest in output dataframe (col 1)
  # write species name here! with up-to-date nest summaries
  if(nrow(nest_summaries[nest_summaries$Nest.ID==nest[1,1],])>0){
    density_data[nrow(density_data),2]<-nest_summaries[nest_summaries$Nest.ID==nest[1,1],'Species']
    density_data[nrow(density_data),3]<-edge_density(graph, loops = FALSE)
    # next line puts number of chambers in current network into column in aa_on_density
    density_data[nrow(density_data),4]<-length(unique(nest[grepl("C", nest$V4),4]))
  }
  else{
    density_data[nrow(density_data),3]<-"Not in Nest Summaries!"
  }
}
# name cols
names(density_data)<-c("Nest.ID","Species","Density","Num.Chambers")
# remove nests that don't have density measurements
density_data<-filter(density_data, Density !="Not in Nest Summaries!")

# create object w/ individual widths for all species w/ > 2 nests
nets_w_reps<-density_data %>%
  group_by(Species) %>%
  filter(n() > 2)
nets_w_reps<-as.data.frame(nets_w_reps)
nets_w_reps$Density<-as.numeric(nets_w_reps$Density)

## Make density_by_species w/ mean density for each species - for phylo analysis
density_by_species<-merge(
  aggregate(nets_w_reps$Density, list(nets_w_reps$Species), FUN=mean),
  aggregate(nets_w_reps$Density, list(nets_w_reps$Species), FUN=cv),
  by="Group.1")
# add "genus" column to be paired w/ tree
density_by_species$genus<-sapply(strsplit(as.character(density_by_species$Group.1), "\\ "), "[[", 1) 
# add column for number of nests included
for (i in 1:nrow(density_by_species)) {
  density_by_species[i,"num.nest"]<-nrow(filter(nets_w_reps, Species==density_by_species[i,1]))
}

# Name columns
names(density_by_species)[1]<-"Species"
names(density_by_species)[2]<-"Density"
names(density_by_species)[3]<-"CV"

################# Number of Chambers
# Remove density as we'll include species w/o networks
num.cham_w_reps<-nets_w_reps[,-3]

# Add non-networked nests w/ num.cham known
networkless_num.cham<-read.csv("C:/Users/seano/Desktop/Projects/ComparativeAnalysisNests/data-analysis/num-cham_no-net_upload.csv", header=FALSE)
# name cols
names(networkless_num.cham)<-c("Nest.ID","Species","Num.Chambers")

# Combine non-networked nests and num.cham taken from networks
num.cham_w_reps<-rbind(num.cham_w_reps, networkless_num.cham)

# Atta is an outlier - remove it
num.cham_w_reps<-num.cham_w_reps[!(num.cham_w_reps$Species=="Atta laevigata"),]

## Make density_by_species w/ mean num.cham for each species - for phylo analysis
num.cham_by_species<-merge(
  aggregate(num.cham_w_reps$Num.Chambers, list(num.cham_w_reps$Species), FUN=mean),
  aggregate(num.cham_w_reps$Num.Chambers, list(num.cham_w_reps$Species), FUN=cv),
  by="Group.1")

# add "genus" column to be paired w/ tree
num.cham_by_species$genus<-sapply(strsplit(as.character(num.cham_by_species$Group.1), "\\ "), "[[", 1) 
# add column for number of nests included
for (i in 1:nrow(num.cham_by_species)) {
  num.cham_by_species[i,"num.nest"]<-nrow(filter(num.cham_w_reps, Species==num.cham_by_species[i,1]))
}

# Name columns
names(num.cham_by_species)[1]<-"Species"
names(num.cham_by_species)[2]<-"Num.Cham"
names(num.cham_by_species)[3]<-"CV"



################# Nest Depth
## Prepare object w/ depth for each nest
# filter just rows w/ a measurement
nest_depths<-nest_summaries[nest_summaries$Nest.Depth...m.!="n/a"&
                              nest_summaries$Nest.Depth...m.!="ask"&
                              nest_summaries$Nest.Depth...m.!="n"&
                              nest_summaries$Nest.Depth...m.!="",
                            c(1,2,9)]
nest_depths$Nest.Depth...m.<-as.numeric(nest_depths$Nest.Depth...m.)
names(nest_depths)[3]<-"Nest.Depth"
# filter out species w/ <3 nests 
nest_depths<-nest_depths %>%
  group_by(Species) %>%
  filter(n() > 2)
nest_depths<-as.data.frame(nest_depths)

## Make nest.depth_by_species w/ mean nest depth for each species - for phylo analysis
nest.depth_by_species<-merge(
  aggregate(nest_depths$Nest.Depth, list(nest_depths$Species), FUN=mean),
  aggregate(nest_depths$Nest.Depth, list(nest_depths$Species), FUN=cv),
  by="Group.1")

# add "genus" column to be paired w/ tree
nest.depth_by_species$genus<-sapply(strsplit(as.character(nest.depth_by_species$Group.1), "\\ "), "[[", 1) 
# add column for number of nests included
for (i in 1:nrow(nest.depth_by_species)) {
  nest.depth_by_species[i,"num.nest"]<-nrow(filter(nest_depths, Species==nest.depth_by_species[i,1]))
}

# Name columns
names(nest.depth_by_species)[1]<-"Species"
names(nest.depth_by_species)[2]<-"Nest.Depth"
names(nest.depth_by_species)[3]<-"CV"

################# Nest Depth / Weber's Length
## combine mean webers length with nest_depths df
nest_depths<-merge(mean_webers, nest_depths, by="Species")

# add empty columns to put in scaled nest depth
nest_depths[,'scaled.depth']<-NA

# loop to divide nest depth by webers length for each row
for(i in 1:nrow(nest_depths)){
  nest_depths[i,'scaled.depth']<-nest_depths[i,"Nest.Depth"]/nest_depths[i,"Mean_Webers"]
}

## Make scaled.nest.depth_by_species w/ mean nest depth / WL for each species - for phylo analysis
scaled.nest.depth_by_species<-merge(
  aggregate(nest_depths$scaled.depth, list(nest_depths$Species), FUN=mean),
  aggregate(nest_depths$scaled.depth, list(nest_depths$Species), FUN=cv),
  by="Group.1")

# add "genus" column to be paired w/ tree
scaled.nest.depth_by_species$genus<-sapply(strsplit(as.character(scaled.nest.depth_by_species$Group.1), "\\ "), "[[", 1) 
# add column for number of nests included
for (i in 1:nrow(scaled.nest.depth_by_species)) {
  scaled.nest.depth_by_species[i,"num.nest"]<-nrow(filter(nest_depths, Species==scaled.nest.depth_by_species[i,1]))
}

# Name columns
names(scaled.nest.depth_by_species)[1]<-"Species"
names(scaled.nest.depth_by_species)[2]<-"scaled.depth"
names(scaled.nest.depth_by_species)[3]<-"CV"


####################################### DIFFERENCES AMONG SPECIES ####################################
### Generate big by-species summary of measures
## Combine species summary objects
# list summaries to be joined
summary.dfs<-list(widths_by_species,
                  density_by_species,
                  num.cham_by_species,
                  nest.depth_by_species,
                  scaled.nest.depth_by_species)
# join them
all.features_by_species<-summary.dfs %>% reduce(full_join, by="Species")

# add col to all.features_by_species w/ num.nest - max num.nest for that species of all features
all_num.nest<-all.features_by_species[,c("num.nest.x","num.nest.y","num.nest.x.x","num.nest.y.y","num.nest")]
all.features_by_species$max.num.nest<-apply(all_num.nest,1,max,na.rm=T)

# Remove redundant/extraneous columns
species.summary_table<-all.features_by_species[,c("Species",
                                                  "webers.width","ww.cv","width.by.head","hw.cv",
                                                  "Density","CV.x",
                                                  "Num.Cham","CV.y",
                                                  "Nest.Depth","CV.x.x",
                                                  "scaled.depth","CV.y.y",
                                                  "max.num.nest")]
# Name cols properly
names(species.summary_table)<-c("Species",
                                "ECW / Weber's Length","ECW / Weber's Length CV",
                                "ECW / Head Width","ECW / Head Width CV",
                                "Network Density","Network Density CV",
                                "Number of Chambers","Number of Chambers CV",
                                "Nest Depth","Nest Depth CV",
                                "Scaled Nest Depth","Scaled Nest Depth CV",
                                "Number of Nests")
# rearrange order of cols for consistency
species.summary_table<-species.summary_table[,c(1,14,10,11,12,13,8,9,6,7,2,3,4,5)]
# round numbers to 3 digits
species.summary_table<-cbind(species.summary_table[,1], round(species.summary_table[,-1], digits = 3))
names(species.summary_table)[1]<-"Species"

# make into Gtable to display
species.summary_table.print<-tableGrob(species.summary_table,
                                       theme=ttheme_minimal(),
                                       rows=NULL)
species.summary_table.print <- gtable_add_grob(species.summary_table.print,
                                               grobs = rectGrob(gp = gpar(fill = NA, lwd = 1)),
                                               t = 1, b = nrow(species.summary_table.print), l = 1, r = ncol(species.summary_table.print))
species.summary_table.print <- gtable_add_grob(species.summary_table.print,
                                               grobs = rectGrob(gp = gpar(fill = NA, lwd = 1)),
                                               t = 2, b = nrow(species.summary_table.print), l = 1, r = ncol(species.summary_table.print))

# Kruskal-Wallis on nest.depth by species
nest.depth_kw <- kruskal.test(Nest.Depth ~ Species,
                              data = nest_depths)

# Kruskal-Wallis on nest depth / weber's length by species
scaled.nest.depth_kw <- kruskal.test(scaled.depth ~ Species,
                                     data = nest_depths)

# Kruskal-Wallis on num.chambers by species
num.cham_kw <- kruskal.test(Num.Chambers ~ Species,
                            data = num.cham_w_reps)

# Kruskal-Wallis on density
density_kw <- kruskal.test(Density ~ Species,
                           data = nets_w_reps)

# ANOVA on widths scaled by head width
hw_ec_anova <- aov(log(Width_by_Head) ~ Species,
                   data = ec_w_reps)

# ANOVA on widths scaled by weber's length
ww_ec_anova <- aov(log(Width_by_Webers) ~ Species,
                   data = ec_w_reps)

### Report tests of differences among species
## Make table for features tested by ANOVA
hw.aov_stats<-c(summary(hw_ec_anova)[[1]][["Df"]][[1]],
                summary(hw_ec_anova)[[1]][["Sum Sq"]][[1]],
                summary(hw_ec_anova)[[1]][["F value"]][[1]],
                summary(hw_ec_anova)[[1]][["Pr(>F)"]][[1]])

ww.aov_stats<-c(summary(ww_ec_anova)[[1]][["Df"]][[1]],
                summary(ww_ec_anova)[[1]][["Sum Sq"]][[1]],
                summary(ww_ec_anova)[[1]][["F value"]][[1]],
                summary(ww_ec_anova)[[1]][["Pr(>F)"]][[1]])

aov_stats<-rbind(hw.aov_stats, ww.aov_stats)

aov_stats<-as.data.frame(aov_stats)

aov_stats$feature.name<-c("ECW / Head Width","ECW / Weber's Length")

names(aov_stats)<-c("Df","Sum Sq","F value","p-value","Feature Name")

# round to 3 digits
aov_stats$`Sum Sq` <- round(aov_stats$`Sum Sq`, digits = 3)
aov_stats$`F value` <- round(aov_stats$`F value`, digits = 3)
aov_stats$`p-value` <- round(aov_stats$`p-value`, digits = 3)

# replace really small numbers with "<0.0001"
for (i in 1:nrow(aov_stats)) {
  if (as.numeric(aov_stats$`p-value`[i]) < 0.0001) {
    aov_stats$`p-value`[i] <- "<0.0001"
  }
}

aov_stats<-aov_stats[,c(5,1,2,3,4)]

aov.species_table<-tableGrob(aov_stats,
                             theme=ttheme_minimal(),
                             rows=NULL)

aov.species_table <- gtable_add_grob(aov.species_table,
                                     grobs = rectGrob(gp = gpar(fill = NA, lwd = 1)),
                                     t = 1, b = nrow(aov.species_table), l = 1, r = ncol(aov.species_table))

aov.species_table <- gtable_add_grob(aov.species_table,
                                     grobs = rectGrob(gp = gpar(fill = NA, lwd = 1)),
                                     t = 2, b = nrow(aov.species_table), l = 1, r = ncol(aov.species_table))

## Make table for features tested w/ K-W
kw_stats<-as.data.frame(rbind(c(density_kw$parameter,density_kw$statistic,density_kw$p.value),
                              c(num.cham_kw$parameter,num.cham_kw$statistic,num.cham_kw$p.value),
                              c(nest.depth_kw$parameter,nest.depth_kw$statistic,nest.depth_kw$p.value),
                              c(scaled.nest.depth_kw$parameter,scaled.nest.depth_kw$statistic,scaled.nest.depth_kw$p.value)))

kw_stats<-as.data.frame(rbind(c(nest.depth_kw$parameter,nest.depth_kw$statistic,nest.depth_kw$p.value),
                              c(scaled.nest.depth_kw$parameter,scaled.nest.depth_kw$statistic,scaled.nest.depth_kw$p.value),
                              c(num.cham_kw$parameter,num.cham_kw$statistic,num.cham_kw$p.value),
                              c(density_kw$parameter,density_kw$statistic,density_kw$p.value)))

# make col w/ feature names
kw_stats$feature.name<-c("Nest Depth","Scaled Nest Depth","Number of Chambers","Network Density")

# name col w/ p-values
names(kw_stats)[3]<-"p-value"
names(kw_stats)[4]<-"Feature Name"

# round KW chi-squareds to 3 digits
kw_stats$`Kruskal-Wallis chi-squared`<-round(kw_stats$`Kruskal-Wallis chi-squared`, digits = 3)

# round p-values to 3 digits, replace tiny numbers w/ "<0.0001"
for (i in 1:nrow(kw_stats)) {
  if (as.numeric(kw_stats$'p-value'[i]) < 0.0001) {
    kw_stats$'p-value'[i] <- "<0.0001"
  }
  else{
    kw_stats$'p-value'[i] <- round(as.numeric(kw_stats$'p-value'[i]), digits = 3)
  }
}

kw_stats<-kw_stats[,c(4,1,2,3)]

kw.species_table<-tableGrob(kw_stats,
                            theme=ttheme_minimal(),
                            rows=NULL)

kw.species_table <- gtable_add_grob(kw.species_table,
                                    grobs = rectGrob(gp = gpar(fill = NA, lwd = 1)),
                                    t = 1, b = nrow(kw.species_table), l = 1, r = ncol(kw.species_table))

kw.species_table <- gtable_add_grob(kw.species_table,
                                    grobs = rectGrob(gp = gpar(fill = NA, lwd = 1)),
                                    t = 2, b = nrow(kw.species_table), l = 1, r = ncol(kw.species_table))


##################################### PHYLOGENETIC SIGNAL / HEATMAP ###################################

######################### ECWs
## Load tree
# Load tree from Blanchard and Moreau 2017 - https://doi.org/10.1111/evo.13117, Fig S1
# read.nexus() must be called on file path
ant.tree.moreau<-read.nexus("C:/Users/seano/Desktop/Projects/Milleretal2021/Miller_mats/MoreauTree2016")

## Prune tree to the genera represented in ECW analysis
# tips to remain
tips<-c("Acromyrmex_versicolor",
        "Aphaenogaster_occidentalis_NW",
        "Camponotus_maritimus",
        "Diacamma_rugosum",
        "Dinoponera_australis",
        "Dorymyrmex_bicolor",
        "Ectatomma_opaciventre",
        "Mycetagroicus_triangularis",
        "Mycetarotes_acutus",
        "Mycetophylax_conformis",
        "Odontomachus_coquereli",
        "Pachycondyla_harpax",
        "Pheidole_longispinosa",
        "Sericomyrmex_Sp")
# new tree with only represented genera
pruned.tree.widths<-drop.tip(ant.tree.moreau,ant.tree.moreau$tip.label[-match(tips, ant.tree.moreau$tip.label)])

## Rename tip labels to create genus level tree
for (i in 1:length(pruned.tree.widths$tip.label)) {
  split.tips<-strsplit(pruned.tree.widths$tip.label[i],"_")
  genus.tip<-split.tips[[1]]
  pruned.tree.widths$tip.label[i]<-genus.tip[1]
}

## "Main" subset for illustrative purposes will contain representative of genus w/ most nests (1st by alph if ties)
width_rep_subset_main<-subset(widths_by_species, Species!="Aphaenogaster ashmeadi" &
                                Species!="Aphaenogaster treatae" &
                                Species!="Dorymyrmex bureni" &
                                Species!="Ectatomma brunneum" &
                                Species!="Ectatomma edentatum" &
                                Species!="Ectatomma opaciventre" &
                                Species!="Mycetarotes acutus")
# rownames need to match tree tips (genus)
rownames(width_rep_subset_main)<-width_rep_subset_main$genus

## Save species and ECW cols w/ genus rownames for heatmap
width.for.heatmap<-width_rep_subset_main[,-c(1,4,5)]

## Test main subset for phylosig
wl.trait<-setNames(width_rep_subset_main[,2],rownames(width_rep_subset_main))
hw.trait<-setNames(width_rep_subset_main[,3],rownames(width_rep_subset_main))
## test for phylogenetic signal in ECW/WL
# w/ lambda
main_wl.lambda<-phylosig(pruned.tree.widths, wl.trait, method = "lambda", test=T)
# w/ K
main_ww.K<-phylosig(pruned.tree.widths, wl.trait, method = "K", test=T)
## test for phylogenetic signal in ECW/HW
# w/ lambda
main_hw.lambda<-phylosig(pruned.tree.widths, hw.trait, method = "lambda", test=T)
# w/ K
main_hw.K<-phylosig(pruned.tree.widths, hw.trait, method = "K", test=T)

## Write results
# make df to take ECW/WL results
main_ww.res<-data.frame(feature.name = character(),
                        lambda=numeric(),
                        lambda.p=numeric(),
                        K=numeric(),
                        K.p=numeric())
# write
main_ww.res[1,"feature.name"]<-"ECW / WL"
main_ww.res[1,"lambda"]<-main_wl.lambda$lambda
main_ww.res[1,"lambda.p"]<-main_wl.lambda$P
main_ww.res[1,"K"]<-main_ww.K$K
main_ww.res[1,"K.p"]<-main_ww.K$P

# make df to take ECW/HW results
main_hw.res<-data.frame(feature.name=character(),
                        lambda=numeric(),
                        lambda.p=numeric(),
                        K=numeric(),
                        K.p=numeric())
# write
main_hw.res[1,"feature.name"]<-"ECW / HW"
main_hw.res[1,"lambda"]<-main_hw.lambda$lambda
main_hw.res[1,"lambda.p"]<-main_hw.lambda$P
main_hw.res[1,"K"]<-main_hw.K$K
main_hw.res[1,"K.p"]<-main_hw.K$P

### Get phylogenetic signal for every combination of species (1 / genus)
## Make object with 1 row for each combination of species
# make object with just each species and its genus
sp_gn<-widths_by_species %>% dplyr::select(Species,genus)

# make df w/ species+genus for repeated genera
# empty df to hold app rows
rep_sp_gn<-data.frame(Species=character(),
                      genus=character())
# loop fills ^obj w/ rows whose "genus" is in sp_gn >1x
for (i in 1:nrow(sp_gn)) {
  if (nrow(sp_gn[sp_gn$genus==sp_gn[i,2],])>1) {
    rep_sp_gn<-rep_sp_gn %>% add_row(sp_gn[i,]) 
  }
}
# make df w/ species+genus for species in every combo (only one species per genus)
# empty df to hold app rows
single_sp_gn<-data.frame(Species=character(),
                         genus=character())
# loop fills ^ obj w/ rows whose "genus" is in sp_gn 1x
for (i in 1:nrow(sp_gn)) {
  if (nrow(sp_gn[sp_gn$genus==sp_gn[i,2],])==1) {
    single_sp_gn<-single_sp_gn %>% add_row(sp_gn[i,]) 
  }
}

tbl=table(rep_sp_gn$genus)
all_sp_poss = matrix(NA, ncol=length(unique(rep_sp_gn$genus)))
unq_gn=unique(rep_sp_gn$genus)

for (i in 1:tbl[1]){
  for (j in 1:tbl[2]){
    for (k in 1:tbl[3]){
      for (m in 1:tbl[4]){
        gni=unq_gn[1]
        spi=rep_sp_gn$Species[rep_sp_gn$genus==gni][i]
        gnj=unq_gn[2]
        spj=rep_sp_gn$Species[rep_sp_gn$genus==gnj][j]
        gnk=unq_gn[3]
        spk=rep_sp_gn$Species[rep_sp_gn$genus==gnk][k]
        gnm=unq_gn[4]
        spm=rep_sp_gn$Species[rep_sp_gn$genus==gnm][m]
        all_sp_poss=rbind(all_sp_poss,c(spi,spj,spk,spm))
      } 
    }
  }
}
#remove NA row at top
all_sp_poss<-all_sp_poss[-1,]
#make empty matrix to be filled w/ all combos
all_width_subs<-matrix(data=NA,nrow=48,ncol=14)
# combine single rep species (single_sp_gn) w/ all_sp_poss to have matrix w/ each combo (all_width_subs)
for (i in 1:nrow(all_sp_poss)) {
  all_width_subs[i,]<-c(all_sp_poss[i,], single_sp_gn[,1])
}

## Use each row in a phylosig call to get lambda and K
# make empty object to fill w/ results
widths.res<-data.frame(ww.lambda=numeric(),
                       ww.lambda.p=numeric(),
                       ww.K=numeric(),
                       ww.K.p=numeric(),
                       hw.lambda=numeric(),
                       hw.lambda.p=numeric(),
                       hw.K=numeric(),
                       hw.K.p=numeric())

for (i in 1:nrow(all_width_subs)) {
  ## save ith subset to df
  width_subi<-widths_by_species %>% filter(Species %in% all_width_subs[i,])
  rownames(width_subi)<-width_subi$genus 
  ## make objects for phylosig call
  wl.trait<-setNames(width_subi[,2],rownames(width_subi)) 
  hw.trait<-setNames(width_subi[,3],rownames(width_subi)) 
  ## test for phylogenetic signal
  # w/ lambda
  ww.lambda<-phylosig(pruned.tree.widths, wl.trait, method = "lambda", test=T)
  hw.lambda<-phylosig(pruned.tree.widths, hw.trait, method = "lambda", test=T)
  # w/ K
  ww.K<-phylosig(pruned.tree.widths, wl.trait, method = "K", test=T)
  hw.K<-phylosig(pruned.tree.widths, hw.trait, method = "K", test=T)
  ## write to output object widths.res
  widths.res[i,"ww.lambda"]<-ww.lambda$lambda
  widths.res[i,"ww.lambda.p"]<-ww.lambda$P
  widths.res[i,"ww.K"]<-ww.K$K
  widths.res[i,"ww.K.p"]<-ww.K$P
  widths.res[i,"hw.lambda"]<-hw.lambda$lambda
  widths.res[i,"hw.lambda.p"]<-hw.lambda$P
  widths.res[i,"hw.K"]<-hw.K$K
  widths.res[i,"hw.K.p"]<-hw.K$P
}

######################### Density
## Re-prune tree to the genera in Density analysis
tips<-c("Acromyrmex_versicolor",
        "Aphaenogaster_occidentalis_NW",
        "Camponotus_maritimus",
        "Diacamma_rugosum",
        "Dorymyrmex_bicolor",
        "Ectatomma_opaciventre",
        "Formica_moki",
        "Mycetagroicus_triangularis",
        "Mycetarotes_acutus",
        "Mycetophylax_conformis",
        "Odontomachus_coquereli",
        "Pachycondyla_harpax",
        "Pheidole_longispinosa",
        "Prenolepis_imparis",
        "Sericomyrmex_Sp")
pruned.tree.density<-drop.tip(ant.tree.moreau,ant.tree.moreau$tip.label[-match(tips, ant.tree.moreau$tip.label)])

# Rename tip labels to create genus level tree
for (i in 1:length(pruned.tree.density$tip.label)) {
  split.tips<-strsplit(pruned.tree.density$tip.label[i],"_")
  genus.tip<-split.tips[[1]]
  pruned.tree.density$tip.label[i]<-genus.tip[1]
}

# "Main" subset for illustrative purposes will contain representative of genus w/ most nests (1st by alph if ties)
density_subset_main<-subset(density_by_species, Species!="Aphaenogaster ashmeadi" &
                              Species!="Aphaenogaster treatae" &
                              Species!="Dorymyrmex bureni" &
                              Species!="Ectatomma ruidum" &
                              Species!="Formica archboldi" &
                              Species!="Mycetarotes acutus")
rownames(density_subset_main)<-density_subset_main$genus 

## Save species and Density cols w/ genus rownames for heatmap
density.for.heatmap<-density_subset_main[,-c(1,4,5)]

## Test main subset for phylosig
density.trait<-setNames(density_subset_main[,2],rownames(density_subset_main)) 
## test for phylogenetic signal
# w/ lambda
main_density.lambda<-phylosig(pruned.tree.density, density.trait, method = "lambda", test=T)
# w/ K
main_density.K<-phylosig(pruned.tree.density, density.trait, method = "K", test=T)

## Save phylogenetic signal res to main_density.res to be combined w/ other features later
# start df
main_density.res<-data.frame(feature.name = character(),
                             lambda=numeric(),
                             lambda.p=numeric(),
                             K=numeric(),
                             K.p=numeric())
# Write results
main_density.res[1,"feature.name"]<-"Density"
main_density.res[1,"lambda"]<-main_density.lambda$lambda
main_density.res[1,"lambda.p"]<-main_density.lambda$P
main_density.res[1,"K"]<-main_density.K$K
main_density.res[1,"K.p"]<-main_density.K$P

### Get phylogenetic signal for every combination of species (1 / genus)
## Make object with 1 row for each combination of species
# make object with just each species and its genus
sp_gn<-density_by_species %>% dplyr::select(Species,genus)

# make df w/ species+genus for repeated genera
# empty df to hold app rows
rep_sp_gn<-data.frame(Species=character(),
                      genus=character())
# loop fills ^obj w/ rows whose "genus" is in sp_gn >1x
for (i in 1:nrow(sp_gn)) {
  if (nrow(sp_gn[sp_gn$genus==sp_gn[i,2],])>1) {
    rep_sp_gn<-rep_sp_gn %>% add_row(sp_gn[i,]) 
  }
}
# make df w/ species+genus for species in every combo (only one species per genus)
# empty df to hold app rows
single_sp_gn<-data.frame(Species=character(),
                         genus=character())
# loop fills ^ obj w/ rows whose "genus" is in sp_gn 1x
for (i in 1:nrow(sp_gn)) {
  if (nrow(sp_gn[sp_gn$genus==sp_gn[i,2],])==1) {
    single_sp_gn<-single_sp_gn %>% add_row(sp_gn[i,]) 
  }
}

tbl=table(rep_sp_gn$genus)
all_sp_poss = matrix(NA, ncol=length(unique(rep_sp_gn$genus)))
unq_gn=unique(rep_sp_gn$genus)

for (i in 1:tbl[1]){
  for (j in 1:tbl[2]){
    for (k in 1:tbl[3]){
      for (m in 1:tbl[4]){
        for (n in 1:tbl[5]) {
          gni=unq_gn[1]
          spi=rep_sp_gn$Species[rep_sp_gn$genus==gni][i]
          gnj=unq_gn[2]
          spj=rep_sp_gn$Species[rep_sp_gn$genus==gnj][j]
          gnk=unq_gn[3]
          spk=rep_sp_gn$Species[rep_sp_gn$genus==gnk][k]
          gnm=unq_gn[4]
          spm=rep_sp_gn$Species[rep_sp_gn$genus==gnm][m]
          gnn=unq_gn[5]
          spn=rep_sp_gn$Species[rep_sp_gn$genus==gnn][n]
          all_sp_poss=rbind(all_sp_poss,c(spi,spj,spk,spm,spn))
        } 
      }
    }
  }
}
#remove NA row at top
all_sp_poss<-all_sp_poss[-1,]
#make empty matrix to be filled w/ all combos
all_density_subs<-matrix(data=NA,nrow=48,ncol=15)
# combine single rep species (single_sp_gn) w/ all_sp_poss to have matrix w/ each combo (all_aa_subs)
for (i in 1:nrow(all_sp_poss)) {
  all_density_subs[i,]<-c(all_sp_poss[i,], single_sp_gn[,1])
}

## Use each row in a phylosig call to get results for every combination of rep species
# Start df
density.res<-data.frame(density.lambda=numeric(),
                        density.lambda.p=numeric(),
                        density.K=numeric(),
                        density.K.p=numeric())
# Write results for each subset
for (i in 1:nrow(all_density_subs)) {
  ## save ith subset to df
  density_subi<-density_by_species %>% filter(Species %in% all_density_subs[i,])
  rownames(density_subi)<-density_subi$genus 
  ## make object for phylosig call
  density.trait<-setNames(density_subi[,2],rownames(density_subi)) 
  ## test for phylogenetic signal
  # w/ lambda
  density.lambda<-phylosig(pruned.tree.density, density.trait, method = "lambda", test=T)
  # w/ K
  density.K<-phylosig(pruned.tree.density, density.trait, method = "K", test=T)
  ## write to output object density.res
  density.res[i,"density.lambda"]<-density.lambda$lambda
  density.res[i,"density.lambda.p"]<-density.lambda$P
  density.res[i,"density.K"]<-density.K$K
  density.res[i,"density.K.p"]<-density.K$P
}

#################### Number of Chambers
# Re-prune tree to the genera in num.cham analysis
tips<-c("Acromyrmex_versicolor",
        "Aphaenogaster_occidentalis_NW",
        "Camponotus_maritimus",
        "Diacamma_rugosum",
        "Dinoponera_australis",
        "Dorymyrmex_bicolor",
        "Ectatomma_opaciventre",
        "Formica_moki",
        "Mycetagroicus_triangularis",
        "Mycetarotes_acutus",
        "Mycetophylax_conformis",
        "Mycocepurus_goeldii",
        "Odontomachus_coquereli",
        "Pachycondyla_harpax",
        "Pheidole_longispinosa",
        "Prenolepis_imparis",
        "Sericomyrmex_Sp")
pruned.tree.num.cham<-drop.tip(ant.tree.moreau,ant.tree.moreau$tip.label[-match(tips, ant.tree.moreau$tip.label)])

#Rename tip labels to create genus level tree
for (i in 1:length(pruned.tree.num.cham$tip.label)) {
  split.tips<-strsplit(pruned.tree.num.cham$tip.label[i],"_")
  genus.tip<-split.tips[[1]]
  pruned.tree.num.cham$tip.label[i]<-genus.tip[1]
}

# "Main" subset for illustrative purposes will contain representative of genus w/ most nests (1st by alph if ties)
num.cham_subset_main<-subset(num.cham_by_species,
                             Species!="Acromyrmex landolti" &
                               Species!="Acromyrmex rugosus" &
                               Species!="Aphaenogaster ashmeadi" &
                               Species!="Aphaenogaster treatae" &
                               Species!="Dorymyrmex bureni" &
                               Species!="Ectatomma edentatum" &
                               Species!="Ectatomma opaciventre" &
                               Species!="Ectatomma ruidum" &
                               Species!="Formica archboldi" &
                               Species!="Mycetarotes acutus" &
                               Species!="Mycocepurus goeldii")

rownames(num.cham_subset_main)<-num.cham_subset_main$genus 

## Save species and num.cham cols w/ genus rownames for heatmap
num.cham.for.heatmap<-num.cham_subset_main[,-c(1,4,5)]

## Test main subset for phylosig
num.cham.trait<-setNames(num.cham_subset_main[,2],rownames(num.cham_subset_main)) 
## test for phylogenetic signal
# w/ lambda
main_num.cham.lambda<-phylosig(pruned.tree.num.cham, num.cham.trait, method = "lambda", test=T)
# w/ K
main_num.cham.K<-phylosig(pruned.tree.num.cham, num.cham.trait, method = "K", test=T)

## Save main subset results to be combined w/ other features later
main_num.cham.res<-data.frame(feature.name=character(),
                              lambda=numeric(),
                              lambda.p=numeric(),
                              K=numeric(),
                              K.p=numeric())
# enter res into cols
main_num.cham.res[1,"feature.name"]<-"# Chambers"
main_num.cham.res[1,"lambda"]<-main_num.cham.lambda$lambda
main_num.cham.res[1,"lambda.p"]<-main_num.cham.lambda$P
main_num.cham.res[1,"K"]<-main_num.cham.K$K
main_num.cham.res[1,"K.p"]<-main_num.cham.K$P

### Get phylogenetic signal for every combination of species (1 / genus)
## Make object with 1 row for each combination of species
# make object with just each species and its genus
sp_gn<-num.cham_by_species %>% dplyr::select(Species,genus)

# make df w/ species+genus for repeated genera
# empty df to hold app rows
rep_sp_gn<-data.frame(Species=character(),
                      genus=character())
# loop fills ^obj w/ rows whose "genus" is in sp_gn >1x
for (i in 1:nrow(sp_gn)) {
  if (nrow(sp_gn[sp_gn$genus==sp_gn[i,2],])>1) {
    rep_sp_gn<-rep_sp_gn %>% add_row(sp_gn[i,]) 
  }
}
# make df w/ species+genus for species in every combo (only one species per genus)
# empty df to hold app rows
single_sp_gn<-data.frame(Species=character(),
                         genus=character())
# loop fills ^ obj w/ rows whose "genus" is in sp_gn 1x
for (i in 1:nrow(sp_gn)) {
  if (nrow(sp_gn[sp_gn$genus==sp_gn[i,2],])==1) {
    single_sp_gn<-single_sp_gn %>% add_row(sp_gn[i,]) 
  }
}
# make objects needed to generate every combination of subsets
tbl=table(rep_sp_gn$genus)
all_sp_poss = matrix(NA, ncol=length(unique(rep_sp_gn$genus)))
unq_gn=unique(rep_sp_gn$genus)
# save every combination to all_sp_poss
for (i in 1:tbl[1]){
  for (j in 1:tbl[2]){
    for (k in 1:tbl[3]){
      for (m in 1:tbl[4]){
        for (n in 1:tbl[5]){
          for (o in 1:tbl[6]){
            for (p in 1:tbl[7]) {
              gni=unq_gn[1]
              spi=rep_sp_gn$Species[rep_sp_gn$genus==gni][i]
              gnj=unq_gn[2]
              spj=rep_sp_gn$Species[rep_sp_gn$genus==gnj][j]
              gnk=unq_gn[3]
              spk=rep_sp_gn$Species[rep_sp_gn$genus==gnk][k]
              gnm=unq_gn[4]
              spm=rep_sp_gn$Species[rep_sp_gn$genus==gnm][m]
              gnn=unq_gn[5]
              spn=rep_sp_gn$Species[rep_sp_gn$genus==gnn][n]
              gno=unq_gn[6]
              spo=rep_sp_gn$Species[rep_sp_gn$genus==gno][o]
              gnp=unq_gn[7]
              spp=rep_sp_gn$Species[rep_sp_gn$genus==gnp][p]
              all_sp_poss=rbind(all_sp_poss,c(spi,spj,spk,spm,spn,spo,spp))
            }
          }
        }
      } 
    }
  }
}
#remove NA row at top
all_sp_poss<-all_sp_poss[-1,]
#make empty matrix to be filled w/ all combos
all_num.cham_subs<-matrix(data=NA,nrow=576,ncol=17)
# combine single rep species (single_sp_gn) w/ all_sp_poss to have matrix w/ each combo (all_num.cham_subs)
for (i in 1:nrow(all_sp_poss)) {
  all_num.cham_subs[i,]<-c(all_sp_poss[i,], single_sp_gn[,1])
}

## Use each row in a phylosig call to get lambda and K
# make df to fill w/ results
num.cham.res<-data.frame(num.cham.lambda=numeric(),
                         num.cham.lambda.p=numeric(),
                         num.cham.K=numeric(),
                         num.cham.K.p=numeric())
# save results from each combination of representative species
for (i in 1:nrow(all_num.cham_subs)) {
  ## save ith subset to df
  num.cham_subi<-num.cham_by_species %>% filter(Species %in% all_num.cham_subs[i,])
  rownames(num.cham_subi)<-num.cham_subi$genus 
  ## make object for phylosig call
  num.cham.trait<-setNames(num.cham_subi[,2],rownames(num.cham_subi)) 
  ## test for phylogenetic signal
  # w/ lambda
  num.cham.lambda<-phylosig(pruned.tree.num.cham, num.cham.trait, method = "lambda", test=T)
  # w/ K
  num.cham.K<-phylosig(pruned.tree.num.cham, num.cham.trait, method = "K", test=T)
  ## write to output object num.cham.res
  num.cham.res[i,"num.cham.lambda"]<-num.cham.lambda$lambda
  num.cham.res[i,"num.cham.lambda.p"]<-num.cham.lambda$P
  num.cham.res[i,"num.cham.K"]<-num.cham.K$K
  num.cham.res[i,"num.cham.K.p"]<-num.cham.K$P
}

#################### Nest Depth
# Re-prune tree to the genera in nest depth analysis
tips<-c("Acromyrmex_versicolor",
        "Aphaenogaster_occidentalis_NW",
        "Atta_texana",
        "Camponotus_maritimus",
        "Dinoponera_australis",
        "Dorymyrmex_bicolor",
        "Ectatomma_opaciventre",
        "Formica_moki",
        "Mycetagroicus_triangularis",
        "Mycetarotes_acutus",
        "Mycetophylax_conformis",
        "Odontomachus_coquereli",
        "Pachycondyla_harpax",
        "Prenolepis_imparis")
pruned.tree.nest.depth<-drop.tip(ant.tree.moreau,ant.tree.moreau$tip.label[-match(tips, ant.tree.moreau$tip.label)])

#Rename tip labels to create genus level tree
for (i in 1:length(pruned.tree.nest.depth$tip.label)) {
  split.tips<-strsplit(pruned.tree.nest.depth$tip.label[i],"_")
  genus.tip<-split.tips[[1]]
  pruned.tree.nest.depth$tip.label[i]<-genus.tip[1]
}

# "Main" subset for illustrative purposes will contain representative of genus w/ most nests (1st by alph if ties)
nest.depth_subset_main<-subset(nest.depth_by_species,
                               Species!="Acromyrmex landolti" &
                                 Species!="Acromyrmex rugosus" &
                                 Species!="Aphaenogaster ashmeadi" &
                                 Species!="Aphaenogaster treatae" &
                                 Species!="Atta laevigata" &
                                 Species!="Dorymyrmex bureni" &
                                 Species!="Ectatomma opaciventre" &
                                 Species!="Formica archboldi" &
                                 Species!="Mycetarotes acutus" &
                                 Species!="Odontomachus chelifer")

rownames(nest.depth_subset_main)<-nest.depth_subset_main$genus 

## Save species and nest depth cols w/ genus rownames for heatmap
depth.for.heatmap<-nest.depth_subset_main[,-c(1,4,5)]

## Test main subset for phylosig
nest.depth.trait<-setNames(nest.depth_subset_main[,2],rownames(nest.depth_subset_main)) 
## test for phylogenetic signal
# w/ lambda
main_nest.depth.lambda<-phylosig(pruned.tree.nest.depth, nest.depth.trait, method = "lambda", test=T)
# w/ K
main_nest.depth.K<-phylosig(pruned.tree.nest.depth, nest.depth.trait, method = "K", test=T)

## Save main subset results to be combined w/ other features later
main_nest.depth.res<-data.frame(feature.name=character(),
                                lambda=numeric(),
                                lambda.p=numeric(),
                                K=numeric(),
                                K.p=numeric())
# write res to each col
main_nest.depth.res[1,"feature.name"]<-"Nest Depth"
main_nest.depth.res[1,"lambda"]<-main_nest.depth.lambda$lambda
main_nest.depth.res[1,"lambda.p"]<-main_nest.depth.lambda$P
main_nest.depth.res[1,"K"]<-main_nest.depth.K$K
main_nest.depth.res[1,"K.p"]<-main_nest.depth.K$P

### Get phylogenetic signal for every combination of species (1 / genus)
## Make object with 1 row for each combination of species
# make object with just each species and its genus
sp_gn<-nest.depth_by_species %>% dplyr::select(Species,genus)

# make df w/ species+genus for repeated genera
# empty df to hold app rows
rep_sp_gn<-data.frame(Species=character(),
                      genus=character())
# loop fills ^obj w/ rows whose "genus" is in sp_gn >1x
for (i in 1:nrow(sp_gn)) {
  if (nrow(sp_gn[sp_gn$genus==sp_gn[i,2],])>1) {
    rep_sp_gn<-rep_sp_gn %>% add_row(sp_gn[i,]) 
  }
}
# make df w/ species+genus for species in every combo (only one species per genus)
# empty df to hold app rows
single_sp_gn<-data.frame(Species=character(),
                         genus=character())
# loop fills ^ obj w/ rows whose "genus" is in sp_gn 1x
for (i in 1:nrow(sp_gn)) {
  if (nrow(sp_gn[sp_gn$genus==sp_gn[i,2],])==1) {
    single_sp_gn<-single_sp_gn %>% add_row(sp_gn[i,]) 
  }
}

tbl=table(rep_sp_gn$genus)
all_sp_poss = matrix(NA, ncol=length(unique(rep_sp_gn$genus)))
unq_gn=unique(rep_sp_gn$genus)

for (i in 1:tbl[1]){
  for (j in 1:tbl[2]){
    for (k in 1:tbl[3]){
      for (m in 1:tbl[4]){
        for (n in 1:tbl[5]){
          for (o in 1:tbl[6]){
            for (p in 1:tbl[7]) {
              gni=unq_gn[1]
              spi=rep_sp_gn$Species[rep_sp_gn$genus==gni][i]
              gnj=unq_gn[2]
              spj=rep_sp_gn$Species[rep_sp_gn$genus==gnj][j]
              gnk=unq_gn[3]
              spk=rep_sp_gn$Species[rep_sp_gn$genus==gnk][k]
              gnm=unq_gn[4]
              spm=rep_sp_gn$Species[rep_sp_gn$genus==gnm][m]
              gnn=unq_gn[5]
              spn=rep_sp_gn$Species[rep_sp_gn$genus==gnn][n]
              gno=unq_gn[6]
              spo=rep_sp_gn$Species[rep_sp_gn$genus==gno][o]
              gnp=unq_gn[7]
              spp=rep_sp_gn$Species[rep_sp_gn$genus==gnp][p]
              all_sp_poss=rbind(all_sp_poss,c(spi,spj,spk,spm,spn,spo,spp))
            }
          }
        }
      } 
    }
  }
}
#remove NA row at top
all_sp_poss<-all_sp_poss[-1,]
#make empty matrix to be filled w/ all combos
all_nest.depth_subs<-matrix(data=NA,nrow=288,ncol=14)
# combine single rep species (single_sp_gn) w/ all_sp_poss to have matrix w/ each combo (all_aa_subs)
for (i in 1:nrow(all_sp_poss)) {
  all_nest.depth_subs[i,]<-c(all_sp_poss[i,], single_sp_gn[,1])
}

## Use each row in a phylosig call to get lambda and K
nest.depth.res<-data.frame(nest.depth.lambda=numeric(),
                           nest.depth.lambda.p=numeric(),
                           nest.depth.K=numeric(),
                           nest.depth.K.p=numeric())

for (i in 1:nrow(all_nest.depth_subs)) {
  ## save ith subset to df
  nest.depth_subi<-nest.depth_by_species %>% filter(Species %in% all_nest.depth_subs[i,])
  rownames(nest.depth_subi)<-nest.depth_subi$genus 
  ## make object for phylosig call
  nest.depth.trait<-setNames(nest.depth_subi[,2],rownames(nest.depth_subi)) 
  ## test for phylogenetic signal
  # w/ lambda
  nest.depth.lambda<-phylosig(pruned.tree.nest.depth, nest.depth.trait, method = "lambda", test=T)
  # w/ K
  nest.depth.K<-phylosig(pruned.tree.nest.depth, nest.depth.trait, method = "K", test=T)
  ## write to output object nest.depth.res
  nest.depth.res[i,"nest.depth.lambda"]<-nest.depth.lambda$lambda
  nest.depth.res[i,"nest.depth.lambda.p"]<-nest.depth.lambda$P
  nest.depth.res[i,"nest.depth.K"]<-nest.depth.K$K
  nest.depth.res[i,"nest.depth.K.p"]<-nest.depth.K$P
}


#################### Nest Depth / Weber's Length
# "Main" subset for illustrative purposes will contain representative of genus w/ most nests (1st by alph if ties)
scaled.nest.depth_subset_main<-subset(scaled.nest.depth_by_species,
                                      Species!="Acromyrmex landolti" &
                                        Species!="Acromyrmex rugosus" &
                                        Species!="Aphaenogaster ashmeadi" &
                                        Species!="Aphaenogaster treatae" &
                                        Species!="Atta laevigata" &
                                        Species!="Dorymyrmex bureni" &
                                        Species!="Ectatomma opaciventre" &
                                        Species!="Formica archboldi" &
                                        Species!="Mycetarotes acutus" &
                                        Species!="Odontomachus chelifer")
# rownames need to be same as tree tips (genera) 
rownames(scaled.nest.depth_subset_main)<-scaled.nest.depth_subset_main$genus 

## Save species and scaled nest depth cols w/ genus rownames for heatmap
scaled.depth.for.heatmap<-scaled.nest.depth_subset_main[,-c(1,4,5)]

## Test main subset for phylosig
scaled.nest.depth.trait<-setNames(scaled.nest.depth_subset_main[,2],rownames(scaled.nest.depth_subset_main)) 
## test for phylogenetic signal
# w/ lambda
main_scaled.nest.depth.lambda<-phylosig(pruned.tree.nest.depth, scaled.nest.depth.trait, method = "lambda", test=T)
# w/ K
main_scaled.nest.depth.K<-phylosig(pruned.tree.nest.depth, scaled.nest.depth.trait, method = "K", test=T)

## Save main subset results to be combined w/ other features later
main_scaled.nest.depth.res<-data.frame(feature.name=character(),
                                       lambda=numeric(),
                                       lambda.p=numeric(),
                                       K=numeric(),
                                       K.p=numeric())
# write res to cols
main_scaled.nest.depth.res[1,"feature.name"]<-"Nest Depth / WL"
main_scaled.nest.depth.res[1,"lambda"]<-main_scaled.nest.depth.lambda$lambda
main_scaled.nest.depth.res[1,"lambda.p"]<-main_scaled.nest.depth.lambda$P
main_scaled.nest.depth.res[1,"K"]<-main_scaled.nest.depth.K$K
main_scaled.nest.depth.res[1,"K.p"]<-main_scaled.nest.depth.K$P

## Use each row in a phylosig call to get lambda and K
# make df to fill w/ res from each combination of representative species
scaled.nest.depth.res<-data.frame(scaled.nest.depth.lambda=numeric(),
                                  scaled.nest.depth.lambda.p=numeric(),
                                  scaled.nest.depth.K=numeric(),
                                  scaled.nest.depth.K.p=numeric())
# write res to df
for (i in 1:nrow(all_nest.depth_subs)) {
  ## save ith subset to df
  scaled.nest.depth_subi<-scaled.nest.depth_by_species %>% filter(Species %in% all_nest.depth_subs[i,])
  rownames(scaled.nest.depth_subi)<-scaled.nest.depth_subi$genus 
  ## make object for phylosig call
  scaled.nest.depth.trait<-setNames(scaled.nest.depth_subi[,2],rownames(scaled.nest.depth_subi)) 
  ## test for phylogenetic signal
  # w/ lambda
  scaled.nest.depth.lambda<-phylosig(pruned.tree.nest.depth, scaled.nest.depth.trait, method = "lambda", test=T)
  # w/ K
  scaled.nest.depth.K<-phylosig(pruned.tree.nest.depth, scaled.nest.depth.trait, method = "K", test=T)
  ## write to output object nest.depth.res
  scaled.nest.depth.res[i,"scaled.nest.depth.lambda"]<-scaled.nest.depth.lambda$lambda
  scaled.nest.depth.res[i,"scaled.nest.depth.lambda.p"]<-scaled.nest.depth.lambda$P
  scaled.nest.depth.res[i,"scaled.nest.depth.K"]<-scaled.nest.depth.K$K
  scaled.nest.depth.res[i,"scaled.nest.depth.K.p"]<-scaled.nest.depth.K$P
}

######################### HEATMAP
### Generate genus-level summary of measures for heatmap
## Combine .for.heatmap objects
# give each heatmap a genus col out of its rownames
width.for.heatmap<-tibble::rownames_to_column(width.for.heatmap, "genus")
density.for.heatmap<-tibble::rownames_to_column(density.for.heatmap, "genus")
num.cham.for.heatmap<-tibble::rownames_to_column(num.cham.for.heatmap, "genus")
depth.for.heatmap<-tibble::rownames_to_column(depth.for.heatmap, "genus")
scaled.depth.for.heatmap<-tibble::rownames_to_column(scaled.depth.for.heatmap, "genus")

# list .for.heatmaps to be joined
heatmap.dfs<-list(width.for.heatmap,
                  density.for.heatmap,
                  num.cham.for.heatmap,
                  depth.for.heatmap,
                  scaled.depth.for.heatmap)
# join them
all.for.heatmap<-heatmap.dfs %>% reduce(full_join, by="genus")

# set rownames (to be paired w/ tree tips) and col names
rownames(all.for.heatmap)<-all.for.heatmap$genus 
names(all.for.heatmap)[1]<-c("Genus",
                             "ECW / WL","ECW / HW",
                             "ECW / WL CV","ECW / HW CV",
                             "Density","Density CV",
                             "Chamber #","Chamber # CV",
                             "Nest Depth","Nest Depth CV",
                             "Nest Depth / WL","Nest Depth / WL CV")
# separate CVs and measures for two maps
features.for.heatmap<-all.for.heatmap[,c(2,3,6,8,10,12)]
cvs.for.heatmap<-all.for.heatmap[,c(4,5,7,9,11,13)]

#Prune tree to the genera in heatmap - all sp w/ >2 nests measured for at least 1 feature
tips<-c("Acromyrmex_versicolor",
        "Aphaenogaster_occidentalis_NW",
        "Atta_texana",
        "Camponotus_maritimus",
        "Diacamma_rugosum",
        "Dinoponera_australis",
        "Dorymyrmex_bicolor",
        "Ectatomma_opaciventre",
        "Formica_moki",
        "Mycetagroicus_triangularis",
        "Mycetarotes_acutus",
        "Mycetophylax_conformis",
        "Mycocepurus_goeldii",
        "Odontomachus_coquereli",
        "Pachycondyla_harpax",
        "Pheidole_longispinosa",
        "Prenolepis_imparis",
        "Sericomyrmex_Sp")
pruned.tree.all<-drop.tip(ant.tree.moreau,ant.tree.moreau$tip.label[-match(tips, ant.tree.moreau$tip.label)])

#Rename tip labels to create genus level tree
for (i in 1:length(pruned.tree.all$tip.label)) {
  split.tips<-strsplit(pruned.tree.all$tip.label[i],"_")
  genus.tip<-split.tips[[1]]
  pruned.tree.all$tip.label[i]<-genus.tip[1]
}

########################### PHYLOGENETIC SIGNAL SUMMARY
### Report tests of phylogenetic signal on main subset for each feature
## Put them together in a table
phylosig_all.features<-rbind(main_ww.res,main_hw.res,
                             main_density.res,
                             main_num.cham.res,
                             main_nest.depth.res,main_scaled.nest.depth.res)
# make number cols numeric to allow rounding
phylosig_all.features[,c(2:5)]<-phylosig_all.features[,as.numeric(c(2:5))]
# round off extra digits
phylosig_all.features$lambda <- round(phylosig_all.features$lambda, digits = 3)
phylosig_all.features$lambda.p <- round(phylosig_all.features$lambda.p, digits = 3)
phylosig_all.features$K <- round(phylosig_all.features$K, digits = 3)
# replace really small numbers with "<0.0001"
for (i in 1:nrow(phylosig_all.features)) {
  if (as.numeric(phylosig_all.features$'lambda'[i]) < 0.0001) {
    phylosig_all.features$'lambda'[i] <- "<0.0001"
  }
}
# set column names to something presentable
colnames(phylosig_all.features)[colnames(phylosig_all.features) == 'feature.name'] <- "Feature Name"
colnames(phylosig_all.features)[colnames(phylosig_all.features) == 'lambda'] <- "Pagel's Lambda"
colnames(phylosig_all.features)[colnames(phylosig_all.features) == 'lambda.p'] <- "p-value"
colnames(phylosig_all.features)[colnames(phylosig_all.features) == 'K'] <- "Blomberg's K"
colnames(phylosig_all.features)[colnames(phylosig_all.features) == 'K.p'] <- "p-value"
# fix order of features to align w/ figs and text
phylosig_all.features<-phylosig_all.features[c(5,6,4,3,2,1),]

# split into main text table and supp mat table (lambda for main, K for supp mat)
phylosig_all.features_main<-phylosig_all.features[,c(1,2,3)]
phylosig_all.features_supp<-phylosig_all.features[,c(1,4,5)]

## prepare main for plotting
# turn into gtable for customizble plotting
phylosig_table_main<-tableGrob(phylosig_all.features_main,
                               theme=ttheme_minimal(),
                               rows=NULL)
# make lines pretty
phylosig_table_main <- gtable_add_grob(phylosig_table_main,
                                       grobs = rectGrob(gp = gpar(fill = NA, lwd = 1)),
                                       t = 1, b = nrow(phylosig_table_main), l = 1, r = ncol(phylosig_table_main))
phylosig_table_main <- gtable_add_grob(phylosig_table_main,
                                       grobs = rectGrob(gp = gpar(fill = NA, lwd = 1)),
                                       t = 2, b = nrow(phylosig_table_main), l = 1, r = ncol(phylosig_table_main))

## prepare supp mat table for plotting
# turn into gtable for customizble plotting
phylosig_table_supp<-tableGrob(phylosig_all.features_supp,
                               theme=ttheme_minimal(),
                               rows=NULL)
# make lines pretty
phylosig_table_supp <- gtable_add_grob(phylosig_table_supp,
                                       grobs = rectGrob(gp = gpar(fill = NA, lwd = 1)),
                                       t = 1, b = nrow(phylosig_table_supp), l = 1, r = ncol(phylosig_table_supp))
phylosig_table_supp <- gtable_add_grob(phylosig_table_supp,
                                       grobs = rectGrob(gp = gpar(fill = NA, lwd = 1)),
                                       t = 2, b = nrow(phylosig_table_supp), l = 1, r = ncol(phylosig_table_supp))

### Summarize phylogenetic signal for ALL subsets
## Ready for combining
# split widths
ww.widths.res<-widths.res[,c(1:4)]
hw.widths.res<-widths.res[,c(5:8)]
# add feature name
ww.widths.res$feature<-"ww.width"
hw.widths.res$feature<-"hw.width"
density.res$feature<-"density"
num.cham.res$feature<-"num.cham"
nest.depth.res$feature<-"nest.depth"
scaled.nest.depth.res$feature<-"scaled.nest.depth"
## make all objects have same col names for rbind
names(ww.widths.res)<-c("lambda","lambda p-value","K","K p-value","feature")
names(hw.widths.res)<-c("lambda","lambda p-value","K","K p-value","feature")
names(density.res)<-c("lambda","lambda p-value","K","K p-value","feature")
names(num.cham.res)<-c("lambda","lambda p-value","K","K p-value","feature")
names(nest.depth.res)<-c("lambda","lambda p-value","K","K p-value","feature")
names(scaled.nest.depth.res)<-c("lambda","lambda p-value","K","K p-value","feature")

## rbind all subset res
all.subset.res<-rbind(ww.widths.res,
                      hw.widths.res,
                      density.res,
                      num.cham.res,
                      nest.depth.res,
                      scaled.nest.depth.res)

## arrange features as factors to display in chosen order
all.subset.res$feature<-as.factor(all.subset.res$feature)
all.subset.res$feature<-factor(all.subset.res$feature, levels = c("nest.depth","scaled.nest.depth",
                                                                  "num.cham","density",
                                                                  "hw.width","ww.width"))


################################################ FORAGING V. NEST FEATURES ###########################################
## First, add strat to objects w/ each nest measured for each feature
# Make new objects for foraging analysis
ec_foraging<-ec_w_reps
density_foraging<-nets_w_reps
num.cham_foraging<-num.cham_w_reps
depths_foraging<-nest_depths

### Add foraging strat column to each object
## E. chamber widths
# associate foraging strategy with each nest
for (i in 1:nrow(ec_foraging)){
  species_name <- ec_foraging[i,"Species"]
  if (sum(foraging_strats$Species.Name==species_name) > 0) {
    ec_foraging$foraging.strat[i] <- foraging_strats[foraging_strats$Species.Name==species_name,
                                                     "Foraging.Strategy"]
  }
  else {
    ec_foraging$foraging.strat[i] <- NA
  }
}
# remove species w/o foraging strat assigned
ec_foraging <- ec_foraging[!(ec_foraging$foraging.strat==""),]
ec_foraging<-na.omit(ec_foraging)

## Density
# associate foraging strategy with each nest
for (i in 1:nrow(density_foraging)){
  species_name <- density_foraging[i,"Species"]
  if (sum(foraging_strats$Species==species_name) > 0) {
    density_foraging[i,"foraging.strat"] <- foraging_strats[foraging_strats$Species==species_name,
                                                            "Foraging.Strategy"]
  }
  else {
    density_foraging$foraging.strat[i] <- NA
  }
}
# remove species w/o foraging strat assigned
density_foraging<-density_foraging[!(density_foraging$foraging.strat==""),]
density_foraging<-na.omit(density_foraging)

## Num.Cham
# associate foraging strategy with each nest
for (i in 1:nrow(num.cham_foraging)){
  species_name <- num.cham_foraging[i,"Species"]
  if (sum(foraging_strats$Species==species_name) > 0) {
    num.cham_foraging[i,"foraging.strat"] <- foraging_strats[foraging_strats$Species==species_name,
                                                             "Foraging.Strategy"]
  }
  else {
    num.cham_foraging$foraging.strat[i] <- NA
  }
}
# remove species w/o foraging strat assigned
num.cham_foraging<-num.cham_foraging[!(num.cham_foraging$foraging.strat==""),]
num.cham_foraging<-na.omit(num.cham_foraging)
# remove Atta laevigata bc it's an outlier
num.cham_foraging<-num.cham_foraging[!(num.cham_foraging$Species=="Atta laevigata"),]


## Depths
# associate foraging strategy with each nest
for (i in 1:nrow(depths_foraging)){
  species_name <- depths_foraging[i,"Species"]
  if (sum(foraging_strats$Species==species_name) > 0) {
    depths_foraging[i,"foraging.strat"] <- foraging_strats[foraging_strats$Species==species_name,
                                                           "Foraging.Strategy"]
  }
  else {
    depths_foraging$foraging.strat[i] <- NA
  }
}
# remove species w/o foraging strat assigned
depths_foraging<-depths_foraging[!(depths_foraging$foraging.strat==""),]
depths_foraging<-na.omit(depths_foraging)


#### Test features vs. foraging strat
### W/ individual nests
# Run ANOVA on ECW/WL w/ foraging strat and species as "effects"
ecw.wl_foraging_aov<-aov(Width_by_Webers ~ foraging.strat + Species, data = ec_foraging)

# Run ANOVA on ECW/WL w/ foraging strat and species as "effects"
ecw.hw_foraging_aov<-aov(Width_by_Head ~ foraging.strat + Species, data = ec_foraging)

# Kruskal-Wallis on density w/ foraging strat as effect
density_foraging_kw <- kruskal.test(Density ~ foraging.strat,
                                    data = density_foraging)

# Kruskal-Wallis on num.cham w/ foraging strat as effect
num.cham_foraging_kw <- kruskal.test(Num.Chambers ~ foraging.strat,
                                     data = num.cham_foraging)

# Kruskal-Wallis on unscaled depth w/ foraging strat as effect
nest.depth_foraging_kw <- kruskal.test(Nest.Depth ~ foraging.strat,
                                       data = depths_foraging)

# Kruskal-Wallis on unscaled depth w/ foraging strat as effect
scaled.depth_foraging_kw <- kruskal.test(scaled.depth ~ foraging.strat,
                                         data = depths_foraging)

### Summarize results
## Report results summary of features v. foraging strats in tables
## Make table for features tested by ANOVA
hw.aov.for_stats<-c(summary(ecw.hw_foraging_aov)[[1]][["Df"]][[1]],
                    summary(ecw.hw_foraging_aov)[[1]][["Sum Sq"]][[1]],
                    summary(ecw.hw_foraging_aov)[[1]][["F value"]][[1]],
                    summary(ecw.hw_foraging_aov)[[1]][["Pr(>F)"]][[1]])

ww.aov.for_stats<-c(summary(ecw.wl_foraging_aov)[[1]][["Df"]][[1]],
                    summary(ecw.wl_foraging_aov)[[1]][["Sum Sq"]][[1]],
                    summary(ecw.wl_foraging_aov)[[1]][["F value"]][[1]],
                    summary(ecw.wl_foraging_aov)[[1]][["Pr(>F)"]][[1]])

aov.for_stats<-rbind(hw.aov.for_stats, ww.aov.for_stats)

aov.for_stats<-as.data.frame(aov.for_stats)

aov.for_stats$feature.name<-c("ECW / Head Width","ECW / Weber's Length")

# give cols names for table
names(aov.for_stats)<-c("Df","Sum Sq","F value","p-value","Feature Name")

# round to 3 digits
aov.for_stats$`Sum Sq` <- round(aov.for_stats$`Sum Sq`, digits = 3)
aov.for_stats$`F value` <- round(aov.for_stats$`F value`, digits = 3)
aov.for_stats$`p-value` <- round(aov.for_stats$`p-value`, digits = 3)

# replace really small numbers with "<0.0001"
for (i in 1:nrow(aov.for_stats)) {
  if (as.numeric(aov.for_stats$`p-value`[i]) < 0.0001) {
    aov.for_stats$`p-value`[i] <- "<0.0001"
  }
}

# fix order of rows
aov.for_stats<-aov.for_stats[,c(5,1,2,3,4)]

aov.for_table<-tableGrob(aov.for_stats,
                         theme=ttheme_minimal(),
                         rows=NULL)

aov.for_table <- gtable_add_grob(aov.for_table,
                                 grobs = rectGrob(gp = gpar(fill = NA, lwd = 1)),
                                 t = 1, b = nrow(aov.for_table), l = 1, r = ncol(aov.for_table))

aov.for_table <- gtable_add_grob(aov.for_table,
                                 grobs = rectGrob(gp = gpar(fill = NA, lwd = 1)),
                                 t = 2, b = nrow(aov.for_table), l = 1, r = ncol(aov.for_table))


## Make table for features v. foraging strats tested w/ K-W

kw.for_stats<-as.data.frame(rbind(c(nest.depth_foraging_kw$parameter,nest.depth_foraging_kw$statistic,nest.depth_foraging_kw$p.value),
                                  c(scaled.depth_foraging_kw$parameter,scaled.depth_foraging_kw$statistic,scaled.depth_foraging_kw$p.value),
                                  c(num.cham_foraging_kw$parameter,num.cham_foraging_kw$statistic,num.cham_foraging_kw$p.value),
                                  c(density_foraging_kw$parameter,density_foraging_kw$statistic,density_foraging_kw$p.value)))

# make col w/ feature names
kw.for_stats$feature.names<-c("Nest Depth",
                              "Scaled Nest Depth",
                              "Number of Chambers",
                              "Network Density")

# name col w/ p-values
names(kw.for_stats)[3]<-"p-value"
names(kw.for_stats)[4]<-"Feature Name"

# round KW chi-squareds to 3 digits
kw.for_stats$`Kruskal-Wallis chi-squared`<-round(kw.for_stats$`Kruskal-Wallis chi-squared`, digits = 3)

# round p-values to 3 digits, replace tiny numbers w/ "<0.0001"
for (i in 1:nrow(kw.for_stats)) {
  if (as.numeric(kw.for_stats$'p-value'[i]) < 0.0001) {
    kw.for_stats$'p-value'[i] <- "<0.0001"
  }
  else{
    kw.for_stats$'p-value'[i] <- round(as.numeric(kw.for_stats$'p-value'[i]), digits = 3)
  }
}

#fix order of rows
kw.for_stats<-kw.for_stats[,c(4,1,2,3)]

kw.for_table<-tableGrob(kw.for_stats,
                        theme=ttheme_minimal(),
                        rows=NULL)

kw.for_table <- gtable_add_grob(kw.for_table,
                                grobs = rectGrob(gp = gpar(fill = NA, lwd = 1)),
                                t = 1, b = nrow(kw.for_table), l = 1, r = ncol(kw.for_table))

kw.for_table <- gtable_add_grob(kw.for_table,
                                grobs = rectGrob(gp = gpar(fill = NA, lwd = 1)),
                                t = 2, b = nrow(kw.for_table), l = 1, r = ncol(kw.for_table))


### Post hoc testing differences among foraging groups
# Dunn test for nest depth / WL
dunnTest(scaled.depth ~ foraging.strat,
         data = depths_foraging)
# Dunn test for unscaled nest depth
dunnTest(Nest.Depth ~ foraging.strat,
         data = depths_foraging)
# Dunn test for density
dunnTest(Density ~ foraging.strat,
         data = density_foraging)
# Dunn test for number of chambers
dunnTest(Num.Chambers ~ foraging.strat,
         data = num.cham_foraging)
# Tukey test for ECW / HW
TukeyHSD(ecw.hw_foraging_aov)
# Tukey test for ECW / WL
TukeyHSD(ecw.wl_foraging_aov)

### Prepare for visualization
# modify objects to order foraging strats
ec_foraging$foraging.strat<-factor(ec_foraging$foraging.strat, levels = c("Solitary Foraging",
                                                                          "Group Recruitment","Stable Trail","Mass Recruitment"))
density_foraging$foraging.strat<-factor(density_foraging$foraging.strat, levels = c("Solitary Foraging",
                                                                                    "Group Recruitment","Stable Trail","Mass Recruitment"))
num.cham_foraging$foraging.strat<-factor(num.cham_foraging$foraging.strat, levels = c("Solitary Foraging",
                                                                                      "Group Recruitment","Stable Trail","Mass Recruitment"))
depths_foraging$foraging.strat<-factor(depths_foraging$foraging.strat, levels = c("Solitary Foraging",
                                                                                  "Group Recruitment","Stable Trail","Mass Recruitment"))



############################################### OUTPUTS #################################################
density_v_num.cham<-plot(Density ~ Num.Chambers, data=nets_w_reps)


dev.off()
par(mar=c(5.1, 4.1, 0, 2.1))
features_heatmap<-phylo.heatmap(mar = c(0.1,0.1,0,0.1),
                                pruned.tree.all,features.for.heatmap,
                                standardize=TRUE,
                                lwd=3,
                                pts=FALSE,
                                legend=TRUE,
                                split=c(0.6,0.4),
                                ylim = c(-0.18,1.23))


dev.off()
par(mar=c(5.1, 4.1, 0, 2.1))
cvs_heatmap<-phylo.heatmap(mar = c(0.1,0.1,0,0.1),
                           pruned.tree.all,
                           cvs.for.heatmap,
                           standardize=TRUE,
                           lwd=3,
                           pts=FALSE,
                           legend=TRUE,
                           split=c(0.6,0.4),
                           ylim = c(-0.18,1.23))

# make presentable
dev.off()
grid.draw(species.summary_table.print)


### Report tests of differences among species
dev.off()
grid.draw(aov.species_table)

# print table to Plots
dev.off()
grid.draw(kw.species_table)

# display main table in Plots
dev.off()
grid.draw(phylosig_table_main)

# display main table in Plots
dev.off()
grid.draw(phylosig_table_supp)


### Pagel's Lambda for main text
## plot lambda values
dev.off()
par(mfrow = c(2,1),
    mar = c(4.5,5.5,2,5))
boxplot(lambda ~ feature, all.subset.res,
        ylab = "Lambda",
        names = c("Nest Depth","Nest Depth / WL","Chamber #","Density","ECW / HW","ECW / WL"),
        xlab = "Nest Feature",
        col = "white",
        ylim = c(0,1.25),
        cex.axis = 1,
        cex.lab = 1.5,
        las = 1)
# add n - number of subsets
text(x= 1, y= max(all.subset.res$lambda[all.subset.res$feature == "nest.depth"])+0.05, labels= "n = 288", adj = 0.4, cex = 1)
text(x= 2, y= max(all.subset.res$lambda[all.subset.res$feature == "scaled.nest.depth"])+0.05, labels= "n = 288", adj = 0.4, cex = 1)
text(x= 3, y= max(all.subset.res$lambda[all.subset.res$feature == "num.cham"])+0.05, labels= "n = 576", adj = 0.4, cex = 1)
text(x= 4, y= max(all.subset.res$lambda[all.subset.res$feature == "density"])+0.05, labels= "n = 48", adj = 0.4, cex = 1)
text(x= 5, y= max(all.subset.res$lambda[all.subset.res$feature == "hw.width"])+0.05, labels= "n = 48", adj = 0.4, cex = 1)
text(x= 6, y= max(all.subset.res$lambda[all.subset.res$feature == "ww.width"])+0.05, labels= "n = 48", adj = 0.4, cex = 1)
# add * for representative "main" subset value
text(x= 6, y= 0.497, labels= "*", col = "blue", adj = 0.5, cex = 2)
text(x= 5, y= 0, labels= "*", col = "blue", adj = 0.5, cex = 2)
text(x= 4, y= 0.141, labels= "*", col = "blue", adj = 0.5, cex = 2)
text(x= 3, y= 0.115, labels= "*", col = "blue", adj = 0.5, cex = 2)
text(x= 2, y= 0, labels= "*", col = "blue", adj = 0.5, cex = 2)
text(x= 1, y= 0, labels= "*", col = "blue", adj = 0.5, cex = 2)

## plot lambda p-values
boxplot(all.subset.res$`lambda p-value` ~ feature, all.subset.res,
        ylab = "Likelihood ratio test p-value",
        names = c("Nest Depth","Nest Depth / WL","Chamber #","Density","ECW / HW","ECW / WL"),
        xlab = "Nest Feature",
        col = "white",
        ylim = c(0,1),
        cex.axis = 1,
        cex.lab = 1.5,
        las = 1)
# add * for representative "main" subset value
text(x= 6, y= 0.199, labels= "*", col = "blue", adj = 0.5, cex = 2)
text(x= 5, y= 1, labels= "*", col = "blue", adj = 0.5, cex = 2)
text(x= 4, y= 0.851, labels= "*", col = "blue", adj = 0.5, cex = 2)
text(x= 3, y= 0.060, labels= "*", col = "blue", adj = 0.5, cex = 2)
text(x= 2, y= 1, labels= "*", col = "blue", adj = 0.5, cex = 2)
text(x= 1, y= 1, labels= "*", col = "blue", adj = 0.5, cex = 2)

### Blomberg's K for supp mat
dev.off()
par(mfrow = c(2,1),
    mar = c(4.5,5.5,2,5))
boxplot(K ~ feature, all.subset.res,
        ylab = "Blomberg's K",
        names = c("Nest Depth","Nest Depth / WL","Chamber #","Density","ECW / HW","ECW / WL"),
        xlab = "Nest Feature",
        col = "white",
        ylim = c(0,1.25),
        cex.axis = 1,
        cex.lab = 1.5,
        las = 1)
# add n - number of subsets
text(x= 1, y= max(all.subset.res$lambda[all.subset.res$feature == "nest.depth"])+0.05, labels= "n = 288", adj = 0.4, cex = 1)
text(x= 2, y= max(all.subset.res$lambda[all.subset.res$feature == "scaled.nest.depth"])+0.05, labels= "n = 288", adj = 0.4, cex = 1)
text(x= 3, y= max(all.subset.res$lambda[all.subset.res$feature == "density"])+0.05, labels= "n = 48", adj = 0.4, cex = 1)
text(x= 4, y= max(all.subset.res$lambda[all.subset.res$feature == "num.cham"])+0.05, labels= "n = 576", adj = 0.4, cex = 1)
text(x= 5, y= max(all.subset.res$lambda[all.subset.res$feature == "hw.width"])+0.05, labels= "n = 48", adj = 0.4, cex = 1)
text(x= 6, y= max(all.subset.res$lambda[all.subset.res$feature == "ww.width"])+0.05, labels= "n = 48", adj = 0.4, cex = 1)
# add * for representative "main" subset value
text(x= 6, y= 0.797, labels= "*", col = "blue", adj = 0.5, cex = 2)
text(x= 5, y= 0.469, labels= "*", col = "blue", adj = 0.5, cex = 2)
text(x= 4, y= 0.895, labels= "*", col = "blue", adj = 0.5, cex = 2)
text(x= 3, y= 0.646, labels= "*", col = "blue", adj = 0.5, cex = 2)
text(x= 2, y= 0.555, labels= "*", col = "blue", adj = 0.5, cex = 2)
text(x= 1, y= 0.492, labels= "*", col = "blue", adj = 0.5, cex = 2)

## plot lambda p-values
boxplot(all.subset.res$`K p-value` ~ feature, all.subset.res,
        ylab = "p-value",
        names = c("Nest Depth","Nest Depth / WL","Chamber #","Density","ECW / HW","ECW / WL"),
        xlab = "Nest Feature",
        col = "white",
        ylim = c(0,1),
        cex.axis = 1,
        cex.lab = 1.5,
        las = 1)
# add * for representative "main" subset value
text(x= 6, y= 0.089, labels= "*", col = "blue", adj = 0.5, cex = 2)
text(x= 5, y= 0.674, labels= "*", col = "blue", adj = 0.5, cex = 2)
text(x= 4, y= 0.067, labels= "*", col = "blue", adj = 0.5, cex = 2)
text(x= 3, y= 0.302, labels= "*", col = "blue", adj = 0.5, cex = 2)
text(x= 2, y= 0.485, labels= "*", col = "blue", adj = 0.5, cex = 2)
text(x= 1, y= 0.537, labels= "*", col = "blue", adj = 0.5, cex = 2)








dev.off()
grid.draw(aov.for_table)

# print table to Plots
dev.off()
grid.draw(kw.for_table)


## Set multi-panel output
dev.off()
par(mfrow=c(3,2),
    mar=c(5,6,4,2)+0.1)

## Unscaled nest depth
nest.depth_boxplot<-boxplot(Nest.Depth ~ foraging.strat, data = depths_foraging,
                            ylab = "Nest Depth", xlab = "",
                            #names = c("Solitary\nForaging","Group\nRecruitment","Stable\nTrail","Mass\nRecruitment"),
                            col = c("#9965b9","#6cb966","#b86565","#6598b9"),
                            border = c("#5c3278","#327732","#773332","#355474"),
                            outcol = c("#5c3278","#327732","#773332","#355474"),
                            whiskcol = c("#9965b9","#6cb966","#b86565","#6598b9"),
                            whisklty = 1,
                            ylim = c(0,4),
                            cex.axis = 1.5,
                            cex.lab = 2.1,
                            pars = list(boxwex = 0.65, staplewex = 0.38, outwex = 0.7),
                            las=1)
title("A. Nest Depth", adj = 0, cex.main = 2.5)
text(x= 1, y= 3.9, labels= "A", adj = 0.4, cex = 1.3)
text(x= 2, y= 3.9, labels= "A", adj = 0.4, cex = 1.3)
text(x= 3, y= 3.9, labels= "B", adj = 0.4, cex = 1.3)
text(x= 4, y= 3.9, labels= "B", adj = 0.4, cex = 1.3)

## Nest depth scaled by Weber's Length
scaled.depth_boxplot<-boxplot(scaled.depth ~ foraging.strat, data = depths_foraging,
                              ylab = "Nest Depth / Weber's Length", xlab = "",
                              col = c("#9965b9","#6cb966","#b86565","#6598b9"),
                              border = c("#5c3278","#327732","#773332","#355474"),
                              outcol = c("#5c3278","#327732","#773332","#355474"),
                              whiskcol = c("#9965b9","#6cb966","#b86565","#6598b9"),
                              whisklty = 1,
                              ylim = c(0,3.5),
                              cex.axis = 1.5,
                              cex.lab = 2.1,
                              pars = list(boxwex = 0.65, staplewex = 0.38, outwex = 0.7),
                              las=1)
title("B. Nest Depth / Weber's Length", adj = 0, cex.main = 2.5)
text(x= 1, y= 3.5, labels= "A", adj = 0.4, cex = 1.3)
text(x= 2, y= 3.5, labels= "A", adj = 0.4, cex = 1.3)
text(x= 3, y= 3.5, labels= "B", adj = 0.4, cex = 1.3)
text(x= 4, y= 3.5, labels= "C", adj = 0.4, cex = 1.3)

## Density
density_boxplot<-boxplot(Density ~ foraging.strat, data = density_foraging,
                         ylab = "Network Density", xlab = "",
                         col = c("#9965b9","#6cb966","#b86565","#6598b9"),
                         border = c("#5c3278","#327732","#773332","#355474"),
                         outcol = c("#5c3278","#327732","#773332","#355474"),
                         whiskcol = c("#9965b9","#6cb966","#b86565","#6598b9"),
                         whisklty = 1,
                         ylim = c(0,1.2),
                         cex.axis = 1.5,
                         cex.lab = 2.1,
                         pars = list(boxwex = 0.65, staplewex = 0.38, outwex = 0.7),
                         las=1)
title("C. Network Density", adj = 0, cex.main = 2.5)
text(x= 1, y= 1.15, labels= "A", adj = 0.4, cex = 1.3)
text(x= 2, y= 1.15, labels= "AB", adj = 0.4, cex = 1.3)
text(x= 3, y= 1.15, labels= "B", adj = 0.4, cex = 1.3)
text(x= 4, y= 1.15, labels= "B", adj = 0.4, cex = 1.3)

## Number of chambers
num.cham_boxplot<-boxplot(Num.Chambers ~ foraging.strat, data = num.cham_foraging,
                          ylab = "Number of Chambers", xlab = "",
                          col = c("#9965b9","#6cb966","#b86565","#6598b9"),
                          border = c("#5c3278","#327732","#773332","#355474"),
                          outcol = c("#5c3278","#327732","#773332","#355474"),
                          whiskcol = c("#9965b9","#6cb966","#b86565","#6598b9"),
                          whisklty = 1,
                          cex.axis = 1.5,
                          cex.lab = 2.1,
                          pars = list(boxwex = 0.65, staplewex = 0.38, outwex = 0.7),
                          las=1)
title("D. Number of Chambers", adj = 0, cex.main = 2.5)

## E. chamber width scaled by head width
ecw.hw_boxplot<-boxplot(Width_by_Head ~ foraging.strat, data = ec_foraging,
                        ylab = "Entrance Chamber Width / Head Width", xlab = "",
                        col = c("#9965b9","#6cb966","#b86565","#6598b9"),
                        border = c("#5c3278","#327732","#773332","#355474"),
                        outcol = c("#5c3278","#327732","#773332","#355474"),
                        whiskcol = c("#9965b9","#6cb966","#b86565","#6598b9"),
                        whisklty = 1,
                        ylim = c(0,18),
                        cex.axis = 1.5,
                        cex.lab = 2.1,
                        pars = list(boxwex = 0.65, staplewex = 0.38, outwex = 0.7),
                        las=1)
title("E. Entrance Chamber Width / Head Width", adj = 0, cex.main = 2.5)
text(x= 1, y= 17.5, labels= "A", adj = 0.4, cex = 1.3)
text(x= 2, y= 17.5, labels= "B", adj = 0.4, cex = 1.3)
text(x= 3, y= 17.5, labels= "A", adj = 0.4, cex = 1.3)
text(x= 4, y= 17.5, labels= "A", adj = 0.4, cex = 1.3)

## E. chamber width scaled by Weber's length
ecw.wl_boxplot<-boxplot(Width_by_Webers ~ foraging.strat, data = ec_foraging,
                        ylab = "Entrance Chamber Width / Weber's Length", xlab = "",
                        col = c("#9965b9","#6cb966","#b86565","#6598b9"),
                        border = c("#5c3278","#327732","#773332","#355474"),
                        outcol = c("#5c3278","#327732","#773332","#355474"),
                        whiskcol = c("#9965b9","#6cb966","#b86565","#6598b9"),
                        whisklty = 1,
                        ylim = c(0,12),
                        cex.axis = 1.5,
                        cex.lab = 2.1,
                        pars = list(boxwex = 0.65, staplewex = 0.38, outwex = 0.7),
                        las=1)
title("F. Entrance Chamber Width / Weber's Length", adj = 0, cex.main = 2.5)
text(x= 1, y= 11.5, labels= "A", adj = 0.4, cex = 1.3)
text(x= 2, y= 11.5, labels= "B", adj = 0.4, cex = 1.3)
text(x= 3, y= 11.5, labels= "A", adj = 0.4, cex = 1.3)
text(x= 4, y= 11.5, labels= "B", adj = 0.4, cex = 1.3)

