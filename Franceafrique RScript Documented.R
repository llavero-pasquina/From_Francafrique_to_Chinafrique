#Marcel Llavero-Pasquina. 2024. Contact details: marcel.llavero@uab.cat / marcelllaveropasquina@gmail.com
#R Script including data analysis and graphs for the paper:
#Cantoni, R., Llavero-Pasquina, M., Apostolopoulou, E., Gerber, J.F., Martínez-Alier, J.,
#Bond, P. (2024). From Françafrique to Chinafrique? Ecologically unequal exchange and environmental 
#conflicts in extractivist practices in Africa. In World Development
#The script is made available under the CC-BY-NC-SA 4.0 license. Please cite the 
#paper above if you use or modify the script

####LIBRARIES####
library(ggplot2)
library("rmatio")
library(dplyr)

options(scipen=999)

####FUNCTIONS####

add.n <- function(factor){
  factor <- as.factor(factor)
  if(length(grep("n = ", levels(factor))) >1){}else{
    for(i in 1:length(levels(factor))){
      levels(factor)[i] <- paste0(levels(factor)[i], " (n = ", table(factor)[i], ")")
    }
  }
  return(factor)
}

chi.table<- function(Xfactor, Yfactor){
  #Xfactor:
  #The independent variable to compare in the Chi.square (ie. cases$Multinationals or cases$State)
  #Yfactor:
  #The dependent variable to compare (ie. cases$Project.Status)
  Xfactor <- as.factor(Xfactor)
  if(is.list(Yfactor)){
    Yfactor$compile <- rep(NA,length.out = nrow(Yfactor)) 
    for(i in 1:nrow(Yfactor)){
      compilation <- NULL
      for(j in 2:(which(colnames(Yfactor) %in% "Multinationals")-1)){
        if(Yfactor[i,j] == 1|Yfactor[i,j] == "V"){
          compilation <- c(compilation,colnames(Yfactor)[j]) 
        }else{}
      }
      Yfactor$compile[i] <- paste(compilation, collapse = "\n")
    }
    Ytrans <- Yfactor[which(Xfactor == levels(Xfactor)[1]),]
    Ytrans <- as.data.frame(table(unlist(strsplit(Ytrans$compile,"\n"))))
    heading <- c("factor",levels(Xfactor)[1])
    for(i in 2:nlevels(Xfactor)){
      Ytrans1 <- Yfactor$compile[which(Xfactor == levels(Xfactor)[i])]
      Ytrans <- cbind.data.frame(Ytrans,as.data.frame(table(unlist(strsplit(Ytrans1,"\n"))))[,-1])
      heading <- c(heading, levels(Xfactor)[i])
    }
    #Add total n for each outcome
    levels(Ytrans$Var1) <- paste0(levels(Ytrans$Var1)," (n = ",table(unlist(strsplit(Yfactor$compile,"\n"))),")")
    colnames(Ytrans) <- heading
    chitable <- Ytrans
    
  }else{
    if(length(grep("\n",Yfactor)) > 0){
      Ytrans <- Yfactor[which(Xfactor == levels(Xfactor)[1])]
      Ytrans <- as.data.frame(table(unlist(strsplit(Ytrans,"\n"))))
      Ytrans$XFactor <- levels(Xfactor)[1] 
      heading <- c("factor",paste(levels(Xfactor)[1]))
      for(i in 2:nlevels(Xfactor)){
        Ytrans1 <- Yfactor[which(Xfactor == levels(Xfactor)[i])]
        Ytrans1 <- as.data.frame(table(unlist(strsplit(Ytrans1,"\n"))))
        Ytrans1$XFactor <- levels(Xfactor)[i]
        Ytrans <- rbind(Ytrans,Ytrans1)
        heading <- c(heading,paste(levels(Xfactor)[i]))
      }
      Xfac <- unique(Ytrans$Var1)
      data <- NULL
      XFactor <- NULL
      for(i in 1:length(Xfac)){
        num<-NULL
        for(j in 1:nlevels(Xfactor)){
          if(length(which((Ytrans$XFactor %in% levels(Xfactor)[j]) & (Ytrans$Var1 %in% levels(Xfac)[i])))>0){
            num1 <- Ytrans$Freq[which((Ytrans$XFactor %in% levels(Xfactor)[j]) & (Ytrans$Var1 %in% levels(Xfac)[i]))]
          } else {num1 <- 0}
          num <- c(num,num1)
          XFactor <- c(XFactor, levels(Xfactor)[j])
        }
        data <- rbind(data,num)
      }
      YFactor <- Xfac
      Ytrans <- cbind.data.frame(YFactor,data)
      #Add total n for each outcome
      for(i in 1:nlevels(Ytrans$YFactor)){
        levels(Ytrans$YFactor)[i] <- paste0(levels(Ytrans$YFactor)[i]," (n = ",sum(Ytrans[i,2:ncol(Ytrans)]),")")
      }
      colnames(Ytrans) <- heading
      chitable <- Ytrans
    }else{Yfactor <- as.factor(Yfactor)
    Yfactor <- add.n(Yfactor)
    chitable <- as.data.frame(levels(Yfactor))
    heading <- "factor"
    #deparse(substitute(Xfactor)) 
    for(i in 1:nlevels(Xfactor)){
      chitable <- cbind.data.frame(chitable,as.integer(table(Yfactor[which(Xfactor %in% levels(Xfactor)[i])])))
      heading <- c(heading, levels(Xfactor)[i])
    }
    colnames(chitable) <- heading}}
  
  chitable$Ratio <- as.numeric(rep(NA,length.out = nrow(chitable)))
  
  for(i in 1:nrow(chitable)){
    chitable$Ratio[i] <- chitable[i,2]/sum(chitable[i,3:ncol(chitable)],na.rm =T)
  }
  
  chisq <- chisq.test(chitable[,c(-1,-ncol(chitable))])
  chitable <- cbind(chitable,chisq$residuals)
  
  residuals <- rep(NA,length.out = nrow(chitable))
  for(i in 1:nlevels(Xfactor)){
    res<-NULL
    for(j in 1:nrow(chitable)){
      res <- c(res,as.numeric(table(Xfactor)[i]) - chitable[j,i+1])
    }
    residuals <- cbind.data.frame(residuals,res)
  }
  residuals <- residuals[,-1]
  
  chitable$pvalue <- rep(NA,length.out = nrow(chitable))
  
  for(i in 1:nrow(chitable)){
    Pretest <- NULL
    for(j in 1:nlevels(Xfactor)){
      Sample <- as.numeric(c(chitable[i,j+1],residuals[i,j]))
      Pretest <- rbind.data.frame(Pretest,Sample)
    }
    Test <- chisq.test(Pretest)
    chitable$pvalue[i] <- Test$p.value
  }
  expected <- chisq$expected
  colnames(expected) <- gsub("\\b$"," \\(expected\\)",colnames(expected))
  chitable <- cbind(chitable,expected)
  relative_abundance <- chitable[,1:nlevels(Xfactor)+1]/chitable[,(ncol(chitable)-nlevels(Xfactor)+1):ncol(chitable)]
  colnames(relative_abundance) <- gsub("\\b$"," \\(relative_abundance\\)",colnames(relative_abundance))
  chitable <- cbind(chitable,relative_abundance)
  chitable <- chitable[order(chitable$Ratio, decreasing = TRUE), ]
  return(chitable)
}

chi.plot <- function(chitable,Xfactors, style = "proportion"){
  Xfactors <- as.factor(Xfactors)
  total <- length(Xfactors) 
  #Xfactors <- Xfactors[-which(Xfactors %in% levels(Xfactors)[nlevels(Xfactors)])]
  chitable$Yfactor <- factor(chitable$factor, levels = c(as.character(chitable$factor)))
  data <- NULL
  Xfactor <- NULL
  yintercept <- NULL
  a <- 0
  #n <- (ncol(chitable)-4)/2
  n <- nlevels(Xfactors)
  if(style %in% c("bars","points"))
  {chitable <- chitable[order(chitable$pvalue, decreasing = TRUE), ]
  for(i in 1:n){
    data <- c(data,chitable[,(ncol(chitable)-n-1+i)])
    Xfactor <- c(Xfactor,rep(colnames(chitable)[i+1],nlevels(chitable$Yfactor)))
    a <- a+table(Xfactors)[i]
  }
  Yfactor <- rep(chitable$Yfactor,n)
  chitable_data <- cbind.data.frame(Yfactor,data)
  chitable_data$Xfactor <- Xfactor
  significant <- NULL
  for(i in 1:nrow(chitable)){
    significant <- c(significant,if(chitable$pvalue[i] > 0.05) {"no"} else {"yes"})
  }
  chitable_data$significant <- significant
  chitable_data$Xfactor <- as.factor(chitable_data$Xfactor)
  for(i in 1:nrow(chitable_data)){
    if(chitable_data$data[i] == 0){chitable_data$data[i] <-NA}
    else{
      chitable_data$data[i] <- log2(chitable_data$data[i])
    }
  }
  if(style == "bars")
  {chiplot<-ggplot(chitable_data, aes(x = Yfactor, y = data, fill=Xfactor)) +
    geom_col(position=position_dodge()) +
    labs(x = NULL, y = NULL)  +
    theme(legend.position = "bottom") +
    scale_fill_manual(name = paste(Xfactors),values=c("#1e818b","#e2c537","#3BB449","#AF3456","#ED2224","#734C20","#111111","#9B989B","#BE85BA"))+
    coord_flip()+
    #gghighlight(chitable_data$significant %in% "yes",
    #            unhighlighted_params =  list(fill = NULL, alpha = 0.5))+
    geom_hline(yintercept = 0, size=0.5)}
  if(style == "points")
  {a <- max(na.omit(abs(chitable_data$data)))
  chiplot<-ggplot(chitable_data, aes(x = Yfactor, y = data, colour=Xfactor)) +
    geom_point(size = 6) +
    labs(x = NULL, y = NULL)  +
    theme(legend.position = "bottom") +
    scale_colour_manual(values=c("#1e818b","#e2c537","#3BB449","#AF3456","#ED2224","#734C20","#111111","#9B989B","#BE85BA"))+
    coord_flip()+
    gghighlight(chitable_data$significant %in% "yes",
                unhighlighted_params =  list(colour = NULL, alpha = 0.3))+
    geom_hline(yintercept = 0, size=0.5)+
    ylim(c(-a,a))
  #geom_text(hjust=0.5, vjust=0.5, colour="#ffffff", fontface = "bold")
  }
  chiplot
  }
  else
  {for(i in 1:n){
    data <- c(data,chitable[,i+1])
    Xfactor <- c(Xfactor,rep(colnames(chitable)[i+1],nlevels(chitable$Yfactor)))
    a <- a+table(Xfactors)[i]
    yintercept <- c(yintercept,a/total)
  }
    yintercept <- 1-yintercept
    Yfactor <- rep(chitable$Yfactor,n)
    chitable_data <- cbind.data.frame(Yfactor,data)
    chitable_data$Xfactor <- Xfactor
    significant <- NULL
    for(i in 1:nrow(chitable)){
      significant <- c(significant,if(chitable$pvalue[i] > 0.05) {"no"} else {"yes"})
    }
    chitable_data$significant <- significant
    
    chiplot<-ggplot(chitable_data, aes(x = Yfactor, y = data, fill=Xfactor)) +
      geom_col(position="fill",) +
      labs(x = NULL, y = NULL)  +
      theme(legend.position = "bottom") +
      scale_fill_manual(name = paste(Xfactors),values=c("#1e818b","#e2c537","#F6EB13","#7A752F","#734C20","#111111","#9B989B","#BE85BA","#3BB449","#ED2224","#AF3456"))+
      coord_flip()+
      gghighlight(chitable_data$significant %in% "yes",
                  unhighlighted_params =  list(fill = NULL, alpha = 0.5))+
      geom_hline(yintercept = yintercept, size=0.5)}
  
  return(chiplot)
}

#### 1. Import Data ####

setwd("C:/Users/Usuario/Desktop/Marcel/EJAtlas/Writings/Franceafrique/DataAnalysis")

#Data from the EJAtlas - data can be requested at 
#https://ejatlas.org/backoffice/cms/en/data-use-policy/
#For information on companies refer to
#Llavero-Pasquina, M. (2024). A Global Analysis of A Global Analysis of 
#Multinational Corporations' role in Environmental Conflicts 
cases<-read.csv("cases.csv", sep = ",")
Regions <- read.csv("Regions.csv")
Summary <- read.csv("Summarywithattributes.csv", sep = ",")
company_conflicts<-read.csv("company_conflicts.csv", sep = ",")

#Data from EORA EEMRIO database
#Free academic licenses available. Download link: https://worldmrio.com/footprints/
#Wiedmann, T. O., Schandl, H., Lenzen, M., Moran, D., Suh, S., West, J., & Kanemoto, K. 
#(2015). The material footprint of nations. 
#Proceedings of the national academy of sciences, 112(20), 6271-6276.
bilateral <- read.mat("bilateraltrade.mat") #2GB file!

#Data from United Nations World Population Prospects - Accessed: https://www.macrotrends.net/countries/CHN/china/population
Population <- read.csv("population.csv")

#Data from Atlas of Economic Complexity - Accessed: https://oec.world/en/profile/country/chn
OEC <- read.csv("OEC.csv", dec = ",")

#### 2. Select African cases ####

#Assigns a world region based on the case's country information
Regions$Region <- as.character(Regions$Region)
cases$Region <- as.numeric(rep(NA,nrow(cases)))
for(i in 1:nrow(cases)){
  cases$Region[i] <- if(cases$Country[i] %in% "") {"ND"} else {
    Regions$Region[which(Regions$Country %in% cases$Country[i])]
  }
}

#Selects conflicts ocurring in Africa
cases$Region <- as.factor(cases$Region)
levels(cases$Region) <- c("Africa","China","East Asia","Europe","Latin America and Caribbean","Middle East","North America","Pacific","Russia and Central Asia","South Asia","South East Asia")
Africa <- cases[which(cases$Region %in% "Africa"),]

#### 3. Assign subregion & former French colony categories ####

#Assigns the African Union subregion categories based on the case's country information
for(i in 1:nrow(Africa)){
  Africa$Subregion[i] <- Regions$Africa[which(Regions$Country %in% Africa$Country[i])]
  Africa$Colony[i] <- Regions$Colony[which(Regions$Country %in% Africa$Country[i])]
}

#### 4. Identify cases with French and Chinese Companies ####

#For each of the conflicts in Africa, it searches whether the involved companies are documented as originally French or Chinese
Africa$France <- as.character(rep(NA,nrow(Africa)))
Africa$China <- as.character(rep(NA,nrow(Africa)))
for(i in 1:nrow(Africa)){
  if(Africa$Conflict.Id[i] %in% company_conflicts$conflict_id){
    Comp <- as.numeric(company_conflicts[which(company_conflicts$conflict_id %in% Africa$Conflict.Id[i]),seq(from = 2, length.out = 44, by = 4)], na.rm = TRUE)
    Comp <- as.character(Comp[!is.na(Comp)])
    Ctry <- NULL
    for(j in 1:length(Comp)){
      Ctry <- c(Ctry,Summary$Country[grep(paste("\\<",Comp[j],"\\>", sep = ""),Summary$companyIDs)])
    }
  if("France" %in% Ctry)
  {Africa$France[i] <- "yes"} else {Africa$France[i] <- "no"}
  if("China" %in% Ctry)
  {Africa$China[i] <- "yes"} else {Africa$China[i] <- "no"}
}}
rm(i,j,Comp,Ctry)

#Selects conflicts with French and Chinese companies respectively
France <- Africa[which(Africa$France %in% "yes"),]
China <- Africa[which(Africa$China %in% "yes"),]

France$Case[which(France$Conflict.Id %in% China$Conflict.Id)]

#Creates a category to define whether there are companies from France, China or both involved in a given conflict
Africa$frique <- as.character(rep("ND",nrow(Africa)))
for(i in 1:nrow(Africa)){
    if(Africa$Conflict.Id[i] %in% France$Conflict.Id)
    {Africa$frique[i] <- "French"} else {}
    if(Africa$Conflict.Id[i] %in% China$Conflict.Id)
    {Africa$frique[i] <- "Chinese"} else {}
    if(Africa$Conflict.Id[i] %in% France$Conflict.Id & Africa$Conflict.Id[i] %in% China$Conflict.Id)
    {Africa$frique[i] <- "FrenchChinese"} else {}
}

#Selects only the cases with French or Chinese companies involved
AfricaFC <- Africa[which(Africa$frique %in% c("French","Chinese","FrenchChinese")),]

#Duplicates entries containing both French and Chinese companies for graph-making and statistics
FC <- AfricaFC[which(AfricaFC$frique %in% "FrenchChinese"),]
AfricaFC$frique <- gsub("FrenchChinese","French",AfricaFC$frique)
AfricaFC <- rbind.data.frame(AfricaFC,FC)
AfricaFC$frique <- gsub("FrenchChinese","Chinese",AfricaFC$frique)
table(AfricaFC$frique)
rm(FC)

#### 5. Make graphs ####

#Rewords the category names
AfricaFC$First.level.category <- as.factor(AfricaFC$First.level.category)
levels(AfricaFC$First.level.category) <- gsub(" \\(Forests, Agriculture, Fisheries and Livestock Management\\)","",levels(AfricaFC$First.level.category))
levels(Africa$First.level.category) <- gsub(" \\(Forests, Agriculture, Fisheries and Livestock Management\\)","",levels(AfricaFC$First.level.category))

#Plots the distribution of cases by conflict category for conflicts with French or Chinese companies
Category <- ggplot(AfricaFC, aes(x= AfricaFC$frique, fill = AfricaFC$First.level.category))+
  geom_bar(position = "fill",)+
  xlab("Origin of companies")+
  ylab(NULL)+
  scale_fill_manual(name = "Category", values=c("#3BB449","#734C20","#111111","#ED2224","#9B989B","#F47C20","#F6EB13","#7A752F","#51C2EF"))
Category

#Plots the distribution of cases by conflict category for all African conflicts
CategoryAfrica <- ggplot(Africa, aes(x= Africa$Region, fill = Africa$First.level.category))+
  geom_bar(position = "fill",)+
  xlab("Origin of companies")+
  ylab(NULL)+
  scale_fill_manual(name = "Category", values=c("#3BB449","#734C20","#111111","#ED2224","#9B989B","#F47C20","#F6EB13","#BE85BA","#7A752F","#51C2EF"))
CategoryAfrica

#Plots the distribution of cases by African subregions for conflicts with French or Chinese companies
Subregion <- ggplot(AfricaFC, aes(x= AfricaFC$frique, fill = AfricaFC$Subregion))+
  geom_bar(position = "fill",)+
  xlab("Origin of companies")+
  ylab(NULL)+
  scale_fill_manual(name = "Subregion", values=c("#3BB449","#734C20","#111111","#ED2224","#9B989B","#F47C20","#F6EB13","#7A752F","#51C2EF"))
Subregion

#Plots the distribution of cases by African subregions for all African conflicts
SubregionAfrica <- ggplot(Africa, aes(x= Africa$Region, fill = Africa$Subregion))+
  geom_bar(position = "fill",)+
  xlab("Origin of companies")+
  ylab(NULL)+
  scale_fill_manual(name = "Subregion", values=c("#3BB449","#734C20","#111111","#ED2224","#9B989B","#F47C20","#F6EB13","#7A752F","#51C2EF"))
SubregionAfrica

#Plots the distribution of cases by French former colony for conflicts with French or Chinese companies
Colony <- ggplot(AfricaFC, aes(x= AfricaFC$frique, fill = AfricaFC$Colony))+
  geom_bar(position = "fill",)+
  xlab("Origin of companies")+
  ylab(NULL)+
  scale_fill_manual(name = "Subregion", values=c("#3BB449","#734C20","#111111","#ED2224","#9B989B","#F47C20","#F6EB13","#7A752F","#51C2EF"))
Colony

#Plots the distribution of cases by French former colony for all African conflicts
ColonyAfrica <- ggplot(Africa, aes(x= Africa$Region, fill = Africa$Colony))+
  geom_bar(position = "fill",)+
  xlab("Origin of companies")+
  ylab(NULL)+
  scale_fill_manual(name = "Subregion", values=c("#3BB449","#734C20","#111111","#ED2224","#9B989B","#F47C20","#F6EB13","#7A752F","#51C2EF"))
ColonyAfrica

#### 7. Chi-squared ####

#Calculates the chi-squared expected values, residuals and p-value for the distribution
#of category, subregion and former French colony between conflicts with Chinese and French company involvement
CatChi <- chi.table(AfricaFC$frique,AfricaFC$First.level.category)
RegChi <- chi.table(AfricaFC$frique,AfricaFC$Subregion)
ColChi  <- chi.table(AfricaFC$frique,AfricaFC$Colony)

#### 8. Calculate EUE ####

#Selects the raw material indicators of interest
indicators<-bilateral$metadata$indicators
indicators1<-unlist(indicators)
indicators1 <- indicators1[seq(from=2, to=252,by=2)]
selind <- c(1,grep("Raw material",indicators1))

#Selects the countries of interest
a3s <- bilateral$metadata$a3s
a3s <- unlist(a3s)
countrynames <- unlist(bilateral$metadata$countrynames)
years <- unlist(bilateral$metadata$years)
cbind(a3s,countrynames)
Africa <- grep("DZA|AGO|BEN|BWA|GFA|BDI|CMR|CPV|CAF|TCD|COG|CIV|COD|DJI|EGY|ERI|ETH|GAB|
                 GMB|GHA|GIN|KEN|LSO|LBR|LBY|MDG|MWI|MLI|MRT|MAR|MOZ|NAM|NER|NGA|RWA|STP|SEN|SLE|SOM
                 |ZAF|SDS|SUD|SWZ|TGO|UGA|TZA|ZMB|ZWE",a3s)
Colonial <- grep("FRA|CHN",a3s)
bil <-bilateral$bilateraltrade
mydata <- bil[,selind,,]

#Repackages the array into a data frame to facilitate subsequent plotting
values <- as.vector(matrix(aperm(mydata,c(1,2,3,4))))
Year <- rep(years, each = 1, length.out = length(values))
Indicator <- rep(indicators1[selind], each = length(years), length.out=length(values))
Exp <- rep(countrynames, each = length(years)*length(selind), length.out = length(values))
Imp <- rep(countrynames, each = length(years)*length(selind)*length(countrynames), length.out = length(values))
df <- cbind.data.frame(Year,Indicator,Exp,Imp,values)

rm(Exp,Imp,Indicator,mydata,values,Year)

#### 9. Graphs EUE - Material trade ####

#Adds up raw material export data from Africa to France and China since 1989 
df2 <- df[which(df$Year > 1989),]
df2 <- df2[which(df2$Exp %in% c(countrynames[Africa])),]
df2 <- df2[which(df2$Imp %in% countrynames[Colonial]),]
df2 <- df2 %>% group_by(Year,Indicator,Imp) %>% summarise(total = sum(values))
df2$Indicator <- as.factor(df2$Indicator)
levels(df2$Indicator)

#Plots African exports to China and France by type of raw material
df3 <- df2[which(df2$Indicator %in% levels(df2$Indicator)[3:6]),]
Export_bars_item <- ggplot(df3, aes(x=df3$Year, y=df3$total, fill = df3$Indicator))+
  geom_bar(stat = "identity")+
  scale_fill_manual(name = "Indicator", values = c("#734C20","#9B989B","#111111","#F47C20"))+
  facet_wrap(df3$Imp~., scales = "free_y", ncol = 1)+
  ylab(NULL)+
  xlab("Year")+
  ylim(0,400000000)
Export_bars_item

#Adds up raw material import data from Africa to France and China since 1989 
df5 <- df[which(df$Year > 1989),]
df5 <- df5[which(df5$Exp %in% c(countrynames[Colonial])),]
df5 <- df5[which(df5$Imp %in% countrynames[Africa]),]
df5 <- df5 %>% group_by(Year,Indicator,Exp) %>% summarise(total = sum(values))
df5$Indicator <- as.factor(df5$Indicator)
levels(df5$Indicator)

#Plots African imports from China and France by type of raw material
df6 <- df5[which(df5$Indicator %in% levels(df5$Indicator)[3:6]),]
Import_bars_item <- ggplot(df6, aes(x=df3$Year, y=df6$total, fill = df3$Indicator))+
  geom_bar(stat = "identity")+
  scale_fill_manual(name = "Indicator", values = c("#734C20","#9B989B","#111111","#F47C20"))+
  facet_wrap(df6$Exp~., scales = "free_y", ncol = 1)+
  ylab(NULL)+
  xlab("Year")+
  ylim(0,400000000)
Import_bars_item

#### 10. Normalise by population ####

Population$China <- as.numeric(gsub(",","",Population$China))
Population$France <- as.numeric(gsub(",","",Population$France))

#Divides raw material exports across categories per population of China and France
df10 <- df2[which(df2$Year < 2025),]
for(i in 1:nrow(df10)){
  if(df10$Imp[i] %in% "China"){
    df10$percapita[i] <- df10$total[i]/Population$China[which(Population$Year %in% df10$Year[i])]  
  }else
  {df10$percapita[i] <- df10$total[i]/Population$France[which(Population$Year %in% df10$Year[i])]}
}

#Divides raw material imports across categories per population of China and France
df11<- df5[which(df5$Year < 2025),]
for(i in 1:nrow(df11)){
  if(df11$Exp[i] %in% "China"){
    df11$percapita[i] <- df11$total[i]/Population$China[which(Population$Year %in% df11$Year[i])]  
  }else
  {df11$percapita[i] <- df11$total[i]/Population$France[which(Population$Year %in% df11$Year[i])]}
}

#Calculates net raw materials exports from Africa to China and France
df11$net_percapita <- df11$percapita - df10$percapita
df11$Indicator <- as.factor(df11$Indicator)

#Aggregates up raw material net flows across categories
df12 <- df11[which(df11$Indicator %in% levels(df11$Indicator)[c(3:6)]),]
df12 <- df12 %>% group_by(Year,Exp) %>% summarise(total = sum(net_percapita))

#Plots net exports from Africa per capita
Export_net_percapita <- ggplot(df12, aes(x=df12$Year, y=df12$total, colour=df12$Exp))+
  geom_line()+
  geom_point()+
  scale_colour_manual(name = "Importer", values = c("#de0008","#1034ad"))+
  ylab(NULL)+
  xlab("Year")
Export_net_percapita

#### 11. Calculate dependency ####

#Calculates percentage of raw materials consumed in China originating in China, Africa, and rest of the world
df13 <- df[which(df$Year == 2023),]
df13 <- df13[which(df13$Imp %in% "China"),]
df13$Indicator <- as.factor(df13$Indicator)
df13 <- df13[which(df13$Indicator %in% levels(df13$Indicator)[c(3:6)]),]
df13 <- df13 %>% group_by(Year,Exp,Imp) %>% summarise(total = sum(values))
df13$total[which(df13$Exp %in% "China")]
sum(df13$total[which(df13$Exp %in% c(countrynames[Africa]))])
sum(df13$total)

#Calculates percentage of raw materials consumed in France originating in China, Africa, and rest of the world
df14 <- df[which(df$Year == 2023),]
df14 <- df14[which(df14$Imp %in% "France"),]
df14$Indicator <- as.factor(df14$Indicator)
df14 <- df14[which(df14$Indicator %in% levels(df14$Indicator)[c(3:6)]),]
df14 <- df14 %>% group_by(Year,Exp,Imp) %>% summarise(total = sum(values))
df14$total[which(df14$Exp %in% "France")]
sum(df14$total[which(df14$Exp %in% c(countrynames[Africa]))])
sum(df14$total)

#### 12. Calculate trade balance per capita ####

#Selects relevant data rows
OEC <- OEC[,c(-4,-7)]

#Subtracts imports from exports and divides by the population of China
OEC$NetChina <- OEC$Imports.to.China..bn.USD.-OEC$Exports.from.China..bn.USD.
OEC$NetChinapercapita <- OEC$NetChina/Population$China[30:4]

#Subtracts imports from exports and divides by the population of France
OEC$NetFrance <- OEC$Imports.to.France..bn.USD.-OEC$Exports.from.France..bn.USD.
OEC$NetFrancepercapita <- OEC$NetFrance/Population$France[30:4]

#Repackages a dataframe for subsequent plotting
OEC <- cbind.data.frame(rep(OEC$Year,2),c(rep("China", nrow(OEC)),rep("France", nrow(OEC))), c(OEC$NetChinapercapita,OEC$NetFrancepercapita))
colnames(OEC) <- c("Year","Country","Value")

#Plots the net trade balance per capita between Africa and China, and between Africa and France
OEC_plot <- ggplot(OEC, aes(x=OEC$Year, y=OEC$Value, colour=OEC$Country))+
  geom_point()+
  geom_smooth() #Adds a smoothing line for visualisation
OEC_plot

#### 13. Save graphs and data ####

ggsave(Category, filename = "Category.svg")
ggsave(CategoryAfrica, filename = "CategoryAfrica.svg")
ggsave(Subregion, filename = "Subregion.svg")
ggsave(SubregionAfrica, filename = "SubregionAfrica.svg")
ggsave(Colony, filename = "Colony.svg")
ggsave(ColonyAfrica, filename = "ColonyAfrica.svg")

ggsave(Export_bars_item, filename = "Export_bars_item.svg")
ggsave(Import_bars_item, filename = "Import_bars_item.svg")
ggsave(Export_net_percapita, filename = "Export_net_percapita.svg")
ggsave(OEC_plot, filename = "OEC_plot.svg")
#All graphs are subsequently formatted in Adobe Illustrator to generate the publication figures

Export <- cbind.data.frame(AfricaFC$Conflict.Id,AfricaFC$Case,AfricaFC$Country,AfricaFC$Subregion,AfricaFC$frique,AfricaFC$First.level.category,AfricaFC$Commodity,AfricaFC$Intensity.of.Conflict,AfricaFC$Project.Status,AfricaFC$Company)
write.csv(Export, "Cases Françafrique-Chinafrique.csv")
write.csv(Africa, "Africa.csv")
#These files are used as input for the map in QGIS

write.csv(df, "EUE.csv")