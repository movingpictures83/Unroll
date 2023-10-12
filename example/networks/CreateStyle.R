#########################################################################################################################
#########################################################################################################################
###### PROJECT:        DBNs
###### NAME:           CreateStyle.R
###### AUTHOR:         Daniel Ruiz-Perez, PhD Student
###### AFFILIATION:    Florida International University
###### 
###### DESCRIPTION:    This file creates the style.xml file needed to visualize it in Cytoscape. It takes as input a network
######                , a file with the mean abundance and a base style file to update.
#########################################################################################################################
#########################################################################################################################


library(scales)
library(stringr)
options("scipen"=100, "digits"=4)

folder = "C:/Users/danir/Google Drive/FIU/RESEARCH/DBNs/Heterogeneous/Visualization/Grant/" #"Demo" #"finalFiguresBootscore" #"finalFiguresBootscore"
setwd(folder)

#nameOfNetwork = "DemoFigure.graphml"#list.files(pattern =".*\\.graphml")
nameOfNetwork = "human_ibd_microbiota_metabolites+genes_dbn_sample_alignment_sr14d_dbnIntraBoot.graphml"#list.files(pattern =".*\\.graphml")
nameOfNetwork = paste(folder,nameOfNetwork,sep="/")
abundances = unlist(read.csv('infant_oral_microbiota_dbn_sample_average.txt',sep="\t"))


#New infant gut
names = c("s__Escherichia_coli",
          "s__Haemophilus_parainfluenzae",
          "s__Neisseria_subflava",
          "s__Veillonella_dispar",
          "s__Akkermansia_muciniphila",
          "s__Rothia_mucilaginosa",
          "s__Prevotella_melaninogenica",
          "s__Lactobacillus_iners")


##non-taxa for infant
attr = c("")

clinical = c("Week sample obtained")


all = c(names, attr, clinical)

f = readChar("styleBase.xml", file.info("styleBase.xml")$size)

split = strsplit(f, "<dependency name=\"nodeSizeLocked\" value=\"true\"/>", fixed = T)
before = unlist(split)[1]
after = unlist(split)[2]





################## NODE COLOR -> timepoint
aux ="\n            <visualProperty name=\"NODE_FILL_COLOR\" default=\"#ff9232\">\n                <discreteMapping attributeType=\"string\" attributeName=\"name\">"
for (i in 1:length(all)){
  aux = paste(aux,"\n                    <discreteMappingEntry value=\"#4d93a8\" attributeValue=\"", paste(all[i],"_ti",sep=""),"\"/>",sep="")
}
aux = paste(aux,"\n                </discreteMapping>\n            </visualProperty>" ,sep="")
f = paste(before,"<dependency name=\"nodeSizeLocked\" value=\"true\"/>",aux,after,sep= "")




############### SIZE -> incoming edges
network = readChar(nameOfNetwork, file.info(nameOfNetwork)$size)
network= unlist(strsplit(x =network, "\r\n"))
network= unlist(strsplit(x =network, "\n"))
maxAcum = 1
for (i in 1:length(names))
  maxAcum = max(maxAcum, length(grep(pattern = paste(".*target=\"",names[i],"_ti+1","\".*",sep=""), x = network, ignore.case = T)))

aux =paste(aux,"\n            <visualProperty name=\"NODE_SIZE\" default=\"40\">\n                <discreteMapping attributeType=\"string\" attributeName=\"name\">",sep="")
for (i in 1:length(names)){
  aux = paste(aux,"\n                    <discreteMappingEntry value=\"",(20*length(grep(pattern = paste(".*target=\"",names[i],"_ti","\".*",sep=""), x = network, ignore.case = T))/maxAcum+40),"\" attributeValue=\"",names[i],"_ti","\"/>",sep="")
  aux = paste(aux,"\n                    <discreteMappingEntry value=\"",(20*length(grep(pattern = paste(".*target=\"",names[i],"_ti\\+1","\".*",sep=""), x = network, ignore.case = T))/maxAcum+40),"\" attributeValue=\"",names[i],"_ti+1","\"/>",sep="")
}
aux = paste(aux,"\n                </discreteMapping>\n            </visualProperty>" ,sep="")
f = paste(before,"<dependency name=\"nodeSizeLocked\" value=\"true\"/>",aux,after,sep= "")



############### TRANSPARENCY -> abundance
abundancesScaled = rescale(abundances,c(90,255))
abundancesScaled = c(rep(200,length(attr)),abundancesScaled)
aux =paste(aux,"\n            <visualProperty name=\"NODE_TRANSPARENCY\" default=\"255\">\n                <discreteMapping attributeType=\"string\" attributeName=\"name\">",sep="")
for (i in 1:length(names)){
  aux = paste(aux,"\n                    <discreteMappingEntry value=\"",round(abundancesScaled[i],1),"\" attributeValue=\"",names[i],"_ti","\"/>",sep="")
  aux = paste(aux,"\n                    <discreteMappingEntry value=\"",round(abundancesScaled[i],1),"\" attributeValue=\"",names[i],"_ti+1","\"/>",sep="")
}
aux = paste(aux,"\n                </discreteMapping>\n            </visualProperty>" ,sep="")
f = paste(before,"<dependency name=\"nodeSizeLocked\" value=\"true\"/>",aux,after,sep= "")



############### SHAPE -> type of data

aux =paste(aux,"\n            <visualProperty name=\"NODE_SHAPE\" default=\"ELLIPSE\">\n                <discreteMapping attributeType=\"string\" attributeName=\"name\">",sep="")
for (i in 1:length(attr)){
  aux = paste(aux,"\n                    <discreteMappingEntry value=\"DIAMOND\" attributeValue=\"",attr[i],"_ti","\"/>",sep="")
  aux = paste(aux,"\n                    <discreteMappingEntry value=\"DIAMOND\" attributeValue=\"",attr[i],"_ti+1","\"/>",sep="")
}
for (i in 1:length(clinical)){
  aux = paste(aux,"\n                    <discreteMappingEntry value=\"TRIANGLE\" attributeValue=\"",clinical[i],"_ti","\"/>",sep="")
  aux = paste(aux,"\n                    <discreteMappingEntry value=\"TRIANGLE\" attributeValue=\"",clinical[i],"_ti+1","\"/>",sep="")
}
aux = paste(aux,"\n                </discreteMapping>\n            </visualProperty>" ,sep="")
f = paste(before,"<dependency name=\"nodeSizeLocked\" value=\"true\"/>",aux,after,sep= "")


split = strsplit(f, "<dependency name=\"arrowColorMatchesEdge\" value=\"false\"/>", fixed = T)
before = unlist(split)[1]
after = unlist(split)[2]



################## EDGE LINE TYPE -> intra-inter
aux = "\n            <visualProperty name=\"EDGE_LINE_TYPE\" default=\"SOLID\">\n                <discreteMapping attributeType=\"string\" attributeName=\"shared name\">"
for (i in 1:length(all)){
   for (j in 1:length(all)){
      aux = paste(aux,"\n                    <discreteMappingEntry value=\"EQUAL_DASH\" attributeValue=\"",all[i],"_ti+1 (-) ", all[j],"_ti+1","\"/>",sep="")
    aux = paste(aux,"\n                    <discreteMappingEntry value=\"EQUAL_DASH\" attributeValue=\"",all[i],"_ti (-) ", all[j],"_ti","\"/>",sep="")
  }
}
aux = paste(aux,"\n                </discreteMapping>\n            </visualProperty>" ,sep="")


################## EDGE PAINT
# aux =paste(aux,"\n            <visualProperty name=\"EDGE_UNSELECTED_PAINT\" default=\"#CC0033\">\n                <discreteMapping attributeType=\"string\" attributeName=\"name\">",sep="")
# for (i in 1:length(names)){
#     aux = paste(aux,"\n                    <discreteMappingEntry value=\"#FF9999\" attributeValue=\"",names[i],"_ti (-) ", names[i],"_ti+1","\"/>",sep="")
#   }
# 
# aux = paste(aux,"\n                </discreteMapping>\n            </visualProperty>" ,sep="")


# ##### Make the self loops more transparent
# aux =paste(aux,"\n            <visualProperty name=\"EDGE_TRANSPARENCY\" default=\"255\">\n                <discreteMapping attributeType=\"string\" attributeName=\"name\">",sep="")
# for (i in 1:length(names)){
#     aux = paste(aux,"\n                    <discreteMappingEntry value=\"70\" attributeValue=\"",names[i],"_ti (-) ", names[i],"_ti+1","\"/>",sep="")
#   }
# 
# aux = paste(aux,"\n                </discreteMapping>\n            </visualProperty>" ,sep="")
# 









################ Now transparency is based on the boot score
aux =paste(aux,"\n            <visualProperty name=\"EDGE_TRANSPARENCY\" default=\"255\">\n                <discreteMapping attributeType=\"string\" attributeName=\"name\">",sep="")
# for (i in 1:length(all)){
#   aux = paste(aux,"\n                    <discreteMappingEntry value=\"70\" attributeValue=\"",all[i],"_ti (-) ", all[i],"_ti+1","\"/>",sep="")
# }

aux = paste(aux,"\n                </discreteMapping>\n            </visualProperty>" ,sep="")


# 
# ################## EDGE TRANSPARENCY
# weights = str_extract(network[grep("key=\"key_bootScore",network)], "\\-*\\d+\\.+\\d+")
# weights = weights[!is.na(weights)]
# maxAbsWeight =  max(abs(as.numeric(weights)))
# minAbsWeight =  min(abs(as.numeric(weights)))
# maxAbsWeight = 150 #when we hide intra edges we need this
# 
# ######## For edge coefficient
# aux =paste(aux,"\n<visualProperty name=\"EDGE_TRANSPARENCY\" default=\"2.0\">\n")
# aux =paste(aux,"  <continuousMapping attributeType=\"float\" attributeName=\"bootScore\">\n")
# aux =paste(aux,"  <continuousMappingPoint lesserValue=\"130.0\" greaterValue=\"130.0\" equalValue=\"130.0\" attrValue=\"",minAbsWeight,"\"/>\n",sep="")
# aux =paste(aux,"  <continuousMappingPoint lesserValue=\"130.0\" greaterValue=\"255.0\" equalValue=\"255.0\" attrValue=\"",maxAbsWeight,"\"/>\n",sep="")
# aux =paste(aux,"  </continuousMapping>\n")
# aux =paste(aux,"  </visualProperty>\n")



################## EDGE WIDTH
weights = str_extract(network[grep("key=\"key_weight",network)], "\\-*\\d+\\.+\\d+")
weights = weights[!is.na(weights)]

maxAbsWeight =  max(abs(as.numeric(weights)))
medianAbsWeight =  median(abs(as.numeric(weights)))
# maxAbsWeight = 150 #when we hide intra edges we need this



######## For edge coefficient normalized
aux =paste(aux,"\n<visualProperty name=\"EDGE_WIDTH\" default=\"2.0\">\n")
aux =paste(aux,"  <continuousMapping attributeType=\"float\" attributeName=\"weight\">\n")
aux =paste(aux,"  <continuousMappingPoint lesserValue=\"15.0\" greaterValue=\"15.0\" equalValue=\"15.0\" attrValue=\"-",1,"\"/>\n",sep="")
aux =paste(aux,"  <continuousMappingPoint lesserValue=\"1.0\" greaterValue=\"1.0\" equalValue=\"1.0\" attrValue=\"0.0\"/>\n")
aux =paste(aux,"  <continuousMappingPoint lesserValue=\"15.0\" greaterValue=\"15.0\" equalValue=\"15.0\" attrValue=\"",1,"\"/>\n",sep="")
aux =paste(aux,"  </continuousMapping>\n")
aux =paste(aux,"  </visualProperty>\n")



# ######## For edge coefficient
# aux =paste(aux,"\n<visualProperty name=\"EDGE_WIDTH\" default=\"2.0\">\n")
# aux =paste(aux,"  <continuousMapping attributeType=\"float\" attributeName=\"weight\">\n")
# aux =paste(aux,"  <continuousMappingPoint lesserValue=\"10.0\" greaterValue=\"15.0\" equalValue=\"15.0\" attrValue=\"-",maxAbsWeight,"\"/>\n",sep="")
# aux =paste(aux,"  <continuousMappingPoint lesserValue=\"1.0\" greaterValue=\"10.0\" equalValue=\"10.0\" attrValue=\"-",medianAbsWeight,"\"/>\n",sep="")
# aux =paste(aux,"  <continuousMappingPoint lesserValue=\"1.0\" greaterValue=\"1.0\" equalValue=\"1.0\" attrValue=\"0.0\"/>\n")
# aux =paste(aux,"  <continuousMappingPoint lesserValue=\"1.0\" greaterValue=\"10.0\" equalValue=\"10.0\" attrValue=\"",medianAbsWeight,"\"/>\n",sep="")
# aux =paste(aux,"  <continuousMappingPoint lesserValue=\"10.0\" greaterValue=\"15.0\" equalValue=\"15.0\" attrValue=\"",maxAbsWeight,"\"/>\n",sep="")
# aux =paste(aux,"  </continuousMapping>\n")
# aux =paste(aux,"  </visualProperty>\n")
  
########### For edge confidence
# aux =paste(aux,"<visualProperty name=\"EDGE_WIDTH\" default=\"2.0\">\n")
# aux =paste(aux,"  <continuousMapping attributeType=\"float\" attributeName=\"weight\">\n")
# aux =paste(aux,"  <continuousMappingPoint lesserValue=\"16.0\" greaterValue=\"20.0\" equalValue=\"20.0\" attrValue=\"-",maxAbsWeight,"\"/>\n",sep="")
# aux =paste(aux,"  <continuousMappingPoint lesserValue=\"8.0\" greaterValue=\"16.0\" equalValue=\"16.0\" attrValue=\"-",80.0,"\"/>\n",sep="")
# aux =paste(aux,"  <continuousMappingPoint lesserValue=\"8.0\" greaterValue=\"8.0\" equalValue=\"8.0\" attrValue=\"-",20.0,"\"/>\n",sep="")
# aux =paste(aux,"  <continuousMappingPoint lesserValue=\"1.0\" greaterValue=\"1.0\" equalValue=\"1.0\" attrValue=\"0.0\"/>\n")
# aux =paste(aux,"  <continuousMappingPoint lesserValue=\"8.0\" greaterValue=\"8.0\" equalValue=\"8.0\" attrValue=\"",20.0,"\"/>\n",sep="")
# aux =paste(aux,"  <continuousMappingPoint lesserValue=\"8.0\" greaterValue=\"16.0\" equalValue=\"16.0\" attrValue=\"",80.0,"\"/>\n",sep="")
# aux =paste(aux,"  <continuousMappingPoint lesserValue=\"16.0\" greaterValue=\"20.0\" equalValue=\"20.0\" attrValue=\"",maxAbsWeight,"\"/>\n",sep="")
# aux =paste(aux,"  </continuousMapping>\n")
# aux =paste(aux,"  </visualProperty>\n")


fileConn<-file("style.xml")
writeLines(paste(before,"<dependency name=\"arrowColorMatchesEdge\" value=\"false\"/>",aux,after,sep=""), fileConn)
close(fileConn)
