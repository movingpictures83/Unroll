#
#
# ADD BOOTSTRAP SCORE OUTPUT OR FILTERING?
#
#
#


library(scales)
library(stringr)
library(dplyr)
options("scipen"=100, "digits"=4)

readNetwork <- function(nameOfNetwork, bootscoreName){

  network = tryCatch({
    readChar(nameOfNetwork, file.info(nameOfNetwork)$size)
  }, warning = function(e) {
    ""
  }, error = function(e) {
    ""
  })
  if (network == "")
    return (c())
  
  
  # bootscore is d1
  
  
  #READ NETWORK
  network =  gsub("\\(", "", gsub("\\)", "",gsub("\\[", "", gsub("\\]","",gsub("t0", "ti", gsub("tn", "ti+1", tolower(unlist(strsplit(x =network, "\r?\n")))))))))
  network = gsub(paste("\n<data key=\"",bootscoreName,"\">[0-9]*\\.[0-9]*</data>\n",sep=""),"",network)
  newNetwork = network[grep(pattern = paste(".*target=\".*_ti(\\+1)?\">","",sep=""), x = network, ignore.case = T)]
  newNetwork = gsub("_ti\\+1","",gsub("\\(", "", gsub("\\)", "", gsub("\\[", "", gsub("\\]", "", gsub("\">", "", gsub("<edge id=\"e?\\d*\" source=\"", "", newNetwork)))))))
  newNetwork = gsub("_ti","",newNetwork)
  
  # newNetwork = newNetwork[setdiff(1:length(newNetwork),grep(".*hg.__.*",newNetwork))]

  edgeBootscore = network[grep(paste("<data key=\\\"",bootscoreName,".*",sep=""),network)]
  edgeBootscore = as.numeric(gsub("\"", "", gsub(paste("<data key=\\\"",bootscoreName,"\\\">",sep=""), "", gsub("</data>", "", edgeBootscore)))) 
  
  
  #Extract the information from the network
  pairs = c()
  for (i in 1:length(newNetwork)){
    edge =  c(trimws(strsplit(newNetwork[i],"\" target=\"")[[1]], which ="both"))
    pairs = rbind(pairs,edge)
  }
  pairs = cbind(pairs,edgeBootscore)
  
  return(pairs)
}

getUnrolledTGMTpairs <- function(pairs){
  unrolledpairs = c()
  #get the taxa ti+1
  taxa = unique(as.character(pairs[grep("s__.*",pairs[,1]),1]))
  for (t in taxa){
    #for each taxa, get the genes it expresses
    tinteraction = unique(as.character(pairs[grep(t,pairs[,1]),2]))
    genes =  unique(as.character(tinteraction[grep("g__.*",tinteraction)]))
    for (g in genes){
      # for each gene, get the metabolites it produces
      ginteraction = unique(as.character(pairs[grep(g,pairs[,1]),2]))
      metabolites =  unique(as.character(ginteraction[grep("m__.*",ginteraction)]))
      bootscoreTG = pairs[intersect(grep(t,pairs[,1]), grep(g,pairs[,2]))[1],3]
      #now try to see if any of these metabolites itneract with the taxa
      for (m in metabolites){
        minteraction = unique(as.character(pairs[grep(m,pairs[,1]),2]))
        taxainfluenced = unique(as.character(minteraction[grep("s__.*",minteraction)]))
        bootscoreGM = pairs[intersect(grep(g,pairs[,1]), grep(m,pairs[,2]))[1],3]
        for (ti in taxainfluenced){
          # print(paste(t,"->",g,"->",m,"->",ti))
          # write(paste(".*",t,".*",ti), outputFile,  append=TRUE)
          # write(paste(t,"->",g,"->",m,"->",ti), outputFile,  append=TRUE)
          
          #Get the bootscore of each
          bootscoreMTi = pairs[intersect(grep(m,pairs[,1]), grep(ti,pairs[,2]))[1],3]
          unrolledpairs = rbind(unrolledpairs,c(t,"",g,m,ti,bootscoreTG,bootscoreGM,bootscoreMTi))
          
        }
      }
    }
  }
  if (length(unrolledpairs) > 0)
    colnames(unrolledpairs) = c("t","","g","m","ti","bTG","bGM","bMTi")
  return(unrolledpairs)
}

getUnrolledTGTpairs <- function(pairs){
  unrolledpairs = c()
  #get the taxa ti+1
  taxa = unique(as.character(pairs[grep("s__.*",pairs[,1]),1]))
  for (t in taxa){
    #for each taxa, get the genes it expresses
    tinteraction = unique(as.character(pairs[grep(t,pairs[,1]),2]))
    genes =  unique(as.character(tinteraction[grep("g__.*",tinteraction)]))
    for (g in genes){
      # for each gene, get the metabolites it produces
      ginteraction = unique(as.character(pairs[grep(g,pairs[,1]),2]))
      taxainfluenced =  unique(as.character(ginteraction[grep("s__.*",ginteraction)]))
      bootscoreTG = pairs[intersect(grep(t,pairs[,1]), grep(g,pairs[,2]))[1],3]
      #now try to see if any of these metabolites itneract with the taxa
      for (ti in taxainfluenced){
        # write(paste(t,"->",g,"->",ti), outputFile,  append=TRUE)
        # unrolledpairs = rbind(unrolledpairs,c(t,g,"",ti))
        
        #Get the bootscore of each
        bootscoreGTi = pairs[intersect(grep(g,pairs[,1]), grep(ti,pairs[,2]))[1],3]
        unrolledpairs = rbind(unrolledpairs,c(t,"",g,"",ti,bootscoreTG,"",bootscoreGTi))
        
      }
    }
  }
  if (length(unrolledpairs) > 0)
    colnames(unrolledpairs) = c("t","","g","m","ti","bTG","","bGTi")
  return(unrolledpairs)
}

getUnrolledTMTpairs <- function(pairs){
  unrolledpairs = c()
  #get the taxa ti+1
  taxa = unique(as.character(pairs[grep("s__.*",pairs[,1]),1]))
  for (t in taxa){
    #for each taxa, get the genes it expresses
    tinteraction = unique(as.character(pairs[grep(t,pairs[,1]),2]))
    metabolites =  unique(as.character(tinteraction[grep("m__.*",tinteraction)]))
    for (m in metabolites){
      # for each gene, get the metabolites it produces
      minteraction = unique(as.character(pairs[grep(m,pairs[,1]),2]))
      taxainfluenced =  unique(as.character(minteraction[grep("s__.*",minteraction)]))
      bootscoreTM = pairs[intersect(grep(t,pairs[,1]), grep(m,pairs[,2]))[1],3]
      #now try to see if any of these metabolites itneract with the taxa
      for (ti in taxainfluenced){
        # write(paste(t,"->",m,"->",ti), outputFile,  append=TRUE)
        # unrolledpairs = rbind(unrolledpairs,c(t,"",m,ti))
        
        #Get the bootscore of each
        bootscoreMTi = pairs[intersect(grep(m,pairs[,1]), grep(ti,pairs[,2]))[1],3]
        unrolledpairs = rbind(unrolledpairs,c(t,"","",m,ti,"",bootscoreTM,bootscoreMTi))
      }
    }
  }
  if (length(unrolledpairs) > 0)
    colnames(unrolledpairs) = c("t","","g","m","ti","","bTM","bMTi")
  return(unrolledpairs)
}

getUnrolledTTTpairs <- function(pairs){
  unrolledpairs = c()
  #get the taxa ti+1
  taxa = unique(as.character(pairs[grep("s__.*",pairs[,1]),1]))
  for (t in taxa){
    #for each taxa, get the genes it expresses
    tinteraction = unique(as.character(pairs[grep(t,pairs[,1]),2]))
    middleTaxa =  unique(as.character(tinteraction[grep("s__.*",tinteraction)]))
    for (tx in middleTaxa){
      # for each gene, get the metabolites it produces
      middleSinteraction = unique(as.character(pairs[grep(tx,pairs[,1]),2]))
      finalTaxaInfluenced =  unique(as.character(middleSinteraction[grep("s__.*",middleSinteraction)]))
      bootscoreTti = pairs[intersect(grep(t,pairs[,1]), grep(tx,pairs[,2]))[1],3]
      #now try to see if any of these metabolites itneract with the taxa
      for (ti in finalTaxaInfluenced){
        #Get the bootscore of each
        bootscoreTTiTj = pairs[intersect(grep(tx,pairs[,1]), grep(ti,pairs[,2]))[1],3]
        unrolledpairs = rbind(unrolledpairs,c(t,tx,"",m,ti,bootscoreTti,"",bootscoreTTiTj))
      }
    }
  }
  # unrolledpairs = rbind(unrolledpairs,c(t,g,m,ti,bootscoreTG,bootscoreGM,bootscoreMTi))
  
  if (length(unrolledpairs) > 0)
    colnames(unrolledpairs) = c("t","tx","g","m","ti","","bTM","bMTi")
  return(unrolledpairs)
}

intersectUnrolledAndOriginal <- function(unrolledpairs,pairsTT){
  #now try to do the itnereseciton between pairsTT and unrolledpairs
  intersected = c()
  for (i in 1:nrow(pairsTT)){
      indexOfverifiedinteraction = intersect(grep(pairsTT[i,1],unrolledpairs[,1]), grep(pairsTT[i,2],unrolledpairs[,5]))[1]
      if (!is.na(indexOfverifiedinteraction)){
        # print(unrolledpairs[indexOfverifiedinteraction,])
        intersected = rbind(intersected,c(unrolledpairs[indexOfverifiedinteraction,],pairsTT[i,3]))
      }
    # }
  }
  return(intersected)
}

getfinalTGTM <- function(intersectedAll,TGTchains,TMTchains){
  finalTGTM = c()
  for (i in 1:nrow(intersectedAll)){
    tgm = intersectedAll[i,]
    #Check if this interaction is also in TGT
    isInTGT = intersect(which(TGTchains[,1] %in% tgm[1]),intersect(which(TGTchains[,5] %in% tgm[5]),which(TGTchains[,3] %in% tgm[3])))
    #Check if this interaction is also in TMT
    isInTMT = intersect(which(TMTchains[,1] %in% tgm[1]),intersect(which(TMTchains[,5] %in% tgm[5]),which(TMTchains[,4] %in% tgm[4])))
    finalTGTM = rbind(finalTGTM,c(tgm,length(isInTGT)!=0,length(isInTMT)!=0))
  }
  colnames(finalTGTM) = c("T","Tx","G","M","Ti", "T->(G|Tx)_Bs","(T|G)->M_Bs","(M|G|Tx)->Ti_Bs" ,"T->Ti_Bs","ChainInTTxTi","ChainInTGMTi","ChainInTGTi","ChainInTMTi")
  
  # Don't remove duplicates, since they are triangles, so independent interactions
  ## finalTGTM= removeDuplicates(finalTGTM)

  finalTGTM = finalTGTM[order( finalTGTM[,11], finalTGTM[,12],finalTGTM[,13], finalTGTM[,10], finalTGTM[,1], finalTGTM[,2], finalTGTM[,3], finalTGTM[,4], finalTGTM[,5],decreasing=T),]
  # finalTGTM = finalTGTM[duplicated(finalTGTM), ]
  return(finalTGTM)
}

removeDuplicates <- function(d){
  remove = c()
  for (i in 1:(nrow(d)-1)){
    for (j in (i+1):nrow(d)){
      if (all(d[i,-c(6,7,8,9)] == d[j,-c(6,7,8,9)])){
        remove = c(remove,i)
        break
      }
    }
  }
  if (length(remove) != 0){
    d = d[-remove,]
  }
  return(d)
}

getNumberOfInteractionTypes <- function(pairs){
  
  tg = length(intersect(grep("s__.*",pairs[,1]),grep("g__.*",pairs[,2])))
  gm = length(intersect(grep("g__.*",pairs[,1]),grep("m__.*",pairs[,2])))
  mt = length(intersect(grep("m__.*",pairs[,1]),grep("s__.*",pairs[,2])))
  tt = length(intersect(grep("s__.*",pairs[,1]),grep("s__.*",pairs[,2])))
  mm = length(intersect(grep("m__.*",pairs[,1]),grep("m__.*",pairs[,2])))
  gg = length(intersect(grep("g__.*",pairs[,1]),grep("g__.*",pairs[,2])))
  gt = length(intersect(grep("g__.*",pairs[,1]),grep("s__.*",pairs[,2])))
  
  stats = c(gt,tg,gm,mt,tt,mm,gg)
  names(stats) = c("#g->t","#t->g","#g->m","#m->t","#t->t","#m->m","#g->g")
  
  return(stats)
  # return(paste("Number of t->g: ", tg, ", g->m: ", gm, ", m->t: ", mt, ", t->t: " , tt, ", m->m: ", mm, ", g->g: ", gg,sep=""))
}



main <- function(folder,TGMName,TGName,TMName,TName,TTName,outputFile, method){
  if (method == "Tigramite_"){
    bootscoreName = "d1"
  }else if (method == "PyCausal_"){
    bootscoreName = "d0"
  }else{
    bootscoreName = "key_bootscore"
  }
  
  
  pairsTGMT = readNetwork(paste(folder,TGMName,sep="/"),bootscoreName)
  pairsTT = readNetwork(paste(folder,TName,sep="/"),bootscoreName)
  pairsTGT = readNetwork(paste(folder,TGName,sep="/"),bootscoreName)
  pairsTMT = readNetwork(paste(folder,TMName,sep="/"),bootscoreName)
  pairsTTT = readNetwork(paste(folder,TTName,sep="/"),bootscoreName)
  
  interactionTypespairsTT = getNumberOfInteractionTypes(pairsTT)
  interactionTypespairsTGMT = getNumberOfInteractionTypes(pairsTGMT)
  interactionTypespairsTGT = getNumberOfInteractionTypes(pairsTGT)
  interactionTypespairsTMT = getNumberOfInteractionTypes(pairsTMT)
  interactionTypespairspairsTTT = getNumberOfInteractionTypes(pairsTTT)
  print("pairsTT")
  print(interactionTypespairsTT)
  print("pairsTGMT")
  print(interactionTypespairsTGMT)
  print("pairsTGT")
  print(interactionTypespairsTGT)
  print("pairsTMT")
  print(interactionTypespairsTMT)
  print("pairsTTT")
  print(interactionTypespairspairsTTT)
  
  TGMTchains = getUnrolledTGMTpairs(pairsTGMT)
  TGTchains = getUnrolledTGTpairs(pairsTGT)
  TMTchains = getUnrolledTMTpairs(pairsTMT)
  TTTchains = getUnrolledTTTpairs(pairsTTT)
  
  intersectedTGTM = intersectUnrolledAndOriginal(TGMTchains,pairsTT)
  intersectedTGT = intersectUnrolledAndOriginal(TGTchains,pairsTT)
  intersectedTMT = intersectUnrolledAndOriginal(TMTchains,pairsTT)
  intersectedTTT = intersectUnrolledAndOriginal(TTTchains,pairsTT)

  if (is.null(intersectedTGTM) && is.null(intersectedTGT) && is.null(intersectedTMT) &&is.null(intersectedTTT)){
    finalAll = NULL
  }else{
 
    intersectedAll= c()
    if (!is.null(intersectedTGTM)){
      intersectedAll = rbind(intersectedAll,cbind(intersectedTGTM,rep("FALSE",nrow(intersectedTGTM)),rep("TRUE",nrow(intersectedTGTM))))
    }
    if (!is.null(intersectedTGT)){
      intersectedAll = rbind(intersectedAll,cbind(intersectedTGT,rep("FALSE",nrow(intersectedTGT)),rep("FALSE",nrow(intersectedTGT))))
    }
    if (!is.null(intersectedTMT)){
      intersectedAll = rbind(intersectedAll,cbind(intersectedTMT,rep("FALSE",nrow(intersectedTMT)),rep("FALSE", nrow(intersectedTMT))))
    }
    if (!is.null(intersectedTTT)){
      intersectedAll = rbind(intersectedAll,cbind(intersectedTTT,rep("TRUE",nrow(intersectedTTT)),rep("FALSE", nrow(intersectedTTT))))
    }
      finalAll = getfinalTGTM(intersectedAll,TGTchains,TMTchains)
  }
  if (!is.null(finalAll)){
    write.table(finalAll, outputFile, sep=",", append=F, row.names = F,quote=F, col.names = T)
  }
  fractionUnrolledTGMT =if (is.null(intersectedTGTM)) 0 else (nrow(intersectedTGTM))/nrow(pairsTT)
  fractionUnrolledTGT = if (is.null(intersectedTGT)) 0 else (nrow(intersectedTGT))/nrow(pairsTT)
  fractionUnrolledTMT = if (is.null(intersectedTMT)) 0 else (nrow(intersectedTMT))/nrow(pairsTT)
  fractionUnrolledTTT = if (is.null(intersectedTTT)) 0 else (nrow(intersectedTTT))/nrow(pairsTTT)
  stats = c( if (is.null(TTTchains)) 0 else nrow(TTTchains),
              if (is.null(TGTchains)) 0 else nrow(TGTchains),
              if (is.null(TMTchains)) 0 else nrow(TMTchains),
              if (is.null(TGMTchains)) 0 else nrow(TGMTchains),
              if (is.null(intersectedTTT)) 0 else nrow(intersectedTTT),
              if (is.null(intersectedTGTM)) 0 else nrow(intersectedTGTM),
              if (is.null(intersectedTGT)) 0 else nrow(intersectedTGT),
              if (is.null(intersectedTMT)) 0 else nrow(intersectedTMT),
              fractionUnrolledTTT,fractionUnrolledTGMT,fractionUnrolledTGT,fractionUnrolledTMT)

  #TGTchains: Number of chains t1->g1->t2 in the TG network.
  #intersectedTGT: How many from  #TGTchains are actually an unrolling. This means that the interaction t1->g1->t2 appears in the form t1->t2 in the TT network. This is expanded into the individual interactions list that I generate in another file 
  #%UnrolledTGT: simply dividing  intersectedTGT /#TGTchains
    
  firstStatsNames = c("#TTTchains","#TGTchains","#TMTchains","#TGMTchains","intersectedTTT","intersectedTGTM","intersectedTGT","intersectedTMT","%UnrolledTTT","%UnrolledTGMT","%UnrolledTGT","%UnrolledTMT")
  #Augment the stasts with the number of each interaction type for each network
  stats = c(interactionTypespairsTT,interactionTypespairsTGMT,interactionTypespairsTGT,interactionTypespairsTMT,interactionTypespairspairsTTT,stats)
  names(stats) = c(paste("TT:",names(interactionTypespairsTT)),paste("TGMT:",names(interactionTypespairsTGMT)),paste("TGT:",names(interactionTypespairsTGT)),
                   paste("TMT:",names(interactionTypespairsTMT)),paste("TTT:",names(interactionTypespairspairsTTT)),firstStatsNames)
  
  
  stats = c(TGMName,stats)
  names(stats)[1]  = "Name"
  
  return(stats)
}
  


folder = "C:/Users/danir/Google Drive/RESEARCH/Causality/unrolling/networks/DBN"
TGMName = "/ibd_t_g_m_alignment_sr14d_dbnIntraBoot.graphml"
TGName = "/ibd_t_g_noalignment_sr14d_dbnIntraBoot.graphml"
TMName = "/ibd_t_m_alignment_sr14d_dbnIntraBoot.graphml"
TName = "/ibd_t_noalignment_sr14d_dbnIntraBoot.graphml"
outputFile = "finalTGTM.csv"
# finalAll = main(folder,TGMName,TGName,TMName,TName,outputFile,"DBN")
# print(finalAll)



#Do a for loop over here and store the overal info of the network and overal unrolling info and output some statistics






folder = "." 
method = "Tigramite_"; dataset = "filtered_ibd"; t = "_t"; g = "_g"; m = "_m"; alignment = "_noalignment"; sr = "_sr7d"
bs = "_bs0.1"; a = "_a0.01"; p = "_ptreshold0.05"; tau = "_T1"; test = "_test_ParCorr"; matrix = ".*"


# matrix = "HEGTM_Skeleton.*"
statsAll = c()
srList = c("_sr7d","_sr14d")
# for (method in c("PyCausal_","Tigramite_")){
for (method in c("PyCausal_","Tigramite_","")){
    if (method == "Tigramite_"){
    bs = "_bs0.1"
    testList = c("_test_ParCorr","_test_GPDC","_test_CMIknn")
    aList = c("_a0.5","_a0.1","_a0.01","_a0.001","_a0.0001")
    scoreList = c("")
    matrix = ".*"
    p = "_ptreshold0.05"
    nboots = ""
    tau = "_T1"
    sr = "_sr7d"
  }else if (method == "PyCausal_") {
    bs = ""
    scoreList = c("_scoreDiscreteMixedBicScore","_scoreConditionalGaussianBicScore","_scoreDegenerateGaussianBicScore","_scoreFisherZScore","_scoreMNLRBicScore","_scoreMVPBicScore","_scorePeterScore","_scoreSemBicScore","_scoreSemBicScoreDeterministic")
    testList = c("_testMultinomialLogisticRegressionWald","_testConditionalGaussianLRT", "DegenerateGaussianLRT","_testFisherZ","_testKci","_testMNLRLRT","_testPositiveCorr","_testSemBicTest")
    
    scoreList = c("_scoreFisherZScore")
    testList = c("_testPositiveCorr")
    aList = c("_a0.5","_a0.1","_a0.01","_a0.001","_a0.0001")
    p = ""
    sr = "_sr7d"
    nboots = "_nboots5"
    matrix = ".*"
    tau = ""
  }else{
    # DBN
    testList = c("")
    aList = c("")
    scoreList = c("")
    bs = ""
    p = ""
    tau = ""
    test =""
    sr = "_sr14d"
    matrix ="_dbnIntraBoot"
    nboots = ""
  }
  for (alignment in c("_noalignment","_alignment")){
    for (a in aList){
      for (test in testList){
        for (score in scoreList){
          for (sr in srList){
              TGMName = list.files(path = folder, pattern = paste(method,dataset,t,g,m,alignment,sr,bs,a,nboots,p,tau,score,test,matrix,sep=""))
              TGName = list.files(path = folder, pattern = paste(method,dataset,t,g,alignment,sr,bs,a,nboots,p,tau,score,test,matrix,sep=""))
              TMName = list.files(path = folder, pattern = paste(method,dataset,t,m,alignment,sr,bs,a,nboots,p,tau,score,test,matrix,sep=""))
              TName = list.files(path = folder, pattern = paste(method,dataset,t,alignment,sr,bs,a,nboots,p,tau,score,test,matrix,sep=""))
              TTName = list.files(path = folder, pattern = paste(method,dataset,t,t,alignment,sr,bs,a,nboots,p,tau,score,test,matrix,sep=""))
              
            
            # Tigramite_filtered_ibd_t_g_m_noalignment_sr7dTrim_bs0.1_a0.5_ptreshold0.05_T1_test_ParCorrHEGTM_Skeleton.matrix
            # print("Tigramite_filtered_ibd_t_g_m_noalignment_sr7d_bs0.1_a0.5_ptreshold0.05_T1_test_ParCorrHEGTM_Skeleton.matrix")
             # print(paste(method,dataset,t,g,m,alignment,sr,bs,a,nboots,p,tau,score,test,matrix,sep=""))
             print(paste(method,dataset,t,t,alignment,sr,bs,a,nboots,p,tau,score,test,matrix,sep=""))
            
            if (length(TGMName) != 1){
              print("Aborting, TGMName is:")
              print(TGMName)
              next
            }
            if (length(TGName) != 1){
              print("Aborting, TGName is:")
              print(TGName)
              next
            }
            if (length(TMName) != 1){
              print("Aborting, TMName is:")
              print(TMName)
              next
            }
            if (length(TName) != 1){
              print("Aborting, TName is:")
              print(TName)
              next
            }
            
            outputFile = paste(method,dataset,t,g,m,alignment,sr,bs,a,p,tau,test,".csv",sep="")
            
            stats = main(folder = folder,TGMName = TGMName, TGName = TGName, TMName = TMName, TName= TName,TTName = TTName, method = method, outputFile = outputFile)
            
            print(stats)
            statsAll = rbind(statsAll,stats)
          }
        }
      }
    }
  }
}
print(statsAll)
write.table(statsAll, "statsAllAugmented.csv", sep=",", append=F, row.names = F,quote=F, col.names = T)












