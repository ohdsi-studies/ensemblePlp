getAnalysisEnsembles <- function(
  analysis = 'Analysis_1', 
  databaseDetails = list(
    list(name = 'ccae'),
    ist(name = 'mdcd'),
    ist(name = 'mdcr'),
    ist(name = 'optumClaims'),
    ist(name = 'optumEhr')
  ),
  modelLocation = 'D:/ensembleDatabase',
  ensembleLocation = 'D:/ensembleDatabase/revisionsEnsemble',
  getDetails = T
  ){
  
  #get auc/sim/demo details
  if(getDetails){
  getDetails(
    databaseDetails = databaseDetails, 
    modelLocation = modelLocation, 
    analysis = analysis,  
    ensembleLocation = ensembleLocation
  )
  }
  
  aucs <- loadAucs(ensembleLocation, analysis)
  ages <- loadAges(ensembleLocation, analysis)
  
  databases <- unlist(lapply(databaseDetails, function(x) x$name))
  
  for(valData in databases){
    
    # get all predictions for valData
    predictionWithAge <- getPrediction(
      analysis = analysis, 
      databases = databases, 
      valData = valData,
      modelLocation = modelLocation
      )
    
    # predictions only
    prediction <- predictionWithAge[,colnames(predictionWithAge) %in% c('rowId','outcomeCount', databases)]

 
    databasesInc <- databases[databases != valData]
    
    prediction <- prediction[, c('rowId','outcomeCount', databasesInc)]
    predictionWithAge <- predictionWithAge[, c('rowId','outcomeCount','ageYear', databasesInc)]
    
    # make sure these are ordered the same as databasesInc
    aucsInc <- aucs[databasesInc]
    agesInc <- ages[databasesInc]
    
      # get individual database performance:
      #==============
      for(tdb in databasesInc){
        
        risk <- data.frame(
          value = prediction[,tdb],
          rowId = prediction[,'rowId'],
          outcomeCount = prediction[,'outcomeCount']
          )
        attr(risk, "metaData")$modelType <- "binary"
        level1 <- list(
          auc = PatientLevelPrediction::computeAuc(risk, confidenceInterval = T), 
          #calDF = PatientLevelPrediction::getCalibration(risk, 100),
          calInLarge = PatientLevelPrediction:::calibrationInLarge(risk),
          calGradient = PatientLevelPrediction:::calibrationWeak(risk)
        )
        
        saveRDS(level1, file.path(ensembleLocation,paste0(analysis,'_',valData,'_',tdb)))
      }
      
      
      
      # Ensembles
      #==============
      
      meanRisk <- data.frame(
        rowId = prediction[,'rowId'],
        value = apply(prediction[,databasesInc], 1, mean),
        outcomeCount = prediction[,'outcomeCount']
        )
      
      attr(meanRisk, "metaData")$modelType <- "binary"
      
      meanEnsemble <- list(
        auc = PatientLevelPrediction::computeAuc(meanRisk, confidenceInterval = T), 
        ###calDF = PatientLevelPrediction::getCalibration(meanRisk, 100),
        calInLarge = PatientLevelPrediction:::calibrationInLarge(meanRisk),
        calGradient = PatientLevelPrediction:::calibrationWeak(meanRisk)
        )
      saveRDS(meanEnsemble, file.path(ensembleLocation,paste0(analysis,'_',valData,'_mean')))
      
      
# AUC fusion
        tauc <- as.matrix(abs(aucsInc-0.5))
        aucRisk1 <- data.frame(
          rowId = prediction[,'rowId'],
          outcomeCount = prediction[,'outcomeCount'],
          value = as.matrix(prediction[,databasesInc])%*%tauc/sum(tauc)
        )
        attr(aucRisk1, "metaData")$modelType <- "binary"
  
      aucEnsemble1 <- list(
        auc = PatientLevelPrediction::computeAuc(aucRisk1, confidenceInterval = T), 
        #calDF = PatientLevelPrediction::getCalibration(aucRisk1, 100),
        calInLarge = PatientLevelPrediction:::calibrationInLarge(aucRisk1),
        calGradient = PatientLevelPrediction:::calibrationWeak(aucRisk1)
      )
      saveRDS(aucEnsemble1, file.path(ensembleLocation,paste0(analysis,'_',valData,'_auc1')))
      
      
# AUC 2
      tauc <- as.matrix(aucsInc-0.5)
      aucRisk2 <- data.frame(
        rowId = prediction[,'rowId'],
        outcomeCount = prediction[,'outcomeCount'],
        value = as.matrix(prediction[,databasesInc])%*%tauc
      )
      attr(aucRisk2, "metaData")$modelType <- "binary"
        
      aucEnsemble2 <- list(
        auc = PatientLevelPrediction::computeAuc(aucRisk2, confidenceInterval = T), 
        #calDF = PatientLevelPrediction::getCalibration(aucRisk2, 100),
        calInLarge = PatientLevelPrediction:::calibrationInLarge(aucRisk2),
        calGradient = PatientLevelPrediction:::calibrationWeak(aucRisk2)
      )
      saveRDS(aucEnsemble2, file.path(ensembleLocation,paste0(paste0(analysis,'_',valData,'_auc2'))))
        
        
        
        
        
        # covariateSummary
        
        # first get database models and restrict sum to predictors, use cosine sim
        covSum = read.csv(file.path(ensembleLocation,paste0(analysis,'_covSum_',valData,'.csv')))
        covSum <- covSum[,c('covariateId','CovariateMean')]
        
        sims <- c()
        for(database in databasesInc){
          dat = read.csv(file.path(ensembleLocation, paste0(analysis,'_covSum_',database,'.csv')))
          dat <- dat[dat$covariateValue!=0,c('covariateId','CovariateMean')]
          
          datSim <- merge(dat, covSum, by='covariateId', all.x=T)
          datSim[is.na(datSim)] <- 0
          sims <- c(sims,sum(datSim$CovariateMean.x*datSim$CovariateMean.y)/sqrt(sum(datSim$CovariateMean.x^2)*sum(datSim$CovariateMean.y^2)))
        }
        names(sims) <- databasesInc
        

          tsim <- as.matrix(sims)
          simRisk <- data.frame(
            rowId = prediction[,'rowId'],
            outcomeCount = prediction[,'outcomeCount'],
            value = as.matrix(prediction[,databasesInc])%*%t(t(tsim))/sum(tsim))
          
          attr(simRisk, "metaData")$modelType <- "binary"

        simEnsemble <- list(
          auc = PatientLevelPrediction::computeAuc(simRisk, confidenceInterval = T), 
          #calDF = PatientLevelPrediction::getCalibration(simRisk, 100),
          calInLarge = PatientLevelPrediction:::calibrationInLarge(simRisk),
          calGradient = PatientLevelPrediction:::calibrationWeak(simRisk)
          )
        saveRDS(simEnsemble, file.path(ensembleLocation,paste0(analysis,'_',valData,'_sim')))
        
        
        
        
        
# age ensemble
          ageDiff <- as.matrix(1/(1+abs(agesInc - rep(ages[valData],4))))
          ageRisk <- data.frame(
            rowId = prediction[,'rowId'],
            outcomeCount = prediction[,'outcomeCount'],
            value = as.matrix(prediction[,databasesInc])%*%ageDiff/sum(ageDiff)
            )
          attr(ageRisk, "metaData")$modelType <- "binary"

        
        ageEnsemble <- list(
          auc = PatientLevelPrediction::computeAuc(ageRisk, confidenceInterval = T), 
          #calDF = PatientLevelPrediction::getCalibration(ageRisk, 100),
          calInLarge = PatientLevelPrediction:::calibrationInLarge(ageRisk),
          calGradient = PatientLevelPrediction:::calibrationWeak(ageRisk)
          )
        saveRDS(ageEnsemble, file.path(ensembleLocation,paste0(analysis,'_',valData,'_age')))
        
        
        
        
        
        
        # mixture of experts
        meAgeVal <- lapply(1:nrow(predictionWithAge), function(i) predictionWithAge[i,3+which.min(abs(agesInc-predictionWithAge$ageYear[i]))])
        predictionWithAge$value <- unlist(meAgeVal)
        attr(predictionWithAge, "metaData") <- list(modelType = "binary") 
        
        ageME <- list(
          auc = PatientLevelPrediction::computeAuc(predictionWithAge, confidenceInterval = T), 
          #calDF = PatientLevelPrediction::getCalibration(predictionWithAge, 100),
          calInLarge = PatientLevelPrediction:::calibrationInLarge(predictionWithAge),
          calGradient = PatientLevelPrediction:::calibrationWeak(predictionWithAge)
        )
        saveRDS(ageME, file.path(ensembleLocation,paste0(analysis,'_',valData,'_ageME')))
        
        
        
        
        
        
        # super learner:
        dataF <- prediction[,-1]
        colnames(dataF)[1] <- 'outcome'
        lrMod <- glm(outcome ~ ., data = dataF, family = "binomial")
        slRisk <- data.frame(value = predict(lrMod, prediction[,-c(1:2)], type = "response"),
          rowId = prediction[,1],
          outcomeCount = prediction[,2])
        attr(slRisk, "metaData")$modelType <- "binary"
        
        slEnsemble <- list(
          auc = PatientLevelPrediction::computeAuc(slRisk, confidenceInterval = T), 
          #calDF = PatientLevelPrediction::getCalibration(slRisk, 100),
          calInLarge = PatientLevelPrediction:::calibrationInLarge(slRisk),
          calGradient = PatientLevelPrediction:::calibrationWeak(slRisk),
          mod = lrMod
          )
        saveRDS(slEnsemble, file.path(ensembleLocation,paste0(analysis,'_',valData,'_slAll')))
        
        ind <- sample(nrow(dataF), 1000)
        lrMod <- glm(outcome ~ ., data = dataF[ind,], family = "binomial")
        slRisk <- data.frame(value = predict(lrMod, prediction[-ind,-c(1:2)], type = "response"),
          rowId = prediction[-ind,1],
          outcomeCount = prediction[-ind,2])
        attr(slRisk, "metaData")$modelType <- "binary"
        
        slEnsemble <- list(
          auc = PatientLevelPrediction::computeAuc(slRisk, confidenceInterval = T), 
          #calDF = PatientLevelPrediction::getCalibration(slRisk, 100),
          calInLarge = PatientLevelPrediction:::calibrationInLarge(slRisk),
          calGradient = PatientLevelPrediction:::calibrationWeak(slRisk),
          mod = lrMod
          )
        saveRDS(slEnsemble, file.path(ensembleLocation,paste0(analysis,'_',valData,'_sl1000')))
        
        ind <- sample(nrow(dataF), 10000)
        lrMod <- glm(outcome ~ ., data = dataF[ind,], family = "binomial")
        slRisk <- data.frame(value = predict(lrMod, prediction[-ind,-c(1:2)], type = "response"),
          rowId = prediction[-ind,1],
          outcomeCount = prediction[-ind,2])
        attr(slRisk, "metaData")$modelType <- "binary"
        
        slEnsemble <- list(auc = PatientLevelPrediction::computeAuc(slRisk, confidenceInterval = T), 
          #calDF = PatientLevelPrediction::getCalibration(slRisk, 100),
          calInLarge = PatientLevelPrediction:::calibrationInLarge(slRisk),
          calGradient = PatientLevelPrediction:::calibrationWeak(slRisk),
          mod = lrMod)
        saveRDS(slEnsemble, file.path(ensembleLocation,paste0(analysis,'_',valData,'_sl10000')))
        
  } # val data
}



# Helpers

getDetails <- function(
  databaseDetails, 
  modelLocation = 'D:/ensembleDatabase', 
  analysis,  
  ensembleLocation = 'D:/ensembleDatabase/emsemblePred'
){
  
  databases <- unlist(lapply(databaseDetails, function(x) x$name))
  
  for(i in 1:length(databases)){
    database <- databases[i]
    
    #model[[i]] <- PatientLevelPrediction::loadPlpResult(file.path(modelLocation, databases[i], analysis, 'plpResult'))
    #inputSetting <- model[[i]]$inputSetting
    #performanceEvaluation <- model[[i]]$performanceEvaluation
    #covariateSummary <- model[[i]]$covariateSummary
    inputSetting <- readRDS(file.path(modelLocation, databases[i], analysis, 'plpResult', 'inputSetting.rds'))
    performanceEvaluation <- readRDS(file.path(modelLocation, databases[i], analysis, 'plpResult', 'performanceEvaluation.rds'))
    covariateSummary <- readRDS(file.path(modelLocation, databases[i], analysis, 'plpResult', 'covariateSummary.rds'))
    pop<- readRDS(file.path(modelLocation,database,paste0('StudyPop_L1_T8931_O',inputSetting$populationSettings$outcomeId,'_P1.rds')))
    
    demoInfo <- pop %>% dplyr::select(ageYear, gender) %>% summarise(meanAge=mean(ageYear), 
      maleFrac = mean(gender==8507))
    write.csv(demoInfo, file.path(ensembleLocation, paste0(analysis,'_demo_',database,'.csv')))
    
    res <- as.data.frame(performanceEvaluation$evaluationStatistics)
    write.csv(res$Value[res$Eval=='test' & res$Metric=='AUC.auc'], file.path(ensembleLocation, paste0(analysis,'_auc_',database,'.csv')))
    
    write.csv(covariateSummary,file.path(ensembleLocation, paste0(analysis,'_covSum_',database,'.csv')))
  }
  
  
}

loadAges <- function(resultLocation, analysis){
  agefiles <- dir(resultLocation, pattern = paste0(analysis,'_demo_'))
  
  ages <- data.frame(
    
    unlist(
      lapply(agefiles, 
        function(x){
          read.csv(file.path(resultLocation,x))[,2]
        }
      )
    )
  )
  
  ages <- t(ages)
  colnames(ages) <- gsub('.csv','',gsub(paste0(analysis,'_demo_'), '',agefiles ))
  rownames(ages) <- NULL
  names(ages) <- colnames(ages)
  return(ages)
}

loadAucs <- function(resultLocation, analysis){
  aucfiles <- dir(resultLocation, pattern = paste0(analysis,'_auc_'))
  
  aucs <- data.frame(
    
    unlist(
      lapply(aucfiles, 
        function(x){
          read.csv(file.path(resultLocation,x))[,2]
        }
      )
    )
  )
  
  aucs <- t(aucs)
  colnames(aucs) <- gsub('.csv','',gsub(paste0(analysis,'_auc_'), '',aucfiles ))
  rownames(aucs) <- NULL
  names(aucs) <- colnames(aucs)
  return(aucs)
}



getPrediction <- function(analysis, databases, valData, modelLocation ){
  
  predWithAge <- NULL 
  for(ensembleDB in databases[databases!=valData]){
    
    # val results:
    res <- readRDS(file.path(modelLocation,ensembleDB,'Validation',valData, analysis ,'validationResult.rds'))
    
    if(is.null(predWithAge)){
      predWithAge <- res$prediction[,c('rowId','subjectId', 'cohortStartDate', 'ageYear', 'outcomeCount', 'value')]
      colnames(predWithAge)[6] <- ensembleDB
    } else{
      predt <- res$prediction[,c('subjectId', 'cohortStartDate', 'value')]
      colnames(predt)[3] <- ensembleDB
      predWithAge <- merge(predWithAge, predt, by = c('subjectId', 'cohortStartDate'))
    }
    
  }
  
  return(predWithAge)
  
}




