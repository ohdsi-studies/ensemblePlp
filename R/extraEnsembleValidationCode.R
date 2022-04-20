# function to validate the models across other databases 
validateModels <- function(
  databaseDetails = list(
    list(
      name = 'ccae',
      connectionDetailList = list(),
      cdmDatabaseSchema = '',
      tempEmulationSchema = '',
      cohortDatabaseSchema = '',
      cohortTable = '',
      outcomeDatabaseSchema = '',
      outcomeTable = ''
      ),
    ist(name = 'mdcd'),
    ist(name = 'mdcr'),
    ist(name = 'optumClaims'),
    ist(name = 'optumEhr')
  ),
  modelLocation = 'D:/ensembleDatabase',
  cohortTable = 'databaseEnsembleCohort'
){
  
  databases <- unlist(lapply(databaseDetails, function(x) x$name))


for(ind in 1:length(databases)){
  
  modelDatabase <- databases[ind]

  for(valDatabase in databases[databases!=modelDatabase]){
    
    connectionDetails <- do.call(
      DatabaseConnector::createConnectionDetails, 
      databaseDetails[[which(databases == valDatabase)]]$connectionDetailList
      )
    
    for(analysisId in  1:21){
      
      res <- PatientLevelPrediction::loadPlpResult(file.path(modelLocation,modelDatabase,paste0('Analysis_',analysisId), 'PlpResult'))
      pop <- readRDS(file.path(modelLocation,valDatabase,paste0('StudyPop_L1_T8931_O',res$model$trainDetails$outcomeId,'_P1.rds')))

      plpData <- PatientLevelPrediction::getPlpData(
        databaseDetails = PatientLevelPrediction::createDatabaseDetails(
          connectionDetails = connectionDetails,
          cdmDatabaseSchema = databaseDetails[[which(databases == valDatabase)]]$cdmDatabaseSchema,
          cdmDatabaseName = valDatabase,
          tempEmulationSchema = databaseDetails[[which(databases == valDatabase)]]$tempEmulationSchema,
          cohortDatabaseSchema = databaseDetails[[which(databases == valDatabase)]]$cohortDatabaseSchema,
          cohortTable = databaseDetails[[which(databases == valDatabase)]]$cohortTable,
          outcomeDatabaseSchema = databaseDetails[[which(databases == valDatabase)]]$outcomeDatabaseSchema,
          outcomeTable = databaseDetails[[which(databases == valDatabase)]]$outcomeTable,
          cohortId = res$model$trainDetails$cohortId,
          outcomeIds = res$model$trainDetails$outcomeId
            ), 
        covariateSettings =  res$model$trainDetails$covariateSettings, 
        restrictPlpDataSettings = res$model$trainDetails$getPlpSettings
      )
      
      cohort <- plpData$cohorts
      md <- attr(pop, 'metaData')
      newPop <- merge(pop[,colnames(pop)!='rowId'], cohort[,c('rowId', 'subjectId', 'cohortStartDate')], by = c('subjectId', 'cohortStartDate'))
      attr(newPop, 'metaData') <- md
      
      exVal <- list(
        prediction = PatientLevelPrediction::predictPlp(
          plpModel =  res$model, 
          plpData = plpData, 
          population = newPop
        )
      )

      if(!dir.exists(file.path(modelLocation,modelDatabase, 'validation',valDatabase,paste0('Analysis_',analysisId)))){
        dir.create(file.path(modelLocation,modelDatabase, 'validation',valDatabase,paste0('Analysis_',analysisId)),recursive = T )
      }
      saveRDS(exVal, file.path(modelLocation,modelDatabase, 'validation',valDatabase,paste0('Analysis_',analysisId), 'validationResult.rds' ))
    }
  }
}
  
}
