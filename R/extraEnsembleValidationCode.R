# externally validate the models in the new data:
options(andromedaTempFolder = "D:/andromedaTemp")

databases <- c('ccae','mdcd','mdcr','optumEhr','optumClaims')

for(modelDatabase in databases){

  for(valDatabase in databases[databases!=modelDatabase]){
    
    database <- valDatabase
    # USER INPUTS
    #=======================
    # The folder where the study intermediate and result files will be written:
    
    connectionDetails <- DatabaseConnector::createConnectionDetails(dbms = dbms, 
                                                                    user = user, 
                                                                    password = pw, 
                                                                    server = server, 
                                                                    port = port,
                                                                    extraSettings = extraSettings)
    
    # Add the database containing the OMOP CDM data
    cdmDatabaseSchema <- cdmDatabaseSchema
    # Add a sharebale name for the cdmDatabaseSchema database
    databaseName <- cdmDatabaseName
    # Add a database with read/write access as this is where the cohorts will be generated
    cohortDatabaseSchema <- cohortDatabaseSchema
    
    oracleTempSchema <- NULL
    
    # table name where the cohorts will be generated
    cohortTable <- 'databaseEnsembleCohort'
    
    
    for(analysisId in  1:21){
      
      model <- PatientLevelPrediction::loadPlpResult(file.path('D:/ensembleDatabase',modelDatabase,paste0('Analysis_',analysisId), 'PlpResult'))
      pop <- readRDS(file.path('D:/ensembleDatabase',valDatabase,paste0('StudyPop_L1_T8931_O',model$inputSetting$populationSettings$outcomeId,'_P1.rds')))

      plpData <- PatientLevelPrediction::similarPlpData(plpModel = model$model, 
                                                        newConnectionDetails = connectionDetails, 
                                                        newCdmDatabaseSchema = cdmDatabaseSchema, 
                                                        newCohortDatabaseSchema = cohortDatabaseSchema, 
                                                        newCohortTable = cohortTable, 
                                                        newCohortId = 8931, 
                                                        newOutcomeDatabaseSchema = cohortDatabaseSchema, 
                                                        newOutcomeTable = cohortTable, 
                                                        newOutcomeId = model$inputSetting$populationSettings$outcomeId, 
                                                        createPopulation = T)
      
      
      cohort <- plpData$plpData$cohorts
      md <- attr(pop, 'metaData')
      newPop <- merge(pop[,colnames(pop)!='rowId'], cohort[,c('rowId', 'subjectId', 'cohortStartDate')], by = c('subjectId', 'cohortStartDate'))
      attr(newPop, 'metaData') <- md
      exVal <- PatientLevelPrediction::applyModel(population = newPop, plpData = plpData$plpData, plpModel = model$model)
      
      
      if(!dir.exists(file.path('D:/ensembleDatabase',modelDatabase, 'validation',valDatabase,paste0('Analysis_',analysisId)))){
        dir.create(file.path('D:/ensembleDatabase',modelDatabase, 'validation',valDatabase,paste0('Analysis_',analysisId)),recursive = T )
      }
      saveRDS(exVal, file.path('D:/ensembleDatabase',modelDatabase, 'validation',valDatabase,paste0('Analysis_',analysisId), 'validationResult.rds' ))
    }
  }
}
