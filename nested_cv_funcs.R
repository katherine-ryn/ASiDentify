trainCV <- function(data, K = 5, indep_var, standardized_col, feature_col) {
  
  modelList <- list()
  prList <- list()
  testList <- list()
  coefList <- list()
  aucList <- list()
  aucprList <- list()
  
  switch(indep_var,
         autism = y <- data$`Autism susceptibility`)
  
  # Assign folds evenly using the modulus operator
  set.seed(3019)
  fold0 <- sample.int(sum(y == 0)) %% K
  fold1 <- sample.int(sum(y == 1)) %% K
  foldid <- numeric(length(y))
  foldid[y == 0] <- fold0
  foldid[y == 1] <- fold1
  foldid <- foldid + 1 
  
  for (i in 1:K) {
    # Creating training and test data sets
    trainingData <- data[which(foldid != i), ]
    testData <- data[which(foldid == i), ]
    
    # Standardize the training set
    standardizedTrain <- trainingData %>%
      dplyr::mutate_at(standardized_col, ~(scale(.) %>% as.vector))
    
    # Standardize the test set
    standardizedTest <- testData %>%
      dplyr::mutate_at(standardized_col, ~(scale(.) %>% as.vector))
    
    # Adding the test set to testList
    testList[[i]] <- standardizedTest
    
    # Find lambda - size of the penalty
    x <- as.matrix(standardizedTrain[, feature_col]) 
    switch(indep_var,
           autism = y <- standardizedTrain$`Autism susceptibility`)
    
    fraction_0 <- rep(1 - sum(y == 0) / length(y), sum(y == 0))
    fraction_1 <- rep(1 - sum(y == 1) / length(y), sum(y == 1))
    # assign that value to a "weights" vector
    weights <- numeric(length(y))
    weights[y == 0] <- fraction_0
    weights[y == 1] <- fraction_1
    
    nfold <- 5
    # assign folds evenly using the modulus operator
    fold0_inner <- sample.int(sum(y == 0)) %% nfold
    fold1_inner <- sample.int(sum(y == 1)) %% nfold
    foldid_inner <- numeric(length(y))
    foldid_inner[y == 0] <- fold0_inner
    foldid_inner[y == 1] <- fold1_inner
    foldid_inner <- foldid_inner + 1
    
    set.seed(1028)
    cvfit <- cv.glmnet(x, y, alpha = 1, family = "binomial", nfolds = nfold, foldid = foldid_inner,
                       standardize = TRUE, type.measure = "default", weights = weights)
    modelList[[i]] <- cvfit
    
    coefs <- as.matrix(coef(cvfit, s = cvfit$lambda.1se))
    ix <- which(abs(coefs[,1]) > 0) #keep the coefs which are >0
    
    coefList[[i]] <- coefs[ix, 1, drop = FALSE]
    
    set.seed(2382)
    pred <- predict(cvfit, newx = as.matrix(standardizedTest[, feature_col]),  
                    type = 'response')[,1]
    
    standardizedTest$predictions <- pred
    predicted <- standardizedTest %>% arrange(desc(predictions))
    
    prList[[i]] <- predicted
    
    switch(indep_var,
           autism = pred_labels <- predicted$`Autism susceptibility`)
    pred_ROCR <- prediction(predicted$predictions, predicted$`Autism susceptibility`)
    auc <- performance(pred_ROCR, measure = "auc")
    auc <- auc@y.values[[1]]
    aucList[[i]] <- auc
    
    aucpr <- ROCR::performance(pred_ROCR, "aucpr")
    aucpr <- aucpr@y.values[[1]]
    aucprList[[i]] <- aucpr
    
  }
  
  results <- list(models = modelList,
                  predictions = prList,
                  testSets = testList,
                  betaCoefs = coefList,
                  aucs = aucList,
                  aucprs = aucprList)
  
  class(results) <- "trainCV"
  
  return(results)
}
