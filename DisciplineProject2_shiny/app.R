#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(GEOquery)  ## go to library(devtools)
library(R.utils)
library(reshape2)
library(ggplot2)
library(limma)
library(dplyr)
library(caret)

clinical_outcome <-getGEO("GSE120396")
clinical_outcome <- clinical_outcome$GSE120396_series_matrix.txt.gz


# Output 
ui <- fluidPage(

    # Application title
    titlePanel("Classifying Rejection Predictions on Genomic Data"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            selectInput(
              inputId = "classifiers",
              label = "Choose a classifier to display:",
              choices = c("kNN" = "kNN","Decision Tree" = "DT","Random Forest" = "RF", "Linear Support Vector Machine" = "SVM", "Naive Bayes" = "NB", "Multilayer Perception" = "MLP"),
              selected = "kNN"
            ),
            sliderInput("fold",
                        "Number of Folds for Cross Validation:",
                        value = 5,
                        min = 1,
                        max = 30),
            
            sliderInput("repeats",
                        "Number of Repeats:",
                        value = 25,
                        min = 1,
                        max = 50),
            
            sliderInput("pca",
                        "Number of Principal Components:",
                        value = 15,
                        min = 1,
                        max = 88),
            
          actionButton(inputId = "update", label = "Update View")
        ),

        # Show a plot of the generated distribution
        mainPanel(
          plotOutput("classifiersPlot")
        )   
    )
)


# Define server logic required 
# Reactive
server <- function(input, output) {
  
  clinical_outcome <-getGEO("GSE120396")
  clinical_outcome <- clinical_outcome$GSE120396_series_matrix.txt.gz
  
  rejection_status  <- clinical_outcome$characteristics_ch1.1
  rejection_status <- unlist(lapply( strsplit(as.character(rejection_status), ": " ) , `[[` , 2)  )
  
  # Note: please change this dir to point to the folder where your dataset is
  datadir = "data/GSE120396_RAW/"
  
  # Read in the files
  fileNames <- list.files(datadir)
  
  # Read in all 88 files to make a table
  fileNames = list.files(datadir)
  gse = c()
  for(i in 1:length(fileNames)){
    temptable <- read.delim(file.path(datadir, fileNames[i]), header=TRUE)
    gse <- cbind(gse, temptable[,2])
    colnames(gse)[i] <- colnames(temptable)[2]
  }
  gse_pca <- prcomp(t(gse), scale. = TRUE)
  pca_data = gse_pca$x
  largevar = apply(pca_data, 1, var)
  ind = which(largevar > quantile(largevar, 0.9))
  
  
  data_input <- eventReactive(input$update, {
    switch(input$classifiers,
           "kNN" = "kNN",
           "DT" = "DT",
           "RF" = "RF",
           "SVM" = "SVM",
           "NB" = "NB",
           "MLP" = "MLP")
    return(input$classifiers)
  })
  
  pca_input = reactive({
    pca = input$pca
    return(pca)
  })
  
  output$classifiersPlot <- renderPlot({
    
    X = as.matrix(pca_data[,1:pca_input()])
    y = rejection_status
    
    ml_Data = data.frame(X)
    ml_Data$y = y
    
    classifier <- data_input()
    fold <- input$fold
    repeats <- input$repeats
    
    fitControl = trainControl(method = "repeatedcv", number = fold, repeats = repeats) 
    
    if (classifier == "kNN"){
      knnFit1 = train(factor(y)~., data = ml_Data, method = "knn", trControl = fitControl, tuneLength = 10)
      kNN_d = as.data.frame(knnFit1$resample$Accuracy)
      kNN_d %>% ggplot(aes(x = factor(0), y = `knnFit1$resample$Accuracy`)) + geom_boxplot() + 
        labs(title="Accuracy using kNN" , y = "Accuracy") + 
        theme_classic()+ theme(axis.title.x=element_blank(),
                               axis.text.x=element_blank(),
                               axis.ticks.x=element_blank())
    }
    
    
    else if (classifier == "DT"){
      rpartFit1 = train(factor(y)~., data = ml_Data, method = "rpart", trControl = fitControl)
      DT_d = as.data.frame(rpartFit1$resample$Accuracy)
      ggplot(DT_d, aes(x = factor(0),y = `rpartFit1$resample$Accuracy`)) + geom_boxplot() + 
        labs(title="Accuracy using Decision Trees" , y = "Accuracy") + 
        theme_classic()+ theme(axis.title.x=element_blank(),
                               axis.text.x=element_blank(),
                               axis.ticks.x=element_blank())
      
    }
    
    else if (classifier == "RF"){
      rfFit1 = train(factor(y)~., data = ml_Data, method = "rf", trControl = fitControl)
      RF_d = as.data.frame(rfFit1$resample$Accuracy)
      ggplot(RF_d, aes(x = factor(0),y = `rfFit1$resample$Accuracy`)) + geom_boxplot() + 
        labs(title="Accuracy using Random Forest" , y = "Accuracy") + 
        theme_classic()+ theme(axis.title.x=element_blank(),
                               axis.text.x=element_blank(),
                               axis.ticks.x=element_blank())
    }
    
    else if (classifier == "SVM"){
      SVMFit1 = train(factor(y)~., data = ml_Data, method = "svmLinear", trControl = fitControl)
      SVM_d = as.data.frame(SVMFit1$resample$Accuracy)
      ggplot(SVM_d, aes(x = factor(0),y = `SVMFit1$resample$Accuracy`)) + geom_boxplot() + 
        labs(title="Accuracy using Linear Support Vector Machine" , y = "Accuracy") + 
        theme_classic()+ theme(axis.title.x=element_blank(),
                               axis.text.x=element_blank(),
                               axis.ticks.x=element_blank())
    }
    
    else if (classifier == "NB"){
      nbFit1 = train(factor(y)~., data = ml_Data, method = "naive_bayes", trControl = fitControl)
      NB_d = as.data.frame(nbFit1$resample$Accuracy)
      ggplot(NB_d, aes(x = factor(0),y = `nbFit1$resample$Accuracy`)) + geom_boxplot() + 
        labs(title="Accuracy using Naive Bayes" , y = "Accuracy") + 
        theme_classic()+ theme(axis.title.x=element_blank(),
                               axis.text.x=element_blank(),
                               axis.ticks.x=element_blank())
    }
    
    else if (classifier == "MLP"){
      mlpFit1 = train(factor(y)~., data = ml_Data, method = "mlp", trControl = fitControl, preProc = c("center", "scale"))
      MLP_d = as.data.frame(mlpFit1$resample$Accuracy)
      ggplot(MLP_d, aes(x = factor(0),y = `mlpFit1$resample$Accuracy`)) + geom_boxplot() + 
        labs(title="Accuracy using Multilayer Perception" , y = "Accuracy") + 
        theme_classic()+ theme(axis.title.x=element_blank(),
                               axis.text.x=element_blank(),
                               axis.ticks.x=element_blank())
    }
    
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
