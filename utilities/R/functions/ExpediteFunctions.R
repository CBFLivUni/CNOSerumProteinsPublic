
# Preface ----

# The functions of this script are copyrighted to the author Eva Caamano Gutierrez


# Libraries ----
library(ggpubr)
library(variancePartition)
library(readxl)
library(limma)
library(glmnet)
library(randomForest)
library(caret)
library(varSelRF)
library(VennDiagram)
library(limma)
library(psych)
library(reshape)
library(RColorBrewer)
library(readr)
library(ggrepel)
library(verification)
library(pROC)

#variables needed


#functions needed ####



#PCA function
do_PCA_Plot_Nov16<-function(data, groups, useLabels=F, labels = "", pcs=c(1,2), type='scores', scale=T, legendName="Treatment",colores=NULL){
  
  
  if(is.null(colores)){
    colores<-c(1:length(unique(groups)))
    colores<-rainbow(length(unique(groups)))
    #colores<-c("olivedrab","orange","steelblue","darkgrey")
  }
  
  # INPUTS:
  #
  #  data - data.frame or matrix
  #   - data to analyse with variables in columns and samples in rows
  #  groups - factor
  #   - a grouping factor that determines the colours of the points in scores plots
  #  useLabels (optional) - boolean - default=FALSE
  #   - if TRUE will draw labels next to each point
  #  labels (optional) - default=""
  #   - labels to be drawn next to each point. If useLabels=TRUE and labels is empty will use rownames(data) as labels.
  #  pcs (optional) - a vector of 2 integers - default=c(1,2)
  #   - principal components to be plotted
  #  type (optional) - string - default='scores'
  #   - which plot to produce at the end. Possible: 'scores', 'loadings', 'varAcc'.
  #  scale (optional) - boolean - default=TRUE
  #   - determines if variables are to be scaled (recommended)
  #  legendName (optional) - string - default='Groups'
  #   - sets the name for the legend of colours (set by groups)
  #
  # OUTPUTS:
  # a ggplot object. If not assigned to a variable will show the plot.
  
  #here we remove the columns with 0s (this allows the scaling to work if chosen and if there are columns with 0 values)
  
  
  data<-as.matrix(data)
  
  data<-data[,!apply(data,2,function(x) all(x==0))]
  
  data<-data[,!apply(data,2,function(x) any(is.na(x)))]
  
  
  if(scale){
    pc<-prcomp(data, scale = T)
  } else {
    pc <- prcomp(data)
  }
  
  #View(pc$x)
  if(type=='scores'){
    if(useLabels & length(labels) != nrow(data)){
      print("Warning: The labels not given or given incorrectly. Using rownames.")
      labels <- rownames(data)
    }
    
    pcdf<-data.frame(pc1=pc$x[,pcs[1]], pc2=pc$x[,pcs[2]])
    
    if(useLabels) pcdf$labels<-labels
    
    perc_accounted<-summary(pc)$importance[2,pcs]*100
    #print(perc_accounted)
    
    label_offset_x <- 0.035 * (range(pcdf$pc1)[2] - range(pcdf$pc1)[1])
    label_offset_y <- 0.035 * (range(pcdf$pc2)[2] - range(pcdf$pc2)[1])
    
    .e <- environment()
    p <- ggplot(data=pcdf, aes(x=pc1, y=pc2), environment=.e) +
      geom_point(aes(fill=groups),colour="black",size=5.5,pch=21)+
      scale_fill_manual(values=colores,name=legendName)
    
    if(useLabels)  p <- p + geom_text_repel(aes(label = labels))# geom_text(aes(x=pc1+label_offset_x, y=pc2+label_offset_y, label=labels))
    
    p <- p+
      
      xlab(paste("PC",pcs[1], " (", round(perc_accounted[1],2), "%)", sep=""))+
      ylab(paste("PC",pcs[2], " (", round(perc_accounted[2],2), "%)", sep=""))+
      
      theme_bw(base_size=20)+
      theme(legend.position="bottom")
    
    
    p
    
  } else if(type=='loadings'){
    
    if(useLabels & length(labels) != nrow(pc$rotation)){
      print("Warning: loadings labels not given or given incorrectly. Using the column names.")
      labels <- colnames(data)
    }
    
    pcdf<-data.frame(load1=pc$rotation[,pcs[1]], load2=pc$rotation[,pcs[2]], var=labels)
    
    #label_offset_x <- 0.035 * (range(pcdf$load1)[2] - range(pcdf$load1)[1])
    #label_offset_y <- 0.035 * (range(pcdf$load2)[2] - range(pcdf$load2)[1])
    
    .e <- environment()
    
    p <- ggplot(data=pcdf, aes(x=load1, y=load2), environment=.e) + geom_point()
    
    if(useLabels) p <- p + geom_text_repel(aes(x=load1,y=load2),label=labels)
    
    # geom_text(aes(x=load1+label_offset_x, y=load2+label_offset_y, label=labels))
    
    
    
    #geom_text(aes(x=load1+label_offset_x, y=load2+label_offset_y, label=labels))
    
    p <- p+
      xlab(paste("Loadings for PC",pcs[1],sep=""))+
      ylab(paste("Loadings for PC",pcs[2],sep=""))+
      ggtitle("PCA loadings plot")+
      theme_bw(base_size=20)
    p
    
  } else if(type=='varAcc'){
    perc_accounted <- (pc$sdev/sum(pc$sdev)*100)
    perc_with_cumsum <- data.frame(pc = as.factor(1:length(perc_accounted)),
                                   perc_acc = perc_accounted,
                                   perc_cumsum = cumsum(perc_accounted))
    p<- ggplot(data = perc_with_cumsum, aes(x=pc, y=perc_cumsum))+
      geom_bar(stat='identity', col='black', fill='white')+
      geom_hline(yintercept = 95, col='red')+
      geom_hline(yintercept = 0, col='black')+
      xlab('PC')+
      ylab('% Variance')+
      ggtitle('% Variance accounted for by principle components')+
      theme_bw()
    print(p)
    #}else if(type=='biplot'){
    #ggsave(PCbiplot(pc),"blah.pdf")
  }else if(type=='biplot'){
    
    p<- ggbiplot(pc,choices = pcs,scale=3,obs.scale=T,var.scale=T,groups=groups,ellipse=F,
                 circle=F,labels.size=3,varname.size = 5,varname.adjust = 1.5,)+
      geom_point(size=3,aes(colour=groups),alpha=1)+
      #scale_colour_discrete(name="",values=colores)+
      scale_fill_manual(values=c("firebrick1","dodgerblue","black"))+
      #geom_point(size=3,aes(colour=groups),alpha=0.3)+
      
      #circle=T,labels=groups)+scale_colour_discrete(name="")+
      theme_bw(base_size=20)+
      theme(legend.position="none")
    
    
  }else {
    cat(sprintf("\nError: no type %s", type))
  }
  #print(p)
  return(p)
}


split_train_test = function(grp, split=split){
  grpU = unique(grp)
  out = lapply(grpU, function(x) {
    n.choose = round(split * sum(grp == x))
    sample(which(grp==x), n.choose)
  })
  do.call(c, out)
}


#function to split data in train/test----
prepData<-function(data,factor_var,split=0.8,setSeed="reproduce"){
  # for reproducibility
  if(setSeed=="reproduce"){
    set.seed(123)
    
  }
  
  
  #train_idx <- sample(seq_len(nrow(data)), size = floor(split*nrow(data)))
  train_idx <- split_train_test(grp = factor_var,split = split)
  train_data <- data[train_idx, ]
  test_data <- data[-train_idx, ]
  factor_var_train<-factor(factor_var[train_idx])
  factor_var_test<-factor(factor_var[-train_idx])
  #print(factor_var_train)
  #modified the function a posteriori to also return the indeces as that was a mistake on my part!
  out<-list("train"=as.matrix(train_data),"test"=as.matrix(test_data),"factor_train"=factor_var_train,"factor_test"=factor_var_test,"idx"=train_idx)
  
}


#function to check overlap in Univariate Stats----
getProtThre<-function(x, t=0.05,pvalCol="adj.pvalue",nameCol="somaID",AveExpr=5,logFC=0.5){
  #x is a df with the results from the univar statistics
  # t is the threshold to select significance
  #pvalCol is the name of the adj pvalue column in x
  # AveExpr- denotes if a cutoff for AveExpr is needed
  # logFC - cut off for FC
  if(nameCol=="rownames"){
    namesout<-rownames(x[x[,pvalCol]<t& x[,"AveExpr"]>AveExpr & abs(x[,"logFC"])>logFC,])
  }else{
    #namesout<-x[x[,pvalCol]<t,nameCol]
    namesout<-x[x[,pvalCol]<t & x[,"AveExpr"]>AveExpr & abs(x[,"logFC"])>logFC,nameCol]
  }
  return(namesout)
}

#function to do boxplots of selected variables


boxplotSubset<-function(data, selVars,group){
  
  #note that data is the dataset with variables in columns and samples in rows
  #selVars is the vector of names of interest
  #group is the colour group to use given in the same order than the data
  
  dataNew<-data.frame("Group"=group,data[,colnames(data)%in%selVars])
  dataNew_melt<-melt(data=dataNew,id.vars="Group")
  
  #dataNew_melt$Group<-factor(somaNew_melt$Group,ordered=T)
  
  
  pp<-ggplot(dataNew_melt)+
    geom_boxplot(aes(x=Group,y=value,fill=Group))+
    facet_wrap(~variable, scales="free_y",)+
    theme_bw()+
    xlab("")+ylab("log2 abundance")+
    scale_fill_manual(values=rainbow(length(unique(group)))) +
    theme(axis.text.x=element_text(angle=45,vjust = 0.6,hjust=0.6),legend.position = "bottom",legend.title = element_blank())
  
  print(pp)
  
  #return(pp)
}





#function to do different splits of data and store them on a list ----
#the idea is to apply this function to the train data to get slightly different slices of data
applyPrepData <- function(data_train, prepData=prepData,factor_var,split=0.9, n=10) {
  result <- list()  # Create an empty list to store the iterations
  
  for (i in 1:n) {
    iteration <- prepData(data_train,factor_var,split,setSeed = "notReproduce")  # Call the prepData function
    result[[i]] <- iteration  # Store the iteration in the result list
  }
  
  return(result)
}


#function to do multimodal Lasso selection iteratively a number of times and count the number of times each variable is selected ----
multimodalLasso<-function(train_data=train_data,factor_var_train=factor_var_train,n_runs=100,min_appearances=25){
  
  #train_data: proportion of the data to do the variable selection with
  # factor_var_train - the multimodal factor
  # n_runs - number of times we run Lasso
  # min_appearances =25 - minimum frequency to keep a var
  
  # Initialize variable appearance counter
  var_counts <- data.frame("names"=colnames(train_data), "counts"=rep(0, ncol(train_data)) )
  
  #start a list to store results
  out<-list()
  
  #find lambda
  lambda_seq <- rev(10^seq(2.25, -2.5, by = -.1))
  v_output <- cv.glmnet(train_data, 
                        factor_var_train,
                        alpha = 1, 
                        lambda = lambda_seq,
                        nfolds = 5,
                        family="multinomial",
                        type.multinomial = "ungrouped")
  best_lam <- v_output$lambda.min
  
  # Run Lasso variable selection on training data n_runs times
  for (i in 1:n_runs) {
    # Sample data for Lasso run, slightly different samples in each run
    set.seed(i)
    train_xid <- sample(seq_len(nrow(train_data)), size = floor(0.9*nrow(train_data)))
    train_x <- train_data[train_xid, ]
    factor_var_train_x<-factor_var_train[train_xid]
    
    # Fit Lasso model
    lasso_best <- glmnet(train_x, 
                         factor_var_train_x, 
                         alpha = 1, 
                         lambda = best_lam,
                         family="multinomial",
                         type.multinomial = "ungrouped")
    #pred <- predict(lasso_best, s = best_lam, newx = x_test)
    
    #we select the variables and put them in a matrix
    varsLass<-coef(object = lasso_best,s = best_lam);
    varsLassOut<-do.call(cbind, lapply(varsLass, as.matrix));
    colnames(varsLassOut)<-names(varsLass)
    
    varsLassOut_no0<-varsLassOut[!apply(varsLassOut,1,function(x) all(x==0)),]
    #note that in this multimodal approach too many are selected as it is one against the others
    
    
    # Update variable appearance counter
    selected_vars<-rownames(varsLassOut_no0)
    var_counts[var_counts[,1]%in%selected_vars,2] <- var_counts[var_counts[,1]%in%selected_vars,2]+1
    
    #store vars in a list just in case
    out[[i]]<-varsLassOut
  }
  
  # Select variables that appear at least min_appearances times
  final_selected_vars <- var_counts[var_counts[,2]>min_appearances,]
  
  xout<-list("finalSelVars"=final_selected_vars,"allVarsFreq"=var_counts ,"allRes"=out)
  
  return(xout)
  
}

#function to find variables correlated to the Lasso calculations done above----
LassoCorrelated<-function(data,selVars,thres=0.8){
  
  correDat<-cor(data)
  correDat<-correDat[,colnames(correDat)%in%selVars]
  #make the diagonal zero because if not it will be selected
  diag(correDat)<-rep(0,ncol(correDat))
  varsToKeep<-apply(X = correDat,MARGIN = 1,FUN = function(x){max(abs(x))})
  varsToKeep<-names(varsToKeep[varsToKeep>thres])
  varsToKeep<-varsToKeep[!varsToKeep%in%selVars]
  if(length(varsToKeep)==0){
    cat("There were zero proteins highly correlated")
    
  }else{
    correDat<-correDat[varsToKeep,]
    
    out<-list(correDatSel=correDat,correVars=rownames(correDat))
    return(out)
  }
}


#function to do 2-class Lasso selection iteratively a number of times and count the number of times each variable is selected ----

TwoClassLasso<-function(train_data=train_data,factor_var_train=factor_var_train,n_runs=100,min_appearances=25){
  
  #train_data: proportion of the data to do the variable selection with
  # factor_var_train - the multimodal factor
  # n_runs - number of times we run Lasso
  # min_appearances =25 - minimum frequency to keep a var
  
  # Initialize variable appearance counter
  var_counts <- data.frame("names"=colnames(train_data), "counts"=rep(0, ncol(train_data)) )
  
  #start a list to store results
  out<-list()
  
  #find lambda
  lambda_seq <- rev(10^seq(2.25, -2.5, by = -.1))
  v_output <- cv.glmnet(train_data, 
                        factor_var_train,
                        alpha = 1, 
                        lambda = lambda_seq,
                        nfolds = 5, 
                        family = "binomial")
  best_lam <- v_output$lambda.min
  
  # Run Lasso variable selection on training data n_runs times
  for (i in 1:n_runs) {
    # Sample data for Lasso run, slightly different samples in each run
    set.seed(i)
    train_xid <- sample(seq_len(nrow(train_data)), size = floor(0.9*nrow(train_data)))
    train_x <- train_data[train_xid, ]
    factor_var_train_x<-factor_var_train[train_xid]
    
    # Fit Lasso model
    lasso_best <- glmnet(train_x, factor_var_train_x, alpha = 1, lambda = best_lam,family="binomial")
    #pred <- predict(lasso_best, s = best_lam, newx = x_test)
    
    #we select the variables and put them in a matrix
    varsLass<-coef(object = lasso_best,s = best_lam);
    
    varsLass<-as.matrix(varsLass)
    #select the ones that are not zero penalised
    varsLass_no0<-varsLass[varsLass[,1]!=0,];
    #remove intercept
    varsLass_no0<-varsLass_no0[-1]
    
    #update our counter
    var_counts[var_counts[,1]%in%names(varsLass_no0),2]<-var_counts[var_counts[,1]%in%names(varsLass_no0),2]+1
    
    
    #store vars in a list just in case
    out[[i]]<-varsLass
  }
  
  # Select variables that appear at least min_appearances times
  final_selected_vars <- var_counts[var_counts[,2]>min_appearances,]
  
  xout<-list("finalSelVars"=final_selected_vars,"allVarsFreq"=var_counts ,"allRes"=out)
  
  return(xout)
  
}



#function to do backwards elimination selection iteratively a number of times and count the number of times each variable is selected
BackSelec_wTimes<-function(train_data=train_data,factor_var_train=factor_var_train,n_runs=100,min_appearances=25){
  
  #train_data: proportion of the data to do the variable selection with
  # factor_var_train - the multimodal factor
  # n_runs - number of times we run Lasso
  # min_appearances =25 - minimum frequency to keep a var
  
  # Initialize variable appearance counter
  var_counts <- data.frame("names"=colnames(train_data), "counts"=rep(0, ncol(train_data)) )
  
  #start a list to store results
  out<-list()
  
  # Run Lasso variable selection on training data n_runs times
  for (i in 1:n_runs) {
    # Sample data for Lasso run, slightly different samples in each run
    set.seed(i)
    train_xid <- sample(seq_len(nrow(train_data)), size = floor(0.9*nrow(train_data)))
    train_x <- train_data[train_xid, ]
    factor_var_train_x<-factor_var_train[train_xid]
    
    # Fit Lasso model
    rfModel_backselVars<-varSelRF(xdata=train_x,
                                  Class=factor_var_train_x, 
                                  ntree=500, 
                                  ntreeIterat=500,
                                  vars.drop.num = NULL,
                                  vars.drop.frac = 0.1)
    
    
    
    #we select the variables and put them in a matrix
    varsBack<-rfModel_backselVars$selected.vars
    
    
    # Update variable appearance counter
    
    var_counts[var_counts[,1]%in%varsBack,2] <- var_counts[var_counts[,1]%in%varsBack,2]+1
    
    #store vars in a list just in case
    out[[i]]<-varsBack
  }
  
  # Select variables that appear at least min_appearances times
  final_selected_vars <- var_counts[var_counts[,2]>min_appearances,]
  
  xout<-list("finalSelVars"=final_selected_vars,"allVarsFreq"=var_counts ,"allRes"=out)
  
  return(xout)
  
}


#function to apply one of the variable selection methods to different splits of data and store them as well as some statistics ----

applyMultiVarSelectionForCrossVal_withModel <- function(inputListwSplits,
                                                        functionSelection,
                                                        metaData,
                                                        posCase=NULL,
                                                        runLASSOSelection=TRUE,
                                                        thresMod=90,
                                                        ntree=2000,
                                                        tuneParams=FALSE,
                                                        covs2Select=NULL,
                                                        paramStepSize=2,
                                                        predModel="random forest", 
                                                        minNodeSizes=NULL,
                                                        nRFTrees=NULL,
                                                        maxDepths=NULL,
                                                        mTrys=NULL,
                                                        splitRules=NULL,
                                                        sampleFrac=NULL,
                                                        eta=NULL,
                                                        gamma=NULL,
                                                        subSample=NULL,
                                                        colSample_byTree=NULL,
                                                        lambda=NULL,
                                                        alpha=NULL,
                                                        maxDepth=NULL,
                                                        minChildWeight=NULL,
                                                        nRounds=NULL,
                                                        evalMetric=NULL,
                                                        ...) {
  
  #inputListwSplits assumes an input that will have n elements with each element containing
  # train (train data), test (test data), factor_train, factor_test and idx
  
  #functionSelection is a function to undertake variable selection a number of times. The outputs of these functions
  #must be a list that will have at least one matrix called finalSelVars with 2 columns (probes and frequencies)
  
  #ntree is for the RF model
  
  # If tuning parameters
  if (tuneParams) {
    
    # If use random forest
    if (predModel == "random forest") {
      
      # Get number of selected variables
      nVarSelected <- ncol(inputListwSplits[[1]]$train)
      
      # Define all levels
      allLevels <- unique(c(sapply(inputListwSplits, function(splitSet) splitSet$factor_test),
                            sapply(inputListwSplits, function(splitSet) splitSet$factor_train)))
      
      # Get mtrys as 1 to the double the square root of input (twice recommended mtry)
      if (is.null(mTrys)) { mTrys <- seq(2, floor(sqrt(nVarSelected))*2, paramStepSize) }
      
      # Decide if using hellinger (only works for binary)
      if (is.null(splitRules)) {
        
        # Define split rules, not hellinger if binary
        splitRules <- c("gini",
                        "extratrees",
                        ifelse(length(allLevels) > 2, NA, "hellinger"))
        splitRules <- splitRules[!is.na(splitRules)]
        
      }
      
      # If minimum node sizes
      if (is.null(minNodeSizes)) { minNodeSizes <- seq(1, ceiling(nVarSelected/2), paramStepSize) }
      
      # If number of trees
      if (is.null(nRFTrees)) { nRFTrees <- seq(1000, 3000, 1000) }
      
      # If max depths
      if (is.null(maxDepths)) { maxDepths <- c(0, seq(3, 9, paramStepSize)) }
      
      # Define combinations of mtry/RF
      comboGrid <- expand.grid(nRFTrees=nRFTrees,
                               MaxDepth=maxDepths,
                               mTry=mTrys,
                               SplitRule=splitRules,
                               MinNodeSize=minNodeSizes,
                               SampleFrac=sampleFrac,
                               stringsAsFactors=FALSE)
      
      # Report number of models
      message(paste0("// Tuning models across ",nrow(comboGrid), " configurations"))

        } else { stop("Model choice invalid") }
    
  }
  
  # Create an empty list to store the iterations
  result <- list()
  
  # Loop over splits
  for (i in seq_len(length(inputListwSplits))) {
    
    # Report
    message(paste0("// On data split = ",i))
    set.seed(1)
    
    # Get train/test data and labels
    train<-inputListwSplits[[i]]$train
    test<-inputListwSplits[[i]]$test
    factor_train<-inputListwSplits[[i]]$factor_train
    factor_test<-inputListwSplits[[i]]$factor_test
    
    # If want to run initial round of feature selection
    if (runLASSOSelection) {
      
      VarSelec<-functionSelection(as.matrix(train[,!colnames(train) %in% covs2Select]),factor_train,100,25)
      #100 is the number of runs
      # 25 is the min of variables to keep
      #print(VarSelec$)
      
      #fit a model with those Vars (at a given selection thresh) and report the stats
      varsForModel<-VarSelec$finalSelVars[VarSelec$finalSelVars[,2]>as.numeric(thresMod),1]
      
      # If we have covariates to select for
      if (!is.null(covs2Select)) { varsForModel <- union(varsForModel, covs2Select) }
      # if (i == 4) { browser() }
      
    } else {
      
      # Just take columns
      varsForModel<- colnames(train)
      
      # Create dummy VarSelec
      VarSelec <- data.frame(allVarsFreq=NA,
                             finalSelVars=varsForModel)
      
    }
    #print(varsForModel)
    if(length(varsForModel)<2){
      next
    }
    
    # If using random forest
    if (predModel == "random forest") {
      
      # Predict model with random forest
      modelRes <- RandomForestWithStats(train         = train[,colnames(train)%in%varsForModel,drop=FALSE],
                                        test          = test[,colnames(test)%in%varsForModel,drop=FALSE],
                                        metaData      = metaData,
                                        factor_train  = as.factor(factor_train),
                                        factor_test   = as.factor(factor_test),
                                        posCase       = posCase,
                                        covs2Select   = covs2Select,
                                        tuneParams    = tuneParams,
                                        comboGrid     = comboGrid,
                                        nRounds       = nRounds,
                                        splitIdx=i,
                                        ...)
      
      # Else if xgboost
    } 
    
    # Calculate MCC
    rfMCC <- mltools::mcc(as.character(modelRes$Predictions$predictions[names(factor_test)]),
                          as.character(factor_test))
    
    # Append to output
    result[[i]] <- list("data"=inputListwSplits[[i]],
                        "SampleConditions"=list(Train=factor_train,
                                                Test=factor_test),
                        "varSel"=varsForModel,
                        "VarsModelling"=varsForModel,
                        "rfmodelRes"=modelRes,
                        "rfMCC"=rfMCC,
                        "allSelVars"=VarSelec$allVarsFreq,
                        "finalSelVars"=VarSelec$finalSelVars,
                        "lassoModel"=VarSelec$out,
                        "RFModel"=modelRes$Model)
    #in the line above there is a problem and i should have reported as well all the variables seleced by adding as output VarSelec$allVarsFreq or at least VarSelec$finalSelVars, now added as an after thought
    
  }
  
  return(result)
}


RandomForestWithStats <- function(train,
                                  test,
                                  factor_train,
                                  factor_test,
                                  metaData,
                                  splitIdx=NA,
                                  predModel="random forest", 
                                  posCase=NULL,
                                  covs2Select=NULL,
                                  tuneParams=FALSE,
                                  comboGrid=NULL,
                                  ...){
  require(ranger)
  
  # Get all levels
  allLevels <- unique(c(factor_train, factor_test)) # Need this as sometimes levels missing in either
  
  # If not tuning params, just run
  if (!tuneParams) {
    
    # Set seed
    set.seed(splitIdx)
    
    # Run model
    Model <- ranger(data = cbind(Condition=factor_train[rownames(train)], data.frame(train)),
                    formula = as.formula(paste0(" Condition ~ ", paste(c(colnames(train)), collapse=" + "))),
                    always.split.variables = covs2Select,
                    ...) # Pass specified hyperparameters to model
    
    # Predict on test data
    Predictions <- predict(object = Model,
                           data = cbind(Condition=factor_test[rownames(test)],
                                        data.frame(test))) # CHANGED factor_test FROM factor_train
    
    # Add names to predictions
    names(Predictions$predictions) <- names(factor_test)
    
    # Generate confusion matrix
    ConfMat_Stats <- caret::confusionMatrix(factor(Predictions$predictions,
                                                   levels=allLevels),
                                            factor(factor_test,
                                                   levels=allLevels))
    
    if(length(levels(factor_train))>2){
      
      ConfMat_Stats_byClass<-data.frame("Class"=strsplit2(rownames(ConfMat_Stats$byClass),": ")[,2]
                                        ,ConfMat_Stats$byClass)
      
    }else{
      
      ConfMat_Stats_byClass<-paste("Not by class as 2-level, positive class is",ConfMat_Stats$positive)
      ConfMat_Stats_byClass<-ConfMat_Stats$byClass
      
    }
    
    # Return
    return(list("Model"=Model,
                "PredStats_byClass"=ConfMat_Stats_byClass,
                "PredStatsOverall"=ConfMat_Stats$overall,
                "Predictions"=Predictions))
    
    # Else
  } else {
    
    # Initialise model vectors
    models <- vector(mode="list",
                     length=nrow(comboGrid))
    confStats <- vector(mode="list",
                        length=nrow(comboGrid))
    confStatsByClass <- vector(mode="list",
                               length=nrow(comboGrid))
    
    # Loop over combinations
    for (i in seq_along(models)) {
      
      set.seed(i)
      
      # Run model
      Model <- ranger(data = cbind(Condition=factor_train[rownames(train)], data.frame(train)),
                      formula = as.formula(paste0(" Condition ~ ", paste(c(colnames(train)), collapse=" + "))),
                      num.trees = as.numeric(comboGrid[i,]$nRFTrees),
                      mtry = as.numeric(comboGrid[i,]$mTry),
                      splitrule = as.character(comboGrid[i,]$SplitRule),
                      max.depth = as.numeric(comboGrid[i,]$MaxDepth),
                      min.node.size = as.numeric(comboGrid[i,]$MinNodeSize),
                      sample.fraction = as.numeric(comboGrid[i,]$SampleFrac),
                      always.split.variables = covs2Select)
      
      # Predict
      Predictions <- predict(object = Model,
                             data = cbind(Condition=factor_test[rownames(test)],
                                          data.frame(test)), # CHANGED factor_test FROM factor_train
                             predict.all=TRUE)
      
      # If null
      if (is.null(Predictions)) {
        
        # Set model to NULL
        models[[i]] <- Predictions
        
        # Skip to next
        next
        
      }
      
      # Get prediction data
      predData <- Predictions$predictions
      
      # If multi-sample test
      if (length(factor_test) > 1) {
        
        # Initialise probability matrix
        probM <- matrix(NA,
                        nrow=nrow(predData),
                        ncol=length(allLevels))
        rownames(probM) <- names(predData)
        
        # Loop over samples
        for (j in seq_len(nrow(probM))) {
          
          # Set factor names
          perTreeVote <- as.character(allLevels[match(predData[j,], as.numeric(as.factor(allLevels)))])
          probM[j,] <- sapply(as.character(allLevels), function(class) 
            sum( perTreeVote == class ) / length(perTreeVote))
          
        }
        
        # Make dataframe
        probM <- data.frame(probM)
        
        # Get proportion of trees voted for by each class
        colnames(probM) <- as.character(allLevels)
        
        # Get predictions and max predicted probability
        predictions <- sapply(seq_len(nrow(probM)), function(i) names(which.max(probM[i,])))
        maxPredProb <- sapply(seq_len(nrow(probM)), function(i) base::max(probM[i,]))
        names(predictions) <- names(factor_test)
        names(maxPredProb) <- names(factor_test)
        
        # If binary
        if (length(allLevels) == 2) {
          
          # Calc briers
          brierScore <- Brier(probabilities=unlist(probM[posCase]),
                              truth=as.character(factor_test),
                              positive=posCase,
                              negative=setdiff(allLevels, posCase))
          
          # Calculate per-sample Brier's
          perSmplBrier <- sapply(seq_along(factor_test), function(i) Brier(probabilities=unlist(probM[i,posCase]),
                                                                           truth=as.character(factor_test[i]),
                                                                           positive=posCase,
                                                                           negative=setdiff(allLevels, posCase)))
          names(perSmplBrier) <- names(factor_test)
          
          # Calculate AUC ROC
          rocVal <- suppressMessages(roc(as.numeric(factor(factor_test, levels=c("Other", "CNO"))) - 1,
                                         as.numeric(factor(predictions, levels=c("Other", "CNO"))) - 1))
          aucVal <- pROC::auc(rocVal)
          
          # MCC val
          mccVal <- mltools::mcc(factor(factor_test, levels=c("Other", "CNO")),
                                 factor(predictions, levels=c("Other", "CNO")))
          
          # Per sample AUC/ROC is AUC/ROC
          perCondROCAUC <- lapply(seq_along(allLevels), function(i) list(ROC=rocVal,
                                                                         AUC=aucVal))
          names(perCondROCAUC) <- allLevels
          
          # If multi-group
        } else {
          
          # if (splitIdx == 7) { browser() }
          
          # Calc briers
          brierScore <- multiclass.Brier(probM,
                                         factor_test)
          
          # Calculate per-sample Brier's
          # Loop over test factors
          perSmplBrier <- sapply(seq_along(factor_test), function(i) multiclass.Brier(probM[i,],
                                                                                      factor_test[i]))
          names(perSmplBrier) <- names(factor_test)
          
          # Calculate AUC ROC
          rocVal <- multiclass.roc(response=factor_test,
                                   predictor=probM)
          aucVal <- pROC::auc(rocVal)
          
          # Loop over conditions and get per-condition ROC/AUC
          perCondROCAUC <- lapply(seq_along(allLevels), function(i) {
            
            # Get per-condition test
            factor_test <- factor(factor_test, levels=allLevels) # Ensure all factors present when testing
            perCondActual <- as.numeric(factor(ifelse(factor_test == allLevels[[i]], allLevels[[i]], "Other"), 
                                               levels=c("Other", allLevels[[i]]))) - 1
            perCondPredict <- as.numeric(factor(ifelse(predictions == allLevels[[i]], allLevels[[i]], "Other"), 
                                                levels=c("Other", allLevels[[i]]))) - 1
            
            # MCC
            mccVal <- mltools::mcc(perCondActual,
                                   perCondPredict)
            
            # If either class predicted
            if (length(unique(perCondActual)) > 1) {
              
              # Calculate AUC ROC
              rocVal <- suppressMessages(roc(perCondActual,
                                             perCondPredict))
              aucVal <- pROC::auc(rocVal)
              
              # Else if just one class predicted
            } else {
              
              # Else set NAs
              rocVal <- NA
              aucVal <- NA
              
            }
            
            # Return
            return(list(ROC=rocVal,
                        AUC=aucVal))
            
          })
          names(perCondROCAUC) <- allLevels
          
        }
        
        
        # If single sample test
      } else {
        
        # Map numeric factor to factor level names
        predData <- allLevels[match(predData,
                                    as.numeric(as.factor(allLevels)))]
        
        # Get proportion of trees voted for by each class
        probM <- data.frame(lapply(allLevels, function(class) sum( predData == class ) / comboGrid[i,]$nRFTrees))
        colnames(probM) <- as.character(allLevels)
        
        # Get predictions and max predicted probability
        predictions <- sapply(seq_len(nrow(probM)), function(i) names(which.max(probM[i,])))
        maxPredProb <- sapply(seq_len(nrow(probM)), function(i) base::max(probM[i,]))
        names(predictions) <- names(factor_test)
        names(maxPredProb) <- names(factor_test)
        
        # Get matrix of true scores
        trueM <- rep(0, length(allLevels))
        names(trueM) <- allLevels
        trueM[names(trueM) == factor_test] <- 1
        
        # Sort
        trueM <- trueM[names(probM)]
        
        # IF binary
        if (length(allLevels) == 2) {
          
          # Calc briers
          brierScore <- Brier(probabilities=unlist(probM)[posCase],
                              truth=as.character(factor_test),
                              positive=posCase,
                              negative=setdiff(allLevels, posCase))
          
          # Per-sample brier is this brier
          perSmplBrier <- brierScore
          names(perSmplBrier) <- names(factor_test)
          
          # Calculate AUC ROC
          rocVal <- suppressMessages(roc(as.numeric(factor(factor_test, levels=c("Other", "CNO"))) - 1,
                                         as.numeric(factor(predictions, levels=c("Other", "CNO"))) - 1))
          aucVal <- pROC::auc(rocVal)
          
          # Per sample AUC/ROC is AUC/ROC
          perCondROCAUC <- lapply(seq_along(allLevels), function(i) list(ROC=rocVal,
                                                                         AUC=aucVal))
          names(perCondROCAUC) <- allLevels
          
          # If multi-group
        } else {
          
          # Calc briers
          brierScore <- multiclass.Brier(as.matrix(probM),
                                         names(trueM[trueM == 1]))
          
          # Per-sample brier is this brier
          perSmplBrier <- brierScore
          names(perSmplBrier) <- names(factor_test)
          
          # Calculate AUC ROC
          rocVal <- multiclass.roc(response=factor_test,
                                   predictor=probM)
          aucVal <- pROC::auc(rocVal)
          
          # Loop over conditions and get per-condition ROC/AUC
          perCondROCAUC <- lapply(seq_along(allLevels), function(i) {
            
            # Get per-condition test
            perCondActual <- as.numeric(factor(ifelse(factor_test == allLevels[[i]], allLevels[[i]], "Other"), 
                                               levels=c("Other", allLevels[[i]]))) - 1
            perCondPredict <- as.numeric(factor(ifelse(predictions == allLevels[[i]], allLevels[[i]], "Other"), 
                                                levels=c("Other", allLevels[[i]]))) - 1
            
            # Calculate AUC ROC
            rocVal <- suppressMessages(roc(perCondActual,
                                           perCondPredict))
            aucVal <- pROC::auc(rocVal)
            
            # Return
            return(list(ROC=rocVal,
                        AUC=aucVal))
            
          })
          names(perCondROCAUC) <- allLevels
          
        }
        
      }
      
      # Calculate MCC
      rfMCC <- mltools::mcc(as.character(predictions[names(factor_test)]),
                            as.character(factor_test))
      
      # Get actual
      rownames(metaData) <- metaData$Sample_ID
      actualClasses <- c(metaData[names(predictions),]$Condition)
      names(actualClasses) <- names(predictions)
      
      # Add to output
      models[[i]] <- list(Params=comboGrid[i,],
                          MCC=rfMCC,
                          ROC=rocVal,
                          AUC=aucVal,
                          PerCondRoc=sapply(perCondROCAUC, function(cond) cond$ROC,
                                            simplify = FALSE, USE.NAMES = TRUE),
                          PerCondAUC=sapply(perCondROCAUC, function(cond) cond$AUC,
                                            simplify = FALSE, USE.NAMES = TRUE),
                          BrierScore=brierScore,
                          PerSampleBrier=perSmplBrier,
                          Predictions=Predictions,
                          Actual=actualClasses,
                          PredProbM=probM,
                          PredList=predictions,
                          maxPredProb=maxPredProb)
      
      # Remove model
      rm(Model)
      
      # Generate confusion matrix
      ConfMat_Stats <- caret::confusionMatrix(factor(predictions,
                                                     levels=allLevels),
                                              factor(factor_test,
                                                     levels=allLevels))
      
      # Get by-class confidence stats
      if(length(levels(factor_train))>2){
        
        ConfMat_Stats_byClass<-data.frame("Class"=strsplit2(rownames(ConfMat_Stats$byClass),": ")[,2],
                                          ConfMat_Stats$byClass)
        
      }else{
        
        ConfMat_Stats_byClass<-paste("Not by class as 2-level, positive class is",ConfMat_Stats$positive)
        ConfMat_Stats_byClass<-ConfMat_Stats$byClass
        
      }
      
      # Append confidence stats
      confStats[[i]] <- ConfMat_Stats
      confStatsByClass[[i]] <- ConfMat_Stats_byClass
      
    }
    
    # Return
    return(list("Models"=models,
                "PredStats_byClass"=confStats,
                "PredStatsOverall"=confStatsByClass,
                "Params"=comboGrid))
    
  }
  
}

MLConfMatResultsAgregg<-function(ListOfresults,levels=c(1,2)){
  #to this function we give it a list where the results of the modelling are stored
  #from the previous functions
  
  #careful not tested
  
  
  listOfStats<-list()
  for(i in 1:length(ListOfresults)){
    
    listOfStats[[i]]<-ListOfresults[[i]]$rfmodelRes$PredStats_byClass
  }
  
  AllStatsByClass<-do.call(rbind, listOfStats)
  
  statis<-describeBy(AllStatsByClass[,-1],AllStatsByClass[,1],mat = F)
  
  return(statis)
}




#function to apply s ----

applyModelsToDataSplits <- function(inputListwSplits, varsForModel ,ntree=2000) {
  
  #inputListwSplits assumes an input that will have n elements with each element containing
  # train (train data), test (test data), factor_train, factor_test and idx
  
  #varsForModel is a vector with the variables selected from previous steps
  
  #ntree is for the RF model
  
  result <- list()  # Create an empty list to store the iterations
  
  for (i in 1:length(inputListwSplits)) {
    train<-as.matrix(inputListwSplits[[i]]$train)
    test<-as.matrix(inputListwSplits[[i]]$test)
    factor_train<-inputListwSplits[[i]]$factor_train
    factor_test<-inputListwSplits[[i]]$factor_test
    
    
    rfmodelRes<-RandomForestWithStats(train=train[,colnames(train)%in%varsForModel],test[,colnames(test)%in%varsForModel],factor_train,factor_test,ntree=ntree)
    
    
    result[[i]] <- list("data"=inputListwSplits[[i]], "varSel"=varsForModel,"rfmodelRes"=rfmodelRes)
    
  }
  
  
  
  
  return(result)
}


#bring together variables
concatenateVarSelec<-function(res){
  out<-c()
  for(i in 1:length(res)){
    out<-c(out,res[[i]]$varSel)
  }
  return(out)
}

#performance evaluation
modelMetrExtract<-function(res){
  out<-data.frame()
  vars<-data.frame()
  
  for(i in 1:length(res)){
    
    if(is.null(res[[i]])){
      next
    }else{
      
      out<-rbind(out,res[[i]]$rfmodelRes$PredStats_byClass)
      vars<-rbind(vars,paste(res[[i]]$varSel,collapse = " "))
    }
  }
  colnames(out)<-names(res[[1]]$rfmodelRes$PredStats_byClass)
  finalRes<-data.frame("vars"=vars,out)
  colnames(finalRes)[1]<-"vars"
  
  return(finalRes)
}




fit_LR = function(data_tr, labels_tr, data_ts){
  
  data_tr$labels = labels_tr
  
  model = glm(labels~. ,data = data_tr ,family='binomial')
  preds = predict(model , newdata = data_ts ,type = 'response' )
  preds
  
}