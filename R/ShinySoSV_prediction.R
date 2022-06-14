#' ShinySoSV prediction
#'
#' This function predict the SV calling performance of SV callers and their combinations
#'
#' @param candidate_callers names of callers
#' @param newdata data frame of variables with which to predict
#' @param performance should be any of "sensitivity", "precision" or "F1 score"
#' @param callset individual caller, pairwise union or intersection
#' @return data frame of predicted performance
#' @export
ShinySoSV_prediction <- function(Candidate_callers, newdata, performance, callset){
  model_name1 <- paste0(c("sen", "pre_off", "F1_score")[c("sensitivity", "precision", "F1_score") %in% performance])
  model_name2 <- paste0(c("", "_UnionIntersect", "_UnionIntersect")[c("individual", "union", "intersection") %in% callset])

  model_name <- apply(expand.grid(model_name1, model_name2), 1, paste, collapse="")
  for(i in c(1:length(model_name))){
    #load(paste0("./Shiny-SoSV/data/","gam",model_name[i],"_callers.RData"))
    load(paste0("sysdata.rda"))
  }

  combine_SV_SVcaller <- c()
  for(i in c(1:length(callset))){
    if(callset[i] %in% c("union","intersection")){
      for (i in c(1:length(candidate_callers))){
        combine_SV_SVcaller <- c(combine_SV_SVcaller, paste0(candidate_callers[i], candidate_callers[!(c(1:length(candidate_callers)) %in% i)],c("Union","Intersect")[c("union", "intersection") %in% callset]))
      }
    }else{
      combine_SV_SVcaller <- c(combine_SV_SVcaller, candidate_callers)
    }
  }

  df_prediction <- c()
  for(i in c(1: length(combine_SV_SVcaller))){
    for(j in c(1: length(model_name1))){
      prediction <- predict(eval(parse(text = paste0("gam", model_name1[j],"_", combine_SV_SVcaller[i]))), newdata, type = "response",se.fit = T,unconditional = TRUE)$fit
      df_prediction <- data.frame(cbind(df_prediction, prediction))
      tmp <- paste0(c("sensitivity", "precision", "F1_score")[c("sen", "pre_off", "F1_score") %in% model_name1[j]])
      colnames(df_prediction)[ncol(df_prediction)] <- paste0("fit_",tmp,"_", combine_SV_SVcaller[i])
    }
  }
  df_prediction <- cbind(newdata, df_prediction)
  write.csv(df_prediction,file = "./Shiny-SoSV_prediction.csv",row.names = FALSE)

  return(df_prediction)
}
