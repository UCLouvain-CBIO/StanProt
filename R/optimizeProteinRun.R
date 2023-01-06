##' @export
##' 
##' @author Philippe Hauchamps
optimizeProteinRun <- function(nPeptPerProt, nMissingDataPerProt, tasks){
  nProt <- length(nPeptPerProt)
  if(nProt < 1) stop("nPeptPerProt array empty")
  if(length(nMissingDataPerProt) != nProt) stop("Array size mismatch")
  protDifficultyOrder <- order(nPeptPerProt, nMissingDataPerProt, decreasing = TRUE)
  
  integerTaskSize <- nProt %/% tasks
  moduloTaskSize <- nProt %% tasks
  
  taskSizes <- rep(integerTaskSize, tasks)
  if(moduloTaskSize > 0){
    for(i in 1:moduloTaskSize){
      taskSizes[i] <- taskSizes[i] + 1
    }
  }
  
  tasksProts <- list()
  for(i in 1:tasks){
    tasksProts[[i]] <- rep(0., taskSizes[i])
  }
  # nRemainingTasks <- tasks
  # currentTaskIndex <- 1
  # currentProtIndex <- 1
  # while(nRemainingTasks > 0){
  #   allocatedTasks <- min(c(nRemainingTasks, cores))
  #   allocatedProts <- 
  #     sum(taskSizes[currentTaskIndex:(currentTaskIndex + allocatedTasks - 1)])
  #   remainingProts <- allocatedProts
  #   currentProtIndexWithinTask <- rep(1., allocatedTasks)
  #   while(remainingProts > 0){
  #     # do once direct order
  #     for(i in 1:allocatedTasks){
  #       if(currentProtIndexWithinTask[i] <= taskSizes[currentTaskIndex + i -1]){
  #         tasksProts[[currentTaskIndex + i -1]][currentProtIndexWithinTask[i]] <- protDifficultyOrder[currentProtIndex]
  #         currentProtIndex <- currentProtIndex + 1
  #         currentProtIndexWithinTask[i] <- currentProtIndexWithinTask[i] + 1
  #         remainingProts <- remainingProts - 1
  #       }
  #     }
  #     # do once reverse order
  #     for(i in seq(allocatedTasks,1,by = -1)){
  #       if(currentProtIndexWithinTask[i] <= taskSizes[currentTaskIndex + i -1]){
  #         tasksProts[[currentTaskIndex + i - 1]][currentProtIndexWithinTask[i]] <- protDifficultyOrder[currentProtIndex]
  #         currentProtIndex <- currentProtIndex + 1
  #         currentProtIndexWithinTask[i] <- currentProtIndexWithinTask[i] + 1
  #         remainingProts <- remainingProts - 1
  #       }
  #     }
  #     
  #   }
  #   nRemainingTasks <- nRemainingTasks - allocatedTasks
  #   currentTaskIndex <- currentTaskIndex + allocatedTasks
  # }
  
  currentProtIndex <- 1
  remainingProts <- nProt
  currentProtIndexWithinTask <- rep(1., tasks)
  while(remainingProts > 0){
    # do once direct order
    for(i in 1:tasks){
      if(currentProtIndexWithinTask[i] <= taskSizes[i]){
        tasksProts[[i]][currentProtIndexWithinTask[i]] <- protDifficultyOrder[currentProtIndex]
        currentProtIndex <- currentProtIndex + 1
        currentProtIndexWithinTask[i] <- currentProtIndexWithinTask[i] + 1
        remainingProts <- remainingProts - 1
      }
    }
    # do once reverse order
    for(i in seq(tasks,1,by = -1)){
      if(currentProtIndexWithinTask[i] <= taskSizes[i]){
        tasksProts[[i]][currentProtIndexWithinTask[i]] <- protDifficultyOrder[currentProtIndex]
        currentProtIndex <- currentProtIndex + 1
        currentProtIndexWithinTask[i] <- currentProtIndexWithinTask[i] + 1
        remainingProts <- remainingProts - 1
      }
    }
  }
    
  res <- list()
  res$runOrder <- unlist(tasksProts)
  res$runIndex <- order(res$runOrder)
  res$protDifficultyOrder <- protDifficultyOrder
  res
}
