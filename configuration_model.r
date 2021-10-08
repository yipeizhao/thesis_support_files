#include necessary libraries, and setting up the environment
library("reticulate")
library("igraph")
library("plotly")
#fill the file path with your python path
use_python("/Home/anaconda3/bin")

#function1 return the edges list to build the network in networkx
#function2 return the graph itself for tespecing in the R environment
configuration_model <- function(series){
  g <- sample_degseq(series, method="vl")
  return(as_edgelist(g))
}
return_graph<- function(series){
  g <- sample_degseq(series, method="vl")
  return(g)
}

#python functions
#you need to install all packages in the following python package to run the script
#how to install python packages in R:
#https://rstudio.github.io/reticulate/articles/python_packages.html
source_python("R_functions.py")

#Building the network in python
n = 50
gamma = 3
gamma_list = seq(gamma-0.3,gamma+0.3,0.01)
mari_results = integer(length(gamma_list))
OdC_results = integer(length(gamma_list))
C1espec_results = integer(length(gamma_list))
graphs = list()
for(niter in 1:50){
for(i in 1:length(gamma_list)){
  series = power_law_series(n,gamma_list[i])
  if(sum(series)%%2 !=0){
    series[1]=series[1]+1
  }
  edge_list = configuration_model(series)
  g = construct_network(edge_list)
  graphs = append(graphs,g)
  mari_results[i] = mari_results[i]+MAri(g)
  OdC_results[i] = OdC_results[i]+OdC(g)
  C1espec_results[i] = C1espec_results[i]+C1espec(g)
}
  }


mari_results = mari_results/niter
OdC_results = OdC_results/niter
C1espec_results = C1espec_results/niter

dev.new(3,5)
plot(gamma_list, C1espec_results, col="black",pch="*","l", lty=1, ylim=c(0,0.5), ylab="Complexity",xlab="Gamma" )
lines(gamma_list, OdC_results, col="blue",lty=2)

lines(gamma_list, mari_results, col="magenta", lty=4)

legend(3.1,0.45,legend=c(expression(paste('C'["1e,spec"])),"OdC","MAri"),
       col=c("black","blue","red"),
       lty=c(1,2,3), ncol=1)
