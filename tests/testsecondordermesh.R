source('../R/datasets.R',chdir=T)
source('../R/fdaPDE.objects.R')
source("../R/fdaPDE.smoothing.R")
source("../R/fdaPDE.checkParameters.R")
source("../R/fdaPDE.smoothing_CPP.R")
source("../R/fdaPDE.smoothing_R.R")
source("../R/fdaPDE.locator.R")
source("../R/mesh.2D.R")
source("../R/zzz.R")
source("../R/fdaPDE.smoothing.manifold_CPP.R")
source("../R/fdaPDE.plot.mesh_R.R")
dyn.load("../src/fdaPDE.so") 


V = read.table(file="V_caramella2.csv",header=F,sep=",")
T = read.table(file="T_caramella2.csv",header=F,sep=",")
mymesh=list(nodes=V,triangles=T,nnodes=nrow(V),ntriangles=nrow(T))

times=50
MSE_caramella=numeric(times)

for(j in 1:times){
  ########## Generating Test data ############
  
  a1 = rnorm(1,mean = 1, sd = 1)
  a2 = rnorm(1,mean = 1, sd = 1)
  a3 = rnorm(1,mean = 1, sd = 1)
  
  nnodes=mymesh$nnodes
  func_evaluation = numeric(nnodes)
  
  for (i in 1:nnodes){
    func_evaluation[i] = a1* sin(2*pi*mymesh$nodes[i,1]) +  a2* sin(2*pi*mymesh$nodes[i,2]) +  a3*sin(2*pi*mymesh$nodes[i,3]) +1
  }
  
  data=func_evaluation+rnorm(nnodes,mean=0,sd=0.5)
  
  
  mesh <- create.surface.mesh(mymesh$nodes, mymesh$triangles, order=2)
  
  FEMbasis <- create.FEM.basis(mesh)
  
  lambda=c(0.00375)
  output_CPP =smooth.FEM.basis(observations = data, 
                               FEMbasis = FEMbasis, lambda = lambda,
                               CPP_CODE = TRUE) 
  
  MSE_caramella[j] = sum((as.vector(output_CPP$fit.FEM$coeff)-func_evaluation)^2)/nnodes
}

write.table(MSE_caramella,file="MSE_caramella2order.csv",sep=",")
