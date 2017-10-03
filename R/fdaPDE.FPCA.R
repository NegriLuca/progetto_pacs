#observation is a matrix with nrow=number of location points (or mesh nodes) and n cols=n of observations

smooth.FEM.FPCA<-function(locations = NULL, observations, FEMbasis, lambda, BC = NULL, GCV = FALSE, CPP_CODE = TRUE, nPC)
{
 covariates=NULL
 if(class(FEMbasis$mesh) == "MESH2D"){
 	ndim = 2
 	mydim = 2
 }else if(class(FEMbasis$mesh) == "MESH.2.5D"){
 	ndim = 3
 	mydim = 2
 }else{
 	stop('Unknown mesh class')
 }
 
##################### Checking parameters, sizes and conversion ################################

  checkSmoothingParameters(locations, observations, FEMbasis, lambda, covariates, BC, GCV, CPP_CODE, PDE_parameters_constant = NULL, PDE_parameters_func = NULL) 
  ## Coverting to format for internal usage
  if(!is.null(locations))
    locations = as.matrix(locations)
  observations = as.matrix(observations)
  lambda = as.matrix(lambda)
  if(!is.null(BC))
  {
    BC$BC_indices = as.matrix(BC$BC_indices)
    BC$BC_values = as.matrix(BC$BC_values)
  }
#MODIFICARE PER SURFACE_MESH	  
  checkSmoothingParametersSize(locations, observations, FEMbasis, lambda, covariates, BC, GCV, CPP_CODE, PDE_parameters_constant = NULL, PDE_parameters_func = NULL, ndim, mydim,FPCA=TRUE)
	  ################## End checking parameters, sizes and conversion #############################

  if(class(FEMbasis$mesh) == 'MESH2D'){	  
  	bigsol = NULL
	print('C++ Code Execution')
	bigsol = CPP_smooth.FEM.FPCA(locations, observations, FEMbasis, lambda,
	ndim, mydim, BC, GCV,nPC)
	numnodes = nrow(FEMbasis$mesh$nodes)
	  
  } else if(class(FEMbasis$mesh) == 'MESH.2.5D'){

	  bigsol = NULL  
	  print('C++ Code Execution')
	  bigsol = CPP_smooth.manifold.FEM.FPCA(locations, observations, FEMbasis$mesh,
	  lambda, ndim, mydim, BC, GCV,nPC)
	  numnodes = FEMbasis$mesh$nnodes
  }
  return(bigsol)
  }
