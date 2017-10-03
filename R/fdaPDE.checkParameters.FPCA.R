checkSmoothingParametersSizeFPCA<-function(locations = NULL, observations, FEMbasis, lambda, covariates = NULL, BC = NULL, GCV = FALSE, CPP_CODE = TRUE, PDE_parameters_constant = NULL, PDE_parameters_func = NULL, ndim, mydim)
{
  #################### Parameter Check #########################
  if(nrow(observations) < 1)
    stop("'observations' must contain at least one element")
  if(is.null(locations))
  {
    if(class(FEMbasis$mesh) == "MESH2D"){
    	if(nrow(observations) > nrow(FEMbasis$mesh$nodes))
     	 stop("Size of 'observations' is larger then the size of 'nodes' in the mesh")
    }else if(class(FEMbasis$mesh) == "MESH2.5D"){
    	if(nrow(observations) > FEMbasis$mesh$nnodes)
     	 stop("Size of 'observations' is larger then the size of 'nodes' in the mesh")
    }
  }
  if(!is.null(locations))
  {
    if(ncol(locations) != ndim)
      stop("'locations' and the mesh points have incompatible size;")
    if(nrow(locations) != nrow(observations))
      stop("'locations' and 'observations' have incompatible size;")
  }
  if(ncol(lambda) != 1)
    stop("'lambda' must be a column vector")
  if(nrow(lambda) < 1)
    stop("'lambda' must contain at least one element")
  if(!is.null(covariates))
  {
    if(nrow(covariates) != nrow(observations))
      stop("'covariates' and 'observations' have incompatible size;")
  }
  
  if(!is.null(BC))
  {
    if(ncol(BC$BC_indices) != 1)
      stop("'BC_indices' must be a column vector")
    if(ncol(BC$BC_values) != 1)
      stop("'BC_values' must be a column vector")
    if(nrow(BC$BC_indices) != nrow(BC$BC_values))
      stop("'BC_indices' and 'BC_values' have incompatible size;")
    if(class(FEMbasis$mesh) == "MESH2D"){
	    if(sum(BC$BC_indices>nrow(nrow(FEMbasis$mesh$nodes))) > 0)
	      stop("At least one index in 'BC_indices' larger then the number of 'nodes' in the mesh;")
    }else if((class(FEMbasis$mesh) == "MESH2D")){
    	if(sum(BC$BC_indices>FEMbasis$mesh$nnodes) > 0)
	      stop("At least one index in 'BC_indices' larger then the number of 'nodes' in the mesh;")
    	}
  }
  
  if(!is.null(PDE_parameters_constant))
  {
    if(!all.equal(dim(PDE_parameters_constant$K), c(2,2)))
      stop("'K' in 'PDE_parameters must be a 2x2 matrix")
    if(!all.equal(dim(PDE_parameters_constant$b), c(2,1)))
      stop("'b' in 'PDE_parameters must be a column vector of size 2")
    if(!all.equal(dim(PDE_parameters_constant$c), c(1,1)))
      stop("'c' in 'PDE_parameters must be a double")
  }
  
  if(!is.null(PDE_parameters_func))
  {
    
    n_test_points = min(nrow(FEMbasis$mesh$nodes), 5)
    test_points = FEMbasis$mesh$nodes[1:n_test_points, ]
    
    try_K_func = PDE_parameters_func$K(test_points)
    try_b_func = PDE_parameters_func$b(test_points)
    try_c_func = PDE_parameters_func$c(test_points)
    try_u_func = PDE_parameters_func$u(test_points)
    
    if(!is.numeric(try_K_func))
      stop("Test on function 'K' in 'PDE_parameters' not passed; output is not numeric")
    if(!all.equal(dim(try_K_func), c(2,2,n_test_points)) )
      stop("Test on function 'K' in 'PDE_parameters' not passed; wrong size of the output")
    
    if(!is.numeric(try_b_func))
      stop("Test on function 'b' in 'PDE_parameters' not passed; output is not numeric")
    if(!all.equal(dim(try_b_func), c(2,n_test_points)))
      stop("Test on function 'b' in 'PDE_parameters' not passed; wrong size of the output")
    
    if(!is.numeric(try_c_func))
      stop("Test on function 'c' in 'PDE_parameters' not passed; output is not numeric")
    if(length(try_c_func) != n_test_points)
      stop("Test on function 'c' in 'PDE_parameters' not passed; wrong size of the output")
    
    if(!is.numeric(try_u_func))
      stop("Test on function 'u' in 'PDE_parameters' not passed; output is not numeric")
    if(length(try_u_func) != n_test_points)
      stop("Test on function 'u' in 'PDE_parameters' not passed; wrong size of the output")
  }
  
}
