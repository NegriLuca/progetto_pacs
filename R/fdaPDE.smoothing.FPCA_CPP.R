CPP_smooth.FEM.basis<-function(locations, observations, FEMbasis, lambda, ndim, mydim, BC = NULL, GCV,nPC)
{ 
  # Indexes in C++ starts from 0, in R from 1, opportune transformation
  ##TO BE CHANGED SOON: LOW PERFORMANCES, IMPLIES COPY OF PARAMETERS
  FEMbasis$mesh$triangles = FEMbasis$mesh$triangles - 1
  FEMbasis$mesh$edges = FEMbasis$mesh$edges - 1
  FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] = FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] - 1
  
  if(is.null(locations))
  {
    locations<-matrix(nrow = 0, ncol = ndim)
  }
  
  if(is.null(BC$BC_indices))
  { 
    BC$BC_indices<-vector(length=0)
  }else
  { #the conversion of these indices is done in c++ now
    BC$BC_indices<-as.vector(BC$BC_indices)
  }
  
  if(is.null(BC$BC_values))
  {
    BC$BC_values<-vector(length=0)
  }else
  {
    BC$BC_values<-as.vector(BC$BC_values)
  }
  
  ## Set propr type for correct C++ reading
  locations <- as.matrix(locations)
  storage.mode(locations) <- "double"
  storage.mode(FEMbasis$mesh$points) <- "double"
  storage.mode(FEMbasis$mesh$triangles) <- "integer"
  storage.mode(FEMbasis$mesh$edges) <- "integer"
  storage.mode(FEMbasis$mesh$neighbors) <- "integer"
  storage.mode(FEMbasis$order) <- "integer"
  storage.mode(lambda)<- "double"
  storage.mode(ndim)<- "integer"
  storage.mode(mydim)<- "integer"
  nPC<-as.integer(nPC)
  storage.mode(nPC)<-"integer"
  storage.mode(BC$BC_indices)<- "integer"
  storage.mode(BC$BC_values)<-"double"
  
  GCV = as.integer(GCV)
  storage.mode(GCV)<-"integer"
  
  ## Call C++ function
  bigsol <- .Call("FPCA_Laplace", locations, observations, FEMbasis$mesh, 
                  FEMbasis$order,mydim,ndim, lambda,
                  BC$BC_indices, BC$BC_values, GCV,nPC,
                  package = "fdaPDE")
  
  ## Reset them correctly
  #fdobj$basis$params$mesh$triangles = fdobj$basis$params$mesh$triangles + 1
  #fdobj$basis$params$mesh$edges = fdobj$basis$params$mesh$edges + 1
  #fdobj$basis$params$mesh$neighbors = fdobj$basis$params$mesh$neighbors + 1
  return(bigsol)
}



CPP_smooth.manifold.FEM.basis<-function(locations, observations, mesh, lambda, ndim, mydim, BC = NULL, GCV,nPC)
{
  # Indexes in C++ starts from 0, in R from 1, opportune transformation
  # This is done in C++ now to optimize speed

  if(is.null(locations))
  {
    locations<-matrix(nrow = 0, ncol = ndim)
  }
  
  if(is.null(BC$BC_indices))
  {
    BC$BC_indices<-vector(length=0)
  }else
  { 
    BC$BC_indices<-as.vector(BC$BC_indices)
  
  }
  
  if(is.null(BC$BC_values))
  {
    BC$BC_values<-vector(length=0)
  }else
  {
    BC$BC_values<-as.vector(BC$BC_values)
  }
  
  ## Set propr type for correct C++ reading
  locations <- as.matrix(locations)
  storage.mode(locations) <- "double"
  observations <- as.matrix(observations)
  storage.mode(observations) <- "double"
  storage.mode(mesh$order) <- "integer"
  storage.mode(mesh$nnodes) <- "integer"
  storage.mode(mesh$ntriangles) <- "integer"
  storage.mode(mesh$nodes) <- "double"
  storage.mode(mesh$triangles) <- "integer"
  storage.mode(lambda) <- "double"
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"
  nPC<-as.integer(nPC)
  storage.mode(nPC)<-"integer"
  storage.mode(BC$BC_indices) <- "integer"
  storage.mode(BC$BC_values)  <- "double"
  GCV = as.integer(GCV)
  storage.mode(GCV)<-"integer"
  
  ## Call C++ function
  bigsol <- .Call("FPCA_Laplace", locations, observations, mesh, 
                  mesh$order, mydim, ndim, lambda,
                  BC$BC_indices, BC$BC_values, GCV, nPC,
                  package = "fdaPDE")

  return(bigsol)
}

