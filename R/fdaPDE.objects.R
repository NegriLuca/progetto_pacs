#' Create a FEM basis
#'
#' @param mesh A \code{MESH2D} or \code{SURFACE_MESH}  object representing the domain triangulation. See \link{create.MESH.2D}, \code{\1ink{create.surface.mesh}}.
#' @return A  \code{FEMbasis} object. This contains the \code{mesh}, along with some additional quantities:
#' if \code{class(mesh) == MESH2D}
#' 	\item{\code{order}}{Either "1" or "2". Order of the Finite Element basis.}
#' 	\item{\code{nbasis}}{Scalar. The number of basis.}
#' 	\item{\code{detJ}}{The determinant of the transformation from the nodes of the reference triangle to the nodes of the i-th triangle; this coincides with the double of the area of the i-th triangle.}
#' 	\item{\code{transf}}{A three-dimensional array such that  \code{transf[i,,]} is the 2-by-2 matrix that transforms the nodes of the reference triangle to the nodes of the i-th triangle.}
#' 	\item{\code{metric}}{A three-dimensional array such that \code{metric[i,,]} is the 2-by-2 matrix \cr
#' 	\code{transf[i,,]^{-1}*transf[i,,]^{-T}}. This matrix is used for the computation
#' of the integrals over the elements of the mesh.}
#' if \code{class(mesh) == SURFACE_MESH}
#' 	\item{\code{order}}{Either "1" or "2". Order of the Finite Element basis.}
#' 	\item{\code{nbasis}}{Scalar. The number of basis.}
#' @description Sets up a Finite Element basis. It requires a triangular mesh, a \code{MESH2D} object, as input.
#' The basis' functions are globally continuos surfaces, that are polynomials once restricted to a triangle in the mesh.
#' Linear if (\code{order = 1}) in the input \code{mesh} and quadratic if (\code{order = 2}) in the input \code{mesh}
#' Finite Element are currently implemented.
#' @usage create.FEM.basis(mesh)
#' @seealso \code{\link{create.MESH.2D}}
#' @examples
#' ## Creates a simple triangulated domain with a concavity; this is a MESH2D object
#' mesh<-create.MESH.2D(nodes=rbind(c(0, 0), c(0, 1), c(0.5, 0.5), c(1, 1), c(1, 0)),
#' segments=rbind(c(1, 2), c(2, 3), c(3, 4), c(4, 5), c(5, 1)), order=1)
#' ## Plot it
#' plot(mesh)
#' ## Creates the basis
#' FEMbasis = create.FEM.basis(mesh)

create.FEM.basis = function(mesh)
{
  if (class(mesh)=="MESH2D"){

	  #  The number of basis functions corresponds to the number of vertices
	  #  for order = 1, and to vertices plus edge midpoints for order = 2

	  nbasis = dim(mesh$nodes)[[1]]
	  eleProp = R_elementProperties(mesh)

	  #eleProp = NULL
	  #if(CPP_CODE == FALSE)
	  #{
	  #  eleProp = R_elementProperties(mesh)
	  #}

	  FEMbasis = list(mesh = mesh, order = as.integer(mesh$order), nbasis = nbasis, detJ=eleProp$detJ, transf = eleProp$transf, metric = eleProp$metric)
	  class(FEMbasis) = "FEMbasis"

	  FEMbasis
  } else if (class(mesh) == "MESH.2.5D"){

  	  FEMbasis = list(mesh = mesh, order = as.integer(mesh$order),nbasis = mesh$nnodes)
  	  class(FEMbasis) = "FEMbasis"
  	  FEMbasis
  	}
}


#' Define a surface or spatial field by a Finite Element basis expansion
#'
#' @param coeff A vector or a matrix containing the coefficients for the Finite Element basis expansion. The number of rows
#' (or the vector's length) corresponds to the number of basis in \code{FEMbasis}.
#' The number of columns corresponds to the number of functional replicates.
#' @param FEMbasis A \code{FEMbasis} object defining the Finite Element basis, created by \link{create.FEM.basis}.
#' @description This function defines a FEM object. This is not usualled called directly by users.
#' @usage FEM(coeff,FEMbasis)
#' @return An \code{FEM} object. This contains a list with components \code{coeff} and \code{FEMbasis}.
#' @examples
#' ## Upload a triangular mesh and plot it
#' data("mesh.2D.rectangular")
#' plot(mesh.2D.rectangular)
#' ## Create a linear Finite Element basis
#' FEMbasis = create.FEM.basis(mesh.2D.rectangular)
#' ## Define a sinusoidal function as expansion of this basis and plot it
#' coeff <- sin(mesh.2D.rectangular$nodes[,1])*cos(mesh.2D.rectangular$nodes[,2])
#' FEM_object<- FEM(coeff, FEMbasis)
#' plot(FEM_object)

FEM<-function(coeff,FEMbasis)
{
  if (is.null(coeff))
    stop("coeff required;  is NULL.")
  if (is.null(FEMbasis))
    stop("FEMbasis required;  is NULL.")
  if(class(FEMbasis) != "FEMbasis")
    stop("FEMbasis not of class 'FEMbasis'")
  coeff = as.matrix(coeff)
  if(nrow(coeff) != FEMbasis$nbasis)
    stop("Number of row of 'coeff' different from number of basis")

  fclass = NULL
  fclass = list(coeff=coeff, FEMbasis=FEMbasis)
  class(fclass)<-"FEM"
  return(fclass)
}

#' Plot a \code{FEM} object
#'
#' @param x A \code{FEM} object.
#' @param num_refinements A natural number specifying how many bisections should by applied to each triangular element for
#' plotting purposes. This functionality is useful where a discretization with 2nd order Finite Element is applied.
#' @param ... Arguments representing graphical options to be passed to \link[rgl]{plot3d}.
#' @description Three-dimensional plot of a \code{FEM} object, generated by \code{FEM} or returned by \code{smooth.FEM.basis}, \code{smooth.FEM.PDE.basis} or
#' \code{smooth.FEM.PDE.sv.basis}.
#' @usage \method{plot}{FEM}(x, num_refinements, ...)
#' @seealso \code{\link{image.FEM}}
#' @examples
#' ## Upload a triangular mesh and plot it
#' data("mesh.2D.rectangular")
#' plot(mesh.2D.rectangular)
#' ## Create a linear Finite Element basis
#' FEMbasis = create.FEM.basis(mesh.2D.rectangular)
#' ## Define a sinusoidal function as expansion of this basis and plot it
#' coeff <- sin(mesh.2D.rectangular$nodes[,1])*cos(mesh.2D.rectangular$nodes[,2])
#' FEM_object<- FEM(coeff, FEMbasis)
#' plot(FEM_object)

plot.FEM = function(x, num_refinements = NULL, ...)
{
  if(x$FEMbasis$order == 1)
  {
    R_plot.ORD1.FEM(x, ...)
  }else{
    R_plot.ORDN.FEM(x, num_refinements, ...)
  }
}

#' Image Plot of a FEM object
#'
#' @param x A \code{FEM} object.
#' @param num_refinements A natural number specifying how many bisections should by applied to each triangular element for
#' plotting purposes. This functionality is useful where a discretization with 2nd order Finite Element is applied.
#' @param ... Arguments representing  graphical options to be passed to \link[rgl]{plot3d}.
#' @description Image plot of a \code{FEM} object, generated by the function \code{FEM} or returned by \code{smooth.FEM.basis}, \code{smooth.FEM.PDE.basis} or
#' \code{smooth.FEM.PDE.sv.basis} can be visualized through an image plot.
#' @usage \method{image}{FEM}(x, num_refinements, ...)
#' @seealso \code{\link{plot.FEM}}
#' @examples
#' ## Upload a triangular mesh and plot it
#' data("mesh.2D.rectangular")
#' plot(mesh.2D.rectangular)
#' ## Create a linear Finite Element basis
#' FEMbasis = create.FEM.basis(mesh.2D.rectangular)
#' ## Define a sinusoidal function as expansion of this basis and plot it
#' coeff <- sin(mesh.2D.rectangular$nodes[,1])*cos(mesh.2D.rectangular$nodes[,2])
#' FEM_object<- FEM(coeff, FEMbasis)
#' image(FEM_object)
image.FEM = function(x, num_refinements = NULL, ...)
{
  if(x$FEMbasis$order == 1)
  {
    R_image.ORD1.FEM(x, ...)
  }else{
    R_image.ORDN.FEM(x, num_refinements, ...)
  }
}

#' Create a \code{MESH2.5D} object from the connectivty matrix and nodes locations
#'
#' @param nodes A nnodes x 3 matrix specifying the locations of each node
#' @param triangles A ntriangles x 3*order matrix specifying the indices of the nodes in each triangle
#' @param order{Either "1" or "2". Order of the Finite Element basis default is order = 1
#' @return A \code{SURFACE_MESH} object
#' @examples
#' #read the matrix nodes and triangles from file
#' nodes = read.table(file="mynodes.csv",header=F,sep=",")
#' triangles = read.table(file="mytriangles.csv",header=F,sep=",")
#' mesh = create.MESH.2.5D(nodes,triangles)

create.MESH.2.5D<- function(nodes, triangles, order = 1)
{
  nnodes = dim(nodes)[1]

  ntriangles = dim(triangles)[1]

  if(dim(triangles)[2]!= 3*order){
    if (order==1)
      stop("The matrix 'triangles' has the wrong number of columns. See second.order.mesh(...)")
  	stop("The matrix 'triangles' has wrong number of columns. Should be 3*order \n")
  	}
  out = list(nnodes=nnodes, ntriangles=ntriangles, nodes=c(t(nodes)), triangles = c(t(triangles)), order=as.integer(order))

  class(out)<-"MESH.2.5D"

  return(out)
}

#' Double the order of a fist order Finite Element mesh by adding middle points to each side of the triangles in the triangulation
#' @param mesh an object of class 'SURFACE_MESH' is the starting mesh of order 1
#' @param bc A vector specifying the indices of the nodes on which boundary conditions are applied
#' @return if no boundaries conditions are passed
#' \item{\code{mesh}} An object of class 'SURFACE_MESH' with the mesh of order 2. Otherwhise a \code{\list} with parameters:
#' \item{\code{mesh}} An object of class 'SURFACE_MESH' with the mesh of order 2.
#' \item{\code{bc_index}} An update of the vector specifying the indices of the nodes on which boundary conditions are applied
#' @examples
#' data(hub) #loading mesh hub, order=1
#' hub_order2 = second.order.MESH2.5D(hub)

second.order.MESH2.5D<-function(mesh,bc=NULL){
  if(class(mesh) != 'MESH.2.5D'){
    stop('This method is implemented only for a mesh of class SURFACE_MESH')
  }else if(mesh$order != 1){
    stop('The object mesh must have order = 1')
  }else{
    toll=1e-5
    T = matrix(mesh$triangles,nrow=mesh$ntriangles,ncol=3, byrow = TRUE)
    V = matrix(mesh$nodes, nrow = mesh$nnodes, ncol= 3, byrow = TRUE)
    T <- cbind(T, matrix(0,nrow=nrow(T),ncol=3))
    nnodes=nrow(V)
    index=nrow(V)
    points = V[T[1,],]
    midpoints<-rbind((points[2,]+points[3,])/2,(points[1,]+points[3,])/2, (points[1,]+points[2,])/2);
    if(!is.null(bc)){
      isBC<-c( any(bc==T[1,2]) & any(bc==T[1,3]),
               any(bc==T[1,1]) & any(bc==T[1,3]),
               any(bc==T[1,2]) & any(bc==T[1,1]))
    }

    for (side in 1:3){
      point<-midpoints[side,]
      index<-index+1;
      V<-rbind(V,point)
      T[1,3+side]<-index;

      if(!is.null(bc)&&isBC[side]==1){
        bc<-c(bc,index)
      }

    }

    for (i in 2:nrow(T)){
      points = V[T[i,],]
      midpoints<-rbind((points[2,]+points[3,])/2,(points[1,]+points[3,])/2, (points[1,]+points[2,])/2);
      if(!is.null(bc)){
        isBC<-c( any(bc==T[i,2]) & any(bc==T[i,3]),
                 any(bc==T[i,1]) & any(bc==T[i,3]),
                 any(bc==T[i,2]) & any(bc==T[i,1]))
      }

      for (side in 1:3){
        point<-midpoints[side,]
        isthere<-apply(V[(nnodes+1):nrow(V),], 1, function(x) identical(as.vector(x), point))
        loc = which(isthere)
        if(length(loc)>0){
          loc = loc+nnodes
          T[i,3+side]<-loc[1]
        }else{
          index<-index+1;
          V<-rbind(V,point)
          T[i,3+side]<-index;

          if(!is.null(bc)&&isBC[side]==1){
            bc<-c(bc,index)
          }
        }
      }
    }
  }
  if(is.null(bc)){
    out = list(nnodes=nrow(V), ntriangles=nrow(T), nodes=c(t(V)), triangles = c(t(T)), order=2)
    class(out)<-"MESH.2.5D"
    return(out)
  }else{
    out = list(nnodes=nrow(V), ntriangles=nrow(T), nodes=c(t(V)), triangles = c(t(T)), order=2)
    class(out)<-"MESH.2.5D"
    retlist = list(mesh = out, bc_index=bc)
    return(retlist)
  }
}
