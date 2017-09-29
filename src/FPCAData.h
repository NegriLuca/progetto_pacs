#ifndef __FPCADATA_HPP__
#define __FPCADATA_HPP__

#include "fdaPDE.h"
#include "mesh_objects.h"
#include "param_functors.h"

//!  An IO handler class for objects passed from R
/*!
 * This class, given the data from R, convert them in a C++ format, offering a
 * series of method for their access, so isolating the more the possible the specific
 * code for R/C++ data conversion.
*/
class  FPCAData{
	private:

		// Design matrix pointer and dimensions
		std::vector<Point> locations_;
		
		
		//Design matrix
		MatrixXr observations_;
		UInt n_;
		UInt p_;
		std::vector<UInt> observations_indices_;
		bool locations_by_nodes_;

		//Other parameters
		UInt order_;
		std::vector<Real> lambda_;

		std::vector<Real> bc_values_;
		std::vector<UInt> bc_indices_;
		
		//Number of Principal Components to compute
		UInt nPC_;
		

		//bool inputType;
		bool DOF_;

	public:
		//! A basic version of the constructor.

		/*!
			It initializes the object storing the R given objects. This is the simplest of the two possible interfaces with R
			\param Robservations an R-vector containing the values of the observations.
			\param Rdesmat an R-matrix containing the design matrix for the regression.
			\param Rorder an R-integer containing the order of the approximating basis.
			\param Rlambda an R-double containing the penalization term of the empirical evidence respect to the prior one.
			\param Rbindex an R-integer containing the indexes of the nodes the user want to apply a Dirichlet Condition,
					the other are automatically considered in Neumann Condition.
			\param Rbvalues an R-double containing the value to impose for the Dirichlet condition, on the indexes specified in Rbindex
		*/


		//! A complete version of the constructor.
		/*!
			It initializes the object storing the R given objects. This is the simplest of the two possible interfaces with R
			\param Robservations an R-vector containing the values of the observations.
			\param Rdesmat an R-matrix containing the design matrix for the regression.
			\param Rorder an R-integer containing the order of the approximating basis.
			\param Rlambda an R-double containing the penalization term of the empirical evidence respect to the prior one.
			\param Rbindex an R-integer vector containing the indexes of the nodes the user want to apply a Dirichlet Condition,
					the other are automatically considered in Neumann Condition.
			\param Rbvalues an R-double vector containing the value to impose for the Dirichlet condition, on the indexes specified in Rbindex
			\param Rc an R-double that contains the coefficient of the REACTION term
			\param Rbeta an R-double 2-dim vector that contains the coefficients for the TRANSPORT coefficients.
			\param RK an R-double 2X2 matrix containing the coefficients for a anisotropic DIFFUSION term.
			\param (UNSUPPORTED put it zero) Ru an R-double vector of length #triangles that contaiins the forcing term integrals.
		*/
		FPCAData(){};

		explicit FPCAData(std::vector<Point>& locations, MatrixXr& observations, UInt order, std::vector<Real> lambda, std::vector<UInt>& bc_indices, std::vector<Real>& bc_values,UInt nPC, bool DOF);


		void printObservations(std::ostream & out) const;
		void printLocations(std::ostream & out) const;

		//! A method returning a reference to the observations vector
		inline MatrixXr const & getObservations() const {return observations_;}
		//! A method returning the number of location points
		inline UInt const getNumberofLocations() const {return observations_.rows();}
		//! A method returning the number of observations
		inline UInt const getNumberofObservations() const {return observations_.cols();}
		//! A method returning the locations of the observations
		inline std::vector<Point> const & getLocations() const {return locations_;}
		inline bool isLocationsByNodes() const {return locations_by_nodes_;}
		inline bool computeDOF() const {return DOF_;}
		inline std::vector<UInt> const & getObservationsIndices() const {return observations_indices_;}
		//! A method returning the the penalization term
		inline std::vector<Real> const & getLambda() const {return lambda_;}
		//! A method returning the input order
		inline UInt const getOrder() const {return order_;}
		//! A method returning the indexes of the nodes for which is needed to apply Dirichlet Conditions
		inline std::vector<UInt> const & getDirichletIndices() const {return bc_indices_;}
		//! A method returning the values to apply for Dirichlet Conditions
		inline std::vector<Real> const & getDirichletValues() const {return bc_values_;}
		//! A method returning the number of Principal Components to compute
		inline UInt const getNPC() const {return nPC_;}
};

#include "FPCAData_imp.h"

#endif
