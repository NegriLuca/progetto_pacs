#ifndef __MIXEDFEREGRESSION_FPCA_SPEC_HPP__
#define __MIXEDFEREGRESSION_FPCA_SPEC_HPP__

template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class MixedFERegression<FPCAData,Integrator,ORDER,mydim,ndim>{
	private:
		const MeshHandler<ORDER,mydim,ndim> &mesh_;
		FPCAData& fPCAData_;
		std::vector<coeff> tripletsData_;

		//SpMat NWblock_;
		SpMat psi_;
		SpMat DMat_;
		SpMat AMat_;
		SpMat MMat_;

		MatrixXr Q_;
		MatrixXr H_;

		SpMat _coeffmatrix;       //!A Eigen::VectorXr: Stores the system right hand side.
		VectorXr _b;			  //!A Eigen::VectorXr : Stores the system solution
		std::vector<VectorXr> _solution;
		std::vector<Real> _dof;

		void setPsi();
		void setQ();
		void setH();

		void buildCoeffMatrix(const SpMat& DMat,  const SpMat& AMat,  const SpMat& MMat, SpMat& coeffmatrix);

	public:
		//!A Constructor.
		MixedFERegression(const MeshHandler<ORDER,mydim,ndim>& mesh, FPCAData& fPCAData):mesh_(mesh), fPCAData_(fPCAData){};
		
		//!A Destructor
		//~Model(){};
		
		 //! A normal member taking six arguments. Build the linear system.
		 /*!
		 * This member builds the system matrix and system right hand side, storing them in _coeffmatr and _b
		\param L is a const reference Eigen::SparseMatrix<> : NW block of the system matrix .
		\param opMat is a const reference Eigen::SparseMatrix<> : NE block of the system matrix .
		\param opMat2 is a const reference Eigen::SparseMatrix<> : SW block of the system matrix .
		\param mass is a const reference Eigen::SparseMatrix<> : SE block of the system matrix .
		\the method modifies _coeffmatrix and _b.
		*/
//		void build( SpMat& L,  SpMat& opMat,  SpMat& opMat2,  SpMat& mass,
//						  const VectorXr& rightside, const VectorXr& forcing_term );

		//
		//		|DMat | AMat^T  |
		//		|AMat | MMat	|


		
		void smoothLaplace();
		
		void smoothFPCA();
		
		//void smoothEllipticPDE();
		//void smoothEllipticPDESpaceVarying();

		 //! A template member for the system resolution.
         /*! the template P is the solutor, possible choices are: SpLu, SpQR, SpCholesky,SpConjGrad.
          *  the solution is stored in _solution		
		*/

		template<typename P>
		void solve(UInt output_index, SpMat& coeffmatrix);
		//! A inline member that returns a VectorXr, returns the whole _solution. 
		inline std::vector<VectorXr> const & getSolution() const{return _solution;};
		inline std::vector<Real> const & getDOF() const{return _dof;};
		//Real eval_sol(MeshTria const & mesh,VectorXr Point p);
		//! A member for printing the solution.
		//void printSolution(std::ostream & out) {out<<_solution; out<<"dim"<<_solution.rows()<<"x"<<_solution.cols();};
		//apply dirichlet boundary condition (first attempt)
		//void applyDir(int* bcindex, double* bcvalues, int n_b, UInt order);
		 //! A normal member taking two arguments: Dirichlet Boundary condition
		 /*!
		  * This member applies Dirichlet boundary conditions on the linear system with the penalization method.
		\param bcindex is a const reference to vector<int> : the global indexes of the nodes to which the boundary condition has to be applied.
		\param bcvalues is a const reference to vector<double> : the values of the boundary conditions relative to bcindex.
		\the method modifies _coeffmatrix and _b
		*/		
		void addDirichletBC(const vector<int>& bcindex, const vector<Real>& bcvalues,  SpMat& coeffmatrix);
		void getDataMatrix(SpMat& DMat);
		void getDataMatrixByIndices(SpMat& DMat);
		void getRightHandData(VectorXr& rightHandData);
		void computeDegreesOfFreedom(UInt output_index, SpMat& coeffmatrix);
};

#include "mixedFERegression_FPCA_spec_imp.h"

#endif
