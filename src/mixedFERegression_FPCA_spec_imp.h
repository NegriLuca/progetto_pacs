#ifndef __MIXEDFEREGRESSION_FPCA_SPEC_IMP_HPP__
#define __MIXEDFEREGRESSION_FPCA_SPEC_IMP_HPP__

#include <iostream>


template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegression<FPCAData,Integrator,ORDER, mydim, ndim>::buildCoeffMatrix(const SpMat& DMat,  const SpMat& AMat,  const SpMat& MMat, SpMat& coeffmatrix)
{
	//I reserve the exact memory for the nonzero entries of each row of the coeffmatrix for boosting performance
	//_coeffmatrix.setFromTriplets(tripletA.begin(),tripletA.end());

	UInt nnodes = mesh_.num_nodes();

	std::vector<coeff> tripletAll;
	tripletAll.reserve(DMat.nonZeros() + 2*AMat.nonZeros() + MMat.nonZeros());

	for (int k=0; k<DMat.outerSize(); ++k)
	  for (SpMat::InnerIterator it(DMat,k); it; ++it)
	  {
		  tripletAll.push_back(coeff(it.row(), it.col(),it.value()));
	  }
	for (int k=0; k<MMat.outerSize(); ++k)
	  for (SpMat::InnerIterator it(MMat,k); it; ++it)
	  {
		  tripletAll.push_back(coeff(it.row()+nnodes, it.col()+nnodes,it.value()));
	  }
	for (int k=0; k<AMat.outerSize(); ++k)
	  for (SpMat::InnerIterator it(AMat,k); it; ++it)
	  {
		  tripletAll.push_back(coeff(it.col(), it.row()+nnodes,it.value()));
	  }
	for (int k=0; k<AMat.outerSize(); ++k)
	  for (SpMat::InnerIterator it(AMat,k); it; ++it)
	  {
		  tripletAll.push_back(coeff(it.row()+nnodes, it.col(), it.value()));
	  }

	coeffmatrix.setZero();
	coeffmatrix.resize(2*nnodes,2*nnodes);
	coeffmatrix.setFromTriplets(tripletAll.begin(),tripletAll.end());
	coeffmatrix.makeCompressed();

}

template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegression<FPCAData,Integrator,ORDER, mydim, ndim>::addDirichletBC(const vector<int>& bcindex, const vector<Real>& bcvalues, SpMat& coeffmatrix)
{
	UInt id1,id3;

	UInt nnodes = mesh_.num_nodes();

	const std::vector<UInt>& bc_indices = fPCAData_.getDirichletIndices();
	const std::vector<Real>& bc_values = fPCAData_.getDirichletValues();
	UInt nbc_indices = bc_indices.size();

	Real pen=10e20;

	for( auto i=0; i<nbc_indices; i++)
	 {
			id1=bcindex[i];
			id3=id1+nnodes;

			//_coeffmatrix.prune([id1,id3](int i, int j, Real) { return (i!=id1 && i!=id3); });

			coeffmatrix.coeffRef(id1,id1)=pen;
			coeffmatrix.coeffRef(id3,id3)=pen;


			_b(id1)+=bc_values[i]*pen;
			_b(id3)=0;
	 }

	coeffmatrix.makeCompressed();

}

//construct NW block of the system matrix in Ramsay with covariates format
//!! Depends on setPsi and setQ
template<typename Integrator, UInt ORDER,  UInt mydim, UInt ndim>
void MixedFERegression<FPCAData,Integrator,ORDER,mydim,ndim>::getDataMatrix(SpMat& DMat)
{
		UInt nnodes = mesh_.num_nodes();
		//UInt nlocations = regressionData_.getNumberofObservations();

		DMat.resize(nnodes,nnodes);
		DMat = psi_.transpose()*psi_;

}

template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegression<FPCAData,Integrator,ORDER,mydim,ndim>::getDataMatrixByIndices(SpMat& DMat)
{
		UInt nnodes = mesh_.num_nodes();
		UInt nlocations = fPCAData_.getNumberofObservations();

		DMat.resize(nnodes,nnodes);

			DMat.reserve(1);
			for (auto i = 0; i<nlocations; ++i)
			{
				auto index = fPCAData_.getObservationsIndices()[i];
				DMat.insert(index,index) = 1;
			}
}

template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegression<FPCAData,Integrator,ORDER,mydim,ndim>::setPsi(){

		//std::cout<<"Data Matrix Computation by Basis Evaluation.."<<std::endl;
		UInt nnodes = mesh_.num_nodes();
		UInt nlocations = fPCAData_.getNumberofObservations();
		Real eps = 2.2204e-016,
			 tolerance = 100 * eps;

		//cout<<"Nodes number "<<nnodes<<"Locations number "<<nlocations<<endl;

		//std::vector<coeff> entries;
		//entries.resize((ORDER * 3)*nlocations);


		psi_.resize(nlocations, nnodes);
		//psi_.reserve(Eigen::VectorXi::Constant(nlocations,ORDER*3));

		Triangle<ORDER*3,mydim,ndim> tri_activated;
		Eigen::Matrix<Real,ORDER * 3,1> coefficients;

		Real evaluator;

		for(UInt i=0; i<nlocations;i++)
		{
			tri_activated = mesh_.findLocationNaive(fPCAData_.getLocations()[i]);
			if(tri_activated.getId() == Identifier::NVAL)
			{
				#ifdef R_VERSION_
				Rprintf("ERROR: Point %d is not in the domain, remove point and re-perform smoothing\n", i+1);
				#else
				std::cout << "ERROR: Point " << i+1 <<" is not in the domain\n";
				#endif
			}else
			{
				for(UInt node = 0; node < ORDER*3 ; ++node)
				{
					coefficients = Eigen::Matrix<Real,ORDER * 3,1>::Zero();
					coefficients(node) = 1; //Activates only current base
					evaluator = evaluate_point<ORDER,mydim,ndim>(tri_activated, fPCAData_.getLocations()[i], coefficients);
					psi_.insert(i, tri_activated[node].getId()) = evaluator;
				}
			}
		}

		psi_.prune(tolerance);
		psi_.makeCompressed();

}

template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegression<FPCAData,Integrator,ORDER,mydim,ndim>::setQ()
{
	//std::cout<<"Computing Orthogonal Space Projection Matrix"<<std::endl;
	Q_.resize(H_.rows(),H_.cols());
	Q_ = -H_;
	for (int i=0; i<H_.rows();++i)
	{
		Q_(i,i) += 1;
	}

}

template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegression<FPCAData,Integrator,ORDER,mydim,ndim>::setH()
{
	//std::cout<<"Computing Projection Matrix"<<std::endl;
	//UInt nnodes = mesh_.num_nodes();
	UInt nlocations = fPCAData_.getNumberofObservations();

	//regressionData_.printCovariates(std::cout);
	MatrixXr W;
	//std::cout<<"W "<< W <<std::endl;
	//total number of mesh nodes
	//UInt nnodes = mesh_.num_nodes();
	if(fPCAData_.isLocationsByNodes())
	{
		MatrixXr W_reduced(fPCAData_.getNumberofObservations(), W.cols());
		for (auto i=0; i<nlocations;++i)
		{
			auto index_i = fPCAData_.getObservationsIndices()[i];
			for (auto j=0; j<W.cols();++j)
			{
				W_reduced(i,j) = W(index_i,j);
			}
		}
		W = W_reduced;
	}


	MatrixXr WTW(W.transpose()*W);

	H_=W*WTW.ldlt().solve(W.transpose()); // using cholesky LDLT decomposition for computing hat matrix

}

template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegression<FPCAData,Integrator,ORDER,mydim,ndim>::getRightHandData(VectorXr& rightHandData)
{
	UInt nnodes = mesh_.num_nodes();
	UInt nlocations = fPCAData_.getNumberofObservations();
	//rightHandData.resize(nnodes);
	rightHandData = VectorXr::Zero(nnodes);

		if(fPCAData_.isLocationsByNodes())
		{

			for (auto i=0; i<nlocations;++i)
			{
				auto index_i = fPCAData_.getObservationsIndices()[i];
				rightHandData(index_i) = fPCAData_.getDataForRegression()[i];
			}
		}
		else
		{
			rightHandData=psi_.transpose()*fPCAData_.getObservations();
		}
}


template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegression<FPCAData,Integrator,ORDER,mydim,ndim>::computeDegreesOfFreedom(UInt output_index, SpMat& coeffmatrix)
{
	UInt nnodes = mesh_.num_nodes();
	UInt nlocations = fPCAData_.getNumberofObservations();

	Eigen::SparseLU<SpMat> solver;
	solver.compute(coeffmatrix);
	SpMat I(coeffmatrix.rows(),coeffmatrix.cols());
	I.setIdentity();
	SpMat coeff_inv = solver.solve(I);


	Real degrees=0;

		if(fPCAData_.isLocationsByNodes())
		{
			VectorXr d = coeff_inv.diagonal();

			for (auto i=0; i<nlocations;++i)
			{
				auto index_i = fPCAData_.getObservationsIndices()[i];
				degrees+=d(index_i);
			}
		}
		else
		{
			MatrixXr An(coeff_inv.topLeftCorner(nnodes, nnodes));
			MatrixXr S = psi_*An*psi_.transpose();
			for (auto i=0; i<nlocations;++i)
			{
				degrees+=S(i,i);
			}
		}



	std::cout<<"TRACE "<<degrees<<std::endl;

	_dof[output_index] = degrees;
}



//Implementation kept from Sangalli et al
template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegression<FPCAData,Integrator,ORDER,mydim,ndim>::smoothLaplace()
{
	//std::cout<<"Laplace Penalization - Order: "<<ORDER<<std::endl;

	//UInt ndata=regressionData_.getObservations().size();
	UInt nnodes=mesh_.num_nodes();

	FiniteElement<Integrator, ORDER,mydim,ndim> fe;

	typedef EOExpr<Mass> ETMass;
	typedef EOExpr<Stiff> ETStiff;

	Mass EMass;
	Stiff EStiff;

	ETMass mass(EMass);
	ETStiff stiff(EStiff);

	if(!fPCAData_.isLocationsByNodes())
	{
		setPsi();
	}

    if(!fPCAData_.isLocationsByNodes())
    {
    	getDataMatrix(DMat_);
    }
    else
    {
    	getDataMatrixByIndices(DMat_);
    }
    //std::cout<<"Block Data"<<DMat_<<std::endl;


    Assembler::operKernel(stiff, mesh_, fe, AMat_);
    Assembler::operKernel(mass, mesh_, fe, MMat_);

    VectorXr rightHandData;
    getRightHandData(rightHandData);
    _b = VectorXr::Zero(2*nnodes);
    _b.topRows(nnodes)=rightHandData;

    _solution.resize(fPCAData_.getLambda().size());
    _dof.resize(fPCAData_.getLambda().size());

#pragma omp parallel for
    for(UInt i = 0; i<fPCAData_.getLambda().size(); ++i)
	{
    	//build(tripletsData_,(-regressionData_.getLambda())*stiff, (-regressionData_.getLambda())*mass, righthand, forcing);

    	Real lambda = fPCAData_.getLambda()[i];
    	SpMat AMat_lambda = (-lambda)*AMat_;
    	SpMat MMat_lambda = (-lambda)*MMat_;

    	//this->buildCoeffMatrix(DMat_, AMat_lambda, MMat_lambda);

	SpMat coeffmatrix_lambda;
	buildCoeffMatrix(DMat_, AMat_lambda, MMat_lambda, coeffmatrix_lambda);


    	//Appling border conditions if necessary
    	if(fPCAData_.getDirichletIndices().size() != 0){
    		addDirichletBC(fPCAData_.getDirichletIndices(), fPCAData_.getDirichletValues(), coeffmatrix_lambda);}

    	//prova.solveSystem<SpConjGrad>();
    	this-> template solve<SpLU>(i, coeffmatrix_lambda);
	if(fPCAData_.computeDOF()){
		computeDegreesOfFreedom(i, coeffmatrix_lambda);
    	}else
    		_dof[i] = -1;

	}

}

/*
//Implementation kept from Sangalli et al
//Templatization in this method is only partial! For the moment it only works in
//the case mydim=ndim=2
template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegression<FPCAData,Integrator,ORDER,mydim,ndim>::smoothEllipticPDE()
{

if(mydim!=2 || ndim !=2){

	#ifdef R_VERSION_
		Rprintf("ERROR: these dimensions are not yet implemented, for the moment smoothEllipticPDE is available for mydim=ndim=2");
	#else
		std::cout << "ERROR: these dimensions are not yet implemented, for the moment smoothEllipticPDE is available for mydim=ndim=2\n";
	#endif

}else{
	//std::cout<<"Elliptic PDE Penalization - Order: "<<ORDER<<std::endl;

	//UInt ndata=regressionData_.getObservations().size();
	UInt nnodes=mesh_.num_nodes();

	FiniteElement<Integrator, ORDER,mydim,ndim> fe;

	typedef EOExpr<Mass> ETMass;
	typedef EOExpr<Stiff> ETStiff;
	typedef EOExpr<Grad> ETGrad;

	Mass EMass;
	Stiff EStiff;
	Grad EGrad;

	ETMass mass(EMass);
	ETStiff stiff(EStiff);
	ETGrad grad(EGrad);

	if(!regressionData_.isLocationsByNodes())
	{
		//std::cout<<"HERE";
		setPsi();
	}

	if(!regressionData_.getCovariates().rows() == 0)
	{
		setH();
		setQ();
	}

    if(!regressionData_.isLocationsByNodes())
    {
    	getDataMatrix(DMat_);
    }
    else
    {
    	getDataMatrixByIndices(DMat_);
    }
    //std::cout<<"Block Data"<<DMat<<std::endl;



    const Real& c = regressionData_.getC();
    const Eigen::Matrix<Real,2,2>& K = regressionData_.getK();
    const Eigen::Matrix<Real,2,1>& beta = regressionData_.getBeta();
    Assembler::operKernel(c*mass+stiff[K]+dot(beta,grad), mesh_, fe, AMat_);
    Assembler::operKernel(mass, mesh_, fe, MMat_);

    VectorXr rightHandData;
    getRightHandData(rightHandData);
    _b = VectorXr::Zero(2*nnodes);
    _b.topRows(nnodes)=rightHandData;
    //std::cout<<"b vector"<<_b;

    _solution.resize(regressionData_.getLambda().size());
    _dof.resize(regressionData_.getLambda().size());
#pragma omp parallel for
    for(UInt i = 0; i<regressionData_.getLambda().size(); ++i)
	{
    	//build(tripletsData_,(-regressionData_.getLambda())*stiff, (-regressionData_.getLambda())*mass, righthand, forcing);

    	Real lambda = regressionData_.getLambda()[i];
    	SpMat AMat_lambda = (-lambda)*AMat_;
    	SpMat MMat_lambda = (-lambda)*MMat_;
    	//this->buildCoeffMatrix(DMat_, AMat_lambda, MMat_lambda);

	SpMat coeffmatrix_lambda;
	buildCoeffMatrix(DMat_, AMat_lambda, MMat_lambda, coeffmatrix_lambda);
    	//std::cout<<"AMat"<<std::endl<<_coeffmatrix;


    	//Appling border conditions if necessary
    	if(regressionData_.getDirichletIndices().size() != 0)
    		addDirichletBC(regressionData_.getDirichletIndices(), regressionData_.getDirichletValues(), coeffmatrix_lambda);

    	//prova.solveSystem<SpConjGrad>();
    	this->template solve<SpLU>(i, coeffmatrix_lambda);
    	if(regressionData_.computeDOF())
    		computeDegreesOfFreedom(i,coeffmatrix_lambda);
    	else
    		_dof[i] = -1;

	}

}

}


//Implementation kept from Sangalli et al
//Templatization in this method is only partial! For the moment it only works in
//the case mydim=ndim=2
template<typename InputHandler, typename Integrator, UInt ORDER,UInt mydim, UInt ndim>
void MixedFERegression<InputHandler,Integrator,ORDER, mydim, ndim>::smoothEllipticPDESpaceVarying()
{

if(mydim!=2 || ndim !=2){

	#ifdef R_VERSION_
		Rprintf("ERROR: these dimensions are not yet implemented, for the moment SpaceVarying is available for mydim=ndim=2");
	#else
		std::cout << "ERROR: these dimensions are not yet implemented, for the moment SpaceVarying is available for mydim=ndim=2\n";
	#endif

}else{
	//std::cout<<"Space-varying Coefficient - Elliptic PDE Penalization - Order: "<<ORDER<<std::endl;
		//UInt ndata=regressionData_.getObservations().size();
		UInt nnodes=mesh_.num_nodes();

		FiniteElement<Integrator, ORDER,mydim,ndim> fe;

		typedef EOExpr<Mass> ETMass;
		typedef EOExpr<Stiff> ETStiff;
		typedef EOExpr<Grad> ETGrad;

		Mass EMass;
		Stiff EStiff;
		Grad EGrad;

		ETMass mass(EMass);
		ETStiff stiff(EStiff);
		ETGrad grad(EGrad);

	if(!regressionData_.isLocationsByNodes())
	{
		//std::cout<<"HERE";
		setPsi();
	}

	if(!regressionData_.getCovariates().rows() == 0)
	{
		setH();
		setQ();
	}

    if(!regressionData_.isLocationsByNodes())
    {
    	getDataMatrix(DMat_);
    }
    else
    {
    	getDataMatrixByIndices(DMat_);
    }
    //std::cout<<"Block Data"<<DMat_<<std::endl;

    const Reaction& c = regressionData_.getC();
    const Diffusivity& K = regressionData_.getK();
    const Advection& beta = regressionData_.getBeta();
    Assembler::operKernel(c*mass+stiff[K]+dot(beta,grad), mesh_, fe, AMat_);
    Assembler::operKernel(mass, mesh_, fe, MMat_);

    const ForcingTerm& u = regressionData_.getU();
    //for(auto i=0;i<18;i++) std::cout<<u(i)<<std::endl;
    VectorXr forcingTerm;
    Assembler::forcingTerm(mesh_,fe, u, forcingTerm);

    VectorXr rightHandData;
    getRightHandData(rightHandData);
    //_b.resize(2*nnodes);
    _b = VectorXr::Zero(2*nnodes);
    _b.topRows(nnodes)=rightHandData;

    _solution.resize(regressionData_.getLambda().size());
    _dof.resize(regressionData_.getLambda().size());
#pragma omp parallel for
    for(UInt i = 0; i<regressionData_.getLambda().size(); ++i)
	{
    	//build(tripletsData_,(-regressionData_.getLambda())*stiff, (-regressionData_.getLambda())*mass, righthand, forcing);

    	Real lambda = regressionData_.getLambda()[i];
    	SpMat AMat_lambda = (-lambda)*AMat_;
    	SpMat MMat_lambda = (-lambda)*MMat_;
    	//this->buildCoeffMatrix(DMat_, AMat_lambda, MMat_lambda);

	SpMat coeffmatrix_lambda;
	buildCoeffMatrix(DMat_, AMat_lambda, MMat_lambda, coeffmatrix_lambda);

        //std::cout<<"Forcing Term "<<std::cout<<forcingTerm<<"END";
        _b.bottomRows(nnodes)=lambda*forcingTerm;
        //std::cout<<"b vector"<<_b;
    	//std::cout<<"AMat"<<std::endl<<_coeffmatrix;


    	//Appling border conditions if necessary
    	if(regressionData_.getDirichletIndices().size() != 0)
    		addDirichletBC(regressionData_.getDirichletIndices(), regressionData_.getDirichletValues(), coeffmatrix_lambda);

    	//prova.solveSystem<SpConjGrad>();

    	//std::cout<< _coeffmatrix;
    	this-> template solve<SpLU>(i, coeffmatrix_lambda);
    	if(regressionData_.computeDOF())
    		computeDegreesOfFreedom(i,coeffmatrix_lambda);
    	else
    		_dof[i] = -1;

	}
}
}
*/

//solve sparse system with P method

template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
template <typename P>
void MixedFERegression<FPCAData,Integrator,ORDER,mydim,ndim>::solve(UInt output_index, SpMat& coeffmatrix)
{
	//std::cout<<this->_coeffmatrix;
	this->_solution[output_index].resize(coeffmatrix.rows());
	P::solve(coeffmatrix,this->_b,this->_solution[output_index]);
}

#endif
