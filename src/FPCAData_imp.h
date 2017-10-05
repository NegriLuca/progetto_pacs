#ifndef __FPCADATA_IMP_HPP__
#define __FPCADATA_IMP_HPP__

FPCAData::FPCAData(std::vector<Point>& locations, MatrixXr& observations, UInt order, std::vector<Real> lambda , std::vector<UInt>& bc_indices, std::vector<Real>& bc_values,UInt nPC, bool DOF):
					locations_(locations), observations_(observations), order_(order), lambda_(lambda),
					bc_values_(bc_values), bc_indices_(bc_indices),nPC_(nPC), DOF_(DOF)
{
	if(locations_.size()==0)
	{
		locations_by_nodes_= true;
		for(int i = 0; i<observations_.rows();++i) observations_indices_.push_back(i);
	}
	else
	{
		locations_by_nodes_= false;
	}
	//Initialize loadings vector
	Eigen::JacobiSVD<MatrixXr> svd(observations_.transpose(),Eigen::ComputeThinU|Eigen::ComputeThinV);
	loadings_=svd.matrixV().col(0);
}

#ifdef R_VERSION_
FPCAData::FPCAData(SEXP Rlocations, SEXP Robservations, SEXP Rorder, SEXP Rlambda,
	SEXP RBCIndices, SEXP RBCValues,SEXP RnPC, SEXP DOF)
{
	setLocations(Rlocations);
	setObservations(Robservations);

	//Initialize loadings vector
	Eigen::JacobiSVD<MatrixXr> svd(observations_.transpose(),Eigen::ComputeThinU|Eigen::ComputeThinV);
	loadings_=svd.matrixV().col(0);
	
	nPC_ = INTEGER(RnPC)[0];
	order_ =  INTEGER(Rorder)[0];
	DOF_ = INTEGER(DOF)[0];
	UInt length_indexes = Rf_length(RBCIndices);
    //for (UInt i = 0; i<length_indexes; ++i)  bc_indices_.push_back(INTEGER(RBCIndices)[i]);
    //for (UInt i = 0; i<length_indexes; ++i)  bc_values_.push_back(REAL(RBCValues)[i]);
	bc_indices_.assign(INTEGER(RBCIndices), INTEGER(RBCIndices) +  length_indexes);
    //conversion between R indices and c++ indices
	std::for_each(bc_indices_.begin(), bc_indices_.end(), [](int& i){i-=1;});
	bc_values_.assign(REAL(RBCValues),REAL(RBCValues) + Rf_length(RBCIndices));

    UInt length_lambda = Rf_length(Rlambda);
    for (UInt i = 0; i<length_lambda; ++i)  lambda_.push_back(REAL(Rlambda)[i]);

}


void FPCAData::setObservations(SEXP Robservations)
{
	n_ = INTEGER(Rf_getAttrib(Robservations, R_DimSymbol))[0];
	p_ = INTEGER(Rf_getAttrib(Robservations, R_DimSymbol))[1];
	observations_.resize(n_,p_);
	observations_indices_.reserve(n_);
	
	for(auto i=0; i<n_; ++i)
	{
		for(auto j=0; j<p_ ; ++j)
		{
			observations_(i,j)=REAL(Robservations)[i+ n_*j];
		}
	}

	if(locations_.size() == 0)
	{
		locations_by_nodes_ = true;
		for(auto i=0;i<n_;++i) observations_indices_.push_back(i);
	}
	else
	{
		locations_by_nodes_ = false;
	}

	//std::cout<<"Observations #"<<observations_.size()<<std::endl<<observations_<<std::endl;
	//for(auto i=0;i<observations_indices_.size();++i)	std::cout<<observations_indices_[i]<<std::endl;
}

void FPCAData::setLocations(SEXP Rlocations)
{
	n_ = INTEGER(Rf_getAttrib(Rlocations, R_DimSymbol))[0];
	if(n_>0){
		int ndim = INTEGER(Rf_getAttrib(Rlocations, R_DimSymbol))[1];

	  if (ndim == 2){
			for(auto i=0; i<n_; ++i)
			{
				locations_.emplace_back(REAL(Rlocations)[i+ n_*0],REAL(Rlocations)[i+ n_*1]);
			}
		}else{
			for(auto i=0; i<n_; ++i)
			{
				locations_.emplace_back(REAL(Rlocations)[i+ n_*0],REAL(Rlocations)[i+ n_*1],REAL(Rlocations)[i+ n_*2]);
			}
		}
	}
}

#endif


void FPCAData::printObservations(std::ostream & out) const
{

	for(auto i=0;i<observations_.rows(); i++)
	{
		for(auto j=0;j<observations_.cols();j++)
		{
		out<<observations_(i,j)<<"\t";
		}
		out<<std::endl;
	}
}


void FPCAData::printLocations(std::ostream & out) const
{

	for(std::vector<Point>::size_type i=0;i<locations_.size(); i++)
	{
		locations_[i].print(out);
		//std::cout<<std::endl;
	}
}

void FPCAData::setScores()
{
	scores_=observations_.transpose()*loadings_;
	scores_=scores_/scores_.norm();
}

void FPCAData::setDataForRegression()
{
	ObservationData_=observations_*scores_;
}

void FPCAData::setLoadings(SpMat psi, VectorXr f_sol)
{
	loadings_=psi.transpose()*f_sol;
}

#endif

