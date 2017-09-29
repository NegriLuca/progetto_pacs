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
}

void RegressionData::printObservations(std::ostream & out) const
{

	for(auto i=0;i<observations_.cols(); i++)
	{
		for(auto j=0;j<observations_.rows();j++)
		{
		out<<i<<"\t"<<observations_(j,i)<<std::endl;
		}
	}
}


void RegressionData::printLocations(std::ostream & out) const
{

	for(std::vector<Point>::size_type i=0;i<locations_.size(); i++)
	{
		locations_[i].print(out);
		//std::cout<<std::endl;
	}
}

#endif

