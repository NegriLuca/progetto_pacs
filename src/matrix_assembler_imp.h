#ifndef MATRIX_ASSEMBLER_IMP_H_
#define MATRIX_ASSEMBLER_IMP_H_


template<UInt ORDER, typename Integrator, typename A>
void Assembler::operKernel(EOExpr<A> oper,const MeshHandler<ORDER,2,2>& mesh,
	                     FiniteElement<Integrator, ORDER,2,2>& fe, SpMat& OpMat)
{
	Real eps = 2.2204e-016,
		 tolerance = 10 * eps;
	std::vector<coeff> triplets;


  	for(auto t=0; t<mesh.num_triangles(); t++)
  	{
		fe.updateElement(mesh.getTriangle(t));

		// Vector of vertices indices (link local to global indexing system)
		std::vector<UInt> identifiers;
		identifiers.resize(ORDER*3);
		for( auto q=0; q<ORDER*3; q++)
			identifiers[q]=mesh.getTriangle(t)[q].id();

		//localM=localMassMatrix(currentelem);
		for(int i = 0; i < 3*ORDER; i++)
		{
			for(int j = 0; j < 3*ORDER; j++)
			{
				Real s=0;

				for(int l = 0;l < Integrator::NNODES; l++)
				{
					s += oper(fe,i,j,l) * fe.getDet() * fe.getAreaReference()* Integrator::WEIGHTS[l];
				}
			  triplets.push_back(coeff(identifiers[i],identifiers[j],s));
			}
		}

	}

  	UInt nnodes = mesh.num_nodes();
  	OpMat.resize(nnodes, nnodes);
	OpMat.setFromTriplets(triplets.begin(),triplets.end());
	OpMat.prune(tolerance);
}

template<UInt ORDER, typename Integrator>
void Assembler::forcingTerm(const MeshHandler<ORDER,2,2>& mesh,
	                     FiniteElement<Integrator, ORDER,2,2>& fe, const ForcingTerm& u, VectorXr& forcingTerm)
{

	forcingTerm = VectorXr::Zero(mesh.num_nodes());

  	for(auto t=0; t<mesh.num_triangles(); t++)
  	{
		fe.updateElement(mesh.getTriangle(t));

		// Vector of vertices indices (link local to global indexing system)
		std::vector<UInt> identifiers;
				identifiers.resize(ORDER*3);

		for( auto q=0; q<ORDER*3; q++)
			identifiers[q]=mesh.getTriangle(t)[q].id();


		//localM=localMassMatrix(currentelem);
		for(int i = 0; i < 3*ORDER; i++)
		{
			Real s=0;

			for(int iq = 0;iq < Integrator::NNODES; iq++)
			{
				UInt globalIndex = fe.getGlobalIndex(iq);
				s +=  fe.phiMaster(i,iq)* u(globalIndex) * fe.getDet() * fe.getAreaReference()* Integrator::WEIGHTS[iq];//(*)
			}
			forcingTerm[identifiers[i]] += s;
		}

	}
}


//! Surface mesh implementation

template<UInt ORDER, typename Integrator, typename A>
void Assembler::operKernel(EOExpr<A> oper,const MeshHandler<ORDER,2,3>& mesh,
	                     FiniteElement<Integrator, ORDER,2,3>& fe, SpMat& OpMat)
{
	Real eps = 2.2204e-016,
		 tolerance = 10 * eps;
	std::vector<coeff> triplets;


  	for(auto t=0; t<mesh.num_triangles(); t++)
  	{
		fe.updateElement(mesh.getTriangle(t));

		// Vector of vertices indices (link local to global indexing system)
		std::vector<UInt> identifiers;
		identifiers.resize(ORDER*3);
		for( auto q=0; q<ORDER*3; q++)
			identifiers[q]=mesh.getTriangle(t)[q].id();

		//localM=localMassMatrix(currentelem);
		for(int i = 0; i < 3*ORDER; i++)
		{
			for(int j = 0; j < 3*ORDER; j++)
			{
				Real s=0;

				for(int l = 0;l < Integrator::NNODES; l++)
				{
					s += oper(fe,i,j,l) * std::sqrt(fe.getDet()) * fe.getAreaReference()* Integrator::WEIGHTS[l];
				}
			  triplets.push_back(coeff(identifiers[i],identifiers[j],s));
			}
		}

	}

  	UInt nnodes = mesh.num_nodes();
  	OpMat.resize(nnodes, nnodes);
	OpMat.setFromTriplets(triplets.begin(),triplets.end());
	OpMat.prune(tolerance);
}



template<UInt ORDER, typename Integrator>
void Assembler::forcingTerm(const MeshHandler<ORDER,2,3>& mesh,
	                     FiniteElement<Integrator, ORDER,2,3>& fe, const ForcingTerm& u, VectorXr& forcingTerm)
{

	forcingTerm = VectorXr::Zero(mesh.num_nodes());

  	for(auto t=0; t<mesh.num_triangles(); t++)
  	{
		fe.updateElement(mesh.getTriangle(t));

		// Vector of vertices indices (link local to global indexing system)
		std::vector<UInt> identifiers;
				identifiers.resize(ORDER*3);

		for( auto q=0; q<ORDER*3; q++)
			identifiers[q]=mesh.getTriangle(t)[q].id();


		//localM=localMassMatrix(currentelem);
		for(int i = 0; i < 3*ORDER; i++)
		{
			Real s=0;

			for(int iq = 0;iq < Integrator::NNODES; iq++)
			{
				UInt globalIndex = fe.getGlobalIndex(iq);
				s +=  fe.phiMaster(i,iq)* u(globalIndex) * std::sqrt(fe.getDet()) * fe.getAreaReference()* Integrator::WEIGHTS[iq];//(*)
			}
			forcingTerm[identifiers[i]] += s;
		}

	}

}




#endif
