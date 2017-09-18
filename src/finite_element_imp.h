#ifndef __FINITE_ELEMENT_IMP_HPP__
#define __FINITE_ELEMENT_IMP_HPP__

#include<iostream>

//! Implementazione EF mydim=2, ndim=2


template <class Integrator, UInt ORDER>
FiniteElement<Integrator, ORDER,2,2>::FiniteElement()
{
	//Set the properties of the reference element
	std::vector<Point> reference_nodes;
	reference_nodes.push_back(Point(0,0));
	reference_nodes.push_back(Point(1,0));
	reference_nodes.push_back(Point(0,1));

	reference_ = Triangle<ORDER*3,2,2> (Id(0), reference_nodes);

	//How it will be used, it does not depend on J^-1 -> set one time
	setPhiMaster();
	setPhiDerMaster();
}


template <class Integrator, UInt ORDER>
void FiniteElement<Integrator, ORDER,2,2>::updateElement(Triangle<ORDER*3,2,2> t)
{
	t_ = t;

	//it does depend on J^-1 -> set for each triangle
	setInvTrJPhiDerMaster();

}

template <class Integrator, UInt ORDER>
void FiniteElement<Integrator, ORDER,2,2>::setPhiMaster()
{
	Eigen::Matrix<Real,3*ORDER,1> coefficients;
	for (auto i=0; i < 3*ORDER; i++)
	{
		coefficients = MatrixXr::Zero(3*ORDER,1);
		coefficients(i) = 1;
		for (auto iq=0; iq < Integrator::NNODES; iq++)
		{
			Real phi = evaluate_point<ORDER>(reference_,Integrator::NODES[iq],coefficients);
			phiMapMaster_(i,iq) = phi;
		}
	}
}

template <class Integrator, UInt ORDER>
void FiniteElement<Integrator, ORDER,2,2>::setPhiDerMaster()
{
	Eigen::Matrix<Real,3*ORDER,1> coefficients;
	Eigen::Matrix<Real,2,1> der;
	Eigen::Matrix<Real,2,1> der_transf;

	for (auto i=0; i < 3*ORDER; i++)
	{
		coefficients = MatrixXr::Zero(3*ORDER,1);
		coefficients(i) = 1;
		for (auto iq=0; iq < Integrator::NNODES; iq++)
		{
			der = evaluate_der_point<ORDER>(reference_,Integrator::NODES[iq],coefficients);
			// we need J^(-1) nabla( phi)
			//der_transf = t_.getM_invJ()*der;
			phiDerMapMaster_(i,iq*2) = der[0];
			phiDerMapMaster_(i,iq*2+1) = der[1];
		}
	}
}

template <class Integrator, UInt ORDER>
Real FiniteElement<Integrator, ORDER,2,2>::phiMaster(UInt i, UInt iq) const
{
	return phiMapMaster_(i, iq);
}

template <class Integrator, UInt ORDER>
Real FiniteElement<Integrator, ORDER,2,2>::phiDerMaster(UInt i, UInt ic, UInt iq) const
{
	return phiDerMapMaster_(i, iq*2 + ic);
}

template <class Integrator, UInt ORDER>
Real FiniteElement<Integrator, ORDER,2,2>::invTrJPhiDerMaster(UInt i, UInt ic, UInt iq) const
{
	return invTrJPhiDerMapMaster_(i, iq*2 + ic);
}

template <class Integrator, UInt ORDER>
void FiniteElement<Integrator, ORDER,2,2>::setInvTrJPhiDerMaster()
{
	Eigen::Matrix<Real,2,1> der;
	Eigen::Matrix<Real,2,1> der_transf;

	for (auto i=0; i < 3*ORDER; i++)
	{
		for (auto iq=0; iq < Integrator::NNODES; iq++)
		{
			der[0] = phiDerMaster(i, 0, iq);
			der[1] = phiDerMaster(i, 1, iq);
			// we need J^(-T) nabla( phi)
			der_transf = t_.getM_invJ().transpose()*der;
			invTrJPhiDerMapMaster_(i,iq*2) = der_transf[0];
			invTrJPhiDerMapMaster_(i,iq*2+1) = der_transf[1];
		}
	}
}



//! Implementazione EF mydim=2, ndim=3


template <class Integrator, UInt ORDER>
FiniteElement<Integrator, ORDER,2,3>::FiniteElement()
{
	//Set the properties of the reference element
	std::vector<Point> reference_nodes;
	reference_nodes.push_back(Point(0,0));
	reference_nodes.push_back(Point(1,0));
	reference_nodes.push_back(Point(0,1));

	reference_ = Triangle<ORDER*3,2,2> (Id(0), reference_nodes);

	//How it will be used, it does not depend on J^-1 -> set one time
	setPhiMaster();
	setPhiDerMaster();
}


template <class Integrator, UInt ORDER>
void FiniteElement<Integrator, ORDER,2,3>::updateElement(Triangle<ORDER*3,2,3> t)
{
	t_ = t;

	//it does depend on J^-1 -> set for each triangle
	metric_ = t.getMetric();

}

template <class Integrator, UInt ORDER>
void FiniteElement<Integrator, ORDER,2,3>::setPhiMaster()
{
	Eigen::Matrix<Real,3*ORDER,1> coefficients;
	for (auto i=0; i < 3*ORDER; i++)
	{
		coefficients = MatrixXr::Zero(3*ORDER,1);
		coefficients(i) = 1;
		for (auto iq=0; iq < Integrator::NNODES; iq++)
		{
			Real phi = evaluate_point<ORDER>(reference_,Integrator::NODES[iq],coefficients);
			phiMapMaster_(i,iq) = phi;
		}
	}
}

template <class Integrator, UInt ORDER>
void FiniteElement<Integrator, ORDER,2,3>::setPhiDerMaster()
{
	Eigen::Matrix<Real,3*ORDER,1> coefficients;
	Eigen::Matrix<Real,2,1> der;
	Eigen::Matrix<Real,2,1> der_transf;

	for (auto i=0; i < 3*ORDER; i++)
	{
		coefficients = MatrixXr::Zero(3*ORDER,1);
		coefficients(i) = 1;
		for (auto iq=0; iq < Integrator::NNODES; iq++)
		{
			der = evaluate_der_point<ORDER,2,2>(reference_,Integrator::NODES[iq],coefficients);
			// we need J^(-1) nabla( phi)
			//der_transf = t_.getM_invJ()*der;
			phiDerMapMaster_(i,iq*2) = der[0];
			phiDerMapMaster_(i,iq*2+1) = der[1];
		}
	}
}

template <class Integrator, UInt ORDER>
Real FiniteElement<Integrator, ORDER,2,3>::phiMaster(UInt i, UInt iq) const
{
	return phiMapMaster_(i, iq);
}

template <class Integrator, UInt ORDER>
Real FiniteElement<Integrator, ORDER,2,3>::phiDerMaster(UInt i, UInt ic, UInt iq) const
{
	return phiDerMapMaster_(i, iq*2 + ic);
}

/*template <class Integrator, UInt ORDER>
Real FiniteElement<Integrator, ORDER,2,2>::invTrJPhiDerMaster(UInt i, UInt ic, UInt iq) const
{
	return invTrJPhiDerMapMaster_(i, iq*2 + ic);
}

template <class Integrator, UInt ORDER>
void FiniteElement<Integrator, ORDER,2,2>::setInvTrJPhiDerMaster()
{
	Eigen::Matrix<Real,2,1> der;
	Eigen::Matrix<Real,2,1> der_transf;

	for (auto i=0; i < 3*ORDER; i++)
	{
		for (auto iq=0; iq < Integrator::NNODES; iq++)
		{
			der[0] = phiDerMaster(i, 0, iq);
			der[1] = phiDerMaster(i, 1, iq);
			// we need J^(-T) nabla( phi)
			der_transf = t_.getM_invJ().transpose()*der;
			invTrJPhiDerMapMaster_(i,iq*2) = der_transf[0];
			invTrJPhiDerMapMaster_(i,iq*2+1) = der_transf[1];
		}
	}
}*/

#endif
