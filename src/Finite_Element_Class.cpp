#include "Finite_Element_Class.h"


namespace FINITE_ELEMENT
{
	
	
	unsigned const &FiniteElement::NodesPerElement() const{ return Npe; }
	
	
	LinearElement::LinearElement(): FiniteElement(2) {};
	
	Eigen::RowVectorXd LinearElement::ShFunc(Eigen::MatrixXd &qp)
	{
		double xi = qp(0);
		Eigen::RowVectorXd Narr(2);
		
		Narr(0) = 0.5 * (1 - xi); 	Narr(1) = 0.5 * (1 + xi);
		
		return Narr;
	}
	
	Eigen::MatrixXd LinearElement::ShFuncGrad(Eigen::MatrixXd &qp)
	{
		double xi = qp(0);
		Eigen::MatrixXd NGradArr(1,2);
		
		NGradArr(0) = -0.5; 	NGradArr(1) = 0.5;
		
		return NGradArr;
	}
	
	QuadraticElement::QuadraticElement(): FiniteElement(3) {};
	
	Eigen::RowVectorXd QuadraticElement::ShFunc(Eigen::MatrixXd &qp)
	{
		double xi = qp(0);
		Eigen::RowVectorXd Narr(3);
		
		Narr(0) = -0.5 * xi * (1.0 - xi);  Narr(1) = 0.5 * xi * (1.0 + xi);  Narr(2) = (1.0 - xi) * (1.0 + xi);
		
		return Narr;
	}
	
	Eigen::MatrixXd QuadraticElement::ShFuncGrad(Eigen::MatrixXd &qp)
	{
		double xi = qp(0);
		Eigen::MatrixXd NGradArr(1,3);
		
		NGradArr(0) = -0.5 + xi;  NGradArr(1) = 0.5 + xi;  NGradArr(2) = -2.0 * xi;
		
		return NGradArr;
	}


	BilinearElement::BilinearElement(): FiniteElement(4) {};
	
	Eigen::RowVectorXd BilinearElement::ShFunc(Eigen::MatrixXd &qp) 
	{
		double xi = qp(0), eta = qp(1);
		Eigen::RowVectorXd Narr(4);
		
		Narr(0) = 0.25 * (1.0 - xi)*(1.0 - eta); 		Narr(1) = 0.25 * (1.0 + xi)*(1.0 - eta);
		Narr(2) = 0.25 * (1.0 + xi)*(1.0 + eta); 		Narr(3) = 0.25 * (1.0 - xi)*(1.0 + eta);
		
		return Narr;
	}
	
	Eigen::MatrixXd BilinearElement::ShFuncGrad(Eigen::MatrixXd &qp)
	{
		double xi = qp(0), eta = qp(1);
		
		Eigen::MatrixXd NGradArr(2,4);
		
		NGradArr(0,0) = -0.25 *(1. - eta); 			NGradArr(0,1) =  0.25 *(1. - eta);
		NGradArr(0,2) =  0.25 *(1. + eta); 			NGradArr(0,3) = -0.25 *(1. + eta);
		
		NGradArr(1,0) = -0.25 *(1. - xi); 			NGradArr(1,1) = -0.25 *(1. + xi);
		NGradArr(1,2) =  0.25 *(1. + xi); 			NGradArr(1,3) =  0.25 *(1. - xi);
		
		
		return  NGradArr ;
	}
	
	
	QuadraticSerendipElement::QuadraticSerendipElement(): FiniteElement(8) {};
	
	Eigen::RowVectorXd QuadraticSerendipElement::ShFunc(Eigen::MatrixXd &qp) 
	{
		double xi = qp(0), eta = qp(1);
		Eigen::RowVectorXd Narr(8);
		
		Narr(0) = 0.5 * (1.0 - xi) * (1.0 - eta * eta); 		Narr(1) = 0.5 * (1.0 - xi * xi) * (1.0 - eta);
		Narr(2) = 0.5 * (1.0 + xi) * (1.0 - eta * eta); 		Narr(3) = 0.5 * (1.0 - xi * xi) * (1.0 + eta);
		
		Narr(4) = 0.25 * (1.0 - xi) * (1.0 - eta) * (-xi - eta - 1.0);   Narr(5) = 0.25 * (1.0 + xi) * (1.0 - eta) * ( xi - eta - 1.0);
		Narr(6) = 0.25 * (1.0 + xi) * (1.0 + eta) * ( xi + eta - 1.0);   Narr(7) = 0.25 * (1.0 - xi) * (1.0 + eta) * (-xi + eta - 1.0);
		
		return Narr;
	}
	
	Eigen::MatrixXd QuadraticSerendipElement::ShFuncGrad(Eigen::MatrixXd &qp)
	{
		double xi = qp(0), eta = qp(1);
		
		Eigen::MatrixXd NGradArr(2,8);
		
		NGradArr(0,0) = -0.5 * (1.0 - eta * eta); 	
		NGradArr(0,1) = -1.0 * (1.0 - eta) * xi;
		NGradArr(0,2) =  0.5 * (1.0 - eta * eta);
		NGradArr(0,3) = -1.0 * (1.0 + eta) * xi;
		NGradArr(0,4) = -0.25 * (1.0 - eta) * (1.0 - xi) - 0.25 * (1.0 - eta) * (-1.0 - eta - xi);
		NGradArr(0,5) =  0.25 * (1.0 - eta) * (1.0 + xi) + 0.25 * (1.0 - eta) * (-1.0 - eta + xi);
		NGradArr(0,6) =  0.25 * (1.0 + eta) * (1.0 + xi) + 0.25 * (1.0 + eta) * (-1.0 + eta + xi);
		NGradArr(0,7) = -0.25 * (1.0 + eta) * (1.0 - xi) - 0.25 * (1.0 + eta) * (-1.0 + eta - xi);
		
		
		
		NGradArr(1,0) = -1.0 * eta * (1.0 - xi);
		NGradArr(1,1) = -0.5 * (1.0 - xi * xi);
		NGradArr(1,2) = -1.0 * eta * (1.0 + xi);			
		NGradArr(1,3) =  0.5 * (1.0 - xi * xi);
		NGradArr(1,4) = -0.25 * (1.0 - eta) * (1.0 - xi) - 0.25 * (1.0 - xi) * (-1.0 - eta - xi);
		NGradArr(1,5) = -0.25 * (1.0 - eta) * (1.0 + xi) - 0.25 * (1.0 + xi) * (-1.0 - eta + xi);
		NGradArr(1,6) =  0.25 * (1.0 + eta) * (1.0 + xi) + 0.25 * (1.0 + xi) * (-1.0 + eta + xi);
		NGradArr(1,7) =  0.25 * (1.0 + eta) * (1.0 - xi) + 0.25 * (1.0 - xi) * (-1.0 + eta - xi);
		
		return  NGradArr ;
	}
	
}
