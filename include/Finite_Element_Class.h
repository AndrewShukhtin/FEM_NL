#ifndef __FINITE_ELEMENT_CLASS_H__
#define __FINITE_ELEMENT_CLASS_H__

#include<Eigen/Dense>
#include<vector>

namespace FINITE_ELEMENT 
{
	class FiniteElement 
	{
	public:
		
		FiniteElement(unsigned const Npe,unsigned const Dim): Npe(Npe), Dim(Dim) {};
		
		virtual unsigned const &NodesPerElement() const;	
		
		virtual unsigned const &Dimension() const;
		
		virtual Eigen::RowVectorXd ShFunc(Eigen::MatrixXd &qp)  = 0;
		
		virtual Eigen::MatrixXd ShFuncGrad(Eigen::MatrixXd &qp) = 0;
		
	protected:
		
		unsigned const Npe; // Количество узлов на элементе
		unsigned const Dim; // Размерность элемента
		
	};	
	
	
	class LinearElement: public FiniteElement
	{
	public:
		LinearElement();
		
		Eigen::RowVectorXd ShFunc(Eigen::MatrixXd &qp);
		
		Eigen::MatrixXd ShFuncGrad(Eigen::MatrixXd &qp);
				
	};
	
	class QuadraticElement: public FiniteElement
	{
	public:
		QuadraticElement();
		
		Eigen::RowVectorXd ShFunc(Eigen::MatrixXd &qp);
		
		Eigen::MatrixXd ShFuncGrad(Eigen::MatrixXd &qp);
				
	};
	
	
	
	class BilinearElement: public FiniteElement
	{
	public:
		BilinearElement();
		
		Eigen::RowVectorXd ShFunc(Eigen::MatrixXd &qp);
		
		Eigen::MatrixXd ShFuncGrad(Eigen::MatrixXd &qp);
				
	};
	
	class QuadraticSerendipElement: public FiniteElement
	{
	public:
		QuadraticSerendipElement();
		
		Eigen::RowVectorXd ShFunc(Eigen::MatrixXd &qp);
		
		Eigen::MatrixXd ShFuncGrad(Eigen::MatrixXd &qp);
				
	};
}




#endif
