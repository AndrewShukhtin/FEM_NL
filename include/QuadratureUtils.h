#ifndef __QUADRATURE_UTILS_H__
#define __QUADRATURE_UTILS_H__

#include "Finite_Element_Class.h"

namespace FINITE_ELEMENT 
{
	class QuadratureUtils
	{
	public:
		QuadratureUtils(FiniteElement *pE);
	
		const Eigen::RowVectorXd &Weights() const;
	
		const Eigen::MatrixXd &QuadPoint() const;
		
		const std::vector<Eigen::RowVectorXd> &Narr() const;
	
		const std::vector<Eigen::MatrixXd> &NGradArr() const;
		
		const int NumberOfQP() const;
	
		~QuadratureUtils();
	
	private:
		
		FiniteElement *pE;
		
		Eigen::RowVectorXd _Weights;
		Eigen::MatrixXd _QuadPoint;
		
		std::vector<Eigen::RowVectorXd> _Narr;
		std::vector<Eigen::MatrixXd> _NGradArr;
		
		void SetQuadratureRule(unsigned const &Npe);
	
		void SetShFVal(FiniteElement *pE);

	};

}
#endif
