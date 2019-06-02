#include "QuadratureUtils.h"

namespace FINITE_ELEMENT 
{
	QuadratureUtils::QuadratureUtils(FiniteElement *pE, unsigned const Order): pE(pE), Order(Order)
	{		
		SetQuadratureRule();
		SetShFVal();
	}
	
	const Eigen::RowVectorXd &QuadratureUtils::Weights() const { return _Weights; }

	const Eigen::MatrixXd &QuadratureUtils::QuadPoint() const { return _QuadPoint; }
	
	const std::vector<Eigen::RowVectorXd> &QuadratureUtils::Narr() const { return _Narr; }

	const std::vector<Eigen::MatrixXd> &QuadratureUtils::NGradArr() const { return _NGradArr; }
	
	const int QuadratureUtils::NumberOfQP() const {return _Weights.size(); }
	
	
	QuadratureUtils::~QuadratureUtils()
	{
		pE = nullptr;
	};
	

	void QuadratureUtils::SetQuadratureRule()
	{
		unsigned Dim = pE->Dimension();
		
		if(Dim == 1)
		{	
			switch (Order)
			{
				case 1:
				{
					_Weights.resize(1); _Weights(0) = 2.0;
					
					_QuadPoint.resize(1,1);
					
					_QuadPoint(0,0) = 0.0;
					
					break;
					
				}
				case 3:
				{
					
					_Weights.resize(2); _Weights(0) = 1.0; _Weights(1) = 1.0; 
					
					_QuadPoint.resize(2,1);
					
					_QuadPoint(0,0) = -1.0/sqrt(3.0);
					_QuadPoint(1,0) =  1.0/sqrt(3.0);
					

					break;
					
				}
				default:
					break;
			}
		}
		
		if(Dim == 2)
		{	switch (Order)
			{	
				case 3:
				{
					_Weights.resize(4); 
					_Weights(0) = 1.0; _Weights(1) = 1.0;
					_Weights(2) = 1.0; _Weights(3) = 1.0;
					
					_QuadPoint.resize(4,2);
					
					_QuadPoint(0,0) = -1.0/sqrt(3.0); _QuadPoint(0,1) = -1.0/sqrt(3.0);
					_QuadPoint(1,0) =  1.0/sqrt(3.0); _QuadPoint(1,1) = -1.0/sqrt(3.0);
					_QuadPoint(2,0) =  1.0/sqrt(3.0); _QuadPoint(2,1) =  1.0/sqrt(3.0);
					_QuadPoint(3,0) = -1.0/sqrt(3.0); _QuadPoint(3,1) =  1.0/sqrt(3.0);
					
					break;
				}
				
				case 5:
				{
					double q = sqrt(0.6);
					_Weights.resize(9); 
					_QuadPoint.resize(9,2);
					
					if(pE->NodesPerElement() == 4)
					{	
						
						_Weights(0) = 25.0/81.0;
						_Weights(1) = 25.0/81.0;
						_Weights(2) = 25.0/81.0;
						_Weights(3) = 25.0/81.0;
						_Weights(4) = 40.0/81.0;
						_Weights(5) = 40.0/81.0;
						_Weights(6) = 40.0/81.0;
						_Weights(7) = 40.0/81.0;
						_Weights(8) = 64.0/81.0;
						

						_QuadPoint(0,0) =  -q;   _QuadPoint(0,1) =  -q;
						_QuadPoint(1,0) =   q;   _QuadPoint(1,1) =  -q;
						_QuadPoint(2,0) =   q;   _QuadPoint(2,1) =   q;
						_QuadPoint(3,0) =  -q;   _QuadPoint(3,1) =   q;
						_QuadPoint(4,0) =  -q;   _QuadPoint(4,1) = 0.0;
						_QuadPoint(5,0) = 0.0;   _QuadPoint(5,1) =  -q;
						_QuadPoint(6,0) =   q;   _QuadPoint(6,1) = 0.0;
						_QuadPoint(7,0) = 0.0;   _QuadPoint(7,1) =   q;
						_QuadPoint(8,0) = 0.0;   _QuadPoint(8,1) = 0.0;
						
					} else if(pE->NodesPerElement() == 8)
					{

						_Weights(0) = 40.0/81.0;
						_Weights(1) = 40.0/81.0;
						_Weights(2) = 40.0/81.0;
						_Weights(3) = 40.0/81.0;
						_Weights(4) = 25.0/81.0;
						_Weights(5) = 25.0/81.0;
						_Weights(6) = 25.0/81.0;
						_Weights(7) = 25.0/81.0;
						_Weights(8) = 64.0/81.0;
						

						_QuadPoint(0,0) =  -q;   _QuadPoint(0,1) = 0.0;
						_QuadPoint(1,0) = 0.0;   _QuadPoint(1,1) =  -q;
						_QuadPoint(2,0) =   q;   _QuadPoint(2,1) = 0.0;
						_QuadPoint(3,0) = 0.0;   _QuadPoint(3,1) =   q;
						_QuadPoint(4,0) =  -q;   _QuadPoint(4,1) =  -q;
						_QuadPoint(5,0) =   q;   _QuadPoint(5,1) =  -q;
						_QuadPoint(6,0) =   q;   _QuadPoint(6,1) =   q;
						_QuadPoint(7,0) =  -q;   _QuadPoint(7,1) =   q;
						_QuadPoint(8,0) = 0.0;   _QuadPoint(8,1) = 0.0;
					}
					break;
				}
				
				default:
				{
					break;
				}
			}
	}
	};
	
	void QuadratureUtils::SetShFVal()
	{				
				_Narr.resize(_Weights.size());
				_NGradArr.resize(_Weights.size());
				
				for(size_t i = 0; i < _Weights.size(); ++i)
				{
					
					Eigen::MatrixXd QP = _QuadPoint.row(i);
					
					_Narr[i] = pE->ShFunc(QP);
					
					_NGradArr[i] = pE->ShFuncGrad(QP);
				}
				

	}
		
	
		
}
