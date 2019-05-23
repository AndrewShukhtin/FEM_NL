#include <iostream>
#include "QuadratureUtils.h"
#include "ISOFEMSOL.h"
#include "MESH.h"




int main(int argc, char *argv[])
{
	assert(argc > 1);
	
	unsigned N = std::stoi(argv[2]);
	
	double p1 = std::stod(argv[3]);
	
	double L = std::stod(argv[4]);
	
	const auto Load = [](Eigen::RowVector2d &X)->Eigen::Vector2d 
	{
		Eigen::Vector2d V; 
		V(0) = -10.0; 
		V(1) = 0.0;  
		return V;
		
	};
	
	{
		ISOFEMSOL Task(argv[1], 1.0);	
	
		Task.Static_Analysis(Load);
	}
	
	if(argc > 2)
	{
		
		const auto phi1 =[](const double &r, const double &L)->double 
		{
			double res = r*r/(L * L);

			res = exp(-res) / (L * L * M_PI);

			return res;
		};
		
		const auto phi2 =[](const double &r, const double &L)->double 
		{
			double res = 1.0 - r/L;

			res = (res + fabs(res)) * 1.5 / (L * L * M_PI);

			return res;
		};
		
		
		ISOFEMSOL Task(argv[1], p1, L);
		
		switch(N)
		{
			case 1:
				Task.NonLocal_Static_Analysis(Load, phi1);
				break;
			case 2:
				Task.NonLocal_Static_Analysis(Load, phi2);
				break;
			default:
				break;
				
		}
	
 	}
	return 0;
}
