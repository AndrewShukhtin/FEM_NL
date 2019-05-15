#include <iostream>
#include "QuadratureUtils.h"
#include "ISOFEMSOL.h"
#include "MESH.h"




int main(int argc, char *argv[])
{
	assert(argc > 1);
	
	{
		ISOFEMSOL Task(argv[1], 1.0);	
	
		Task.Static_Analysis([](Eigen::RowVector2d &X)->Eigen::Vector2d {Eigen::Vector2d V; V(0) = -10.0; V(1) = 0.0;  return V;});
	}
	{
		ISOFEMSOL Task(argv[1], 0.8, 100.0);
		Task.NonLocal_Static_Analysis(
		[](Eigen::RowVector2d &X)->Eigen::Vector2d 
		{
			Eigen::Vector2d V; 
			V(0) = -10.0; V(1) = 0.0;  
			return V;
			
			},
			
		[](Eigen::RowVector2d &X, Eigen::RowVector2d &Y, const double &R)->double 
		{
			Eigen::Vector2d V = X - Y; 
			
			double l = (R - 4.0);
			
			double r = V.norm();
			
			double res = 1.0 - r/(l);
			
			//std::cout << "\nres 1 = " << res << std::endl;
			
			res = (res + fabs(res)) * 1.5 / (l * l * 3.141592653589793238);
			
			//std::cout << "\nres 2 = " << res << std::endl;
			
			return res;
			});
 	}
	return 0;
}
