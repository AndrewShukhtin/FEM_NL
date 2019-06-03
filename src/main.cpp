#include <iostream>
#include "QuadratureUtils.h"
#include "ISOFEM.h"
#include "FMMISOFEM.h"
#include "MESH.h"




int main(int argc, char *argv[])
{

	assert(argc > 1);
	
	Eigen::initParallel();
	
	omp_set_num_threads(4);
	Eigen::setNbThreads(4);
	
	std::cout<<"Eigen Version: " << EIGEN_WORLD_VERSION << "." << EIGEN_MAJOR_VERSION << "." << EIGEN_MINOR_VERSION<< std::endl;
	
	const auto Load = [](Eigen::RowVector2d &X)->Eigen::Vector2d 
	{
		Eigen::Vector2d V; 
		V(0) = -10.0; 
		V(1) = 0.0;  
		return V;
		
	};
	
	if (argc == 2)
	{
		ISOFEM Task(argv[1]);	
	
		Task.StaticAnalysis(Load);
	}
	
	if(argc > 2)
	{
		const auto A1 =[](const double r, const double L)->double 
		{
			double l = 0.8 * L;
			double res = r*r/(l * l);

			res = exp(-res) / (l * l * M_PI);

			return res;
		};
		
		const auto A2 =[](const double r, const double L)->double 
		{
			double res = 1.0 - r/L;

			res = (res + fabs(res)) * 1.5 / (L * L * M_PI);

			return res;
		};
		
				
		
		if(argc == 5)
		{
			int N = std::stoi(argv[2]);
	
			double p1 = std::stod(argv[3]);
		
			double L = std::stod(argv[4]);
			
			ISOFEM Task(argv[1], p1, L);
						
			switch(N)
			{
				case 1:
					Task.NonLocalStaticAnalysis(Load, A1);
					break;
				case 2:
					Task.NonLocalStaticAnalysis(Load, A2);
					break;
				default:
					break;
					
			}
		}
		
		if(argc == 6)
		{
			int N = std::stoi(argv[2]);
	
			double p1 = std::stod(argv[3]);
		
			double L = std::stod(argv[4]);
			
			std::string str = argv[5];
			
			if(str != "-fmm")
			{	
				int HighOrder = std::stoi(argv[5]);
				
				ISOFEM Task(argv[1], p1, L);
							
				switch(N)
				{
					case 1:
						Task.NonLocalStaticAnalysis(Load, A1);
						break;
					case 2:
						Task.NonLocalStaticAnalysis(Load, A2);
						break;
					default:
						break;
						
				}
			}else if(str == "-fmm")
			{
				FMMISOFEM Task(argv[1], p1, L);
				
				switch(N)
				{
					case 1:
						Task.FMMStaticAnalysis(Load, A1);
						break;
					case 2:
						Task.FMMStaticAnalysis(Load, A2);
						break;
					default:
						break;
				}
				
			}else
			{
				std::cout << "Error! Wrong paramaters!\n";
			}
		}
		
		if(argc == 7)
		{	
			int N = std::stoi(argv[2]);
	
			double p1 = std::stod(argv[3]);
		
			double L = std::stod(argv[4]);
		
			int HighOrder = std::stoi(argv[6]);
			
			FMMISOFEM Task(argv[1], p1, L, HighOrder);
	
			switch(N)
			{
				case 1:
					Task.FMMStaticAnalysis(Load, A1);
					break;
				case 2:
					Task.FMMStaticAnalysis(Load, A2);
					break;
				default:
					break;
					
			}
		}
	
 	}
 	
	return 0;
}
