#ifndef __ISOFEM_H__
#define __ISOFEM_H__

#include "MESH.h"
#include "QuadratureUtils.h"

#include <map>
#include <Eigen/Sparse>
#include <Eigen/Core>
#include "omp.h"
#include <memory>
#include <functional>
#include <algorithm>
#include <iostream>
#include <chrono>
#include <stdlib.h>
#include <iomanip>

using Triplet = Eigen::Triplet<double>;

using NdxArray = std::vector<std::vector<Eigen::MatrixXd>>;

using JArray = std::vector<std::vector<Eigen::Matrix2d>>;

using JacArray = std::vector<std::vector<double>>;



class ISOFEM
{
public:
	ISOFEM(std::string);
	
	ISOFEM(std::string, double, double);
	
	ISOFEM(std::string, double, double, int);
	
	void StaticAnalysis(const std::function<Eigen::Vector2d(Eigen::RowVector2d &)>&);
	
	void NonLocalStaticAnalysis(const std::function<Eigen::Vector2d(Eigen::RowVector2d &)>  &,
								const std::function<double(const double, const double)> &);
	
	const double &SystemEnergy() const;	
	
protected:
	
	Eigen::SparseMatrix<double> K;
	Eigen::MatrixXd	Ke,        // матирца жесткости элемента
					BT,        // оператор связи деформаций и перемещений (транспонированный)
					B;
	Eigen::VectorXd F;         // глобальный вектор правой части
	Eigen::VectorXd Fe;        // вектор правой части элемента
	Eigen::Matrix2d Jmatr;     // матрица Якоби в двумерном случае
	Eigen::RowVector2d Jvec;   // матрица Якоби в одномернмом случае
	Eigen::Matrix3d D;         // матрица Гука
	
	Eigen::VectorXd U;         // глобальный вектор Перемещений
	Eigen::VectorXd EpsiXX,  EpsiYY,  EpsiXY,  EpsiZZ;
	Eigen::VectorXd SigmaXX, SigmaYY, SigmaXY, SigmaVM;
	
	std::map<int, int> cN;
	
	std::vector<std::vector<int>> RnnArr;
	
	double p1, p2;             // параметры вклада локальных и нелокальных эффектов
	
	double L;                  // радиус влияния 
	
	double DeformationEnergy = 0.0;
	
	int NumberOfNodes, NumberOfElements, NodesPerElement, UxBoundNumber, 
		UyBoundNumber, UxUyBoundNumber, LoadBoundNumber;
		
	int HighOrder;
	
	FINITE_ELEMENT::MESH MESH;
	
	void ConstructRightPart(const std::function<Eigen::Vector2d(Eigen::RowVector2d &)> &);

	void ConstructStiffMatr(std::vector<Triplet> &);
	
	void ConstructStiffMatr(std::vector<Triplet> &, NdxArray &, JArray &, JacArray &);
	
	void ConstructStiffMatr(std::vector<Triplet> &, const NdxArray &, const JArray &, const JacArray &, 
							const std::function<double(const double, const double)> &);
	
	void ApplyingConstraints();
	
	std::vector<int> UniqueNodes(const Eigen::MatrixXi &);
	
	inline void ComputeB(const Eigen::MatrixXd &);
	
 	inline void ComputeKe(const int &, const Eigen::RowVectorXd &, const Eigen::MatrixXd &,
						  const std::vector<Eigen::MatrixXd> &,          Eigen::MatrixXd &);
	
	inline void AddKeToTriplet(std::vector<Triplet> &, const std::vector<int> &);
	
	void ComputeStrains();
	
	void ComputeStress();
	
	void ComputeStress(const JacArray &, const NdxArray &, const std::function<double(const double, const double)> &);
	
	void SetFiniteElement(const unsigned, std::shared_ptr<FINITE_ELEMENT::FiniteElement> &, unsigned &);
	
	unsigned SetFiniteElement(const unsigned, std::shared_ptr<FINITE_ELEMENT::FiniteElement> &);
	
	std::map<int, int> Counter(const Eigen::MatrixXi &);
	
	void NaiveRnnSearch(std::vector<std::vector<int>> &, const double &);
	
	virtual void WriteToVTK(std::string);
	
	
};

	


#endif
