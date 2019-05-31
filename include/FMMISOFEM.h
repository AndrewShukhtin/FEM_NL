#ifndef __FMMISOFEM_H__
#define __FMMISOFEM_H__

#include "ISOFEM.h"


class FMMISOFEM : public ISOFEM
{
public:
	FMMISOFEM(std::string _filename, double _p1, double R): ISOFEM(_filename, _p1, R) {};
	
	FMMISOFEM(std::string _filename, double _p1, double R, int _HighOrder): ISOFEM(_filename, _p1, R, _HighOrder) {};
	
	void FMMStaticAnalysis(std::function<Eigen::Vector2d(Eigen::RowVector2d &)>, 
						std::function<double(const double &, const double &)>);
	
private:	
	
	void ConstructStiffMatr(std::vector<Triplet> &, std::vector<Eigen::MatrixXd> &, std::vector<Eigen::MatrixXd> &, Eigen::MatrixXd &, NdxArray &, JArray &, JacArray &);
	void ConstructStiffMatr(std::vector<Triplet> &, const std::vector<Eigen::MatrixXd> &, const std::vector<Eigen::MatrixXd> &, const Eigen::MatrixXd &, const NdxArray &, const JArray &, const JacArray &, std::function<double(const double &, const double &)>);
	void WriteToVTK(std::string _filename) override;
	
};

#endif

