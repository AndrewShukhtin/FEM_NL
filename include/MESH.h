#ifndef __MESH_H__
#define __MESH_H__

#include "Finite_Element_Class.h"
#include "QuadratureUtils.h"
#include <fstream>
#include <string>
#include <iostream>

namespace FINITE_ELEMENT
{
	struct Bound
	{
		std::string tag;
		int type;
		Eigen::MatrixXi Elements;
	};
	
	class MESH
	{
	public:
		MESH();
		
		void ReadMesh(const std::string);
		
		const std::string &FileName() const;
		
		const Eigen::MatrixXd &Nodes() const;
		
		const Eigen::MatrixXi &Elements() const;
		
		const std::vector<Bound> &Bounds() const;
		
		const int NumberOfBounds();
		
		const int NodesPerElement();
		
		const int NumberOfNodes();
		
		const int NumberOfElements();
		
		~MESH();
	private:
		std::ifstream meshfile;
		std::string _filename;
		Eigen::MatrixXd _Nodes;
		Eigen::MatrixXi _Elements;
		std::vector<Bound> _Bounds;
	};

	std::ostream &operator <<(std::ostream &os, const std::vector<Bound> &Bounds);

	
}

#endif
