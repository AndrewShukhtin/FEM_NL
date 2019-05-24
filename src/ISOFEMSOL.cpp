#include "ISOFEMSOL.h"


ISOFEMSOL::ISOFEMSOL(std::string _filename)
{
	MESH.ReadMesh(_filename);
	
	NumberOfElements = MESH.NumberOfElements();
	
	NumberOfNodes = MESH.NumberOfNodes();
	
	p1 = 1.0;
	p2 = 1.0 - p1;
	
	auto &Bounds = MESH.Bounds();
	
	HighOrder = UxBoundNumber = UyBoundNumber = UxUyBoundNumber = LoadBoundNumber = -1;
	
	std::string Ux("Ux");
	std::string Uy("Uy");	
	std::string UxUy("UxUy");
	std::string Load("Load");
	
	for (int i = 0; i < MESH.NumberOfBounds(); ++i)
	{
		if(!Ux.compare(Bounds[i].tag)) UxBoundNumber = i;
		if(!Uy.compare(Bounds[i].tag)) UyBoundNumber = i;
		if(!UxUy.compare(Bounds[i].tag)) UxUyBoundNumber = i;
		if(!Load.compare(Bounds[i].tag)) LoadBoundNumber = i;
	}
	
	double E = 210000.0, nu = 0.3;
	
	D << 1.0,   nu,                       0.0,
				nu,      1.0,             0.0,
				0.0,     0.0,  0.5*(1.0 - nu);
							
	D*=E/(1.0 - nu*nu);

};




ISOFEMSOL::ISOFEMSOL(std::string _filename, double _p1, double _R):p1(_p1), L(_R)
{
	MESH.ReadMesh(_filename);
	
	NumberOfElements = MESH.NumberOfElements();
	
	NumberOfNodes = MESH.NumberOfNodes();
	
	p2 = 1.0 - p1;
	
	auto &Bounds = MESH.Bounds();
	
	HighOrder = UxBoundNumber = UyBoundNumber = UxUyBoundNumber = LoadBoundNumber = -1;
	
	std::string Ux("Ux");
	std::string Uy("Uy");	
	std::string UxUy("UxUy");
	std::string Load("Load");
	
	for (int i = 0; i < MESH.NumberOfBounds(); ++i)
	{
		if(!Ux.compare(Bounds[i].tag)) UxBoundNumber = i;
		if(!Uy.compare(Bounds[i].tag)) UyBoundNumber = i;
		if(!UxUy.compare(Bounds[i].tag)) UxUyBoundNumber = i;
		if(!Load.compare(Bounds[i].tag)) LoadBoundNumber = i;
	}
	
	double E = 210000.0, nu = 0.3;
	
	D << 1.0,   nu,                       0.0,
				nu,      1.0,             0.0,
				0.0,     0.0,  0.5*(1.0 - nu);
							
	D*=E/(1.0 - nu*nu);	
}


ISOFEMSOL::ISOFEMSOL(std::string _filename, double _p1, double _R, int _HighOrder):p1(_p1), L(_R), HighOrder( _HighOrder)
{
	MESH.ReadMesh(_filename);
	
	NumberOfElements = MESH.NumberOfElements();
	
	NumberOfNodes = MESH.NumberOfNodes();
	
	p2 = 1.0 - p1;
	
	auto &Bounds = MESH.Bounds();
	
	UxBoundNumber = UyBoundNumber = UxUyBoundNumber = LoadBoundNumber = -1;
	
	std::string Ux("Ux");
	std::string Uy("Uy");	
	std::string UxUy("UxUy");
	std::string Load("Load");
	
	for (int i = 0; i < MESH.NumberOfBounds(); ++i)
	{
		if(!Ux.compare(Bounds[i].tag)) UxBoundNumber = i;
		if(!Uy.compare(Bounds[i].tag)) UyBoundNumber = i;
		if(!UxUy.compare(Bounds[i].tag)) UxUyBoundNumber = i;
		if(!Load.compare(Bounds[i].tag)) LoadBoundNumber = i;
	}
	
	double E = 210000.0, nu = 0.3;
	
	D << 1.0,   nu,                       0.0,
				nu,      1.0,             0.0,
				0.0,     0.0,  0.5*(1.0 - nu);
							
	D*=E/(1.0 - nu*nu);	
}



void ISOFEMSOL::Static_Analysis(std::function<Eigen::Vector2d(Eigen::RowVector2d &)> load)
{
	auto T1 = std::chrono::high_resolution_clock::now();
	
	auto t1 = std::chrono::high_resolution_clock::now();
	ConstructRightPart(load);
// 	std::cout << "\nF:\n" << F << std::endl;
	auto t2 = std::chrono::high_resolution_clock::now();
	std::cout << "ConstructRightPart time = "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
              << " milliseconds\n";
	
	t1 = std::chrono::high_resolution_clock::now();	
	
	{
	std::vector<Triplet> TripletList;
	
	ConstructStiffMatr(TripletList);
	t2 = std::chrono::high_resolution_clock::now();
	std::cout << "ConstructStiffMatr time = "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
              << " milliseconds\n";
	
	TripletList.shrink_to_fit();
	
	t1 = std::chrono::high_resolution_clock::now();
	K.setFromTriplets(TripletList.begin(), TripletList.end());
	t2 = std::chrono::high_resolution_clock::now();
	std::cout << "SetFromTriplets time = "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
              << " milliseconds\n";
	}
	//std::cout << "\nK\n:" << K << std::endl;
	
	t1 = std::chrono::high_resolution_clock::now();
	ApplyingConstraints();
	t2 = std::chrono::high_resolution_clock::now();
	std::cout << "ApplyingConstraints time = "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
              << " milliseconds\n";
			  
	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower> Solver;
	
	//std::cout << "\nK\n:" << K << std::endl;
	
	Solver.compute(K);	
	
	t1 = std::chrono::high_resolution_clock::now();
	U = Solver.solve(F);
	
	t2 = std::chrono::high_resolution_clock::now();
	std::cout << "System solve time = "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
              << " milliseconds\n";
	
	
	t1 = std::chrono::high_resolution_clock::now();
	cN = Counter(MESH.Elements());
	ComputeStrains();
	t2 = std::chrono::high_resolution_clock::now();
	std::cout << "ComputeStrains time = "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
              << " milliseconds\n";
	
	t1 = std::chrono::high_resolution_clock::now();
	ComputeStress();
	t2 = std::chrono::high_resolution_clock::now();
	std::cout << "ComputeStress time = "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
              << " milliseconds\n";
	
	t1 = std::chrono::high_resolution_clock::now();
	WriteToVTK(MESH.FileName());
	t2 = std::chrono::high_resolution_clock::now();
	std::cout << "WriteToVTK time = "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
              << " milliseconds\n";
			  
	auto T2 = std::chrono::high_resolution_clock::now();
	
	std::cout << "Total time = "
              << std::chrono::duration_cast<std::chrono::seconds>(T2-T1).count()
              << " seconds\n\n";

};

void ISOFEMSOL::NonLocal_Static_Analysis(std::function<Eigen::Vector2d(Eigen::RowVector2d &)> load, 
		std::function<double(const double &, const double &)> phi)
{
	auto T1 = std::chrono::high_resolution_clock::now();
	
	auto t1 = std::chrono::high_resolution_clock::now();
	ConstructRightPart(load);
	auto t2 = std::chrono::high_resolution_clock::now();
	std::cout << "ConstructRightPart time = "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
              << " milliseconds\n";
	
	t1 = std::chrono::high_resolution_clock::now();
	
	NdxArray NdxArr;  NdxArr.resize(MESH.NumberOfElements());
	JArray JArr;      JArr.resize(MESH.NumberOfElements());
	JacArray JacArr;  JacArr.resize(MESH.NumberOfElements());
	
	{
	std::vector<Triplet> TripletList;
	
	ConstructStiffMatr(TripletList, NdxArr, JArr, JacArr);
	
	ConstructStiffMatr(TripletList, NdxArr, JArr, JacArr,phi);
	t2 = std::chrono::high_resolution_clock::now();
	std::cout << "ConstructStiffMatr time = "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
              << " milliseconds\n";
	
	TripletList.shrink_to_fit();
	
	t1 = std::chrono::high_resolution_clock::now();
	K.setFromTriplets(TripletList.begin(), TripletList.end());
	
	t2 = std::chrono::high_resolution_clock::now();
	std::cout << "SetFromTriplets time = "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
              << " milliseconds\n";
	}
	t1 = std::chrono::high_resolution_clock::now();
	ApplyingConstraints();
	t2 = std::chrono::high_resolution_clock::now();
	std::cout << "ApplyingConstraints time = "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
              << " milliseconds\n";
			  
	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,Eigen::Lower> Solver;
	
	Solver.compute(K);	
	
	t1 = std::chrono::high_resolution_clock::now();
	U = Solver.solve(F);
	t2 = std::chrono::high_resolution_clock::now();
	std::cout << "System solve time = "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
              << " milliseconds\n";
	
// 	std::cout <<"\nU:\n"<< U <<std::endl;
	
	t1 = std::chrono::high_resolution_clock::now();
	cN = Counter(MESH.Elements());
	ComputeStrains();
	t2 = std::chrono::high_resolution_clock::now();
	std::cout << "ComputeStrains time = "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
              << " milliseconds\n";
	
	t1 = std::chrono::high_resolution_clock::now();
	ComputeStress(JacArr,phi);
	t2 = std::chrono::high_resolution_clock::now();
	std::cout << "ComputeStress time = "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
              << " milliseconds\n";
	
	t1 = std::chrono::high_resolution_clock::now();
	WriteToVTK(MESH.FileName());
	t2 = std::chrono::high_resolution_clock::now();
	std::cout << "WriteToVTK time = "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
              << " milliseconds\n";	
	
	auto T2 = std::chrono::high_resolution_clock::now();
	
	std::cout << "Total time = "
              << std::chrono::duration_cast<std::chrono::seconds>(T2-T1).count()
              << " seconds\n\n";
}


void ISOFEMSOL::ConstructStiffMatr(std::vector<Triplet> &TripletList)
{
	
	std::shared_ptr<FINITE_ELEMENT::FiniteElement> pFE;
	
	unsigned Order;
	
	switch (MESH.NodesPerElement())
	{
		case 4: 
		{
			pFE.reset(new FINITE_ELEMENT::BilinearElement);
			Order = 3;
			break;
		}
		case 8:
		{
			pFE.reset(new FINITE_ELEMENT::QuadraticSerendipElement);
			Order = 5;
			break;
		}
	}
	
	NodesPerElement = pFE->NodesPerElement();
	
	NumberOfElements = MESH.NumberOfElements();
	
	
	const FINITE_ELEMENT::QuadratureUtils QuadUtil(pFE.get(),Order);
	
	
	K.resize(NumberOfNodes * 2, NumberOfNodes * 2);
	
	Ke.resize(NodesPerElement * 2, NodesPerElement * 2);
	
	
	Eigen::MatrixXd Ndx = Eigen::MatrixXd::Zero(2, NodesPerElement); // производные функции форм в глобаной ск
	
	Eigen::RowVectorXi ElementNodesNumbers(NodesPerElement);
	
	BT = Eigen::MatrixXd::Zero(NodesPerElement * 2, 3);
	
	B = Eigen::MatrixXd::Zero(3, NodesPerElement * 2);
	
	
	int NumberOfTriplets = 0;
	

	NumberOfTriplets = NumberOfElements * 4 * NodesPerElement * NodesPerElement;

	
	TripletList.reserve(NumberOfTriplets);
	
	Eigen::MatrixXd ElementNodesCoord(NodesPerElement, 2);
	
	const auto &Nodes = MESH.Nodes();
	const auto &Elements = MESH.Elements();
	
	std::vector<int> IdX(NodesPerElement * 2);
	
	for(int e = 0; e < NumberOfElements; ++e)
	{
		ElementNodesNumbers = Elements.row(e);
		
		Ke = Eigen::MatrixXd::Zero(NodesPerElement * 2, NodesPerElement * 2);
		
		for(int j = 0; j < NodesPerElement; ++j)
		{		
			ElementNodesCoord.row(j) = Nodes.row(ElementNodesNumbers(j));
			IdX[j * 2] = ElementNodesNumbers(j) * 2;
			IdX[j * 2 + 1] = ElementNodesNumbers(j) * 2 + 1; 
		}	
		
		ComputeKe(QuadUtil.NumberOfQP(), QuadUtil.Weights(), ElementNodesCoord, QuadUtil.NGradArr(), Ndx);
		
		AddKeToTriplet(TripletList,IdX);
			
	}

}


void ISOFEMSOL::ConstructStiffMatr(std::vector<Triplet> &TripletList, NdxArray &NdxArr, JArray &JArr, JacArray &JacArr)
{
	
	std::shared_ptr<FINITE_ELEMENT::FiniteElement> pFE;
	
	unsigned Order;
	
	switch (MESH.NodesPerElement())
	{
		case 4: 
		{
			pFE.reset(new FINITE_ELEMENT::BilinearElement);
			Order = 3;
			break;
		}
		case 8:
		{
			Order = 5;
			pFE.reset(new FINITE_ELEMENT::QuadraticSerendipElement);
			break;
		}
	}
	
	NodesPerElement = pFE->NodesPerElement();
	
	NumberOfElements = MESH.NumberOfElements();
	
	
	const FINITE_ELEMENT::QuadratureUtils QuadUtil(pFE.get(),Order);
	
	
	K.resize(NumberOfNodes * 2, NumberOfNodes * 2);
	
	Ke.resize(NodesPerElement * 2, NodesPerElement * 2);
	
	
	Eigen::MatrixXd Ndx = Eigen::MatrixXd::Zero(2, NodesPerElement); // производные функции форм в глобаной ск
	
	Eigen::RowVectorXi ElementNodesNumbers(NodesPerElement);
	
	BT = Eigen::MatrixXd::Zero(NodesPerElement * 2, 3);
	
	B = Eigen::MatrixXd::Zero(3, NodesPerElement * 2);
	
	
	int NumberOfTriplets = 0;
	
	NaiveRnnSearch(RnnArr,L);
		
	NumberOfTriplets = NumberOfElements * 4 * NodesPerElement * NodesPerElement;
		
	for (size_t ei = 0; ei < NumberOfElements; ++ei)
	{
		NumberOfTriplets += RnnArr[ei].size() * 4 * NodesPerElement * NodesPerElement;
	}

	const auto NumberOfQP = QuadUtil.NumberOfQP();
	
	TripletList.reserve(NumberOfTriplets);
	
	Eigen::MatrixXd ElementNodesCoord(NodesPerElement, 2);
	
	const auto &Nodes = MESH.Nodes();
	const auto &Elements = MESH.Elements();
	
	const auto &Narr = QuadUtil.Narr();
	const auto &NGradArr = QuadUtil.NGradArr();
	
	const auto &Weights = QuadUtil.Weights();
	
	double Jac = 0.0;
	

	std::vector<int> IdX(NodesPerElement * 2);
	
	for(size_t e = 0; e < NumberOfElements; ++e)
	{
		ElementNodesNumbers = Elements.row(e);
		
		Ke = Eigen::MatrixXd::Zero(NodesPerElement * 2, NodesPerElement * 2);
		
		NdxArr[e].resize(NumberOfQP);
		JArr[e].resize(NumberOfQP);
		JacArr[e].resize(NumberOfQP);
		
		for(size_t j = 0; j < NodesPerElement; ++j)
		{		
			ElementNodesCoord.row(j) = Nodes.row(ElementNodesNumbers(j));
			IdX[j * 2] = ElementNodesNumbers(j) * 2;
			IdX[j * 2 + 1] = ElementNodesNumbers(j) * 2 + 1; 
		}	
		

		for (size_t q = 0; q < NumberOfQP; ++q)
		{
			Jmatr = NGradArr[q] * ElementNodesCoord;
			
			JArr[e][q] = Jmatr;
			
			Ndx = Jmatr.inverse() * NGradArr[q];
			
			NdxArr[e][q] = Ndx;
			
			Jac = Jmatr.determinant();
			
			JacArr[e][q] = Jac;
			
			ComputeB(Ndx);
		
			Ke += Weights(q) * BT * D * B * Jac;
		}
	
		Ke*=p1;

		for(size_t j = 0; j < NodesPerElement * 2; ++j)		
			for(size_t k = 0; k < NodesPerElement * 2; ++k)
			{
				if(IdX[k] <= IdX[j])
				TripletList.push_back(Triplet(IdX[j], IdX[k], Ke(j, k)));
			}
	}
}


void ISOFEMSOL::ConstructStiffMatr(std::vector<Triplet> &TripletList, const NdxArray &NdxArr, const JArray &JArr, const JacArray &JacArr, 
								   std::function<double(const double &, const double &)> phi)
{
	std::shared_ptr<FINITE_ELEMENT::FiniteElement> pFE;
	
	unsigned Order;
	switch (MESH.NodesPerElement())
	{
		case 4: 
		{
			pFE.reset(new FINITE_ELEMENT::BilinearElement);
			Order = 3;
			break;
		}
		case 8:
		{
			pFE.reset(new FINITE_ELEMENT::QuadraticSerendipElement);
			Order = 5;
			break;
		}
	}
	
	if(HighOrder!=-1) Order = HighOrder;
	
	NodesPerElement = pFE->NodesPerElement();
	
	NumberOfElements = MESH.NumberOfElements();
	
	const FINITE_ELEMENT::QuadratureUtils QuadUtil(pFE.get(), Order);
	
	Eigen::MatrixXd Ndxi = Eigen::MatrixXd::Zero(2, NodesPerElement); // производные функции форм в глобаной ск
	Eigen::MatrixXd Ndxj = Eigen::MatrixXd::Zero(2, NodesPerElement); // производные функции форм в глобаной ск
	
	Eigen::RowVectorXi ElementNodesNumbersi(NodesPerElement);
	Eigen::RowVectorXi ElementNodesNumbersj(NodesPerElement);
	
	BT = Eigen::MatrixXd::Zero(NodesPerElement * 2, 3);
	
	B = Eigen::MatrixXd::Zero(3, NodesPerElement * 2);
	
	Eigen::MatrixXd ElementNodesCoordi(NodesPerElement, 2);
	Eigen::MatrixXd ElementNodesCoordj(NodesPerElement, 2);
	
	Eigen::RowVector2d X, Y, R;
	
	double r = 0.0;
	
	const auto &Nodes = MESH.Nodes();
	const auto &Elements = MESH.Elements();
	
	const auto &Narr = QuadUtil.Narr();
	const auto &NGradArr = QuadUtil.NGradArr();
	
	const auto &Weights = QuadUtil.Weights();
	
	const auto NumberOfQP = QuadUtil.NumberOfQP();
	
	std::vector<Eigen::MatrixXd> NdxiArr(NumberOfQP), NdxjArr(NumberOfQP);
	
	for(size_t i = 0; i < NumberOfQP; ++i)
	{
		NdxiArr[i] = Eigen::MatrixXd::Zero(2, NodesPerElement);
		NdxjArr[i] = Eigen::MatrixXd::Zero(2, NodesPerElement);
	}
		
	std::vector<Eigen::Matrix2d> JArri, JArrj; 
	
	JArri.resize(NumberOfQP); JArrj.resize(NumberOfQP);
	
	std::vector<double> Jaci(NumberOfQP), Jacj(NumberOfQP);
	
	std::vector<int> IdX1(NodesPerElement * 2), IdX2(NodesPerElement * 2);

	for(int ei = 0; ei < NumberOfElements; ++ei)
	{
		ElementNodesNumbersi = Elements.row(ei);

		for(int j = 0; j < NodesPerElement; ++j)
		{		
			ElementNodesCoordi.row(j) = Nodes.row(ElementNodesNumbersi(j));
			IdX1[j * 2] = ElementNodesNumbersi(j) * 2;
			IdX1[j * 2 + 1] = ElementNodesNumbersi(j) * 2 + 1;
		}
		
		JArri = JArr[ei];
		
		NdxiArr = NdxArr[ei];
		
		Jaci = JacArr[ei];
		
		for(int ej = 0; ej < RnnArr[ei].size(); ++ej)
		{	
			Ke = Eigen::MatrixXd::Zero(NodesPerElement * 2,NodesPerElement * 2);
			
			int Ej = RnnArr[ei][ej];
			
			ElementNodesNumbersj = Elements.row(Ej);
			
			JArrj = JArr[Ej];
		
			NdxjArr = NdxArr[Ej];
		
			Jacj = JacArr[Ej];
			
			for(int j = 0; j < NodesPerElement; ++j)
			{
				ElementNodesCoordj.row(j) = Nodes.row(ElementNodesNumbersj(j));
				IdX2[j * 2] = ElementNodesNumbersj(j) * 2;
				IdX2[j * 2 + 1] = ElementNodesNumbersj(j) * 2 + 1;
			}
			
			for(int qi = 0; qi < NumberOfQP; ++qi)
			{	
				
				Ndxi = NdxiArr[qi];
				
				for(int k = 0; k < NodesPerElement; ++k)
				{
					BT(k * 2, 0) = Ndxi(0,k);
					BT(k * 2 + 1, 1) = Ndxi(1,k);
					BT(k * 2, 2) = Ndxi(1,k);    BT(k * 2 + 1, 2) = Ndxi(0,k);
				}
				
				X = Narr[qi] * ElementNodesCoordi;
				
				for(int qj = 0; qj < QuadUtil.NumberOfQP(); ++qj)
				{
					Ndxj = NdxjArr[qj];
					

					for(int k = 0; k < NodesPerElement; ++k)
					{
						B(0, k * 2) = Ndxj(0,k);
						B(1,k * 2 + 1) = Ndxj(1,k);
						B(2,k * 2) = Ndxj(1,k);    B(2,k * 2 + 1) = Ndxj(0,k);
					}					
					
					
					Y = Narr[qj] * ElementNodesCoordj;
					
					R = X - Y;
					
					r = R.norm();
					
					Ke += Weights(qi) * Weights(qj) * BT * D * B * phi(r,L) *  Jaci[qi] * Jacj[qj];
				}
			}	
			
			Ke *=p2;
				
			for(int j = 0; j < NodesPerElement * 2; ++j)		
				for(int k = 0; k < NodesPerElement * 2; ++k)
				{
					if(IdX2[k] <= IdX1[j])
					TripletList.push_back(Triplet(IdX1[j], IdX2[k], Ke(j ,k)));
				}
		
		}	
			

	}
	
}


void ISOFEMSOL::NaiveRnnSearch(std::vector<std::vector<int>> &RnnArr, const double &L)
{
	const auto &Elements = MESH.Elements();
	const auto &Nodes = MESH.Nodes();
	
	Eigen::MatrixXd COE(MESH.NumberOfElements(),2);
	
	Eigen::MatrixXd ElementNodesCoord(4, 2);
	
	Eigen::RowVectorXi ElementNodesNumbers(4);
	
	for (size_t e = 0; e < MESH.NumberOfElements(); ++e)
	{
		for(size_t i = 0; i < 4; ++i)
			ElementNodesNumbers(i) = Elements(e,i);
		
		for(size_t i = 0; i < 4; ++i)
			ElementNodesCoord.row(i) = Nodes.row(ElementNodesNumbers(i));
		
		COE.row(e) = 0.25 * (ElementNodesCoord.row(0) + ElementNodesCoord.row(1) + ElementNodesCoord.row(2) + ElementNodesCoord.row(3));	
	}
	
	RnnArr.resize(MESH.NumberOfElements());
	
	double dist = 0.0;
	Eigen::RowVector2d distV;
	
	for(size_t ei = 0; ei < MESH.NumberOfElements(); ++ei)
	{
		RnnArr[ei].reserve(MESH.NumberOfElements() /2);
		
		for(size_t ej = 0; ej < MESH.NumberOfElements(); ++ej)
		{
			distV = COE.row(ei) - COE.row(ej);
			dist = distV.norm();
			
			if (dist <= 1.05 * L)
				RnnArr[ei].push_back(ej);
		}
		RnnArr[ei].shrink_to_fit();
	}
	


}


void ISOFEMSOL::ConstructRightPart(std::function<Eigen::Vector2d(Eigen::RowVector2d &)> load)
{
	std::shared_ptr<FINITE_ELEMENT::FiniteElement> pFE;
	
	auto &Bounds = MESH.Bounds();
	
	auto &LoadBound = Bounds[LoadBoundNumber];
	
	unsigned Order;
	switch (LoadBound.Elements.cols())
	{
		case 2:
		{
			pFE.reset(new FINITE_ELEMENT::LinearElement);
			Order = 1;
			break;
		}
		case 3:
		{
			pFE.reset(new FINITE_ELEMENT::QuadraticElement);
			Order = 3;
			break;
		}
	}
	
	NodesPerElement = pFE->NodesPerElement();
	
	NumberOfElements = LoadBound.Elements.rows();
	
	const FINITE_ELEMENT::QuadratureUtils QuadUtil(pFE.get(),Order);
	
	Fe.resize(NodesPerElement * 2);
	
	F.resize(NumberOfNodes * 2);
	
	F = Eigen::VectorXd::Zero(NumberOfNodes * 2);
	
	Eigen::RowVectorXi ElementNodesNumbers(NodesPerElement);
	
	Eigen::MatrixXd ElementNodesCoord(NodesPerElement,2);
	
	Eigen::MatrixXd N = Eigen::MatrixXd::Zero(2 * NodesPerElement, 2);
	
	Eigen::MatrixXd Pn(2, NodesPerElement);
	
	Eigen::Vector2d P;
	
	Eigen::RowVector2d tmp;
	
	const auto &Weights = QuadUtil.Weights();
	const auto &Narr = QuadUtil.Narr();
	const auto &NGradArr = QuadUtil.NGradArr();
	
	
	for (int e = 0; e < NumberOfElements; ++e)		
	{
		Fe = Eigen::VectorXd::Zero(NodesPerElement * 2);
		
		ElementNodesNumbers = LoadBound.Elements.row(e);
		
		for(int j = 0; j <  NodesPerElement; ++j)
		{
			ElementNodesCoord.row(j) = MESH.Nodes().row(ElementNodesNumbers(j));
			
			tmp = ElementNodesCoord.row(j);
			
			Pn.col(j) = load(tmp);
			
		}		
		
		
		for (int j = 0; j < QuadUtil.NumberOfQP(); ++j)
		{
			Jvec = NGradArr[j] * ElementNodesCoord;
			
			P =  Pn*Narr[j].transpose();
		
			for(int k = 0; k < NodesPerElement; ++k)
			{
				N(2 * k, 0) = Narr[j](k);
				N(2 * k + 1, 1) = Narr[j](k);
			}
		

			Fe += Weights(j) * N * P * Jvec.norm();			
			
		}
		
		
		for(int j = 0; j < NodesPerElement; ++j)
		{
			F(ElementNodesNumbers(j) * 2) += Fe(2 * j);
			F(ElementNodesNumbers(j) * 2 + 1) += Fe(2 * j + 1);
		}
	}
	
}



void ISOFEMSOL::ApplyingConstraints()
{
	auto &Bounds = MESH.Bounds();
	std::vector<int> constraintIdx;
	constraintIdx.reserve(1000);
	
// 	printf("\nConstructRightPart\nNumberOfNodes = %d, NumberOfElements = %d, NodesPerElement = %d, UxBoundNumber = %d, \nUyBoundNumber = %d, UxUyBoundNumber = %d, LoadBoundNumber = %d\n", NumberOfNodes, NumberOfElements, NodesPerElement, UxBoundNumber,UyBoundNumber, UxUyBoundNumber, LoadBoundNumber );
	
	
	if(UxUyBoundNumber !=-1)
	{
		auto &fixedBound = Bounds[UxUyBoundNumber];
	
		std::vector<int> fixedNodes = UniqueNodes(fixedBound.Elements);
	
		
		
		//std::cout << "fixedBound.Elements:\n" <<fixedBound.Elements <<std::endl;
	
		for(int i = 0; i < fixedNodes.size(); ++i)
		{
			constraintIdx.push_back(2*fixedNodes[i]);
			constraintIdx.push_back(2*fixedNodes[i] + 1);
		}
		

	} 
	
	if(UxBoundNumber !=-1)
	{
		auto &fixedBound = Bounds[UxBoundNumber];
	
		std::vector<int> fixedNodes = UniqueNodes(fixedBound.Elements);
	
		//constraintIdx.reserve(fixedNodes.size() * 2);
	
	
		for(size_t i = 0; i < fixedNodes.size(); ++i)
		{
			constraintIdx.push_back(2*fixedNodes[i]);
		}		
				
	}
	
	if(UyBoundNumber !=-1)
	{
		auto &fixedBound = Bounds[UyBoundNumber];
	
		std::vector<int> fixedNodes = UniqueNodes(fixedBound.Elements);
	
		//constraintIdx.reserve(fixedNodes.size() * 2);
	
	
		for(size_t i = 0; i < fixedNodes.size(); ++i)
		{
			constraintIdx.push_back(2*fixedNodes[i] + 1);
		}		
				
		
	}
	
	constraintIdx.shrink_to_fit();
		
	for(size_t i = 0; i < constraintIdx.size(); ++i)
		F(constraintIdx[i]) = 0.0;
	
	for(size_t i = 0; i < K.outerSize(); ++i)
		for(Eigen::SparseMatrix<double>::InnerIterator it(K,i); it; ++it)
			for(std::vector<int>::iterator idit = constraintIdx.begin(); idit!= constraintIdx.end(); ++idit)
			{
				if (it.row() == *idit || it.col() == *idit )
				{
					it.valueRef() = it.row() > it.col() ? 0.0 : 0.0;
				}
			}
}

std::vector<int> ISOFEMSOL::UniqueNodes(const Eigen::MatrixXi &Elements)
{
	std::vector<int> unique_arr; unique_arr.reserve(Elements.size());

	for (int i = 0; i < Elements.rows(); ++i)
		for (int j = 0; j < Elements.cols(); ++j)
			unique_arr.push_back(Elements(i,j));

	std::sort(unique_arr.begin(), unique_arr.end());

	auto last = std::unique(unique_arr.begin(), unique_arr.end());

	unique_arr.erase(last, unique_arr.end());
	
	return unique_arr;
}

inline void ISOFEMSOL::ComputeB(const Eigen::MatrixXd &Ndx)
{
	for(int k = 0; k < NodesPerElement; ++k)
		{
			B(0, k * 2) = Ndx(0,k);
			B(1,k * 2 + 1) = Ndx(1,k);
			B(2,k * 2) = Ndx(1,k);    B(2, k * 2 + 1) = Ndx(0,k);
			
			BT(k * 2,0) = Ndx(0,k);
			BT(k * 2 + 1,1) = Ndx(1,k);
			BT(k * 2,2) = Ndx(1,k);    BT(k * 2 + 1, 2) = Ndx(0,k);
		}
}


inline void ISOFEMSOL::ComputeKe(const int &NumberOfQP, const Eigen::RowVectorXd &Weights, const Eigen::MatrixXd &ElementNodesCoord, const std::vector<Eigen::MatrixXd> &NGradArr, Eigen::MatrixXd &Ndx)
{
	for (int j = 0; j < NumberOfQP; ++j)
	{
		Jmatr = NGradArr[j] * ElementNodesCoord;
		
		Ndx = Jmatr.inverse() * NGradArr[j];
			
		ComputeB(Ndx);
		
		Ke += Weights(j) * BT * D * B * Jmatr.determinant();
	}
	
	Ke*=p1;
}

inline void ISOFEMSOL::AddKeToTriplet(std::vector<Triplet> &TripletList, const std::vector<int> &IdX)
{
	for(int j = 0; j < NodesPerElement * 2; ++j)		
		for(int k = 0; k < NodesPerElement * 2; ++k)
		{
			if(IdX[k] <= IdX[j])
			TripletList.push_back(Triplet(IdX[j], IdX[k], Ke(j, k)));
		}

}

void ISOFEMSOL::ComputeStrains()
{
	
	EpsiXX.resize(NumberOfNodes); EpsiXY.resize(NumberOfNodes); EpsiYY.resize(NumberOfNodes); 
	
	
	EpsiXX =  Eigen::VectorXd::Zero(NumberOfNodes);
	EpsiXY =  Eigen::VectorXd::Zero(NumberOfNodes);
	EpsiYY =  Eigen::VectorXd::Zero(NumberOfNodes);
	
	
	Eigen::VectorXd EpsiXX_Element = Eigen::VectorXd::Zero(NodesPerElement);

	Eigen::VectorXd EpsiYY_Element = Eigen::VectorXd::Zero(NodesPerElement);	
	
	Eigen::VectorXd EpsiXY_Element = Eigen::VectorXd::Zero(NodesPerElement);	
	
	
	
	Eigen::VectorXd EpsiXX_QP = Eigen::VectorXd::Zero(NodesPerElement);

	Eigen::VectorXd EpsiYY_QP = Eigen::VectorXd::Zero(NodesPerElement);	
	
	Eigen::VectorXd EpsiXY_QP = Eigen::VectorXd::Zero(NodesPerElement);
	
	
	
	Eigen::VectorXd Ue =  Eigen::VectorXd::Zero(NodesPerElement * 2);
	
	Eigen::VectorXd Epsi =  Eigen::VectorXd::Zero(3);
	
	
	
	std::shared_ptr<FINITE_ELEMENT::FiniteElement> pFE;
	
	unsigned Order;
	
	switch (MESH.NodesPerElement())
	{
		case 4:
		{
			pFE.reset(new FINITE_ELEMENT::BilinearElement);
			Order = 3;
			break;
		}
		case 8:
		{
			pFE.reset(new FINITE_ELEMENT::QuadraticSerendipElement);
			Order = 5;
			break;
		}
	}
	
	const FINITE_ELEMENT::QuadratureUtils QuadUtil(pFE.get(), Order);
	
	Eigen::MatrixXd NQP(NodesPerElement,NodesPerElement);
	
	NodesPerElement = pFE->NodesPerElement();
	
	Eigen::MatrixXd Ndx(2, NodesPerElement);
	
	Eigen::MatrixXd ElementNodesCoord(NodesPerElement, 2);
	
	Eigen::RowVectorXi ElementNodesNumbers(NodesPerElement);
	
	
	
	
	const auto &Nodes = MESH.Nodes();
	const auto &Elements = MESH.Elements();
	
	const auto &Narr = QuadUtil.Narr();
	const auto &NGradArr = QuadUtil.NGradArr();
	
	const auto &Weights = QuadUtil.Weights();
	
	
	for(int q = 0; q < NodesPerElement; ++q)
		NQP.row(q) =  Narr[q];
	

	
	Eigen::MatrixXd b = Eigen::MatrixXd::Identity(NodesPerElement,NodesPerElement);
	
	Eigen::MatrixXd NQPI = NQP.fullPivLu().solve(b);
	
	
	
	for(int e = 0; e < NumberOfElements; ++e)
	{		
		B = Eigen::MatrixXd::Zero(3, NodesPerElement * 2);
		
		ElementNodesNumbers = Elements.row(e);
		
		
		for(int j = 0; j < NodesPerElement; ++j)
		{
			Ue(j * 2) = U(ElementNodesNumbers(j) * 2);
			Ue(j * 2 + 1) = U(ElementNodesNumbers(j) * 2 + 1);
		}
			
		
		for(int j = 0; j < NodesPerElement; ++j)		
			ElementNodesCoord.row(j) = Nodes.row(ElementNodesNumbers(j));
		
		
		for (int j = 0; j < NodesPerElement; ++j)
		{
			Jmatr = NGradArr[j] * ElementNodesCoord;
			
			Ndx = Jmatr.inverse() * NGradArr[j];
			
			for(int k = 0; k < NodesPerElement; ++k)
			{
				B(0, k * 2) = Ndx(0, k);
				B(1, k * 2 + 1) = Ndx(1, k);
				B(2, k * 2) = Ndx(1, k);    B(2, k * 2 + 1) = Ndx(0, k);
			}
			
			Epsi =  B * Ue;
			
			EpsiXX_QP(j) = Epsi(0);
			EpsiYY_QP(j) = Epsi(1);
			EpsiXY_QP(j) = Epsi(2);
			
			
		}	
		
		EpsiXX_Element = NQPI * EpsiXX_QP;
		EpsiYY_Element = NQPI * EpsiYY_QP;
		EpsiXY_Element = 0.5 * NQPI * EpsiXY_QP;
		
		for (int j = 0; j < NodesPerElement; ++j)
		{
			
			EpsiXX(ElementNodesNumbers(j)) += EpsiXX_Element(j);
			EpsiYY(ElementNodesNumbers(j)) += EpsiYY_Element(j);
			EpsiXY(ElementNodesNumbers(j)) += EpsiXY_Element(j);
			
		}	
		
		
	}	
	
	
	for (auto count = cN.begin(); count != cN.end(); ++count) 
	{
		
		EpsiXX(count->first)/=count->second;
		EpsiXY(count->first)/=count->second;
		EpsiYY(count->first)/=count->second;
		
	}

}



std::map<int, int> ISOFEMSOL::Counter(const Eigen::MatrixXi &Elements)
{
	std::vector<int> vals; vals.reserve(NodesPerElement*NumberOfElements);
	
	for(size_t i = 0; i < NumberOfElements; ++i)
		for(size_t j = 0; j < NodesPerElement; ++j)
			vals.push_back(Elements(i,j));
	
	vals.shrink_to_fit();
	std::map<int, int> rv;

    for (auto val = vals.begin(); val != vals.end(); ++val) 
		rv[*val]++;

	return rv;
}



void ISOFEMSOL::ComputeStress()
{
	SigmaXX.resize(NumberOfNodes); SigmaYY.resize(NumberOfNodes); SigmaXY.resize(NumberOfNodes);
	
	SigmaXX =  Eigen::VectorXd::Zero(NumberOfNodes);
	SigmaXY =  Eigen::VectorXd::Zero(NumberOfNodes);
	SigmaYY =  Eigen::VectorXd::Zero(NumberOfNodes);
	
	Eigen::VectorXd Sigma =  Eigen::VectorXd::Zero(3);
	Eigen::VectorXd Epsi  =  Eigen::VectorXd::Zero(3);
	
	Eigen::RowVectorXi ElementNodesNumbers(NodesPerElement);
	
	const auto &Elements = MESH.Elements();
	
	for(int e = 0; e < NumberOfElements; ++e)
	{		
		ElementNodesNumbers = Elements.row(e);
		
		for (int j = 0; j < NodesPerElement; ++j)
		{
			
			Epsi(0) = EpsiXX(ElementNodesNumbers(j));
			Epsi(1) = EpsiYY(ElementNodesNumbers(j));
			Epsi(2) = EpsiXY(ElementNodesNumbers(j));
			
			Sigma = p1 * D * Epsi;
			
			SigmaXX(ElementNodesNumbers(j)) += Sigma(0);
			SigmaYY(ElementNodesNumbers(j)) += Sigma(1);
			SigmaXY(ElementNodesNumbers(j)) += 0.5*Sigma(2);
			
		}	
		
		
	}	
	
	
	for (auto count = cN.begin(); count != cN.end(); ++count) 
	{		
		SigmaXX(count->first)/=count->second;
		SigmaXY(count->first)/=count->second;
		SigmaYY(count->first)/=count->second;
	}
	
	auto &Bounds = MESH.Bounds();
	
	auto &LoadBound = Bounds[LoadBoundNumber];
	
}


void ISOFEMSOL::ComputeStress(const JacArray &JacArr ,std::function<double(const double&, const double &)> phi)
{
	SigmaXX.resize(NumberOfNodes); SigmaYY.resize(NumberOfNodes); SigmaXY.resize(NumberOfNodes);
	
	SigmaXX =  Eigen::VectorXd::Zero(NumberOfNodes);
	SigmaXY =  Eigen::VectorXd::Zero(NumberOfNodes);
	SigmaYY =  Eigen::VectorXd::Zero(NumberOfNodes);
	
	Eigen::VectorXd Sigma =  Eigen::VectorXd::Zero(3);
	Eigen::VectorXd Epsi  =  Eigen::VectorXd::Zero(3);
	
	std::shared_ptr<FINITE_ELEMENT::FiniteElement> pFE;
	
	unsigned Order;
	
	switch (MESH.NodesPerElement())
	{
		case 4: 
		{
			pFE.reset(new FINITE_ELEMENT::BilinearElement);
			Order = 3;
			break;
		}
		case 8:
		{
			pFE.reset(new FINITE_ELEMENT::QuadraticSerendipElement);
			Order = 5;
			break;
		}
	}
	
	if(HighOrder!=-1) Order = HighOrder;
	
	const FINITE_ELEMENT::QuadratureUtils QuadUtil(pFE.get(), Order);
	
	
	NodesPerElement = pFE->NodesPerElement();
	
	NumberOfElements = MESH.NumberOfElements();
	
	Eigen::MatrixXd NQP(NodesPerElement,NodesPerElement);
	
		
	Eigen::VectorXd EpsiXX_Element = Eigen::VectorXd::Zero(NodesPerElement);

	Eigen::VectorXd EpsiYY_Element = Eigen::VectorXd::Zero(NodesPerElement);	
	
	Eigen::VectorXd EpsiXY_Element = Eigen::VectorXd::Zero(NodesPerElement);	
	
	
	Eigen::VectorXd SigmaXX_Element = Eigen::VectorXd::Zero(NodesPerElement);

	Eigen::VectorXd SigmaYY_Element = Eigen::VectorXd::Zero(NodesPerElement);	
	
	Eigen::VectorXd SigmaXY_Element = Eigen::VectorXd::Zero(NodesPerElement);	
	
	
	Eigen::VectorXd SigmaXX_QP = Eigen::VectorXd::Zero(NodesPerElement);

	Eigen::VectorXd SigmaYY_QP = Eigen::VectorXd::Zero(NodesPerElement);	
	
	Eigen::VectorXd SigmaXY_QP = Eigen::VectorXd::Zero(NodesPerElement);

	
	Eigen::RowVectorXi ElementNodesNumbersi(NodesPerElement);
	Eigen::RowVectorXi ElementNodesNumbersj(NodesPerElement);
		
	Eigen::MatrixXd ElementNodesCoordi(NodesPerElement, 2);
	Eigen::MatrixXd ElementNodesCoordj(NodesPerElement, 2);
	
	Eigen::RowVector2d X, Y, R;
	
	double r;
	
	const auto &Nodes = MESH.Nodes();
	const auto &Elements = MESH.Elements();
	
	const auto &Narr = QuadUtil.Narr();
	const auto &NGradArr = QuadUtil.NGradArr();
	
	const auto &Weights = QuadUtil.Weights();
	

	for(int q = 0; q < NodesPerElement; ++q)
		NQP.row(q) = Weights(q) * Narr[q];

	
	Eigen::MatrixXd b = Eigen::MatrixXd::Identity(NodesPerElement,NodesPerElement);
	
	Eigen::MatrixXd NQPI = NQP.fullPivLu().solve(b);
	
	const auto NumberOfQP = QuadUtil.NumberOfQP();
	
	std::vector<double> Jac(NumberOfQP);
	
	for(int ei = 0; ei < NumberOfElements; ++ei)
	{		
		
		ElementNodesNumbersi = Elements.row(ei);
		
		for(int j = 0; j < NodesPerElement; ++j)		
			ElementNodesCoordi.row(j) = Nodes.row(ElementNodesNumbersi(j));
	
		
		for (int j = 0; j < NodesPerElement; ++j)
		{
			Epsi(0) = EpsiXX(ElementNodesNumbersi(j));
			Epsi(1) = EpsiYY(ElementNodesNumbersi(j));
			Epsi(2) = EpsiXY(ElementNodesNumbersi(j));
			
			Sigma = p1 * D * Epsi;
			
			SigmaXX(ElementNodesNumbersi(j)) += Sigma(0);
			SigmaYY(ElementNodesNumbersi(j)) += Sigma(1);
			SigmaXY(ElementNodesNumbersi(j)) += 0.5 * Sigma(2);
			
		}	
		
		SigmaXX_Element = Eigen::VectorXd::Zero(NodesPerElement);

		SigmaYY_Element = Eigen::VectorXd::Zero(NodesPerElement);	
	
		SigmaXY_Element = Eigen::VectorXd::Zero(NodesPerElement);
		
		for(int ej = 0; ej < RnnArr[ei].size(); ++ej)
		{
			ElementNodesNumbersj = Elements.row(RnnArr[ei][ej]);
			
			Jac = JacArr[RnnArr[ei][ej]];
			
			for(int j = 0; j < NodesPerElement; ++j)
				ElementNodesCoordj.row(j) = Nodes.row(ElementNodesNumbersj(j));
			
			for (int j = 0; j < NodesPerElement; ++j)
			{
				EpsiXX_Element(j) = EpsiXX(ElementNodesNumbersj(j));
				EpsiYY_Element(j) = EpsiYY(ElementNodesNumbersj(j));
				EpsiXY_Element(j) = EpsiXY(ElementNodesNumbersj(j));
			}
			
		

			for(int q = 0; q < NodesPerElement; ++q)
			{
				X = Narr[q] * ElementNodesCoordi;
				
				Y = Narr[q] * ElementNodesCoordj;
				
				R = X - Y;
				
				r = R.norm();
				
				Epsi(0) = Narr[q] * EpsiXX_Element;
				Epsi(1) = Narr[q] * EpsiYY_Element;
				Epsi(2) = Narr[q] * EpsiXY_Element;
				
				
				Sigma = Weights(q) * phi(r,L) * D * Epsi * Jac[q];
				
				SigmaXX_QP(q) = Sigma(0);
				SigmaYY_QP(q) = Sigma(1);
				SigmaXY_QP(q) = Sigma(2);
				
			}

		
			
			SigmaXX_Element += NQPI * SigmaXX_QP;
			
			SigmaYY_Element += NQPI * SigmaYY_QP;
			
			SigmaXY_Element += NQPI * SigmaXY_QP;
			
	}
			
		for (int j = 0; j < NodesPerElement; ++j)
		{
			
			SigmaXX(ElementNodesNumbersi(j)) += p2 * SigmaXX_Element(j);
			SigmaYY(ElementNodesNumbersi(j)) += p2 * SigmaYY_Element(j);
			SigmaXY(ElementNodesNumbersi(j)) += 0.5 * p2 * SigmaXY_Element(j);
			
		}	
		
	}	
	
	for (auto count = cN.begin(); count != cN.end(); ++count) 
	{		
		SigmaXX(count->first)/=count->second;
		SigmaXY(count->first)/=count->second;
		SigmaYY(count->first)/=count->second;		
	}
	
}


void ISOFEMSOL::WriteToVTK(std::string _filename)
{
	
	_filename.erase(_filename.end()-4, _filename.end());
	
	std::string path ("VTK/");
	
	path += _filename;
	
	std::string command = "mkdir -p " + path;
	
	system(command.c_str());
	
	_filename = "VTK/" + _filename  + "/"+_filename + "_" + std::to_string(p1) +".vtk";
	
	std::cout << _filename << std::endl;
	
	std::ofstream vtkfile;
	vtkfile.open(_filename);
	
	vtkfile << "# vtk DataFile Version 4.2" << std::endl;
	vtkfile << "Displacement, Stresses and Strains" << std::endl;
	vtkfile << "ASCII" << std::endl;
	vtkfile << "DATASET UNSTRUCTURED_GRID" << std::endl;
	vtkfile << "POINTS " << MESH.NumberOfNodes() <<" double" << std::endl;
	
	const auto &Nodes = MESH.Nodes();
	
	for(int i = 0; i < MESH.NumberOfNodes(); ++i)
	{
		vtkfile << Nodes(i,0) << " " << Nodes(i,1) << " " << 0.0 << std::endl;
	}
	
	const auto &Elements = MESH.Elements();
	
	if(MESH.NodesPerElement() == 4)
	{
		vtkfile << "CELLS " << MESH.NumberOfElements() << " " << MESH.NumberOfElements()*5 << std::endl;
		
		for(int i = 0; i < MESH.NumberOfElements(); ++i)
		{
			vtkfile << 4 << " ";
			vtkfile << Elements(i,0) << " " << Elements(i,1) << " " <<Elements(i,2) << " " << Elements(i,3) << std::endl;
		}
		
		vtkfile << "CELL_TYPES " << MESH.NumberOfElements() << std::endl;
		
		for(int i = 0; i < MESH.NumberOfElements(); ++i)
		{
			vtkfile << 9 << std::endl;
		}
	}else if(MESH.NodesPerElement() == 8)
	{
		vtkfile << "CELLS " << MESH.NumberOfElements() << " " << MESH.NumberOfElements()*9 << std::endl;
		
		for(int i = 0; i < MESH.NumberOfElements(); ++i)
		{
			vtkfile << 8 << " ";
// 			vtkfile << Elements(i,7) << " " << Elements(i,4) << " " <<Elements(i,5) << " " << Elements(i,6) << " ";
// 			vtkfile << Elements(i,3) << " " << Elements(i,0) << " " <<Elements(i,1) << " " << Elements(i,2) << std::endl;
			
			vtkfile << Elements(i,0) << " " << Elements(i,1) << " " <<Elements(i,2) << " " << Elements(i,3) << " ";
			vtkfile << Elements(i,4) << " " << Elements(i,5) << " " <<Elements(i,6) << " " << Elements(i,7) << std::endl;
		}
		
		vtkfile << "CELL_TYPES " << MESH.NumberOfElements() << std::endl;
		
		for(int i = 0; i < MESH.NumberOfElements(); ++i)
		{
			vtkfile << 23 << std::endl;
		}
				
	}
	
	vtkfile << "POINT_DATA " << MESH.NumberOfNodes() << std::endl;
	vtkfile << "SCALARS U_X double " << 1 << std::endl;
	vtkfile << "LOOKUP_TABLE default" << std::endl;
	
	for(int i = 0; i < MESH.NumberOfNodes(); ++i)
		vtkfile << U(2 * i) << std::endl;
	
	vtkfile << "SCALARS U_Y double " << 1 << std::endl;
	vtkfile << "LOOKUP_TABLE default" << std::endl;
	
	for(int i = 0; i < MESH.NumberOfNodes(); ++i)
		vtkfile << U(2 * i + 1) << std::endl;
	
	vtkfile << "SCALARS EpsiXX double " << 1 << std::endl;
	vtkfile << "LOOKUP_TABLE default" << std::endl;
	
	for(int i = 0; i < MESH.NumberOfNodes(); ++i)
		vtkfile << EpsiXX(i) << std::endl;
	
	vtkfile << "SCALARS EpsiYY double " << 1 << std::endl;
	vtkfile << "LOOKUP_TABLE default" << std::endl;
	
	for(int i = 0; i < MESH.NumberOfNodes(); ++i)
		vtkfile << EpsiYY(i) << std::endl;
	
	vtkfile << "SCALARS EpsiXY double " << 1 << std::endl;
	vtkfile << "LOOKUP_TABLE default" << std::endl;
	
	for(int i = 0; i < MESH.NumberOfNodes(); ++i)
		vtkfile << EpsiXY(i) << std::endl;
	
	
	vtkfile << "SCALARS SigmaXX double " << 1 << std::endl;
	vtkfile << "LOOKUP_TABLE default" << std::endl;
	
	for(int i = 0; i < MESH.NumberOfNodes(); ++i)
		vtkfile << SigmaXX(i) << std::endl;
	
	vtkfile << "SCALARS SigmaYY double " << 1 << std::endl;
	vtkfile << "LOOKUP_TABLE default" << std::endl;
	
	for(int i = 0; i < MESH.NumberOfNodes(); ++i)
		vtkfile << SigmaYY(i) << std::endl;
	
	vtkfile << "SCALARS SigmaXY double " << 1 << std::endl;
	vtkfile << "LOOKUP_TABLE default" << std::endl;
	
	for(int i = 0; i < MESH.NumberOfNodes(); ++i)
		vtkfile << SigmaXY(i) << std::endl;
	
	vtkfile.close();
}

