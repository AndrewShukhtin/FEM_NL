#include "FMMISOFEM.h"


void FMMISOFEM::FMMStaticAnalysis(std::function<Eigen::Vector2d(Eigen::RowVector2d &)> load, 
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
	{
		NdxArray NdxArr;  NdxArr.resize(MESH.NumberOfElements());
		JArray JArr;      JArr.resize(MESH.NumberOfElements());
		JacArray JacArr;  JacArr.resize(MESH.NumberOfElements());
		
		{
			std::vector<Triplet> TripletList;
			
			std::vector<Eigen::MatrixXd> KeArr; KeArr.resize(MESH.NumberOfElements());
			
			Eigen::MatrixXd COE(2,MESH.NumberOfElements());
			
			std::vector<Eigen::MatrixXd> GaussNodes;  GaussNodes.resize(MESH.NumberOfElements());
			
			ConstructStiffMatr(TripletList, KeArr, GaussNodes, COE, NdxArr, JArr, JacArr);
			
			ConstructStiffMatr(TripletList, KeArr, GaussNodes, COE, NdxArr, JArr, JacArr, phi);
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
	}
	
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





void FMMISOFEM::ConstructStiffMatr(std::vector<Triplet> &TripletList, std::vector<Eigen::MatrixXd> &KeArr, std::vector<Eigen::MatrixXd> &GaussNodes, Eigen::MatrixXd &COE , NdxArray &NdxArr, JArray &JArr, JacArray &JacArr)
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
	
	if(HighOrder!=-1) Order = HighOrder;
	
	const FINITE_ELEMENT::QuadratureUtils QuadUtilNL(pFE.get(),Order);
	
	
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
	const auto NumberOfQPNL = QuadUtilNL.NumberOfQP();
	
	TripletList.reserve(NumberOfTriplets);
	
	Eigen::MatrixXd ElementNodesCoord(NodesPerElement, 2);
	
	const auto &Nodes = MESH.Nodes();
	const auto &Elements = MESH.Elements();
	
	const auto &Narr = QuadUtil.Narr();
	const auto &NGradArr = QuadUtil.NGradArr();
	
	const auto &Weights = QuadUtil.Weights();
	
	
	const auto &NarrNL = QuadUtilNL.Narr();
	const auto &NGradArrNL = QuadUtilNL.NGradArr();
	
	const auto &WeightsNL = QuadUtilNL.Weights();
	
	Eigen::Vector2d X , Y , R;
	
	Eigen::MatrixXd ElementGaussNodes(2, NumberOfQPNL);
	
	double Jac = 0.0;

	std::vector<int> IdX(NodesPerElement * 2);
	
	for(size_t e = 0; e < NumberOfElements; ++e)
	{
		ElementNodesNumbers = Elements.row(e);
		
		Ke = Eigen::MatrixXd::Zero(NodesPerElement * 2, NodesPerElement * 2);
		
		NdxArr[e].resize(NumberOfQPNL);
		JArr[e].resize(NumberOfQPNL);
		JacArr[e].resize(NumberOfQPNL);
		
		
		GaussNodes[e].resize(2, NumberOfQPNL);
		
		for(size_t j = 0; j < NodesPerElement; ++j)
		{		
			ElementNodesCoord.row(j) = Nodes.row(ElementNodesNumbers(j));
			IdX[j * 2] = ElementNodesNumbers(j) * 2;
			IdX[j * 2 + 1] = ElementNodesNumbers(j) * 2 + 1; 
		}	
		
		COE(0,e) = 0.25 * (ElementNodesCoord(0, 0) + ElementNodesCoord(1, 0) + ElementNodesCoord(2, 0) + ElementNodesCoord(3, 0));
		COE(1,e) = 0.25 * (ElementNodesCoord(0, 1) + ElementNodesCoord(1, 1) + ElementNodesCoord(2, 1) + ElementNodesCoord(3, 1));
		
		
		for (size_t q = 0; q < NumberOfQP; ++q)
		{
			Jmatr = NGradArr[q] * ElementNodesCoord;
			
			Ndx = Jmatr.inverse() * NGradArr[q];

			Jac = Jmatr.determinant();
			
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
			
		Ke = Eigen::MatrixXd::Zero(NodesPerElement * 2, NodesPerElement * 2);
		
		Y = COE.col(e);
		
		for (size_t q = 0; q < NumberOfQPNL; ++q)
		{
			
			X = Narr[q] * ElementNodesCoord;
			
			ElementGaussNodes.col(q) = X;

			Jmatr = NGradArrNL[q] * ElementNodesCoord;
			
			JArr[e][q] = Jmatr;
			
			Ndx = Jmatr.inverse() * NGradArrNL[q];
			
			NdxArr[e][q] = Ndx;
			
			Jac = Jmatr.determinant();
			
			JacArr[e][q] = Jac;
			
			R = X - Y;
			
			ComputeB(Ndx);
		
			Ke += WeightsNL(q) * BT * D * B * Jac * R(0) * R(1) / (L * L);
		}
		
		GaussNodes[e] = ElementGaussNodes;
		
		KeArr[e] = Ke;

	}
}


void FMMISOFEM::ConstructStiffMatr(std::vector<Triplet> &TripletList, const std::vector<Eigen::MatrixXd> &KeArr, const std::vector<Eigen::MatrixXd> &GaussNodes, const Eigen::MatrixXd &COE, 
						const NdxArray &NdxArr, const JArray &JArr, const JacArray &JacArr, std::function<double(const double &, const double &)> phi)
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
	
	Eigen::MatrixXd Kej = Eigen::MatrixXd::Zero(NodesPerElement * 2,NodesPerElement * 2);
	
	Eigen::MatrixXd ElementGaussNodes(2,NumberOfQP);

	for(size_t ei = 0; ei < NumberOfElements; ++ei)
	{
		ElementNodesNumbersi = Elements.row(ei);

		for(size_t j = 0; j < NodesPerElement; ++j)
		{		
			ElementNodesCoordi.row(j) = Nodes.row(ElementNodesNumbersi(j));
			IdX1[j * 2] = ElementNodesNumbersi(j) * 2;
			IdX1[j * 2 + 1] = ElementNodesNumbersi(j) * 2 + 1;
		}
		
		for(size_t ej = 0; ej < RnnArr[ei].size(); ++ej)
		{	
			Ke = Eigen::MatrixXd::Zero(NodesPerElement * 2,NodesPerElement * 2);
			
			size_t Ej = RnnArr[ei][ej];
			
			ElementNodesNumbersj = Elements.row(Ej);
		
			Jacj = JacArr[Ej];
			
			for(int j = 0; j < NodesPerElement; ++j)
			{
				IdX2[j * 2] = ElementNodesNumbersj(j) * 2;
				IdX2[j * 2 + 1] = ElementNodesNumbersj(j) * 2 + 1;
			}
			
				
			Y = COE.col(Ej);
			
			Kej = KeArr[Ej];
			
			ElementGaussNodes = GaussNodes[ej];
			
			for(size_t qj = 0; qj < QuadUtil.NumberOfQP(); ++qj)
			{
				X = ElementGaussNodes.col(qj);
										
				Y = Narr[qj] * ElementNodesCoordj;
				
				R = X - Y;
				
				r = R.norm();
				
				Ke += Weights(qj) * Kej* phi(r,L) * Jacj[qj] * 4.0 * R(0) * R(1) /(L * L);
			}
				
			Ke *=p2;
				
			for(size_t j = 0; j < NodesPerElement * 2; ++j)		
				for(size_t k = 0; k < NodesPerElement * 2; ++k)
				{
					if(IdX2[k] <= IdX1[j])
					TripletList.push_back(Triplet(IdX1[j], IdX2[k], Ke(j ,k)));
				}
	}
			
			

	}
}
