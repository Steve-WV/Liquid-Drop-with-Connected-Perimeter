#ifndef __PDE_HPP__
#define __PDE_HPP__

#include <string>
#include "mesh.hpp"
#include <thread>
#include <mutex>



#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <Eigen/SPQRSupport>
#include <Eigen/UmfPackSupport>


class ImplicitPDESystem {
public:
	
  typedef Eigen::VectorXd Vec;
  typedef Eigen::SparseMatrix<double, Eigen::ColMajor> SpMat;
  typedef Eigen::Triplet<double> T;
  
  // Define different solvers
  typedef Eigen::SimplicialLDLT<SpMat> SpSolver; 
  //typedef Eigen::ConjugateGradient<SpMat, Eigen::Lower|Eigen::Upper> cg;
 typedef Eigen::ConjugateGradient<SpMat, Eigen::Lower|Eigen::Upper> CGSolver;
 //typedef Eigen::BiCGSTAB<SpMat, Eigen::IncompleteLUT<double> >  BCGST;
 typedef Eigen::LeastSquaresConjugateGradient<SpMat> BCGST; 
 typedef Eigen::SPQR<SpMat> SpSolver2;
 typedef Eigen::UmfPackLU<SpMat> SpSolver3;
 typedef Eigen::SimplicialLLT<SpMat> SPDsolver;
 typedef Eigen::BiCGSTAB<SpMat> Bicgstable;
 //typedef Eigen::BiCGSTAB<SpMat,_IdentityPreconditioner> Bicgstable1;
    
  inline ImplicitPDESystem( Mesh* mesh ) : m(mesh) { }
    
  void Init(std::string, bool);
  void InitU();
  void Step(int );
  void Solve();
  void Export(std::string) const;

private:
	void VLoop(double*, Mesh::TriIt, Mesh:: TriIt);
	void ELoop( double*, double*, Mesh::TriIt,Mesh::TriIt); 
	void CalcV();
	void CalcE(); 
	void PrepKM();
	void PrepNonLin(); 
	void KMLoop( SpMat*, SpMat*,SpMat*,SpMat*,SpMat*, SpMat*, SpMat*, Vec*, Mesh::TriIt, Mesh::TriIt);
	void Non_linLoop( SpMat*,Mesh::TriIt, Mesh::TriIt);
	void CalcF();
	void ForceLoop(Mesh::TriIt, Mesh::TriIt);
	
	
	void eta_vec(Vec);
	
	typedef double (* function_arg)(double arg);
	void CF_thread(double*, Vec*, function_arg, function_arg, function_arg, function_arg );
	void CalcDC();
	
	
	
  // Access to the basis functions, etc.
  Mesh* m;
    
  // basis functions (that are not constrained) and elements
  int N_dof, N_tri;
  
    double v, energy,mas,Area, energy1;
  
  Vec U, U_Db,U1,V;
  Vec F, B,F1,F2,F3,F_neumann,F_N;
  SpMat L, K, M,S1,S2,S3,S4,S_block,S_block1,S_B,S_neumann,S_N;
  
   Vec Fc; 
  double E_C; 
  
  //solver
  SpSolver solver;
  SpSolver3 solver3;
  CGSolver cg; 
  SPDsolver psolver;
  Bicgstable Bicg;
  BCGST bst;
 // par lin_solver; 
   
  // mutexes and threading
  mutable std::vector<std::mutex> idx_mutex;
  mutable std::mutex mut;
  mutable std::mutex mut1;
  std::vector<std::pair<Mesh::TriIt, Mesh::TriIt> > tri_th;
  std :: vector<std::pair<Mesh::TriIt,Mesh::TriIt>> tri_tl;
 
  // general params
   const int numThreads1 = 1; // number of threads
  const int numThreads = 1; // number of threads
 
    
};


#endif
