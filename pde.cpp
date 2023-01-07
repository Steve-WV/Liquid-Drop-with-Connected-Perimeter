#include <iostream>
#include<sstream>
#include <fstream>
#include "pde.hpp"
#include "mesh.hpp"
#include "vertex.hpp"
#include "tmv.hpp"

using namespace std; 
#include "boost/tokenizer.hpp"

#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include "pde.hpp"


double const c0 = 1.0; 
const double eps = 0.004;
const double eps0 = 0.004;
const double tau=0.5*1e-5; 

const double sigma = 0.5*(eps/eps0);
const double sigma1 = 0.5; 
const double delta = 1.0; 

const int N_conn =-4000;
const double kappa = 0.001; 


inline double W( double s ) {
	
return 1.0/4.0 * sqr(s)*sqr((1-s));

}

inline double DW( double s ) {

return 	(1/2.0) * s *(sqr(s-1)+ 2* s+(s-1));

}

inline double D2W( double s ) {

return 	 9*sqr(s)-8*s+1;
}


void ImplicitPDESystem::InitU() {
	
	auto vi = m->VertexBegin(), ve = m->VertexEnd();
	for (; vi!=ve; ++vi) {
		
	        double x = (*vi)->x();
	        double y = (*vi)->y();
		
		double r = std::sqrt(sqr(x)+sqr(y));
		double theta = atan2(x,y);

       	   	(*vi)->u() = r < 0.01+0.35*cos(2.0*theta) ? 1.0 : 0.0 ;
       	   	
	}
}


void ImplicitPDESystem::Init(std::string config_file, bool vis_only=false) {
	
	/*
	assert( (config_file!="") );
	if (config_file != "") {
		std::cerr << "Loading config file" << std::endl;
		LoadOptions(config_file);
	}*/
	
	std::cerr << "Setting boundary... ";
	auto bvi = m->BdryVertexBegin(), bve = m->BdryVertexEnd();
	for (; bvi!=bve; ++bvi) {
		(*bvi)->MarkDirBdry();
	}
	std::cerr << "done; ";
	
	// now we can index...
	std::cerr << "Indexing " << m->VertexSize() << " vertices, ";
	auto vi = m->VertexBegin(), ve = m->VertexEnd();
	N_dof = 0;
	// count basis functions
	for (; vi!=ve; ++vi) ++N_dof;  //if ( !(*vi)->DirBdry() ) We need all the points on the boundary here
	idx_mutex = std::vector<std::mutex>(N_dof);
		
	std::cerr << N_dof << " dof, " << m->VertexSize() - N_dof << " on Dirichlet boundary... ";
    
	int idx = 0, idx_bdry = 0;
	vi = m->VertexBegin();
	for (; vi!=ve; ++vi) {
	//if ( !(*vi)->DirBdry() )
		// {
			(*vi)->SetIndex(idx);
			++idx;
		//} else 	{
		//	(*vi)->SetIndex(idx_bdry+N_dof);
		//	++idx_bdry;
		//}
	}
	std::cerr << "done; ";
	Area=0.0;	
	std::cerr << "Preparing " << m->TriSize() << " triangles... ";
	N_tri = 0;
	Mesh::TriIt ti = m->TriBegin(), te = m->TriEnd();
	for (; ti!=te; ++ti) {
		(*ti)->Init_area();
		Area+= (*ti)->Init_area();
		(*ti)->Init_grad_s();
		(*ti)->dijk().length = sqrt(sqrt(3.0)*(*ti)->Area()); // i don't like this here, but we need the length of a tet
		++N_tri;
	}
		
		
	std::cerr << "Preparing " << numThreads1 << " threads... ";
	size_t elpth = ceil( (double)N_tri/numThreads);
	ti = m->TriBegin();
	for ( int i=0; i<numThreads; ++i ) {
		std::pair< Mesh::TriIt, Mesh::TriIt > bounds;
		bounds.first = ti;
		int j = 0;
		while (j<elpth && ti!=te) {++j; ++ti;}
		bounds.second = ti;
		tri_th.push_back(bounds);
	}
#ifdef EIGEN_HAS_OPENMP
	std::cerr << "Eigen is running in parallel... ";
	Eigen::setNbThreads(numThreads1);
	Eigen::initParallel();
	omp_set_num_threads(numThreads1); 
	
	
#endif 
	std::cerr << "ok; ";

	std::cerr << "zustandsvektor... ";
	U = Vec(2*N_dof);
	V = Vec(N_dof+1); 
	//U_Db = Vec(m->VertexSize()-N_dof);
	
	vi = m->VertexBegin();
	for (; vi!=ve; ++vi) {
		//if (!(*vi)->DirBdry()) 
		(*vi)->Attach( &U[(*vi)->Index()]);
		(*vi)-> Attach1( &V[(*vi) -> Index()]); 
		//else (*vi)->Attach( &U_Db[(*vi)->Index()-N_dof] );
	}
	std::cerr << "done." << std::endl;
	
	// now that we have a zustandsvector attached, we can set the values of u
	// and possible init matrices
	InitU();
	
	CalcV();
	mas=v; 
	std::cerr << "Int u of initial condition = " << v << std::endl;
	std::cerr << "Area of the domain Omega = " << Area << std::endl;
	std::cerr << "Local energy density = " << sigma << std::endl;
	
	PrepKM();
		
}


void ImplicitPDESystem::VLoop(double* v, Mesh::TriIt ti, Mesh::TriIt te ) {
	double u[NumIntPts];
	for (; ti!=te; ++ti) {
		(*ti)->Calc_u( u );
		for (int k = 0; k<NumIntPts; ++k) {
			*v += u[k] * GaussWeights[k]*(*ti)->Area();
		}
	}
}


void ImplicitPDESystem::CalcV() {
	std::vector<double> th_v;
	for ( int j = 0; j < numThreads; ++j ) {
		th_v.push_back( 0.0 );
	}
	// distribute elements on threads
	std::vector<std::thread> threads;
	for ( int j = 0; j < numThreads; ++j )
		threads.push_back( std::thread( &ImplicitPDESystem::VLoop, this, &th_v[j], tri_th[j].first, tri_th[j].second)  );

	// join threads
	for (auto &thread : threads) thread.join();
	v = 0.0;
	// add up all volume from all threads
	
	for (int j = 0; j < numThreads; ++j) v += (th_v[j]/Area);
}


// Attaching the block matrix 
void ImplicitPDESystem::KMLoop ( SpMat* th_k, SpMat* th_m, SpMat* th_k1, SpMat* th_m1,SpMat* th_mass,SpMat* th_neumann, SpMat* th_l,  Vec* th_b, Mesh::TriIt ti, Mesh::TriIt te ) {
	
	std::vector<T> vals_k, vals_m,vals_k1,vals_m1, vals_mass, vals_neumann,vals_l; 
	*th_b = Vec::Zero(N_dof);
	double u[NumIntPts];
	double sj[NumIntPts];
	double si[NumIntPts];
	
	
	Vec2 grad_sj, grad_si;
	
	for (; ti!=te; ++ti) {
	(*ti)->Calc_u( u );
	
		for (int bj=0; bj<3; ++bj) {
		(*ti)->Calc_s(bj, sj);
			
			Vertex* vert_j = (*ti)->v(static_cast<VertexName>(bj));
			int idx_j = vert_j->Index();
			
			(*ti)->Calc_grad_s(bj, grad_sj);
			
			for (int bi=0; bi<3; ++bi) {
			(*ti)->Calc_s(bi, si);
				
				Vertex* vert_i = (*ti)->v(static_cast<VertexName>(bi));
				int idx_i = vert_i->Index();
					
				(*ti)->Calc_grad_s(bi, grad_si);
				
				
				// scale or do not scale the time in the gradient flow by eps
				// Here not scaled	
				double k = (tau)*(dot(grad_sj, grad_si) * (*ti)->Area());
				double m = (1+(sigma*tau))*((bi == bj ? 1.0/6.0 : 1.0/12.0) * (*ti)->Area());
				double k1 = - (delta*sqr(eps)/c0)*((dot(grad_sj, grad_si))* (*ti) -> Area()); 
				double m1 = (bi == bj ? 1.0/6.0 : 1.0/12.0) * (*ti)->Area(); 
				double NM= (1/tau)*k; 
				 
				 
				vals_k.push_back( T(idx_i, idx_j+N_dof, k) );
				vals_m.push_back( T(idx_i, idx_j, m) );
				vals_k1.push_back(T(idx_i+N_dof,idx_j,k1)) ;  
				vals_m1.push_back(T(idx_i+N_dof,idx_j+N_dof,m1)); 
				vals_mass.push_back(T(idx_i,idx_j,m1)); 
				vals_neumann.push_back(T(idx_i,idx_j,NM)); 
					
				}
				
				double l= (1.0/3.0)* (*ti)->Area(); 
				vals_l.push_back( T(N_dof, idx_j, l) );
				vals_l.push_back( T(idx_j, N_dof, l) );			
				(*th_b)[idx_j] += (1/3.0)* (*ti)-> Area();  
			
			}
		}
	
	(*th_k).setFromTriplets( vals_k.begin(), vals_k.end() );	
	(*th_m).setFromTriplets( vals_m.begin(), vals_m.end() );
	(*th_k1).setFromTriplets( vals_k1.begin(), vals_k1.end() );
	(*th_m1).setFromTriplets( vals_m1.begin(), vals_m1.end() );
	(*th_mass).setFromTriplets( vals_mass.begin(), vals_mass.end() );
	(*th_neumann).setFromTriplets(vals_neumann.begin(), vals_neumann.end()); 
	(*th_l).setFromTriplets(vals_l.begin(), vals_l.end()); 
	
}

void ImplicitPDESystem::PrepKM() {
	std::vector<SpMat> th_k, th_m,th_k1,th_m1,th_mass,th_neumann,th_l;
	std::vector<Vec> th_b;

	for ( int j = 0; j < numThreads; ++j ) {
		th_k.push_back( SpMat(2*N_dof, 2*N_dof) );
		th_m.push_back( SpMat(2*N_dof, 2*N_dof) );
		th_k1.push_back( SpMat(2*N_dof, 2*N_dof) );
		th_m1.push_back( SpMat(2*N_dof, 2*N_dof) );
		th_mass.push_back( SpMat(N_dof,N_dof) ); 
		th_neumann.push_back ( SpMat(N_dof+1,N_dof+1)); 
		th_l.push_back(SpMat(N_dof+1,N_dof+1)); 
		th_b.push_back( Vec() );
	}
  
	std::vector<std::thread> threads;
	for ( int j = 0; j < numThreads; ++j ) {
	
		threads.push_back( std::thread( &ImplicitPDESystem::KMLoop,this, &th_k[j],&th_m[j],&th_k1[j],&th_m1[j],&th_mass[j], &th_neumann[j],&th_l[j],&th_b[j],tri_th[j].first, tri_th[j].second));
														  
	}																																			 													 
	// join threads
	for (auto &thread : threads) thread.join();
  
	//K = SpMat(N_dof,N_dof);
	M = SpMat(N_dof,N_dof);
	S_block= SpMat(2*N_dof,2*N_dof);
	S_neumann = SpMat(N_dof+1,N_dof+1); 
	L= SpMat(N_dof+1,N_dof+1); 
	B = Vec::Zero(N_dof);
	for (int j = 0; j < numThreads; ++j) {
	
	S_block+= th_k[j];
	S_block+= th_m[j];
	S_block+= th_k1[j];
	S_block+= th_m1[j];
	S_neumann += th_neumann[j]; 
	L += th_l[j]; 
	//K += th_k[j];
	M += th_mass[j];
	B += th_b[j];
	
	}
	
	//S_neumann.coeffRef(N_dof,N_dof)=1.0; 	
}


void ImplicitPDESystem::Non_linLoop ( SpMat* th_phu,Mesh::TriIt ti, Mesh::TriIt te ) {
	
	std::vector<T> vals_phu;
	double u[NumIntPts];
	double sj[NumIntPts];
	double si[NumIntPts];
	
	
	
	Vec2 grad_sj, grad_si;
	
	for (; ti!=te; ++ti) {
	(*ti)->Calc_u( u );
	
		for (int bj=0; bj<3; ++bj) {
		(*ti)->Calc_s(bj, sj);
			
			Vertex* vert_j = (*ti)->v(static_cast<VertexName>(bj));
			int idx_j = vert_j->Index();
			
			(*ti)->Calc_grad_s(bj, grad_sj);
			
			for (int bi=0; bi<3; ++bi) {
			(*ti)->Calc_s(bi, si);
				
				Vertex* vert_i = (*ti)->v(static_cast<VertexName>(bi));
				int idx_i = vert_i->Index();
					
				
				double phu= 0.0;
				
				for (int k = 0; k<NumIntPts; ++k) {

					// here we have a few ways to treat the nonlinearity W
					//double temp = D2W(u[k]) * sj[k]*si[k];
					//phu += temp * GaussWeights[k]*(*ti)->Area();
					//double temp =(1/c0) * ((sqr(u[k]))*sj[k]*si[k]);
					//double temp =(1/c0*eps) * ((sqr(u[k]))*sj[k]*si[k]);
					 double temp = delta*(0.5/(c0))*(2*sqr(u[k])-3*u[k]+1) * sj[k] *si[k];
					        phu += temp * GaussWeights[k]*(*ti)->Area();
				
						}
								
								//phu*=(1/eps);  
								vals_phu.push_back(T(idx_i+N_dof,idx_j,phu)) ;  
				
						}
					}
			}
	
				(*th_phu).setFromTriplets( vals_phu.begin(), vals_phu.end() );
}


void ImplicitPDESystem::PrepNonLin() {
	std::vector<SpMat> th_phu;
	
	for ( int j = 0; j < numThreads; ++j ) {
		th_phu.push_back( SpMat(2*N_dof, 2*N_dof) );
		
					}
  
	std::vector<std::thread> threads;
	for ( int j = 0; j < numThreads; ++j ) 
		threads.push_back( std::thread( &ImplicitPDESystem::Non_linLoop, this, &th_phu[j], tri_th[j].first, tri_th[j].second)  );
	// join threads
	for (auto &thread : threads) thread.join();
  
	S_block1= SpMat(2*N_dof,2*N_dof);
	for (int j = 0; j < numThreads; ++j) {
	
	S_block1+= th_phu[j];
	
	}
}	
	
	
void ImplicitPDESystem::ForceLoop(Mesh::TriIt ti, Mesh::TriIt te) {
	
	double sj[NumIntPts];
	double u[NumIntPts];
	Vec2 p[NumIntPts];
	
	for (; ti!=te; ++ti ) {
		(*ti)->Calc_Coord( p );
		(*ti)->Calc_u( u );
		
		for (int bj=0; bj<3; ++bj) {
			Vertex* vert_j = (*ti)->v(static_cast<VertexName>(bj));
			
			(*ti)->Calc_s(bj, sj);
				
			double f = 0.0;
			double f1 = 0.0;
			for (int k = 0; k<NumIntPts; ++k) {
						
				double temp = DW(u[k]) * sj[k];
						
				f += temp * GaussWeights[k]*(*ti)->Area();
				
				double temp1 = -D2W(u[k])*u[k]*sj[k];
				
				f1 += temp * GaussWeights[k]*(*ti)->Area();
			}
				
			//if (!vert_j->DirBdry()) {
				int idx = vert_j->Index();
				std::lock_guard<std::mutex> lock( idx_mutex[idx] );
				F[idx] += f;
				F[idx]+= f1; 
				//lock_u.unlock();
			//}

		}
	}
		
} 


void ImplicitPDESystem::CalcF() {
	
	F = Vec::Zero(N_dof);
	std::vector<std::thread> threads;
	for ( int j = 0; j < numThreads; ++j ) 
		threads.push_back( std::thread( &ImplicitPDESystem::ForceLoop, this, tri_th[j].first, tri_th[j].second)  );
 
	// join threads
	for (auto &thread : threads) thread.join();
  
}

void ImplicitPDESystem::ELoop(double* e, double* e1,  Mesh::TriIt ti, Mesh::TriIt te ) {
	double u[NumIntPts];
	Vec2 grad_u, grad_v;
	//double s=0.0;
	for (; ti!=te; ++ti) {
		(*ti)->Calc_u( u );
		(*ti)->Calc_grad_u( grad_u );
		(*ti)->Calc_grad_v( grad_v );

		for (int k = 0; k<NumIntPts; ++k) {
			*e += delta * W(u[k]) * GaussWeights[k]*(*ti)->Area();
			//*e += (delta/eps) * W(u[k]) * GaussWeights[k]*(*ti)->Area();
		}
		*e += ( (delta*sqr(eps))/(2.0) * grad_u.normsqr()  ) * (*ti)->Area();
		*e1 +=( (sigma1/(2.0))* grad_v.normsqr())* (*ti)-> Area();   
	}
}

void ImplicitPDESystem::CalcE() {
	std::vector<double> th_e, th_e1;
	for ( int j = 0; j < numThreads; ++j ) {
		th_e.push_back( 0.0 );
		th_e1.push_back( 0.0); 
	}
	// distribute elements on threads
	std::vector<std::thread> threads1;
	for ( int j = 0; j < numThreads; ++j )
		threads1.push_back( std::thread( &ImplicitPDESystem::ELoop,
	this, &th_e[j],&th_e1[j], tri_th[j].first, tri_th[j].second)  );

	// join threads
	for (auto &thread : threads1) thread.join();
	energy = 0.0;
	energy1=0.0; 
	// add up all energies from all threads
	
	for (int j = 0; j < numThreads; ++j)
	{
			energy += th_e[j];
			energy1 += th_e1[j]; 
	}
}
	

void ImplicitPDESystem::Step(int step) {


std::cerr << "Step." << std::endl;

double t1,t2; 
t1= omp_get_wtime();
Eigen::setNbThreads(numThreads);
Eigen::initParallel();
omp_set_num_threads(numThreads); 
	
PrepNonLin();

bst.setMaxIterations(0.005*N_dof); 
bst.setTolerance(1e-12); 
	
S_B=S_block-S_block1;
	
bst.compute(S_B);

U1=Vec::Zero(N_dof);
F1 =Vec::Zero(N_dof); 
U1+=U.head(N_dof);
	 
//F1+= M*U1+((sigma*tau*mas))*B; // right hand side if we scale the time by eps 
F1+= M*U1+(sigma*tau*mas)*B;

// for certain treatments of the nonlinearity we need CalcF() or not	
//CalcF();

F2=Vec::Zero(N_dof);
//F2 += (1/sqr(eps))*F;
	
// Include disconnectedness penalty or not. Remove /*, if you like do penalize disconnectedness	
	/*if (step>N_conn) {
		CalcDC();
		//Fc *= -(kappa*tau)/sqr(eps);
		//Fc *= (kappa*tau)/(eps); 
		Fc *= (kappa*tau)/(sqr(eps)); 
		//Fc= Fc.cwiseAbs();
	std::cerr<< " " << Fc.minCoeff() << " " << Fc.maxCoeff() << std::endl; 
	std::cerr<< " " << Fc(N_dof) << std::endl; 
		//F += Fc.cwiseAbs();
		F2+= Fc; 
			} */
		
Vec F3(2*N_dof) ;	
F3<< F1,F2 ; 
U= bst.solve(F3);
CalcV();

t2 = omp_get_wtime();
	
std::cout<< "time needed=" << t2-t1<< std::endl; 
	
}
	

void ImplicitPDESystem::Solve() {
	
	int k = 0;
	double t11,t22; 
	Export("out" + std::to_string(k) + ".vtk");
	while (k< 5.0*(1e6)) {
		Step(k);
		++k;
		
		// save energies and vtks
		if (k%5000 ==0)
		{
		
			CalcV();
			std::cerr << "Int u of initial condition = " << v << std::endl;
			S_N= S_neumann+L;	
			solver.compute(S_N);
			F_neumann = Vec::Zero(N_dof); 
			F_neumann += M*U1-mas*B; 
		
			Vec F_N(N_dof+1); 
			F_N << F_neumann, 0.0; 
			V= solver.solve(F_N);  
			CalcE(); 	
		 	Export("out" + std::to_string(k) + ".vtk");
		 
			std::fstream dataablage; 
			dataablage.open("daten.dat", ios::app); 
			dataablage << "Energy perimeter ="<< energy<< std::endl; 
			dataablage << "Non local energy ="<< energy1<< std::endl; 
	 		dataablage.close(); 
		}
	
	}
}
