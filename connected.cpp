#include <iostream>
#include <unordered_set>
#include <queue>
#include <algorithm>
#include "boost/heap/fibonacci_heap.hpp"
#include "boost/timer/timer.hpp"
#include "utils.hpp"
#include "pde.hpp"

typedef double (* function_arg)(double arg);

const double distance = 0.068; // how far away from pure phase do we start looking (should be ~sqrt(eps))
//const double distance = 0.35; // how far away from pure phase do we start looking (should be ~sqrt(eps))
//const double distance = 0.298; // how far away from pure phase do we start looking (should be ~sqrt(eps))
//const double distance = 0.35; // how far away from pure phase do we start looking (should be ~sqrt(eps))
const double ns_con = 1.0-distance; //lower boundary for phi's support
const double m_norm = 2.0*3.0/cube(1.0-ns_con);  //ensures that integral phi from ns_con to one is one

double Wt_upper(double s) {
	if (s<ns_con)
		return 0.0;
	else
		return 0.5*sqr(s-ns_con)*m_norm;
}
double DWt_upper(double s) {
	if (s < ns_con)
		return 0.0;
	else
		return (s-ns_con)*m_norm;
}

inline double phi_upper(double s){
	if (s>ns_con)
		return 0.0;
	else{
		return 0.5*sqr(s-ns_con);
	}
}
double Dphi_upper(double s) {
	if (s>ns_con)
		return 0.0;
	else{
		return s-ns_con;
	}
}

double Wt_lower(double s) {	return Wt_upper(1-s); }
double DWt_lower(double s) { return -DWt_upper(1-s); }
double phi_lower(double s) { return phi_upper(1-s); }
double Dphi_lower(double s) { return -Dphi_upper(1-s); }

struct dist_info {
	int i,j; // component index
	double d; // distance between component i and j
	std::vector<Tri*> arc; // the geodesic curve from i to j
};


void ImplicitPDESystem::CalcDC() {
	Vec fc_upper = Vec::Zero(N_dof); // remove init in production!
	Vec fc_lower = Vec::Zero(N_dof);
	double ec_upper, ec_lower;
  
  	std::cerr << "Upper... ";
	std::thread th_upper( &ImplicitPDESystem::CF_thread, this, &ec_upper, &fc_upper, Wt_upper, DWt_upper, phi_upper, Dphi_upper);
	th_upper.join();
	std::cerr << "Lower... ";
 	std::thread th_lower( &ImplicitPDESystem::CF_thread, this, &ec_lower, &fc_lower, Wt_lower, DWt_lower, phi_lower, Dphi_lower);
	th_lower.join();
    	std::cerr << std::endl;
	//Fc =  fc_lower;
	//E_C = ec_lower;
	
	//Fc =  fc_upper;
	//E_C = ec_upper;
	
	Fc = fc_upper+fc_lower; 
	E_C = ec_upper+ec_lower; 			
}

std::vector<Tri*> FindComponent( Tri* start, std::unordered_set<Tri*> &itris, function_arg phi ) {
	assert( phi((*start).dijk().u_av) < m_deps );
	
	std::vector<Tri*> component;
	std::queue<Tri*> q;
	q.push(start);
	(*start).dijk().closed = true;
	component.push_back(start);
	itris.erase(start);
	while (!q.empty()) {
		Tri* t = q.front(); q.pop();
		for (int j = 0; j<3; ++j) {
			if (!(*t).Across(static_cast<VertexName>(j))) continue;
			Tri* nb = (*t).Across(static_cast<VertexName>(j));
			if ( !(*nb).dijk().closed && phi((*nb).dijk().u_av) < m_deps ) {
				q.push(nb);
				(*nb).dijk().closed = true;
				if ( (*nb).dijk().in_supp ) {
					component.push_back(nb);
					itris.erase(nb);
				}
			}
		}
		
	}
	return component;
}

void Dijkstra(Tri* start, Mesh* m, int sz_supp, function_arg phi ) {
	// first reset some parameters
	auto ti = m->TriBegin(), te = m->TriEnd();
	for (; ti!=te; ++ti) {
		(*ti)->dijk().dist = m_dinf;
		(*ti)->dijk().pre = NULL;
	}
	// our heap data structure
	boost::heap::fibonacci_heap <Tri*, boost::heap::compare<Dijk::tri_dist_larger> > q;
	
	//no_of_ielements=N_e;
	start->dijk().dist = 0.0;
	start->dijk().handle = q.push(start);
    
	while ( sz_supp > 0 ) {
		Tri* active = q.top(); q.pop();
    
		if ( active->dijk().in_supp ) --sz_supp; //if this was in the support, we're done with it now.
  
		for( int i=0; i<3; ++i ) {
			Tri* nb = active->Across( static_cast<VertexName>(i) );
			if ( !nb ) continue;
 
			double ew =  0.5*( phi(active->dijk().u_av)*active->dijk().length + phi(nb->dijk().u_av)*nb->dijk().length );
			
			if (  nb->dijk().dist == m_dinf ) {
				nb->dijk().dist = active->dijk().dist + ew;
				nb->dijk().pre = active;
				nb->dijk().handle = q.push(nb);
			} else if (nb->dijk().dist > active->dijk().dist + ew ) {
				nb->dijk().dist = active->dijk().dist + ew;
				nb->dijk().pre = active;
				q.increase(nb->dijk().handle);
			}
		}
	}
}
  
std::vector<Tri*> MakeGeodesic(Tri* start) {
	std::vector<Tri*> arc;
	Tri* t = start;
	while (t) {
		if ( !(t->dijk().in_supp) ) arc.push_back(t);
		t = t->dijk().pre;
	}
	return arc;
}

void ImplicitPDESystem::CF_thread(double* e_c, Vec* f, function_arg Wt, function_arg DWt, function_arg phi, function_arg Dphi) {
	std::unordered_set<Tri*> supp_Wt;
	auto ti = m->TriBegin(), te = m->TriEnd();
	for (; ti!=te; ++ti) {
		double avg = 0.0;
		for (int j=0; j<3; ++j) avg += (*ti)->v(static_cast<VertexName>(j))->u(); // for piecewise linear this is correct
		avg /= 3.0; // we just average the nodal values.
		(*ti)->dijk().u_av = avg; // save the average for this element somewhere
		(*ti)->dijk().closed = false; // tri not seen for component search
		//(*ti)->which_c = -1;
		if ( Wt(avg)>m_deps ) {
			supp_Wt.insert(*ti); // add to the set of interesting elements
			(*ti)->dijk().in_supp = true;
		} else (*ti)->dijk().in_supp = false;
	}
	int sz_supp = supp_Wt.size();
	
	// now separate into connected components
	int num_cc = 0;
	std::vector< std::vector<Tri*> > components;
	while ( !supp_Wt.empty() ) {
		std::vector<Tri*> new_component = FindComponent( *(supp_Wt.begin()), supp_Wt, phi);
		++num_cc; 
		components.push_back(new_component);
	}
	std::cerr << num_cc << ( (num_cc>1) ? " components... " : " component... " );
	
	// if there is only one connected component, we're done here.
	if (num_cc < 2) {
		*f = Vec::Zero(N_dof);
		*e_c = 0.0;
		return; 
	} 
	
	// once the components are done, we can integrate Wt.
	std::vector<double> Wt_sum;
	for (auto c : components ) {
		double s = 0.0;
		for (auto t : c ) {
			s += Wt( (*t).dijk().u_av ) * (*t).Area();
		}
		Wt_sum.push_back(s);
	}

	std::vector<dist_info> geodesics;
	for (int i=0; i<num_cc-1; ++i) {
		Dijkstra( components[i][0], m, sz_supp, phi );
		for (int j=i+1; j<num_cc; ++j) {
		    dist_info di;
			di.i = i; di.j = j;
			di.d = (components[j][0])->dijk().dist;
			di.arc = MakeGeodesic(components[j][0]);
			geodesics.push_back(di);
		}
	}

	*e_c = 0.0;
	*f = Vec::Zero(N_dof);
	for (auto g : geodesics) {
		int i = g.i, j = g.j;
		*e_c += 2.0*Wt_sum[i]*Wt_sum[j]*g.d;
		
		for (auto t : components[j]) {
			for (int k=0; k<3; ++k) {
				Vertex* v = t->v(static_cast<VertexName>(k));
				if (v->DirBdry()) continue;
				int idx = v->Index();
				
				(*f)[idx] += 2.0*Wt_sum[i]*g.d * DWt(t->dijk().u_av) * 1.0/3.0 * t->Area(); // integral over a basis function is 1/3 * Vol
			}
		}
		for (auto t : components[i]) {
			for (int k=0; k<3; ++k) {
				Vertex* v = t->v(static_cast<VertexName>(k));
				if (v->DirBdry()) continue;
				int idx = v->Index();
			
				(*f)[idx] += 2.0*Wt_sum[j]*g.d * DWt(t->dijk().u_av) * 1.0/3.0 * t->Area(); // integral over a basis function is 1/3
			}
		}
		
		for (auto t : g.arc) {
			for (int k=0; k<3; ++k) {
				Vertex* v = t->v(static_cast<VertexName>(k));
				if (v->DirBdry()) continue;
				int idx = v->Index();
			
				(*f)[idx] += 2.0*Wt_sum[i]*Wt_sum[j] * Dphi(t->dijk().u_av) * 1.0/3.0 * t->dijk().length; // integral over a basis function is 1/3
			}
		}
		
	}
	
	
	/*
	// check component assignment.
	int j = 0;
	for (auto c : components ) { 
		for (auto t : c ) (*t).which_c = j;
		++j;
	}
	for (auto g : geodesics) {
		for (auto t : g.arc) t->which_c = -10.0;	
	}*/
	
}

