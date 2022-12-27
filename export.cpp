#include <iomanip>
#include <fstream>
#include <unordered_map>
#include "pde.hpp"
#include "mesh.hpp"
#include "vertex.hpp"
#include "tmv.hpp"
/*
#ifdef HAVE_VTK
#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkTetra.h>
#include <vtkCellArray.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPointData.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkFloatArray.h>

#include <chrono>

void ImplicitPDESystem::Export(std::string fname) const 
{ 
	std::lock_guard<std::mutex> lock(mut);
	
	auto start = std::chrono::steady_clock::now();
	
	auto points = vtkSmartPointer<vtkPoints>::New();
	points->SetNumberOfPoints(m->VertexSize());
	Mesh::const_VertexIt vi = m->VertexBegin(), ve = m->VertexEnd(); 
	float p[3]; 
	for (int j=0; vi!=ve; ++vi) {
		p[0] = (*vi)->x();
		p[1] = (*vi)->y();
		p[2] = (*vi)->z();
		points->SetPoint(j++, p);
	}
	
	auto stop = std::chrono::steady_clock::now();
	auto diff =  stop - start;
	start = stop;
	std::cerr << std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl;
  
	// we need some indexes
	std::unordered_map<Vertex*, int> Vt_list_idx;
	vi = m->VertexBegin(); 
	for (int j = 0; vi!=ve; ++vi) Vt_list_idx[*vi] = j++;
	
	stop = std::chrono::steady_clock::now();
	diff =  stop - start;
	start = stop;
	std::cerr << std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl;
 	
	auto Conn = vtkSmartPointer<vtkIdTypeArray>::New();
	Conn->SetNumberOfValues(5*m->TetSize());
	
	auto ti = m->TetBegin(), te = m->TetEnd(); 
	for (int j = 0; ti!=te; ++ti) {
		Conn->SetValue( j++, 4 );
		Conn->SetValue( j++, Vt_list_idx[(*ti)->a()] );
		Conn->SetValue( j++, Vt_list_idx[(*ti)->b()] );
		Conn->SetValue( j++, Vt_list_idx[(*ti)->c()] );
		Conn->SetValue( j++, Vt_list_idx[(*ti)->d()] );
	}

	auto Tets = vtkSmartPointer<vtkCellArray>::New();
	Tets->SetCells(m->TetSize(),Conn);
	
	stop = std::chrono::steady_clock::now();
	diff =  stop - start;
	start = stop;
	std::cerr << std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl;
 
	auto unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	unstructuredGrid->SetPoints(points);
	unstructuredGrid->SetCells(VTK_TETRA, Tets);
	
	stop = std::chrono::steady_clock::now();
	diff =  stop - start;
	start = stop;
	std::cerr << std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl;
  
	auto u = vtkSmartPointer<vtkFloatArray>::New();
	u->SetNumberOfComponents(1);
	u->SetName("u");
	u->SetNumberOfValues(m->VertexSize());
	vi = m->VertexBegin(); 
	for (int j = 0; vi!=ve; ++vi) u->SetTuple1(j++, (*vi)->u());
	unstructuredGrid->GetPointData()->AddArray(u);
	
	auto v = vtkSmartPointer<vtkFloatArray>::New();
	v->SetNumberOfComponents(1);
	v->SetName("v");
	v->SetNumberOfValues(m->VertexSize());
	vi = m->VertexBegin(); 
	for (int j = 0; vi!=ve; ++vi) v->SetTuple1(j++, (*vi)->v());
	unstructuredGrid->GetPointData()->AddArray(v);
	
	auto H = vtkSmartPointer<vtkFloatArray>::New();
	H->SetNumberOfComponents(1);
	H->SetName("H");
	H->SetNumberOfValues(m->VertexSize());
	vi = m->VertexBegin(); 
	for (int j = 0; vi!=ve; ++vi) {
		if (W((*vi)->u())>0.1) H->SetTuple1(j++, -1.0/2.0*1.0/sqrt(2.0*W((*vi)->u()))*(*vi)->v() + H0((*vi)->Coord()) );
		else H->SetTuple1(j++, vtkMath::Nan());
	}
	unstructuredGrid->GetPointData()->AddArray(H);
	
	
	
  	stop = std::chrono::steady_clock::now();
  	diff =  stop - start;
  	start = stop;
  	std::cerr << std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl;
  

	// Write file
	auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
	writer->SetDataModeToAppended();
	writer->EncodeAppendedDataOff();
	std::string fn = fname+".vtu";
	writer->SetFileName(fn.c_str());
#if VTK_MAJOR_VERSION <= 5
	writer->SetInput(unstructuredGrid);
#else
	writer->SetInputData(unstructuredGrid);
#endif
	writer->Write();
	
	stop = std::chrono::steady_clock::now();
	diff =  stop - start;
	start = stop;
	std::cerr << std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl;
  
}
#endif
#ifndef HAVE_VTK */
void ImplicitPDESystem::Export(std::string fname) const {
	
	std::lock_guard<std::mutex> lock(mut);
	std::ofstream os;
	os.open(fname);
	os << std::setprecision(16);
	os << "# vtk DataFile Version 2.0" << std::endl;
	os << "s2d output data" << std::endl;
	os << "ASCII" << std::endl;
	os << "DATASET POLYDATA" << std::endl;
	
	os << "POINTS " << m->VertexSize() << " float" << std::endl;
	Mesh::const_VertexIt vi = m->VertexBegin(), ve = m->VertexEnd();
	for (; vi!=ve; ++vi)
		os << (*vi)->x() << " " << (*vi)->y() << " " << "0.0" << std::endl;

	// we need some indexes
	std::unordered_map<Vertex*, int> Vt_list_idx;
	vi = m->VertexBegin(); int idx = 0;
	for (; vi!=ve; ++vi) {
		Vt_list_idx[*vi] = idx;
		++idx;
	}

	os << "POLYGONS " << m->TriSize() << " " << m->TriSize()*4 << std::endl;
	Mesh::const_TriIt ti = m->TriBegin(), te = m->TriEnd();
	for (; ti!=te; ++ti)
		os << "3 " << Vt_list_idx[(*ti)->a()] << " " << Vt_list_idx[(*ti)->b()] << " "
			<< Vt_list_idx[(*ti)->c()] << " " << std::endl;
    
	os << "POINT_DATA " << m->VertexSize() << std::endl;
	os << "SCALARS u float" << std::endl;
	os << "LOOKUP_TABLE default" << std::endl;
	vi = m->VertexBegin();
	for (; vi!=ve; ++vi)
		os << (*vi)->u() << std::endl;
    
	os << "CELL_DATA " << m->TriSize() << std::endl;
	os << "SCALARS Tri_area float" << std::endl;
	os << "LOOKUP_TABLE default" << std::endl;
	ti = m->TriBegin();
	for (; ti!=te; ++ti)
		os << (*ti)->Area() << std::endl;
	os.close();
    
}
/*
#endif
*/
