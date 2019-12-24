#include <valarray>
#include <iostream>
#include <ostream>
#include "Program.hpp"
#include <limits>
#include <array>

#include "DimensionType.hpp"
#include "ReferenceShape.hpp"

#include "ReferenceShapeNode.hpp"
#include "ReferenceShapeEdge.hpp"
#include "ReferenceShapeFaceTriangle.hpp"
#include "ReferenceShapeFaceQuadrilateral.hpp"
#include "ReferenceShapeVolumeTetrahedron.hpp"
#include "ReferenceShapeVolumePyramid.hpp"
#include "ReferenceShapeVolumeWedge.hpp"
#include "ReferenceShapeVolumeHexahedron.hpp"

#include "Config.hpp"
#include "Point.hpp"
#include "CellToNodes.hpp"

#include "Input/Medit.hpp"

#include "Output/Medit.hpp"
#include "Output/Vtk.hpp"
#include "AHF/Mesh.hpp"

class MeshBuilder3D
{
public:
  template <VolumeType::enum_t 	_volume> static AHF::Mesh3D * Transfinite(const int_t 		nx_,
									  const int_t 		ny_,
									  const int_t 		nz_)
  {
    AHF::MeshTopology3D * topology = BuildTransfiniteTopology<_volume>	(nx_,
									 ny_,
									 nz_);
    
    AHF::MeshGeometry * geometry   = BuildTransfiniteGeometry<_volume>	(topology,
									 nx_,
									 ny_,
									 nz_);
    return new AHF::Mesh3D(topology,geometry);
  };
  
private: template <VolumeType::enum_t 	_volume> static AHF::MeshTopology3D * BuildTransfiniteTopology	(const int_t 		nx_,
													 const int_t 		ny_,
													 const int_t 		nz_);

  
private: template <VolumeType::enum_t 	_volume> static AHF::MeshGeometry * BuildTransfiniteGeometry	(AHF::MeshTopology3D * 	topology_,
													 const int_t 		nx_,
													 const int_t 		ny_,
													 const int_t 		nz_);  
};



template<>
AHF::MeshGeometry * MeshBuilder3D::BuildTransfiniteGeometry<VolumeType::Hexahedron>(AHF::MeshTopology3D * 	topology_,
										    const int_t 		nx_,
										    const int_t 		ny_,
										    const int_t 		nz_)
{  
  AHF::MeshGeometry* geometry = new AHF::MeshGeometry(topology_);
  using pts_t = typename AHF::MeshGeometry::pts_t;
  const int_t
    nxm1	= nx_-1,
    nym1	= ny_-1,
    nzm1	= nz_-1;

  static constexpr const double s_one(1.0);
  static constexpr const double s_two(2.0);
  static constexpr const double s_negOne(-1.0);
  
  const auto dx 	= s_two/double(nxm1);
  const auto dy 	= s_two/double(nym1);
  const auto dz 	= s_two/double(nzm1);
  
  int_t nodeIndex	= 0;
  
  auto z = s_negOne;
  for (int_t iz=0;iz<nzm1;++iz)
    {
      auto x = s_negOne;
      for (int_t ix=0;ix<nxm1;++ix)
	{
	  
	  { auto y = s_negOne;
	    for (int_t iy=0;iy<nym1;++iy)
	      {
		(*geometry)[nodeIndex++] = pts_t(x,y,z,1000);
		y += dy;
	      } }
	  
	  (*geometry)[nodeIndex++] = pts_t(x,s_one,z,1000);
	  x += dx;
	} 
      
      {
	auto y = s_negOne;
	for (int_t iy=0;iy<nym1;++iy)
	  {
	    (*geometry)[nodeIndex++] = pts_t(s_one,y,z,1000);
	    y += dy;
	  }	
	(*geometry)[nodeIndex++] = pts_t(s_one,s_one,z,1000);
      }
      
      z += dz;
    }
  
  {
    auto x = -1.0;
    for (int_t ix=0;ix<nxm1;++ix)
      {
	{ auto y = s_negOne;
	  for (int_t iy=0;iy<nym1;++iy)
	    {
	      (*geometry)[nodeIndex++] = pts_t(x,y,s_one,1000);
	      y += dy;
	    } }
	
	(*geometry)[nodeIndex++] = pts_t(x,s_one,s_one,1000);
	x += dx;	
      }
    
    { auto y = s_negOne;
      for (int_t iy=0;iy<nym1;++iy)
	{
	  (*geometry)[nodeIndex++] = pts_t(s_one,y,s_one,1000);
	  y += dy;
	}
      (*geometry)[nodeIndex++] = pts_t(s_one,s_one,s_one,1000); }
  }

#define access(_i,_j,_k) (_k)*nx_*ny_ + (_i)*ny_ + (_j)

  (*geometry)[access(0,0,0)].SetCod(1);
  (*geometry)[access(nxm1,0,0)].SetCod(2);
  (*geometry)[access(nxm1,nym1,0)].SetCod(3);
  (*geometry)[access(0,nym1,0)].SetCod(4);

  (*geometry)[access(0,0,nzm1)].SetCod(5);
  (*geometry)[access(nxm1,0,nzm1)].SetCod(6);
  (*geometry)[access(nxm1,nym1,nzm1)].SetCod(7);
  (*geometry)[access(0,nym1,nzm1)].SetCod(8);


  for (int_t ix = 1;ix < nxm1;++ix)
    {
      (*geometry)[access(ix,0,0)].SetCod(12);
    }

  for (int_t ix = 1;ix < nxm1;++ix)
    {
      (*geometry)[access(ix,nym1,0)].SetCod(34);
    }

  for (int_t ix = 1;ix < nxm1;++ix)
    {
      (*geometry)[access(ix,0,nzm1)].SetCod(56);
    }

  for (int_t ix = 1;ix < nxm1;++ix)
    {
      (*geometry)[access(ix,nym1,nzm1)].SetCod(78);
    }

  for (int_t iy = 1;iy < nym1;++iy)
    {
      (*geometry)[access(0,iy,0)].SetCod(14);
    }

  for (int_t iy = 1;iy < nym1;++iy)
    {
      (*geometry)[access(nxm1,iy,0)].SetCod(23);
    }

  for (int_t iy = 1;iy < nym1;++iy)
    {
      (*geometry)[access(0,iy,nzm1)].SetCod(58);
    }

  for (int_t iy = 1;iy < nym1;++iy)
    {
      (*geometry)[access(nxm1,iy,nzm1)].SetCod(67);
    }


  for (int_t iz = 1;iz < nzm1;++iz)
    {
      (*geometry)[access(0,0,iz)].SetCod(15);
    }
  for (int_t iz = 1;iz < nzm1;++iz)
    {
      (*geometry)[access(nxm1,0,iz)].SetCod(26);
    }
  for (int_t iz = 1;iz < nzm1;++iz)
    {
      (*geometry)[access(nxm1,nym1,iz)].SetCod(37);
    }
  for (int_t iz = 1;iz < nzm1;++iz)
    {
      (*geometry)[access(0,nym1,iz)].SetCod(48);
    }

    for (int_t ix = 1;ix < nxm1;++ix)
    {
      for (int_t iy = 1;iy < nym1;++iy)
	{
	  (*geometry)[access(ix,iy,0)].SetCod(101);
	}
    }  

    for (int_t ix = 1;ix < nxm1;++ix)
    {
      for (int_t iy = 1;iy < nym1;++iy)
	{
	  (*geometry)[access(ix,iy,nzm1)].SetCod(102);
	}
    }  



    for (int_t iz = 1;iz < nzm1;++iz)
    {
      for (int_t iy = 1;iy < nym1;++iy)
	{
	  (*geometry)[access(0,iy,iz)].SetCod(103);
	}
    }  
    
    for (int_t iz = 1;iz < nzm1;++iz)
      {
      for (int_t iy = 1;iy < nym1;++iy)
	{
	  (*geometry)[access(nxm1,iy,iz)].SetCod(104);
	}
    }  

    for (int_t iz = 1;iz < nzm1;++iz)
    {
      for (int_t ix = 1;ix < nxm1;++ix)
	{
	  (*geometry)[access(ix,0,iz)].SetCod(105);
	}
    }  
    
    for (int_t iz = 1;iz < nzm1;++iz)
      {
      for (int_t ix = 1;ix < nxm1;++ix)
	{
	  (*geometry)[access(ix,nym1,iz)].SetCod(106);
	}
    }  

    
    
    
#undef access
  
  return geometry;

};

template<>
AHF::MeshGeometry * MeshBuilder3D::BuildTransfiniteGeometry<VolumeType::Wedge>(AHF::MeshTopology3D * 	topology_,
									       const int_t 		nx_,
									       const int_t 		ny_,
									       const int_t 		nz_)
{
  return MeshBuilder3D::BuildTransfiniteGeometry<VolumeType::Hexahedron>(topology_,
									 nx_,
									 ny_,
									 nz_);

};

template<>
AHF::MeshGeometry * MeshBuilder3D::BuildTransfiniteGeometry<VolumeType::Tetrahedron>(AHF::MeshTopology3D * 	topology_,
										     const int_t 		nx_,
										     const int_t 		ny_,
										     const int_t 		nz_)
{
  return MeshBuilder3D::BuildTransfiniteGeometry<VolumeType::Hexahedron>(topology_,
									 nx_,
									 ny_,
									 nz_);

};


template <>
AHF::MeshTopology3D * MeshBuilder3D::BuildTransfiniteTopology<VolumeType::Hexahedron>(const int_t 		nx_,
										      const int_t 		ny_,
										      const int_t 		nz_)
{

  const int_t
    nx 		= nx_,
    ny 		= ny_,
    nxm1	= nx_-1,
    nym1	= ny_-1,
    nzm1	= nz_-1;
  
  int_t numCellsPerKind[VolumeType::NumKinds]{0};
  numCellsPerKind[VolumeType::Hexahedron] = nxm1*nym1*nzm1;
  
  AHF::MeshTopology3D * topology = new AHF::MeshTopology3D(nx_*ny_*nz_,
							   numCellsPerKind);

  int_t cncelm[8];
  int_t ith = 0;
#define _mflatIndex(_i,_j,_k) (_k)*nx*ny + (_i)*ny + (_j)
  for (int_t ielm=0;ielm<nxm1;++ielm)
    {
      for (int_t jelm=0;jelm<nym1;++jelm)
	{
	  for (int_t kelm=0;kelm<nzm1;++kelm)
	    {
	      cncelm[0] = _mflatIndex(ielm+1,jelm+1,kelm);
	      cncelm[1] = _mflatIndex(ielm,jelm+1,kelm);
	      cncelm[2] = _mflatIndex(ielm,jelm,kelm);
	      cncelm[3] = _mflatIndex(ielm+1,jelm,kelm);
	      cncelm[4] = _mflatIndex(ielm+1,jelm+1,kelm+1);
	      cncelm[5] = _mflatIndex(ielm,jelm+1,kelm+1);
	      cncelm[6] = _mflatIndex(ielm,jelm,kelm+1);
	      cncelm[7] = _mflatIndex(ielm+1,jelm,kelm+1);
	      topology->SetCellToNodes<VolumeType::Hexahedron>(ith++,
							       cncelm);
	    }
	}
    }
#undef _mflatIndex
  return topology;
};


template <>
AHF::MeshTopology3D * MeshBuilder3D::BuildTransfiniteTopology<VolumeType::Wedge>(const int_t 		nx_,
										 const int_t 		ny_,
										 const int_t 		nz_)
{
  const int_t
    nx 		= nx_,
    ny 		= ny_,
    nxm1	= nx_-1,
    nym1	= ny_-1,
    nzm1	= nz_-1;
  
  int_t numCellsPerKind[VolumeType::NumKinds]{0};
  numCellsPerKind[VolumeType::Wedge] = nxm1*nym1*nzm1 * 2;
  
  AHF::MeshTopology3D * topology = new AHF::MeshTopology3D(nx_*ny_*nz_,
							   numCellsPerKind);

  int_t cncelm[6];
  int_t ith = 0;

#define access(_i,_j,_k) (_k)*nx*ny + (_i)*ny + (_j)

  { int_t kelm;
    for (kelm=0;kelm<nzm1;++kelm)
      {
	
	{ int_t ielm;
	  for (ielm=0;ielm<nxm1;++ielm)
	    {
	      if (ielm%2==0)
		{
		  { int_t jelm;
		    for (jelm=0;jelm<nym1;++jelm)
		      {
			if (jelm%2==0)
			  {
			    cncelm[0] = access(ielm,jelm,kelm);
			    cncelm[1] = access(ielm+1,jelm,kelm);
			    cncelm[2] = access(ielm+1,jelm+1,kelm);
			    cncelm[3] = access(ielm,jelm,kelm+1);
			    cncelm[4] = access(ielm+1,jelm,kelm+1);
			    cncelm[5] = access(ielm+1,jelm+1,kelm+1);
			    
			    topology->SetCellToNodes<VolumeType::Wedge>(ith++,
									cncelm);
			    
			    cncelm[0] = access(ielm,jelm,kelm);
			    cncelm[1] = access(ielm+1,jelm+1,kelm);
			    cncelm[2] = access(ielm,jelm+1,kelm);
			    cncelm[3] = access(ielm,jelm,kelm+1);
			    cncelm[4] = access(ielm+1,jelm+1,kelm+1);
			    cncelm[5] = access(ielm,jelm+1,kelm+1);
			    topology->SetCellToNodes<VolumeType::Wedge>(ith++,
									cncelm);
			  }
			else
			  {
			    cncelm[0] = access(ielm,jelm,kelm);
			    cncelm[1] = access(ielm+1,jelm,kelm);
			    cncelm[2] = access(ielm,jelm+1,kelm);
			    cncelm[3] = access(ielm,jelm,kelm+1);
			    cncelm[4] = access(ielm+1,jelm,kelm+1);
			    cncelm[5] = access(ielm,jelm+1,kelm+1);
			    topology->SetCellToNodes<VolumeType::Wedge>(ith++,
									cncelm);

			    cncelm[0] = access(ielm+1,jelm,kelm);
			    cncelm[1] = access(ielm+1,jelm+1,kelm);
			    cncelm[2] = access(ielm,jelm+1,kelm);
			    cncelm[3] = access(ielm+1,jelm,kelm+1);
			    cncelm[4] = access(ielm+1,jelm+1,kelm+1);
			    cncelm[5] = access(ielm,jelm+1,kelm+1);
			    topology->SetCellToNodes<VolumeType::Wedge>(ith++,
									cncelm);

			  }
		      } }
		}
	      else
		{
		  { int_t jelm;
		    for (jelm=0;jelm<nym1;++jelm)
		      {
			if (jelm%2==0)
			  {
			    cncelm[0] = access(ielm,jelm,kelm);
			    cncelm[1] = access(ielm+1,jelm,kelm);
			    cncelm[2] = access(ielm,jelm+1,kelm);
			    cncelm[3] = access(ielm,jelm,kelm+1);
			    cncelm[4] = access(ielm+1,jelm,kelm+1);
			    cncelm[5] = access(ielm,jelm+1,kelm+1);
			    topology->SetCellToNodes<VolumeType::Wedge>(ith++,
									cncelm);
			    cncelm[0] = access(ielm+1,jelm,kelm);
			    cncelm[1] = access(ielm+1,jelm+1,kelm);
			    cncelm[2] = access(ielm,jelm+1,kelm);
			    cncelm[3] = access(ielm+1,jelm,kelm+1);
			    cncelm[4] = access(ielm+1,jelm+1,kelm+1);
			    cncelm[5] = access(ielm,jelm+1,kelm+1);
			    topology->SetCellToNodes<VolumeType::Wedge>(ith++,
									cncelm);
			  }
			else
			  {
			    cncelm[0] = access(ielm,jelm,kelm);
			    cncelm[1] = access(ielm+1,jelm,kelm);
			    cncelm[2] = access(ielm+1,jelm+1,kelm);
			    cncelm[3] = access(ielm,jelm,kelm+1);
			    cncelm[4] = access(ielm+1,jelm,kelm+1);
			    cncelm[5] = access(ielm+1,jelm+1,kelm+1);
			    topology->SetCellToNodes<VolumeType::Wedge>(ith++,
									cncelm);
			    cncelm[0] = access(ielm,jelm,kelm);
			    cncelm[1] = access(ielm+1,jelm+1,kelm);
			    cncelm[2] = access(ielm,jelm+1,kelm);
			    cncelm[3] = access(ielm,jelm,kelm+1);
			    cncelm[4] = access(ielm+1,jelm+1,kelm+1);
			    cncelm[5] = access(ielm,jelm+1,kelm+1);
			    topology->SetCellToNodes<VolumeType::Wedge>(ith++,
									cncelm);
				    
			  }
		      } }
		}
	    } }
      } }
	  
  if (nym1%2!=0)
    {
#if 0
      printf("coin superieur gauche\n");
#endif
      { int_t ielm = 0;
	int_t jelm = nym1-1;
	{ int_t kelm;
	  for (kelm=0;kelm<nzm1;++kelm)
	    {
			    
	      cncelm[0] = access(ielm,jelm,kelm);
	      cncelm[1] = access(ielm+1,jelm,kelm);
	      cncelm[2] = access(ielm,jelm+1,kelm);
	      cncelm[3] = access(ielm,jelm,kelm+1);
	      cncelm[4] = access(ielm+1,jelm,kelm+1);
	      cncelm[5] = access(ielm,jelm+1,kelm+1);
	      topology->SetCellToNodes<VolumeType::Wedge>(kelm * nxm1*nym1*2 + 2*(nym1-1),
							  cncelm);
	      
	      cncelm[0] = access(ielm+1,jelm,kelm);
	      cncelm[1] = access(ielm+1,jelm+1,kelm);
	      cncelm[2] = access(ielm,jelm+1,kelm);
	      cncelm[3] = access(ielm+1,jelm,kelm+1);
	      cncelm[4] = access(ielm+1,jelm+1,kelm+1);
	      cncelm[5] = access(ielm,jelm+1,kelm+1);
	      topology->SetCellToNodes<VolumeType::Wedge>(kelm * nxm1*nym1*2+2*(nym1-1)+1,
							  cncelm);
	      

	    } }
      } 
		    
      if (nxm1%2==0)
	{
#if 0
	  printf("coin superieur droit\n");
#endif

	  { int_t ielm = nxm1-1;
	    int_t jelm = nym1-1;

	    { int_t kelm;
	      for (kelm=0;kelm<nzm1;++kelm)
		{

		  cncelm[0] = access(ielm,jelm,kelm);
		  cncelm[1] = access(ielm+1,jelm,kelm);
		  cncelm[2] = access(ielm+1,jelm+1,kelm);
		  cncelm[3] = access(ielm,jelm,kelm+1);
		  cncelm[4] = access(ielm+1,jelm,kelm+1);
		  cncelm[5] = access(ielm+1,jelm+1,kelm+1);
		  topology->SetCellToNodes<VolumeType::Wedge>(kelm * nxm1*nym1*2+2*(nym1*nxm1-1),
							  cncelm);
		  cncelm[0] = access(ielm,jelm,kelm);
		  cncelm[1] = access(ielm+1,jelm+1,kelm);
		  cncelm[2] = access(ielm,jelm+1,kelm);
		  cncelm[3] = access(ielm,jelm,kelm+1);
		  cncelm[4] = access(ielm+1,jelm+1,kelm+1);
		  cncelm[5] = access(ielm,jelm+1,kelm+1);
		  topology->SetCellToNodes<VolumeType::Wedge>(kelm * nxm1*nym1*2+2*(nym1*nxm1-1)+1,
							      cncelm);
		  
		} }
	  }
	}
      else
	{
#if 0
	  printf("coin inferieur droit\n");
#endif
	  { int_t ielm = nxm1-1;
	    int_t jelm = 0;
	    { int_t kelm;
	      for (kelm=0;kelm<nzm1;++kelm)
		{
		  cncelm[0] = access(ielm,jelm,kelm);
		  cncelm[1] = access(ielm+1,jelm,kelm);
		  cncelm[2] = access(ielm,jelm+1,kelm);
		  cncelm[3] = access(ielm,jelm,kelm+1);
		  cncelm[4] = access(ielm+1,jelm,kelm+1);
		  cncelm[5] = access(ielm,jelm+1,kelm+1);
		  topology->SetCellToNodes<VolumeType::Wedge>(kelm * nxm1*nym1*2+2*(nym1*nxm1-nym1),
							      cncelm);
		  		

		  cncelm[0] = access(ielm+1,jelm,kelm);
		  cncelm[1] = access(ielm+1,jelm+1,kelm);
		  cncelm[2] = access(ielm,jelm+1,kelm);
		  cncelm[3] = access(ielm+1,jelm,kelm+1);
		  cncelm[4] = access(ielm+1,jelm+1,kelm+1);
		  cncelm[5] = access(ielm,jelm+1,kelm+1);
		topology->SetCellToNodes<VolumeType::Wedge>(kelm * nxm1*nym1*2+2*(nym1*nxm1-nym1)+1,
							      cncelm);
		  		

		} }
	  }
			
	}
    }
  else if (nxm1%2!=0)
    {
#if 0
      printf("coin superieur droit\n");
#endif
      { int_t ielm = nxm1-1;
	int_t jelm = nym1-1;
	{ int_t kelm;
	  for (kelm=0;kelm<nzm1;++kelm)
	    {

	      cncelm[0] = access(ielm,jelm,kelm);
	      cncelm[1] = access(ielm+1,jelm,kelm);
	      cncelm[2] = access(ielm+1,jelm+1,kelm);
	      cncelm[3] = access(ielm,jelm,kelm+1);
	      cncelm[4] = access(ielm+1,jelm,kelm+1);
	      cncelm[5] = access(ielm+1,jelm+1,kelm+1);
	      topology->SetCellToNodes<VolumeType::Wedge>(kelm * nxm1*nym1*2+2*(nym1*nxm1-1),
							  cncelm);
	      

	      cncelm[0] = access(ielm,jelm,kelm);
	      cncelm[1] = access(ielm+1,jelm+1,kelm);
	      cncelm[2] = access(ielm,jelm+1,kelm);
	      cncelm[3] = access(ielm,jelm,kelm+1);
	      cncelm[4] = access(ielm+1,jelm+1,kelm+1);
	      cncelm[5] = access(ielm,jelm+1,kelm+1);
	      topology->SetCellToNodes<VolumeType::Wedge>(kelm * nxm1*nym1*2+2*(nym1*nxm1-1)+1,
							  cncelm);


	    } }
      }
#if 0
      printf("coin inferieur droit\n");
#endif
      { int_t ielm = nxm1-1;
	int_t jelm = 0;

	  for (int_t kelm=0;kelm<nzm1;++kelm)
	    {
	      cncelm[0] = access(ielm,jelm,kelm);
	      cncelm[1] = access(ielm+1,jelm,kelm);
	      cncelm[2] = access(ielm,jelm+1,kelm);
	      cncelm[3] = access(ielm,jelm,kelm+1);
	      cncelm[4] = access(ielm+1,jelm,kelm+1);
	      cncelm[5] = access(ielm,jelm+1,kelm+1);
	      topology->SetCellToNodes<VolumeType::Wedge>(kelm * nxm1*nym1*2+2*(nym1*nxm1-nym1),
							  cncelm);

	      cncelm[0] = access(ielm+1,jelm,kelm);
	      cncelm[1] = access(ielm+1,jelm+1,kelm);
	      cncelm[2] = access(ielm,jelm+1,kelm);
	      cncelm[3] = access(ielm+1,jelm,kelm+1);
	      cncelm[4] = access(ielm+1,jelm+1,kelm+1);
	      cncelm[5] = access(ielm,jelm+1,kelm+1);
	      topology->SetCellToNodes<VolumeType::Wedge>(kelm * nxm1*nym1*2+2*(nym1*nxm1-nym1)+1,
							  cncelm);


	    } 
      }
    }
  
#undef access

  return topology;
};





static int_t addTetraMin(const int_t* 			cnc_,
			 const int_t 			ith_,
			 const int_t 			dec_,
			 const int 			sign_,
			 AHF::MeshTopology3D* 		topology_)
{
  int_t ith = ith_;
  int_t cncelm[4];
  int_t imin=0;
  int_t imax=0;
  int_t i=0;

  if (cnc_[1]<cnc_[imin])
    imin=1;
  if (cnc_[2]<cnc_[imin])
    imin=2;
  if (cnc_[1]>cnc_[imax])
    imax=1;
  if (cnc_[2]>cnc_[imax])
    imax=2;
  if ( (i==imin) || (i==imax) )
    ++i;
  if ( (i==imin) || (i==imax) )
    ++i;

  { int_t j;
    j=imin;
    imin=imax;
    imax = j; }

  cncelm[0] 	= cnc_[0];
  cncelm[1] 	= cnc_[1];
  cncelm[2] 	= cnc_[2];
  cncelm[3] 	= cnc_[imin]+dec_;
  topology_->SetCellToNodes<VolumeType::Tetrahedron>(ith++,
							  cncelm);

  cncelm[0] 	= cnc_[0]+dec_;
  cncelm[2] 	= cnc_[1]+dec_;
  cncelm[1] 	= cnc_[2]+dec_;
  cncelm[3] 	= cnc_[imax];
  topology_->SetCellToNodes<VolumeType::Tetrahedron>(ith++,
							  cncelm);

  if (sign_>0)
    {
      cncelm[0] 	= cnc_[i];
      cncelm[1] 	= cnc_[imax];
      cncelm[2] 	= cnc_[i]+dec_;
      cncelm[3] 	= cnc_[imin]+dec_;
    }
  else
    {
      cncelm[0] 	= cnc_[i];
      cncelm[2] 	= cnc_[imax];
      cncelm[1] 	= cnc_[i]+dec_;
      cncelm[3] 	= cnc_[imin]+dec_;
    }
  topology_->SetCellToNodes<VolumeType::Tetrahedron>(ith++,
							  cncelm);
  return ith;
}



static int_t addTetraNoMin(const int_t* 			cnc_,
			   const int_t 			ith_,
			   const int_t 			dec_,
			   const int		sign_,
			   AHF::MeshTopology3D* 		topology_)
{
  int_t ith = ith_;
  int_t cncelm[4];
  int_t imin=0;
  int_t imax=0;
  int_t i=0;

  if (cnc_[1]<cnc_[imin])
    imin=1;
  if (cnc_[2]<cnc_[imin])
    imin=2;
  if (cnc_[1]>cnc_[imax])
    imax=1;
  if (cnc_[2]>cnc_[imax])
    imax=2;
  if ( (i==imin) || (i==imax) )
    ++i;
  if ( (i==imin) || (i==imax) )
    ++i;

  { int_t j;
    j=imin;
    imin=imax;
    imax = j; }

  cncelm[0] 	= cnc_[0];
  cncelm[1] 	= cnc_[1];
  cncelm[2] 	= cnc_[2];
  cncelm[3] 	= cnc_[imin]+dec_;
  topology_->SetCellToNodes<VolumeType::Tetrahedron>(ith++,
							  cncelm);

  cncelm[0] 	= cnc_[0]+dec_;
  cncelm[2] 	= cnc_[1]+dec_;
  cncelm[1] 	= cnc_[2]+dec_;
  cncelm[3] 	= cnc_[i];
  topology_->SetCellToNodes<VolumeType::Tetrahedron>(ith++,
							  cncelm);

  if (sign_>0)
    {
      cncelm[0] 	= cnc_[imax];
      cncelm[1] 	= cnc_[i];
      cncelm[2] 	= cnc_[imin]+dec_;
      cncelm[3] 	= cnc_[imax]+dec_;
    }
  else
    {
      cncelm[0] 	= cnc_[imax];
      cncelm[2] 	= cnc_[i];
      cncelm[1] 	= cnc_[imin]+dec_;
      cncelm[3] 	= cnc_[imax]+dec_;
    }
  topology_->SetCellToNodes<VolumeType::Tetrahedron>(ith++,
							  cncelm);
  return ith;
}


static int_t addTetraNoMax(const int_t* 		cnc_,
			   const int_t 			ith_,
			   const int_t 			dec_,
			   const int			sign_,
			   AHF::MeshTopology3D* 	topology_)
{
  int_t ith = ith_;
  int_t cncelm[4];
  int_t imin=0;
  int_t imax=0;
  int_t i=0;

  if (cnc_[1]<cnc_[imin])
    imin=1;
  if (cnc_[2]<cnc_[imin])
    imin=2;
  if (cnc_[1]>cnc_[imax])
    imax=1;
  if (cnc_[2]>cnc_[imax])
    imax=2;
  if ( (i==imin) || (i==imax) )
    ++i;
  if ( (i==imin) || (i==imax) )
    ++i;

  cncelm[0] 	= cnc_[0];
  cncelm[1] 	= cnc_[1];
  cncelm[2] 	= cnc_[2];
  cncelm[3] 	= cnc_[i]+dec_;
  topology_->SetCellToNodes<VolumeType::Tetrahedron>(ith++,
							  cncelm);

  cncelm[0] 	= cnc_[0]+dec_;
  cncelm[2] 	= cnc_[1]+dec_;
  cncelm[1] 	= cnc_[2]+dec_;
  cncelm[3] 	= cnc_[imax];
  topology_->SetCellToNodes<VolumeType::Tetrahedron>(ith++,
							  cncelm);

  if (sign_>0)
    {
      cncelm[0] 	= cnc_[imin];
      cncelm[1] 	= cnc_[imax];
      cncelm[2] 	= cnc_[imin]+dec_;
      cncelm[3] 	= cnc_[i]+dec_;
    }
  else
    {
      cncelm[0] 	= cnc_[imin];
      cncelm[2] 	= cnc_[imax];
      cncelm[1] 	= cnc_[imin]+dec_;
      cncelm[3] 	= cnc_[i]+dec_;
    }
  topology_->SetCellToNodes<VolumeType::Tetrahedron>(ith++,
							  cncelm);
  return ith;
}




static int_t addTetraMax(const int_t* 			cnc_,
		     const int_t 			ith_,
		     const int_t 			dec_,
		     const int 			sign_,
		AHF::MeshTopology3D* 		topology_)
{
  int_t ith = ith_;
  int_t cncelm[4];
  int_t imin=0;
  int_t imax=0;
  int_t i=0;
  if (cnc_[1]>cnc_[imax])
    imax=1;
  if (cnc_[2]>cnc_[imax])
    imax=2;
  if (cnc_[1]<cnc_[imin])
    imin=1;
  if (cnc_[2]<cnc_[imin])
    imin=2;
  if ( (i==imin) || (i==imax) )
    ++i;
  if ( (i==imin) || (i==imax) )
    ++i;
  cncelm[0] 	= cnc_[0];
  cncelm[1] 	= cnc_[1];
  cncelm[2] 	= cnc_[2];
  cncelm[3] 	= cnc_[imin]+dec_;
  topology_->SetCellToNodes<VolumeType::Tetrahedron>(ith++,
							  cncelm);

  cncelm[0] 	= cnc_[0]+dec_;
  cncelm[2] 	= cnc_[1]+dec_;
  cncelm[1] 	= cnc_[2]+dec_;
  cncelm[3] 	= cnc_[imax];
  topology_->SetCellToNodes<VolumeType::Tetrahedron>(ith++,
							  cncelm);
  if (sign_>0)
    {
      cncelm[0] 	= cnc_[i];
      cncelm[1] 	= cnc_[imax];
      cncelm[2] 	= cnc_[i]+dec_;
      cncelm[3] 	= cnc_[imin]+dec_;
    }
  else
    {
      cncelm[0] 	= cnc_[i];
      cncelm[2] 	= cnc_[imax];
      cncelm[1] 	= cnc_[i]+dec_;
      cncelm[3] 	= cnc_[imin]+dec_;
    }
  topology_->SetCellToNodes<VolumeType::Tetrahedron>(ith++,
							  cncelm);
  return ith;
}

template <>
AHF::MeshTopology3D * MeshBuilder3D::BuildTransfiniteTopology<VolumeType::Tetrahedron>(const int_t 		nx_,
										       const int_t 		ny_,
										       const int_t 		nz_)
{
  const int_t
    nx 		= nx_,
    ny 		= ny_,
    nxm1	= nx_-1,
    nym1	= ny_-1,
    nzm1	= nz_-1;
  
  int_t numCellsPerKind[VolumeType::NumKinds]{0};
  numCellsPerKind[VolumeType::Tetrahedron] = nxm1*nym1*nzm1 * 6;  
  AHF::MeshTopology3D * topology = new AHF::MeshTopology3D(nx_*ny_*nz_,
							   numCellsPerKind);  
  int_t ith = 0;  
  int_t cnctria[3];
#define access(_i,_j,_k) (_k)*nx*ny + (_i)*ny + (_j)
  for (int_t kelm=0;kelm<nzm1;++kelm)
    {
      if (kelm%2==0)
	{		    
	  for (int_t ielm=0;ielm<nxm1;++ielm)
	    {
	      if (ielm%2==0)
		{
		  for (int_t jelm=0;jelm<nym1;++jelm)
		    {
		      if (jelm%2==0)
			{
			  cnctria[0] = access(ielm,jelm,kelm);
			  cnctria[1] = access(ielm+1,jelm,kelm);
			  cnctria[2] = access(ielm+1,jelm+1,kelm);					
			  ith =  addTetraMin(cnctria,ith,nx*ny,1,topology);
			  cnctria[0] = access(ielm,jelm,kelm);
			  cnctria[1] = access(ielm+1,jelm+1,kelm);
			  cnctria[2] = access(ielm,jelm+1,kelm);
			  ith =  addTetraMin(cnctria,ith,nx*ny,-1,topology);
			}
		      else
			{
			  cnctria[0] = access(ielm,jelm,kelm);
			  cnctria[1] = access(ielm+1,jelm,kelm);
			  cnctria[2] = access(ielm,jelm+1,kelm);
			  ith 	=  addTetraNoMin(cnctria,ith,nx*ny,-1,topology);
			  cnctria[0] = access(ielm+1,jelm,kelm);
			  cnctria[1] = access(ielm+1,jelm+1,kelm);
			  cnctria[2] = access(ielm,jelm+1,kelm);
			  ith =  addTetraMin(cnctria,ith,nx*ny,1,topology);					

			}
		    } 
		}
	      else
		{

		  for (int_t jelm=0;jelm<nym1;++jelm)
		    {
		      if (jelm%2!=0)
			{
			  cnctria[0] = access(ielm,jelm,kelm);
			  cnctria[1] = access(ielm+1,jelm,kelm);
			  cnctria[2] = access(ielm+1,jelm+1,kelm);
					
			  ith =  addTetraNoMax(cnctria,ith,nx*ny,1,topology);
					
			  cnctria[0] = access(ielm,jelm,kelm);
			  cnctria[1] = access(ielm+1,jelm+1,kelm);
			  cnctria[2] = access(ielm,jelm+1,kelm);
			  ith =  addTetraNoMax(cnctria,ith,nx*ny,-1,topology);

			}
		      else
			{
					
			  cnctria[0] = access(ielm,jelm,kelm);
			  cnctria[1] = access(ielm+1,jelm,kelm);
			  cnctria[2] = access(ielm,jelm+1,kelm);
			  ith =  addTetraNoMax(cnctria,ith,nx*ny,-1,topology);

			  cnctria[0] = access(ielm+1,jelm,kelm);
			  cnctria[1] = access(ielm+1,jelm+1,kelm);
			  cnctria[2] = access(ielm,jelm+1,kelm);
			  ith =  addTetraNoMin(cnctria,ith,nx*ny,1,topology);					
			}
		    } 
		}
	    } 
	}
      else
	{


	  for (int_t ielm=0;ielm<nxm1;++ielm)
	    {
	      if (ielm%2==0)
		{

		  for (int_t jelm=0;jelm<nym1;++jelm)
		    {
		      if (jelm%2==0)
			{
			  cnctria[0] = access(ielm,jelm,kelm);
			  cnctria[1] = access(ielm+1,jelm,kelm);
			  cnctria[2] = access(ielm+1,jelm+1,kelm);					
			  ith =  addTetraMax(cnctria,ith,nx*ny,-1,topology);
			  cnctria[0] = access(ielm,jelm,kelm);
			  cnctria[1] = access(ielm+1,jelm+1,kelm);
			  cnctria[2] = access(ielm,jelm+1,kelm);
			  ith =  addTetraMax(cnctria,ith,nx*ny,1,topology);
			}
		      else
			{
			  cnctria[0] = access(ielm,jelm,kelm);
			  cnctria[1] = access(ielm+1,jelm,kelm);
			  cnctria[2] = access(ielm,jelm+1,kelm);
			  ith 	=  addTetraNoMax(cnctria,ith,nx*ny,-1,topology);
			  cnctria[0] = access(ielm+1,jelm,kelm);
			  cnctria[1] = access(ielm+1,jelm+1,kelm);
			  cnctria[2] = access(ielm,jelm+1,kelm);
			  ith =  addTetraMax(cnctria,ith,nx*ny,-1,topology);					

			}
		    } 
		}
	      else
		{
		      
		  for (int_t jelm=0;jelm<nym1;++jelm)
		    {
		      if (jelm%2!=0)
			{
			  cnctria[0] = access(ielm,jelm,kelm);
			  cnctria[1] = access(ielm+1,jelm,kelm);
			  cnctria[2] = access(ielm+1,jelm+1,kelm);
					
			  ith =  addTetraNoMin(cnctria,ith,nx*ny,1,topology);
					
			  cnctria[0] = access(ielm,jelm,kelm);
			  cnctria[1] = access(ielm+1,jelm+1,kelm);
			  cnctria[2] = access(ielm,jelm+1,kelm);
			  ith =  addTetraNoMin(cnctria,ith,nx*ny,-1,topology);

			}
		      else
			{
					
			  cnctria[0] = access(ielm,jelm,kelm);
			  cnctria[1] = access(ielm+1,jelm,kelm);
			  cnctria[2] = access(ielm,jelm+1,kelm);
			  ith =  addTetraNoMin(cnctria,ith,nx*ny,-1,topology);

			  cnctria[0] = access(ielm+1,jelm,kelm);
			  cnctria[1] = access(ielm+1,jelm+1,kelm);
			  cnctria[2] = access(ielm,jelm+1,kelm);
			  ith =  addTetraNoMax(cnctria,ith,nx*ny,1,topology);					
			}
		    } 
		}
	    }

	}
    } 
#undef  access

  return topology;
};







class MnsMakeMesh3D : public Program
{
  int_t m_nxyz[3];
  std::string m_nameVolume;
public:
  
  MnsMakeMesh3D(int 		argc,
		char * 	argv[]) : Program(argc,argv,true)
  {
    this->AddOption(new Option<int_t>(&this->m_nxyz[0],
				      "-x",
				      7,
				      true,
				      "The number of points in X-direction"));
    
    this->AddOption(new Option<int_t>(&this->m_nxyz[1],
				      "-y",
				      7,
				      true,
				      "The number of points in Y-direction"));
    this->AddOption(new Option<int_t>(&this->m_nxyz[2],
				      "-z",
				      7,
				      true,
				      "The number of points in Z-direction"));
    
    this->AddOption(new Option<std::string>(&m_nameVolume,
					    "-e",
					    std::string("Hexahedron"),
					    true,
					    "The volume kind."));    
  };
  
  virtual void Main()
  {
    if (this->HasVerbose())
      {
	std::cout  << "element "  << m_nameVolume << std::endl;
	std::cout  << "nx " 	 << m_nxyz[0] 	<< std::endl;
	std::cout  << "ny " 	 << m_nxyz[1] 	<< std::endl;
	std::cout  << "nz " 	 << m_nxyz[2] 	<< std::endl;
	std::cout  << "ofilename " << this->GetOfilename() << std::endl;
      }

    AHF::Mesh3D * mesh;
    if ("Hexahedron" == m_nameVolume)
      {
	mesh = MeshBuilder3D::Transfinite<VolumeType::Hexahedron>(m_nxyz[0],
								  m_nxyz[1],
								  m_nxyz[2]);
      }
    else if ("Tetrahedron" == m_nameVolume)
      {
	mesh = MeshBuilder3D::Transfinite<VolumeType::Tetrahedron>(m_nxyz[0],
								   m_nxyz[1],
								   m_nxyz[2]);
      }
    else
      {
	mesh = MeshBuilder3D::Transfinite<VolumeType::Wedge>(m_nxyz[0],
							     m_nxyz[1],
							     m_nxyz[2]);
      }

    
    {
      Output::Medit outputMedit(this->GetOfilename().c_str());
      outputMedit << *mesh;
    }
    
    {
      Output::Vtk::Writer outputVtkWriter(this->GetOfilename().c_str());
      outputVtkWriter << *mesh;
    }
    
    delete [] mesh;    
  };
};


int main(int 	argc_, 
	 char* argv_[])
{
  MnsMakeMesh3D application(argc_,argv_);
  application.Run();
  return 0;
}
