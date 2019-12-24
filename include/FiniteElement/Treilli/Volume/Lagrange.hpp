#pragma once

#include "CRTP.hpp"
#include "VolumeType.hpp"

namespace FiniteElement
{
  namespace Treilli
  {
    //
    // Declaration of the class Lagrange
    //
    namespace Volume
    {
      template <VolumeType::enum_t _volumeShape,unsigned int _degree> class Lagrange;
    };

    //
    // Specialization of the traits CRTP
    //
    template <VolumeType::enum_t _volumeShape,unsigned int _degree> struct Traits_CRTP< Volume::Lagrange<_volumeShape,_degree> >
    {
      typedef VolumeType::enum_t eTypeShape;
    };

    
    namespace Volume
    {
      
      template <VolumeType::enum_t _volumeShape,unsigned int _degree> struct Traits_Lagrange
      {
      };
      
      template <unsigned int _degree> struct Traits_Lagrange<VolumeType::Tetrahedron,_degree>
      {
	static constexpr const unsigned int NumNodes 		= ( (_degree+1)*(_degree+2)*(_degree+3))/6;
	static constexpr const unsigned int NumSubElements 	= _degree*_degree*_degree;
	static constexpr const unsigned int Dimension 		= 3;    
	static constexpr const unsigned int Degree 		= _degree;    
	static constexpr const unsigned int NumNodesInVolume	= 4;
      };
      
      template <unsigned int _degree> struct Traits_Lagrange<Mesh::Topology::HEXAHEDRON,_degree>
      {
	static constexpr const unsigned int NumNodes 		= (_degree+1)*(_degree+1)*(_degree+1);
	static constexpr const unsigned int NumSubElements 	= _degree*_degree*_degree;
	static constexpr const unsigned int Dimension 		= 3;    
	static constexpr const unsigned int Degree 		= _degree;    
	static constexpr const unsigned int NumNodesInVolume	= 8;
      };
      
      template <VolumeType::enum_t _volumeShape,unsigned int _degree> struct LagrangeUtils
      {	
      };
      
      template <unsigned int _degree> struct LagrangeUtils<Mesh::Topology::HEXAHEDRON,_degree>
      {
	using Traits = Traits_Lagrange<Mesh::Topology::HEXAHEDRON,_degree>;		
	using tSubVolumesToNodes = std::array< std::array<unsigned int, Traits::NumNodesInVolume >, Traits::NumSubElements >;
	using tIntegerNodesCoordinates = std::array< std::array<unsigned int, Traits::Dimension >, Traits::NumNodes >;
	using tLobattoNodesCoordinates = std::array< std::array<double, Traits::Dimension >, Traits::NumNodes >;

	static inline void ComputeSubcnc(tSubVolumesToNodes& subCellsToNodes_) noexcept
	{
#define _dec(_i,_j,_k) (_degree+1)*(_degree+1)*(_k) + (_degree+1) * (_i) + (_j)
	  unsigned int subCellIndex = 0;
	  for (unsigned int k=0;k<_degree;++k)
	    {
	      for (unsigned int i=0;i<_degree;i++)
		{
		  for (unsigned int j=0;j<_degree;j++)
		    {		      
		      subCellsToNodes_[subCellIndex][0] = _dec(i+1,j+1,k);
		      subCellsToNodes_[subCellIndex][1] = _dec(i,j+1,k);
		      subCellsToNodes_[subCellIndex][2] = _dec(i,j,k);
		      subCellsToNodes_[subCellIndex][3] = _dec(i+1,j,k);
		      subCellsToNodes_[subCellIndex][4] = _dec(i+1,j+1,k+1);
		      subCellsToNodes_[subCellIndex][5] = _dec(i,j+1,k+1);
		      subCellsToNodes_[subCellIndex][6] = _dec(i,j,k+1);
		      subCellsToNodes_[subCellIndex][7] = _dec(i+1,j,k+1);
		      ++subCellIndex;
		    } 
		}
	    }
#undef _dec
	};
	
	static inline void ComputeCoordinates2(tIntegerNodesCoordinates& icoo_) noexcept
	{
	  unsigned int nodeIndex = 0;
	  for (unsigned int k=0;k<_degree+1;++k)
	    {
	      for (unsigned int i=0;i<_degree+1;++i)
		{		
		  for (unsigned int j=0;j<_degree+1;j++)
		    {
		      const unsigned int nodeIndex = (k*(_degree+1)*(_degree+1)+i*(_degree+1)+j);
		      icoo_[nodeIndex][0] = i;
		      icoo_[nodeIndex][1] = j;
		      icoo_[nodeIndex][2] = k;
		    }
		} 
	    }
	};
	
	
	static inline void ComputeCoordinates(tIntegerNodesCoordinates& icoo_) noexcept
	{
	  unsigned int nodeIndex = 0;
	  icoo_[nodeIndex][0] = _degree;	  
	  icoo_[nodeIndex][1] = _degree;	
	  icoo_[nodeIndex][2] = 0;	
	  ++nodeIndex;
	  
	  icoo_[nodeIndex][0] = 0;	  
	  icoo_[nodeIndex][1] = _degree;	
	  icoo_[nodeIndex][2] = 0;	
	  ++nodeIndex;
	  
	  icoo_[nodeIndex][0] = 0;	  
	  icoo_[nodeIndex][1] = 0;	
	  icoo_[nodeIndex][2] = 0;	
	  ++nodeIndex;
	  
	  icoo_[nodeIndex][0] = _degree;	  
	  icoo_[nodeIndex][1] = 0;	
	  icoo_[nodeIndex][2] = 0;	
	  ++nodeIndex;
	  

	  icoo_[nodeIndex][0] = _degree;	  
	  icoo_[nodeIndex][1] = _degree;	
	  icoo_[nodeIndex][2] = _degree;	
	  ++nodeIndex;

	  icoo_[nodeIndex][0] = 0;	  
	  icoo_[nodeIndex][1] = _degree;	
	  icoo_[nodeIndex][2] = _degree;	
	  ++nodeIndex;

	  icoo_[nodeIndex][0] = 0;	  
	  icoo_[nodeIndex][1] = 0;	
	  icoo_[nodeIndex][2] = _degree;	
	  ++nodeIndex;

	  icoo_[nodeIndex][0] = _degree;	  
	  icoo_[nodeIndex][1] = 0;	
	  icoo_[nodeIndex][2] = _degree;	
	  ++nodeIndex;

	  static const unsigned int hexaedge_cnc[]  = { 0,1,
							1,2,
							2,3,
							3,0,
							4,5,
							5,6,
							6,7,
							7,4,
							0,4,
							1,5,
							2,6,
							3,7 }; 
	  /* deja reference dans eVolume */
	  static const unsigned int hexaface_cnc[]  = {0,3,2,1,
						       4,5,6,7,				 
						       0,1,5,4,
						       1,2,6,5,				 
						       2,3,7,6,
						       3,0,4,7 };
	
	  static const I 	refhexa_icoo[] =  { 1,1,0,
						    0,1,0,
						    0,0,0,
						    1,0,0,
						    1,1,1,
						    0,1,1,
						    0,0,1,
						    1,0,1};

	  { unsigned int reedge[8];
	    for (unsigned int iedge=0;iedge<12;++iedge)
	      {
		reedge[0] = refhexa_icoo[3*hexaedge_cnc[iedge * 2 + 0]+0];
		reedge[1] = refhexa_icoo[3*hexaedge_cnc[iedge * 2 + 0]+1];
		reedge[2] = refhexa_icoo[3*hexaedge_cnc[iedge * 2 + 0]+2];

		reedge[3] = refhexa_icoo[3*hexaedge_cnc[iedge * 2 + 1]+0];
		reedge[4] = refhexa_icoo[3*hexaedge_cnc[iedge * 2 + 1]+1];
		reedge[5] = refhexa_icoo[3*hexaedge_cnc[iedge * 2 + 1]+2];
		
		for (unsigned int i=0;i<_degree-1;++i)
		  {
		    const unsigned int l1 		= (i+1);
		    const unsigned int l0 		= _degree-(i+1);
		    icoo_[nodeIndex][0] 		= reedge[0] * l0 + reedge[3] * l1;
		    icoo_[nodeIndex][1] 		= reedge[1] * l0 + reedge[4] * l1;
		    icoo_[nodeIndex][2] 		= reedge[2] * l0 + reedge[5] * l1;
		    ++nodeIndex;
		  } 
	      }
	  }


	  { unsigned int reface[32];
	    for (unsigned int iface=0;iface<6;++iface)
	      {
		reface[0] = refhexa_icoo[3*hexaface_cnc[iface * 4 + 0]+0];
		reface[1] = refhexa_icoo[3*hexaface_cnc[iface * 4 + 0]+1];
		reface[2] = refhexa_icoo[3*hexaface_cnc[iface * 4 + 0]+2];

		reface[3] = refhexa_icoo[3*hexaface_cnc[iface * 4 + 1]+0];
		reface[4] = refhexa_icoo[3*hexaface_cnc[iface * 4 + 1]+1];
		reface[5] = refhexa_icoo[3*hexaface_cnc[iface * 4 + 1]+2];

		reface[6] = refhexa_icoo[3*hexaface_cnc[iface * 4 + 2]+0];
		reface[7] = refhexa_icoo[3*hexaface_cnc[iface * 4 + 2]+1];
		reface[8] = refhexa_icoo[3*hexaface_cnc[iface * 4 + 2]+2];

		reface[9]  = refhexa_icoo[3*hexaface_cnc[iface * 4 + 3]+0];
		reface[10] = refhexa_icoo[3*hexaface_cnc[iface * 4 + 3]+1];
		reface[11] = refhexa_icoo[3*hexaface_cnc[iface * 4 + 3]+2];

		for (unsigned int i=0;i<_degree-1;++i)
		  {
		    for (unsigned int j=0;j<_degree-1;++j)
		      {
			
			const unsigned int l0 			= (i+1)*(j+1);
			const unsigned int l1 			= (_degree-(i+1))*(j+1);
			const unsigned int l2 			= (_degree-(i+1))*(_degree-(j+1));
			const unsigned int l3 			= (_degree-(j+1))*(i+1);
			
			const unsigned int xi 			= reface[0] * l0 + reface[3] * l1 + reface[6] * l2 + reface[9] *l3;
			const unsigned int xj 			= reface[1] * l0 + reface[4] * l1 + reface[7] * l2 + reface[10]*l3;
			const unsigned int xk 			= reface[2] * l0 + reface[5] * l1 + reface[8] * l2 + reface[11]*l3;
			
			icoo_[nodeIndex][0] 	= xi/_degree;	  
			icoo_[nodeIndex][1] 	= xj/_degree;	
			icoo_[nodeIndex][2] 	= xk/_degree;
			++nodeIndex;
		      }
		  } 
	      }
	  }
	

	  /* unsigned intNSunsigned intDE */
	  for (unsigned int k=0;k<_degree-1;++k)
	    {	    
	      for (unsigned int i=0;i<_degree-1;++i)
		{
		  for (unsigned int j=0;j<_degree-1;j++)
		    {			  
		      icoo_[nodeIndex][0] = i+1;
		      icoo_[nodeIndex][1] = j+1;	
		      icoo_[nodeIndex][2] = k+1;
		      ++nodeIndex;
		    }
		}
	    }
	};
	  
      };





      template <unsigned int _degree> struct LagrangeUtils<VolumeType::Tetrahedron,_degree>
      {
	using Traits = Traits_Lagrange<VolumeType::Tetrahedron,_degree>;		
	using tSubVolumesToNodes = std::array< std::array<unsigned int, Traits::NumNodesInVolume >, Traits::NumSubElements >;
	using tIntegerNodesCoordinates = std::array< std::array<unsigned int, Traits::Dimension >, Traits::NumNodes >;
	using tLobattoNodesCoordinates = std::array< std::array<double, Traits::Dimension >, Traits::NumNodes >;
	using tTrianglesToNodes = std::array< std::array<unsigned int, 3 >, (_degree+1)*(_degree+1) >;

	
	static void TreilliVolume_compute_tetra(const unsigned int		nframe_,
						const tTrianglesToNodes& 	ecnci_,						
						tSubVolumesToNodes& 		cnci_)
	{
	  static constexpr const unsigned int s_numTriangles = (_degree+1)*(_degree+1);
	  unsigned int count	= 0;
	  count 		= 0;
	  unsigned int nb 			= nframe_-1;
	  unsigned int dec0 		= 0;
	  unsigned int dec1 		= ((nb+1)*(nb+2))/2;
	  unsigned int N = s_numTriangles;
	  for (unsigned int irot=0;irot<nb;++irot)
	    {
	      /*	fprintf(stderr,"%d %d\n",dec0,dec1);*/
	      N=(nb-irot-1)*(nb-irot-1);

	      for (unsigned int i=0;i<N;++i)
		{
		    
		  unsigned int mx = ecnci_[i][0];
		    
		  for (unsigned int j=1;j<3;++j)
		    {
		      if (ecnci_[i][j]<mx)		    
			{
			  mx = ecnci_[i][j];
			}
		    } 
		    
		  cnci_[count][0]=dec0+ecnci_[i][0];
		  cnci_[count][1]=dec0+ecnci_[i][1];
		  cnci_[count][2]=dec0+ecnci_[i][2];
		  cnci_[count][3]=dec1+mx;
		  ++count;
		    
		  unsigned int mx2 = ecnci_[i][0];
		    
		  for (unsigned int j=1;j<3;++j)
		    {
		      if (ecnci_[i][j]>mx2)		    
			{
			  mx2 = ecnci_[i][j];
			}
		    } 
		    
		  cnci_[count][0]=dec1+ecnci_[i][0];
		  cnci_[count][2]=dec1+ecnci_[i][1];
		  cnci_[count][1]=dec1+ecnci_[i][2];
		  cnci_[count][3]=dec0+mx2;
		  ++count;
		} 
	


	      for (unsigned int i=0;i<N;++i)
		{
		  unsigned int j1,j2,j3,mx1 = ecnci_[i][0];
		  j1=0;
		  j2=0;
		  j3=0;
		  unsigned int j;
		  for (j=1;j<3;++j)
		    {
		      if (mx1>ecnci_[i][j])
			{
			  mx1 = ecnci_[i][j];
			  j1=j;
			}
		    }
		  I mx2 = ecnci_[i][0];
		  for (j=1;j<3;++j)
		    {
		      if (mx2<ecnci_[i][j])
			{
			  mx2 = ecnci_[i][j];
			  j2=j;
			}
		    }
		  if ((j1==0)&&(j2==1))
		    {
		      j3 = 2;
		    }
		  else if ((j1==1)&&(j2==2))
		    {
		      j3 = 0;
		    }
		  else if ((j1==0)&&(j2==2))
		    {
		      j3 = 1;
		    }
		  else if ((j1==1)&&(j2==0))
		    {
		      j3 = 2;
		    }
		  else if ((j1==2)&&(j2==1))
		    {
		      j3 = 0;
		    }
		  else if ((j1==2)&&(j2==0))
		    {
		      j3 = 1;
		    }

		  cnci_[count][0]=dec0+ecnci_[i][j3];
		  cnci_[count][1]=dec0+mx2;
		  cnci_[count][2]=dec1+ecnci_[i][j3];
		  cnci_[count][3]=dec1+mx1;
		  if (cnci_[count][0]!=cnci_[count][1]-1)
		    {
		      cnci_[count][0]=dec0+ecnci_[i][j3];
		      cnci_[count][2]=dec0+mx2;
		      cnci_[count][1]=dec1+ecnci_[i][j3];
		      cnci_[count][3]=dec1+mx1;		  
		    }
#if 0
		  printf("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA "ifmt" "ifmt" "ifmt" "ifmt"\n",cnci_[count][0],cnci_[count][1],cnci_[count][2],cnci_[count][3]);
#endif
		    
		  ++count;
		} 


	      unsigned int n = N;
	      N=(nb-irot)*(nb-irot);
		
	      for (unsigned int i=n;i<N;++i)
		{
		  unsigned int mx = ecnci_[i][0];

		  for (unsigned int j=1;j<3;++j)
		    {
		      if (ecnci_[i][j]<mx)		    
			{
			  mx = ecnci_[i][j];
			}
		    } 
		  cnci_[count][0]=dec0+ecnci_[i][0];
		  cnci_[count][1]=dec0+ecnci_[i][1];
		  cnci_[count][2]=dec0+ecnci_[i][2];
		  cnci_[count][3]=dec1+mx;
		  ++count;
		} 

	      for (unsigned int i=n+1;i<N;i+=2)
		{
		  unsigned int j1,j2,j3,mx1 = ecnci_[i][0];
		  j1=0;
		  j2=0;
		  j3=0;
		  unsigned int j;
		  for (j=1;j<3;++j)
		    {
		      if (mx1>ecnci_[i][j])
			{
			  mx1 = ecnci_[i][j];
			  j1=j;
			}
		    }
		  unsigned int mx2 = ecnci_[i][0];
		  for (j=1;j<3;++j)
		    {
		      if (mx2<ecnci_[i][j])
			{
			  mx2 = ecnci_[i][j];
			  j2=j;
			}
		    }
		  if ((j1==0)&&(j2==1))
		    j3 = 2;
		  else if ((j1==1)&&(j2==2))
		    j3 = 0;
		  else if ((j1==0)&&(j2==2))
		    j3 = 1;
		  else if ((j1==1)&&(j2==0))
		    j3 = 2;
		  else if ((j1==2)&&(j2==1))
		    j3 = 0;
		  else if ((j1==2)&&(j2==0))
		    j3 = 1;
		  cnci_[count][0]=dec0+ecnci_[i][j3];
		  cnci_[count][1]=dec0+mx2;
		  cnci_[count][2]=dec1+ecnci_[i][j3];
		  cnci_[count][3]=dec1+mx1;
		  if (cnci_[count][0]!=cnci_[count][1]-1)
		    {
		
		      cnci_[count][0]=dec0+ecnci_[i][j3];
		      cnci_[count][2]=dec0+mx2;
		      cnci_[count][1]=dec1+ecnci_[i][j3];
		      cnci_[count][3]=dec1+mx1;
		
		    }
		  ++count;
		} 
	    
	    
	      dec0 = dec1;
	      dec1 = dec1+((nb-irot)*(nb-irot+1))/2;
	    }


	};
	
	
	static inline void ComputeSubcnc(tSubVolumesToNodes& subCellsToNodes_) noexcept
	{

	  tTrianglesToNodes icnctria;
#define _dec(_i,_j) ((  (_i)+(_j) + 1 )*( (_i)+(_j) ))/2+(_i)
	  unsigned int triangleIndex = 0;
	  for (unsigned int j=0;j<(_degree+1);j++)
	    {
	      
	      for (unsigned int i=0;i<j;i++)
		{
		  icnctria[triangleIndex][0] 	= _dec(i,j-i);
		  icnctria[triangleIndex][1] 	= _dec(i+1,j-i);
		  icnctria[triangleIndex][2] 	= _dec(i,j+1-i); 
		  ++triangleIndex;
		  
		  icnctria[triangleIndex][0] 	= _dec(i,j-i);
		  icnctria[triangleIndex][1] 	= _dec(i+1,j-1-i);
		  icnctria[triangleIndex][2] 	= _dec(i+1,j-i);
		  ++triangleIndex;		  
		}
	      
	      icnctria[triangleIndex][0]	=  _dec(j,0);
	      icnctria[triangleIndex][1]	=  _dec(j+1,0);
	      icnctria[triangleIndex][2]	=  _dec(j,1);	      
	      ++triangleIndex;
	    } 
#undef _dec	
	
	  TreilliVolume_compute_tetra(_degree+1,				      
				      icnctria,
				      subCellsToNodes_);
	};
	
	static inline void ComputeCoordinates2(tIntegerNodesCoordinates& icoo_) noexcept
	{
	  unsigned int nodeIndex = 0;
	  for (unsigned int i=0;i<_degree+1;i++)
	    {

	      for (unsigned int j=0;j<=i;j++)
		{
		  icoo_[nodeIndex][0] = j;
		  icoo_[nodeIndex][1] = i-j;
		  icoo_[nodeIndex][2] = 0;
		  ++nodeIndex;
		} 
	    }
	  
	  { I N     = ((_degree+1)*(_degree+1+1))/2-(_degree+1);
	    for (unsigned int irot=1;irot<_degree+1;++irot)
	      {
		  for (unsigned int i=0;i<N;++i)
		    {
		      icoo_[nodeIndex][0] = icoo_[i][0];
		      icoo_[nodeIndex][1] = icoo_[i][1];
		      icoo_[nodeIndex][2] = irot;
		      ++nodeIndex;
		    } 
		N = N-(_degree+1-irot);
	      } }
	};
	
	
	static inline void ComputeCoordinates(tIntegerNodesCoordinates& icoo_) noexcept
	{
	  unsigned int next = 0;
	  static constexpr const I reftetra_icoo[] = {0,0,0,
						      1,0,0,
						      0,1,0,
						      0,0,1};

	  const I 	n 		= (_degree>0)? _degree-1 : 0;
	  icoo_[next][0] = 0;	  
	  icoo_[next][1] = 0;	
	  icoo_[next][2] = 0;	
	  ++next;
	  icoo_[next][0] = _degree;	  
	  icoo_[next][1] = 0;	
	  icoo_[next][2] = 0;	
	  ++next;
	  icoo_[next][0] = 0;	  
	  icoo_[next][1] = _degree;	
	  icoo_[next][2] = 0;	
	  ++next;
	  icoo_[next][0] = 0;	  
	  icoo_[next][1] = 0;	
	  icoo_[next][2] = _degree;	
	  ++next;
	
	  static const I tetraedge_cnc[] = { 1,2,
					     2,0,
					     0,1,
					     2,3,
					     0,3,
					     1,3};			
	
	  static const I tetraface_cnc[] = {1,2,3,
					    2,0,3,
					    0,1,3,
					    0,2,1};
	  
	  { I reedge[8];
	    for (unsigned int iedge=0;iedge<6;++iedge)
	      {
		reedge[0] = reftetra_icoo[3*tetraedge_cnc[iedge * 2 + 0]+0];
		reedge[1] = reftetra_icoo[3*tetraedge_cnc[iedge * 2 + 0]+1];
		reedge[2] = reftetra_icoo[3*tetraedge_cnc[iedge * 2 + 0]+2];
		reedge[3] = reftetra_icoo[3*tetraedge_cnc[iedge * 2 + 1]+0];
		reedge[4] = reftetra_icoo[3*tetraedge_cnc[iedge * 2 + 1]+1];
		reedge[5] = reftetra_icoo[3*tetraedge_cnc[iedge * 2 + 1]+2];
		for (unsigned int i=0;i<n;++i)
		  {
		    const unsigned int l1 = (i+1);
		    const unsigned int l0 = _degree-(i+1);
		    std::cout << "ddddddd " << l0 << " " << l1 << std::endl;
		    icoo_[next][0] = reedge[0] * l0 + reedge[3] * l1;
		    icoo_[next][1] = reedge[1] * l0 + reedge[4] * l1;
		    icoo_[next][2] = reedge[2] * l0 + reedge[5] * l1;			    
		    ++next;
		  } 
	      } } 
	  

	  if (_degree>2)
	    {
	  { I reface[32];
	    for (unsigned int iface=0;iface<4;++iface)
	      {
		  reface[0] = reftetra_icoo[3*tetraface_cnc[iface * 3 + 0]+0];
		  reface[1] = reftetra_icoo[3*tetraface_cnc[iface * 3 + 0]+1];
		  reface[2] = reftetra_icoo[3*tetraface_cnc[iface * 3 + 0]+2];

		  reface[3] = reftetra_icoo[3*tetraface_cnc[iface * 3 + 1]+0];
		  reface[4] = reftetra_icoo[3*tetraface_cnc[iface * 3 + 1]+1];
		  reface[5] = reftetra_icoo[3*tetraface_cnc[iface * 3 + 1]+2];

		  reface[6] = reftetra_icoo[3*tetraface_cnc[iface * 3 + 2]+0];
		  reface[7] = reftetra_icoo[3*tetraface_cnc[iface * 3 + 2]+1];
		  reface[8] = reftetra_icoo[3*tetraface_cnc[iface * 3 + 2]+2];
		  
		  for (unsigned int i=0;i<n;++i)
		    {
		      for (unsigned int j=0;j<n-(i+1);++j)
			{
			  const unsigned int l0 		= _degree-( (i+1) + (j+1) );
			  const unsigned int l1 		= (i+1);
			  const unsigned int l2 		= (j+1);
			  //			  std::cout << l0 << " " << l1 <<  " " << l2 <<  " " << std::endl;
			  
			  
			  icoo_[next][0]	= reface[0] * l0 + reface[3] * l1 + reface[6] * l2;
			  icoo_[next][1] = reface[1] * l0 + reface[4] * l1 + reface[7] * l2;
			  icoo_[next][2] = reface[2] * l0 + reface[5] * l1 + reface[8] * l2;
			  //std::cout << icoo_[next][0] << " " << icoo_[next][1] <<  " " << icoo_[next][2] <<  " " << std::endl;
			  ++next;
			} 
		    }
	      		  
		} } 
	    }
	  /* INSIDE TETRA */

	  for (unsigned int  k=0;k<n-1;++k)
	    {
	      for (unsigned int  i=0;i<n-(k+1);++i)
		{		    
		  for (unsigned int j=0;j<n-(i+1)-(k+1);j++)
		    {			  
		      icoo_[next][0] = i+1;
		      icoo_[next][1] = j+1;	
		      icoo_[next][2] = k+1;	
		      ++next;
		    } 
		} 
	    }
	  
	};
	  
      };

    template <VolumeType::enum_t _volumeShape,unsigned int _degree> class Lagrange : public CRTP<Lagrange<_volumeShape,_degree> > 
    {
    protected:
      typedef CRTP<Lagrange<_volumeShape,_degree> > 	BaseClass;
      typedef typename BaseClass::eTypeShape eTypeShape;
      
      using Traits 			= Traits_Lagrange<_volumeShape,_degree>;		
      using tSubVolumesToNodes 		= std::array< std::array<unsigned int, Traits::NumNodesInVolume >, Traits::NumSubElements >;
      using tIntegerNodesCoordinates 	= std::array< std::array<unsigned int, Traits::Dimension >, Traits::NumNodes >;
      using tLobattoNodesCoordinates 	= std::array< std::array<double, Traits::Dimension >, Traits::NumNodes >;

      
      static constexpr const unsigned int s_numNodes 		= Traits::NumNodes;
      static constexpr const unsigned int s_numElements 	= Traits::NumSubElements;
      static constexpr const unsigned int s_dimension 		= Traits::Dimension;
      static constexpr const unsigned int s_degree 		= Traits::Degree;
      static constexpr const unsigned int s_numNodesInVolume	= Traits::NumNodesInVolume;

      tSubVolumesToNodes m_subcnc;
      tIntegerNodesCoordinates m_icoo;

    public:            

      inline unsigned int 	GetDegree() 		const noexcept { return s_degree; };
      inline unsigned int 	GetDimension() 		const noexcept { return s_dimension; };
      inline unsigned int 	GetNumNodes() 		const noexcept { return s_numNodes; };    
      inline unsigned int 	GetNumSubElements() 	const noexcept { return s_numElements; };
      inline eTypeShape 	GetShape()		const noexcept { return _volumeShape; };
      inline unsigned int 	GetNumNodesInCell() 	const noexcept { return s_numNodesInVolume; };    
      
      template <typename _float_type> inline	_float_type GetCoordinate(const unsigned int&nodeIndex_,
									  const unsigned int&dimensionIndex_) const
      {
	static const _float_type idegree = _float_type(1.0) / _float_type(_degree);
#ifndef NDEBUG
	Debug::IsInRange(__TRACE__,nodeIndex_,(unsigned int)0,this->s_numNodes-1);
	Debug::IsInRange(__TRACE__,dimensionIndex_,(unsigned int)0,this->s_dimension-1);
#endif
	return _float_type(this->m_icoo[nodeIndex_][dimensionIndex_]) * idegree;
      };
      
      inline unsigned int GetNodeIndex(const unsigned int&subElementIndex_,
				       const unsigned int&localNodeIndex_) const
      {
#ifndef NDEBUG
	Debug::IsInRange(__TRACE__,subElementIndex_,(unsigned int)0,this->s_numElements-1);
	Debug::IsInRange(__TRACE__,localNodeIndex_,(unsigned int)0,this->s_numNodesInVolume-1);
#endif
	return m_subcnc[subElementIndex_][localNodeIndex_];
      };
      
      inline Lagrange()
      {		

	// COMPUTE COORDINATES
	LagrangeUtils<_volumeShape,_degree>::ComputeCoordinates(this->m_icoo);
	{
	  std::array<unsigned int,(_degree+1)*(_degree+1)*(_degree+1)> perm;
	  tIntegerNodesCoordinates ilagr;
	  
	  LagrangeUtils<_volumeShape,_degree>::ComputeCoordinates2(ilagr);
	  
	  for (unsigned int i=0;i<s_numNodes;++i)
	    {	
	      perm[(_degree+1)*(_degree+1)*this->m_icoo[i][2] + this->m_icoo[i][0]*(_degree+1) + this->m_icoo[i][1]] = 1+i;
	    } 
	  
	  LagrangeUtils<_volumeShape,_degree>::ComputeSubcnc(this->m_subcnc); 
	  for (unsigned int subElementIndex=0;subElementIndex<s_numElements;++subElementIndex)
	    {
	      for (unsigned int iv=0;iv<s_numNodesInVolume;++iv)
		{	  
		  const unsigned int l 				= m_subcnc[subElementIndex][iv];
		  const unsigned int i 				= ilagr[l][0];
		  const unsigned int j 				= ilagr[l][1];
		  const unsigned int k 				= ilagr[l][2];		  
		  m_subcnc[subElementIndex][iv] = perm[(_degree+1)*(_degree+1)*k+(_degree+1)*i+j]-1;
		}
	    } 	
	  

	}

      };
      
      inline ~Lagrange()
      {
      };


    protected:

      
    };

  };
  };
};
