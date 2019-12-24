// #include "Quadrature/Edge/Legendre.hpp"

#include "CRTP.hpp"
#include "FaceType.hpp"
#include "VolumeType.hpp"


int tetlobatto(int n,double* w,double *rst)
{
  int m = n-1;
  int i,j,k,l;
  int next = 0;
  for (i=0;i<m+1;++i)
    {
      for (j=0;j<m+1-i;++j)
	{
	  l = m-i-j;
	  rst[3*next+0] = (1.0+2.0*w[i]-w[j]-w[l]) / 3.0;
	  rst[3*next+1] = (1.0+2.0*w[j]-w[i]-w[l]) / 3.0;
	  rst[3*next+2] = 0.0;
	  ++next;
	}
    }

  for (j=1;j<=m;++j)
    {
      for (k=2;k<=m+2-j;++k)
	{
	  l = m+3-j-k;
	  rst[3*next+0] = 0.0;
	  rst[3*next+1] = (1.0+2.0*w[j-1]-w[k-1]-w[l-1]) / 3.0;
	  rst[3*next+2] = (1.0+2.0*w[k-1]-w[j-1]-w[l-1]) / 3.0;
	  ++next;

	}
    }

  for (i=2;i<=m;++i)
    {
      for (k=2;k<=m+2-i;++k)
	{
	  l = m+3-i-k;

	  rst[3*next+0] = (1.0+2.0*w[i-1]-w[k-1]-w[l-1]) / 3.0;
	  rst[3*next+1] = 0.0;
	  rst[3*next+2] = (1.0+2.0*w[k-1]-w[i-1]-w[l-1]) / 3.0;
	  
	  ++next;
	}
    }

  for (i=2;i<=m;++i)
    {
      for (j=2;j<=m+1-i;++j)
	{
	  l = m+3-i-j;
	  
	  double xi = (1.0+2.0*w[i-1]-w[j-1]-w[l-1]) / 3.0;
	  double eta = (1.0+2.0*w[j-1]-w[i-1]-w[l-1]) / 3.0;
	  rst[3*next+0] = xi;
	  rst[3*next+1] = eta;
	  rst[3*next+2] = 1.0-xi-eta;
	  ++next;
	  
	}
    }

  for (i=2;i<=m;++i)
    {
      for (j=2;j<=m+1-i;++j)
	{
	  for (k=2;k<=m+2-i-j;++k)
	    {
	      l = m+4-i-j-k;
	      rst[3*next+0] = (1.0+3.0*w[i-1]-w[j-1]-w[k-1]-w[l-1]) / 4.0;
	      rst[3*next+1] = (1.0+3.0*w[j-1]-w[i-1]-w[k-1]-w[l-1]) / 4.0;
	      rst[3*next+2] = (1.0+3.0*w[k-1]-w[i-1]-w[j-1]-w[l-1]) / 4.0;
	      ++next;
	    }
	}
    }
  return next;  
}




template <typename impl_t,
	  VolumeType::enum_t 	_volumeShape,
	  unsigned int 	_degree>
class treilli3d_t;

template <VolumeType::enum_t _volumeShape,unsigned int _degree> struct treilli3d_traits_t
{
};

template <unsigned int _degree> struct treilli3d_traits_t<VolumeType::Tetrahedron,_degree>
{
  static constexpr const unsigned int nnodes 		= ( (_degree+1)*(_degree+2)*(_degree+3))/6;
  static constexpr const unsigned int NumSubElements 	= _degree*_degree*_degree;
  static constexpr const unsigned int dim 		= 3;    
  static constexpr const unsigned int Degree 		= _degree;    
  static constexpr const unsigned int nnodesInVolume	= 4;
};

template <unsigned int _degree> struct treilli3d_traits_t<VolumeType::Hexahedron,_degree>
{
  static constexpr const unsigned int nnodes 		= (_degree+1)*(_degree+1)*(_degree+1);
  static constexpr const unsigned int NumSubElements 	= _degree*_degree*_degree;
  static constexpr const unsigned int dim 		= 3;    
  static constexpr const unsigned int Degree 		= _degree;    
  static constexpr const unsigned int nnodesInVolume	= 8;
};

template <unsigned int _degree> struct treilli3d_traits_t<VolumeType::Wedge,_degree>
{
  ;
  static constexpr const unsigned int nnodes 		= ( ( (_degree+1) *(_degree+2) ) /2 ) * (_degree+1);
  static constexpr const unsigned int NumSubElements 	= _degree*_degree*_degree;
  static constexpr const unsigned int dim 		= 3;    
  static constexpr const unsigned int Degree 		= _degree;    
  static constexpr const unsigned int nnodesInVolume	= 6;
};


template <VolumeType::enum_t _volumeShape,unsigned int _degree> struct treilli3d_utils
{	
};

template <unsigned int _degree> struct treilli3d_utils<VolumeType::Hexahedron,_degree>
{
  using traits_t 		= treilli3d_traits_t<VolumeType::Hexahedron,_degree>;		
  using tSubVolumesToNodes 	= std::array< std::array<unsigned int, traits_t::nnodesInVolume >, traits_t::NumSubElements >;
  using nodes_int_t = std::array< std::array<unsigned int, traits_t::dim >, traits_t::nnodes >;
  using nodes_real_t = std::array< std::array<double, traits_t::dim >, traits_t::nnodes >;
  
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
	
	static inline void ComputeCoordinates2(nodes_int_t& icoo_) noexcept
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
	
	
	static inline void ComputeCoordinates(nodes_int_t& icoo_) noexcept
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
	
	  static const unsigned int 	refhexa_icoo[] =  { 1,1,0,
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



      template <unsigned int _degree> struct treilli3d_utils<VolumeType::Tetrahedron,_degree>
      {
	using traits_t = treilli3d_traits_t<VolumeType::Tetrahedron,_degree>;		
	using tSubVolumesToNodes = std::array< std::array<unsigned int, traits_t::nnodesInVolume >, traits_t::NumSubElements >;
	using nodes_int_t = std::array< std::array<unsigned int, traits_t::dim >, traits_t::nnodes >;
	using nodes_real_t = std::array< std::array<double, traits_t::dim >, traits_t::nnodes >;
	using tTrianglesToNodes = std::array< std::array<unsigned int, 3 >, (_degree+1)*(_degree+1) >;

#if 0

	template<typename real_t> static void LobattoPoints(const std::array<real_t,_degree+1>& 	w_,														nodes_real_t& 		p_) noexcept
	{


	  
	  static constexpr const unsigned int n_ 		= _degree + 1;
	  static constexpr const unsigned int numPointsOnEdge 	= (_degree > 1) ? _degree - 1 : 0;
	  static constexpr const real_t zero(0.0);
	  static constexpr const real_t one(1.0);
	  static constexpr const real_t two(2.0);
	  static constexpr const real_t three(3.0);
	  static constexpr const real_t six(6.0);
	  
	  unsigned int startedge0 = 3;
	  unsigned int startedge1 = 3 + numPointsOnEdge + numPointsOnEdge-1;  
	  unsigned int startedge2 = 3 + numPointsOnEdge*2 + numPointsOnEdge-1;
	  unsigned int startInterior = 3*numPointsOnEdge+3;
	  
	  p_[0][0] 	= zero;
	  p_[0][1] 	= zero;
	  p_[1][0] 	= one;
	  p_[1][1] 	= zero;
	  p_[2][0] 	= zero;
	  p_[2][1] 	= one;

	  //
	  // Third edge
	  //
	  for (unsigned int j=1;j<_degree;++j)
	    {
	      const auto wj 	= w_[j];
	      const auto wk 	= w_[_degree-j];
	      p_[startedge2][0] 	= zero;
	      p_[startedge2][1] 	= (three + two * wj - wk) / six;
	      std::cout << p_[startedge2][1] << std::endl;
	      std::cout << (three + two * wj - wk) / six << std::endl;
	      --startedge2;
	    }
	  
	  //
	  // First edge
	  //
	  for (unsigned int i=1;i<_degree;++i)
	    {
	      const auto wi 	= w_[i];
	      const auto wk 	= w_[_degree - i];
	      p_[startedge0][0] 	= (three + two * wi - wk) / six;
	      p_[startedge0][1] 	= zero;
	      ++startedge0;
	    }
	  
	  startedge0=3;
	  // 
	  // Second edge
	  //
	  for (unsigned int i=1;i<_degree;++i)
	    {
	      p_[startedge1][0] 	= p_[startedge0++][0]; 
	      p_[startedge1][1] 	= p_[++startedge2][1];
	      --startedge1;
	    }
	  //
	  // Interior
	  //
	  for (unsigned int i=1;i<_degree;++i)
	    {
	      const auto wi = w_[i];
	      for (unsigned int j=1;j<_degree-i;++j)
		{
		  const auto wj 	= w_[j];
		  const auto wk 	= w_[_degree-i-j];
		  p_[startInterior][0] 	= ( two* (one + wi) - (wj + wk) ) / six;
		  p_[startInterior][1] 	= ( two* (one + wj) - (wi + wk) ) / six;
		  startInterior++;
		}
	    }
	  
	};
#endif
	
	
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
	
	static inline void ComputeCoordinates2(nodes_int_t& icoo_) noexcept
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
	  
	  { unsigned int N     = ((_degree+1)*(_degree+1+1))/2-(_degree+1);
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
	
	
	static inline void ComputeCoordinates(nodes_int_t& icoo_) noexcept
	{
	  unsigned int next = 0;
	  static constexpr const unsigned int  reftetra_icoo[] = {0,0,0,
						      1,0,0,
						      0,1,0,
						      0,0,1};

	  const unsigned int  	n 		= (_degree>0)? _degree-1 : 0;
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
	
	  static const unsigned int  tetraedge_cnc[] = { 1,2,
					     2,0,
					     0,1,
					     2,3,
					     0,3,
					     1,3};			
	
	  static const unsigned int  tetraface_cnc[] = {1,2,3,
					    2,0,3,
					    0,1,3,
					    0,2,1};
	  
	  { unsigned int  reedge[8];
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
	  { unsigned int  reface[32];
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




      template <unsigned int _degree> struct treilli3d_utils<VolumeType::Wedge,_degree>
      {
	using traits_t = treilli3d_traits_t<VolumeType::Wedge,_degree>;		
	using tSubVolumesToNodes = std::array< std::array<unsigned int, traits_t::nnodesInVolume >, traits_t::NumSubElements >;
	using nodes_int_t = std::array< std::array<unsigned int, traits_t::dim >, traits_t::nnodes >;
	using nodes_real_t = std::array< std::array<double, traits_t::dim >, traits_t::nnodes >;
	using tTrianglesToNodes = std::array< std::array<unsigned int, 3 >, (_degree+1)*(_degree+1) >;


	static inline void ComputeSubcnc(tSubVolumesToNodes& subCellsToNodes_) noexcept
	{
#define _dec(_i,_j,_k) (_k)* ( ((_degree+1)*(_degree+2)) /2) + ((  (_i)+(_j) + 1 )*( (_i)+(_j) ))/2 + (_i)
	  unsigned int subcell = 0;
	  unsigned int k;
	  for (unsigned int k=0;k<_degree;++k)
	    {
	      for (unsigned int j=0;j<_degree;j++)
		{
		  for (unsigned int i=0;i<j;i++)
		    {
		      subCellsToNodes_[subcell][0] 	= _dec(i,j-i,k);
		      subCellsToNodes_[subcell][1] 	= _dec(i+1,j-i,k);
		      subCellsToNodes_[subcell][2] 	= _dec(i,j+1-i,k); 
		      subCellsToNodes_[subcell][3] 	= _dec(i,j-i,k+1);
		      subCellsToNodes_[subcell][4] 	= _dec(i+1,j-i,k+1);
		      subCellsToNodes_[subcell][5] 	= _dec(i,j+1-i,k+1); 
		      ++subcell;
		      
		      subCellsToNodes_[subcell][0] 	= _dec(i,j-i,k);
		      subCellsToNodes_[subcell][1] 	= _dec(i+1,j-1-i,k);
		      subCellsToNodes_[subcell][2] 	= _dec(i+1,j-i,k);
		      subCellsToNodes_[subcell][3] 	= _dec(i,j-i,k+1);
		      subCellsToNodes_[subcell][4] 	= _dec(i+1,j-1-i,k+1);
		      subCellsToNodes_[subcell][5] 	= _dec(i+1,j-i,k+1);
		      ++subcell;
		    } 
		  
		  subCellsToNodes_[subcell][0]		=  _dec(j,0,k);
		  subCellsToNodes_[subcell][1]		=  _dec(j+1,0,k);
		  subCellsToNodes_[subcell][2]		=  _dec(j,1,k);	      
		  subCellsToNodes_[subcell][3]		=  _dec(j,0,k+1);
		  subCellsToNodes_[subcell][4]		=  _dec(j+1,0,k+1);
		  subCellsToNodes_[subcell][5]		=  _dec(j,1,k+1);	      
		  ++subcell;
		  
		}
	    }
#undef _dec

	};
	
	static inline void ComputeCoordinates2(nodes_int_t& icoo_) noexcept
	{
	  /*
	    6 
	    3 7
	    1 4 8 
	    0 2 5 9
	  */
	  unsigned int nodeIndex = 0;
	  for (unsigned int k=0;k<_degree+1;++k)
	    {
	      for (unsigned int i=0;i<_degree+1;i++)
		{
		  for (unsigned int j=0;j<=i;j++)
		    {
		      icoo_[nodeIndex][0] = j;
		      icoo_[nodeIndex][1] = i-j;
		      icoo_[nodeIndex][2] = k;
		      ++nodeIndex;
		    } 
		} 
	    } 
	  
	};

	static inline void ComputeCoordinates(nodes_int_t& icoo_) noexcept
	{
	  const unsigned int 	n 		= (_degree>0)? _degree-1 : 0;
	  unsigned int next = 0;
	  /* SOMMETS */
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
	  icoo_[next][0] = _degree;
	  icoo_[next][1] = 0;
	  icoo_[next][2] = _degree;
	  ++next;
	  icoo_[next][0] = 0;
	  icoo_[next][1] = _degree;
	  icoo_[next][2] = _degree;
	  ++next;

	  /* EDGE 0 */
	    for (unsigned int i=0;i<n;++i)
	      {
		icoo_[next][0] = (i+1);	  
		icoo_[next][1] = 0;	
		icoo_[next][2] = 0;	
		++next;	      
	      } 
	  /* EDGE 1 */
	    for (unsigned int i=0;i<n;++i)
	      {
		icoo_[next][0] = _degree-(i+1);	  
		icoo_[next][1] = (i+1);	
		icoo_[next][2] = 0;	
		++next;	      
	      } 
	  /* EDGE 2 */
	    for (unsigned int i=0;i<n;++i)
	      {
		icoo_[next][0] = 0;	  
		icoo_[next][1] = _degree-(i+1);
		icoo_[next][2] = 0;	
		++next;	      
	      } 


	  /* EDGE 3 */
	    for (unsigned int i=0;i<n;++i)
	      {
		icoo_[next][0] = (i+1);	  
		icoo_[next][1] = 0;	
		icoo_[next][2] = _degree;	
		++next;	      
	      } 
	  /* EDGE 4 */
	    for (unsigned int i=0;i<n;++i)
	      {
		icoo_[next][0] = _degree-(i+1);	  
		icoo_[next][1] = (i+1);	
		icoo_[next][2] = _degree;	
		++next;	      
	      } 
	  /* EDGE 5 */
	    for (unsigned int i=0;i<n;++i)
	      {
		icoo_[next][0] = 0;	  
		icoo_[next][1] = _degree-(i+1);
		icoo_[next][2] = _degree;	
		++next;	      
	      } 



	  /* EDGE 6 */
	    for (unsigned int i=0;i<n;++i)
	      {
		icoo_[next][0] = 0;	  
		icoo_[next][1] = 0;	
		icoo_[next][2] = (i+1);	
		++next;	      
	      } 
	  /* EDGE 7 */
	    for (unsigned int i=0;i<n;++i)
	      {
		icoo_[next][0] = _degree;	  
		icoo_[next][1] = 0;	
		icoo_[next][2] = (i+1);	
		++next;	      
	      } 
	  /* EDGE 8 */
	    for (unsigned int i=0;i<n;++i)
	      {
		icoo_[next][0] = 0;	  
		icoo_[next][1] = _degree;
		icoo_[next][2] = (i+1);	
		++next;	      
	      } 

#if 0
	  static const I wedgetriaface_cnc[] = {0,1,2,
						3,5,4};
	
#endif
	  static const unsigned int wedgequadface_cnc[] = {0,2,5,3,
							   4,5,2,1,
							   0,3,4,1};
	  
	  static const unsigned int 	refwedge_icoo[] =  { 0,0,0,
							     1,0,0,
							     0,1,0,
							     0,0,1,
							     1,0,1,
							     0,1,1};
	  { unsigned int reface[32];
	    { unsigned int iface;
	      for (unsigned int iface=0;iface<3;++iface)
		{
		  reface[0] = refwedge_icoo[3*wedgequadface_cnc[iface * 4 + 0]+0];
		  reface[1] = refwedge_icoo[3*wedgequadface_cnc[iface * 4 + 0]+1];
		  reface[2] = refwedge_icoo[3*wedgequadface_cnc[iface * 4 + 0]+2];
		  
		  reface[3] = refwedge_icoo[3*wedgequadface_cnc[iface * 4 + 1]+0];
		  reface[4] = refwedge_icoo[3*wedgequadface_cnc[iface * 4 + 1]+1];
		  reface[5] = refwedge_icoo[3*wedgequadface_cnc[iface * 4 + 1]+2];
		  
		  reface[6] = refwedge_icoo[3*wedgequadface_cnc[iface * 4 + 2]+0];
		  reface[7] = refwedge_icoo[3*wedgequadface_cnc[iface * 4 + 2]+1];
		  reface[8] = refwedge_icoo[3*wedgequadface_cnc[iface * 4 + 2]+2];
		
		  reface[9]  = refwedge_icoo[3*wedgequadface_cnc[iface * 4 + 3]+0];
		  reface[10] = refwedge_icoo[3*wedgequadface_cnc[iface * 4 + 3]+1];
		  reface[11] = refwedge_icoo[3*wedgequadface_cnc[iface * 4 + 3]+2];

		  { unsigned int i;
		    for (unsigned int i=0;i<n;++i)
		      {
			{ unsigned int j;
			  for (unsigned int j=0;j<n;++j)
			    {
			      const unsigned int l0 			= (i+1)*(j+1);
			      const unsigned int l1 			= (_degree-(i+1))*(j+1);
			      const unsigned int l2 			= (_degree-(i+1))*(_degree-(j+1));
			      const unsigned int l3 			= (_degree-(j+1))*(i+1);
			    
			      const unsigned int xi 			= reface[0] * l0 + reface[3] * l1 + reface[6] * l2 + reface[9] *l3;
			      const unsigned int xj 			= reface[1] * l0 + reface[4] * l1 + reface[7] * l2 + reface[10]*l3;
			      const unsigned int xk 			= reface[2] * l0 + reface[5] * l1 + reface[8] * l2 + reface[11]*l3;
			    
			      icoo_[next][0] 	= xi/_degree;	  
			      icoo_[next][1] 	= xj/_degree;	
			      icoo_[next][2] 	= xk/_degree;	
			      ++next;
			    } }
		      } } 

		} } }

	  { unsigned int iface;
	    for (unsigned int iface=0;iface<2;++iface)
	      {
		{ unsigned int i;
		  for (unsigned int i=0;i<n;++i)
		    {
		      { unsigned int j;
			for (unsigned int j=0;j<n-(i+1);++j)
			  {
			    icoo_[next][0] = (i+1);	  
			    icoo_[next][1] = (j+1);
			    icoo_[next][2] = iface*_degree;
			    ++next;	      
			  } }
		    } }
	      
	      } } 

	  /* interior */
	  { unsigned int k;
	    for (unsigned int k=0;k<n;++k)
	      {
		{ unsigned int i;
		  for (unsigned int i=0;i<n;++i)
		    {
		      { unsigned int j;
			for (unsigned int j=0;j<n-(i+1);++j)
			  {
			    icoo_[next][0] = (i+1);	  
			    icoo_[next][1] = (j+1);
			    icoo_[next][2] = (k+1);
			    ++next;	      
			  } }
		    } }
	      } }




	  
	  
	};
	  
      };








template <typename impl_t, VolumeType::enum_t _volumeShape, unsigned int _degree>
class treilli3d_t : public CRTP< treilli3d_t<impl_t,_volumeShape ,_degree> >
{
protected:  
  using traits_t 		= treilli3d_traits_t<_volumeShape,_degree>;		
  using volumes_to_nodes_t 	= std::array< std::array<unsigned int, traits_t::nnodesInVolume >, traits_t::NumSubElements >;
  using nodes_int_t 		= std::array< std::array<unsigned int, traits_t::dim >, traits_t::nnodes >;
  using nodes_real_t 		= std::array< std::array<double, traits_t::dim >, traits_t::nnodes >;
  
  static constexpr const unsigned int s_numNodes 		= traits_t::nnodes;
  static constexpr const unsigned int s_numElements 		= traits_t::NumSubElements;
  static constexpr const unsigned int s_dimension 		= traits_t::dim;
  static constexpr const unsigned int s_degree 			= traits_t::Degree;
  static constexpr const unsigned int s_numNodesInVolume	= traits_t::nnodesInVolume;

  volumes_to_nodes_t 	m_subcnc;
  
public:

  inline unsigned int 		degree() 	const noexcept { return s_degree; };
  inline unsigned int 		dim() 		const noexcept { return s_dimension; };
  inline unsigned int 		nnodes() 	const noexcept { return s_numNodes; };    
  inline unsigned int 		nsubcells() 	const noexcept { return s_numElements; };
  inline VolumeType::enum_t 	shape()		const noexcept { return _volumeShape; };
  inline unsigned int 		nnodesincell() 	const noexcept { return s_numNodesInVolume; };    
  
  template <typename _float_type> inline	_float_type GetCoordinate(const unsigned int&nodeIndex_,
									  const unsigned int&dimensionIndex_) const
  {
    return static_cast<const impl_t&>(*this).GetCoordinate<_float_type>(nodeIndex_,dimensionIndex_);
  };
  
  inline unsigned int GetNodeIndex(const unsigned int&subElementIndex_,
				   const unsigned int&localNodeIndex_) const
  {
    return m_subcnc[subElementIndex_][localNodeIndex_];
  };
  
  inline treilli3d_t(nodes_int_t& icoo)
  {		
    
    // COMPUTE COORDINATES
    treilli3d_utils<_volumeShape,_degree>::ComputeCoordinates(icoo);
    {
      std::array<unsigned int,(_degree+1)*(_degree+1)*(_degree+1)> perm;
      nodes_int_t ilagr;
      
      treilli3d_utils<_volumeShape,_degree>::ComputeCoordinates2(ilagr);

      for (unsigned int i=0;i<s_numNodes;++i)
	{	
	  perm[(_degree+1)*(_degree+1)*icoo[i][2] + icoo[i][0]*(_degree+1) + icoo[i][1]] = 1+i;
	} 
	  
      treilli3d_utils<_volumeShape,_degree>::ComputeSubcnc(this->m_subcnc); 
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
      
  inline ~treilli3d_t()
  {
  };


protected:

      
};



template <VolumeType::enum_t _volumeShape,unsigned int _degree>
class Uniform3d : public treilli3d_t<Uniform3d<_volumeShape,_degree> ,_volumeShape,_degree>
{  
private: using traits_t = treilli3d_traits_t<_volumeShape,_degree>;
protected: using nodes_int_t = std::array< std::array<unsigned int, traits_t::dim >, traits_t::nnodes >;
private: nodes_int_t m_icoo;
  
public: template <typename _float_type> inline _float_type GetCoordinate(const unsigned int nodeIndex_,
									 const unsigned int dimensionIndex_) const noexcept
  {
    // static const _float_type idegree = _float_type(1.0) / _float_type(_degree);
    // return _float_type(this->m_icoo[nodeIndex_][dimensionIndex_]) * idegree;
    static constexpr const _float_type idegree = _float_type(1.0) / _float_type(_degree);
    //#ifndef NDEBUG
    //	  Debug::IsInRange(__TRACE__,nodeIndex_,(unsigned int)0,this->s_numNodes-1);
    //	  Debug::IsInRange(__TRACE__,dimensionIndex_,(unsigned int)0,this->s_dimension-1);
    //#endif      
    return _float_type(this->m_icoo[nodeIndex_][dimensionIndex_]) * idegree;
  };
            
public: inline Uniform3d() noexcept : treilli3d_t< Uniform3d<_volumeShape,_degree>, _volumeShape,_degree>(this->m_icoo)
  {  	
  };
	
  inline ~Uniform3d() noexcept
  {
  };
    
};

namespace FiniteElement
{
  template <typename impl_t,
	    FaceType::enum_t 	_faceShape,
	    unsigned int 	_degree>
  class treilli2d_t;

  
  template <FaceType::enum_t _faceShape,unsigned int _degree> struct treilli2d_traits_t
  {
  };
  template <unsigned int _degree> struct treilli2d_traits_t<FaceType::Triangle,_degree>
  {
    static constexpr const unsigned int nnodes 		= ( (_degree+1)*(_degree+2))/2;
    static constexpr const unsigned int nsubcells 	= _degree*_degree;
    static constexpr const unsigned int dim	= 2;    
    static constexpr const unsigned int degree 		= _degree;    
    static constexpr const unsigned int nnodesInFace	= 3;
  };
  
  template <unsigned int _degree> struct treilli2d_traits_t<FaceType::Quadrilateral,_degree>
  {
    static constexpr const unsigned int nnodes 		= (_degree+1)*(_degree+1);
    static constexpr const unsigned int nsubcells 	= _degree*_degree;
    static constexpr const unsigned int dim		= 2;    
    static constexpr const unsigned int degree 		= _degree;    
    static constexpr const unsigned int nnodesInFace	= 4;
  };

    
  template <FaceType::enum_t _faceShape,unsigned int _degree> struct treilli2d_utils
  {	
  };
  
  template <unsigned int _degree> struct treilli2d_utils<FaceType::Triangle,_degree>
  {
    using traits_t = treilli2d_traits_t<FaceType::Triangle,_degree>;
    
    using faces_to_nodes_t 		= std::array< std::array<unsigned int, traits_t::nnodesInFace >, traits_t::nsubcells >;
    using nodes_int_t 			= std::array< std::array<unsigned int, traits_t::dim >, traits_t::nnodes >;
    using nodes_real_t 	= std::array< std::array<double, traits_t::dim >, traits_t::nnodes >;
    
    template<typename real_t> static void LobattoPoints(const std::array<real_t,_degree+1>& 	w_,														nodes_real_t& 		p_) noexcept
    {
      
      static constexpr const unsigned int n_ 		= _degree + 1;
	  static constexpr const unsigned int numPointsOnEdge 	= (_degree > 1) ? _degree - 1 : 0;
	  static constexpr const real_t zero(0.0);
	  static constexpr const real_t one(1.0);
	  static constexpr const real_t two(2.0);
	  static constexpr const real_t three(3.0);
	  static constexpr const real_t six(6.0);
	  
	  unsigned int startedge0 = 3;
	  unsigned int startedge1 = 3 + numPointsOnEdge + numPointsOnEdge-1;  
	  unsigned int startedge2 = 3 + numPointsOnEdge*2 + numPointsOnEdge-1;
	  unsigned int startInterior = 3*numPointsOnEdge+3;
	  
	  p_[0][0] 	= zero;
	  p_[0][1] 	= zero;
	  p_[1][0] 	= one;
	  p_[1][1] 	= zero;
	  p_[2][0] 	= zero;
	  p_[2][1] 	= one;

	  //
	  // Third edge
	  //
	  for (unsigned int j=1;j<_degree;++j)
	    {
	      const auto wj 	= w_[j];
	      const auto wk 	= w_[_degree-j];
	      p_[startedge2][0] 	= zero;
	      p_[startedge2][1] 	= (three + two * wj - wk) / six;
	  std::cout << p_[startedge2][1] << std::endl;
	  std::cout << (three + two * wj - wk) / six << std::endl;
	      --startedge2;
	    }
	  
	  //
	  // First edge
	  //
	  for (unsigned int i=1;i<_degree;++i)
	    {
	      const auto wi 	= w_[i];
	      const auto wk 	= w_[_degree - i];
	      p_[startedge0][0] 	= (three + two * wi - wk) / six;
	      p_[startedge0][1] 	= zero;
	      ++startedge0;
	    }
	  
	  startedge0=3;
	  // 
	  // Second edge
	  //
	  for (unsigned int i=1;i<_degree;++i)
	    {
	      p_[startedge1][0] 	= p_[startedge0++][0]; 
	      p_[startedge1][1] 	= p_[++startedge2][1];
	      --startedge1;
	    }
	  //
	  // Interior
	  //
	  for (unsigned int i=1;i<_degree;++i)
	    {
	      const auto wi = w_[i];
	      for (unsigned int j=1;j<_degree-i;++j)
		{
		  const auto wj 	= w_[j];
		  const auto wk 	= w_[_degree-i-j];
		  p_[startInterior][0] 	= ( two* (one + wi) - (wj + wk) ) / six;
		  p_[startInterior][1] 	= ( two* (one + wj) - (wi + wk) ) / six;
		  startInterior++;
		}
	    }
	  
	};

	
	static inline void ComputeSubcnc(faces_to_nodes_t & subcnc_)
	{
	  unsigned int subCellIndex = 0;
#define _dec(_i,_j) (( (_i)+(_j) + 1 )*( (_i)+(_j) ))/2+(_i)
	  for (unsigned int j=0;j<_degree;j++)
	    {
	    
	      for (unsigned int i=0;i<j;i++)
		{
		 
		  subcnc_[subCellIndex][0] = _dec(i,j-i);
		  subcnc_[subCellIndex][1] = _dec(i+1,j-i);
		  subcnc_[subCellIndex][2] = _dec(i,j+1-i); 
		  ++subCellIndex;
		
		  subcnc_[subCellIndex][0] = _dec(i,j-i);
		  subcnc_[subCellIndex][1] = _dec(i+1,j-1-i);
		  subcnc_[subCellIndex][2] = _dec(i+1,j-i);
		  ++subCellIndex;

		}
	    
	      subcnc_[subCellIndex][0] = _dec(j,0);
	      subcnc_[subCellIndex][1] = _dec(j+1,0);
	      subcnc_[subCellIndex][2] = _dec(j,1);
	      ++subCellIndex;
	    }
#undef _dec  
	};
      

	static inline void ComputeCoordinates2(nodes_int_t& icoo_)
	{
	  // COMPUTE GRID     
	  //	6 
	  //	3 7
	  //	1 4 8 
	  //	0 2 5 9
	  unsigned int nodeIndex = 0;
	  for (unsigned int i=0;i<_degree+1;i++)
	    {
	      for (unsigned int j=0;j<=i;j++)
		{
		  icoo_[nodeIndex][0] = j;
		  icoo_[nodeIndex][1] = i-j;
		  ++nodeIndex;
		}
	    }
	};

	
	static inline void ComputeCoordinates(nodes_int_t& icoo_)
	{
	  unsigned int nodeIndex = 0;
	  // VERTEX 0 
	  icoo_[nodeIndex][0] = 0;
	  icoo_[nodeIndex][1] = 0;
	  ++nodeIndex;
	
	  // VERTEX 1 
	  icoo_[nodeIndex][0] = _degree;
	  icoo_[nodeIndex][1] = 0;
	  ++nodeIndex;

	  // VERTEX 2
	  icoo_[nodeIndex][0] = 0;
	  icoo_[nodeIndex][1] = _degree;
	  ++nodeIndex;
	    
	  // EDGE 0 
	  for (unsigned int i=0;i<_degree - 1;++i)
	    {
	      icoo_[nodeIndex][0] = i+1;
	      icoo_[nodeIndex][1] = 0;
	      ++nodeIndex;
	    }
	  // EDGE 1 
	  for (unsigned int i=0;i<_degree-1;++i)
	    {
	      icoo_[nodeIndex][0] = _degree-(i+1); 
	      icoo_[nodeIndex][1] = i+1;
	      ++nodeIndex;
	    }
	  // EDGE 2 
	  for (unsigned int i=0;i<_degree-1;++i)
	    {
	      icoo_[nodeIndex][0] = 0;	  
	      icoo_[nodeIndex][1] = _degree-(i+1);
	      ++nodeIndex;
	    }
	  // INTERIOR
	  for (unsigned int i=0;i<_degree - 1;++i)
	    {
	      for (unsigned int j=0;j<_degree-i-2;++j)
		{
		  icoo_[nodeIndex][0] = i + 1;
		  icoo_[nodeIndex][1] = j + 1;
		  ++nodeIndex;
		}
	    }
	};
	  
      };


      template <unsigned int _degree> struct treilli2d_utils<FaceType::Quadrilateral,_degree>
      {
	using traits_t = treilli2d_traits_t<FaceType::Quadrilateral,_degree>;		
	using faces_to_nodes_t = std::array< std::array<unsigned int, traits_t::nnodesInFace >, traits_t::nsubcells >;
	using nodes_int_t = std::array< std::array<unsigned int, traits_t::dim >, traits_t::nnodes >;
	using nodes_real_t = std::array< std::array<double, traits_t::dim >, traits_t::nnodes >;

	template<typename real_t> static void LobattoPoints(const std::array<real_t,_degree+1>& w_,													    nodes_real_t& 	p_) noexcept
	{
	  static constexpr const unsigned int n_ 		= _degree + 1;
	  static constexpr const unsigned int numPointsOnEdge 	= (_degree > 1) ? _degree - 1 : 0;
	  static constexpr const real_t mone(-1.0);
	  static constexpr const real_t one(1.0);
	  static constexpr const real_t two(2.0);
	  static constexpr const real_t three(3.0);
	  static constexpr const real_t six(6.0);
	  
	  p_[0][0] 	= mone;
	  p_[0][1] 	= mone;
	  p_[1][0] 	= one;
	  p_[1][1] 	= mone;
	  p_[2][0] 	= one;
	  p_[2][1] 	= one;
	  p_[3][0] 	= mone;
	  p_[3][1] 	= one;

	  unsigned int pointIndex = 4;
	  
	  //
	  // First edge
	  //
	  for (unsigned int i=1;i<_degree;++i)
	    {
	      p_[pointIndex][0] 	= w_[i];
	      p_[pointIndex][1] 	= mone;
	      ++pointIndex;
	    }

	  for (unsigned int i=1;i<_degree;++i)
	    {
	      p_[pointIndex][0] 	= one;
	      p_[pointIndex][1] 	= w_[i];
	      ++pointIndex;
	    }

	  for (unsigned int i=1;i<_degree;++i)
	    {
	      p_[pointIndex][0] 	= w_[_degree - i];
	      p_[pointIndex][1] 	= one;
	      ++pointIndex;
	    }

	  for (unsigned int i=1;i<_degree;++i)
	    {
	      p_[pointIndex][0] 	= mone;
	      p_[pointIndex][1] 	= w_[_degree - i];
	      ++pointIndex;
	    }
	  
	  for (unsigned int i=1;i<_degree;++i)
	    {
	      for (unsigned int j=1;j<_degree;j++)
		{			  
		  p_[pointIndex][0] = w_[i];
		  p_[pointIndex][1] = w_[j];	
		  ++pointIndex;
		} 
	    }
	};

	static inline void ComputeSubcnc(faces_to_nodes_t& subcnc_)
	{
	  unsigned int subCellIndex = 0;
	  // COMPUTE CNC
#define _dec(_i,_j)   (_degree+1) * (_i) + (_j) 
	  { 
	    for (unsigned int i=0;i<_degree;i++)
	      {
		for (unsigned int j=0;j<_degree;j++)
		  {		    
		    subcnc_[subCellIndex][0] = _dec( (i+1), (j+1) );
		    subcnc_[subCellIndex][1] = _dec( (i), (j+1) );
		    subcnc_[subCellIndex][2] = _dec( (i), (j) );
		    subcnc_[subCellIndex][3] = _dec( (i+1), (j) );
		    ++subCellIndex;
		  } 
	      }
	  }
#undef _dec
	};
      

	static inline void ComputeCoordinates2(nodes_int_t& icoo_)
	{
	  for (unsigned int i=0;i<_degree+1;++i)
	    {
	      for (unsigned int j=0;j<_degree+1;j++)
		{
		  icoo_[i*(_degree+1)+j][0] = i;	  
		  icoo_[i*(_degree+1)+j][1] = j;	
		}
	    }
	};


	static inline void ComputeCoordinates(nodes_int_t& icoo_)
	{
	  unsigned int nodeIndex = 0;
	
	  icoo_[nodeIndex][0] = 0;	  
	  icoo_[nodeIndex][1] = 0;	
	  ++nodeIndex;
	
	  icoo_[nodeIndex][0] = _degree;	  
	  icoo_[nodeIndex][1] = 0;	
	  ++nodeIndex;
		
	  icoo_[nodeIndex][0] = _degree;	  
	  icoo_[nodeIndex][1] = _degree;	
	  ++nodeIndex;
		
	  icoo_[nodeIndex][0] = 0;	  
	  icoo_[nodeIndex][1] = _degree;	
	  ++nodeIndex;
		
	  /* EDGE 0 */
	  for (unsigned int i=1;i<_degree;++i)
	    {
	      icoo_[nodeIndex][0] = i;	  
	      icoo_[nodeIndex][1] = 0;
	      ++nodeIndex;
	    }
	
	  /* EDGE 1 */
	  for (unsigned int i=1;i<_degree;++i)
	    {
	      icoo_[nodeIndex][0] = _degree;	  
	      icoo_[nodeIndex][1] = i;
	      ++nodeIndex;
	    } 
	
	  /* EDGE 2 */
	  for (unsigned int i=0;i<_degree-1;++i)
	    {
	      icoo_[nodeIndex][0] = _degree - i - 1;	  
	      icoo_[nodeIndex][1] = _degree;	
	      ++nodeIndex;
	    }

	  /* EDGE 3 */
	  for (unsigned int i=0;i<_degree-1;++i)
	    {
	      icoo_[nodeIndex][0] = 0;	  
	      icoo_[nodeIndex][1] = _degree - i - 1;	
	      ++nodeIndex;
	    }
	
	  /* INSIDE */
	  for (unsigned int i=1;i<_degree;++i)
	    {
	      for (unsigned int j=1;j<_degree;j++)
		{			  
		  icoo_[nodeIndex][0] = i;
		  icoo_[nodeIndex][1] = j;	
		  ++nodeIndex;
		} 
	    }
	  
	};

      };

      template <typename impl_t, FaceType::enum_t _faceShape,unsigned int _degree>
      class treilli2d_t : public CRTP<treilli2d_t<impl_t,_faceShape,_degree> > 
      {
      private: using traits_t 	= treilli2d_traits_t<_faceShape,_degree>;
	
      private: using faces_to_nodes_t = std::array< std::array<unsigned int, traits_t::nnodesInFace >, traits_t::nsubcells >;
      protected: using nodes_int_t = std::array< std::array<unsigned int, traits_t::dim >, traits_t::nnodes >;
	
      protected: static constexpr const unsigned int s_numNodes 	= traits_t::nnodes;
      protected: static constexpr const unsigned int s_numElements 	= traits_t::nsubcells;
      protected: static constexpr const unsigned int s_dimension 	= traits_t::dim;
      protected: static constexpr const unsigned int s_degree 		= traits_t::degree;
      protected: static constexpr const unsigned int s_numNodesInFace	= traits_t::nnodesInFace;
	
      public: inline unsigned int 	degree() 		const noexcept { return s_degree; };
      public: inline unsigned int 	dimension() 		const noexcept { return s_dimension; };
      public: inline unsigned int 	nnodes() 		const noexcept { return s_numNodes; };    
      public: inline unsigned int 	nsubcells() 		const noexcept { return s_numElements; };
      public: inline FaceType::enum_t	celltype()		const noexcept { return _faceShape; };
      public: inline unsigned int 	nnodesincell() 		const noexcept { return s_numNodesInFace; };    
	
      protected: faces_to_nodes_t m_subcnc;

      public: template <typename _float_type>
      inline _float_type GetCoordinate(const unsigned int nodeIndex_,
				       const unsigned int dimensionIndex_) const noexcept
	{
	  return static_cast<const impl_t&>(*this).GetCoordinate<_float_type>(nodeIndex_,dimensionIndex_);
	};
	
      public: inline unsigned int GetNodeIndex(const unsigned int&subElementIndex_,
					       const unsigned int&localNodeIndex_) const noexcept
	{
//#ifndef NDEBUG
//	  Debug::IsInRange(__TRACE__,subElementIndex_,(unsigned int)0,this->s_numElements-1);
//	  Debug::IsInRange(__TRACE__,localNodeIndex_,(unsigned int)0,this->s_numNodesInFace-1);
//#endif
	  return this->m_subcnc[subElementIndex_][localNodeIndex_];
	};
	
      protected: inline treilli2d_t(nodes_int_t& icoo) noexcept
	{
	  //
	  // Compute integer coordinates.
	  //
	  treilli2d_utils<_faceShape,_degree>::ComputeCoordinates(icoo);
	  std::array<unsigned int, (_degree+1)*(_degree+1)> perm;	  
	  nodes_int_t  ilagr;
	  {
	  //
	  // Compute integer coordinates.
	  //
	    treilli2d_utils<_faceShape,_degree>::ComputeCoordinates2(ilagr);	  
	    for (unsigned int i=0;i<s_numNodes;++i)
	      {	
		perm[ icoo[i][0] * (_degree+1) + icoo[i][1] ] = i+1;
	      } 
	  
	    treilli2d_utils<_faceShape,_degree>::ComputeSubcnc(this->m_subcnc);
	  
	    for (unsigned int subElementIndex=0;subElementIndex<s_numElements;++subElementIndex)
	      {
		for (unsigned int iv=0;iv<s_numNodesInFace;++iv)
		  {	  
		    const unsigned int l 				= m_subcnc[subElementIndex][iv];
		    const unsigned int i 				= ilagr[l][0];
		    const unsigned int j 				= ilagr[l][1];
		    m_subcnc[subElementIndex][iv] = perm[(_degree+1)*i+j]-1;
		  }
	      }
	  }
	
	};
    
	inline ~treilli2d_t() noexcept
	{
	};
    
      };
    

      template <FaceType::enum_t _faceShape,unsigned int _degree>
      class Uniform : public treilli2d_t<Uniform<_faceShape,_degree> ,_faceShape,_degree>
      {
      
      private: using traits_t = treilli2d_traits_t<_faceShape,_degree>;
      protected: using nodes_int_t = std::array< std::array<unsigned int, traits_t::dim >, traits_t::nnodes >;
      private: nodes_int_t m_icoo;

      public: template <typename _float_type> inline _float_type GetCoordinate(const unsigned int nodeIndex_,
									       const unsigned int dimensionIndex_) const noexcept
	{
	  static constexpr const _float_type idegree = _float_type(1.0) / _float_type(_degree);
//#ifndef NDEBUG
//	  Debug::IsInRange(__TRACE__,nodeIndex_,(unsigned int)0,this->s_numNodes-1);
//	  Debug::IsInRange(__TRACE__,dimensionIndex_,(unsigned int)0,this->s_dimension-1);
//#endif      
	  return _float_type(this->m_icoo[nodeIndex_][dimensionIndex_]) * idegree;
	};
            
      public: inline Uniform() noexcept : treilli2d_t<Uniform<_faceShape,_degree>,_faceShape,_degree>(this->m_icoo)
	{  	
	};
	
	inline ~Uniform() noexcept
	{
	};
    
      };


  template <FaceType::enum_t _faceShape,unsigned int _degree>
  class Generator : public treilli2d_t<Generator<_faceShape,_degree> ,_faceShape,_degree>
  {
    
  private: using traits_t = treilli2d_traits_t<_faceShape,_degree>;
  protected: using nodes_int_t = std::array< std::array<unsigned int, traits_t::dim >, traits_t::nnodes >;
  protected: using nodes_real_t = std::array< std::array<double, traits_t::dim >, traits_t::nnodes >;

  private: nodes_int_t m_icoo;
  private: nodes_real_t m_rcoo;
    
    //!
    //! @brief Get the coordinates.
    //!
  public: template <typename _float_type>
  inline _float_type GetCoordinate(const unsigned int nodeIndex_,
				   const unsigned int dimensionIndex_) const noexcept
    {
      //#ifndef NDEBUG
      //	  Debug::IsInRange(__TRACE__,nodeIndex_,(unsigned int)0,this->s_numNodes-1);
      //	  Debug::IsInRange(__TRACE__,dimensionIndex_,(unsigned int)0,this->s_dimension-1);
      //#endif      
      return this->m_rcoo[nodeIndex_][dimensionIndex_];
    };
    
    
  public: inline Generator(const double * p) noexcept : treilli2d_t<Generator<_faceShape,_degree>,_faceShape,_degree>(this->m_icoo)
    {
      //      Quadrature::Edge::Legendre<double,_degree-1> l;
      std::array<double,_degree+1> w;
      w[0]=-1.0;
      for (int i=0;i<_degree-1;++i)
	{
	  std::cout << "eee" << std::endl;
	  w[1+i] = p[i]; // l.GetPosition(i,0);
	  std::cout << "eee done" << std::endl;
	}	  
      w[_degree]=1.0;
      treilli2d_utils<_faceShape,_degree>::LobattoPoints(w, 
						       this->m_rcoo);
      
    };
    
    inline ~Generator() noexcept
    {
    };
  };
  
#if 0
  
  template <typename impl_t, FaceType::enum_t _faceShape,unsigned int _degree>
  struct Traits_CRTP<treilli2d_t<impl_t,_faceShape,_degree> >
  {
    typedef FaceType::enum_t eTypeShape;
  };
#endif
  
#if 0


      



  


      
    };
  
  };
#endif
};
