#pragma once
#include "Quadrature/Edge/Legendre.hpp"

#include "CRTP.hpp"
#include "Mesh/Topology/eFaceShape.h"
#include <array>
namespace FiniteElement
{
  template <typename impl_t,
	    FaceType::enum_t 	_faceShape,
	    unsigned int 	_degree>
  class treilli_t;
  
  template <typename impl_t, FaceType::enum_t _faceShape,unsigned int _degree>
  struct Traits_CRTP<treilli_t<impl_t,_faceShape,_degree> >
  {
    typedef FaceType::enum_t eTypeShape;
  };
  
  template <FaceType::enum_t _faceShape,unsigned int _degree> struct treilli_traits_t
  {
  };
  
  template <unsigned int _degree> struct treilli_traits_t<FaceType::Triangle,_degree>
  {
    static constexpr const unsigned int nnodes 		= ( (_degree+1)*(_degree+2))/2;
    static constexpr const unsigned int nsubcells 	= _degree*_degree;
    static constexpr const unsigned int dim	= 2;    
    static constexpr const unsigned int degree 		= _degree;    
    static constexpr const unsigned int nnodesInFace	= 3;
  };
  
  template <unsigned int _degree> struct treilli_traits_t<Mesh::Topology::Quadrilateral,_degree>
  {
    static constexpr const unsigned int nnodes 		= (_degree+1)*(_degree+1);
    static constexpr const unsigned int nsubcells 	= _degree*_degree;
    static constexpr const unsigned int dim		= 2;    
    static constexpr const unsigned int degree 		= _degree;    
    static constexpr const unsigned int nnodesInFace	= 4;
  };
  
  template <FaceType::enum_t _faceShape,unsigned int _degree> struct Utils
  {	
  };
    
  template <unsigned int _degree> struct Utils<Mesh::Topology::TRIANGLE,_degree>
  {
    using traits_t 			= treilli_traits_t<Mesh::Topology::TRIANGLE,_degree>;		
    using tSubFacesToNodes 		= std::array< std::array<unsigned int, traits_t::nnodesInFace >, traits_t::nsubcells >;
    using nodes_int_t 		= std::array< std::array<unsigned int, traits_t::dim >, traits_t::nnodes >;
    using tLobattoNodesCoordinates 	= std::array< std::array<double, traits_t::dim >, traits_t::nnodes >;
    
    template<typename real_t> static void LobattoPoints(const std::array<real_t,_degree+1>& w_,											 
							tLobattoNodesCoordinates& 	p_) noexcept
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

	
	static inline void ComputeSubcnc(tSubFacesToNodes & subcnc_)
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


      template <unsigned int _degree> struct Utils<FaceType::Quadrilateral,_degree>
      {
	using traits_t = treilli_traits_t<FaceType::Quadrilateral,_degree>;		
	using tSubFacesToNodes = std::array< std::array<unsigned int, traits_t::nnodesInFace >, traits_t::nsubcells >;
	using nodes_int_t = std::array< std::array<unsigned int, traits_t::dim >, traits_t::nnodes >;
	using tLobattoNodesCoordinates = std::array< std::array<double, traits_t::dim >, traits_t::nnodes >;

	template<typename real_t> static void LobattoPoints(const std::array<real_t,_degree+1>& w_,													    tLobattoNodesCoordinates& 	p_) noexcept
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

	static inline void ComputeSubcnc(tSubFacesToNodes& subcnc_)
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
      class treilli_t : public CRTP<treilli_t<impl_t,_faceShape,_degree> > 
      {
      private: using Impl = impl_t;
      private: using BaseClass = CRTP<treilli_t<impl_t,_faceShape,_degree> >;
      private: using eTypeShape = typename BaseClass::eTypeShape;      
      private: using traits_t = treilli_traits_t<_faceShape,_degree>;
      private: using tSubFacesToNodes = std::array< std::array<unsigned int, traits_t::nnodesInFace >, traits_t::nsubcells >;
      protected: using nodes_int_t = std::array< std::array<unsigned int, traits_t::dim >, traits_t::nnodes >;
      
      protected: static constexpr const unsigned int s_numNodes 	= traits_t::nnodes;
      protected: static constexpr const unsigned int s_numElements 	= traits_t::nsubcells;
      protected: static constexpr const unsigned int s_dimension 	= traits_t::dim;
      protected: static constexpr const unsigned int s_degree 		= traits_t::degree;
      protected: static constexpr const unsigned int s_numNodesInFace	= traits_t::nnodesInFace;

      protected: tSubFacesToNodes m_subcnc;
	
      public: inline unsigned int 	GetDegree() 		const noexcept { return s_degree; };
      public: inline unsigned int 	GetDimension() 		const noexcept { return s_dimension; };
      public: inline unsigned int 	GetNodes() 		const noexcept { return s_numNodes; };    
      public: inline unsigned int 	GetNumSubElements() 	const noexcept { return s_numElements; };
      public: inline FaceType::enum_t	GetShape()		const noexcept { return _faceShape; };
      public: inline unsigned int 	GetnnodesInCell() 	const noexcept { return s_numNodesInFace; };    

      public: template <typename _float_type> inline _float_type GetCoordinate(const unsigned int nodeIndex_,
									       const unsigned int dimensionIndex_) const noexcept
	{
	  return static_cast<const impl_t&>(*this).GetCoordinate<_float_type>(nodeIndex_,dimensionIndex_);
	};
	
      public: inline unsigned int GetNodeIndex(const unsigned int&subElementIndex_,
					       const unsigned int&localNodeIndex_) const noexcept
	{
#ifndef NDEBUG
	  Debug::IsInRange(__TRACE__,subElementIndex_,(unsigned int)0,this->s_numElements-1);
	  Debug::IsInRange(__TRACE__,localNodeIndex_,(unsigned int)0,this->s_numNodesInFace-1);
#endif
	  return this->m_subcnc[subElementIndex_][localNodeIndex_];
	};

	
      protected: inline treilli_t(nodes_int_t& icoo) noexcept
	{
	  //
	  // Compute integer coordinates.
	  //
	  Utils<_faceShape,_degree>::ComputeCoordinates(icoo);
	  std::array<unsigned int, (_degree+1)*(_degree+1)> perm;	  
	  nodes_int_t  ilagr;
	  {
	  //
	  // Compute integer coordinates.
	  //
	    Utils<_faceShape,_degree>::ComputeCoordinates2(ilagr);	  
	    for (unsigned int i=0;i<s_numNodes;++i)
	      {	
		perm[ icoo[i][0] * (_degree+1) + icoo[i][1] ] = i+1;
	      } 
	  
	    Utils<_faceShape,_degree>::ComputeSubcnc(this->m_subcnc);
	  
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
    
	inline ~treilli_t() noexcept
	{
	};
    
      };
    

      




  
    
      template <FaceType::enum_t _faceShape,unsigned int _degree>
      class Uniform : public treilli_t<Uniform<_faceShape,_degree> ,_faceShape,_degree>
      {
      
      private: using traits_t = treilli_traits_t<_faceShape,_degree>;
      protected: using nodes_int_t = std::array< std::array<unsigned int, traits_t::dim >, traits_t::nnodes >;
      private: nodes_int_t m_icoo;

      public: template <typename _float_type> inline _float_type GetCoordinate(const unsigned int nodeIndex_,
									       const unsigned int dimensionIndex_) const noexcept
	{
	  static constexpr const _float_type idegree = _float_type(1.0) / _float_type(_degree);
#ifndef NDEBUG
	  Debug::IsInRange(__TRACE__,nodeIndex_,(unsigned int)0,this->s_numNodes-1);
	  Debug::IsInRange(__TRACE__,dimensionIndex_,(unsigned int)0,this->s_dimension-1);
#endif      
	  return _float_type(this->m_icoo[nodeIndex_][dimensionIndex_]) * idegree;
	};
            
      public: inline Uniform() noexcept : treilli_t<Uniform<_faceShape,_degree>,_faceShape,_degree>(this->m_icoo)
	{  	
	};
	
	inline ~Uniform() noexcept
	{
	};
    
      };





  

      template <FaceType::enum_t _faceShape,unsigned int _degree>
      class Lobatto : public treilli_t<Lobatto<_faceShape,_degree> ,_faceShape,_degree>
      {
      
      private: using traits_t = treilli_traits_t<_faceShape,_degree>;
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
#ifndef NDEBUG
	  Debug::IsInRange(__TRACE__,nodeIndex_,(unsigned int)0,this->s_numNodes-1);
	  Debug::IsInRange(__TRACE__,dimensionIndex_,(unsigned int)0,this->s_dimension-1);
#endif      
	  return this->m_rcoo[nodeIndex_][dimensionIndex_];
	};

	
      public: inline Lobatto() noexcept : treilli_t<Lobatto<_faceShape,_degree>,_faceShape,_degree>(this->m_icoo)
	{
	  Quadrature::Edge::Legendre<double,_degree-1> l;
	  std::array<double,_degree+1> w;
	  w[0]=-1.0;
	  for (int i=0;i<_degree-1;++i)
	    {
	  std::cout << "eee" << std::endl;
	  w[1+i] = l.GetPosition(i,0);
	  std::cout << "eee done" << std::endl;
	    }	  
	  w[_degree]=1.0;
	  Utils<_faceShape,_degree>::LobattoPoints(w, 
						   this->m_rcoo);

	};
	
	inline ~Lobatto() noexcept
	{
	};
      };

      
    };
  
  };

};

