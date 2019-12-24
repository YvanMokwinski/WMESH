#ifndef __HEADER_CRTP_Treilli_h__
#define __HEADER_CRTP_Treilli_h__

namespace FiniteElement
{

  template <typename _derivedClass> struct Traits_CRTP_Treilli {    
    // eTypeShape;
  };


  template <typename _derivedClass> class CRTP_Treilli {    

    // Prohibit Copy
    inline CRTP_Treilli(const CRTP_Treilli<_derivedClass>&);   
    // Curiously recurring template pattern
    inline const _derivedClass & asImp() const { return static_cast<const _derivedClass&>(*this); }    
    
  protected:
    inline CRTP_Treilli(){};

    typedef Traits_CRTP_Treilli<_derivedClass> Traits;
    typedef typename Traits::eTypeShape eTypeShape;
  public:
    
    inline unsigned int 	GetNumNodesInCell() 	const
    {
      CHECK_INTERFACE_IMPLEMENTATION(asImp().GetNumNodesInCell());
      return asImp().GetNumNodesInCell();      
    };

    inline eTypeShape GetShape() const
    {
      CHECK_INTERFACE_IMPLEMENTATION(asImp().GetShape());
      return asImp().GetShape();      
    };
    
    /**
       \brief Get the degree related to a TreilliFace
       @return degree
    */
    inline unsigned int GetDegree() const
    {
      CHECK_INTERFACE_IMPLEMENTATION(asImp().GetDegree());
      return asImp().GetDegree();
    };
    
    /**
       \brief Get the number of nodes
       @return number of nodes
    */
    inline unsigned int GetNumNodes() const
    {
      CHECK_INTERFACE_IMPLEMENTATION(asImp().GetNumNodes());
      return asImp().GetNumNodes();
    };
    
    /**
       \brief Get the number of linear element of a TreilliFace
       @param self_ pointer to the TreilliFace object
       @return number of linear element
    */
    inline unsigned int GetNumSubElements() const
    {
      CHECK_INTERFACE_IMPLEMENTATION(asImp().GetNumSubElements());
      return asImp().GetNumSubElements();
    };
    
    /**
       \brief Get the dimension of the coordinates of the vertices
       @param self_ pointer to the TreilliFace object
       @return dimension
    */
    inline unsigned int GetDimension() const
    {
      CHECK_INTERFACE_IMPLEMENTATION(asImp().GetDimension());
      return asImp().GetDimension();
    };

    /**
       \brief Get the (double) coordinates of the TreilliFace
       @param self_ pointer to the TreilliFace object
       @param coo_ double coordinates
       @param icooff_ pointer to the offset value
    */
    template <typename _float_type> inline _float_type GetCoordinate(const unsigned int&nodeIndex_,
								     const unsigned int&dimensionIndex_) const
    {
      CHECK_INTERFACE_IMPLEMENTATION(asImp().GetCoordinate<_float_type>(nodeIndex_,dimensionIndex_));
      return asImp().GetCoordinate<_float_type>(nodeIndex_,dimensionIndex_);
    };
    
    /**
       \brief Get the connectivity of the TreilliFace
       @param self_ pointer to the TreilliFace object
       @param cnc_ connectivity
       @param icncoff_ pointer to the offset value
    */
    inline unsigned int GetNodeIndex(const unsigned int&subElementIndex_,const unsigned int&localNodeIndex_) const
    {
      CHECK_INTERFACE_IMPLEMENTATION(asImp().GetNodeIndex(subElementIndex_,localNodeIndex_));
      return asImp().GetNodeIndex(subElementIndex_,localNodeIndex_);
    };
    
  };

};

#endif
