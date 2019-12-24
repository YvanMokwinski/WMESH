#pragma once

template <FaceType::enum_t _faceType> struct FaceCompare;

template <> struct FaceCompare<FaceType::Triangle>
{
public: template<typename int_t>
static inline bool AreSame(const int_t thisFaceToNodes_[],
			   const int_t thatFaceToNodes_[]) noexcept
  {
    
    if  ( ( (thisFaceToNodes_[0] == thatFaceToNodes_[2]) && (thisFaceToNodes_[1] == thatFaceToNodes_[1]) && (thisFaceToNodes_[2] == thatFaceToNodes_[0]) )
	  ||
	  ( (thisFaceToNodes_[0] == thatFaceToNodes_[1]) && (thisFaceToNodes_[1] == thatFaceToNodes_[0]) && (thisFaceToNodes_[2] == thatFaceToNodes_[2]) )
	  ||
	  ( (thisFaceToNodes_[0] == thatFaceToNodes_[0]) && (thisFaceToNodes_[1] == thatFaceToNodes_[2]) && (thisFaceToNodes_[2] == thatFaceToNodes_[1]) ) )
      {
	return true;
      }
#if 1
    if  ( ( (thisFaceToNodes_[0] == thatFaceToNodes_[0]) && (thisFaceToNodes_[1] == thatFaceToNodes_[1]) && (thisFaceToNodes_[2] == thatFaceToNodes_[2]) )
	  ||
	  ( (thisFaceToNodes_[0] == thatFaceToNodes_[1]) && (thisFaceToNodes_[1] == thatFaceToNodes_[2]) && (thisFaceToNodes_[2] == thatFaceToNodes_[0]) )
	  ||
	  ( (thisFaceToNodes_[0] == thatFaceToNodes_[2]) && (thisFaceToNodes_[1] == thatFaceToNodes_[0]) && (thisFaceToNodes_[2] == thatFaceToNodes_[1]) ) )
      {
	fprintf(stderr,"bad orientation\n");
	exit(1);
      }
#endif    
    return false;    
  };
  
public: template<typename int_t>
static inline int_t HashCode(const int_t numCells_,
			     const int_t cncface_[]) noexcept
  {
    const int_t h = ( (cncface_[0] < cncface_[1])
		      ? ( (cncface_[0] < cncface_[2]) ? cncface_[0] : cncface_[2] ) 
		      : ( (cncface_[1] < cncface_[2]) ? cncface_[1] : cncface_[2] ) ) % numCells_;
    
    return (h<0)
      ? -h
      : h;
  };

  
public: template <unsigned int _degree,typename _array_integer_t>
static inline void  OrientationPermutation(const char 	orientation_,
					   _array_integer_t&	perm_) noexcept
  {
    static constexpr unsigned int n = _degree-1;
    unsigned int k 		= 0;
    for (unsigned int i=1;i<=n;++i)
      {
	for (unsigned int j=1;j<=n-i;++j)
	  {
	    unsigned int
	      r = 0,
	      s = 0;
	    switch(orientation_)
	      {
	      case -1:
		{
		  r = i;
		  s = j;
		  break;
		}
	      case 1:
		{
		  r = j;
		  s = i;
		  break;
		}
	      case -2:
		{
		  r = j;
		  s = _degree- i - j;
		  break;
		}
	      case 2:
		{
		  r = _degree- i - j;
		  s = j;
		  break;
		}
	      case -3:
		{
		  r = _degree- i - j;
		  s = i;
		  break;
		}
	      case 3:
		{
		  r = i;
		  s = _degree- i - j;
		  break;
		}
	      }
	    
	    perm_[k++] = ((n-1)*n)/2 - ((n-r)*(n-r+1))/2 + (s-1);
	    
	    //   4
	    //   3 8 
	    //   2 7 11
	    //   1 6 10 13
	    //   0 5 9  12 14
#if 0
	    if (perm_[k-1]!=0)
	      {
		printf("################# i %lld \n",i);
		printf("################# j %lld \n",j);
		printf("################# r %lld \n",r);
		printf("################# s %lld \n",s);
	      }
#endif	    
	  }
      }
  };

public: template <typename _array_integer_t> static inline char Orientation(const _array_integer_t& icnc_,
									    const _array_integer_t& jcnc_) noexcept
  {
    const auto i0 = icnc_[0];
    const auto i1 = icnc_[1];
    const auto i2 = icnc_[2];
    
    const auto j0 = jcnc_[0];
    const auto j1 = jcnc_[1];
    const auto j2 = jcnc_[2];
    
    if ( (i0==j0) && (i1==j1) && (i2==j2) )
      {
	return 1;
      }    
    else if ( (i0==j1) && (i1==j2) && (i2==j0) )
      {
	return 2;
      }
    else if ( (i0==j2) && (i1==j0) && (i2==j1) )
      {
	return 3;
      }
    else if ( (i0==j0) && (i1==j2) && (i2==j1) )
      {
	return -1;
      }    
    else if ( (i0==j1) && (i1==j0) && (i2==j2) )
      {
	return -2;
      }
    else if ( (i0==j2) && (i1==j1) && (i2==j0) )
      {
	return -3;
      }
    
    return 0;
  };
  
};




template <> struct FaceCompare<FaceType::Quadrilateral>
{
  
public: template<typename int_t> static inline bool AreSame(const int_t cncface[],
								 const int_t cncface2[]) noexcept
  {    
    unsigned int imin = 0;
    int_t c = cncface[0];
    for (unsigned int i=1;i<4;++i)
      {
	if (c > cncface[i])
	  {
	    imin = i;
	    c = cncface[i];
	  }
      }
    
    unsigned int imin2 = 0;
    int_t c2 = cncface2[0];
    for (unsigned int i=1;i<4;++i)
      {
	if (c2 > cncface2[i])
	  {
	    imin2 = i;
	    c2 = cncface2[i];
	  }
      }

    return ( (c == c2) && (cncface[(imin+2)%4] == cncface2[(imin2+2)%4]) );    
  };
  
public: template<typename int_t>
static inline int_t HashCode(const int_t numCells_,
			     const int_t cncface_[]) noexcept
  {
    int_t h = cncface_[0];
    h = (h < cncface_[1]) ? h : cncface_[1];
    h = (h < cncface_[2]) ? h : cncface_[2];
    h = (h < cncface_[3]) ? h : cncface_[3];    
    h = h % numCells_;
    if (h<0)
      h *= -1;
    return h;
  };
  
  
public: template <unsigned int _degree,typename _array_integer_t>
static inline void  OrientationPermutation(const char 	orientation_,
					   _array_integer_t&	perm_) noexcept
  {    
    static constexpr const unsigned int n = _degree-1;
    unsigned int k = 0;
    for (unsigned int i=1;i<=n;++i)
      {
	for (unsigned int j=1;j<=n;++j)
	  {
	    unsigned int r = 0;
	    unsigned int s = 0;
	    switch(orientation_)
	      {
		/* on a interchange 2 et -2 attends toi a faire de meme avec 3 et 4*/
	      case 1:
		{
		  r = j;
		  s = i;
		  /* (j,i) 
		     0,3,2,1
		  */
		  break;
		}
	      case -2:
		{
		  /* (j,k-i)
		     1,2,3,0
		  */
		  r = j;
		  s = _degree-i;		      
		  break;
		}
	      case -3:
		{
		  /* (k-i,k-j) 
		     2,3,0,1
		  */
		  r = _degree-i;
		  s = _degree-j;
		  break;
		}
	      case -4:
		{
		  r = _degree-j;
		  s = i;
		  /* (k-j,i) 
		     3,0,1,2*/
		  break;
		}
		
	      case -1:
		{
		  /* (i,j)
		     0,1,2,3
		  */
		  r = i;
		  s = j;
		  break;
		}
		
	      case 2:
		{
		  r = _degree-i;
		  s = j;
		  /* (k-i,j)
		     1,0,3,2
		  */
		  break;
		}	      
	      case 3:
		{
		  r = _degree-j;
		  s = _degree-i;
		  /* (k-j,k-i)
		     2,1,0,3
		  */
		  break;
		}
	      case 4:
		{
		  r = i;
		  s = _degree-j;
		  /* (i,k-j)
		     3,2,1,0
		  */
		  break;
		}
		
	      }
	    
	    perm_[k++] = n*(r-1)+(s-1);
#if 0
	    if (perm_[k-1]!=0)
	      {
		printf("################# i %lld \n",i);
		printf("################# j %lld \n",j);
		printf("################# r %lld \n",r);
		printf("################# s %lld \n",s);
	      }
#endif
	  } 	
      }     
  };
  
public:  template <typename _array_integer_t> static inline char Orientation(const _array_integer_t& icnc_,
									     const _array_integer_t& jcnc_) noexcept
  {
#if 0
    static const int_t facequad_orientation[]={0,1,2,3,/*1*/
					   1,2,3,0,/*2*/
					   2,3,0,1,/*3*/
					   3,0,1,2,/*4*/
					   
					   0,3,2,1,/*-1*/
					   1,0,3,2,/*-2*/
					   2,1,0,3,/*-3*/
					   3,2,1,0/*-4*/};
#endif
    const auto i0 = icnc_[0];
    const auto i1 = icnc_[1];
    const auto i2 = icnc_[2];
    const auto i3 = icnc_[3];
    
    const auto j0 = jcnc_[0];
    const auto j1 = jcnc_[1];
    const auto j2 = jcnc_[2];
    const auto j3 = jcnc_[3];
    
    if ( (i0==j0) && (i1==j1) && (i2==j2) && (i3==j3) )
      {
	return 1;
      }    
    else if ( (i0==j1) && (i1==j2) && (i2==j3) && (i3==j0) )
      {
	return 2;
      }
    else if ( (i0==j2) && (i1==j3) && (i2==j0) && (i3==j1) )
      {
	return 3;
      }
    else if ( (i0==j3) && (i1==j0) && (i2==j1) && (i3==j2) )
      {
	return 4;
      }  
    else if ( (i0==j0) && (i1==j3) && (i2==j2) && (i3==j1) )
      {
	return -1;
      }
    else if ( (i0==j1) && (i1==j0) && (i2==j3) && (i3==j2) )
      {
	return -2;
      }
    else if ( (i0==j2) && (i1==j1) && (i2==j0) && (i3==j3) )
      {
	return -3;
      }
    else if ( (i0==j3) && (i1==j2) && (i2==j1) && (i3==j0) )
      {
	return -4;
      }
    
    return 0;
  };

};

using TriangleCompare = FaceCompare<FaceType::Triangle>;
