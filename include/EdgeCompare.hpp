#pragma once

struct EdgeCompare
{
  
public: template<typename int_t>
static inline bool HasPositiveOrientation(const int_t thisEdgeToNodes_[]) noexcept
  {
    return (thisEdgeToNodes_[0] < thisEdgeToNodes_[1]);
  };
  
public: template<typename int_t>
static inline bool AreSame(const int_t thisEdgeToNodes_[],
			   const int_t thatEdgeToNodes_[]) noexcept
  {
    return ( ( (thisEdgeToNodes_[0] == thatEdgeToNodes_[0]) && (thisEdgeToNodes_[1] == thatEdgeToNodes_[1]) )
	     ||
	     ( (thisEdgeToNodes_[1] == thatEdgeToNodes_[0]) && (thisEdgeToNodes_[0] == thatEdgeToNodes_[1]) ) );
    
  };
  
public: template<typename int_t>
static inline int_t HashCode(const int_t numCells_,
			     const int_t thisEdgeToNodes_[]) noexcept
  {
    const int_t h = ( (thisEdgeToNodes_[0] < thisEdgeToNodes_[1]) ? thisEdgeToNodes_[0] : thisEdgeToNodes_[1] ) % numCells_;
    return (h<0)
      ? -h
      : h;
  };
  
};
