#pragma once

template <unsigned int _dim> class Vertices : public std::valarray< std::array<double,_dim> >
{
public: using Vertex = std::array<double,_dim>;
public: inline Vertices(const long unsigned int numVertices_) : std::valarray< std::array<double,_dim> >(numVertices_)
  {    
  };
  
  
public: inline unsigned int GetDim() const noexcept
  {
    return _dim;
  };
  
public:  inline void Get(const I 		ith_,
			 pR 			xyz_) const noexcept
  {
    const auto& v = (*this)[ith_];
    xyz_[0] = v[0];
    xyz_[1] = v[1];
    xyz_[2] = v[2];
  };

};
