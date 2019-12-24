#pragma once
#if 0
template <unsigned int _dimension,typename _real_t> class MeshGeometry;

template <typename _real_t> 
struct Traits_CRTP_MeshGeometry<MeshGeometry<3,_real_t> >
{
public: using pts_t = Point3d<double>;
  // public: using meshtopology_t = MeshTopology3D;  
};
#endif
template <unsigned int _dimension,typename _real_t> class MeshGeometry
{
private: int_t m_numPoints = 0;
private: Point3d<_real_t>* m_points = nullptr;
  
public: static constexpr const unsigned int Dimension = _dimension;
  
public: MeshGeometry(const int_t numPoints_)
  : m_numPoints(numPoints_),
    m_points(new Point3d<_real_t>[m_numPoints])  
  {
  };
  
public: ~MeshGeometry()
  {
    if (nullptr != this->m_points)
      {
	delete [] this->m_points;
	this->m_points = nullptr;
      }
  };

public:  Point3d<_real_t>* GetPoints() noexcept
  {
    return m_points;
  };
  
public: template <typename _array_t> void GetPoint(const int_t pointId_,
						   _array_t& values_) const noexcept
  {
    for (unsigned int i=0;i<_dimension;++i)
      {
	values_[i] = this->m_points[i].m_x[i];
      }
  };
  
};
