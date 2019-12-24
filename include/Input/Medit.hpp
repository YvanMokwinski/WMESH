#pragma once

#include "Config.hpp"
#include "ReferenceShape.hpp"
#include "libmeshb7.h"
#include "CellToNodes.hpp"
namespace Input
{
class Medit
{
public: inline ~Medit()
  {
    GmfCloseMesh(this->inm);  
  };


public: void ReadTopology(VolumeType::enum_t volumeType_,
			  int_t * cellsToNodes_,
			  const int_t	cellsToNodesOffset_,
			  int_t * cellsToCod_,
			  const int_t	cellsToCodOffset_) const noexcept;

public: void ReadTopology(FaceType::enum_t faceType_,
			  int_t * cellsToNodes_,
			  const int_t	cellsToNodesOffset_,
			  int_t * cellsToCod_,
			  const int_t	cellsToCodOffset_) const noexcept;


  
public: template<typename cellToNodes_t> void ReadTopology(cellToNodes_t * cellsToNodes_) const noexcept;
public: template<typename _real_t> void ReadGeometry(Point2d<_real_t> * points_) const noexcept;
public: template<typename _real_t> void ReadGeometry(Point3d<_real_t> * points_) const noexcept;
public: template<typename _real_t> void ReadGeometry(std::valarray<Point3d<_real_t> >& points_) const noexcept;


public: template<typename _int_t> static int GmfType();
  
  
public: Medit(const char * filename_,...)
  {
    static enum GmfKwdCod __eFace_TO_GmfKwCod[ 2 ]	= {GmfTriangles,
							   GmfQuadrilaterals};
    
    static enum GmfKwdCod __eVolume_TO_GmfKwCod[4 ] 	= {GmfTetrahedra,
							   GmfPyramids,
							   GmfPrisms,
							   GmfHexahedra};
    
    char filename[256];
    
    { va_list args;
      va_start(args,filename_);
      vsprintf(filename,filename_,args);
      va_end(args); }
    
    this->inm = GmfOpenMesh(filename,
			    GmfRead,
			    &this->version,
			    &this->dim);
    
    if (!this->inm)
      {
	fprintf(stderr,"wrong file %s\n",filename);
	exit(1);
      }

    
    for(const auto volumeType : VolumeType::All)
      {	
	this->m_numEntities[DimensionType::Volume][volumeType] = GmfStatKwd(this->inm,__eVolume_TO_GmfKwCod[volumeType]);
      }

    for(const auto faceType : FaceType::All)
      {
	this->m_numEntities[DimensionType::Face][faceType] = GmfStatKwd(this->inm,__eFace_TO_GmfKwCod[faceType]);
      }

    this->m_numEntities[DimensionType::Edge][0] = GmfStatKwd(this->inm,GmfEdges);
    this->m_numEntities[DimensionType::Node][0] = GmfStatKwd(this->inm,GmfVertices);
    
    this->m_topologyDimension = DimensionType::Node;
    for(const auto volumeType : VolumeType::All)
      {
	if (this->m_numEntities[DimensionType::Volume][volumeType] > 0)
	  {
	    this->m_topologyDimension = DimensionType::Volume;
	    break;
	  }
      }
    
    if (this->m_topologyDimension == DimensionType::Node)
      {
	for(const auto faceType : FaceType::All)
	  {
	    if (this->m_numEntities[DimensionType::Face][faceType] > 0)
	      {
		this->m_topologyDimension = DimensionType::Face;
		break;
	      }
	  }
    }

    if (this->m_topologyDimension == DimensionType::Node)
      {
	if (this->m_numEntities[DimensionType::Edge][0] > 0)
	  {
	    this->m_topologyDimension = DimensionType::Edge;
	  }
      }

  };

public:template <DimensionType::enum_t _dimensionType>
inline int_t GetNumEntities(const typename ReferenceShapeTraits<_dimensionType>::type_t::enum_t cellType_) const noexcept
  {
    return this->m_numEntities[_dimensionType][cellType_];
  };

public:inline int_t GetNumEdges() const noexcept
  {
    return this->m_numEntities[DimensionType::Edge][0];
  };

public:inline int_t GetNumNodes() const noexcept
  {
    return this->m_numEntities[DimensionType::Node][0];
  };



  
  
public:inline int_t GetNumVolumes(VolumeType::enum_t volumeType_) const noexcept
  {
    return GetNumEntities<DimensionType::Volume>(volumeType_);
  };

public:inline int_t GetNumFaces(FaceType::enum_t faceType_) const noexcept
  {
    return GetNumEntities<DimensionType::Face>(faceType_);
  };

private: int32_t dim;
private: int32_t version;
private: int64_t inm;
private: DimensionType::enum_t m_topologyDimension;
private: int_t m_numEntities[4][4];  
};

template<> int Medit::GmfType<int>()		{return GmfInt;};
template<> int Medit::GmfType<long>()		{return GmfLong;};
template<> int Medit::GmfType<long long>()	{return GmfLong;};
template<> int Medit::GmfType<double>()		{return GmfDouble;};
template<> int Medit::GmfType<float>()		{return GmfFloat;};

template<typename _real_t> void Medit::ReadGeometry(Point2d<_real_t> * points_) const noexcept
{
  const int_t numNodes = this->GetNumNodes();
  GmfGetBlock(this->inm,
	      GmfVertices,
	      1,
	      numNodes,
	      0,
	      NULL,
	      NULL,
	      Medit::GmfType<_real_t>(), &points_[0][0], &points_[numNodes-1][0],
	      Medit::GmfType<_real_t>(), &points_[0][1], &points_[numNodes-1][1],
	      Medit::GmfType<int_t>(), &points_[0].m_cod, &points_[numNodes-1].m_cod); 
};


  template<typename _real_t> void Medit::ReadGeometry(std::valarray<Point3d<_real_t> >& points_) const noexcept
  {
    const int_t numNodes = this->GetNumNodes();
    GmfGetBlock(this->inm,
		GmfVertices,
		1,
		numNodes,
		0,
		NULL,
		NULL,
		Medit::GmfType<_real_t>(), &points_[0][0], &points_[numNodes-1][0],
		Medit::GmfType<_real_t>(), &points_[0][1], &points_[numNodes-1][1],
		Medit::GmfType<_real_t>(), &points_[0][2], &points_[numNodes-1][2],
		Medit::GmfType<int_t>(), &points_[0].m_cod, &points_[numNodes-1].m_cod); 
  };
  
template<typename _real_t> void Medit::ReadGeometry(Point3d<_real_t> * points_) const noexcept
{
  const int_t numNodes = this->GetNumNodes();
  GmfGetBlock(this->inm,
	      GmfVertices,
	      1,
	      numNodes,
	      0,
	      NULL,
	      NULL,
	      Medit::GmfType<_real_t>(), &points_[0][0], &points_[numNodes-1][0],
	      Medit::GmfType<_real_t>(), &points_[0][1], &points_[numNodes-1][1],
	      Medit::GmfType<_real_t>(), &points_[0][2], &points_[numNodes-1][2],
	      Medit::GmfType<int_t>(), &points_[0].m_cod, &points_[numNodes-1].m_cod); 
};
  
  void Medit::ReadTopology(FaceType::enum_t faceType_,
			   int_t * cellsToNodes_,
			   const int_t	cellsToNodesOffset_,
			   int_t * cellsToCod_,
			   const int_t	cellsToCodOffset_) const noexcept
  {
    const int_t numCells = this->GetNumFaces(faceType_);
    switch(faceType_)
      {
      case FaceType::Quadrilateral:
	{
	  GmfGetBlock(this->inm,
		      GmfQuadrilaterals,
		      1,
		      numCells,
		      0,
		      NULL,
		      NULL,
		      Medit::GmfType<int_t>(), &cellsToNodes_[0], &cellsToNodes_[(numCells-1)*cellsToNodesOffset_+0],
		      Medit::GmfType<int_t>(), &cellsToNodes_[1], &cellsToNodes_[(numCells-1)*cellsToNodesOffset_+1],
		      Medit::GmfType<int_t>(), &cellsToNodes_[2], &cellsToNodes_[(numCells-1)*cellsToNodesOffset_+2],
		      Medit::GmfType<int_t>(), &cellsToNodes_[3], &cellsToNodes_[(numCells-1)*cellsToNodesOffset_+3] ,
		      Medit::GmfType<int_t>(),
		      (nullptr != cellsToCod_)
		      ? &cellsToCod_[0]
		      : nullptr,
		      (nullptr != cellsToCod_)
		      ? &cellsToCod_[(numCells-1)*cellsToCodOffset_]
		      : nullptr);	  
	  break;
	}

      case FaceType::Triangle:
	{
	  GmfGetBlock(this->inm,
		      GmfTriangles,
		      1,
		      numCells,
		      0,
		      NULL,
		      NULL,
		      Medit::GmfType<int_t>(), &cellsToNodes_[0], &cellsToNodes_[(numCells-1)*cellsToNodesOffset_+0],
		      Medit::GmfType<int_t>(), &cellsToNodes_[1], &cellsToNodes_[(numCells-1)*cellsToNodesOffset_+1],
		      Medit::GmfType<int_t>(), &cellsToNodes_[2], &cellsToNodes_[(numCells-1)*cellsToNodesOffset_+2],
		      Medit::GmfType<int_t>(),
		      (nullptr != cellsToCod_)
		      ? &cellsToCod_[0]
		      : nullptr,
		      (nullptr != cellsToCod_)
		      ? &cellsToCod_[(numCells-1)*cellsToCodOffset_]
		      : nullptr);

	  break;
	}
      }
    const unsigned int numNodesInCell = FaceType::GetNumNodes(faceType_);
    for (int_t cellIndex = 0;cellIndex < numCells;++cellIndex)
      {
	for (unsigned int localNodeIndex = 0; localNodeIndex < numNodesInCell;++localNodeIndex)
	  {
	    cellsToNodes_[cellIndex * cellsToNodesOffset_ + localNodeIndex] -=1;
	  }
      }
    
  };


  void Medit::ReadTopology(VolumeType::enum_t volumeType_,
			   int_t * cellsToNodes_,
			   const int_t	cellsToNodesOffset_,
			   int_t * cellsToCod_,
			   const int_t	cellsToCodOffset_) const noexcept
  {
    const int_t numCells = this->GetNumVolumes(volumeType_);
    switch(volumeType_)
      {
      case VolumeType::Tetrahedron:
	{
	  GmfGetBlock(this->inm,
		      GmfTetrahedra,
		      1,
		      numCells,
		      0,
		      NULL,
		      NULL,
		      Medit::GmfType<int_t>(), &cellsToNodes_[0], &cellsToNodes_[(numCells-1)*cellsToNodesOffset_+0],
		      Medit::GmfType<int_t>(), &cellsToNodes_[1], &cellsToNodes_[(numCells-1)*cellsToNodesOffset_+1],
		      Medit::GmfType<int_t>(), &cellsToNodes_[2], &cellsToNodes_[(numCells-1)*cellsToNodesOffset_+2],
		      Medit::GmfType<int_t>(), &cellsToNodes_[3], &cellsToNodes_[(numCells-1)*cellsToNodesOffset_+3] ,
		      Medit::GmfType<int_t>(),
		      (nullptr != cellsToCod_)
		      ? &cellsToCod_[0]
		      : nullptr,
		      (nullptr != cellsToCod_)
		      ? &cellsToCod_[(numCells-1)*cellsToCodOffset_]
		      : nullptr);	  
	  break;
	}

      case VolumeType::Pyramid:
	{
	  GmfGetBlock(this->inm,
		      GmfPyramids,
		      1,
		      numCells,
		      0,
		      NULL,
		      NULL,
		      Medit::GmfType<int_t>(), &cellsToNodes_[0], &cellsToNodes_[(numCells-1)*cellsToNodesOffset_+0],
		      Medit::GmfType<int_t>(), &cellsToNodes_[1], &cellsToNodes_[(numCells-1)*cellsToNodesOffset_+1],
		      Medit::GmfType<int_t>(), &cellsToNodes_[2], &cellsToNodes_[(numCells-1)*cellsToNodesOffset_+2],
		      Medit::GmfType<int_t>(), &cellsToNodes_[3], &cellsToNodes_[(numCells-1)*cellsToNodesOffset_+3],
		      Medit::GmfType<int_t>(), &cellsToNodes_[4], &cellsToNodes_[(numCells-1)*cellsToNodesOffset_+4],
		      Medit::GmfType<int_t>(),
		      (nullptr != cellsToCod_)
		      ? &cellsToCod_[0]
		      : nullptr,
		      (nullptr != cellsToCod_)
		      ? &cellsToCod_[(numCells-1)*cellsToCodOffset_]
		      : nullptr);

	  break;
	}

      case VolumeType::Wedge:
	{
	  GmfGetBlock(this->inm,
		      GmfPrisms,
		      1,
		      numCells,
		      0,
		      NULL,
		      NULL,
		      Medit::GmfType<int_t>(), &cellsToNodes_[0], &cellsToNodes_[(numCells-1)*cellsToNodesOffset_+0],
		      Medit::GmfType<int_t>(), &cellsToNodes_[1], &cellsToNodes_[(numCells-1)*cellsToNodesOffset_+1],
		      Medit::GmfType<int_t>(), &cellsToNodes_[2], &cellsToNodes_[(numCells-1)*cellsToNodesOffset_+2],
		      Medit::GmfType<int_t>(), &cellsToNodes_[3], &cellsToNodes_[(numCells-1)*cellsToNodesOffset_+3],
		      Medit::GmfType<int_t>(), &cellsToNodes_[4], &cellsToNodes_[(numCells-1)*cellsToNodesOffset_+4],
		      Medit::GmfType<int_t>(), &cellsToNodes_[5], &cellsToNodes_[(numCells-1)*cellsToNodesOffset_+5],
		      Medit::GmfType<int_t>(),
		      (nullptr != cellsToCod_)
		      ? &cellsToCod_[0]
		      : nullptr,
		      (nullptr != cellsToCod_)
		      ? &cellsToCod_[(numCells-1)*cellsToCodOffset_]
		      : nullptr);
	  break;
	}
	
      case VolumeType::Hexahedron:
	{
	  GmfGetBlock(this->inm,
		      GmfHexahedra,
		      1,
		      numCells,
		      0,
		      NULL,
		      NULL,
		      Medit::GmfType<int_t>(), &cellsToNodes_[0], &cellsToNodes_[(numCells-1)*cellsToNodesOffset_+0],
		      Medit::GmfType<int_t>(), &cellsToNodes_[1], &cellsToNodes_[(numCells-1)*cellsToNodesOffset_+1],
		      Medit::GmfType<int_t>(), &cellsToNodes_[2], &cellsToNodes_[(numCells-1)*cellsToNodesOffset_+2],
		      Medit::GmfType<int_t>(), &cellsToNodes_[3], &cellsToNodes_[(numCells-1)*cellsToNodesOffset_+3],
		      Medit::GmfType<int_t>(), &cellsToNodes_[4], &cellsToNodes_[(numCells-1)*cellsToNodesOffset_+4],
		      Medit::GmfType<int_t>(), &cellsToNodes_[5], &cellsToNodes_[(numCells-1)*cellsToNodesOffset_+5],
		      Medit::GmfType<int_t>(), &cellsToNodes_[6], &cellsToNodes_[(numCells-1)*cellsToNodesOffset_+6],
		      Medit::GmfType<int_t>(), &cellsToNodes_[7], &cellsToNodes_[(numCells-1)*cellsToNodesOffset_+7],
		      Medit::GmfType<int_t>(),
		      (nullptr != cellsToCod_)
		      ? &cellsToCod_[0]
		      : nullptr,
		      (nullptr != cellsToCod_)
		      ? &cellsToCod_[(numCells-1)*cellsToCodOffset_]
		      : nullptr);
	  break;
	}
	
      }

    const unsigned int numNodesInCell = VolumeType::GetNumNodes(volumeType_);
    for (int_t cellIndex = 0;cellIndex < numCells;++cellIndex)
      {
	for (unsigned int localNodeIndex = 0; localNodeIndex < numNodesInCell;++localNodeIndex)
	  {
	    cellsToNodes_[cellIndex * cellsToNodesOffset_ + localNodeIndex] -=1;
	  }
      }
    
  };

  


  
template<> void Medit::ReadTopology<TriangleToNodes>(TriangleToNodes * cellsToNodes_) const noexcept
{
  const int_t numCells = this->GetNumFaces(FaceType::Triangle);
  GmfGetBlock(this->inm,
	      GmfTriangles,
	      1,
	      numCells,
	      0,
	      NULL,
	      NULL,
	      Medit::GmfType<int_t>(), &cellsToNodes_[0].m_nodeIds[0], &cellsToNodes_[numCells-1].m_nodeIds[0],
	      Medit::GmfType<int_t>(), &cellsToNodes_[0].m_nodeIds[1], &cellsToNodes_[numCells-1].m_nodeIds[1],
	      Medit::GmfType<int_t>(), &cellsToNodes_[0].m_nodeIds[2], &cellsToNodes_[numCells-1].m_nodeIds[2],
	      Medit::GmfType<int_t>(), &cellsToNodes_[0].m_cod, &cellsToNodes_[numCells-1].m_cod);

};

template<> void Medit::ReadTopology<QuadrilateralToNodes>(QuadrilateralToNodes * cellsToNodes_) const noexcept
{
  const int_t numCells = this->GetNumFaces(FaceType::Quadrilateral);
  GmfGetBlock(this->inm,
	      GmfQuadrilaterals,
	      1,
	      numCells,
	      0,
	      NULL,
	      NULL,
	      Medit::GmfType<int_t>(), &cellsToNodes_[0].m_nodeIds[0], &cellsToNodes_[numCells-1].m_nodeIds[0],
	      Medit::GmfType<int_t>(), &cellsToNodes_[0].m_nodeIds[1], &cellsToNodes_[numCells-1].m_nodeIds[1],
	      Medit::GmfType<int_t>(), &cellsToNodes_[0].m_nodeIds[2], &cellsToNodes_[numCells-1].m_nodeIds[2],
	      Medit::GmfType<int_t>(), &cellsToNodes_[0].m_nodeIds[3], &cellsToNodes_[numCells-1].m_nodeIds[3],
	      Medit::GmfType<int_t>(), &cellsToNodes_[0].m_cod, &cellsToNodes_[numCells-1].m_cod); 
};

template<> void Medit::ReadTopology<TetrahedronToNodes>(TetrahedronToNodes * cellsToNodes_) const noexcept
{
  const int_t numCells = this->GetNumVolumes(VolumeType::Tetrahedron);
  GmfGetBlock(this->inm,
	      GmfTetrahedra,
	      1,
	      numCells,
	      0,
	      NULL,
	      NULL,
	      Medit::GmfType<int_t>(), &cellsToNodes_[0].m_nodeIds[0], &cellsToNodes_[numCells-1].m_nodeIds[0],
	      Medit::GmfType<int_t>(), &cellsToNodes_[0].m_nodeIds[1], &cellsToNodes_[numCells-1].m_nodeIds[1],
	      Medit::GmfType<int_t>(), &cellsToNodes_[0].m_nodeIds[2], &cellsToNodes_[numCells-1].m_nodeIds[2],
	      Medit::GmfType<int_t>(), &cellsToNodes_[0].m_nodeIds[3], &cellsToNodes_[numCells-1].m_nodeIds[3],
	      Medit::GmfType<int_t>(), &cellsToNodes_[0].m_cod, &cellsToNodes_[numCells-1].m_cod); 
  
  const unsigned int numNodesInCell = VolumeType::GetNumNodes(VolumeType::Tetrahedron);
  for (int_t cellIndex = 0;cellIndex < numCells;++cellIndex)
    {
      for (unsigned int localNodeIndex = 0; localNodeIndex < numNodesInCell;++localNodeIndex)
	{
	  cellsToNodes_[cellIndex].m_nodeIds[localNodeIndex] -=1;
	}
    }
};


template<> void Medit::ReadTopology<HexahedronToNodes>(HexahedronToNodes * cellsToNodes_) const noexcept
{
  const int_t numCells = this->GetNumVolumes(VolumeType::Hexahedron);
  GmfGetBlock(this->inm,
	      GmfHexahedra,
	      1,
	      numCells,
	      0,
	      NULL,
	      NULL,
	      Medit::GmfType<int_t>(), &cellsToNodes_[0].m_nodeIds[0], &cellsToNodes_[numCells-1].m_nodeIds[0],
	      Medit::GmfType<int_t>(), &cellsToNodes_[0].m_nodeIds[1], &cellsToNodes_[numCells-1].m_nodeIds[1],
	      Medit::GmfType<int_t>(), &cellsToNodes_[0].m_nodeIds[2], &cellsToNodes_[numCells-1].m_nodeIds[2],
	      Medit::GmfType<int_t>(), &cellsToNodes_[0].m_nodeIds[3], &cellsToNodes_[numCells-1].m_nodeIds[3],
	      Medit::GmfType<int_t>(), &cellsToNodes_[0].m_nodeIds[4], &cellsToNodes_[numCells-1].m_nodeIds[4],
	      Medit::GmfType<int_t>(), &cellsToNodes_[0].m_nodeIds[5], &cellsToNodes_[numCells-1].m_nodeIds[5],
	      Medit::GmfType<int_t>(), &cellsToNodes_[0].m_nodeIds[6], &cellsToNodes_[numCells-1].m_nodeIds[6],
	      Medit::GmfType<int_t>(), &cellsToNodes_[0].m_nodeIds[7], &cellsToNodes_[numCells-1].m_nodeIds[7],
	      Medit::GmfType<int_t>(), &cellsToNodes_[0].m_cod, &cellsToNodes_[numCells-1].m_cod); 
  const unsigned int numNodesInCell = VolumeType::GetNumNodes(VolumeType::Hexahedron);
  for (int_t cellIndex = 0;cellIndex < numCells;++cellIndex)
    {
      for (unsigned int localNodeIndex = 0; localNodeIndex < numNodesInCell;++localNodeIndex)
	{
	  cellsToNodes_[cellIndex].m_nodeIds[localNodeIndex] -=1;
	}
    }
};



template<> void Medit::ReadTopology<WedgeToNodes>(WedgeToNodes * cellsToNodes_) const noexcept
{
  const int_t numCells = this->GetNumVolumes(VolumeType::Wedge);
  GmfGetBlock(this->inm,
	      GmfPrisms,
	      1,
	      numCells,
	      0,
	      NULL,
	      NULL,
	      Medit::GmfType<int_t>(), &cellsToNodes_[0].m_nodeIds[0], &cellsToNodes_[numCells-1].m_nodeIds[0],
	      Medit::GmfType<int_t>(), &cellsToNodes_[0].m_nodeIds[1], &cellsToNodes_[numCells-1].m_nodeIds[1],
	      Medit::GmfType<int_t>(), &cellsToNodes_[0].m_nodeIds[2], &cellsToNodes_[numCells-1].m_nodeIds[2],
	      Medit::GmfType<int_t>(), &cellsToNodes_[0].m_nodeIds[3], &cellsToNodes_[numCells-1].m_nodeIds[3],
	      Medit::GmfType<int_t>(), &cellsToNodes_[0].m_nodeIds[4], &cellsToNodes_[numCells-1].m_nodeIds[4],
	      Medit::GmfType<int_t>(), &cellsToNodes_[0].m_nodeIds[5], &cellsToNodes_[numCells-1].m_nodeIds[5],
	      Medit::GmfType<int_t>(), &cellsToNodes_[0].m_cod, &cellsToNodes_[numCells-1].m_cod); 
  const unsigned int numNodesInCell = VolumeType::GetNumNodes(VolumeType::Wedge);
  for (int_t cellIndex = 0;cellIndex < numCells;++cellIndex)
    {
      for (unsigned int localNodeIndex = 0; localNodeIndex < numNodesInCell;++localNodeIndex)
	{
	  cellsToNodes_[cellIndex].m_nodeIds[localNodeIndex] -=1;
	}
    }
};


template<> void Medit::ReadTopology<PyramidToNodes>(PyramidToNodes * cellsToNodes_) const noexcept
{
  const int_t numCells = this->GetNumVolumes(VolumeType::Pyramid);
  GmfGetBlock(this->inm,
	      GmfPyramids,
	      1,
	      numCells,
	      0,
	      NULL,
	      NULL,
	      Medit::GmfType<int_t>(), &cellsToNodes_[0].m_nodeIds[0], &cellsToNodes_[numCells-1].m_nodeIds[0],
	      Medit::GmfType<int_t>(), &cellsToNodes_[0].m_nodeIds[1], &cellsToNodes_[numCells-1].m_nodeIds[1],
	      Medit::GmfType<int_t>(), &cellsToNodes_[0].m_nodeIds[2], &cellsToNodes_[numCells-1].m_nodeIds[2],
	      Medit::GmfType<int_t>(), &cellsToNodes_[0].m_nodeIds[3], &cellsToNodes_[numCells-1].m_nodeIds[3],
	      Medit::GmfType<int_t>(), &cellsToNodes_[0].m_nodeIds[4], &cellsToNodes_[numCells-1].m_nodeIds[4],
	      Medit::GmfType<int_t>(), &cellsToNodes_[0].m_cod, &cellsToNodes_[numCells-1].m_cod); 
  const unsigned int numNodesInCell = VolumeType::GetNumNodes(VolumeType::Pyramid);
  for (int_t cellIndex = 0;cellIndex < numCells;++cellIndex)
    {
      for (unsigned int localNodeIndex = 0; localNodeIndex < numNodesInCell;++localNodeIndex)
	{
	  cellsToNodes_[cellIndex].m_nodeIds[localNodeIndex] -=1;
	}
    }
};
};
