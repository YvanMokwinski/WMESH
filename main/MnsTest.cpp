#include <valarray>
#include <iostream>
#include <ostream>
#include <limits>
#include <array>

#include "Program.hpp"

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

//#include "Config.hpp"

#include "FiniteElement/Lagrange.hpp"
#include "FiniteElement/LagrangeNode.hpp"
#include "FiniteElement/LagrangeEdge.hpp"
#include "FiniteElement/LagrangeFaceTriangle.hpp"
#include "FiniteElement/LagrangeFaceQuadrilateral.hpp"
#include "FiniteElement/LagrangeVolumeHexahedron.hpp"
#include "FiniteElement/LagrangeVolumeTetrahedron.hpp"
#include "FiniteElement/LagrangeVolumeWedge.hpp"

#include "EdgeCompare.hpp"
#include "FaceCompare.hpp"

#include "generate_finite_element_space.hpp"

int main(int 	argc_, 
	 char* argv_[])
{
  int_t numCells = 2;
  int_t cnc[] = {0,1,2,3,
		 4,5,6,7};
  int_t cncLd = 4;
  int_t numDofs = 0;
  int_t cncP[1024];
  int_t numDofsPerCell = 20;
  GenerateFiniteElementSpace<3,VolumeType::Tetrahedron>(numCells,
							cnc,
							cncLd,
							&numDofs,
							cncP,
							numDofsPerCell);
  for (int i=0;i<numCells;++i)
    {
      for (int j=0;j<numDofsPerCell;++j)
	{
	  std::cout << " " << cncP[i*numDofsPerCell+j];
	}
      std::cout << std::endl;
    }

  return 0;
}






