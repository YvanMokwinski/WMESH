#include <valarray>
#include <iostream>
#include <ostream>
#include <limits>
#include <array>
#include <stdarg.h>
#include "Point.hpp"
// #include "CellToNodes.hpp"
#include "Input/Medit.hpp"
#include "AHF/Mesh.hpp"



#include "Output/Medit.hpp"
#include "Output/Vtk.hpp"

#include "Program.hpp"

class MnsParrot : public Program
{
  
private: int_t m_degree;
  
public:
  
  MnsParrot(int 	 argc,
			  char * argv[]) : Program(argc, argv, true)
  {
  };
  
  virtual void Main()
  {
    
    std::cout << ReferenceShapeVolumeHexahedron::NumEntities[DimensionType::Node] << std::endl;
    if (this->HasVerbose())
      {
	LogMessage  << "degree "
		    << this->m_degree
		    << std::endl;
	
	LogMessage  << "ofilename "
		    << this->GetOfilename()
		    << std::endl;
      }
    
    if (Mpi::IsMaster())
      {
	AHF::Mesh3D mesh(this->GetInputFilename(0).c_str());	  

	{
	  Output::Medit outputMedit(this->GetOfilename());
	  outputMedit << mesh;
	}
	
	{  
	  Output::Vtk::Writer outputVtk(this->GetOfilename());
	  outputVtk << mesh;
	}
	
      }    
  };
};


int main(int 	argc_, 
	 char* argv_[])
{
  
  MnsParrot application(argc_,argv_);
  application.Run();
  return 0;
}
