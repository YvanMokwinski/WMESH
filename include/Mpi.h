#ifndef __header_Mpi_h__
#define __header_Mpi_h__


#ifdef SIMLAB_MPI
#include "MpiData.h"
#include "mpi.h"
#endif
#define MPI_CHECK(_a) { int _mpi_error = _a; if (_mpi_error != MPI_SUCCESS) {cerr << "MPI failed" << endl; exit(1);} }

class Mpi
{

public:

#ifdef SIMLAB_MPI
  typedef MPI_Request Request;
  typedef MPI_Status Status;
#else
  typedef int Request;
  typedef int Status;
#endif

  
  static int GetCpuName(int cpuId_,char cpuName_[256])
  {
    if (cpuId_==0)
      {
	return sprintf(cpuName_,"MASTER");
      }
    else
      {
	return sprintf(cpuName_,"SLAVE-%d",cpuId_);
      }
  };

  static void Init(int * argc_, char *** argv_)
  {
#ifdef SIMLAB_MPI
    MPI_CHECK(MPI_Init(argc_,argv_));
#endif
  };

  static bool IsMaster()
  {
    return (0==GetCpuId());
  };

  static int GetCpuId()
  {
#ifdef SIMLAB_MPI
    int cpu_id;
    MPI_CHECK( MPI_Comm_rank(MPI_COMM_WORLD,&cpu_id) );
    return cpu_id;
#else
    return 0;
#endif
  };

  static int Rank()
  {
#ifdef SIMLAB_MPI
    int rank;
    MPI_CHECK( MPI_Comm_rank(MPI_COMM_WORLD,&rank) );
    return rank;
#else
    return 0;
#endif
  };

  static int GetNumCpus()
  {
#ifdef SIMLAB_MPI
    int numCpus;
    MPI_CHECK( MPI_Comm_size(MPI_COMM_WORLD,&numCpus) );
    return numCpus;
#else
    return 1;
#endif
  };

  //      call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)

  static std::string GetProcessorName()
  {
#ifdef SIMLAB_MPI
    STR processor_name;    
    int len;
    MPI_CHECK( MPI_Get_processor_name(processor_name,&len) );
    return std::string(processor_name);
#else
    return std::string("undefined");
#endif
  };
  
  static float GetVersion()
  {
#ifdef SIMLAB_MPI
    int mpi_version,mpi_subversion;
    MPI_CHECK( MPI_Get_version(&mpi_version,&mpi_subversion) );
    return ((float)mpi_version) + ((float)mpi_subversion) *0.1;
#else
    return 0.0;
#endif
  };

  static void Barrier()
  {
#ifdef SIMLAB_MPI
    MPI_CHECK( MPI_Barrier(MPI_COMM_WORLD) );
#endif
  };

  static void Finalize()
  {
#ifdef SIMLAB_MPI
    MPI_CHECK( MPI_Finalize() );
#endif
  };


  static void Wait(Request * request_,Status * status_)
  {
#ifdef SIMLAB_MPI
    MPI_Wait(request_,status_);
#endif
  };


  template <typename _type> static void Send(_type * data_,int count_,int destination_,int tag_)
    {
#ifdef SIMLAB_MPI    
      MPI_Send(data_,count_,MpiData<_type>::Datatype,destination_,tag_,MPI_COMM_WORLD);
#endif
    };
    
  template <typename _type> static void Recv(_type * data_,int count_,int source_,int tag_)
    {
#ifdef SIMLAB_MPI    
      MPI_Recv(data_,count_,MpiData<_type>::Datatype,source_,tag_,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
#endif
    };


  template <typename _type> static Request ISend(_type * data_,int count_,int destination_,int tag_)
    {
      Request send_request;      
#ifdef SIMLAB_MPI    
      const int rc = MPI_Isend(data_,count_,MpiData<_type>::Datatype,destination_,tag_,MPI_COMM_WORLD,&send_request);
      if (rc != MPI_SUCCESS) 
	{
	  printf ("MPI_Isend error %d\n", rc);
	  exit(1);
	}      
#endif
      return send_request;
    };
    

    template <typename _type> static Request IRecv(_type * data_,int count_,int source_,int tag_)
    {
      Request recv_request;      
#ifdef SIMLAB_MPI    
      const int rc = MPI_Irecv(data_, count_, MpiData<_type>::Datatype, source_, tag_, MPI_COMM_WORLD, &recv_request);
      if (rc != MPI_SUCCESS) 
	{
	  printf ("MPI_Irecv error %d\n", rc);
	  exit(1);
	}
#endif
      return recv_request;
    };



  
#ifdef SIMLAB_MPI
  template <MPI_Op _op,typename _floatType> static _floatType Allreduce(_floatType local_value)
  {
    _floatType global_value;
    MPI_Allreduce(&local_value, &global_value, 1, MpiData<_floatType>::Datatype, _op, MPI_COMM_WORLD);
    return global_value;
  };
  template <MPI_Op _op,typename _floatType> static _floatType Reduce(_floatType local_value)
  {
    _floatType global_value;
    MPI_Reduce(&local_value, &global_value, 1, MpiData<_floatType>::Datatype, _op, 0,MPI_COMM_WORLD);
    return global_value;
  };
#endif


  template <typename _floatType> static _floatType Allsum(_floatType local_value)
    {
#ifdef SIMLAB_MPI
      return Allreduce<MPI_SUM,_floatType>(local_value);
#else
      return local_value;
#endif
    };

  template <typename _floatType> static _floatType Allmax(_floatType local_value)
    {
#ifdef SIMLAB_MPI
      return Allreduce<MPI_MAX,_floatType>(local_value);
#else
      return local_value;
#endif
    };

  template <typename _floatType> static _floatType Allmin(_floatType local_value)
    {
#ifdef SIMLAB_MPI
      return Allreduce<MPI_MIN,_floatType>(local_value);
#else
      return local_value;
#endif
    };

  template <typename _floatType> static _floatType Allprod(_floatType local_value)
    {
#ifdef SIMLAB_MPI
      return Allreduce<MPI_PROD,_floatType>(local_value);
#else
      return local_value;
#endif
    };

  template <typename _floatType> static _floatType Sum(_floatType local_value)
    {
#ifdef SIMLAB_MPI
      return Reduce<MPI_SUM,_floatType>(local_value);
#else
      return local_value;
#endif
    };

  template <typename _floatType> static _floatType Max(_floatType local_value)
    {
#ifdef SIMLAB_MPI
      return Reduce<MPI_MAX,_floatType>(local_value);
#else
      return local_value;
#endif
    };

  template <typename _floatType> static _floatType Min(_floatType local_value)
    {
#ifdef SIMLAB_MPI
      return Reduce<MPI_MIN,_floatType>(local_value);
#else
      return local_value;
#endif
    };

  template <typename _floatType> static _floatType Prod(_floatType local_value)
    {
#ifdef SIMLAB_MPI
      return Reduce<MPI_PROD,_floatType>(local_value);
#else
      return local_value;
#endif
    };

};

#include "MpiRule.h"


#endif
