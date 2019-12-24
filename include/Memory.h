#ifndef __header_Memory_h__
#define __header_Memory_h__


#include "Monitor.h"
#include "Config.h"
#include "Type.h"

#ifdef __cplusplus
extern "C"
{
#endif

  /**
     \brief Define memory type values
   */
#if ! MNS_DEFINE_ENUM_AS_INTEGER
  typedef enum __eMemoryType { MemoryType_STD=0, /**< Use standard memory */
			       MemoryType_SHARED, /**< Use shared memory */
			       MemoryType_SHAREDMAP /**< Use shared map memory */ } MemoryType;  
#else
/**\brief Use standard memory */
#define MemoryType_STD 		0 
/**\brief Use shared memory */
#define MemoryType_SHARED 	1
 /**\brief Use shared map memory */ 
#define MemoryType_SHAREDMAP 	2
  typedef I MemoryType;
#endif

  /**
     \brief Empty structure to define Memory
  */
  typedef struct
  {
  } Memory,*RESTRICT pMemory;

  /**
     \brief Define a type for constant Memory pointer
  */
  typedef const Memory * RESTRICT 	cst_pMemory;

  
  /**
     \brief Get the size of the Memory
     @param self_ pointer to the Memory object
     @return size of the allocated memory
  */
  size_t 			Memory_get_size			(cst_pMemory 		self_);

  /**
     \brief Get the key related to the Memory
     @param self_ pointer to the Memory object
     @return identification key of memory
  */
  int 				Memory_get_key			(cst_pMemory 		self_);

  /**
     \brief Get the MemoryType value related to the Memory
     @param self_ pointer to the Memory object
     @return MemoryType
  */
  MemoryType			Memory_get_type			(cst_pMemory 		self_);

  /**
     \brief Access to the memory with write permission
     @param self_ pointer to the Memory object
     @return address
  */
  void * RESTRICT		Memory_at			(pMemory		self_);

  /**
     \brief Access to the memory with read only permission
     @param self_ pointer to the Memory object
     @return address
  */
  const void * RESTRICT		cst_Memory_at			(cst_pMemory		self_);

  /**
     \brief Kill the Memory object
     @param self_ pointer to the Memory object
     @return NULL
  */
  pMemory			Memory_kill			(pMemory 		self_);

  /**
     \brief Create a Memory object
     @param ithread_ id of the thread which creates the Memory object
     @param key_ id of the memory
     @param size_ size of the memory
     @param type_ MemoryType of the Memory 
     @return pointer to the Memory object
  */
  pMemory			Memory_new			(const int 		ithread_,
								 const int 		key_,
								 const size_t 		size_,
								 const MemoryType 	type_);  

  /**
     \brief Connect to an existing Memory object
     @param ithread_ id of the thread which connects the Memory object
     @param key_ id of the memory
     @param size_ size of the memory
     @param type_ MemoryType of the Memory 
     @return pointer to the Memory object
  */
  pMemory			Memory_connect			(const int 		ithread_,
								 const int 		key_,
								 const size_t 		size_,
								 const MemoryType 	type_);

  /**
     \brief Disconnect to an existing Memory object
     @param self_ pointer to the Memory object
     @return NULL
  */
  pMemory			Memory_disconnect		(pMemory 		self_);

  /**
     \brief Internal structure of a Memory object
  */
  typedef struct
  {
    size_t 	size; /**< size of the memory */
    int		memkey; /**< key id of the memory */
    int		ithread; /**< id of the associated thread */
    MemoryType	memtype; /**< kind of the memory */
  } MemoryContext,* RESTRICT  pMemoryContext;
  

  /**
     \brief pointer to a constant memory context
  */
  typedef const MemoryContext * RESTRICT 	cst_pMemoryContext;


  /**
     \brief Define a MemoryContext object
     @param self_ pointer to a MemoryContext object
     @param ithread_ id of the thread which defines the Memory object
     @param key_ id of the memory
     @param type_ MemoryType of the Memory 
  */
  void 				MemoryContext_def		(pMemoryContext 	self_,
								 const int 		ithread_,
								 const int 		key_,
								 const MemoryType 	type_);  


  /**
     \brief Create a Memory object from a MemoryContext object
     @param self_ pointer to a MemoryContext object
     @param size_ size of the memory to allocate
     @return pointer to a new Memory object
  */
  pMemory			MemoryContext_new		(cst_pMemoryContext 	self_,
								 const size_t 		size_);



  /**
     \brief Connect to a Memory object from a MemoryContext object
     @param self_ pointer to a MemoryContext object
     @param size_ size of memory
     @return pointer to an existing Memory object
  */
  pMemory			MemoryContext_connect		(cst_pMemoryContext 	self_,const size_t size_);  



#ifdef __cplusplus
};
#endif


#endif
