#ifndef __header_MULTITON_H__
#define __header_MULTITON_H__

#include <map>

template <typename T,typename Key> class Multiton
{
 public:
  static T* 	GetInstance	(const Key& key_)
  {
    typename std::map<Key, T*>::iterator it = instances.find(key_);
    if (it != instances.end()) 
      {
	return (T*)(it->second);
      }
    else
      {
	T* instance = new T(key_);
	instances[key_] = instance;
	return instance;
      }
  };
  static void 	KillInstances	()
  {
    for (typename std::map<Key, T*>::iterator it = instances.begin(); 
	 it != instances.end(); 
	 ++it) 
      {
	delete (*it).second;
      }
  };
 protected:
  static std::map<Key,T*> instances;
 private:
};

template <typename T, typename Key> std::map<Key, T*> Multiton<T,Key>::instances;

#endif
