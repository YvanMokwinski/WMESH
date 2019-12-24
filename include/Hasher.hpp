#pragma once


template <typename _int_t> class Hasher
{
private: _int_t m_size;
private: _int_t*__restrict__ m_link;
public: static constexpr const _int_t s_default_hash = - std::numeric_limits<_int_t>::max();
  
public: inline Hasher(const _int_t size_) noexcept
  : m_size(size_)//,
    //    m_link(new _int_t[size_])
  {
    this->m_link = new _int_t[size_];
    this->Reset();
  };
  
public: inline ~Hasher() noexcept
  {
    if (this->m_link)
      {
	delete[] this->m_link;
	this->m_link = nullptr;
      }
  };
  
public: inline const _int_t&operator[](const _int_t hash_value_) const noexcept
  {
    return this->m_link[hash_value_];
  };
  
public: inline _int_t&operator[](const _int_t hash_value_) noexcept
  {
    return this->m_link[hash_value_];
  };  
  
public: inline void Reset() noexcept
  {
    for (_int_t i=0;i<this->m_size;++i)
      {
	this->m_link[i] = s_default_hash;
      }
  };
  
};
