#pragma once
template <typename T> struct entity_iterator : public std::iterator<std::forward_iterator_tag, T>
{
  T m_val;
public:
  entity_iterator(const T v) : m_val(v) {};
  ~entity_iterator() {};
  inline entity_iterator& operator++() noexcept { ++this->m_val; return *this; };
  inline T    	operator* () const  noexcept  { return this->m_val; };
  inline bool  operator!=(const entity_iterator& rhs) const noexcept { return m_val != rhs.m_val; };
};
