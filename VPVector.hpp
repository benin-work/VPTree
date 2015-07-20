/* 
 * File:   VPVector.hpp
 * Author: V.Karlov
 *
 * Created on 19.04.2010
 */

#ifndef _VPVECTOR_H
#define	_VPVECTOR_H

#include "types.h"

namespace vptree {

template<typename DataType, typename DistanceType>
class VPVector
{
public:
	typedef DataType ValueType;
	typedef std::vector<ValueType> DataArray;
	
public:
	VPVector(const size_t size);
	VPVector(const size_t size, const DataType* data, const uint64_t new_id);
	VPVector(const VPVector& orig);
	VPVector& operator=(const VPVector& orig);
	
	virtual ~VPVector();
	
	const DataArray& get_data() const;

	const size_t get_size() const;

	void set_id(const uint64_t new_id);
	const uint64_t get_id() const;

	DataType& operator[](size_t const n);
	const DataType& operator[](size_t const n) const;
	
	const DistanceType GetDistance(VPVector* vector) const;
	const DistanceType GetDistance(const VPVector& vector) const;

	const DistanceType GetSquaredDistance(VPVector* vector) const;
	const DistanceType GetSquaredDistance(const VPVector& vector) const;

	static const DistanceType GetSquaredDistance(const VPVector& vector1, const VPVector& vector2);

private:
	size_t m_size;
	DataArray m_data;
	uint64_t m_id;
};

template< typename DataType, typename DistanceType>
inline bool operator==(VPVector<DataType, DistanceType> const& A, VPVector<DataType, DistanceType> const& B)
{
	assert(A.get_size() == B.get_size());
	
	for(size_t c_index = 0; c_index < A.get_size(); ++c_index)
		if(A[c_index] != B[c_index])
			return false;

	return true;
}

// VPVector implementation
template< typename DataType, typename DistanceType>
VPVector<DataType, DistanceType>::VPVector(const size_t size)
: m_size(size)
, m_id(0)
{
	m_data.reserve(size);
}

template< typename DataType, typename DistanceType>
VPVector<DataType, DistanceType>::VPVector(const size_t size, const DataType* data, const uint64_t new_id)
: m_id(0)
{
	m_size = size;
	m_id = new_id;
	m_data.reserve(size);

	for(size_t c_index = 0; c_index < m_size; ++c_index)
		m_data.push_back(data[c_index]);
}

template< typename DataType, typename DistanceType>
VPVector<DataType, DistanceType>::VPVector(const VPVector& orig)
{
	if(&orig == this)
		return;

	m_size = orig.m_size;
	m_data = orig.m_data;
	m_id = orig.m_id;
}

template< typename DataType, typename DistanceType>
VPVector<DataType, DistanceType>& VPVector< DataType, DistanceType>::operator=(const VPVector& orig)
{
	if(&orig != this)
	{
		m_size = orig.m_size;
		m_data = orig.m_data;
		m_id = orig.m_id;
	}

	return *this;
}

template< typename DataType, typename DistanceType>
VPVector< DataType, DistanceType>::~VPVector()
{
	m_data.clear();
}

template< typename DataType, typename DistanceType>
const typename VPVector<DataType, DistanceType>::DataArray& VPVector<DataType, DistanceType>::get_data() const
{
	return m_data;
}

template< typename DataType, typename DistanceType>
const size_t VPVector<DataType, DistanceType>::get_size() const
{
	return m_size;
}

template< typename DataType, typename DistanceType>
const uint64_t VPVector<DataType, DistanceType>::get_id() const
{
	return m_id;
}

template< typename DataType, typename DistanceType>
void VPVector<DataType, DistanceType>::set_id(const uint64_t new_id)
{
	m_id = new_id;
}

template< typename DataType, typename DistanceType>
DataType& VPVector<DataType, DistanceType>::operator[](size_t const n)
{
	assert(n < m_size);

	return m_data[n];
}

template< typename DataType, typename DistanceType>
const DataType& VPVector<DataType, DistanceType>::operator[](size_t const n) const
{
	assert(n < m_size);

	return m_data[n];
}

template< typename DataType, typename DistanceType>
const DistanceType VPVector<DataType, DistanceType>::GetDistance(const VPVector& vector) const
{
	assert(get_size() == vector.get_size());
	
	// Expensive method for test
	DistanceType c_result = 0;
	DistanceType c_diff = 0;
	for(size_t c_index = 0; c_index < vector.get_size(); ++c_index)
	{
		c_diff = (m_data[c_index] - vector.m_data[c_index]);
		c_result +=  c_diff * c_diff;
	}

	c_result = sqrt(c_result);

	return c_result;
}

template< typename DataType, typename DistanceType>
const DistanceType VPVector<DataType, DistanceType>::GetDistance(VPVector* vector) const
{
	return this->GetDistance(*vector);
}

template< typename DataType, typename DistanceType>
const DistanceType VPVector<DataType, DistanceType>::GetSquaredDistance(const VPVector& vector) const
{
	assert(get_size() == vector.get_size());

	// Optimize !!
	// Max data crop
	DistanceType c_result = 0;
	DistanceType c_diff = 0;
	for(size_t c_index = 0; c_index < vector.get_size(); ++c_index)
	{
		c_diff = (m_data[c_index] - vector.m_data[c_index]);
		c_result +=  c_diff * c_diff;
	}

	c_result /= m_size;

	return c_result;
}

template< typename DataType, typename DistanceType>
const DistanceType VPVector<DataType, DistanceType>::GetSquaredDistance(VPVector* vector) const
{
	//const VPVector& t_vector(*vector);
	return this->GetDistance(*vector);
}

template< typename DataType, typename DistanceType>
const DistanceType VPVector<DataType, DistanceType>::GetSquaredDistance(const VPVector& vector1, const VPVector& vector2)
{
	return vector1.GetDistance(vector2);
}


} // namespace vptree

#endif	/* _VPVECTOR_H */

