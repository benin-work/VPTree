/* 
 * File:   VPNode.hpp
 * Author: V.Karlov
 *
 * Created on 19.04.2010
 */

#ifndef _VPNODE_H
#define	_VPNODE_H


namespace vptree {

// ValType == VPVector
// m_value == vantage point
template <typename ValType, typename DistanceType>
class VPNode
{
public:
	typedef boost::shared_ptr<VPNode<ValType, DistanceType> > VPNodePtr;
	typedef std::vector<VPNodePtr> ChildList;
	typedef std::vector<DistanceType> DistanceList;
	typedef std::vector<ValType> ObjectsList; // Objects, if this is the leaf node

public:
	VPNode();
	VPNode(const VPNode& orig);

	VPNode(const ValType& value);

	virtual ~VPNode();

	bool get_leaf_node() const;
	void set_leaf_node(const bool leaf_node);
	
	size_t get_branches_count() const;
	void set_branches_count(size_t new_size);

	size_t get_objects_count() const;

	const VPNodePtr& get_parent() const;
	void set_parent(const VPNodePtr& parent);

	void AddChild(const size_t child_pos, const VPNodePtr& child);

	void AddObject(const ValType& new_object);
	bool DeleteObject(const ValType& object);

	const ValType& get_value() const;
	void set_value(const ValType& new_value);	

public:
	DistanceList m_mu_list;		// "mu"  and child lists for vp-tree like
	ChildList m_child_list;		// |chld| MU |chld| MU |chld|
	VPNodePtr m_parent;

	ObjectsList m_objects_list;		// Objects, if this is the leaf node
	
private:
	bool m_leaf_node;	// Is leaf node
	size_t m_branches_count;	// Number of internal branches
	ValType m_value;	// Storage value
};


template <typename ValType, typename DistanceType>
VPNode<ValType, DistanceType>::VPNode()
: m_leaf_node(false)
, m_branches_count(0)
{
}

template <typename ValType, typename DistanceType>
VPNode<ValType, DistanceType>::VPNode(const VPNode& orig)
: m_leaf_node(false)
, m_branches_count(0)
{
}

template <typename ValType, typename DistanceType>
VPNode<ValType, DistanceType>::VPNode(const ValType& value)
: m_leaf_node(true)
, m_branches_count(0)
, m_value(value)
{
	//m_objects_list.push_back(value);
}

template <typename ValType, typename DistanceType>
VPNode<ValType, DistanceType>::~VPNode()
{
}

template <typename ValType, typename DistanceType>
const ValType& VPNode<ValType, DistanceType>::get_value() const
{
	return m_value;
}

template <typename ValType, typename DistanceType>
void VPNode<ValType, DistanceType>::set_value(const ValType& new_value)
{
	m_value = new_value;
}

template <typename ValType, typename DistanceType>
bool VPNode<ValType, DistanceType>::get_leaf_node() const
{
	return m_leaf_node;
}

template <typename ValType, typename DistanceType>
void VPNode<ValType, DistanceType>::set_leaf_node(const bool leaf_node)
{
	m_leaf_node = leaf_node;
}

template <typename ValType, typename DistanceType>
size_t VPNode<ValType, DistanceType>::get_branches_count() const
{
	return m_branches_count;
}

template <typename ValType, typename DistanceType>
void VPNode<ValType, DistanceType>::set_branches_count(size_t new_size)
{
	m_branches_count = new_size;
}

template <typename ValType, typename DistanceType>
size_t VPNode<ValType, DistanceType>::get_objects_count() const
{
	return m_objects_list.size();
}

template <typename ValType, typename DistanceType>
const typename VPNode<ValType, DistanceType>::VPNodePtr& VPNode<ValType, DistanceType>::get_parent() const
{
	return m_parent;
}

template <typename ValType, typename DistanceType>
void VPNode<ValType, DistanceType>::set_parent(const VPNodePtr& parent)
{
	m_parent = parent;
}

template <typename ValType, typename DistanceType>
void VPNode<ValType, DistanceType>::AddChild(const size_t /*child_pos*/, const VPNodePtr& child)
{
	/*
	if(m_child_list.size() <= child_pos)
	{
		m_child_list.reserve(child_pos + 1);
		m_mu_list.reserve(child_pos + 1);
	}
	*/
	//m_child_list[child_pos] = child;
	m_child_list.push_back(child);

	//(*m_child_list.rbegin())->set_parent(shared_from_this());

	if(get_branches_count() != 1)
		m_mu_list.push_back(static_cast<DistanceType>(0));

	++m_branches_count;
}

template <typename ValType, typename DistanceType>
void VPNode<ValType, DistanceType>::AddObject(const ValType& new_object)
{
	m_leaf_node = true;
	m_objects_list.push_back(new_object);
}

template <typename ValType, typename DistanceType>
bool VPNode<ValType, DistanceType>::DeleteObject(const ValType& object)
{
	typename ObjectsList::iterator it_find = std::find(m_objects_list.begin(),
													m_objects_list.end(),
													object);
	if(it_find != m_objects_list.end())
	{
		m_objects_list.erase(it_find);
		return true;
	}
	
	return false;
}



} // namespace vptree

#endif	/* _VPNODE_H */

