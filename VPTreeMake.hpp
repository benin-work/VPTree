/* 
 * File:   VPTreeMake.hpp
 * Author: V.Karlov
 *
 * Created on 19.04.2010
 */

#ifndef _VPTREEMAKE_HPP
#define	_VPTREEMAKE_HPP

#include "VPTree.hpp"

namespace vptree {


template<	typename DataType,		// VPVector data type
			typename ValType,		// VPNode value as VPVector<>
			typename DistanceType,	// Distance type
			typename DistanceObtainer	// Distance type obtainer
>
typename VPTree<DataType, ValType, DistanceType, DistanceObtainer>::VPNodePtr
VPTree<DataType, ValType, DistanceType, DistanceObtainer>::MakeVPTree(const ObjectsList& objects, const VPNodePtr& parent)
{

	if(objects.empty())
		return VPNodePtr(new VPNodeType);

	VPNodePtr new_node(new VPNodeType);

	new_node->set_parent(parent);

	// Set the VP
	new_node->set_value(SelectVP(objects));

	if(objects.size() <= m_non_leaf_branching_factor * m_leaf_branching_factor)
	{
		for(size_t c_pos = 0; c_pos < m_leaf_branching_factor; ++c_pos)
		{
			new_node->AddChild(0, VPNodePtr(new VPNodeType));
			new_node->m_child_list[c_pos]->set_leaf_node(true);
			new_node->m_child_list[c_pos]->set_parent(new_node);
		}

		new_node->m_child_list[0]->m_objects_list.insert(new_node->m_child_list[0]->m_objects_list.begin(), objects.begin()+1, objects.end());

		RedistributeAmongLeafNodes(new_node, *objects.begin());

		return new_node;
	}

	// Init children
	new_node->AddChild(0, VPNodePtr(new VPNodeType));
	new_node->AddChild(0, VPNodePtr(new VPNodeType));

	DistanceType median = Median(new_node->get_value(), objects.begin(), objects.end());
	new_node->m_mu_list[0] = median;

	size_t objects_count = objects.size();
	if(median == 0)
		objects_count = 0;

	bool c_left = false;

	// 60% of size
	size_t reserved_memory = static_cast<size_t>(static_cast<double>(objects_count) * 0.6);
	ObjectsList s_left, s_right;
	s_left.reserve(reserved_memory);
	s_right.reserve(reserved_memory);

	typename ObjectsList::const_iterator it_obj = objects.begin();
	while(it_obj != objects.end())
	{
		DistanceType dist = m_get_distance(new_node->get_value(), *it_obj);
		if(dist < new_node->m_mu_list[0] || (dist == 0  && !c_left))
		{
			s_left.push_back(*it_obj);
			c_left = true;
		}else
		{
			s_right.push_back(*it_obj);
			c_left = false;
		}
		++it_obj;
	}

	size_t left_count = s_left.size();
	size_t right_count = s_right.size();

	// 8( for 2 only now
	new_node->set_branches_count(2);

	VPNodePtr new_node_l(new VPNodeType);
	VPNodePtr new_node_r(new VPNodeType);

	VPNodePtr new_node_lc(new VPNodeType);
	VPNodePtr new_node_rc(new VPNodeType);

	#pragma omp task shared(new_node)
		new_node->m_child_list[0] = MakeVPTree(s_left, new_node);
	#pragma omp task shared(new_node)
		new_node->m_child_list[1] = MakeVPTree(s_right, new_node);

	#pragma omp taskwait

	return new_node;
}

// Obtain the best vantage point for objects set
// now - just take the first data point
template<	typename DataType,		// VPVector data type
			typename ValType,		// VPNode value as VPVector<>
			typename DistanceType,	// Distance type
			typename DistanceObtainer	// Distance type obtainer
>
const ValType& VPTree<DataType, ValType, DistanceType, DistanceObtainer>::SelectVP(const ObjectsList& objects)
{
	assert(!objects.empty());

	return *objects.begin();
}
	

} // namespace vptree

#endif	/* _VPTREEMAKE_HPP */

