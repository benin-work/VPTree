/* 
 * File:   VPTree.hpp
  * Author: V.Karlov
 *
 * Created on 19.04.2010
 */

#ifndef _VPTREE_H
#define	_VPTREE_H

#include "types.h"
#include "VPNode.hpp"
#include "VPVector.hpp"
#include <map>
#include <float.h>

#include <iomanip>
#include <iostream>

using namespace std;

namespace {

	template <typename ValType, typename DistanceObtainer>
	class ValueSorter
	{
	public:
		ValueSorter(const ValType& main_val, const DistanceObtainer& distancer)
		: m_main_val(main_val)
		, m_get_distance(distancer)
		{
		}

		bool operator()(const ValType& val1, const ValType& val2)
		{
			return m_get_distance(m_main_val, val1) < m_get_distance(m_main_val, val2);
		}
	private:
		const ValType& m_main_val;
		const DistanceObtainer& m_get_distance;
	};
}

namespace vptree {

struct TreeStatistic{
	uint64_t distance_count;
	uint64_t distance_threshold_count;
	uint64_t search_jump;

	void clear()
	{
		distance_count = 0;
		distance_threshold_count = 0;
		search_jump = 0;
	}

	TreeStatistic()
	: distance_count(0)
	, distance_threshold_count(0)
	, search_jump(0)
	{
	}
};

template<	typename DataType,		// VPVector data type
			typename ValType,		// VPNode value as VPVector<>
			typename DistanceType,	// Distance type
			typename DistanceObtainer	// Distance type obtainer
>
class VPTree
{
public:
    typedef VPNode<ValType, DistanceType> VPNodeType;
    typedef boost::shared_ptr<VPNodeType> VPNodePtr;

    typedef std::vector<ValType> ObjectsList;	// Objects, at leaf node

    typedef std::pair<DistanceType, ValType> SearchResult;
    typedef std::multimap<DistanceType, ValType> SearchResultMap;

    typedef ValueSorter<ValType, DistanceObtainer> ValueSorterType;
public:

    VPTree(const size_t non_leaf_branching_factor = 2, const size_t leaf_branching_factor = 2)
    : m_root(new VPNodeType())
    , m_non_leaf_branching_factor(non_leaf_branching_factor)
    , m_leaf_branching_factor(leaf_branching_factor)
	, m_out(NULL)
    {
    }


    VPTree(const VPTree& orig)
    : m_non_leaf_branching_factor(3)
    , m_leaf_branching_factor(3)
	, m_out(NULL)
    {
		if(&orig == this)
			return;

		m_root = VPNodePtr(new VPNodeType());
    }

    virtual ~VPTree()
    {
    }

    void set_non_leaf_branching_factor(const size_t branching_factor)
    {
		m_non_leaf_branching_factor = branching_factor;
    }

    void set_leaf_branching_factor(const size_t branching_factor)
    {
		m_leaf_branching_factor = branching_factor;
    }

    const VPNodePtr& get_root() const
    {
		return m_root;
    }

	void set_info_output(std::ostream* stream)
	{
		m_out = stream;
	}

	// Construct the tree from given objects set
	void MakeVPTree(const ObjectsList& objects)
	{
		m_root.reset();

		#pragma omp parallel
		{
			#pragma omp single
				#pragma omp task
					m_root = MakeVPTree(objects, m_root);
		}
	}

    void Insert(const ValType& new_value)
    {
		if(m_out) *m_out << "Insert data to vptree.root" << std::endl;
		Insert(new_value, m_root);

		//draw(m_root, 0);
    }

    void Insert(const ValType& new_value, VPNodePtr& root)
    {
		assert(root.get());
		// 4.  L is the empty root Node
		if(!root->get_branches_count() && !root->get_leaf_node())
		{
				// At first insertion make root as a leaf node
				//InsertSplitRoot(new_value)
				root->AddObject(new_value);
		}else
		{
			// Traverse the tree, choosing the subtree Si, until the L(eaf) node
			// is found

			size_t c_current_node_parent_id = 0;
			VPNodePtr c_parent_node;
			// Go through the tree, searching leaf node
			VPNodePtr c_current_node = c_parent_node = root;
			while(!c_current_node->get_leaf_node() && c_current_node.get())
			{
				c_parent_node = c_current_node;
				// test all distances at node
				for(size_t c_pos = 0; c_pos < c_current_node->m_mu_list.size(); ++c_pos)
				{
					// test new_value with node vantage point
					if( m_get_distance(new_value, c_current_node->get_value()) < c_current_node->m_mu_list[c_pos])
					{
						c_current_node = c_current_node->m_child_list[c_pos];
						c_current_node_parent_id = c_pos;
						break;
					}
				}

				if(c_parent_node == c_current_node )
					c_current_node = *c_current_node->m_child_list.rbegin();
			}

			// Assume c_current_node - Leaf node
			// Have found leaf node, analize ancestros

			// 0. If there is a room at L(eaf) node - insert data
			if(c_current_node->get_objects_count() < m_leaf_branching_factor)
			{
				c_current_node->AddObject(new_value);
				return;
			}

			// Second node - we split the root
			if(c_current_node == root && c_current_node->get_objects_count() >= m_leaf_branching_factor)
			{
				InsertSplitLeafRoot(c_current_node, new_value);
				return;
			}

			// 1. If any sibling leaf node of L(eaf) is not full,
			// redistribute all objects under P(arent), among the leaf nodes
			// Analize sibling nodes
			for(size_t c_pos = 0; c_pos < c_parent_node->get_branches_count(); ++c_pos)
			{
				if(c_parent_node->m_child_list[c_pos]->get_objects_count() < m_leaf_branching_factor)
				{
					RedistributeAmongLeafNodes(c_parent_node, new_value);
					return;
				}
			}


			// 2. If Parent has a room for one more child - split the leaf node
			if (c_parent_node->get_branches_count() < m_non_leaf_branching_factor)
			{
				SplitLeafNode(c_parent_node, c_current_node_parent_id, new_value);
				return;
			}

			// 3.a. Redistribute, among the sibling subtrees
			VPNodePtr c_ancestor = c_parent_node->m_parent;
			if(c_ancestor.get())
			{
				// found an id of full leaf node parent
				size_t c_found_free_subtree_id = m_leaf_branching_factor;
				size_t c_full_subtree_id = m_leaf_branching_factor;
				for(size_t c_anc_pos = 0; c_anc_pos < c_ancestor->get_branches_count(); ++c_anc_pos)
					if(c_ancestor->m_child_list[c_anc_pos] == c_current_node->m_parent)
						c_full_subtree_id = c_anc_pos;

				//assert(c_full_subtree_id != m_leaf_branching_factor);

				if(c_full_subtree_id != m_leaf_branching_factor)
					for(size_t c_anc_pos = 0; c_anc_pos < c_ancestor->get_branches_count(); ++c_anc_pos)
					{
						VPNodePtr c_parent = c_ancestor->m_child_list[c_anc_pos];

						if(c_parent == c_current_node->m_parent)
							continue;

						for(size_t c_par_pos = 0; c_par_pos < c_parent->get_branches_count(); ++c_par_pos)
						{
							if(c_parent->m_child_list[c_par_pos]->get_leaf_node() &&
									c_parent->m_child_list[c_par_pos]->m_objects_list.size() < m_leaf_branching_factor)
							{
								c_found_free_subtree_id = c_anc_pos;
								break;
							}
						}
						if(c_found_free_subtree_id < m_leaf_branching_factor)
						{
							// Found free subtree - redistribute data
							if(c_found_free_subtree_id > c_full_subtree_id)
								RedistributeAmongNonLeafNodes(c_ancestor, c_full_subtree_id, c_found_free_subtree_id, new_value);
							else
								RedistributeAmongNonLeafNodes(c_ancestor, c_found_free_subtree_id, c_full_subtree_id, new_value);

							Insert(new_value, c_ancestor);

							return;
						}
					}
			}

			// 3.b. If Parent-Parent node is not full, spleat non-leaf node
			if(c_current_node->m_parent.get() && c_current_node->m_parent->m_parent.get())
			{
				VPNodePtr c_ancestor = c_current_node->m_parent->m_parent; // A
				VPNodePtr c_current_parent = c_current_node->m_parent; // B

				size_t c_found_free_subtree_id = m_non_leaf_branching_factor;
				for(size_t c_pos = 0; c_pos < c_ancestor->get_branches_count(); ++c_pos)
					if(c_ancestor->m_child_list[c_pos] == c_current_parent)
						c_found_free_subtree_id = c_pos;

				if(c_found_free_subtree_id != m_non_leaf_branching_factor &&
						c_ancestor->get_branches_count() < m_non_leaf_branching_factor)
				{
					SplitNonLeafNode(c_ancestor, c_found_free_subtree_id, new_value);
					return;
				}
			}

			// 4. Cannot find any ancestor, that is not full -> split the root
			// into two new nodes s1 and s2
			InsertSplitRoot(root, new_value);

		}

    }

	void Remove(const ValType& query_value)
	{
		Remove(query_value, m_root);
	}

	size_t get_object_count() const
	{
		return get_object_count(m_root);
	}

    // Nearest Neighbor Search
    void NNSearch(const ValType& query_value, const size_t count, SearchResultMap& result_found)
    {
        DistanceType cq = static_cast<DistanceType>(FLT_MAX);
		//DistanceType cq = 100;

        Search(query_value, count, result_found, m_root, cq);
    }

    void ClearStatistic()
    {
    	m_stat.clear();
    }

    const TreeStatistic& statistic() const
    {
    	return m_stat;
    }

private:

	void Remove(const ValType& query_value, const VPNodePtr& node)
	{
		if(node->get_leaf_node())
			node->DeleteObject(query_value);
		else
		{
			for (size_t c_pos = 0; c_pos < node->get_branches_count(); ++c_pos)
				Remove(query_value, node->m_child_list[c_pos]);
		}
	}

	size_t get_object_count(const VPNodePtr& node) const
	{
		size_t c_count = 0;
		get_object_count(c_count, node);
		return c_count;
	}

	void get_object_count(size_t& obj_count, const VPNodePtr& node) const
	{
		if(node->get_leaf_node())
			obj_count += node->get_objects_count();
		else
		{
			for (size_t c_pos = 0; c_pos < node->get_branches_count(); ++c_pos)
				get_object_count(obj_count, node->m_child_list[c_pos]);
		}
	}

    void Search(const ValType& query_value, const size_t count, SearchResultMap& result_found,
            const VPNodePtr& node, DistanceType& q)
    {
        assert(node.get());

        if(node->get_leaf_node())
        {
            for(size_t c_pos = 0; c_pos < node->m_objects_list.size(); ++c_pos)
            {
            	m_stat.distance_threshold_count++;
				DistanceType c_distance = m_get_distance(query_value, node->m_objects_list[c_pos], q);
                if( c_distance <= q)
                {
					result_found.insert(SearchResult(c_distance, node->m_objects_list[c_pos]));

					while(result_found.size() > count)
					{
						typename SearchResultMap::iterator it_last = result_found.end();
						
						result_found.erase(--it_last);
					}

					if(result_found.size() == count)
						q = (*result_found.rbegin()).first;
                }
            }

        }else
        {
            DistanceType dist = 0; //m_get_distance(node->get_value(), query_value);

			// Search flag
			size_t c_mu_pos = m_non_leaf_branching_factor;
			if(node->m_mu_list.size() == 1)
			{
				c_mu_pos = 0;
				m_stat.distance_threshold_count++;
				dist = m_get_distance(node->get_value(), query_value, node->m_mu_list[c_mu_pos] + q);
			}else
			{
				m_stat.distance_count++;
				dist = m_get_distance(node->get_value(), query_value);
				for(size_t c_pos = 0; c_pos < node->m_mu_list.size() -1 ; ++c_pos)
				{
					if(dist > node->m_mu_list[c_pos] && dist < node->m_mu_list[c_pos + 1] )
					{
						c_mu_pos = c_pos;
						break;
					}
				}
			}
				
			if(c_mu_pos != m_non_leaf_branching_factor)
			{
				DistanceType c_mu = node->m_mu_list[c_mu_pos];
				if(dist < c_mu)
				{
					if(dist < c_mu + q)
					{
						m_stat.search_jump++;
						Search(query_value, count, result_found, node->m_child_list[c_mu_pos], q);
					}
					if(dist >= c_mu - q)
					{
						m_stat.search_jump++;
						Search(query_value, count, result_found, node->m_child_list[c_mu_pos + 1], q);
					}
				}else
				{
					if(dist >= c_mu - q)
					{
						m_stat.search_jump++;
						Search(query_value, count, result_found, node->m_child_list[c_mu_pos + 1], q);
					}
					if(dist < c_mu + q)
					{
						m_stat.search_jump++;
						Search(query_value, count, result_found, node->m_child_list[c_mu_pos], q);
					}
				}
			}
        }

    }

    void InsertSplitLeafRoot(VPNodePtr& root, const ValType& new_value)
    {
		if(m_out) *m_out << "	spit leaf root" << std::endl;
		// Split the root node if root is the leaf
		//
		VPNodePtr s1(new VPNodeType);
		VPNodePtr s2(new VPNodeType);

		// Set vantage point to root
		root->set_value(root->m_objects_list[0]);
		//root->m_objects_list.clear();

		//s1->AddObject(root->get_value());

		root->AddChild(0, s1);
		s1->set_parent(root);
		s1->set_leaf_node(true);

		root->AddChild(1, s2);
		s2->set_parent(root);
		s2->set_leaf_node(true);

		root->set_leaf_node(false);

		for(size_t c_pos = 0; c_pos < root->get_objects_count(); ++c_pos)
			Insert(root->m_objects_list[c_pos], root);

		root->m_objects_list.clear();

		Insert(new_value, root);

		//RedistributeAmongLeafNodes(root, new_value);

		//m_root->m_mu_list[0] = 0;
		//m_root->set_value(new_value); // Set Vantage Point
    }

    void InsertSplitRoot(VPNodePtr& root, const ValType& new_value)
    {
		if(m_out) *m_out << "	split root" << std::endl;
		// Split the root node into 2 new nodes s1 and s2 and insert new data
		// according to the strategy SplitLeafNode() or RedistributeAmongLeafNodes()
		
		VPNodePtr new_root(new VPNodeType);
		VPNodePtr s2_node(new VPNodeType);

		new_root->set_value(root->get_value());
		//new_root->set_value(new_value);
		new_root->AddChild(0, root);
		//new_root->AddChild(0, s2_node);

		root->set_parent(new_root);
		//s2_node->set_parent(new_root);
		//s2_node->set_leaf_node(true);

		root = new_root;

		//Insert(new_value, root);
		//Insert(new_value);

		SplitNonLeafNode(root, 0, new_value);
    }

    // If any sibling leaf node of L(eaf) is not full,
    // redistribute all objects under P(arent), among the leaf nodes
    void RedistributeAmongLeafNodes(const VPNodePtr& parent_node, const ValType& new_value)
    {
		if(m_out) *m_out << "	redistribute among leaf nodes" << std::endl;
		// F - number of leaf nodes under P(arent)
		// F should be greater then 1
		//size_t F = parent_node->m_child_list.size();
		const size_t F = parent_node->get_branches_count();

		VPNodePtr c_node;
		ObjectsList S; // Set of leaf objects + new one;

		// Create Set of whole objects from leaf nodes
		CollectObjects(parent_node, S);
		S.push_back(new_value);

		// Order the objects in S with respect to their distances from P's vantage point
		ValueSorterType val_sorter(parent_node->get_value(), m_get_distance);
		std::sort(S.begin(), S.end(), val_sorter);

		// Devide S into F groups of equal cardinality
		size_t c_whole_count = S.size();
		typename ObjectsList::const_iterator it_obj = S.begin();
		for (size_t c_pos = 0; c_pos < F; ++c_pos)
		{
			size_t c_equal_count = c_whole_count / (F - c_pos);
			c_whole_count -= c_equal_count;

			c_node = parent_node->m_child_list[c_pos];
			c_node->m_objects_list.clear();

			c_node->m_objects_list.insert(c_node->m_objects_list.begin(),
											it_obj, it_obj + c_equal_count);
			c_node->set_leaf_node(true);
			it_obj += c_equal_count;
		}

		// Update the boundary distance values
		for (size_t c_pos = 0; c_pos < F - 1; ++c_pos)
		{
			const ObjectsList& SS1 = parent_node->m_child_list[c_pos]->m_objects_list;
			const ObjectsList& SS2 = parent_node->m_child_list[c_pos + 1]->m_objects_list;

			parent_node->m_mu_list[c_pos] = MedianSumm(SS1, SS2, parent_node->get_value());
		}

    }

    // If L(eaf) node has a P(arent) node and P has room for one more child,
    // split the leaf node L
    void SplitLeafNode(const VPNodePtr& parent_node, const size_t child_id, const ValType& new_value)
    {
		if(m_out) *m_out << "	split leaf node" << std::endl;
		// F - number of leaf nodes under P(arent)
		//
		const size_t F = parent_node->get_branches_count();
		const size_t k = child_id;

		assert(child_id < parent_node->m_child_list.size());

		VPNodePtr c_leaf_node = parent_node->m_child_list[child_id];

		ObjectsList S = c_leaf_node->m_objects_list; // Set of leaf objects + new one
		S.push_back(new_value);

		// Order the objects in S with respect to their distances from P's vantage point
		ValueSorterType val_sorter(parent_node->get_value(), m_get_distance);
		std::sort(S.begin(), S.end(), val_sorter);

		// Divide S into 2 groups of equal cardinality
		VPNodePtr ss1_node(new VPNodeType);
		VPNodePtr ss2_node = c_leaf_node;
		ss2_node->m_objects_list.clear();

		parent_node->AddChild(parent_node->get_branches_count(), ss1_node);
		ss1_node->set_parent(parent_node);

		size_t c_half_count = S.size() / 2;
		for(size_t c_pos = 0; c_pos < S.size(); ++c_pos)
		{
			if(c_pos < c_half_count)
				ss1_node->AddObject(S[c_pos]);
			else
				ss2_node->AddObject(S[c_pos]);
		}

		// insertion/shift process
		for(size_t c_pos = F-2; c_pos >= k; --c_pos)
		{
			parent_node->m_mu_list[c_pos+1] = parent_node->m_mu_list[c_pos];
			if(!c_pos) // !!! hack :(
				break;
		}

		const ObjectsList& SS1 = ss1_node->m_objects_list;
		const ObjectsList& SS2 = ss2_node->m_objects_list;
		parent_node->m_mu_list[k] = MedianSumm(SS1, SS2, parent_node->get_value());

		// !! --c_pos
		for(size_t c_pos = F-1; c_pos >= k+1; --c_pos)
			parent_node->m_child_list[c_pos + 1] = parent_node->m_child_list[c_pos];

		parent_node->m_child_list[k] = ss1_node;
		parent_node->m_child_list[k + 1] = ss2_node;
    }

	// 3.a. Redistribute, among the sibling subtrees
    void RedistributeAmongNonLeafNodes(const VPNodePtr& parent_node, const size_t k_id,
    										const size_t k1_id, const ValType& new_value)
	{
		if(m_out) *m_out << "	redistribute among nodes(subtrees)" << std::endl;
		assert(k_id != k1_id);

		size_t num_k = get_object_count(parent_node->m_child_list[k_id]);
		size_t num_k1 = get_object_count(parent_node->m_child_list[k1_id]);

		size_t average = (num_k + num_k1) / 2;

		if(num_k > num_k1)
		{
			// Create Set of objects from leaf nodes K-th subtree
			ObjectsList S; // Set of leaf objects + new one;
			CollectObjects(parent_node->m_child_list[k_id], S);
			//S.push_back(new_value);
			ValueSorterType val_sorter(parent_node->get_value(), m_get_distance);
			std::sort(S.begin(), S.end(), val_sorter);

			size_t w = num_k - average;

			ObjectsList SS1(S.begin(), S.begin() + num_k - w);
			ObjectsList SS2(S.begin() + num_k - w, S.end());

			SS1.push_back(new_value);

			typename ObjectsList::const_iterator it_obj = SS2.begin();
			for(;it_obj != SS2.end(); ++it_obj)
				Remove(*it_obj, parent_node->m_child_list[k_id]);

			parent_node->m_mu_list[k_id] = MedianSumm(SS1, SS2, parent_node->get_value());

			for(it_obj = SS2.begin(); it_obj != SS2.end(); ++it_obj)
				Insert(*it_obj, parent_node->m_child_list[k1_id]);

		}else
		{
			// Create Set of objects from leaf nodes K-th subtree
			ObjectsList S; // Set of leaf objects + new one;
			CollectObjects(parent_node->m_child_list[k1_id], S);
			//S.push_back(new_value);
			ValueSorterType val_sorter(parent_node->get_value(), m_get_distance);
			std::sort(S.begin(), S.end(), val_sorter);

			size_t w = num_k1 - average;

			ObjectsList SS1(S.begin(), S.begin() + w);
			ObjectsList SS2(S.begin() + w, S.end());
			SS2.push_back(new_value);

			typename ObjectsList::const_iterator it_obj = SS1.begin();
			for(;it_obj != SS1.end(); ++it_obj)
				Remove(*it_obj, parent_node->m_child_list[k1_id]);

			parent_node->m_mu_list[k_id] = MedianSumm(SS1, SS2, parent_node->get_value());

			for(it_obj = SS1.begin(); it_obj != SS1.end(); ++it_obj)
				Insert(*it_obj, parent_node->m_child_list[k_id]);
		}

		num_k = get_object_count(parent_node->m_child_list[k_id]);
		num_k1 = get_object_count(parent_node->m_child_list[k1_id]);

		num_k = 0;
	}

	// 3.b. If Parent-Parent node is not full, spleat non-leaf node
    void SplitNonLeafNode(const VPNodePtr& parent_node, const size_t child_id, const ValType& new_value)
    {
		if(m_out) *m_out << "	split node" << std::endl;
		assert(child_id < parent_node->m_child_list.size());

		const size_t k = child_id;
		const size_t F = parent_node->get_branches_count();

		VPNodePtr c_split_node = parent_node->m_child_list[child_id];
		VPNodePtr c_node;

		ObjectsList S; // Set of leaf objects + new one;

		// Create Set of whole objects from leaf nodes at sub tree
		CollectObjects(c_split_node, S);
		S.push_back(new_value);

		// Order the objects in S with respect to their distances from P's vantage point
		ValueSorterType val_sorter(parent_node->get_value(), m_get_distance);
		std::sort(S.begin(), S.end(), val_sorter);

		// Create free room for instance
		VPNodePtr s2_node(new VPNodeType);
		s2_node->set_parent(parent_node);
		parent_node->AddChild(0, s2_node);

		ObjectsList SS1(S.begin(), S.begin() + S.size() / 2);
		ObjectsList SS2(S.begin() + S.size() / 2, S.end());


		// Shift data at parent node

		if(F > 1)
			for(size_t c_pos = F-2; c_pos >= k; --c_pos)
			{
				parent_node->m_mu_list[c_pos+1] = parent_node->m_mu_list[c_pos];
				if(!c_pos) // !!! hack :(
					break;
			}

		parent_node->m_mu_list[k] = MedianSumm(SS1, SS2, parent_node->get_value());
		for(size_t c_pos = F-1; c_pos >= k+1; --c_pos)
			parent_node->m_child_list[c_pos + 1] = parent_node->m_child_list[c_pos];


		// Construct new vp-tree
		VPNodePtr ss1_node;
		VPNodePtr ss2_node;

		#pragma omp task shared(ss1_node)
			ss1_node = MakeVPTree(SS1, ss1_node);
		#pragma omp task shared(ss2_node)
			ss2_node = MakeVPTree(SS2, ss2_node);

		#pragma omp taskwait

		ss1_node->set_parent(parent_node);
		ss2_node->set_parent(parent_node);

		parent_node->m_child_list[k] = ss1_node;
		parent_node->m_child_list[k+1] = ss2_node;
    }

    // (max{d(v, sj) | sj c SS1} + min{d(v, sj) | sj c SS2}) / 2
    const DistanceType MedianSumm(const ObjectsList& SS1, const ObjectsList& SS2, const ValType& v) const
    {
		DistanceType c_max_distance = 0;
		DistanceType c_min_distance = 0;
		DistanceType c_current_distance = 0;

		if(!SS1.empty())
		{
			// max{d(v, sj) | sj c SSj}
			typename ObjectsList::const_iterator it_obj = SS1.begin();
			c_max_distance = m_get_distance(*it_obj, v);
			++it_obj;
			c_current_distance = c_max_distance;
			while(it_obj != SS1.end())
			{
				c_current_distance = m_get_distance(*it_obj, v);
				if(c_current_distance > c_max_distance)
					c_max_distance = c_current_distance;

				++it_obj;
			}
		}

		if(!SS2.empty())
		{
			// min{d(v, sj) | sj c SSj}
			typename ObjectsList::const_iterator it_obj = SS2.begin();
			c_min_distance = m_get_distance(*it_obj, v);
			++it_obj;
			c_current_distance = c_min_distance;
			while(it_obj != SS2.end())
			{
				c_current_distance = m_get_distance(*it_obj, v);
				if(c_current_distance < c_min_distance)
					c_min_distance = c_current_distance;

				++it_obj;
			}
		}

		return (c_max_distance + c_min_distance) / static_cast<DistanceType>(2);
    }

    // Recursively collect data from subtree, and push them into S
    void CollectObjects(const VPNodePtr& node, ObjectsList& S)
    {
		if(node->get_leaf_node())
			S.insert(S.end(), node->m_objects_list.begin(), node->m_objects_list.end() );
		else
		{
			for (size_t c_pos = 0; c_pos < node->get_branches_count(); ++c_pos)
				CollectObjects(node->m_child_list[c_pos], S);
		}
    }

	// Calc the median value for object set
	DistanceType Median(const ValType& value, const typename ObjectsList::const_iterator it_begin,
												const typename ObjectsList::const_iterator it_end)
	{
		typename ObjectsList::const_iterator it_obj = it_begin;
		DistanceType current_distance = 0;
		size_t count = 0;
		while(it_obj != it_end)
		{
			current_distance += m_get_distance(*it_obj, value);
			++it_obj;
			++count;
		}
		return current_distance / static_cast<DistanceType>(count);
	}

	// Construct the tree from given objects set
	VPNodePtr MakeVPTree(const ObjectsList& objects, const VPNodePtr& parent);

	// Obtain the best vantage point for objects set
	const ValType& SelectVP(const ObjectsList& objects);	
	
private:
    VPNodePtr m_root;
    DistanceObtainer m_get_distance;

    size_t m_non_leaf_branching_factor;
    size_t m_leaf_branching_factor;

	std::ostream* m_out;

	// Statistic
	TreeStatistic m_stat;
};



} // namespace vptree

#include "VPTreeMake.hpp"

#endif	/* _VPTREE_H */

