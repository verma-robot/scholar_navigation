#include <string>
#include <deque>
#include <list>
#include <map>
#include <set>
#include <fstream>

#include <gmapping/utils/stat.h>
#include "gmapping/gridfastslam/scholar_gridslamprocessor.h"

namespace Scholar_GMapping {

using namespace std;

GridSlamProcessor::TNode::TNode(const OrientedPoint& p, double w, TNode* n, unsigned int c)
{
	pose=p;
	weight=w;
	childs=c;
	parent=n;
	reading=0;
	gweight=0;
	if (n)
	{
		n->childs++;
	}
	flag=0;
	accWeight=0;
}


GridSlamProcessor::TNode::~TNode()
{
	if (parent && (--parent->childs)<=0)delete parent;
	assert(!childs);
}



GridSlamProcessor::TNodeVector GridSlamProcessor::getTrajectories() const
{
	
  TNodeVector v;
  TNodeMultimap parentCache;
  TNodeDeque border;
	
	for (ParticleVector::const_iterator it=m_particles.begin(); it!=m_particles.end(); it++)
	{
		TNode* node=it->node;
		while(node)
		{
			node->flag=false;
			node=node->parent;
		}
	}
	
	for (ParticleVector::const_iterator it=m_particles.begin(); it!=m_particles.end(); it++)
	{
		TNode* newnode=new TNode(* (it->node) );
		
		v.push_back(newnode);
		assert(newnode->childs==0);
		if (newnode->parent)
		{
			parentCache.insert(make_pair(newnode->parent, newnode));
			if (! newnode->parent->flag)
			{
				newnode->parent->flag=true;
				border.push_back(newnode->parent);
			}
		}
	}	
	while (! border.empty())
	{
		const TNode* node=border.front();
		border.pop_front();
		if (! node)	continue;
			
		TNode* newnode=new TNode(*node);
		node->flag=false;		
		pair<TNodeMultimap::iterator, TNodeMultimap::iterator> p=parentCache.equal_range(node);
		double childs=0;
		for (TNodeMultimap::iterator it=p.first; it!=p.second; it++)
		{
			assert(it->second->parent==it->first);
			(it->second)->parent=newnode;
			childs++;
		}
		parentCache.erase(p.first, p.second);
		assert(childs==newnode->childs);
		
		if ( node->parent ){
			parentCache.insert(make_pair(node->parent, newnode));
			if(! node->parent->flag)
			{
				border.push_back(node->parent);
				node->parent->flag=true;
			}	
		}
	}
	for (unsigned int i=0; i<v.size(); i++){
		TNode* node= v[i];
		while (node)
		{
			node=node->parent;
		}
	}		
	return v;

}
void  GridSlamProcessor::updateTreeWeights(bool weightsAlreadyNormalized)
{

  if (!weightsAlreadyNormalized)
  {
    normalize();
  }
  resetTree();
  propagateWeights();
  
}

void GridSlamProcessor::resetTree()
{

	for (ParticleVector::iterator it=m_particles.begin(); it!=m_particles.end(); it++){
		TNode* n=it->node;
		while (n){
			n->accWeight=0;
			n->visitCounter=0;
			n=n->parent;
		}
	}
	
}

double propagateWeight(GridSlamProcessor::TNode* n, double weight)
{
	
	if (!n)	return weight;
	double w=0;
	n->visitCounter++;
	n->accWeight+=weight;
	if (n->visitCounter==n->childs)
	{
		w=propagateWeight(n->parent,n->accWeight);
	}
	assert(n->visitCounter<=n->childs);
	return w;
	
}

double GridSlamProcessor::propagateWeights()
{
	
	double lastNodeWeight=0;
	double aw=0;                   

	std::vector<double>::iterator w=m_weights.begin();
	for (ParticleVector::iterator it=m_particles.begin(); it!=m_particles.end(); it++)
	{
		double weight=*w;
		aw+=weight;
		TNode * n=it->node;
		n->accWeight=weight;
		lastNodeWeight+=propagateWeight(n->parent,n->accWeight);
		w++;
	}
	
	if (fabs(aw-1.0) > 0.0001 || fabs(lastNodeWeight-1.0) > 0.0001)
	{
	  cerr << "ERROR: ";
	  cerr << "root->accWeight=" << lastNodeWeight << "    sum_leaf_weights=" << aw << endl;
	  assert(0);         
	}
	return lastNodeWeight;
	
}

};

//END
