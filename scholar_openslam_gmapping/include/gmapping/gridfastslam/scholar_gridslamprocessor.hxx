#ifdef MACOSX
// This is to overcome a possible bug in Apple's GCC.
#define isnan(x) (x==FP_NAN)
#endif

inline void GridSlamProcessor::scanMatch(const double* plainReading)
{
  
  for (ParticleVector::iterator it=m_particles.begin(); it!=m_particles.end(); it++)
  {
    OrientedPoint corrected;
    double score, l, s;
    score=m_matcher.optimize(corrected, *m_map, it->pose, plainReading);//根据map打分

    if (score>m_minimumScore)
    {
      it->pose=corrected;
    } 
    else
    {
    }

    m_matcher.likelihoodAndScore(s, l, *m_map, it->pose, plainReading);//跟map做似然估计

    it->weight+=l;
    it->weightSum+=l;

    m_matcher.invalidateActiveArea();
    m_matcher.computeActiveArea(*m_map, it->pose, plainReading);//根据map计算activearea
    
  }
  
}

inline void GridSlamProcessor::normalize()
{
  
  double gain=1./(m_obsSigmaGain*m_particles.size());
  double lmax= -std::numeric_limits<double>::max();
  for (ParticleVector::iterator it=m_particles.begin(); it!=m_particles.end(); it++)
  {
    lmax=it->weight>lmax?it->weight:lmax;
  }
  
  m_weights.clear();
  double wcum=0;
  m_neff=0;
  for (std::vector<Particle>::iterator it=m_particles.begin(); it!=m_particles.end(); it++)
  {
    m_weights.push_back(exp(gain*(it->weight-lmax)));
    wcum+=m_weights.back();
  }
  
  m_neff=0;
  for (std::vector<double>::iterator it=m_weights.begin(); it!=m_weights.end(); it++)
  {
    *it=*it/wcum;
    double w=*it;
    m_neff+=w*w;
  }
  m_neff=1./m_neff;
  
}

inline bool GridSlamProcessor::resample(const double* plainReading, int adaptSize, const RangeReading* reading)
{
  
  bool hasResampled = false;
  
  TNodeVector oldGeneration;
  for (unsigned int i=0; i<m_particles.size(); i++)
  {
    oldGeneration.push_back(m_particles[i].node);
  }

  if (m_neff<m_resampleThreshold*m_particles.size())
  {		    
    
    uniform_resampler<double, double> resampler;
    m_indexes=resampler.resampleIndexes(m_weights, adaptSize);
    
    ParticleVector temp;
    unsigned int j=0;
    std::vector<unsigned int> deletedParticles;  	
    
    for (unsigned int i=0; i<m_indexes.size(); i++)
    {
      while(j<m_indexes[i])
      {
	      deletedParticles.push_back(j);
	      j++;
			}
      if (j==m_indexes[i])j++;

      Particle & p=m_particles[m_indexes[i]];
      TNode* node=0;
      TNode* oldNode=oldGeneration[m_indexes[i]];
      node=new	TNode(p.pose, 0, oldNode, 0);
      node->reading=reading;
      
      temp.push_back(p);
      temp.back().node=node;
      temp.back().previousIndex=m_indexes[i];
    }

    while(j<m_indexes.size())
    {
      deletedParticles.push_back(j);
      j++;
    }
    for (unsigned int i=0; i<deletedParticles.size(); i++)
    {
      delete m_particles[deletedParticles[i]].node;
      m_particles[deletedParticles[i]].node=0;
    }
    
    m_particles.clear();

    for (ParticleVector::iterator it=temp.begin(); it!=temp.end(); it++)
    {
      it->setWeight(0);
      m_matcher.invalidateActiveArea();
      m_matcher.registerScan(*m_map, it->pose, plainReading);
      m_particles.push_back(*it);
    }
    hasResampled = true;
  } 
  else 
  {
    int index=0;
    
    TNodeVector::iterator node_it=oldGeneration.begin();
    for (ParticleVector::iterator it=m_particles.begin(); it!=m_particles.end(); it++)
    {

      TNode* node=0;
      node=new TNode(it->pose, 0.0, *node_it, 0);
      
      node->reading=reading;
      it->node=node;

      m_matcher.invalidateActiveArea();
      m_matcher.computeActiveArea(*m_map, it->pose, plainReading);

      m_matcher.registerScan(*m_map, it->pose, plainReading);
      it->previousIndex=index;
      index++;
      node_it++;      
    }    
  } 

  //更新地图
 Particle best =    getParticles()[getBestParticleIndex()];  
for(TNode* n = best.node;  n;  n = n->parent)
  {
    if(!n->reading)
    {
      continue;
    }
    m_matcher.invalidateActiveArea();
    m_matcher.computeActiveArea(*m_map, n->pose, &((*n->reading)[0]));
    m_matcher.registerScan(*m_map, n->pose, &((*n->reading)[0]));
  }


  return hasResampled;
  
}
