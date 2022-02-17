#ifndef PARTICLEFILTER_H
#define PARTICLEFILTER_H
#include <stdlib.h>
#include <float.h>
#include <sys/types.h>
#include<vector>
#include<utility>
#include<cmath>
#include<gmapping/utils/gvalues.h>


typedef std::pair<uint,uint> UIntPair;

//BEGIN legacy
template <class Particle, class Numeric>
struct uniform_resampler{
	std::vector<unsigned int> resampleIndexes(const std::vector<Particle> & particles, int nparticles=0) const;
	std::vector<Particle> resample(const std::vector<Particle> & particles, int nparticles=0) const;
	Numeric neff(const std::vector<Particle> & particles) const;
};

template <class Particle, class Numeric>
std::vector<unsigned int> uniform_resampler<Particle, Numeric>:: resampleIndexes(const std::vector<Particle>& particles, int nparticles) const{
	Numeric cweight=0;

	//compute the cumulative weights
	unsigned int n=0;
	for (typename std::vector<Particle>::const_iterator it=particles.begin(); it!=particles.end(); ++it){
		cweight+=(Numeric)*it;
		n++;
	}

	if (nparticles>0)
		n=nparticles;
	
	Numeric interval=cweight/n;

	Numeric target=interval*::drand48();

	cweight=0;
	std::vector<unsigned int> indexes(n);
	n=0;
	unsigned int i=0;
	for (typename std::vector<Particle>::const_iterator it=particles.begin(); it!=particles.end(); ++it, ++i){
		cweight+=(Numeric)* it;
		while(cweight>target){
			indexes[n++]=i;
			target+=interval;
		}
	}
	return indexes;
}

template <class Particle, class Numeric>
std::vector<Particle> uniform_resampler<Particle,Numeric>::resample
  (const typename std::vector<Particle>& particles, int nparticles) const{
	Numeric cweight=0;

	//compute the cumulative weights
	unsigned int n=0;
	for (typename std::vector<Particle>::const_iterator it=particles.begin(); it!=particles.end(); ++it){
		cweight+=(Numeric)*it;
		n++;
	}

	if (nparticles>0)
		n=nparticles;
	
	//weight of the particles after resampling
	double uw=1./n;

	//compute the interval
	Numeric interval=cweight/n;

	//compute the initial target weight
	Numeric target=interval*::drand48();
	//compute the resampled indexes

	cweight=0;
	std::vector<Particle> resampled;
	n=0;
	unsigned int i=0;
	for (typename std::vector<Particle>::const_iterator it=particles.begin(); it!=particles.end(); ++it, ++i){
		cweight+=(Numeric)*it;
		while(cweight>target){
			resampled.push_back(*it);
			resampled.back().setWeight(uw);
			target+=interval;
		}
	}
	return resampled;
}

template <class Particle, class Numeric>
Numeric uniform_resampler<Particle,Numeric>::neff(const std::vector<Particle> & particles) const{
	double cum=0;
	double sum=0;
	for (typename std::vector<Particle>::const_iterator it=particles.begin(); it!=particles.end(); ++it){
		Numeric w=(Numeric)*it;
		cum+=w*w;
		sum+=w;
	}
	return sum*sum/cum;
}



template <class Particle, class EvolutionModel>
struct evolver{
	EvolutionModel evolutionModel;
	void evolve(std::vector<Particle>& particles);
	void evolve(std::vector<Particle>& dest, const std::vector<Particle>& src);
};

template <class Particle, class EvolutionModel>
void evolver<Particle, EvolutionModel>::evolve(std::vector<Particle>& particles){
	for (typename std::vector<Particle>::iterator it=particles.begin(); it!=particles.end(); ++it){
		*it=evolutionModel.evolve(*it);
	}
}

template <class Particle, class EvolutionModel>
void evolver<Particle, EvolutionModel>::evolve(std::vector<Particle>& dest, const std::vector<Particle>& src){
	dest.clear();
	for (typename std::vector<Particle>::const_iterator it=src.begin(); it!=src.end(); ++it)
		dest.push_back(evolutionModel.evolve(*it));
}


template <class Particle, class Numeric, class QualificationModel, class EvolutionModel, class LikelyhoodModel>
struct auxiliary_evolver{
	EvolutionModel evolutionModel;
	QualificationModel qualificationModel;
	LikelyhoodModel likelyhoodModel;
	void evolve(std::vector<Particle>& particles);
	void evolve(std::vector<Particle>& dest, const std::vector<Particle>& src);
};

template <class Particle, class Numeric, class QualificationModel, class EvolutionModel, class LikelyhoodModel>
void auxiliary_evolver<Particle, Numeric, QualificationModel, EvolutionModel, LikelyhoodModel>::evolve
  (std::vector<Particle>&particles){
	std::vector<Numeric> observationWeights(particles.size());
	unsigned int i=0;
	for (typename std::vector<Particle>::const_iterator it=particles.begin(); it!=particles.end(); ++it, i++){
		observationWeights[i]=likelyhoodModel.likelyhood(qualificationModel.evolve(*it));
	}
	uniform_resampler<Numeric, Numeric> resampler;
	std::vector<unsigned int> indexes(resampler.resampleIndexes(observationWeights));
	for (typename std::vector<unsigned int>::const_iterator it=indexes.begin(); it!=indexes.end(); it++){
		Particle & particle=particles[*it];
		particle=evolutionModel.evolve(particle);
		particle.setWeight(likelyhoodModel.likelyhood(particle)/observationWeights[*it]);
	}
}

template <class Particle, class Numeric, class QualificationModel, class EvolutionModel, class LikelyhoodModel>
void auxiliary_evolver<Particle, Numeric, QualificationModel, EvolutionModel, LikelyhoodModel>::evolve
  (std::vector<Particle>& dest, const std::vector<Particle>& src){
	dest.clear();
	std::vector<Numeric> observationWeights(src.size());
	unsigned int i=0;
	for (typename std::vector<Particle>::const_iterator it=src.begin(); it!=src.end(); ++it, i++){
		observationWeights[i]=likelyhoodModel.likelyhood(qualificationModel.evolve(*it));
	}
	uniform_resampler<Numeric, Numeric> resampler;
	std::vector<unsigned int> indexes(resampler.resampleIndexes(observationWeights));
	for (typename std::vector<unsigned int>::const_iterator it=indexes.begin(); it!=indexes.end(); it++){
		Particle & particle=src[*it];
		dest.push_back(evolutionModel.evolve(particle));
		dest.back().weight*=likelyhoodModel.likelyhood(particle)/observationWeights[*it];
	}
}
//END legacy

#endif
