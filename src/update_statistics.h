#ifndef ___HISTOGRAM___
#define ___HISTOGRAM___

#include <cmath>
#include <vector>
#include <complex>
#include <cassert>
#include <math.h>

#include <boost/random.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/discrete_distribution.hpp>

class update_proposer {
public:
    update_proposer() : mode_(0) {};

    update_proposer(int Nv, std::vector<double>& proposal_rate_Nv) :
            Nv_(Nv),mode_(0), proposal_rate_Nv_(proposal_rate_Nv), hist_valid_(Nv, 1.0), hist_invalid_(Nv, 1.0), dist01_(),
            dist_(proposal_rate_Nv.begin(), proposal_rate_Nv.end())
    {}

    void generated_valid_update(int num_vertices) {
        if (mode_!=0) {
            return;
        }
        hist_valid_[num_vertices-1] += 1.0;
    }

    void generated_invalid_update(int num_vertices) {
        if (mode_!=0) {
            return;
        }
        hist_invalid_[num_vertices-1] += 1.0;
    }


    void finish_learning(bool print) {
        if (mode_!=0) {
            throw std::logic_error("No more need for learning!");
        }
        mode_ = 1;

        if (print) {
            std::cout << " Summary of valid and invalid updates proposed in thermalization steps" << std::endl;
            std::cout << " Number of updates proposed (# of vertices in update, valid updates, invalid updates" << std::endl;
            for (int i=0; i<Nv_; ++i) {
                std::cout << "  " << i+1 << " : " << hist_valid_[i] <<  " " << hist_invalid_[i] << std::endl;
            }
        }

        for (int i=1; i<Nv_; ++i) {
            proposal_rate_Nv_[i] *= (proposal_rate_Nv_[i]/proposal_rate_Nv_[0])/(hist_valid_[i]/hist_valid_[0]);
        }
        double sum = std::accumulate(proposal_rate_Nv_.begin(), proposal_rate_Nv_.end(), 0.0);
        for (int i=0; i<Nv_; ++i) {
            proposal_rate_Nv_[i] /= sum;
        }
        if (print) {
            std::cout << " New proposal rates " << std::endl;
            for (int i=0; i<Nv_; ++i) {
                std::cout << "    " << i+1 << " , " << proposal_rate_Nv_[i] << std::endl;
            }
        }
        dist_ = boost::random::discrete_distribution<>(proposal_rate_Nv_.begin(), proposal_rate_Nv_.end());
    }

    int operator() (boost::random::mt19937& random_gen) {
        return gen_Nv(random_gen);
    }

    int gen_Nv(boost::random::mt19937& random_gen) {
        int Nv =  dist_(random_gen)+1;
        return Nv;
    }

private:
    int Nv_;
    int mode_;//0 is learning, 1 is learned
    std::vector<double> hist_valid_, hist_invalid_, proposal_rate_Nv_;
    boost::random::uniform_01<> dist01_;
    boost::random::discrete_distribution<> dist_;
};

class simple_update_statistcs {
public:
    simple_update_statistcs(int Nv) : counter_(), Nv_(Nv) {
        for (int iv=0; iv<Nv_; ++iv) {
            counter_.push_back(std::valarray<double>(0.0, 4));
        }
    }
    void not_valid_state(int iv) {
        assert(iv>=0 && iv<Nv_);
        ++counter_[iv][0];
    }
    void reducible(int iv) {
        assert(iv >= 0 && iv < Nv_);
        ++counter_[iv][1];
    }
    void accepted(int iv) {
        assert(iv>=0 && iv<Nv_);
        ++counter_[iv][2];
    }
    void rejected(int iv) {
        assert(iv>=0 && iv<Nv_);
        ++counter_[iv][3];
    }
    std::valarray<double> get_result(int iv) const {return counter_[iv];}
    void reset() {
        for (int iv=0; iv<Nv_; ++iv) {
            counter_[iv] = 0.0;
        }
    }
private:
    std::vector<std::valarray<double> > counter_;
    int Nv_;
};

//compute a spread of vertices in imaginary time
template<class V>
double compute_spread(typename V::const_iterator v_begin, typename V::const_iterator v_end, double beta) {
   double x=0,y =0;
   const double coeff = 2*M_PI/beta;
   //const int Nv = vertices.size();
   const int Nv = std::distance(v_begin, v_end);
   std::vector<double> angles(Nv);

   //distribute times on a unit circle and compute the math center.
   for (int iv=0; iv<Nv; ++iv) {
       //angles[iv] = vertices[iv].time()*coeff;
       angles[iv] = (v_begin+iv)->time()*coeff;
       x += cos(angles[iv]);
       y += sin(angles[iv]);
   }
   double norm = std::sqrt(x*x+y*y);
   double theta_c;
   if (norm==0) {
       theta_c = 0;
   } else {
       theta_c = atan2(y, x); //-pi<=theta_c<=pi
   }

   double max_dist_n = 0, max_dist_p = 0;
   for (int iv=0; iv<Nv; ++iv) {
       double angle = angles[iv]-theta_c; //-pi<=angle<=2*pi
       if (angle>=M_PI) angle -= 2*M_PI;
       assert(angle<=M_PI && angle>=-M_PI);
       if (angle>=0) {
           max_dist_p = std::max(max_dist_p, std::abs(angle));
       } else {
           max_dist_n = std::max(max_dist_n, std::abs(angle));
       }
   }

   return (max_dist_p+max_dist_n)*beta/(2*M_PI);
}

template<class T>
void rebin(std::valarray<T>& org_array, int nrebin) {
    const int new_len = org_array.size()/nrebin;
    std::valarray<T> new_array(new_len);
    int count=0;
    for (int i=0; i<new_len; ++i) {
        T sum = 0;
        for (int j=0; j<nrebin; ++j) {
            sum += org_array[count];
            ++count;
        }
        new_array[i] = sum;
    }
    org_array.resize(new_len);
    org_array = new_array;
}

template<class T>
void rebin(std::valarray<T>& org_array, double maxval, double maxval_new, int new_len) {
    std::valarray<T> new_array(0.0,new_len);
    const int old_len = org_array.size();

    assert(maxval_new<=maxval);

    for (int i=0; i<old_len; ++i) {
        const double cent = (i+0.5)*(maxval/old_len);
        const int pos = static_cast<int>(std::floor(new_len*cent/maxval_new));
        if (0 <= pos && pos < new_len) {
            new_array[pos] += org_array[i];
        }
        if (pos > new_len) {
            break;
        }
    }
    org_array.resize(new_len);
    org_array = new_array;
}

class scalar_histogram
{
    public:
        scalar_histogram() : num_bins_(0), max_val_(0.0), num_sample_(0), sumval(0.0,0), sumval2(0.0,0), counter(0.0,0) {};
        scalar_histogram(int num_bins_, double max_val) : num_bins_(num_bins_), num_sample_(0), max_val_(max_val), sumval(0.0,num_bins_), sumval2(0.0,num_bins_), counter(0.0,num_bins_){};

        void init(int num_bins_, double max_val) {
        	this->num_bins_ = num_bins_;
        	this->num_sample_ = 0;
        	this->max_val_ = max_val;
        	sumval.resize(num_bins_,0.0);
        	sumval2.resize(num_bins_,0.0);
        	counter.resize(num_bins_,0.0);
        	sumval = 0.0;
        	sumval2 = 0.0;
        	counter = 0.0;
        }

        bool add_sample(double distance, double value) {
            const int pos = static_cast<int>(std::floor(num_bins_*distance/max_val_));
            if (0 <= pos && pos < num_bins_) {
                ++num_sample_;
                sumval[pos] += value;
                sumval2[pos] += value*value;
                ++counter[pos];
                return true;
            } else {
                return false;
            }
        }

        std::valarray<double> get_mean() const {
            std::valarray<double> mean(sumval);
            for (int i=0; i<num_bins_; ++i) {
                mean[i] /= static_cast<double>(counter[i]);
            }
            return mean;
        }

        const std::valarray<double>& get_counter() const {
            return counter;
        }

        const std::valarray<double>& get_sumval() const {
            return sumval;
        }

        //maxdist is not updated if we do not have enough data.
        std::tuple<bool,double> update_cutoff(double cutoff_ratio, double maxdist, double mag=1.2) const {
            assert(cutoff_ratio>=0.0 && cutoff_ratio<=1.0);
            assert(mag>=1.0);
            const int min_count = 100;//for stabilization
            const int ndiv = 4;

            double maxdist_new = maxdist;

            std::valarray<double> counter_tmp = counter;
            std::valarray<double> sumval_tmp = sumval;

            rebin(counter_tmp, max_val_, maxdist, ndiv);
            rebin(sumval_tmp, max_val_, maxdist, ndiv);

            bool flag = false;
            double maxval = -1.0;
            for (int i=0; i<counter_tmp.size(); ++i) {
                flag = flag || (counter_tmp[i]<min_count);
                if (flag) {
                    break;
                }
                maxval = std::max(maxval, sumval_tmp[i]/counter_tmp[i]);
            }
            if (flag || maxval<0) {
            	//do not update maxdist_new
                return std::make_tuple(false,maxdist_new);
            }

            double ratio = (sumval_tmp[ndiv-1]/counter_tmp[ndiv-1])/maxval;
            if (ratio > cutoff_ratio) {
                maxdist_new *= mag;
            } else if (ratio < cutoff_ratio) {
                maxdist_new /= mag;
            }
            return std::make_tuple(true,maxdist_new);
        }

        void reset() {
            num_sample_ = 0;
            for (int i=0; i<num_bins_; ++i) {
                sumval[i] = 0.0;
                sumval2[i] = 0.0;
                counter[i] = 0;
            }
        }

        int num_bins() const {
        	return num_bins_;
        }

        int get_num_sample() const {
            return num_sample_;
        };

    private:
        int num_bins_, num_sample_;
        double max_val_;
        std::valarray<double> sumval, sumval2;
        std::valarray<double> counter;
};

class scalar_histogram_flavors
{
    public:
        scalar_histogram_flavors(int num_bins_, double max_val, int flavors) : flavors(flavors), max_val(max_val), histograms(flavors), num_bins(num_bins_) {
            //for (auto& elem : histograms) {
        	for (int i=0; i<histograms.size(); ++i) {
        		histograms[i].init(num_bins_, max_val);
        	}
        };

        bool add_sample(double distance, double value, int flavor) {
        	assert(flavor>=0 && flavor<flavors);
        	return histograms[flavor].add_sample(distance, value);
        }

        std::valarray<double> get_mean() const {
            std::valarray<double> mean_flavors(num_bins*flavors), mean(num_bins);

            for (int iflavor=0; iflavor<flavors; ++iflavor) {
            	mean = histograms[iflavor].get_mean();
            	for (int ibin=0; ibin<num_bins; ++ibin) {
            		assert(ibin+iflavor*num_bins<mean_flavors.size());
            		mean_flavors[ibin+iflavor*num_bins] = mean[ibin];
            	}
            }
            return mean_flavors;
        }

        std::valarray<double> get_counter() const {
            std::valarray<double> counter_flavors(num_bins*flavors);

            for (int iflavor=0; iflavor<flavors; ++iflavor) {
            	const std::valarray<double>& counter = histograms[iflavor].get_counter();
            	for (int ibin=0; ibin<num_bins; ++ibin) {
            		assert(ibin+iflavor*num_bins<counter_flavors.size());
            		counter_flavors[ibin+iflavor*num_bins] = counter[ibin];
            	}
            }
            return counter_flavors;
        }

        std::valarray<double> get_sumval() const {
            std::valarray<double> sumval_flavors(num_bins*flavors);

            for (int iflavor=0; iflavor<flavors; ++iflavor) {
            	const std::valarray<double>& sumval = histograms[iflavor].get_sumval();
            	for (int ibin=0; ibin<num_bins; ++ibin) {
                    assert(ibin+iflavor*num_bins<sumval_flavors.size());
                    sumval_flavors[ibin+iflavor*num_bins] = sumval[ibin];
            	}
            }
            return sumval_flavors;
        }

        double update_cutoff(double cutoff_ratio, double maxdist, double mag=1.2) const {
        	double maxdist_new = -1.0;
            //for (auto& elem : histograms) {
        	for (int i=0; i<histograms.size(); ++i) {
        		std::tuple<bool,double> r = histograms[i].update_cutoff(cutoff_ratio, maxdist, mag);
        		maxdist_new = std::max(maxdist_new, std::get<1>(r));
        	}
        	assert(maxdist_new>0);
        	return std::max(std::min(max_val, maxdist_new), max_val/num_bins);
        }

        void reset() {
        	//for (auto& elem : histograms) {
            for (int i=0; i<histograms.size(); ++i) {
        		histograms[i].reset();
        	}
        }

        int get_num_bins() const {
        	return num_bins;
        }

    private:
        int flavors, num_bins;
        double max_val;
        std::vector<scalar_histogram> histograms;
};

#endif
