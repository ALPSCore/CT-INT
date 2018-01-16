//
// Created by H. Shinaoka on 2015/12/04.
//

#ifndef IMPSOLVER_WANG_LANDAU_H
#define IMPSOLVER_WANG_LANDAU_H

#include <vector>
#include <stdexcept>
#include <cmath>
#include <algorithm>

#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/collectives.hpp>

class FlatHistogram {
public:
    FlatHistogram(double max_val, double min_val, unsigned int num_bin)
        : max_val_(max_val), min_val_(min_val), num_bin_(num_bin), lambda_(2.7), delta_((max_val_-min_val_)/num_bin_), log_dos_renormalized_(num_bin_, 0.0), log_dos_estimated_(num_bin_, 0.0), count_(0) {}

    double log_dos_renormalized(double value) const {
      if (value>max_val_ || value<min_val_)
        throw std::runtime_error("value is not within the range.");

      return log_dos_renormalized_[to_idx(value)];
    }

    double log_dos_estimated(double value) const {
      if (value>max_val_ || value<min_val_)
        throw std::runtime_error("value is not within the range.");

      return log_dos_estimated_[to_idx(value)];
    }

    double log_dos(double value) const {
      if (value>max_val_ || value<min_val_)
        throw std::runtime_error("value is not within the range.");

      return log_dos_renormalized_[to_idx(value)]+log_dos_estimated_[to_idx(value)];
    }

    void measure(double value) {
      if (value>max_val_ || value<min_val_)
        throw std::runtime_error("value is not within the range.");

      int idx = to_idx(value);
      log_dos_renormalized_[idx] += std::log(lambda_);
      ++count_;

      /*
      if (log_freq_[idx] > 1E+3) {
        const double min_log_dos = *std::min(log_freq_.begin(), log_freq_.end());
        for (std::vector<double>::iterator it= log_freq_.begin(); it != log_freq_.end(); ++it)
          *it -= min_log_dos;
      }
      */
    }

    unsigned int num_bin() const {return num_bin_;}
    unsigned int count() const {return count_;}

    void update_dos() {
      for (int i=0; i<num_bin_; ++i) {
        log_dos_estimated_[i] = log_dos_estimated_[i] + log_dos_renormalized_[i];
        log_dos_renormalized_[i] = 0;
      }
      count_ = 0;
//      lambda_ = std::sqrt(lambda_);
    }

    int to_idx(double value) const {
      if (value>max_val_ || value<min_val_) {
        return -1;
      } else {
        return range(static_cast<int>(value-min_val_)/delta_, num_bin_-1, 0);
      }
    }


private:
    const double max_val_, min_val_;
    const unsigned int num_bin_;
    const double delta_;
    double lambda_;
    unsigned long count_;
    std::vector<double> log_dos_renormalized_, log_dos_estimated_;


    int range(int i, int max, int min) const {
      return std::max(std::min(i, max), min);
    }
};

class FlatHistogramPertOrder {
public:
    typedef unsigned int uint;

    FlatHistogramPertOrder(unsigned int max_val, unsigned int min_val, unsigned long n_therm_steps)
        : max_val_(max_val), min_val_(min_val), num_bin_(max_val+1),
          log_lambda_(init_log_lambda_), log_f_(num_bin_, 0),
          counter_(num_bin_, 0), done_(false), top_index_(0), max_index_(0), has_guess_(true), n_therm_steps_(n_therm_steps) {
      assert(max_val>min_val);

      max_index_ = max_val;

    }

    double weight_ratio(uint value_new, uint value_old) const {
      if (value_new>max_index_)
        return 0;

      //if (value_new==120)
        //std::cout << "debug_ratio " << std::exp(log_weight(value_new)-log_weight(value_old)) << " value_old " << value_old << std::endl;

      return std::exp(log_weight(value_new)-log_weight(value_old));
    }

    double meas_weight(uint value) const {
      if (!done_)
        throw std::runtime_error("meas_measurement is called when done_==false.");

      if (log_f_[value]<std::log(1E-10)) {
        return 0.0;
      } else {
        return std::exp(log_f_[value]);
      }
    }

    void measure(uint value) {
      if (done_ || !has_guess_) return;

      check_range(value);
      if (value>max_index_)
        return;

      log_f_[value] += log_lambda_;
      ++counter_[value];
    }

    /*
    void sync_log_f(const boost::mpi::communicator& comm) {
      //check if the histogram is flat enough
      const double mean = std::accumulate(counter_.begin(), counter_.end(), 0.0)/(max_index_+1);
      const double min = *std::min_element(counter_.begin(), counter_.end());
      if (min<=criterion*mean) {
        std::cout << "Inogoring data..." << std::endl;
        std::fill(log_f_.begin(), log_f_.end(), 0.0);
      }

      std::vector<double> log_f_reduced(max_index_+1,0.0);
      boost::mpi::all_reduce(comm, &log_f_[0], max_index_+1, &log_f_reduced[0], std::plus<double>());
      const double max_log_f = *std::max_element(log_f_reduced.begin(), log_f_reduced.end());
      for (int i=0; i<max_index_+1; ++i)
        log_f_[i] = log_f_reduced[i]-max_log_f;
    }
    */

    boost::tuple<bool,double,double> flat_enough(const boost::mpi::communicator& comm) const {
      const double mean = std::accumulate(counter_.begin(), counter_.end(), 0.0)/(max_index_+1);
      const double min = *std::min_element(counter_.begin(), counter_.end());
      //std::cout << "debug min_pos " <<  std::distance(counter_.begin(),std::min_element(counter_.begin(), counter_.end())) << " node " << comm.rank() << std::endl;
      //double data = (double) (min>criterion*mean);
      return boost::make_tuple(min/mean>criterion, min, mean);
      /*
      double data_rsv;
      std::cout << " size " << counter_.size() << std::endl;
      std::cout << "node = " << comm.rank() << " min/mean " << min/mean << " min " << min << " mean " << mean << std::endl;
      if (comm.rank()==0) {
        for (int i=0; i<max_index_+1; ++i)
          std::cout << "debug " << i << "  " << counter_[i]/mean << " " << counter_[i] << std::endl;
      }
      boost::mpi::all_reduce(comm, &data, 1, &data_rsv, std::plus<double>());
      return std::pair<bool,double>(data/comm.size()>criterion_walker, data/comm.size());
       */
    }

    void update_dos(const boost::mpi::communicator& comm, bool verbose=true) {
      if (done_ || !has_guess_) return;


      //lambda_ = std::sqrt(lambda_);
      log_lambda_ = std::max(0.5*log_lambda_, min_log_lambda_);
      if (verbose) {
        std::cout << " new lambda = " << std::exp(log_lambda_) << std::endl;
        std::cout << " new log_lambda = " << log_lambda_ << std::endl;
      }
      rescale_f();

      //std::vector<double> counter_reduced(max_index_+1,0.0);
      //boost::mpi::all_reduce(comm, &counter_[0], max_index_+1, &counter_reduced[0], std::plus<double>());
      //sync_log_f(comm);

      if (verbose) {
        std::cout << " mean = " << std::accumulate(counter_.begin(), counter_.end(), 0.0)/(max_index_+1) << std::endl;
        std::cout << " max = " << *std::max_element(counter_.begin(), counter_.end())<< std::endl;
        std::cout << " min = " << *std::min_element(counter_.begin(), counter_.end())<< std::endl;
        for (int i=0; i<max_index_+1; ++i)
          std::cout << " counter  " << i << " " << counter_[i] << " " << log_f_[i] << std::endl;
      }

      std::fill(counter_.begin(),counter_.end(),0);
    }

    bool is_learning_done() {
      return done_;
    }

    void finish_learning(const boost::mpi::communicator& comm, bool verbose) {
      done_ = true;

      std::vector<double> log_f_reduced(max_index_+1,0.0);
      boost::mpi::all_reduce(comm, &log_f_[0], max_index_+1, &log_f_reduced[0], std::plus<double>());
      for (int i=0; i<max_index_+1; ++i)
        log_f_[i] = log_f_reduced[i];
      rescale_f();
      std::fill(counter_.begin(),counter_.end(),0);
      if (verbose) {
        std::cout << " log_f for measurement steps " << std::endl;
        for (int i=0; i<max_index_+1; ++i)
          std::cout << " pert_order " << i << " " << log_f_[i] << std::endl;
      }
    }

private:
    const uint max_val_, min_val_, num_bin_;
    const double criterion = 0.8;
    const double init_log_lambda_ = 0.993251773010283; //std::log(2.7);
    const double min_log_lambda_ = 0.000999500333083; //std::log(1.001);
    const unsigned long n_therm_steps_;

    //const double criterion_walker = 0.8;
    double log_lambda_;
    std::vector<double> log_f_, log_dos_guess_;
    std::vector<double> counter_;
    uint top_index_, max_index_;
    bool done_, has_guess_;

    void check_range(uint value) const {
      if (value>max_val_ || value<min_val_)
        throw std::runtime_error("value is not within the range.");
    }

    void rescale_f() {
      const double max_log_f = *std::max_element(log_f_.begin(), log_f_.end());
      for (int i=0; i<max_index_+1; ++i)
        log_f_[i] -= max_log_f;
    }

    double log_weight(uint value) const {
      if (!has_guess_)
        return 0.;

      if (value<=max_index_) {
        return -log_f_[value];
      } else {
        throw std::runtime_error("Not within range!");
      }
    }
};

#endif //IMPSOLVER_WANG_LANDAU_H
