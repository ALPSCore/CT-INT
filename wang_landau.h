//
// Created by H. Shinaoka on 2015/12/04.
//

#ifndef IMPSOLVER_WANG_LANDAU_H
#define IMPSOLVER_WANG_LANDAU_H

#include <vector>
#include <stdexcept>
#include <cmath>
#include <algorithm>

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
      if (log_dos_renormalized_[idx] > 1E+3) {
        const double min_log_dos = *std::min(log_dos_renormalized_.begin(), log_dos_renormalized_.end());
        for (std::vector<double>::iterator it= log_dos_renormalized_.begin(); it != log_dos_renormalized_.end(); ++it)
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

#endif //IMPSOLVER_WANG_LANDAU_H
