#pragma once

#include "types.h"

namespace alps {
    namespace ctint {
        template<class T>
        class operator_time_tmpl
        {
        public:
            operator_time_tmpl() : time_(0), small_idx_(0) {}
            operator_time_tmpl(T time, int small_idx) : time_(time), small_idx_(small_idx) {}

            inline T time() const {return time_;}
            inline int small_index() const {return small_idx_;}

            inline void set_time(T time) { time_ = time;}
            inline void set_small_index(int small_idx) { small_idx_ = small_idx;}

        private:
            T time_;
            int small_idx_;
        };

        typedef operator_time_tmpl<itime_t> operator_time;

        template<class T>
        bool operator<(const operator_time_tmpl<T>& t1, const operator_time_tmpl<T>& t2) {
          if (t1.time()==t2.time()) {
            return (t1.small_index()<t2.small_index());
          } else {
            return (t1.time()<t2.time());
          }
        }

        template<class T>
        bool operator==(const operator_time_tmpl<T>& t1, const operator_time_tmpl<T>& t2) {
          return (t1.time()==t2.time()) && (t1.small_index()==t2.small_index());
        }

/*creation and annihilation operator class*/
        class c_or_cdagger   //represents a creation operator or an annihilation operator
        {
        public:
            c_or_cdagger( const spin_t z,const site_t s, const operator_time t)
            {
              s_ = s;
              z_ = z;
              t_ = t;
            }
            virtual ~c_or_cdagger() {}

            inline const spin_t &flavor() const {return z_;}
            inline spin_t &flavor() {return z_;}
            inline const operator_time &t() const{return t_;}
            inline operator_time &t() {return t_;}
            inline const site_t &s() const {return s_;}
            inline void flavor(spin_t z){z_=z;}
            inline void s(site_t s){s_=s;}


        private:
            site_t s_;      //this vertex's site
            operator_time t_;     //its imaginary time point
            spin_t z_;      //its flavor or flavor

        public:
            void set_time(operator_time time) {
              t_ = time;
            }
        };// creator, annihilator;

        class creator : public c_or_cdagger {
        public:
            using c_or_cdagger::operator=;
            creator(const spin_t z,const site_t s, const operator_time t)
              : c_or_cdagger(z, s, t){};
            creator()
              : c_or_cdagger(-1, -1, operator_time(-1.0,0)){};
            ~creator(){};
        };

        class annihilator : public c_or_cdagger {
        public:
            using c_or_cdagger::operator=;
            annihilator(const spin_t z,const site_t s, const operator_time t)
              : c_or_cdagger(z, s, t){};
            annihilator()
              : c_or_cdagger(-1, -1, operator_time(-1.0,0)){};
            ~annihilator(){};
        };


        inline
        bool operator==(const creator& op1, const creator& op2) {
          return (op1.s()==op2.s()) && (op1.t()==op2.t()) && (op1.flavor()==op2.flavor());
        }

        inline
        bool operator==(const annihilator& op1, const annihilator& op2) {
          return (op1.s()==op2.s()) && (op1.t()==op2.t()) && (op1.flavor()==op2.flavor());
        }

    }
}

