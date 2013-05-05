// pecco -- please enjoy classification with conjunctive features
//  $Id: linear.h 832 2012-05-11 02:09:28Z ynaga $
// Copyright (c) 2008-2012 Naoki Yoshinaga <ynaga@tkl.iis.u-tokyo.ac.jp>
#ifndef LLM_CLASSIFY_H
#define LLM_CLASSIFY_H

#include <sys/stat.h>
#include <cmath>
#include <cassert>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <set>
#include "typedef.h"
#include "timer.h"
#include "classify.h"

// switch off
#undef ABUSE_TRIE

namespace pecco {
  class linear_model : public ClassifierBase <linear_model> {
  private:
    // temporary variables
    std::vector <ny::fv_t>                _fcv;
    std::vector <std::vector <ny::fl_t> > _fw_tmp;
    // internal functions
    void _convertCfstr2Cf (char * &p, ny::fv_t &fv, const bool flag);
    bool _setFtrie ();
    // classification
    double _calcProb (const double * const score, ny::uint target_id) const {
      double sum (0), prob_pos (0);
      for (size_t li = 0; li < _nl; ++li) {
        double prob = std::exp (score[li]);
        sum += prob;
        if (li == target_id) prob_pos = prob;
      }
      return (prob_pos / sum);
    }
  public:
    linear_model (const pecco::option &opt) : ClassifierBase <linear_model> (opt), _fcv (), _fw_tmp () {}
    ~linear_model () {
      for (size_t i = 0; i < _nl; ++i)  delete [] _li2l[i];
      delete [] _fw;
      if (_opt.algo == FST) delete [] _fsw;
    }
    template <bool BINARY>
    void addScore (double * score, const ny::uint pos) const;
    template <bool BINARY>
    void addScore (double * score, const int n, const ny::fl_t * const w) const;
    void baseClassify (const ny::fv_t &fv, double * score, ny::fv_it it, const ny::fv_it &end)
    { _pkeClassify (fv, score, it, end); }
    bool load  (const char * model); // set up model
    void printParam () {
      std::fprintf (stderr, "maximum order of feature combination: %d\n", _d);
      std::fprintf (stderr, "# of active features: %d (%d)\n", _nf, _nunit);
    }
    // classification interface
    void classify (ny::fv_t &fnv, double * score) {
      std::fill (&score[0], &score[_nl], 0);
      TIMER (if (_opt.verbose > 0) _enc_t->startTimer ());
      _convertFv2Fv (fnv);
      TIMER (if (_opt.verbose > 0) _enc_t->stopTimer ());
      switch (_opt.algo) {
        case PKE:
          TIMER (_pke_t->startTimer ());
          if (_d >= 2) _sortFv (fnv); _pkeClassify (fnv, &score[0]);
          TIMER (_pke_t->stopTimer ());
          break;
        case FST:
          TIMER (_fst_t->startTimer ());
          _sortFv (fnv); _fstClassify (fnv, &score[0]);
          TIMER (_fst_t->stopTimer ());
          break;
        default: ny::print_err (HERE "unknown classifier.\n");
      }
    }
    bool is_binary_classification () const { return _nl == 2; }
    bool abuse_trie               () const { return false; }
    bool binClassify (ny::fv_t &fnv, char * target = 0) {
      static std::vector <double> score;
      if (score.empty ()) score.resize (_nl); // initialized later
      classify (fnv, &score[0]);
      const ny::uint li
        = static_cast <ny::uint> (std::max_element (&score[0], &score[_nl]) - &score[0]);
      if (target)
        return std::strcmp (_li2l[li], target) == 0;
      else
        return li == _tli;
    }
    double getProbability (ny::fv_t &fnv, char * target = 0) {
      static std::vector <double> score;
      if (score.empty ()) score.resize (_nl); // initilized later
      ny::uint target_id = _tli;
      if (target) {
        lmap::const_iterator it = _l2li.find (target);
        if (it == _l2li.end ())
          ny::print_err ("unknown label: %s", target);
        target_id = it->second;
      }
      classify (fnv, &score[0]);
      return _calcProb (&score[0], target_id);
    }
    ny::uint getLabel (const double * score) {
      return static_cast <ny::uint> (std::max_element (&score[0], &score[_nl]) - &score[0]);
    }
    void printScore (const ny::uint li, const double * score)
    { std::fprintf (stdout, "%s %f\n", _li2l[li], _calcProb (score, li)); }
  };
  template <>
  inline void linear_model::addScore <BINARY> (double * score, const ny::uint pos) const
  { score[0] += _fw[pos * 2], score[1] += _fw[pos * 2 + 1]; }
  template <>
  inline void linear_model::addScore <MULTI>  (double * score, const ny::uint pos) const {
    const ny::fl_t * const fw = &_fw[static_cast <ny::uint> (pos) * _nl];
    for (ny::uint li = 0; li < _nl; ++li) score[li] += fw[li];
  }
  template <>
  inline void linear_model::addScore <BINARY> (double * score, const int n, const ny::fl_t * const w) const
  { score[0] += w[n], score[1] += w[n + 1]; }
  template <>
  inline void linear_model::addScore <MULTI>  (double * score, const int n, const ny::fl_t * const w) const {
    const ny::fl_t * const fw = &w[static_cast <ny::uint> (n)];
    for (ny::uint li = 0; li < _nl; ++li) score[li] += fw[li];
  }
}

#endif /* LLM_CLASSIFY_H */
