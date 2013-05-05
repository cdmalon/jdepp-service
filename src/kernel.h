// pecco -- please enjoy classification with conjunctive features
//  $Id: kernel.h 827 2012-05-10 14:12:37Z ynaga $
// Copyright (c) 2008-2012 Naoki Yoshinaga <ynaga@tkl.iis.u-tokyo.ac.jp>
#ifndef POLYK_CLASSIFY_H
#define POLYK_CLASSIFY_H

#include <sys/stat.h>
#include <cassert>
#include <cmath>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <set>
#include "typedef.h"
#ifdef USE_HASH
#include <tr1/unordered_map>
#endif
#include "timer.h"
#include "classify.h"

namespace pecco {

  class kernel_model : public ClassifierBase <kernel_model>
  {
    // type alias
    typedef std::vector <ny::uint> ss_t;
    typedef std::vector <ny::uint>::const_iterator ss_it;
  private:
    // kernel parameters; (_s w * x + _r)^_d
    double                 _s;
    double                 _r;
    double *               _b;
    double *               _m0;
    // support vector
    std::vector <ny::fv_t> _sv;   // SV ID -> <alpha, feature vector>
    std::vector <ss_t>     _f2ss; // inverted indices from feature to SV IDs
    std::vector <ny::fl_t> _alph; // weight of SVs (used in PKI)
    // kernel cache
    double *               _polyk;
    double *               _spolyk;
    // temporary variables
    ny::uint *             _dot;  // store feature on/off
    std::vector <bool>     _fbit;
    // polynomial kernel expansion
    const double           _sigma;
    const double           _fratio;
    std::vector <double>   _sigma_pos;
    std::vector <double>   _sigma_neg;
    ny::uint               _processed;
    double                 _coeff[MAX_KERNEL_DEGREE + 1];
    double                 _max_coeff;
    ny::uint               _exnunit;     // # conjunctive features
    ny::uint               _exnunit_cut; // # conjunctive features (pruned)
    // timer
#ifdef USE_TIMER
    ny::Timer *            _pki_t;    // pki
#endif
    // internal functions
    bool _splitInit ();
    bool _precomputeKernel ();
    bool _setFtrie ();
    bool _setPKEcoeff ();
    void _pkePrefixSpan (ny::fv_t &fc, std::vector <ny::fl_t> &fw, const std::vector <std::pair <ny::uint, int> > &proj, std::vector <FeatKey*> &pke_key);
    // classifier
    template <binary_t FLAG>
    void _pkiClassify   (double * score, const ny::fv_it &beg, const ny::fv_it &end) const;
    void _pkiClassify   (const ny::fv_t &fv, double * score) const;
    template <binary_t FLAG>
    void _splitClassify (double * score, ny::fv_it it, const ny::fv_it &beg, const ny::fv_it &end);
    void _splitClassify (const ny::fv_t &fv, double * score, ny::fv_it it, const ny::fv_it &end);
    void _splitClassify (const ny::fv_t &fv, double * score)
    { return _splitClassify (fv, score, fv.begin (), fv.end ()); }
    void _setup_binary_labels ();
  public:
    kernel_model (const pecco::option &opt) : ClassifierBase <kernel_model> (opt), _s (0), _r (0), _b (0), _m0 (0), _sv (), _f2ss (), _alph (), _polyk (0), _spolyk (0), _dot (0), _fbit (), _sigma (ny::strton <double> (_opt.sigma)), _fratio (ny::strton <double> (_opt.fratio)), _sigma_pos (0), _sigma_neg (0), _processed (0), _max_coeff (0), _exnunit (0), _exnunit_cut (0)
#ifdef USE_TIMER
      , _pki_t (_timer_pool.push ("pki", "classify"))
#endif
    {} //  _nl = 1; _li2l.push_back ("+1"); _li2l.push_back ("-1");
    ~kernel_model () {
      delete [] _b;
      delete [] _m0;
      if (_opt.algo != PKE) {
        delete [] _polyk;
        if (_fratio) delete [] _spolyk;
      }
      if (_opt.algo == PKI) delete [] _dot;
#ifndef ABUSE_TRIE
      if (_opt.algo == PKE || _opt.algo == FST) delete [] _fw;
      if (_opt.algo == FST) delete [] _fsw;
#endif
      for (size_t li = 0; li < _li2l.size (); ++li)
        delete [] _li2l[li];
    }
    template <binary_t FLAG>
    void addScore (double * score, const ny::uint pos) const;
    template <binary_t FLAG>
    void addScore (double * score, const int n, const ny::fl_t * const w) const;
    template <binary_t FLAG>
    void addScore (double * score, const ny::uint pos, const double m) const;
    void baseClassify (const ny::fv_t &fv, double * score, ny::fv_it it, const ny::fv_it &end)
    { _splitClassify (fv, score, it, end); }
    bool load (const char * model); // set up model
    void printParam () {
      std::fprintf (stderr, "kernel: (%g * w^T x + %g)^%d\n", _s, _r, _d);
      std::fprintf (stderr, "# support vectors: %d\n", _nunit);
      std::fprintf (stderr, "# active features: %d",   _nf_cut);
      if (_opt.algo == PKE || _opt.algo == FST)
        std::fprintf (stderr, " (%d)", _exnunit_cut);
      std::fprintf (stderr, "\n");
      if (_opt.verbose > 1)
        std::fprintf (stderr, "  # common features: %d\n", _f_r - 1);
    }
    // classification interface
    void classify (ny::fv_t &fv, double * score) {
      for (ny::uint li = 0; li <_nl; ++li)
        score[li] = -_b[li];
      TIMER (if (_opt.verbose > 0) _enc_t->startTimer ());
      _convertFv2Fv (fv);
      TIMER (if (_opt.verbose > 0) _enc_t->stopTimer ());
      switch (_opt.algo) {
        case PKE:
          TIMER (_pke_t->startTimer ());
          if (_d >= 2) _sortFv (fv);
          for (ny::uint li = 0; li <_nl; ++li) score[li] +=_m0[li];
          if (_sigma) _pkeClassify (fv, score); else _splitClassify (fv, score);
          TIMER (_pke_t->stopTimer ());
          break;
        case FST:
          TIMER (_fst_t->startTimer ());
          _sortFv (fv);
          for (ny::uint li = 0; li <_nl; ++li) score[li] +=_m0[li];
          _fstClassify (fv, score);
          TIMER (_fst_t->stopTimer ());
          break;
        case PKI:
          TIMER (_pki_t->startTimer ());
          _pkiClassify (fv, score);
          TIMER (_pki_t->stopTimer ());
          break;
        default: ny::print_err (HERE "unknown classifier.\n");
      }
    }
    bool is_binary_classification () const { return _nl == 1; }
#ifdef ABUSE_TRIE
    bool abuse_trie               () const { return true; }
#else
    bool abuse_trie               () const { return false; }
#endif
    bool binClassify    (ny::fv_t &fnv, char * target = 0) {
      if (is_binary_classification ()) {
        double score = 0; classify (fnv, &score); return score > 0;
      } else {
        static std::vector <double> score;
        if (score.empty ()) score.resize (_nl); // initialized later
        classify (fnv, &score[0]);
        const ny::uint li =
          static_cast <ny::uint> (std::max_element (&score[0], &score[_nl]) - &score[0]);
        if (target)
          return std::strcmp (_li2l[li], target) == 0;
        else
          return li == _tli;
      }
    }
    double getProbability (ny::fv_t &fv) {
      if (! is_binary_classification ())
        ny::print_err (HERE "sorry, probability output is unsupported.\n");
      double score = 0;
      classify (fv, &score);
      return sigmoid (score);
    }
    ny::uint getLabel (const double * score) {
      return is_binary_classification () ?
        (*score >= 0 ? 0 : 1) :
        static_cast <ny::uint> (std::max_element (&score[0], &score[_nl]) - &score[0]);
    }
    void printScore (const ny::uint li, const double * score) {
      std::fprintf (stdout, "%s %f\n", _li2l[li],
                    is_binary_classification () ? *score : score[li]);
    }
    static double sigmoid (double x)
    { return 1.0 / (1.0 + std::exp (-x)); }
  };
  template <>
  inline void kernel_model::addScore <BINARY> (double * score, const ny::uint pos) const
  { *score += _fw[pos]; }
  template <>
  inline void kernel_model::addScore <MULTI>  (double * score, const ny::uint pos) const {
    const ny::fl_t * const fw = &_fw[pos * _nl];
    for (ny::uint li (0); li < _nl; ++li)
      score[li] += fw[li];
  }
  template <>
#ifdef ABUSE_TRIE
  inline void kernel_model::addScore <BINARY> (double * score, const int n, const ny::fl_t * const) const
  { union byte_4 b4 (n); b4.u <<= 1; *score += b4.f; }
#else
  inline void kernel_model::addScore <BINARY> (double * score, const int n, const ny::fl_t * const w) const
  { *score += w[n]; }
#endif
  template <>
  inline void kernel_model::addScore <MULTI>  (double * score, const int n, const ny::fl_t * const w) const {
    const ny::fl_t * const fw = &w[n];
    for (ny::uint li = 0; li < _nl; ++li)
      score[li] += fw[li];
  }
  template <>
  inline void kernel_model::addScore <BINARY> (double * score, const ny::uint pos, const double m) const
  { *score += _alph[pos] * m; }
  template <>
  inline void kernel_model::addScore <MULTI>  (double * score, const ny::uint pos, const double m) const {
     const ny::fl_t * const alph = &_alph[pos * _nl];
    for (ny::uint li (0); li < _nl; ++li)
      score[li] += alph[li] * m;
  }
}
#endif /* POLYK_CLASSIFY_H */
