// opal -- online learning with kernel slicing
//  $Id: pa.h 826 2012-05-10 14:07:00Z ynaga $
// Copyright (c) 2009-2012 Naoki Yoshinaga <ynaga@tkl.iis.u-tokyo.ac.jp>
#ifndef OPAL_PA_H
#define OPAL_PA_H

#include <getopt.h>
#include <stdint.h>
#include <cstdio>
#include <cstdlib>
#include <cstdarg>
#include <cstring>
#include <climits>
#include <cmath>
#include <valarray>
#include <vector>
#include <map>
#include <algorithm>
#include <numeric>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef USE_MT19937
#include <tr1/random>
#endif
#if defined (USE_HASH_TRIE) && ! defined (USE_HASH)
#define USE_HASH
#endif
#ifdef USE_HASH
#include <tr1/unordered_map>
#endif
#if ! defined (USE_HASH_TRIE) && ! defined (USE_MAP_TRIE)
#include "cedar.h"
#endif
#include "timer.h"

#define OPAL_COPYRIGHT  "opal - online learning with kernel slicing\n\
Copyright (c) 2009-2012 Naoki Yoshinaga, All rights reserved.\n\
\n\
Usage: %s [options] train model test\n\
\n\
train   training file       set '-' to skip training\n\
model   model file          set '-' to training/test w/o saving a mode\n\
test    test file           set '-' to skip testing\n\
\n"

#define OPAL_OPT "Optional parameters in training and testing:\n\
  -t, --kernel=TYPE         select type of kernel function\n\
                            * 0 - linear     (w^T * x)\n\
                              1 - polynomial (s^T * x + 1)^d\n\
  -d, --kernel-degree=INT   parameter d in polynomial kernel (0)\n\
  -N, --kernel-splitN=INT   parameter N in kernel splitting (0)\n\
                            (a learner [-N 0] will automatically calibrates N\n\
                             to keep the weight trie within [-T] MiB bytes)\n\
  -T, --max-trie-size=INT   max. RAM size (MiB) for feature weight trie (32)\n\
                            (ignored if kernel-splitN [-N] is specified)\n\
\n"

#define OPAL_OPT_TRAIN0 "Optional parameters in training:\n\
  -m, --model0=FILE         re-train a model trained w/ opal ('-')\n\
  -l, --learner=TYPE        select learning algorithm\n\
                              0 - Perceptron\n\
                              1 - Passive Aggressive    (PA)\n\
                            * 2 - Passive Aggressive I  (PA-I)\n\
                              3 - Passive Aggressive II (PA-II)\n\
  -c, --reg-cost=FLOAT      PA-I/-II aggressiveness parameter C (1.0)\n\
  -i, --iteration=INT       # iterations (10)\n\
  -a, --averaging           average parameters\n\
  -s, --shuffling           shuffle training examples on RAM\n\
  -b, --buffer=TYPE         select type of buffer to read examples\n\
                            * 0 - Use RAM to load examples\n\
                              1 - Use DISK to cache examples\n\
                              2 - Do not cache examples\n\
  -M, --max-examples=INT    max. # examples used in training (0: all)\n"

#ifdef USE_POLYK
#define OPAL_OPT_TRAIN1 "\n"
#else
#define OPAL_OPT_TRAIN1 "  -k, --kernel-slicing      perform kernel slicing\n\
  -p, --pruning-margin      terminate margin computation if unnecessarily\n"
#endif

#ifdef USE_MULTICLASS
#define OPAL_OPT_TRAIN2 "  -C, --max-classes=INT     max. # classes (needed for [-b 2])\n\
\n"
#else
#define OPAL_OPT_TRAIN2 "\n"
#endif

#define OPAL_OPT_TEST "Optional parameters in testing:\n\
  -O, --output=TYPE        select output type of testing\n\
                            * 0 - report accuracy\n\
                              1 - report accuracy per iteration\n\
                              2 - output margins\n\
\n"

#define OPAL_OPT_MISC "Misc.:\n\
  -h, --help               show this help and exit\n"

static const  char * opal_short_options = "t:d:N:T:m:l:c:i:kpasb:M:C:O:h";
static struct option opal_long_options[] = {
  {"kernel",         required_argument, NULL, 't'},
  {"kernel-degree",  required_argument, NULL, 'd'},
  {"kernel-splitN",  required_argument, NULL, 'N'},
  {"max-trie-size",  required_argument, NULL, 'T'},
  {"model0",         required_argument, NULL, 'm'},
  {"learner",        required_argument, NULL, 'l'},
  {"reg-cost",       required_argument, NULL, 'c'},
  {"iteration",      required_argument, NULL, 'i'},
  {"kernel-slicing", required_argument, NULL, 'k'},
  {"pruning-margin", required_argument, NULL, 'p'},
  {"iteration",      required_argument, NULL, 'i'},
  {"averaging",      no_argument,       NULL, 'a'},
  {"shuffing",       no_argument,       NULL, 's'},
  {"buffer",         required_argument, NULL, 'b'},
  {"max-examples",   required_argument, NULL, 'M'},
  {"max-classes",    required_argument, NULL, 'C'},
  {"output",         required_argument, NULL, 'O'},
  {"help",           no_argument,       NULL, 'h'},
  {NULL, 0, NULL, 0}
};

extern char * optarg;
extern int    optind;

namespace opal {
  // global type alias
#ifdef USE_FLOAT
  typedef float   fl_t;
#else
  typedef double  fl_t;
#endif
  typedef unsigned int   uint;
  typedef unsigned char  uchar;
  typedef std::vector <uint>            fv_t;
  typedef fv_t::const_iterator          fv_it;
  typedef fv_t::const_reverse_iterator  fv_rit;
  typedef long  label_t;
#ifdef USE_MULTICLASS
  typedef std::valarray <fl_t>  w_t;
#else
  typedef fl_t                  w_t;
#endif
  // static variables and functions
  static const size_t KEY_SIZE = 8;
  static const size_t BUF_SIZE = 1 << 18;
  static const int    MAX_KERNEL_DEGREE = 4;
  // \sum_{i=0}^{k} nCk (* num_class) features in an array-based pseudo trie
  static const uint   PSEUDO_TRIE_N[MAX_KERNEL_DEGREE] = {0, 21, 11, 8};
  // coefficient for kernel expansion; support s=1 and r=1 case
  // refer a general case to pecco/kernel.cc
  static const fl_t   COEFF[MAX_KERNEL_DEGREE][MAX_KERNEL_DEGREE] =
    {{0,  0,  0,  0},  // 0
     {1,  1,  0,  0},  // 1
     {1,  3,  2,  0},  // 2
     {1,  7, 12,  6}}; // 3
  static inline bool getLine (FILE* &fp, char* &line, size_t &read) {
#ifdef __APPLE__
    if ((line = fgetln (fp, &read)) == NULL) return false;
#else
    static ssize_t read_ = 0; static size_t size = 0; // static helps inlining
    if ((read_ = getline (&line, &size, fp)) == -1) return false;
    read = read_;
#endif
    *(line + read - 1) = '\0';
    return true;
  }
  class byte_encoder {
  private:
    uchar _key[KEY_SIZE];
    uint  _len;
  public:
    byte_encoder () : _key (), _len (0) {}
    byte_encoder (uint i) : _key (), _len (0) { encode (i); }
    uint encode  (uint i, uchar * const key) const {
      uint len = 0;
      for (key[len] = (i & 0x7f); i >>= 7; key[++len] = (i & 0x7f))
        key[len] |= 0x80;
      return ++len;
    }
    void encode (uint i) { _len = encode (i, _key); }
    const char * key  () { return reinterpret_cast <const char *> (&_key[0]); }
    uint         len  () { return _len; }
  };
  static inline void print_err (char const * format, ...) {
    std::va_list argptr;
    va_start (argptr, format);
    vfprintf (stderr, format, argptr);
    va_end   (argptr);
    std::exit (1);
  }
  template <typename T> T strton (const char * s, char ** error)
  { return static_cast <T> (std::strtol (s, error, 10)); }
#ifdef USE_HASH
  // incremental Folwer / Noll / Vo (FNV) Hash (FNV-1)
  static const uint FNV_PRIME = 0x01000193;
  static const uint FNV_BASIS = 0x811c9dc5;
  struct inc_fnv { // http://isthe.com/chongo/tech/comp/fnv/#FNV-param
    size_t operator() (const uint ret, const uint fi) const
    { return (ret * FNV_PRIME) ^ fi; }
  };
#endif
  // map- or unordered_map-based trie implementation
#if defined (USE_HASH_TRIE) || defined (USE_MAP_TRIE)
  typedef uint64_t edge; // 32bit parend id + 32bit integer label
  struct node {
    uint id;
    int  value;
    node (const uint id_) : id (id_), value (0) {};
  };
#endif
  // options
  enum kernel_t  { LINEAR, POLY };
  enum algo_t    { P, PA, PA1, PA2 };
  enum buffer_t  { RAM, DISK, null };
  static const char * algo[] = { "P", "PA", "PA1", "PA2" };
  struct option { // option handler
    enum mode_t { BOTH, TRAIN, TEST };
    const char *com, *train, *model0, *model, *test;
    //
    mutable kernel_t  kernel;   // kernel-type
    mutable uint      d;        // kernel-degree
    mutable uint      N;        // kernel-splitN
    algo_t            algo;     // algorithm
    fl_t              C;        // reg-cost
    uint              iter;     // # iteration
    bool              slicing;
    bool              pruning;
    bool              average;
    bool              shuffle;
    buffer_t          buffer;   // buffer type
    size_t            M;        // max-examples
    size_t            T;        // max-trie-size
    mutable uint      nclass;   // max-classes
    uint              output;
    mode_t            mode;
    bool              shrink;
    option () : com ("--"), train ("-"), model0 ("-"), model ("-"), test ("-"), kernel (LINEAR), d (0), N (0), algo (PA1), C (1.0), iter (10), slicing (false), pruning (false), average (false), shuffle (false), buffer (RAM), M (0), T (32 << 20), nclass (1), output (0), mode (BOTH), shrink (false) {}
    option (int argc, char ** argv) : com (argc ? argv[0] : "--"), train ("-"), model0 ("-"), model ("-"), test ("-"), kernel (LINEAR), d (0), N (0), algo (PA1), C (1.0), iter (10), slicing (false), pruning (false), average (false), shuffle (false), buffer (RAM), M (0), T (32 << 20), nclass (1), output (0), mode (BOTH), shrink (false) {
      set (argc, argv);
    }
    void set (int argc, char ** argv) { // getOpt
      if (argc == 0) return;
      optind = 1;
      while (1) {
        int opt = getopt_long (argc, argv,
                               opal_short_options, opal_long_options, NULL);
        if (opt == -1) break;
        char * err = NULL;
        switch (opt) {
          case 't': kernel  = strton <kernel_t> (optarg, &err); break;
          case 'd': d       = strton <uint> (optarg, &err); break;
          case 'N': N       = strton <uint> (optarg, &err); break;
            // training params
          case 'm': model0  = optarg; break;
          case 'l': algo    = strton <algo_t> (optarg, &err); break;
#ifdef USE_FLOAT
          case 'c': C       = std::strtof (optarg, &err); break;
#else
          case 'c': C       = std::strtod (optarg, &err); break;
#endif
          case 'i': iter    = strton <uint> (optarg, &err);   break;
          case 'a': average = true; break;
          case 's': shuffle = true; break;
          case 'b': buffer  = strton <buffer_t> (optarg, &err);     break;
          case 'M': M       = strton <size_t> (optarg, &err);       break;
          case 'T': T       = strton <size_t> (optarg, &err) << 20; break;
#ifndef USE_POLYK
          case 'k': slicing = true; break;
          case 'p': pruning = true; break;
#endif
#ifdef USE_MULTICLASS
          case 'C': nclass  = strton <uint> (optarg, &err); break;
#endif
            // testing params
          case 'O': output  = strton <uint> (optarg, &err); break;
            // misc
          case 'h': printCredit (); printHelp (); std::exit (0); break;
          default:  printCredit (); std::exit (0);
        }
        if (err && *err)
          print_err ("unrecognized option value: %s\n", optarg);
      }
      // errors & warnings
      if (kernel != POLY && kernel != LINEAR)
        print_err ("unknown kernel fucntion [-t].\n");
      if (algo != P && algo != PA && algo != PA1 && algo != PA2)
        print_err ("unknown learning algorithm [-l].\n");
      if (buffer != RAM && buffer != DISK && buffer != null)
        print_err ("unknown buffering method [-b].\n");
#ifdef USE_MULTICLASS
      if (buffer == null && nclass == 1)
        print_err ("set # classes [-C] when you don't buffer training data.\n");
#else
      if (nclass != 1)
        std::fprintf (stderr, "# classes [-C] is ignored in training a binary classifier.\n");
      nclass = 1;
#endif
      if (iter == 0) print_err ("# iterations  [-i] must be >= 1.\n");
      if (C != 1.0 && (algo == P || algo == PA))
        std::fprintf (stderr, "reg-cost C [-c] is ignored in P [-l 0] and PA [-l 1].\n");
      if (kernel == LINEAR) {
        if (d != 0)
          std::fprintf (stderr, "kernel-degree [-d] is ignored in linear kernel [-t 0].\n");
        d = 0;
      } else {
        if (d == 0 || d >= 4) print_err ("set kernel_degree [-d] to 1-3.\n");
        if (d == 1 && (slicing || pruning)) {
          std::fprintf (stderr, "kernel slicing [-k] (or [-p]) is disabled since it is useless for d=1.\n");
          slicing = pruning = false;
        }
      }
      if (! N) shrink = true, N = UINT_MAX; // enable shrinkage
      if (std::strcmp (com, "--") == 0) return;
      if (argc < optind + 3) {
        printCredit ();
        print_err ("Type `%s --help' for option details.\n", com);
      }
      train = argv[optind];
      model = argv[++optind];
      test  = argv[++optind];
      setMode  (); // induce appropriate mode
    }
    void setMode () {
      if (std::strcmp (train, "-") == 0 && std::strcmp (test, "-") == 0)
        print_err ("specify at least training or test file.\n");
      else if (std::strcmp (test,  "-") == 0) mode = TRAIN;
      else if (std::strcmp (train, "-") == 0) mode = TEST;
      else                                    mode = BOTH;
      if (std::strcmp (model, "-") == 0 && mode != BOTH)
        print_err ("instant mode can be used w/ both train/test files.\n");
      if (mode == TRAIN && output == 1)
        print_err ("per-iteration testing requires test file.\n");
      const char * mode0 [] = {"BOTH", "TRAIN", "TEST"};
      std::fprintf (stderr, "mode: %s\n", mode0[mode]);
    }
    void printCredit () { std::fprintf (stderr, OPAL_COPYRIGHT, com); }
    void printHelp   () { std::fprintf (stderr, OPAL_OPT OPAL_OPT_TRAIN0 OPAL_OPT_TRAIN1 OPAL_OPT_TRAIN2 OPAL_OPT_TEST OPAL_OPT_MISC); }
  };
  // label bag: assuming # classes small; you can use tr1::unordered_map
  class lmap {
  private:
    struct less_charp {
      bool operator () (const char * a, const char * b) const
      { return std::strcmp (a, b) < 0; }
    };
    typedef std::map <const char *, label_t, less_charp>     l2li_t;
    typedef std::vector <std::pair <label_t, const char *> > li2l_t;
    l2li_t _l2li;
    li2l_t _li2l;
  public:
    lmap  () : _l2li (), _li2l () {}
    ~lmap () {
      for (l2li_t::const_iterator it = _l2li.begin (); it != _l2li.end (); ++it)
        delete [] it->first;
    }
    label_t set_id (const char * ys, size_t len = 0) {
      l2li_t::const_iterator it = _l2li.find (ys);
      if (it == _l2li.end ()) {
        if (! len) len = std::strlen (ys);
        char * copy = new char[len + 1];
        std::strcpy (copy, ys);
        long li = static_cast <long> (_li2l.size ());
        it = _l2li.insert (l2li_t::value_type (copy, li)).first;
        _li2l.push_back (li2l_t::value_type (li, copy));
      }
      return it->second;
    }
    label_t get_id (const char * ys) const {
      l2li_t::const_iterator it = _l2li.find (ys);
      return it != _l2li.end () ? it->second : -1;
    }
    const char * get_label (const size_t i) const { return _li2l[i].second; }
    void read (char * p, char * const p_end) { // # labels:
      if (std::strncmp (p, "# labels: ", 10) != 0)
        print_err ("premature label definition.\n");
      p += 9;
      while (++p) {
        char * ys = p; while (p != p_end && *p != ' ') ++p; *p = '\0';
        set_id (ys, static_cast <size_t> (p - ys));
        if (p == p_end) break;
      }
    }
    void write (FILE * fp) {
      std::fprintf (fp, "# labels:");
      for (uint i = 0; i < _li2l.size (); ++i)
        std::fprintf (fp, " %s", _li2l[i].second);
      std::fprintf (fp, "\n");
    }
    uint nclass () const { return static_cast <uint> (_li2l.size ()); }
  };
  // feature bag; pack spare feature indices into ordered dense ones
  //              to minimize memory consumption / trie retrieval cost
  class fmap {
  private:
    typedef std::vector <std::pair <uint, uint> > counter_t;
    std::vector <uint> _fn2fi;  // replace this with unordered_map if too sparse
    std::vector <uint> _fi2fn;
    counter_t          _counter;
  public:
    fmap () : _fn2fi (), _fi2fn (), _counter () { _fi2fn.push_back (0); }
    void revertFv2Fv  (fv_t* fv) const {
      for (fv_t::iterator it = fv->begin (); it != fv->end (); ++it)
        *it = _fi2fn[*it];
      std::sort (fv->begin (), fv->end ());
    }
    void convertFv2Fv (fv_t* fnv) const {
      fv_t::iterator jt = fnv->begin ();
      for (fv_it it = jt; it != fnv->end (); ++it)
        if (uint fi = _fn2fi[*it < _fn2fi.size () ? *it : 0])
          *jt = fi, ++jt;
      fnv->erase (jt, fnv->end ());
      std::sort (fnv->begin (), jt);
    }
    void convertFv2Fv (fv_t* fnv, uint maxF) { // register if needed
      if (maxF >= _fn2fi.size ()) _fn2fi.resize (maxF + 1, 0); // widen
      fv_t::iterator jt (fnv->begin ());
      for (fv_it it = jt; it != fnv->end (); ++it, ++jt) {
        uint &fi = _fn2fi[*it];
        if (fi == 0) // register
          fi = static_cast <uint> (_fi2fn.size ()), _fi2fn.push_back (*it);
        *jt = fi;
      }
      std::sort (fnv->begin (), jt);
    }
    void build () { // build feature mapping; empty element should be removed
      if (_fn2fi.size () < _counter.size ()) // widen if needed
        _fn2fi.resize (_counter.size (), 0);
      std::sort (_counter.rbegin (), _counter.rend ());
      for (counter_t::const_iterator it = _counter.begin ();
           it != _counter.end () && it->first; ++it) {
        uint &fi = _fn2fi[it->second];
        if (fi == 0)
          fi = static_cast <uint> (_fi2fn.size ()), _fi2fn.push_back (it->second);
      }
      counter_t ().swap (_counter);
    }
    void inc_count (const fv_t *fv, uint maxF) { // sorted features assumed
      for (uint i = static_cast <uint> (_counter.size ()); i <= maxF; ++i)
        _counter.push_back (counter_t::value_type (0, i));
      for (fv_it it = fv->begin (); it != fv->end (); ++it)
        ++_counter[*it].first;
    }
    void clear () { _fn2fi.clear (); _fi2fn.clear (); _counter.clear (); } // *
  };
  template <typename T, typename U>
  class ex_base {
  private:
    T* _derived () { return static_cast <T*> (this); }
  protected:
    U     _y;    // label or weight
    fv_t* _vec;  // (binary) feature vector
    ~ex_base () {}
  public:
    ex_base ()                        : _y (0),     _vec (0)       {}
    ex_base (const U &y_, fv_t* vec_) : _y (y_),    _vec (vec_)    {}
    ex_base (const ex_base &ex)       : _y (ex._y), _vec (ex._vec) {}
    ex_base& operator= (const ex_base &ex) { _y = ex._y; _vec = ex._vec; return *this; }
    void set (char * ex, char * const ex_end, fv_t& vec, bool store, lmap * lm, fmap * fm = 0) {
      // set label (weight) of feature vectors (support vectors)
      char * p = ex;
      read (p, lm); // _y
      // set features
      vec.clear ();
      while (p != ex_end) {
#if ! defined (USE_HASH_TRIE) && ! defined (USE_MAP_TRIE)
        static cedar::da <int, -1, -2> s2ibag; // cache std::strtol ()
        const char * q = ++p;
        *ex_end = ':'; while (*p != ':') ++p; *p = '\0';
        int &fi = s2ibag.update (q, static_cast <size_t> (p - q));
        if (! fi) fi = strton <int> (q, &p);
        if (*p)   print_err ("illegal feature index: %s\n", q);
#else
        ++p;
        int fi = strton <int> (p, &p);
#endif
        vec.push_back (static_cast <uint> (fi));
        *ex_end = ' '; while (*p != ' ') ++p;
      }
      _vec = store ? new fv_t (vec) : &vec;
      if (fm) fm->inc_count (_vec, maxF ());
    }
    void read  (char * &p, lmap *lm) { _derived ()->read_impl (p, lm); }
    // interface
    const fv_it  begin   () const { return _vec->begin (); }
    const fv_it  end     () const { return _vec->end   (); }
    const fv_rit rbegin  () const { return _vec->rbegin (); }
    const fv_rit rend    () const { return _vec->rend   (); }
    uint         maxF    () const { return _vec->empty () ? 0 : _vec->back (); }
    fv_t*        body    ()       { return _vec; }
    const fv_t*  getBody () const { return _vec; }
    size_t       getSize () const { return _vec->size (); }
  };
  struct ex_t : public ex_base <ex_t, label_t> {
#ifdef USE_MULTICLASS
    void read_impl  (char * &p, lmap* lm) {
      char * ys = p; while (*p && *p != ' ' && *p != '\t') ++p; *p = '\0';
      _y = lm->set_id (ys, p - ys);
#else
    void read_impl  (char * &p, lmap*) {
      _y = std::strtol (p, &p, 10);
#endif
    }
    void set_ex (const label_t y, fv_t &vec, bool store, fmap* fm = 0) {
      _y = y;
      _vec = store ? new fv_t (vec) : &vec;
      if (fm) fm->inc_count (_vec, maxF ());
    }
    const label_t& getLabel () const { return _y; }
    fl_t           getNorm  () const { return static_cast <fl_t> (_vec->size ()); }
    ex_t (const label_t &y_, fv_t* vec_) : ex_base <ex_t, label_t> (y_, vec_) {}
    ex_t () : ex_base <ex_t, label_t> () {}
  };
  template <typename T, typename U>
  class sv_base : public ex_base <T, U> {
  private:
    T *      _derived  ()          { return static_cast <T*> (this); }
  public:
    const U& getWeight () const    { return this->_y; }
    U&       weight    ()          { return this->_y; }
    void     print     (FILE * fp) {
      _derived ()->print_weight (fp);
      for (fv_it it = this->_vec->begin (); it != this->_vec->end (); ++it)
        std::fprintf (fp, " %d:1", *it);
      std::fprintf (fp, "\n");
    }
  protected:
    sv_base (const U &y_, fv_t* vec_) : ex_base <T, U> (y_, vec_) {}
    sv_base  () : ex_base <T, U> () {}
    ~sv_base () {}
  };
  struct sv_t : public sv_base <sv_t, w_t> {
#ifdef USE_MULTICLASS
    void read_impl (char * &p, lmap* lm) {
      _y.resize (lm->nclass (), 0.0);
      for (size_t i = 0; i < _y.size (); ++p, ++i)
        _y[i] = std::strtod (p, &p);
      --p;
    }
    void print_weight (FILE * fp) const {
      std::fprintf (fp, "%.16g", _y[0]);
      for (size_t i = 1; i < _y.size (); ++i)
        std::fprintf (fp, " %.16g", _y[i]);
    }
#else
    void read_impl    (char * &p, lmap*)
    { _y = static_cast <fl_t> (std::strtod (p, &p)); }
    void print_weight (FILE * fp) const  { std::fprintf (fp, "%.16g", _y); }
#endif
    sv_t (const w_t &y_, fv_t* vec_) : sv_base <sv_t, w_t> (y_, vec_) {}
    sv_t () : sv_base <sv_t, w_t> () {}
  };
  template <template <class> class T, typename U>
  class basic_pool {
  public:
    typedef U elem_t;
    elem_t * init    ()                { return derived ()->init_impl (); }
    void     put     (const elem_t &x) { derived ()->put_impl (x); }
    elem_t * get     ()                { return derived ()->get_impl (); }
    void     setup   (fmap &fm)        { derived ()->setup_impl (fm); }
    virtual void shuffle ()
    { std::fprintf (stderr, "skip shuffling; try ./rand_shuf.\n"); }
    virtual void read (const char* lfn, const size_t M, bool mem, lmap *lm, fmap *fm = 0, const char skip = 0) {
      FILE * fp = std::fopen (lfn, "r");
      if (! fp) print_err ("no such file: %s\n", lfn);
      char buf[BUF_SIZE]; std::setvbuf (fp, &buf[0], _IOFBF, BUF_SIZE);
      char * line = 0;
      size_t read (0), i (0);
      while (getLine (fp, line, read)) {
        if (skip && std::memchr (line, skip, read - 1)) continue;
        if (M && ++i > M) break;
        _x.set (line, line + read - 1, _vec, mem, lm, fm);
        derived ()->put_impl (_x);
      }
      std::fclose (fp);
    }
  protected:
    basic_pool  () : _x (), _vec () {}
    ~basic_pool () {}
    elem_t _x;
    fv_t   _vec;
  private:
    T <U> * derived () { return static_cast <T <U>*> (this); }
    basic_pool (const basic_pool&);
    basic_pool& operator= (const basic_pool&);
  };
  template <typename elem_t>
  class null_pool : public basic_pool <null_pool, elem_t> {
  public:
    null_pool  () : _fp (0), _line (0), _buf (), _read (0), _i (0), _M (0), _lm (0), _fm (0), _s () {}
    virtual ~null_pool () {}
    elem_t * init_impl () {
      _i = 0; _line = 0; _read = 0;
      if (_fp != stdin) std::fseek (_fp, 0, SEEK_SET);
      return get_impl ();
    }
    void read (const char * lfn, const size_t M, bool mem, lmap *lm, fmap *fm = 0, char = 0) {
      _fp = lfn ? std::fopen (lfn, "r") : stdin;  // initialize
      _M = M; _lm = lm; _fm = fm;
      if (! _fp) print_err ("no such file: %s\n", lfn);
      std::setvbuf (_fp, &_buf[0], _IOFBF, BUF_SIZE);
      if (_fm) {
        _i = 0;
        while (getLine (_fp, _line, _read)) {
          if (_M && ++_i > _M) break;
          this->_x.set (_line, _line + _read - 1, this->_vec, mem, lm, fm);
        }
      }
    }
    void     put_impl (const elem_t&) {} // ignore
    elem_t * get_impl () {
      if (_M && ++_i > _M)               return 0;
      if (! getLine (_fp, _line, _read)) return 0;
      this->_x.set (_line, _line + _read - 1, this->_vec, false, _lm, _fm);
      if (_fm) _fm->convertFv2Fv (this->_x.body ());
      return &this->_x;
    }
    void setup_impl (fmap &fm) { _fm = &fm; } // ignore
    FILE *  _fp;
    char *  _line;
    char    _buf[BUF_SIZE];
    size_t  _read;
    size_t  _i;
    size_t  _M;
    lmap *  _lm;
    fmap *  _fm;
  private:
    std::vector <uchar> _s;
  };
  template <typename elem_t>
  class mem_pool : public basic_pool <mem_pool, elem_t> {
  public:
    typedef typename std::vector <elem_t>  data_t;
    mem_pool () : _ex (), _pos (0), _gen () {}
    virtual ~mem_pool () {
      for (typename data_t::iterator it = _ex.begin (); it != _ex.end (); ++it)
        delete it->getBody ();
    }
    elem_t * init_impl () { _pos = 0; return get_impl (); }
    void put_impl (const elem_t &x) { _ex.push_back (x); }
    elem_t *  get_impl () { return (_pos == _ex.size ()) ? 0 : &_ex[_pos++]; }
    void setup_impl (fmap &fm) {
      for (typename data_t::iterator it = _ex.begin (); it != _ex.end (); ++it)
        fm.convertFv2Fv (it->body ());
    }
    void shuffle    () { std::random_shuffle (_ex.begin (), _ex.end (), _gen); }
  private:
    data_t _ex;
    size_t _pos;
    struct rand_ {
#ifdef USE_MT19937
      std::tr1::variate_generator <std::tr1::mt19937,
                                   std::tr1::uniform_int <size_t> > gen;
      rand_ () : gen (std::tr1::mt19937 (), std::tr1::uniform_int <size_t> ()) {}
      size_t operator () (size_t max) { return gen (max); }
#else
      size_t gen () { // Xorshift RNG; http://www.jstatsoft.org/v08/i14/paper
        static size_t x (123456789), y (362436069), z (521288629), w (88675123);
        size_t t = (x ^ (x << 11)); x = y; y = z; z = w;
        return w = (w ^ (w >> 19)) ^ (t ^ (t >> 8));
      }
      size_t operator () (size_t max) { return gen () % max; }
#endif
    } _gen;
  };
  template <typename elem_t>
  class disk_pool : public basic_pool <disk_pool, elem_t> {
  public:
    disk_pool  () : fp (std::tmpfile ()), _s () {}
    virtual ~disk_pool () { std::fclose (fp); }
    elem_t * init_impl () { std::fseek (fp, 0, SEEK_SET); return get_impl (); }
    void put_impl (const elem_t &x) {
      if (_s.size () < x.getSize () * KEY_SIZE)
        _s.resize (x.getSize () * KEY_SIZE);
      const label_t y = x.getLabel (); // prohibit elem_t = sv_t
      uint len (0), prev (0);
      byte_encoder encoder;
      for (fv_it it = x.begin (); it != x.end (); prev = *it, ++it)
        len += encoder.encode (*it - prev, &_s[len]);
      std::fwrite (&y,     sizeof (label_t), 1,   fp);
      std::fwrite (&len,   sizeof (uint),    1,   fp);
      std::fwrite (&_s[0], sizeof (uchar),   len, fp);
    }
    elem_t * get_impl () {
      label_t y   = 0; if (std::fread (&y, sizeof (long), 1, fp) == 0) return 0;
      uint    len = 0; std::fread (&len, sizeof (uint), 1, fp);
      if (_s.size () < len) _s.resize (len);
      std::fread (&_s[0], sizeof (uchar), len, fp);
      this->_vec.clear ();
      for (uint i (0), r (0), b (0), prev (0); i < len; ++i) {
        r += ((static_cast <uint> (_s[i]) & 0x7f) << b);
        if (_s[i] & 0x80) b += 7; else this->_vec.push_back (prev += r), r = b = 0;
      }
      this->_x.set_ex (y, this->_vec, false);
      return &this->_x;
    }
    void setup_impl (fmap &fm) {
      disk_pool <elem_t> tmp;
      for (elem_t * x = this->init (); x; x = this->get ()) {
        fm.convertFv2Fv (x->body ());
        tmp.put (*x);
      }
      std::swap (fp, tmp.fp);
    }
    FILE * fp;
  private:
    std::vector <uchar> _s;
  };
  class Model {
  private:
    // type alias
    typedef std::vector <uint>    ss_t;
    typedef ss_t::const_iterator  ss_it;
    typedef std::vector <w_t>     wv_t;
    typedef std::vector <ss_t>    f2ss_t;
    //
    static const int MAX_TRIAL = 4;  // 4 is good enough;
#ifdef USE_HASH
    struct hash_fv : std::unary_function <const fv_t*, size_t> {
      result_type operator() (argument_type f) const
      { return std::accumulate (f->begin (), f->end (), FNV_BASIS, inc_fnv ()); }
    };
    struct eq_fv : std::binary_function <const fv_t*, const fv_t*, bool> {
      result_type operator() (first_argument_type a, second_argument_type b)
        const { return *a == *b; } // can be faster
    };
    typedef std::tr1::unordered_map <const fv_t*, uint, hash_fv, eq_fv> fv2s_t;
#else
    struct less_fv : std::binary_function <const fv_t*, const fv_t*, bool> {
      result_type operator() (first_argument_type a, second_argument_type b)
        const { return *a < *b; }
    };
    typedef std::map <const fv_t*, uint, less_fv> fv2s_t;
#endif
    // variables
    const option       _opt;
    lmap               _lm;        // label map
    fmap               _fm;        // feature map
    size_t             _nex;       // # processed examples
    uint               _nsv;       // # support vectors
    uint               _nf;        // # features (= max feature id)
    w_t                _m0;        // margin0
    w_t                _bias;      // just keep constant for kernel expansion
    f2ss_t             _f2ss;      // feature id        -> support vector id
    fv2s_t             _fv2s;      // feature vector    -> support vector id
    std::vector <bool> _fbit;      // I don't like vector <bool>, though
    wv_t               _alpha;     // support vector id -> (averaged) alpha
    std::vector <sv_t> _sv;        // support vector id -> support vector
    wv_t               _w;         // feature id        -> weight
    wv_t               _wa;        // feature id        -> (average) weight
    std::vector <fl_t> _polyk;     // poly kernel;        int -> fl_t
    std::vector <fl_t> _spolyk;    // poly kernel sliced; int -> fl_t
#ifdef USE_HASH_TRIE
    typedef std::tr1::unordered_map <edge, node> trie_t;
#endif
#ifdef USE_MAP_TRIE
    typedef std::map <edge, node> trie_t;
#endif
#if ! defined (USE_HASH_TRIE) && ! defined (USE_MAP_TRIE)
    typedef cedar::da <int, -1, -2, MAX_TRIAL, MAX_KERNEL_DEGREE> trie_t;
#endif
    uint               _nbin;  // # feature weight tries
    trie_t **          _ftrie; // feature weight tries
    wv_t *             _fw;
    struct pn_t { // upper- and lower-bounds of weights
      w_t pos; // or max
      w_t neg; // or min
      pn_t (const w_t &pos_, const w_t &neg_) : pos (pos_), neg (neg_) {};
      pn_t& operator= (const pn_t &b) { pos = b.pos; neg = b.neg; return *this; }
    };
    std::vector <pn_t> _f2pn;  // feature id -> weight summation of pos/neg svs
    std::vector <pn_t> _bound;
    struct tpm_t { // TODO: we can mark and sweep older nodes
      uint t;   // time signature
      w_t  tpm; // temporal partial margin
      tpm_t (const uint t_, const w_t &tpm_) : t (t_),   tpm (tpm_)   {};
      tpm_t (const tpm_t &l_)                : t (l_.t), tpm (l_.tpm) {};
      tpm_t& operator= (const tpm_t &l) { t = l.t; tpm = l.tpm; return *this; }
    };
    class ring { // circular buffer
    public:
      uint end;  // current pos in the ring
      uint t;    // # updates
      ring () : end (0), t (0), _M () {}
      const sv_t& operator[] (const size_t index) { return _M[index]; }
      size_t size   () const { return _M.size (); }
      void   resize (const size_t n, const sv_t &sv0)
      { _M.insert (_M.begin () + end, n - _M.size (), sv0); }
      void   push   (const sv_t &s, bool overt = true) {
        if (! _M.empty () && overt) { // circulate to overwrite
          if (end >= static_cast <uint> (_M.size ())) end = 0;
          _M[end] = s;
        } else {
          _M.insert (_M.begin () + end, s);
        }
        ++end; ++t;
      }
    private:
      std::vector <sv_t> _M;
    };
    std::vector <ring>   _f2ss_up;
    trie_t               _pmtrie;
    std::vector <tpm_t>  _tpm;
    std::vector <size_t> _limk; // condtion for reuse
    const fl_t _thresh;
#ifdef USE_TIMER
    ny::TimerPool _timer_pool;
    ny::Timer *   _timer;
#endif
    // uncopyable
    Model (const Model&);
    Model& operator= (const Model&);
    // internal functions
    void _precompute_kernel (const size_t nf) {
      for (uint i = static_cast <uint> (_spolyk.size ()); i <= nf; ++i) {
        const fl_t j (static_cast <fl_t> (i)), d (static_cast <fl_t> (_opt.d));
        _spolyk.push_back (std::pow ((j + 1) + 1, d) - std::pow (j + 1, d));
#ifdef USE_POLYK
        _polyk.push_back  (std::pow (j + 1, d));
#endif
      }
    }
    void _getMargin (w_t &m, const fv_it first, const fv_it last) const {
      m = 0.0;
      for (fv_it it = first; it != last && *it <= _nf; ++it)
        m += _w[*it];
    }
    void _getMarginPoly (w_t &m, const fv_it first, const fv_it last, const label_t y = MAX_NUM_CLASSES) {
      m = 0.0;
#ifdef USE_POLYK
      // polynomial kernel inverted
      std::vector <uint> dot (_nsv, 0);
      for (fv_it it = first; it != last && *it <= _nf; ++it) {
        const ss_t &ss = _f2ss[*it];
        for (ss_it st = ss.begin (); st != ss.end (); ++st) ++dot[*st];
      }
      for (uint i = 0; i < _nsv; ++i)
        m += _sv[i].getWeight () * _polyk[dot[i]];
#else
      if (_opt.slicing)
        while (_limk.size () <= static_cast <size_t> (last - first))
          switch (_opt.d) {
            case 1: _limk.push_back (0); break;
            case 2: _limk.push_back (1); break;
            case 3: _limk.push_back ((_limk.size () >> 1) + (_limk.size () & 0x1)); break;
          }
      _project_ro (m, first, last, y);
#endif
    }
    // update
    void _addTo (const ex_t &x, const w_t &t) {
      if (x.maxF () > _nf) { // widen
        _nf = x.maxF ();
        _w.resize (_nf + 1, _m0);
        if (_opt.average) _wa.resize (_nf + 1, _m0);
      }
      for (fv_it it = x.begin (); it != x.end (); ++it) {
        _w[*it] += t;
        if (_opt.average) _wa[*it] += t * static_cast <fl_t> (_nex);
      }
    }
    template <typename elem_t>
    void _pushTo (const elem_t &x, const w_t &t) {
      // check support vector overlap; a little overhead
      fv2s_t::const_iterator fit =_fv2s.find (x.getBody ());
#ifndef USE_POLYK
      const uint MIN_N
        = _opt.shrink ? (1U << PSEUDO_TRIE_N[_opt.d]) - 1 : _opt.N;
#endif
      if (fit == _fv2s.end ()) { // add support vector
        fit = _fv2s.insert (fv2s_t::value_type (new fv_t (*x.getBody ()), _nsv)).first;
        if (x.maxF () > _nf) {
          _nf = x.maxF ();
          _f2ss.resize (_nf + 1);
          if (_opt.pruning)
            _f2pn.resize (_nf + 1, pn_t (_m0, _m0));
          if (_opt.slicing)
            _f2ss_up.resize (_nf + 1); // processed support vector (signature)
          _fbit.resize (_nf + 1, false);
        }
        _precompute_kernel (x.getSize ());
#ifdef USE_POLYK
        for (fv_rit it = x.rbegin (); it != x.rend (); ++it)
#else
        for (fv_rit it = x.rbegin (); it != x.rend () && *it > MIN_N; ++it)
#endif
          _f2ss[*it].push_back (_nsv);
        _sv.push_back (sv_t (t, const_cast <fv_t *> (fit->first))); ++_nsv;
        if (_opt.average) _alpha.push_back (t * static_cast <fl_t> (_nex));
      } else { // just update
        _sv[fit->second].weight () += t;
        if (_opt.average) _alpha[fit->second] += t * static_cast <fl_t> (_nex);
      }
#ifndef USE_POLYK
      if (_opt.slicing)
        for (fv_it it = x.begin (); it != x.end (); ++it)
          _f2ss_up[*it].push (sv_t (t, const_cast <fv_t *> (fit->first)),
                              _f2ss_up[*it].size () >= _f2ss[*it].size () ||
                              *it <= MIN_N);
      if (_opt.pruning)
        for (fv_it it = x.begin (); it != x.end (); ++it) {
          pn_t &fb = _f2pn[*it];
#ifdef USE_MULTICLASS
          for (uint i = 0; i < _opt.nclass; ++i)
            if (t[i] > 0) fb.pos[i] += t[i]; else fb.neg[i] += t[i];
#else
          if (t > 0) fb.pos += t; else fb.neg += t;
#endif
        }
      _project_rw (t, x.begin (), x.end ()); // explicit conjunctive features
      if (_opt.shrink) _shrink_trie ();
#endif
    }
    bool _reuse_pm (w_t &m_, size_t &pid, const size_t i, const uint fi, const uint prev, tpm_t* &l, const bool projected = true) {
      // retrieve partial margin; TODO: avoid expanding trie when unnecessary
      ring &ss = _f2ss_up[fi];
#if defined (USE_HASH_TRIE) || defined (USE_MAP_TRIE)
      const edge e ((static_cast <uint64_t> (pid) << 32) | fi);
      trie_t::iterator pit = _pmtrie.insert (trie_t::value_type (e, node (static_cast <uint> (_pmtrie.size () + 1)))).first;
      pid = pit->second.id;
      int& n = pit->second.value;
#else
      byte_encoder encoder (fi - prev);
      int& n = _pmtrie.update (encoder.key (), pid, encoder.len (), 0);
#endif
      if (! n) { // new feature sequence
        _tpm.push_back (tpm_t (0, _m0)); n = static_cast <int> (_tpm.size ());
        if (projected && ss.size () < _limk[i])
          ss.resize (_limk[i], sv_t (_m0, 0));
      }
      l = &_tpm[static_cast <size_t> (n - 1)]; // too much
      uint nusv = ss.t - l->t; // # newly updated support vectors
      l->t = ss.t;
      if (! nusv) { m_ += l->tpm; return true; }  // TODO: sweep older weights
      if (nusv > (projected ? _limk[i] : _f2ss[fi].size ())) return false;
      uint p = ss.end >= nusv ? ss.end - nusv :
               static_cast <uint> (ss.size () - (nusv - ss.end));
      if (! ss[p].getBody ()) return false; // not filled
      m_ += l->tpm;
      for (; nusv > 0; --nusv, ++p) { // reuse temporal partial margin
        if (p == ss.size ()) p = 0; // loop; dirty bit
        const sv_t &s = ss[p];
        uint dot_c = 0;
        for (fv_it sit = s.begin (); *sit <= prev; ++sit)
          dot_c += _fbit[*sit];
        m_ += s.getWeight () * _spolyk[dot_c];
      }
      return true;
    }
    void _update_weight (int &id, wv_t &fw, const w_t &t) {
      if (id) fw[static_cast <size_t> (id - 1)] += t; // too much
      else    fw.push_back (t), id = static_cast <int> (fw.size ());
    }
#if defined (USE_HASH_TRIE) || defined (USE_MAP_TRIE)
    void _traverse_rw (const w_t &t, size_t &pid, const uint fi, trie_t * trie, wv_t &fw) {
      const edge e ((static_cast <uint64_t> (pid) << 32) | fi);
      trie_t::iterator tit
        = trie->insert (std::make_pair (e, node (static_cast <uint> (trie->size () + 1)))).first;
      _update_weight (tit->second.value, fw, t);
      pid = tit->second.id; // next node
    }
    bool _traverse_ro (w_t &m, size_t &pid, const uint fi, const trie_t* trie, const wv_t &fw) {
      const edge e ((static_cast <uint64_t> (pid) << 32) | fi);
      trie_t::const_iterator tit = trie->find (e);
      if (tit == trie->end ()) return false;
      m += fw[static_cast <size_t> (tit->second.value - 1)];
      pid = tit->second.id; // next node
      return true;
    }
#else
    void _traverse_rw (const w_t &t, size_t &pid, const uint fi, trie_t* trie, wv_t &fw) {
      byte_encoder encoder (fi);
      _update_weight (trie->update (encoder.key (), pid, encoder.len (), 0), fw, t);
    }
    bool _traverse_ro (w_t &t, size_t &pid, const uint fi, const trie_t* trie, const wv_t &fw) {
      size_t p = 0;
      byte_encoder encoder (fi);
      const int n = trie->traverse (encoder.key (), pid, p, encoder.len ());
      if (n == trie->NO_PATH)  return false;
      if (n != trie->NO_VALUE) t += fw[static_cast <size_t> (n - 1)]; // too much
      return true;
    }
#endif
    void _estimate_bound (const fv_it first, fv_it tail) {
      // compute bound from right-to-left
      const size_t len = static_cast <size_t> (std::distance (first, tail));
      if (_bound.size () < len) _bound.resize (len, pn_t (_m0, _m0));
      pn_t * p = &_bound[len - 1];
      for (p->pos = p->neg = 0;; *(p-1) = *p, --p) {
        const pn_t &fb = _f2pn[*--tail];
        const size_t max_dot // max inner product between common features
          = std::min (_spolyk.size () - 1,
                      static_cast <size_t> (std::distance (first, tail)));
        p->pos += fb.pos * _spolyk[max_dot] + fb.neg * _spolyk[0];
        p->neg += fb.neg * _spolyk[max_dot] + fb.pos * _spolyk[0];
        if (tail == first) break;
      }
    }
    bool _prune (w_t &m, const size_t i, label_t y) {
      const pn_t &b = _bound[i];
      const fl_t thresh = y == MAX_NUM_CLASSES ? 0 : _thresh;
#ifdef USE_MULTICLASS
      if (y == MAX_NUM_CLASSES)
        y = std::max_element (&m[0], &m[_opt.nclass]) - &m[0];
      for (uint j = 0; j < _opt.nclass; ++j)
        if (j != y && (m[y] + b.neg[y]) - (m[j] + b.pos[j]) <= thresh)
          return false;
      for (uint j = 0; j < _opt.nclass; ++j)
        m[j] += j == y ? b.neg[j] : b.pos[j];
#else
      if (y == MAX_NUM_CLASSES)
        y = m >= 0 ? 1 : -1;
      if ((y > 0 ? m + b.neg : - (m + b.pos)) <= thresh)
        return false;
      m += y > 0 ? b.neg : b.pos;
#endif
      return true;
    }
    void _project_rw (const w_t &t, const fv_it first, const fv_it last) {
      _bias += t;
      const uint   N   = std::min ((1U << PSEUDO_TRIE_N[_opt.d]) - 1, _opt.N);
      const fv_it  rit = std::upper_bound (first, last, _opt.N);
      fv_it it = first;
      for (uint bin = PSEUDO_TRIE_N[_opt.d]; it != rit; ++it) {
        if (*it <= N) {
          size_t p[MAX_KERNEL_DEGREE];
          switch (_opt.d) {
            case 1: _w[*it - 1] += t; break;
            case 2:
              _w[p[0] = *it * (*it - 1) / 2] += t;
              for (fv_it jt = first; jt != it; ++jt)
                _w[p[0] + 1 + *jt - 1] += t;
              break;
            case 3:
              _w[p[0] = (*it - 1) * (*it * *it - 2 * *it + 6) / 6] += t;
              for (fv_it jt = first; jt != it; ++jt) {
                _w[p[1] = p[0] + 1 + *jt * (*jt - 1) / 2] += t;
                for (fv_it kt = first; kt != jt; ++kt)
                  _w[p[1] + 1 + *kt - 1] += t;
              }
          }
        } else {
          while (*it >> (bin + 1)) ++bin;
          trie_t * const trie = _ftrie[bin];
          wv_t           &fw  = _fw[bin];
#if defined (USE_HASH_TRIE) || defined (USE_MAP_TRIE)
          size_t p[MAX_KERNEL_DEGREE];
#else
          size_t * const p (&(trie->pid[0])); // traversed node might be moved
#endif
          _traverse_rw (t, p[0] = 0, *it, trie, fw);
          if (_opt.d > 1)
            for (fv_it jt = first; jt != it; ++jt) {
              _traverse_rw (t, p[1] = p[0], *jt, trie, fw);
              if (_opt.d > 2)
                for (fv_it kt = first; kt != jt; ++kt)
                  _traverse_rw (t, p[2] = p[1], *kt, trie, fw);
            }
        }
      }
    }
    void _project_ro (w_t &m, const fv_it first, const fv_it last, const label_t y) {
      // kernel slicing
      const fl_t * const coeff = COEFF[_opt.d];
      m += _bias * coeff[0];
      if (first == last) return; // bug fix
      const uint   N   = std::min ((1U << PSEUDO_TRIE_N[_opt.d]) - 1, _opt.N);
      const fv_it  end = std::upper_bound (first, last, _nf);
      const fv_it  rit = std::upper_bound (first, end,  _opt.N);
      if (_opt.pruning && ! _spolyk.empty ()) _estimate_bound (first, end);
      size_t  pid = 0;
      tpm_t * l   = 0;
      uint   prev = 0;
      fv_it  it   = first;
      size_t p[MAX_KERNEL_DEGREE];
#ifdef USE_MULTICLASS
      static w_t pm, pm1, pm2, pm3;
      if (pm.size () < _opt.nclass) {
        pm.resize  (_opt.nclass, 0);
        pm1.resize (_opt.nclass, 0);
        pm2.resize (_opt.nclass, 0);
        pm3.resize (_opt.nclass, 0);
      }
#else
      w_t pm (0), pm1 (0), pm2 (0), pm3 (0);
#endif
      for (uint bin = PSEUDO_TRIE_N[_opt.d]; it != rit; ++it) { // common part
        const size_t index = static_cast <size_t> (std::distance (first, it));
        if (_opt.pruning && _prune (m, index, y)) { it = end; break; }
        pm = pm1 = pm2 = pm3 = 0.0;
        if (! _opt.slicing || ! _reuse_pm (pm, pid, index, *it, prev, l)) {
          if (*it <= N) {
            switch (_opt.d) {
              case 1: pm1 += _w[*it - 1]; break;
              case 2: pm1 += _w[p[0] = *it * (*it - 1) / 2];
                for (fv_it jt = first; jt != it; ++jt)
                  pm2 += _w[p[0] + 1 + *jt - 1];
                break;
              case 3: pm1 += _w[p[0] = (*it - 1) * (*it * *it - 2 * *it + 6) / 6];
                for (fv_it jt = first; jt != it; ++jt) {
                  pm2 += _w[p[1] = p[0] + 1 + *jt * (*jt - 1) / 2];
                  for (fv_it kt = first; kt != jt; ++kt)
                    pm3 += _w[p[1] + 1 + *kt - 1];
                }
            }
          } else {
            while (*it >> (bin + 1)) ++bin;
            const trie_t * const trie = _ftrie[bin];
            const wv_t           &fw  = _fw[bin];
            if (_traverse_ro (pm1, p[0] = 0, *it, trie, fw) && _opt.d > 1)
              for (fv_it jt = first; jt != it; ++jt)
                if (_traverse_ro (pm2, p[1] = p[0], *jt, trie, fw) && _opt.d > 2)
                  for (fv_it kt = first; kt != jt; ++kt)
                    _traverse_ro (pm3, p[2] = p[1], *kt, trie, fw);
          }
          pm += pm1 * coeff[1] + pm2 * coeff[2] + pm3 * coeff[3];
        }
        if (_opt.slicing) l->tpm = pm;
        m += pm; prev = *it; _fbit[*it] = true;
      }
      for (; it != end; ++it) { // rare part
        const size_t index = static_cast <size_t> (std::distance (first, it));
        if (_opt.pruning && _prune (m, index, y)) { it = end; break; }
        pm = 0.0;
        if (! _opt.slicing || ! _reuse_pm (pm, pid, index, *it, prev, l, false)) {
          const ss_t &ss = _f2ss[*it];
          for (ss_it st = ss.begin (); st != ss.end (); ++st) {
            const sv_t &s    = _sv[*st];
            uint       dot_c = 0;
            for (fv_it sit = s.begin (); *sit <= prev; ++sit)
              dot_c += _fbit[*sit];
            pm += s.getWeight () * _spolyk[dot_c];
          }
        }
        if (_opt.slicing) l->tpm = pm;
        m += pm; prev = *it; _fbit[*it] = true;
      }
      for (fv_it jt = first; jt != it; ++jt) _fbit[*jt] = false; // reset
    }
    void _shrink_trie () {
      size_t rss = 0;
      for (uint i (PSEUDO_TRIE_N[_opt.d]); i < _nbin; ++i) {
#if    defined (USE_HASH_TRIE)
        typedef std::tr1::__detail::_Hash_node <trie_t::value_type, false> _Node;
        rss += sizeof (_Node *) * _ftrie[i]->bucket_count ();
        rss += sizeof (_Node)   * _ftrie[i]->size (); // ignore alignment
#elif defined (USE_MAP_TRIE)
        typedef std::_Rb_tree_node <trie_t::value_type> _Node;
        rss += sizeof (_Node)   * _ftrie[i]->size (); // ignore alignment
#else
        rss += sizeof (trie_t::node)  * _ftrie[i]->capacity ();
        rss += sizeof (trie_t::ninfo) * _ftrie[i]->capacity ();
        rss += sizeof (trie_t::block) * _ftrie[i]->capacity () >> 8;
#endif
        if (rss > _opt.T) {
          std::fprintf (stderr,
                        "shrink splitN: 2^%d-1 (= %ld) => 2^%d-1 (= %d)\n",
                        _nbin, (1UL << _nbin) - 1, i, (1 << i) - 1);
          while (_nbin > i) {
            --_nbin;
            delete _ftrie[_nbin];
            wv_t ().swap (_fw[_nbin]);
          }
          _opt.N = (1U << _nbin) - 1;
          break;
        }
      }
    }
  public:
    Model (const opal::option &opt) :
      _opt (opt), _lm (), _fm (), _nex (0), _nsv (0), _nf (0), _m0 (), _bias (), _f2ss (), _fv2s (), _fbit (), _alpha (), _sv (), _w (), _wa (), _polyk (), _spolyk (), _nbin (static_cast <uint> (std::floor (std::log (opt.N) / std::log (2)) + 1)), _ftrie (new trie_t*[_nbin]), _fw (new wv_t[_nbin] ()), _f2pn (), _bound (), _f2ss_up (), _pmtrie (),_tpm (), _limk (), _thresh (_opt.algo == P ? 0 : 1.0)
#ifdef USE_TIMER
      , _timer_pool (), _timer (0)
#endif
    {
#ifdef USE_MULTICLASS
      init_weight (_opt.nclass);
#endif
      for (uint i = 0; i < _nbin; ++i) _ftrie[i] = new trie_t;
    }
    ~Model () { // clean up
      printStat ();
      for (fv2s_t::iterator it = _fv2s.begin (); it != _fv2s.end (); ++it)
        delete it->first;
      if (_opt.kernel == POLY)
        for (uint i = 0; i < _nbin; ++i) delete _ftrie[i];
      delete [] _ftrie;
      delete [] _fw;
    }
    void init_weight_trie () {
      const uint N = std::min ((1U << PSEUDO_TRIE_N[_opt.d]) - 1, _opt.N);
      switch (_opt.d) {
        case 1: _w.resize (N, _m0);                   break;
        case 2: _w.resize (N * (N + 1) / 2, _m0);     break;
        case 3: _w.resize (N * (N * N + 5) / 6, _m0); break;
      }
    }
#ifdef USE_MULTICLASS
    void init_weight (const long nclass) {
      if (nclass >= MAX_NUM_CLASSES)
        print_err ("set MAX_NUM_CLASSES > %ld\n", nclass);
      _m0.resize  (nclass, 0.0);
      if (_opt.kernel == POLY) _bias.resize (nclass, 0.0);
      _opt.nclass = nclass;
    }
#endif
    // clasification interface
    void printStat () { TIMER (_timer_pool.print ()); }
    void classify (fv_t &fv, w_t * m) {
      if (_opt.kernel == LINEAR)
        _getMargin (*m, fv.begin (), fv.end ());
      else
        _fm.convertFv2Fv (&fv), _getMarginPoly (*m, fv.begin (), fv.end ());
    }
    static double sigmoid (fl_t x)
    { return 1.0 / (1.0 + std::exp (-x)); }
#ifdef USE_MULTICLASS
    double getProbability (fv_t &fv)
    { print_err ("sorry, probability output is unsupported.\n"); return 0; }
    bool binClassify (fv_t &fv, const char * label) {
      static w_t m; if (m.size () < _opt.nclass) m.resize (_opt.nclass, 0.0);
      const label_t y = _lm.get_id (label);
      if (y == -1)
        print_err ("unknown label: %s\n", label);
      classify (fv, &m);
      const label_t y_ = std::max_element (&m[0], &m[_opt.nclass]) - &m[0];
      return y_ == y;
    }
#else
    // sigmoid w/o fitting (a stupid variant of J. Plat 1999)
    double getProbability (fv_t &fv)
    { w_t m = 0; classify (fv, &m); return sigmoid (m); }
    bool binClassify (fv_t &fv)
    { w_t m = 0; classify (fv, &m); return m > 0; }
#endif
    // training interface
#ifdef USE_MULTICLASS
    void set_ex (ex_t &x, const char * ys, fv_t &vec, bool store, bool flag)
    { set_ex (x, _lm.set_id (ys), vec, store, flag); }
#endif
    void set_ex (ex_t &x, const label_t y, fv_t &vec, bool store, bool flag)
    { x.set_ex (y, vec, store, flag ? &_fm : 0); }
    void process_example (ex_t &x, bool online = true) {
      if (_opt.kernel == POLY && online)
        _fm.convertFv2Fv (x.body (), x.maxF ());
#ifdef USE_MULTICLASS
      static w_t m; if (m.size () < _opt.nclass) m.resize (_opt.nclass, 0);
#else
      w_t m = 0;
#endif
      const label_t y = x.getLabel ();
      if (_opt.kernel == LINEAR) _getMargin     (m, x.begin (), x.end ());
      else                       _getMarginPoly (m, x.begin (), x.end (), y);
      ++_nex;
#ifdef USE_MULTICLASS
      label_t y_ = y ? 0 : 1; // most probable label other than y
      for (uint i = 0; i < _opt.nclass; ++i)
        if (i != y && m[i] > m[y_]) y_ = i;
      if (m[y] - m[y_] <= _thresh) {
        w_t &t = m;
        if (_opt.algo == P) {
          t = 0.0, t[y] += 1.0, t[y_] -= 1.0;  // what's a simple
        } else {
          typedef std::vector <std::pair <fl_t, size_t> > loss_t;
          static loss_t ls;
          ls.clear ();
          // compute margin
          for (uint j = 0; j < _opt.nclass; ++j)
            if (j != y) {
              const fl_t lj = std::max (0.0, 1.0 - (m[y] - m[j]));
              if (lj > 0) ls.push_back (loss_t::value_type (lj, j));
            }
          // examine support class
          std::sort (ls.rbegin (), ls.rend ());
          size_t k = 0;
          fl_t   l = 0; // summation of loss
          for (bool is_sc = true; k < ls.size (); l += ls[k].first, ++k) {
            switch (_opt.algo) {
              case PA:  is_sc &= l < (k + 1) * ls[k].first; break;
              case PA1: is_sc &= l < std::min (k * ls[k].first + _opt.C * x.getNorm (), (k + 1) * ls[k].first); break;
              case PA2: is_sc &= l < ((k + 1) * x.getNorm () + k / (2 * _opt.C)) / (x.getNorm () + 1 / (2 * _opt.C)) * ls[k].first; break;
              default: break;
            }
            if (! is_sc) { ls.resize (k); break; }
          } // k = |S|; support class 0 -> k - 1
            // update weights
          fl_t penalty = 0;
          switch (_opt.algo) {
            case PA:  penalty = l / (k + 1); break;
            case PA1: penalty = std::max (l / k - _opt.C * x.getNorm () / k, l / (k + 1)); break;
            case PA2: penalty = (x.getNorm () + 1 / (2 * _opt.C)) / ((k + 1) * x.getNorm () + (k / (2 * _opt.C))) * l; break;
            default: break;
          }
          t = 0.0;
          for (loss_t::iterator it = ls.begin (); it != ls.end (); ++it) {
            const fl_t t_ = std::max (static_cast <fl_t> (0), it->first - penalty) / x.getNorm ();
            t[y] += t_; t[it->second] -= t_;
          }
        }
        if  (_opt.kernel == LINEAR) _addTo (x, t); else _pushTo (x, t);
      }
#else
      m *= static_cast <fl_t> (y);
      if (m <= _thresh) {
        fl_t t = static_cast <fl_t> (y > 0 ? 1 : -1); // tau
        switch (_opt.algo) {
          case P:   break;
          case PA:  t *= (1 - m) / x.getNorm (); break;
          case PA1: t *= std::min (_opt.C, (1 - m) / x.getNorm ()); break;
          case PA2: t *= (1 - m) / (x.getNorm () + 1 / (2 * _opt.C)); break;
        }
        if (_opt.kernel == LINEAR) _addTo (x, t); else _pushTo (x, t);
      }
#endif
    }
    void train_from_file (const char* lfn, const uint iter, const char * tfn = "") {
      switch (_opt.buffer) {
        case RAM:  train <mem_pool  <ex_t> > (lfn, iter, tfn); break;
        case DISK: train <disk_pool <ex_t> > (lfn, iter, tfn); break;
        case null: train <null_pool <ex_t> > (lfn, iter, tfn); break;
        default:   break;
      }
    }
    template <typename Pool>
    void train (const char* lfn, const uint iter, const char * tfn = "") {
      // read examples into pool prior to training
      Pool pool;
      const char * pool_action[] = { "loading", "caching", "preparing" };
      std::fprintf (stderr, "%s examples..", pool_action[_opt.buffer]);
      pool.read (lfn, _opt.M, _opt.buffer == RAM, &_lm, _opt.kernel == LINEAR ? 0 : &_fm);
      std::fprintf (stderr, "done.\n");
      train (pool, iter, tfn); // do
    }
    template <typename Pool>
    void train (Pool &pool, const uint iter, const char * tfn = "") {
#ifdef USE_MULTICLASS
      // # classes has been changed
      if (_opt.nclass < _lm.nclass ()) init_weight (_lm.nclass ());
#endif
      if (_opt.kernel == POLY)
        init_weight_trie (), _fm.build (), pool.setup (_fm);  // feature mapping
      if (_opt.shuffle) pool.shuffle (); // shuffling data
      Pool test_pool;
      if (std::strcmp (tfn, "") != 0) {
        test_pool.read  (tfn, 0, _opt.buffer == RAM, &_lm);
        if (_opt.kernel == POLY) test_pool.setup (_fm);
      }
      for (uint i = 1; i <= iter; ++i) {
        for (typename Pool::elem_t * x = pool.init (); x; x = pool.get ()) {
#ifdef USE_MULTICLASS
          if (_opt.nclass < _lm.nclass ())
            print_err ("set # classes [-C] to a larger value.\n");
#endif
          process_example (*x, false);
        }
        std::fprintf (stderr, "\r%s%s iter=%-3d #ex=%-8ld ",
                      _opt.average ? "Ave " : "", algo[_opt.algo], i, _nex);
        if (_opt.kernel == POLY) std::fprintf (stderr, "#SV=%d; ", _nsv);
        if (std::strcmp (tfn, "")) {
          Model model (_opt);
          if (_opt.kernel == LINEAR)
            model.test (test_pool, _w, _wa, _nex);
          else {
            model.init_weight_trie ();
            model.test (test_pool, _ftrie, _fw, _nex, _sv, _alpha);
          }
        }
      }
    }
    // inherit parameters prior to test
    template <typename Pool>
    void test (Pool &pool, const wv_t &w, const wv_t &wa, const size_t nex) {
      const fl_t n = static_cast <fl_t> (nex + 1);
      wv_t (w).swap (_w); _nf = static_cast <uint> (_w.size () - 1);
      if (_opt.average)
        for (uint fi = 0; fi <= _nf; ++fi) _w[fi] -= wa[fi] / n;
      test (pool);
    }
    template <typename Pool>
    void test (Pool &pool, trie_t** ftrie, wv_t* fw, const size_t nex, const std::vector <sv_t> &sv, const wv_t &alpha) {
      const fl_t n = static_cast <fl_t> (nex + 1);
      std::swap (_ftrie, ftrie); // it's a deal
      for (uint i = PSEUDO_TRIE_N[_opt.d]; i < _nbin; ++i)
        _fw[i].resize (fw[i].size (), _m0);
      for (uint i = 0; i < sv.size (); ++i)
        _pushTo (sv[i], sv[i].getWeight () - (_opt.average ? alpha[i] / n : _m0));
      test (pool);
      std::swap (_ftrie, ftrie);
    }
    template <typename Pool>
    void test (Pool &pool, const uint debug = 0) {
      TIMER (_timer = _timer_pool.push ("classify"));
#ifdef USE_MULTICLASS
      uint corr (0), incorr (0);
#else
      uint pp (0), pn (0), np (0), nn (0);
#endif
#ifdef USE_MULTICLASS
      static w_t m; if (m.size () < _opt.nclass) m.resize (_opt.nclass, 0);
#else
      w_t m = 0;
#endif
      for (typename Pool::elem_t * x = pool.init (); x; x = pool.get ()) {
#ifdef USE_MULTICLASS
        if (_opt.nclass < _lm.nclass ())
          print_err ("set # classes [-C] to a larger value.\n");
#endif
        TIMER (_timer->startTimer ());
        if (_opt.kernel == LINEAR) _getMargin     (m, x->begin (), x->end ());
        else                       _getMarginPoly (m, x->begin (), x->end ());
        TIMER (_timer->stopTimer ());
#ifdef USE_MULTICLASS
        const label_t y_ = std::max_element (&m[0], &m[_opt.nclass]) - &m[0];
        if (y_ == x->getLabel ()) ++corr; else ++incorr;
        if (debug == 2) {
          std::fprintf (stdout, "%s", _lm.get_label (y_));
          for (uint i = 0; i < _opt.nclass; ++i)
            std::fprintf (stdout, " %f", m[i]);
          std::fprintf (stdout, "\n");
        }
#else
        if (m >= 0) if (x->getLabel () > 0) ++pp; else ++np;
        else        if (x->getLabel () > 0) ++pn; else ++nn;
        if (debug == 2)
          std::fprintf (stdout, "%s %f\n", m >= 0 ? "+1" : "-1", m);
#endif
      }
#ifdef USE_MULTICLASS
      std::fprintf (stderr, "acc. %2.3f%% (corr %d) (incorr %d)\n",
                    corr * 100.0 / (corr + incorr), corr, incorr);
#else
      std::fprintf (stderr, "acc. %2.3f%% (pp %d) (pn %d) (np %d) (nn %d)\n",
                    static_cast <fl_t> (pp + nn) * 100.0 / (pp + pn + np + nn),
                    pp, pn, np, nn);
#endif
    }
    void batch () { test_on_file (); }
    void test_on_file (const char * tfn = 0, const uint debug = 0) {
      null_pool <ex_t> pool;
      pool.read (tfn, 0, false, &_lm, 0);
      if (_opt.kernel == POLY) pool.setup (_fm);
      test (pool, debug);
    }
    bool load (const char * mfn) {
      std::fprintf (stderr, "loading..");
      FILE * reader = std::fopen (mfn, "r");
      // examine model type
      if (! reader || std::feof (reader))
        print_err ("cannot read a model: %s\n", mfn);
      char buf[BUF_SIZE]; std::setvbuf (reader, &buf[0], _IOFBF, BUF_SIZE);
      const char flag = static_cast <char> (std::fgetc (reader));
      if (std::fseek (reader, 0, SEEK_SET) != 0) return false;
      char * line = 0;
      size_t read = 0;
      if (flag == 0 || flag == '#') {
#ifdef USE_MULTICLASS
        if (flag == 0)
          print_err ("found a model for single-class; unset USE_MULTICLASS\n");
        if (! getLine (reader, line, read)) return false;
        _lm.read (line, line + read - 1);
        if (_opt.nclass < _lm.nclass ()) init_weight (_lm.nclass ());
#else
        if (flag == '#')
          print_err ("found a model for multi-class; set USE_MULTICLASS\n");
#endif
        _opt.kernel = LINEAR;
        // read a model
        const long offset = std::ftell (reader);
        if (std::fseek (reader, 0, SEEK_END) != 0) return false;
        _nf = static_cast <uint> (std::ftell (reader) - offset) / static_cast <uint> (sizeof (fl_t) * _opt.nclass) - 1;
        if (std::fseek (reader, offset, SEEK_SET) != 0) return false;
        _w.resize (_nf + 1, _m0);
        if (_opt.average) _wa.resize (_nf + 1, _m0);
        size_t nf (0);
        for (size_t i = 0; i <= _nf; ++i) {
#ifdef USE_MULTICLASS
          std::fread (&_w[i][0], sizeof (fl_t), _opt.nclass, reader);
          if (! std::equal (&_w[i][0], &_w[i][_opt.nclass], &_m0[0])) ++nf;
#else
          std::fread (&_w[i], sizeof (fl_t), _opt.nclass, reader);
          if (_w[i]) ++nf;
#endif
        }
        std::fprintf (stderr, "done (# features = %ld).\n", nf);
        std::fclose (reader);
      } else {
        _opt.kernel = POLY;
        if (! getLine (reader, line, read) || std::strncmp (line, "opal", 4))
          return false;
        while (getLine (reader, line, read))
          if (! std::memchr (line, '#', read - 1)) break;
#ifdef USE_MULTICLASS
          else if (char * p = std::strstr (line, "# labels: ")) {
            _lm.read (p, line + read - 1);
            if (_opt.nclass < _lm.nclass ()) init_weight (_lm.nclass ());
          }
#endif
          else if (std::strstr (line, "# kernel parameter -d") != NULL) {
            const uint d = strton <uint> (line, NULL);
            if (_opt.d) {
              if (_opt.d != d)
                print_err ("input kernel_degree [-d] conflicts with %s\n", mfn);
            } else
              _opt.d = d; init_weight_trie ();
          }
        std::fclose (reader);
#ifdef USE_MULTICLASS
        if (_lm.nclass () == 0)
          print_err ("found a model for single-class; unset USE_MULTICLASS\n");
#endif
        mem_pool <sv_t> pool;
        pool.read (mfn, 0UL, true, &_lm, &_fm, '#');
        _fm.build ();
        pool.setup (_fm);
        for (mem_pool <sv_t>::elem_t * x = pool.init (); x; x = pool.get ())
          _pushTo (*x, x->getWeight ());
        size_t ncf (0);
        for (size_t i = 0; i < _w.size (); ++i) {
#ifdef USE_MULTICLASS
          std::fread (&_w[i][0], sizeof (fl_t), _opt.nclass, reader);
          if (! std::equal (&_w[i][0], &_w[i][_opt.nclass], &_m0[0])) ++ncf;
#else
          std::fread (&_w[i], sizeof (fl_t), _opt.nclass, reader);
          if (_w[i]) ++ncf;
#endif
        }
        for (uint i = PSEUDO_TRIE_N[_opt.d]; i < _nbin; ++i)
          ncf += _fw[i].size ();
        std::fprintf (stderr, "done (# explicit features = %ld).\n", ncf);
      }
      return true;
    }
    void average (bool reuse = true) {
      const fl_t n = static_cast <fl_t> (_nex + 1);
      if (_opt.kernel == LINEAR)
        for (uint fi = 0; fi <= _nf; ++fi) _w[fi] -= _wa[fi] / n;
      else
        for (uint i = 0; i < _nsv; ++i)
          if (reuse) _pushTo (_sv[i], - _alpha[i] / n);
          else       _sv[i].weight () -= _alpha[i] / n;
    }
    void save (const char* mfn) { // const
      std::fprintf (stderr, "saving..");
      // write weight vectors
      FILE * writer = std::fopen (mfn, "w");
      if (! writer)
        print_err ("cannot write the model: %s\n", mfn);
      char buf[BUF_SIZE]; std::setvbuf (writer, &buf[0], _IOFBF, BUF_SIZE);
      if (_opt.average) average (false); // averaging
      if (_opt.kernel == LINEAR) { // linear training
#ifdef USE_MULTICLASS
        _lm.write (writer);
        for (size_t i = 0; i <= _nf; ++i)
          std::fwrite (&_w[i][0], sizeof (fl_t), _opt.nclass, writer);
#else
        std::fwrite (&_w[0], sizeof (fl_t), _nf + 1, writer);
#endif
      } else { // kernel; output model in SVM-Light-like format
        std::fprintf (writer, "opal # $Id: pa.h 826 2012-05-10 14:07:00Z ynaga $\n");
#ifdef USE_MULTICLASS
        std::fprintf (writer, "%d ", _lm.nclass ()); _lm.write (writer);
#endif
        std::fprintf (writer, "1 # kernel type\n");
        std::fprintf (writer, "%d # kernel parameter -d\n", _opt.d);
        std::fprintf (writer, "1 # kernel parameter -s\n");
        std::fprintf (writer, "1 # kernel parameter -r\n");
        for (uint i = 0; i < _opt.nclass; ++i) std::fprintf (writer, "0 ");
        std::fprintf (writer, "# threshold b\n"); // for pecco
        // output support vectors w/ alpha
        for (uint i = 0; i < _nsv; ++i) {
          _fm.revertFv2Fv (_sv[i].body ());
          _sv[i].print (writer);
        }
      }
      std::fclose (writer);
      std::fprintf (stderr, "done.\n");
    }
  };
}
#endif /* OPAL_PA_H */
// averaging for retraining; keep averaged vector + nex (reluctant)
// having weights for multi-class in a sparse vector to reduce weights
