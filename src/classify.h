// pecco -- please enjoy classification with conjunctive features
//  $Id: classify.h 827 2012-05-10 14:12:37Z ynaga $
// Copyright (c) 2008-2012 Naoki Yoshinaga <ynaga@tkl.iis.u-tokyo.ac.jp>
#ifndef CLASSIFY_H
#define CLASSIFY_H

#include <sys/stat.h>
#include <getopt.h>
#include <cmath>
#include <cassert>
#include <vector>
#include <sstream>
#include <string>
#include <map>
#include <set>
#include <algorithm>
#include <numeric>
#include "typedef.h"
#ifdef USE_HASH
#include <tr1/unordered_map>
#else
#include <map>
#endif
#include "timer.h"

#define PECCO_COPYRIGHT  "pecco - please enjoy classification w/ conjunctive features\n\
Copyright (c) 2008-2012 Naoki Yoshinaga, All rights reserved.\n\
\n\
Usage: %s [options] model [test]\n\
\n\
model   model file             model\n\
test    test file              test examples; read STDIN if omitted\n\
\n"

#define PECCO_OPT "Optional parameters:\n\
  -t, --classifier=TYPE        select classifier type\n\
                                 0 - PKE | SPLIT\n\
                               * 1 - FST\n\
                                 2 - PKI (kernel model only)\n\
  -e, --event=FILE             examples to obtain feature count for reordering ("")\n\
  -f, --fst-event=FILE         examples to enumerate feature sequence in FST ("")\n\
  -i, --fst-prune-factor=INT   use FST with 2^-INT of feature sequences (0)\n\
  -b, --fst-build-verbose      build FSTs with 2^-i (i > 0) of feature sequences\n\
  -c, --force-compile          force to recompile model\n\
\n"

#define PECCO_OPT_KERNEL "Optional parameters in kernel model:\n\
  -s, --pke-sigma=FLOAT        threshold to feature weight (0)\n\
  -r, --split-ratio=FLOAT      threshold to feature frequency ratio (0)\n\
\n"

#define PECCO_OPT_MISC "Misc.:\n\
  -v, --verbose=INT           verbosity level (1)\n\
  -h, --help                  show this help and exit\n"

static const  char * pecco_short_options = "t:e:f:i:bcs:r:v:h";
static struct option pecco_long_options[] = {
  {"classifier type",   required_argument, NULL, 't'},
  {"event",             required_argument, NULL, 'e'},
  {"fst-event",         required_argument, NULL, 'f'},
  {"fst-prune_ratio",   required_argument, NULL, 'i'},
  {"fst-build-verbose", no_argument,       NULL, 'b'},
  {"force-compile",     no_argument,       NULL, 'c'},
  {"pke-sigma",         required_argument, NULL, 's'},
  {"split-ratio",       required_argument, NULL, 'r'},
  {"verbose",           required_argument, NULL, 'v'},
  {"help",              no_argument,       NULL, 'h'},
  {NULL, 0, NULL, 0}
};

extern char * optarg;
extern int    optind;

namespace pecco {
  enum model_t  { KERNEL, LINEAR };
  enum algo_t   { PKE, FST, PKI };
  enum binary_t { BINARY = true, MULTI = false };
  // options
  template <typename T> T strton (const char * s, char ** error)
  { return static_cast <T> (std::strtol (s, error, 10)); }
  struct option { // option handler
    const char *com, *train, *test, *model;
    std::string   event;
    //
    model_t       type;    // model-type
    algo_t        algo;
    const char *  sigma;
    const char *  fratio;
    size_t        fst_factor;
    bool          fst_verbose;
    bool          force;
    size_t        verbose;
    option () : com ("--"), train (0), test (0), model (0), event (""), type (KERNEL), algo (FST), sigma ("0"), fratio ("0"), fst_factor (0), fst_verbose (false), force (false), verbose (0)  {}
    option (int argc, char ** argv) : com (argc ? argv[0] : "--"), train (0), test (0), model (""), event (""), type (KERNEL), algo (FST), sigma ("0"), fratio ("0"), fst_factor (0), fst_verbose (false), force (false), verbose (0)  {
      set (argc, argv);
    }
    void set (int argc, char ** argv) { // getOpt
      if (argc == 0) return;
      optind = 1;
      double sigma_ (0.0), fratio_ (0.0);
      while (1) {
        int opt = getopt_long (argc, argv,
                               pecco_short_options, pecco_long_options, NULL);
        if (opt == -1) break;
        char * err = NULL;
        switch (opt) {
          case 't': algo        = strton <algo_t>  (optarg, &err); break;
          case 'e': train       = optarg; break;
          case 'f': event       = optarg; break;
          case 'i': fst_factor  = strton <size_t> (optarg, &err);  break;
          case 'b': fst_verbose = true;   break;
          case 'c': force       = true;   break;
#ifdef USE_FLOAT
          case 's': sigma       = optarg; sigma_  = std::strtof (optarg, &err); break;
          case 'r': fratio      = optarg; fratio_ = std::strtof (optarg, &err); break;
#else
          case 's': sigma       = optarg; sigma_  = std::strtod (optarg, &err); break;
          case 'r': fratio      = optarg; fratio_ = std::strtod (optarg, &err); break;
#endif
          case 'v': verbose     = strton <size_t> (optarg, &err);  break;
          case 'h': printCredit (); printHelp (); std::exit (0);   break;
          default:  printCredit (); std::exit (0);
        }
        if (err && *err)
          ny::print_err (HERE "unrecognized option value: %s\n", optarg);
      }
      // errors & warnings
      if (algo != PKE && algo != FST && algo != PKI)
        ny::print_err (HERE "unknown classifier type [-t].\n");
      if (algo == PKI && fratio_)
        ny::print_err (HERE "SPLIT ratio [-r] must be 0 in PKI [-t 2].\n");
      if (std::strcmp (argv[0], "--") == 0) return; // skip
      if (argc < optind + 1) {
        printCredit ();
        ny::print_err ("Type `%s --help' for option details.\n", com);
      }
#ifndef USE_MODEL_SUFFIX
      if (fst_verbose)
        ny::print_err ("[-b] building multiple FSTs are useless since model suffix disabled.\n");
#endif
      model = argv[optind];
      setType ();
      if (type == LINEAR) {
        if (algo == PKI)
          ny::print_err (HERE "PKI [-t 2] is not available for LLM [-t 1].\n");
        if (sigma_)
          std::fprintf (stderr, HERE "PKE sigma [-s] is ignored in LLM.\n");
        if (fratio_)
          std::fprintf (stderr, HERE "SPLIT ratio [-r] is ignored in LLM.\n");
      }
      if (algo == FST && event == "")
        ny::print_err (HERE "FST [-t 1] requires possible examples [-f];\n (you can use the training examples)\n");
      if (++optind != argc) test = argv[optind];
    }
    void setType () {
      FILE * fp = std::fopen (model, "r");
      if (! fp || std::feof (fp))
        ny::print_err (HERE "no model found: %s; train a model first\n", model);
      switch (std::fgetc (fp)) {
        case 'o': // delegate
        case 'T': type = KERNEL;  break;
        default:  type = LINEAR;  break;
      }
      std::fclose (fp);
    }
    void printCredit ()
    { std::fprintf (stderr, PECCO_COPYRIGHT, com); }
    void printHelp ()
    { std::fprintf (stderr, PECCO_OPT PECCO_OPT_KERNEL PECCO_OPT_MISC); }
  };
  // \sum_{i=0}^{k} nCk * num_class is assigned to an array-based pseudo trie
  static const size_t PSEUDO_TRIE_N[] = {0, 21, 11, 8, 6};
  // type alias
  // uchar* -> fl_t  ; conj. feat. -> weight (float)
  typedef ny::TrieKeyBase  <ny::uchar, ny::fl_t> FeatKey;
  typedef ny::TrieKeypLess <ny::uchar, ny::fl_t> FeatKeypLess;
  // uchar* -> fl_t  ; feat. seq.  -> weight (float) or weight ID (int)
  struct FstKey : public FeatKey {
    size_t   weight; // used for node cutoff
    ny::uint count;  // used for node cutoff
    bool     leaf;
    FstKey () : FeatKey (), weight (0), count (0), leaf (false) {}
    FstKey (ny::uchar* k, ny::fl_t* c, size_t l, ny::uint nl)
      : FeatKey (k, c, l, nl), weight (0), count (0), leaf (false) {}
    bool is_prefix (FstKey *a) const { // whether key is a->key's prefix
      if (len > a->len)
        return false;
      else
        for (size_t i = 0; i < len; ++i)
          if (key[i] != a->key[i]) return false;
      return true;
    }
  };
  // feature sequence ordering; see Yoshinaga and Kitsuregawa (EMNLP 2009)
  struct FstKeypLess {
    bool operator () (const FstKey* a, const FstKey* b) const {
      // keep frequent / long keys
      if (a->count * a->weight < b->count * b->weight)
        return true;
      else if (a->count * a->weight > b->count * b->weight)
        return false;
      // keep keys with frequent features
      else return FeatKeypLess () (b, a);
      // return (std::memcmp (a->key, b->key, a->len) > 0);
    }
  };
  // int <-> uint <-> float converter
  union byte_4 {
    int i;
    ny::uint u;
    float f;
    byte_4 (int n)   : i (n) {};
    byte_4 (float b) : f (b) {};
  };
  // build / save da at once
  static inline void build_trie (ny::trie* da,
                                 const std::string name,
                                 const std::string &fn,
                                 std::vector <const char*> &str,
                                 const std::vector <size_t> &len,
                                 const std::vector <ny::trie::result_type> &val,
                                 bool flag, const char* mode = "wb") {
    if (flag) std::fprintf (stderr, " building %s..", name.c_str ());
    if (da->build (str.size (), &str[0], &len[0], &val[0]) != 0 ||
        da->save  (fn.c_str (), mode) != 0 ) {
      ny::print_err (HERE "failed to build %s trie.\n", name.c_str ());
    }
    if (flag) std::fprintf (stderr, "done.\n");
  }
  // check file refreshness
  static inline bool newer (const char * newer, const char * base) {
    struct stat newer_fs, base_fs;
    if (stat (newer, &newer_fs) == -1) return false;
    stat (base,  &base_fs); // need not exist
    return difftime (newer_fs.st_mtime, base_fs.st_mtime) >= 0;
  }
  // bytewise coding (cf. Williams & Zobel, 1999)
  static const size_t KEY_SIZE = 8; // >= log_{2^7} _nf + 1
  class byte_encoder {
  private:
    ny::uchar _key[KEY_SIZE];
    ny::uint  _len;
  public:
    byte_encoder () : _key (), _len (0) {}
    ny::uint encode (ny::uint i, ny::uchar * key) const {
      ny::uint len = 0;
      for (key[len] = (i & 0x7f); i >>= 7; key[++len] = (i & 0x7f))
        key[len] |= 0x80;
      return ++len;
    }
    void encode (ny::uint i) { _len = encode (i, _key); }
    const char *   key    () { return reinterpret_cast <const char *> (&_key[0]); }
    ny::uint       len    () { return _len; }
  };
  template <typename T>
  class ClassifierBase : private ny::Uncopyable {
  protected:
    static const int MAX_KERNEL_DEGREE = 4;
    static const int SORT_MAX = 32;            // 2^SORT_MAX must be > _nf
    static const int BIN_BITS = 6;             // BIN for radix sort
    static const int BIN_ELM  = 1 << BIN_BITS; // MAX for bucket sort
    static const int BIN_MAX  = ((1 << BIN_BITS) - 1);
    // type alias
    typedef ny::map <ny::uint, ny::uint>::type counter_t;
    typedef std::map <char *, ny::uint, ny::pless <char> > lmap; // unordered_map
    // options and classifier parameters
    const option           _opt;
    ny::fv_t               _fv;
    ny::uint               _d;        // degree of feature combination
    ny::uint               _nl;       // # of labels
    // various sizes
    ny::uint               _nf;       // # active feature
    ny::uint               _nf_cut;   // # active feature (pruned)
    ny::uint               _nt;       // # training sample
    ny::uint               _nunit;    // # support vectors / conjunctive features
    // for mapping label to label id
    ny::uint               _tli;      // default target label id
    std::vector <char *>   _li2l;     // label id -> label
    lmap                   _l2li;     // label -> label id
    // for mapping feature numbers to feature indices
    std::vector <ny::uint> _fn2fi;    // feature number -> feature index
    std::vector <ny::uint> _fi2fn;    // feature index  -> feature number
    counter_t              _fncnt;    // feature number -> count
    // double arrays for conjunctive features and fstrie
    ny::trie               _ftrie;    // conjunctive feature -> weight / weight_id
    ny::trie               _fstrie;   // feature sequece     -> weight / weight_id
    // temporary variables
    ny::fl_t *             _fw;       // conjunctive feature id -> weight
    ny::fl_t *             _fsw;      // feature sequence id -> weight
    ny::uint               _f_r;      // min feature index of rare feature (= # common feature)
    ny::uint               _maf;      // max active features per vector;
    // timer
#ifdef USE_TIMER
    ny::TimerPool          _timer_pool;
    ny::Timer *            _enc_t;    // feature mapping
    ny::Timer *            _bound_t;    // fst classify
    ny::Timer *            _model_t;  // compiling/loading  model
    ny::Timer *            _pke_t;    // pke classify
    ny::Timer *            _fst_t;    // fst classify
#endif
    // profiling
#ifdef USE_PROFILING
    size_t                 _all;         // all       (for fst)
    size_t                 _hit;         // hit ratio (for fst)
    size_t                 _lookup;      // lookup    (for pke)
    size_t                 _flen;        // traverse  (for pke)
    size_t                 _traverse;    // traverse  (for pke)
    size_t                 _lookup_sp;   // lookup    (for pke)
    size_t                 _traverse_sp; // traverse  (for pke)
#endif
    T *       _derived ()             { return static_cast <T*>       (this); }
    const T * _derived_const () const { return static_cast <const T*> (this); }
    size_t _cost_fun (size_t m, size_t n);
    // bytewise coding for feature indexes (cf. Williams & Zobel, 1999)
    // various functions for sorting feature (index) vectors
    template <typename Iterator>
    static void _insertion_sort (const Iterator &first, const Iterator &last);
    template <typename value_type, typename SrcIterator, typename OutIterator>
    static void _radix (SrcIterator first, SrcIterator last, OutIterator dest, size_t shift);
    template <size_t SIZE, typename Iterator>
    static void _radix_sort (const Iterator &first, const Iterator &last);
    template <typename Iterator>
    static void _bucket_sort (const Iterator &first, const Iterator &last);
    // feature mapping
    void _convertFvstr2Fv (char * &p, const char * const p_end,
                           ny::fv_t &fv, const bool flag);
    void _convertFv2Fv (ny::fv_t &fnv) const;
    void _sortFv (ny::fv_t &fv) const;
    bool _packingFeatures (std::vector <ny::fv_t> &fvv);
    // feature sequence trie builder
    bool _setFStrie ();
    // classification functions
    template <binary_t FLAG>
    void addScore (double * score, const ny::uint pos) const
    { _derived_const ()->addScore <FLAG> (score, pos); }
    template <binary_t FLAG>
    void addScore (double * score, const int n, const ny::fl_t * const w) const
    { _derived_const ()->addScore <FLAG> (score, n, w); }
    template <int D, binary_t FLAG>
    void _pkePseudoInnerLoop (double * score, ny::fv_it it, const ny::fv_it &beg, const ny::fv_it &end, const ny::uint pos);
    template <int D, binary_t FLAG>
    void _pkeInnerLoop (double * score, ny::fv_it it, const ny::fv_it &beg, const ny::fv_it &end, const size_t pos);
    template <int D>
    void _pkeClassify (double * score, ny::fv_it it, const ny::fv_it &beg, const ny::fv_it &end);
    void _pkeClassify (const ny::fv_t &fv, double * score, ny::fv_it it, const ny::fv_it &end);
    void _pkeClassify (const ny::fv_t &fv, double * score)
    { return _pkeClassify (fv, score, fv.begin (), fv.end ()); }
    template <binary_t FLAG>
    void _fstClassify (double * score, ny::fv_it &cit, const ny::fv_it &end);
    void _fstClassify (const ny::fv_t &fv, double * score);
    void _baseClassify (const ny::fv_t &fv, double * score, ny::fv_it it, const ny::fv_it &end)
    { _derived ()->baseClassify (fv, score, it, end); };
    bool _setOpt (const char* opt_str);
  public:
    ClassifierBase (const pecco::option &opt) :
      _opt (opt), _fv (), _d (0), _nl (0), _nf (0), _nf_cut (0), _nt (0), _nunit (0),  _tli (0), _li2l (), _l2li (), _fn2fi (), _fi2fn (), _fncnt (), _ftrie (), _fstrie (), _fw (0), _fsw (0), _f_r (0), _maf (0)
#ifdef USE_TIMER
      , _timer_pool ((std::string ("pecco profiler (") + _opt.model + ")").c_str ()), _enc_t (_timer_pool.push ("enc")), _model_t (_timer_pool.push ("model")), _pke_t (_timer_pool.push ("pke", "classify")), _fst_t (_timer_pool.push ("fst", "classify"))
#endif
#ifdef USE_PROFILING
      , _all (0), _hit (0), _lookup (0), _flen (0), _traverse (0), _lookup_sp (0), _traverse_sp (0)
#endif
    {}
    // interface
    bool load (const char * model)
    { return _derived ()->load (model); }
    void classify (ny::fv_t &fnv, double * score)
    { _derived ()->classify (fnv, score); }
    bool abuse_trie () const
    { return _derived_const ()->abuse_trie (); }
    bool is_binary_classification () const
    { return _derived_const ()->is_binary_classification (); }
    void printScore (const ny::uint li, const double * score)
    { _derived ()->printScore (li, score); }
    ny::uint getLabel (const double * score)
    { return _derived ()->getLabel (score); }
    void classify  (char * p, double * score);
    void batch     ();
    void printStat ();
  protected:
    ~ClassifierBase () {};
  };
}
#endif /* CLASSIFY_H */
