// J.DepP -- Japanese Dependency Parsers
//  $Id: pdep.h 845 2012-05-17 16:49:17Z ynaga $
// Copyright (c) 2008-2012 Naoki Yoshinaga <ynaga@tkl.iis.u-tokyo.ac.jp>
#ifndef PDEP_H
#define PDEP_H

#include <unistd.h>
#include <fcntl.h>
#include <cstdio>
#include <cmath>
#include <set>
#include <map>
#include <stack>
#include <list>

#define MAX_FILENAME_SIZE 500

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef USE_HASH
#include <tr1/unordered_map>
#endif
#include "typedef.h"
// trie type
#if defined (USE_DARTS) || defined (USE_DARTS_CLONE)
#include <darts.h>
#endif
#ifdef USE_CEDAR
#include "cedar.h"
#endif
#ifdef USE_TIMER
#include "timer.h"
#endif
#if defined (USE_OPAL) || defined (USE_SVM)
#include "kernel.h"
#endif
#ifdef USE_SVM
#include <tinysvm.h>
#endif
#ifdef USE_OPAL
#define USE_PRUNE
#define USE_PMT
#include "pa.h"
#endif
#ifdef USE_MAXENT
#include "maxent.h"
#include "linear.h"
#endif
#include "pecco.h"
#ifdef USE_AS_STANDALONE
#include <mecab.h>
#endif
#define JDEPP_COPYRIGHT  "J.DepP - Japanese Dependency Parser\n\
Copyright (c) 2008-2012 Naoki Yoshinaga\n\
\n\
Usage: %s [options] -- [learner options] -- [chunker classifier options] -- [parser classifier options] < test\n\
\n\
test    test file\n\
\n"

#define JDEPP_OPT0 "Optional parameters in training / testing:\n\
  -t, --type=TYPE             select running mode of J.DepP\n\
                                0 - learn\n\
                              * 1 - parse\n\
                                2 - both\n\
                                3 - cache\n\
  -e, --encoding=TYPE         select encoding of input\n\
                              * 0 - UTF-8\n\
                                1 - EUC-JP\n\
  -c, --corpus=FILE           training corpus in JDEPP format ('train.JDP')\n\
  -m, --model-dir=DIR         model directory ('" JDEPP_DEFAULT_MODEL "')\n\
  -p, --parser=TYPE           select parsing algorithm\n\
                              * 0 - shift reduce\n\
                                1 - cascaded chunking\n\
                                2 - backward\n\
                                3 - tournament\n"

#ifdef USE_AS_STANDALONE
#define JDEPP_OPT1 "  -I, --input-format=TYPE     select type of input format\n\
                              * 0 - RAW sentences\n\
                                1 - + POS / BUNSETSU annotation\n\
                                2 - + DEPENDENCY annotation\n\
\n"
#else
#define JDEPP_OPT1 "  -I, --input-format=TYPE     select type of input format\n\
                              * 0 - POS-tagged sentences\n\
                                1 - + BUNSETSU annotation\n\
                                2 - + DEPENDENCY annotation\n\
\n"
#endif

#if defined (USE_OPAL)
#define OPAL_STR "                              * 0 - OPAL\n"
#else
#define OPAL_STR "                                0 - OPAL   (disabled)\n"
#endif
#if defined (USE_SVM)
#if ! defined (USE_OPAL)
#define SVM_STR "                              * 1 - SVM\n"
#else
#define SVM_STR "                                1 - SVM\n"
#endif
#else
#define SVM_STR "                                1 - SVM    (disabled)\n"
#endif
#if defined (USE_MAXENT)
#if ! defined (USE_OPAL) && ! defined (USE_SVM)
#define MAXENT_STR "                              * 2 - MaxEnt\n"
#else
#define MAXENT_STR "                                2 - MaxEnt\n"
#endif
#else
#define MAXENT_STR "                                2 - MaxEnt (disabled)\n"
#endif
#define JDEPP_OPT_TRAIN "Optional parameters in training:\n\
  -l, --learner=TYPE          select type of learning library\n" OPAL_STR SVM_STR MAXENT_STR "\
  -n, --max-sent=INT          max. # processing sentences (0: all)\n\
\n"

#ifdef USE_AS_STANDALONE
#define JDEPP_OPT_TEST "Optional parameters in testing:\n\
  -d, --mecab-dic=DIR         use MeCab dictionary ('" MECAB_DICT "')\n\
\n"
#else
#define JDEPP_OPT_TEST
#endif

#define JDEPP_OPT_MISC "Misc.:\n\
  -v, --verbose=INT           verbosity level (1)\n\
  -h, --help                  show this help and exit\n"

#define MAXENT_OPT "\nOptions for MaxEnt learners are as follows:\n\
  -d, --degree=INT            maximum degree of feature combinations (<=3)\n\
  -l, --algorithm=INT         select type of optimization algorithm\n\
                              * 0 - SGD-L1\n\
                                1 - OWLQN-L1\n\
                                2 - LBFGS-L2\n\
  -c, --reg-cost=INT          cost of regularization\n\
\n"

#ifdef USE_AS_STANDALONE
static const  char * jdepp_short_options = "t:e:c:m:p:I:b:l:n:d:x:v:h:F";
#else
static const  char * jdepp_short_options = "t:e:c:m:p:I:b:l:n:x:v:h:F";
#endif
static struct option jdepp_long_options[] = {
  {"type",         required_argument, NULL, 't'},
  {"encoding",     required_argument, NULL, 'e'},
  {"corpus",       required_argument, NULL, 'c'},
  {"model-dir",    required_argument, NULL, 'm'},
  {"parser",       required_argument, NULL, 'p'},
  {"input-format", required_argument, NULL, 'I'},
  {"cluster-bits", required_argument, NULL, 'b'},
  {"learner",      required_argument, NULL, 'l'},
  {"max-sent",     required_argument, NULL, 'n'},
#ifdef USE_AS_STANDALONE
  {"mecab-dic",    required_argument, NULL, 'd'},
#endif
  {"xcode",        required_argument, NULL, 'x'},
  {"verbose",      required_argument, NULL, 'v'},
  {"help",         no_argument,       NULL, 'h'},
  {"filelist",     no_argument,       NULL, 'F'},
  {NULL, 0, NULL, 0}
};

static const  char * maxent_short_options = "d:l:c:h";
static struct option maxent_long_options[] = {
  {"degree",       required_argument, NULL, 'd'},
  {"algorithm",    required_argument, NULL, 'l'},
  {"reg-cost",     required_argument, NULL, 'c'},
  {"help",         no_argument,       NULL, 'h'},
  {NULL, 0, NULL, 0}
};

extern char * optarg;
extern int    optind;

#define UTF8_COMMA         "\xE8\xAA\xAD\xE7\x82\xB9" // 読点
#define UTF8_PERIOD        "\xE5\x8F\xA5\xE7\x82\xB9" // 句点
#define UTF8_POST_PARTICLE "\xE5\x8A\xA9\xE8\xA9\x9E" // 助詞

#if defined (USE_JUMAN_POS)
#define UTF8_BRACKET_START "\xE6\x8B\xAC\xE5\xBC\xA7\xE5\xA7\x8B" // 括弧始
#define UTF8_BRACKET_END   "\xE6\x8B\xAC\xE5\xBC\xA7\xE7\xB5\x82" // 括弧終
#define UTF8_SPECIAL       "\xE7\x89\xB9\xE6\xAE\x8A" // 特殊
#define UTF8_SUFFIX        "\xE6\x8E\xA5\xE5\xB0\xBE\xE8\xBE\x9E" // 接尾辞
#elif defined (USE_IPA_POS)
#define UTF8_BRACKET_START "\xE6\x8B\xAC\xE5\xBC\xA7\xE9\x96\x8B" // 括弧開
#define UTF8_BRACKET_END   "\xE6\x8B\xAC\xE5\xBC\xA7\xE9\x96\x89" // 括弧閉
#define UTF8_SPECIAL       "\xE8\xA8\x98\xE5\x8F\xB7" // 記号
#endif

#define EUC_COMMA         "\xC6\xC9\xC5\xC0"
#define EUC_PERIOD        "\xB6\xE7\xC5\xC0"
#define EUC_POST_PARTICLE "\xBD\xF5\xBB\xEC"

#if defined (USE_JUMAN_POS)
#define EUC_BRACKET_START "\xB3\xE7\xB8\xCC\xBB\xCF"
#define EUC_BRACKET_END   "\xB3\xE7\xB8\xCC\xBD\xAA"
#define EUC_SUFFIX        "\xC0\xDC\xC8\xF8\xBC\xAD"
#define EUC_SPECIAL       "\xC6\xC3\xBC\xEC"
#elif defined (USE_IPA_POS)
#define EUC_BRACKET_START "\xB3\xE7\xB8\xCC\xB3\xAB"
#define EUC_BRACKET_END   "\xB3\xE7\xB8\xCC\xCA\xC4"
#define EUC_SPECIAL       "\xB5\xAD\xB9\xE6"
#endif

namespace pdep {
  // type alias
  typedef std::vector <uint64_t> flag_t;
  static const ny::uint FLAG_LEN = sizeof (flag_t::value_type) * 8;
#if   defined (USE_JUMAN_POS)
#if   defined (USE_MECAB)
  enum field_t {SURF, POS1, POS2, TYPE, INFL, FIN, YOMI, NUM_FIELD}; // OTHER, C1, C2, C3, C4, C5, C6, C7, C8, C9, C10, C11, C12, C13, C14, C15, C16, C17, C18, NUM_FIELD};
#elif defined (USE_JUMAN)
  enum field_t {SURF, YOMI, FIN, POS1, POS_ID1, POS2, POS_ID2, TYPE, TYPE_ID, INFL, INFL_ID, NUM_FIELD};
#endif
#elif defined (USE_IPA_POS)
  enum field_t {SURF, POS1, POS2, POS3, POS4, TYPE, INFL, FIN, YOMI, PRON, NUM_FIELD};
#endif
  // static const variables
  template <typename T>
  static void widen (T * &array, const ny::uint &avail, const ny::uint &filled = 0) {
    T * tmp = static_cast <T*> (::operator new (avail * sizeof (T)));
    if (filled) {
      std::memcpy (&tmp[0], &array[0], sizeof (T) * filled);   // delegate
      for (ny::uint i = filled; i < avail; ++i) new (&tmp[i]) T (); // fill
    }
    std::swap (tmp, array);
    if (tmp) ::operator delete (tmp); // don't destruct delegated resource
  }
  template <typename T>
  T strton (const char * s, char ** error)
  { return static_cast <T> (std::strtol (s, error, 10)); }
  //
  enum process_t     { LEARN, PARSE, BOTH, CACHE };
  enum parser_t      { LINEAR, CHUNKING, BACKWARD, TOURNAMENT };
  enum learner_t     { OPAL, SVM, MAXENT };
  enum input_t       { RAW, CHUNK, DEPND };
  enum maxent_algo_t { SGD, OWLQN, LBFGS }; // MaxEnt optimizer
  // The command-line arguments will override the following default parameters
  // default learner parameters
  class option { // option handler
  public:
    const char *com, *train;
    std::string model_dir;
    //
    process_t mode;
    parser_t  parser;
    bool      utf8;
    ny::uint  cbits;
    ny::uint  clen;
    bool      batchmode;
    learner_t learner;
    ny::uint  max_sent;
    input_t   input;
    ny::uint  xcode;
#ifdef USE_AS_STANDALONE
    const char *    mecab_dic;
#endif
    int       verbose;
    int       learner_argc;
    char **   learner_argv;
    int       depnd_argc;
    char **   depnd_argv;
    int       chunk_argc;
    char **   chunk_argv;
    //
    option (int argc, char ** argv) :
      com (argv[0]), train ("train.JDP"),
#ifdef USE_STACKING
      model_dir (JDEPP_DEFAULT_MODEL "_stack"),
#else
      model_dir (JDEPP_DEFAULT_MODEL),
#endif
      mode (PARSE), parser (LINEAR), utf8 (true), cbits (0), clen (0), batchmode (false),
#if   defined (USE_OPAL)
      learner (OPAL),
#elif defined (USE_SVM)
      learner (SVM),
#elif defined (USE_MAXENT)
      learner (MAXENT),
#endif
      max_sent (0), input (RAW), xcode (0),
#ifdef USE_AS_STANDALONE
      mecab_dic (MECAB_DICT),
#endif
      verbose (0), learner_argc (0), learner_argv (0), depnd_argc (0), depnd_argv (0), chunk_argc (0), chunk_argv (0) {
      // getOpt
      if (argc == 0) return;
      optind = 1;
      while (1) {
        int opt = getopt_long (argc, argv,
                               jdepp_short_options, jdepp_long_options, NULL);
        if (opt == -1) break;
        char * err = NULL;
        switch (opt) {
          case 't': mode      = strton <process_t> (optarg, &err); break;
          case 'e': utf8      = std::strtol (optarg, &err, 10) == 0; break;
          case 'c': train     = optarg; break;
          case 'm': model_dir = optarg; break;
          case 'p': parser    = strton <parser_t>  (optarg, &err); break;
          case 'I': input     = strton <input_t>   (optarg, &err); break;
          case 'F': batchmode = true; break;
          case 'b':
            do {
              const ny::uint depth = strton <ny::uint> (optarg, &optarg);
              cbits |= 1 << (depth - 1);
              clen = std::max (clen, depth);
            } while (*optarg++ != '\0');
            break;
            // chech bits
            // training parameters
          case 'l': learner   = strton <learner_t> (optarg, &err); break;
          case 'n': max_sent  = strton <ny::uint> (optarg, &err);  break;
            // misc
#ifdef USE_AS_STANDALONE
          case 'd': mecab_dic = optarg; break;
#endif
          case 'x': xcode     = strton <ny::uint> (optarg, &err); break;
          case 'v': verbose   = strton <int> (optarg, &err); break;
          case 'h': printCredit (); printHelp (); std::exit (0);  break;
          default:  printCredit (); std::exit (0);
        }
        if (err && *err)
          ny::print_err (HERE "unrecognized option value: %s\n", optarg);
      }
      // std::fprintf (stderr, "xcode:");
      // for (ny::uint i (0); i < 8; ++i)
      //   std::fprintf (stderr, " %c", ((xcode >> i) & 0x1) ? '+' : '-');
      // std::fprintf (stderr, "\n");
      // errors & warnings
      if (learner != OPAL && learner != SVM && learner != MAXENT)
        ny::print_err (HERE "unknown learner [-l].\n");
      if (mode != LEARN && mode != PARSE && mode != BOTH && mode != CACHE)
        ny::print_err (HERE "unknown running mode [-t].\n");
      if (parser != LINEAR && parser != CHUNKING && parser != BACKWARD &&
          parser != TOURNAMENT)
        ny::print_err (HERE "unknown parsing algorithm [-p].\n");
      if (input != RAW && input != CHUNK && input != DEPND)
        ny::print_err (HERE "unknown input format [-I].\n");
      struct stat st;
      if (stat (model_dir.c_str (), &st) != 0)
        ny::print_err (HERE "no such directory: %s [-m]\n", model_dir.c_str ());
#ifdef USE_AS_STANDALONE
      if (input == RAW && mecab_dic && stat (mecab_dic, &st) != 0)
          ny::print_err (HERE "no such file or directory: %s\n", mecab_dic);
#endif
      if (input == CHUNK && parser != LINEAR)
        std::fprintf (stderr, HERE "parsing algorithm [-p] is ignored in training a chunker.\n");
      // learner options
      if (std::strcmp (argv[optind - 1], "--") == 0) --optind;
      _set_library_options (optind, argc, argv, learner_argc, learner_argv);
      // classifier options for bunsetsu chunker
      _set_library_options (optind, argc, argv, chunk_argc, chunk_argv);
      // classifier options for dependency parser
      _set_library_options (optind, argc, argv, depnd_argc, depnd_argv);
    }
    void printCredit () { std::fprintf (stderr, JDEPP_COPYRIGHT, com); }
    void printHelp   () { std::fprintf (stderr, JDEPP_OPT0 JDEPP_OPT1 JDEPP_OPT_TRAIN JDEPP_OPT_TEST JDEPP_OPT_MISC); }
  private:
    void _set_library_options (int &i, const int argc, char ** argv,
                               int &libargc, char ** &libargv) {
      if (i < argc) {
        if (std::strcmp (argv[optind], "--") == 0) { // library options
          libargv = &argv[optind];
          libargc = 1;
          while (optind + libargc < argc &&
                 std::strcmp (libargv[libargc], "--") != 0)
            ++libargc;
          i += libargc;
        } else {
          printCredit ();
          ny::print_err ("Type `%s --help' for option details.\n", com);
        }
      }
    }
  };
  struct maxent_option { // option handler
    maxent_algo_t algo;
    ny::uint      degree;
    double        reg_cost;
    //
    maxent_option () : algo (SGD), degree (2), reg_cost (1.0) {}
    maxent_option (int argc, char ** argv) : algo (SGD), degree (2), reg_cost (1.0) {
      set (argc, argv);
    }
    void set (int argc, char ** argv) { // getOpt
      if (argc == 0) return;
      optind = 1;
      while (1) {
        int opt = getopt_long (argc, argv,
                               maxent_short_options, maxent_long_options, NULL);
        if (opt == -1) break;
        char * err = NULL;
        switch (opt) {
          case 'd': degree   = strton <ny::uint> (optarg, &err);      break;
          case 'l': algo     = strton <maxent_algo_t> (optarg, &err); break;
          case 'c': reg_cost = std::strtod (optarg, &err); break;
          case 'h': printHelp (); std::exit (0); break;
          default:  std::exit (0);
        }
        if (err && *err)
          ny::print_err ("unrecognized option value: %s\n", optarg);
      }
      // errors
      if (algo != SGD && algo != OWLQN && algo != LBFGS)
        ny::print_err (HERE "unknown optimization algorithm [-l]");
      if (degree == 0 || degree >= 4)
        ny::print_err (HERE "set degree <= 3 for [-d].\n");

    }
    void printHelp () { std::fprintf (stderr, MAXENT_OPT); }
  };
  // dictionary utility variables
  typedef std::map <const char*, ny::uint, ny::pless <char> > sbag_t;
  class dict_base_t : private ny::Uncopyable {
  private:
    char *    _data_ptr;
    ny::trie  _data;
    ny::uint  _num_lexical_features;
    ny::uint  _num_context_features;
  protected:
    dict_base_t (const char *fn) : _data_ptr (0), _data (), _num_lexical_features (0), _num_context_features (0) {
      FILE * fp = std::fopen (fn, "rb");
      if (std::fread (&_num_lexical_features, sizeof (ny::uint), 1, fp) != 1 ||
          std::fread (&_num_context_features, sizeof (ny::uint), 1, fp) != 1)
        ny::print_err (HERE "broken dic: delete %s\n", fn);
      ++_num_lexical_features; // reserve room for unseen feature
      long offset = std::ftell (fp);
      std::fseek (fp, 0,      SEEK_END);
      size_t size   = static_cast <size_t> (std::ftell (fp) - offset);
      _data_ptr = new char[size];
      std::fseek (fp, offset, SEEK_SET);
      std::fread (_data_ptr, sizeof (char), size, fp);
      _data.set_array (_data_ptr, size / _data.unit_size ());
      std::fclose (fp);
    }
    ~dict_base_t () { delete [] _data_ptr; }
  public:
    ny::uint lookup (const char *key, size_t len = 0) const {
      int n = _data.exactMatchSearch <ny::trie::result_type> (key, len);
      return n >= 0 ? static_cast <ny::uint> (n) : _num_lexical_features - 1;
    }
    ny::uint context_feature_bit_len () const
    { return (_num_context_features - 1) / FLAG_LEN + 1; }
    bool     is_context_feature (const ny::uint id) const
    { return id < _num_context_features; }
    ny::uint num_lexical_features () const { return _num_lexical_features; }
  };
  class dict_t : public dict_base_t {
  public: // surface id aliases
    const ny::uint  comma;
    const ny::uint  period;
    const ny::uint  particle;
    const ny::uint  bracket_start;
    const ny::uint  bracket_end;
    const ny::uint  special;
#ifdef USE_JUMAN_POS
    const ny::uint  suffix;
#endif
    dict_t (const char *fn, bool utf8 = true) :
      dict_base_t (fn), comma (lookup (utf8 ? UTF8_COMMA : EUC_COMMA)), period (lookup (utf8 ? UTF8_PERIOD : EUC_PERIOD)), particle (lookup (utf8 ? UTF8_POST_PARTICLE : EUC_POST_PARTICLE)), bracket_start (lookup (utf8 ? UTF8_BRACKET_START : EUC_BRACKET_START)), bracket_end (lookup (utf8 ? UTF8_BRACKET_END : EUC_BRACKET_END)), special (lookup (utf8 ? UTF8_SPECIAL : EUC_SPECIAL))
#ifdef USE_JUMAN_POS
      , suffix (lookup (utf8 ? UTF8_SUFFIX : EUC_SUFFIX))
#endif
    {}
  };
  class Morph {
  private:
    ny::uint _field[NUM_FIELD];
  public:
    size_t       length;
    const char * surface;
    const char * feature;
    bool         chunk;
    bool         chunk_ref;
    double       chunk_prob;
    Morph () : _field (), length (0), surface (0), feature (0), chunk (false), chunk_ref (false), chunk_prob (0) {}
    Morph (const Morph& m) : _field (), length (m.length), surface (m.surface), feature (m.feature), chunk (m.chunk), chunk_ref (m.chunk_ref), chunk_prob (m.chunk_prob)
    { std::copy (&m._field[0], &m._field[0] + NUM_FIELD, &_field[0]); }
    ~Morph () {}
    void set (char * p, const size_t len, const dict_t * dict, bool flag) {
      surface = p;
      char * const p_end = p + len;
      *p_end = '\0'; while (p != p_end && *p != SURFACE_END) ++p; // juman
      length = static_cast <size_t> (p - surface); feature = ++p;
      chunk_ref = flag;
      set (dict);
    }
    void set (const size_t length_,  const char * surface_, 
              const char * feature_, const dict_t * dict) {
      length  = length_; surface = surface_; feature = feature_;
      set (dict);
    }
    void set (const dict_t * dict) {
      _field[SURF] = dict->lookup (surface, length); // read surface
      // read feature
      ny::uint i = 1;
      for (const char * p (feature), * f (p); i < NUM_FIELD; f = ++p, ++i) {
        while (*p != '\0' && *p != FEATURE_SEP) ++p;
        if (i == POS1 || i == POS2 || i == INFL)
          _field[i] = dict->lookup (f, static_cast <size_t> (p - f));
      }
#ifdef USE_JUMAN_POS
      if (i < NUM_FIELD) {
        std::fwrite (surface, sizeof (char), length, stderr);
        ny::print_err (HERE "# fields, %d, is less than %d.\n", i, NUM_FIELD);
      }
#endif
    }
    ny::uint surf () const { return _field[SURF]; }
    ny::uint pos1 () const { return _field[POS1]; }
    ny::uint pos2 () const { return _field[POS2]; }
    ny::uint infl () const { return _field[INFL]; }
    ny::uint cluster (ny::uint i) const { return _field[NUM_FIELD + i]; }
  };
  class Bunsetsu {
  private:
    void _set_fbits (ny::uint i) {
      fbits[i / FLAG_LEN] |= (static_cast <flag_t::value_type> (1) << (i % FLAG_LEN));
    }
  public:
    ny::uint mpos;     // start position
    ny::uint mlen;     // number of morphs
    flag_t   fbits;
    int      head;     // head_morph offset
    int      tail;     // tail
    ny::uint id;       // bunsetsu id
    int      did;      // dest id (gold)
    int      est_did;  // dest id (estimated)
    int      cand;
    char     dt;
    char     candt;
    int      bracket;  // has bracket
    bool     comma;    // has comma
    bool     period;   // has period
    double   depnd_prob;
    // char dtype;
    Bunsetsu () : mpos (0), mlen (0), fbits (), head (-1), tail (-1), id (0), did (-1), est_did (-1), cand (-1), dt ('D'), candt ('D'), bracket (0), comma (false), period (false), depnd_prob (0) {}
    Bunsetsu (const Bunsetsu &b) : mpos (b.mpos), mlen (b.mlen), fbits (b.fbits), head (b.head), tail (b.tail), id (b.id), did (b.did), est_did (b.est_did), cand (b.cand), dt (b.dt), candt (b.candt), bracket (b.bracket), comma (b.comma), period (b.period), depnd_prob (b.depnd_prob) {}
    ~Bunsetsu () {}
    void clear () { // leave other members be overridden
      std::fill (fbits.begin (), fbits.end (), 0);
      head = tail = est_did = cand = -1; dt = 'D'; candt = 'D';
      mpos = mlen = 0; bracket = 0;
      comma = period = false; depnd_prob = 0;
    }
#ifdef USE_STACKING 
    void set (char * p, const size_t len) { // ex. '* 1 7D'
      const char * p_end (p + len);
      id  = ny::strton <ny::uint> (p + 2, &p);
      did = ny::strton <int> (++p, &p);
      if (++p != p_end) {
        cand  = ny::strton <int> (++p, &p); // ignored in testing
        candt = *p;
      }
#else
    void set (char * p, const size_t) { // ex. '* 1 7D'
      id  = ny::strton <ny::uint> (p + 2, &p); ++p;
      did = ny::strton <int> (p, &p);
      dt  = *p;
#endif
    }
    void set (const ny::uint id_) {
      id  = id_;
      did = -1;
#ifdef USE_STACKING
      cand = -1;
#endif
    }
    bool setup (const dict_t * dict, const Morph * morph) {
      fbits.resize (dict->context_feature_bit_len (), 0);
      if (mlen < 1) return false;
      for (ny::uint i = mlen - 1;; --i) {
        const Morph &m = morph[i];
        if (m.pos1 () == dict->special) {
          if      (m.pos2 () == dict->comma)  comma  |= true; // i == mlen - 1;
          else if (m.pos2 () == dict->period) period |= true; // i == mlen - 1;
          else {
            if      (m.pos2 () == dict->bracket_end)   --bracket;
            else if (m.pos2 () == dict->bracket_start) ++bracket;
          }
        } else {
          if (tail == -1) tail = static_cast <int> (i);
          if (m.pos1 () == dict->particle) { // || m.pos1 () != dict->suffix) {
            if (dict->is_context_feature (m.surf ())) // adding suffix may capture agreement
              _set_fbits (m.surf ());
          } else
            if (head == -1) head = static_cast <int> (i);
        }
        if (i == 0) break;
      }
      if (did < 0 || did > static_cast <int> (id))  return true;
#ifndef NDEBUG
      std::fprintf (stderr, "\tbroken dependency?: %d->%d\n", id, did);
#endif
      return false;
    }
  };
  class Sentence : private ny::Uncopyable {
  private:
    ny::uint   _bavail;
    ny::uint   _mavail;
    char       _res[IOBUF_SIZE]; // save output
    char*      _ptr;             // current position in result buffer
  public:
    FILE*      out_fp;
    Bunsetsu * bun;
    Morph    * morph;
    ny::uint   blen;
    ny::uint   mlen;
    Sentence () : _bavail (1), _mavail (1), _res (), _ptr (&_res[0]), bun (static_cast <Bunsetsu *> (::operator new (_bavail * sizeof (Bunsetsu)))), morph (static_cast <Morph *> (::operator new (_mavail * sizeof (Morph)))), blen (0), mlen (0) {
      for (ny::uint i = 0; i < _bavail; ++i) new (&bun[i])   Bunsetsu ();
      for (ny::uint i = 0; i < _mavail; ++i) new (&morph[i]) Morph    ();
    }
    ~Sentence () {
      for (ny::uint i = 0; i < _bavail; ++i) bun[i].~Bunsetsu ();
      ::operator delete (bun);
      for (ny::uint i = 0; i < _mavail; ++i) morph[i].~Morph ();
      ::operator delete (morph);
    }
    void clear () {
      if (blen) do { bun[--blen].clear (); } while (blen);
      mlen = 0;
    };
    void setHeader (char * cs, const size_t len) // skip comment
    { if (_ptr == &_res[0]) { std::memcpy (_ptr, cs, len); _ptr += len; } }
    void setup (const dict_t *dict) { // complete information
      for (ny::uint i = 0; i < blen; ++i) {
        bun[i].mlen = (i == blen - 1 ? mlen : bun[i + 1].mpos) - bun[i].mpos;
        bun[i].setup (dict, &morph[bun[i].mpos]);
      }
    }
    void addBunsetsu (char * cs, const size_t len, const ny::uint mpos) {
      if (blen == _bavail) { _bavail <<= 1; widen (bun, _bavail, blen); }
      bun[blen].mpos = mpos;
      bun[blen].set (cs, len);
      ++blen;
    }
    void addBunsetsu (const ny::uint mpos) {
      if (blen == _bavail) { _bavail <<= 1; widen (bun, _bavail, blen); }
      bun[blen].mpos = mpos;
      bun[blen].set (blen);
      ++blen;
    }
    void addMorph (char * cs, const size_t len, const dict_t * dict,
                   bool flag = false) {
      if (mlen == _mavail) { _mavail <<= 1; widen (morph, _mavail, mlen); }
      morph[mlen].set (cs, len, dict, flag);
      ++mlen;
    }
    void addMorph (const size_t length, const char * surface,
                   const char * feature, const dict_t * dict) {
      if (mlen == _mavail) { _mavail <<= 1; widen (morph, _mavail, mlen); }
      morph[mlen].set (length, surface, feature, dict);
      ++mlen;
    }
    void print (const input_t in) {
      if (mlen) {
        for (ny::uint i = 0; i < blen; ++i) {
          const Bunsetsu &b = bun[i];
          if (in == DEPND)
#ifdef USE_STACKING
            _ptr += std::sprintf (_ptr, "* %d %dD@%f %d%c #%d%c\n", i, b.est_did, b.depnd_prob, b.did, b.dt, b.cand, b.candt);
#else
            _ptr += std::sprintf (_ptr, "* %d %dD@%f %d%c\n", i, b.est_did, b.depnd_prob, b.did, b.dt);
#endif
          else
#ifdef USE_STACKING
            _ptr += std::sprintf (_ptr, "* %d %dD #%d%c\n", i, b.est_did, b.cand, b.candt);
#else
            _ptr += std::sprintf (_ptr, "* %d %dD\n", i, b.est_did);
#endif
          const Morph  * ms = &morph[b.mpos];
          const ny::uint mlen = b.mlen;
          for (ny::uint j = 0; j < mlen; ++j) {
            const Morph &m = ms[j];
            std::memcpy (_ptr, m.surface, m.length); _ptr += m.length;
            if (in == CHUNK)
              _ptr += std::sprintf (_ptr, "%c%s\t%c@%f %c\n",
                                    SURFACE_END, m.feature, m.chunk ? 'B' : 'I',
                                    m.chunk_prob, m.chunk_ref ? 'B' : 'I');
            else
              _ptr += std::sprintf (_ptr, "%c%s\n", SURFACE_END, m.feature);
          }
        }
        std::memcpy (_ptr, "EOS\n", 4); _ptr += 4;
        char * check = &_res[0];
        while (check < _ptr)
          check += fwrite (check, 1, static_cast <size_t> (_ptr - check), out_fp);
      }
      _ptr = &_res[0];
    }
  };
  template <typename T>
  struct stat_base {
    ny::uint snum;  // # sentence
    ny::uint scorr; // # sentence correctly recognized
    stat_base () : snum (0), scorr (0) {}
    void print () {
      if (! snum) return;
      static_cast <T*> (this)->print_impl ();
      std::fprintf (stderr, "acc. (complete)\t%.4f (%5d/%5d)\n\n",
                    scorr * 1.0 / snum, scorr, snum);
    }
  protected:
    ~stat_base () {}
  };
  struct chunk_stat : public stat_base <chunk_stat> {
    ny::uint pp;  // # chunks correctly recognized
    ny::uint np;  // # chunks incorrectly recognized
    ny::uint pn;  // # chunks missed
    chunk_stat () : pp (0), np (0), pn (0) {}
    void print_impl () {
      std::fprintf (stderr, "J.DepP performance statistics (chunk):\n");
      const double prec (pp * 1.0 / (pp + np)), rec (pp * 1.0 / (pp + pn));
      std::fprintf (stderr, "precision\t%.4f (%5d/%5d)\n", prec, pp, pp + np);
      std::fprintf (stderr, "recall   \t%.4f (%5d/%5d)\n", rec,  pp, pp + pn);
      std::fprintf (stderr, "f1       \t%.4f\n", 2 * prec * rec / (prec + rec));
    }
  };
  struct depnd_stat : public stat_base <depnd_stat> {
    ny::uint bnum;   // # dependencies
    ny::uint bcorr;  // # dependencies correctly recognized
    depnd_stat () : bnum (0), bcorr (0) {}
    void print_impl () {
      std::fprintf (stderr, "J.DepP performance statistics (depnd):\n");
      std::fprintf (stderr, "acc. (partial)\t%.4f (%5d/%5d)\n",
                    bcorr * 1.0 / bnum, bcorr, bnum);
    }
  };
  struct chunk_info {
    ny::uint  id;
    chunk_info (ny::uint id_) : id (id_) {};
    // bool done; // deterministic
    // chunk_info (ny::uint id_, bool done_) : id (id_), done (done_) {};
  };
  class parser : private ny::Uncopyable {
  private:
    option           _opt;
    pecco::option    _pecco_opt;
    pecco::pecco *   _pecco;
    pecco::pecco *   _pecco_chunk;
    pecco::pecco *   _pecco_depnd;
#ifdef USE_OPAL
    opal::option     _opal_opt;
    opal::Model *    _opal;
    opal::mem_pool <opal::ex_t> _ex_pool;
#endif
#ifdef USE_SVM
    TinySVM::Param   _tiny_param;
    TinySVM::Model * _tinysvm;
#endif
#ifdef USE_MAXENT
    maxent_option    _maxent_opt;
    ME_Model *       _libme;
#endif
    // variables to handle events
    Sentence *       _s;
    dict_t *         _dict;
    ny::uint         _fi; // offset of feature index
    ny::fv_t         _fv; // feature vector of the example
    flag_t           _context_feature_bits;
    chunk_stat       _chunk_stat;
    depnd_stat       _depnd_stat;
    FILE *           _writer;
    // timer
#ifdef USE_TIMER
    ny::TimerPool    _timer_pool;
    ny::Timer *      _io_t;
    ny::Timer *      _dict_t;
    ny::Timer *      _preproc_t;
    ny::Timer *      _chunk_t;
    ny::Timer *      _depnd_t;
    ny::Timer *      _classify_t;
#endif
    void _print_ex (const bool flag) const {
      std::fprintf (_writer, "%c1", flag ? '+' : '-');
      for (ny::fv_it it = _fv.begin (); it != _fv.end (); ++it)
        std::fprintf (_writer, " %d:1", *it);
      std::fprintf (_writer, "\n");
    }
    // feature vector generators
    void _add_boolean_feature  (const bool flag)__attribute__((always_inline));
    void _add_boolean_feature  (const bool flag, const bool flag_)__attribute__((always_inline));
    void _add_string_feature   (const ny::uint &id)__attribute__((always_inline));
    void _add_string_feature   (const ny::uint &id, bool flag)__attribute__((always_inline));
    void _add_local_feature    (const ny::uint i, const Bunsetsu &bi, const ny::uint h)__attribute__((always_inline));
    void _add_global_feature   (const ny::uint i, const ny::uint j, const Bunsetsu &bi, const Bunsetsu &bj);
    void _add_cluster_feature  (const Morph &m);
    void _add_cluster_feature  (const Bunsetsu &b, const Morph * ms);
    void _add_morpheme_feature (const ny::uint idx, const int offset);
    void _add_morpheme_feature (const ny::uint idx, const int offset, const bool flag);
    void _add_lexical_feature  (const Bunsetsu &b);
    void _add_lexical_feature  (const Bunsetsu &b, bool flag);
    void _add_context_feature  (const ny::uint from, const ny::uint to, const ny::uint j);
    void _event_gen_from_tuple (const ny::uint i);
    void _event_gen_from_tuple (const ny::uint i, const ny::uint j);
    void _event_gen_from_tuple (const ny::uint i, const ny::uint j, const ny::uint k);
    void _setMorphDic ();
    void _registerMorph (char * cs, const size_t &len, sbag_t &sbag,
                         std::set <ny::uint> &context_feature_ids);
    void _learn ();
#if defined (USE_OPAL) || defined (USE_MAXENT)
    void _processSample (const bool flag);
#endif
    template <const process_t MODE> void _analyze ();
    template <const process_t MODE> void _chunk   ();
    template <const process_t MODE> void _parse   ();
    template <const process_t MODE> void _parseLinear     ();
    template <const process_t MODE> void _parseChunking   ();
    template <const process_t MODE> void _parseBackward   ();
    template <const process_t MODE> void _parseTournament ();
    template <const input_t INPUT>  void _collectStat ();
    void _switch_classifier  (const input_t in);
    void _setup_learner      ();
    void _cleanup_learner    ();
    void _setup_classifier   (const input_t in, int argc, char ** argv);
    void _cleanup_classifier (const input_t in);
#ifdef USE_MAXENT
    void _project (ME_Sample &ms) const;
#endif
  public:
    parser (const option &opt) :
      _opt (opt), _pecco_opt (), _pecco_chunk (), _pecco_depnd (),
#ifdef USE_OPAL
      _opal_opt (), _opal (0), _ex_pool (),
#endif
#ifdef USE_SVM
      _tiny_param (), _tinysvm (0),
#endif
#ifdef USE_MAXENT
      _maxent_opt (), _libme (0),
#endif
      _s (0), _dict (0), _fi (1), _fv (), _context_feature_bits (0), _chunk_stat (), _depnd_stat (), _writer (0)
#ifdef USE_TIMER
      , _timer_pool ("J.DepP profiler"), _io_t (_timer_pool.push ("io")), _dict_t (_timer_pool.push ("dict")), _preproc_t (_timer_pool.push ("preproc", "sent.")), _chunk_t (_timer_pool.push ("chunk", "sent.")), _depnd_t (_timer_pool.push ("depnd", "sent.")), _classify_t (_timer_pool.push ("classify"))
#endif
    {}
    ~parser () { delete _dict; }
    // interface
    void      run         ();
    Sentence* setSentence () { return _s =  new Sentence (); }
  };
}
#endif /* PDEP_H */
