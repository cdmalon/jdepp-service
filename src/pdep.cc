// J.DepP -- Japanese Dependency Parsers
//  $Id: pdep.cc 845 2012-05-17 16:49:17Z ynaga $
// Copyright (c) 2008-2012 Naoki Yoshinaga <ynaga@tkl.iis.u-tokyo.ac.jp>
#include "pdep.h"

namespace pdep {

  static const char * input0[] = { "raw", "chunk", "depnd" };
  
  void parser::_add_context_feature (const ny::uint from, const ny::uint to, const ny::uint j = 0) {
    // assuming _context_feature_bits to be cleared
    for (ny::uint k = from; k < to; ++k) { // merge feature bits;
      if (j)
        while (_s->bun[k].est_did != -1 && _s->bun[k].est_did != static_cast <int> (j))
          k = static_cast <ny::uint> (_s->bun[k].est_did); // only competitors
      if (k >= to) break;
      for (size_t b = 0; b < _context_feature_bits.size (); ++b)
        _context_feature_bits[b] |= _s->bun[k].fbits[b];
    }
    // set context features by twiddling bits
    //   cf. http://www-cs-faculty.stanford.edu/~uno/fasc1a.ps.gz (p. 8)
    for (flag_t::iterator it = _context_feature_bits.begin ();
         it != _context_feature_bits.end (); _fi += FLAG_LEN, ++it)
      while (*it) {
        pecco::byte_4 b (static_cast <float> (*it & - *it)); // pull rightmost 1
        _fv.push_back (_fi + (b.u >> 23) - 127); // pick it via floating number
        *it &= (*it - 1); // unset rightmost 1
      }
  }
  inline void parser::_add_boolean_feature (const bool flag)
  { _fv.push_back (flag ? _fi : _fi + 1); _fi += 2; }
  inline void parser::_add_boolean_feature (const bool flag, const bool flag_)
  { _fv.push_back (flag ? _fi : (flag_ ? _fi + 1 : _fi + 2)); _fi += 3; }
  inline void parser::_add_local_feature (const ny::uint i, const Bunsetsu &bi, const ny::uint h) {
    // mod: bos  / parenthesis / punctuation
    _add_boolean_feature (i == h);
    _add_boolean_feature (bi.bracket > 0, bi.bracket < 0);
    _add_boolean_feature (bi.period);
    _add_boolean_feature (bi.comma);
    // (optional) particles
    _add_context_feature (i, i + 1);
  }
  void parser::_add_global_feature (const ny::uint i, const ny::uint j, const Bunsetsu &bi, const Bunsetsu &bj) {
    // activate global features
    int bracket_in = 0;
    for (ny::uint k = i + 1; k < j; ++k)
      bracket_in += _s->bun[k].bracket;
    // mod-head: distance / parenthesis / punctuation / particle
    _add_boolean_feature (j - i == 1, j - i < 6);          //
    _add_boolean_feature (bracket_in > 0, bracket_in < 0); // , flag
    _add_context_feature (i + 1, j, j);
#ifdef USE_STACKING
    _add_boolean_feature (bi.cand == -1);
    _add_boolean_feature (bi.cand ==  j);
#endif
  }
  inline void parser::_add_string_feature (const ny::uint &id)
  { _fv.push_back (_fi + id); _fi += _dict->num_lexical_features (); }
  inline void parser::_add_string_feature (const ny::uint &id, const bool flag)
  { if (flag) _fv.push_back (_fi + id); _fi += _dict->num_lexical_features (); }
  //
  void parser::_add_cluster_feature (const Morph &m) {
    for (ny::uint i = 0; i < _opt.clen; ++i)
      _add_string_feature (m.cluster (i), (1 << i) & _opt.cbits);
  }
  void parser::_add_cluster_feature (const Bunsetsu &b, const Morph * ms) {
    for (ny::uint i = 0; i < _opt.clen; ++i) {
      if ((1 << i) & _opt.cbits) {
        const ny::fv_it::difference_type pos = _fv.end () - _fv.begin ();
        for (ny::uint j = 0; j < b.mlen; ++i)
          _fv.push_back (_fi + ms[j].cluster (i));
        std::sort (_fv.begin () + pos, _fv.end ());
        _fv.erase (std::unique (_fv.begin () + pos, _fv.end ()), _fv.end ());
        _fi += _dict->num_lexical_features ();
      }
    }
  }
  // set features in a particlular morpheme
  void parser::_add_morpheme_feature (const ny::uint idx, const int offset) {
    const Morph &m = _s->morph[idx + static_cast <ny::uint> (offset)];
    _add_string_feature (m.surf ());
    _add_string_feature (m.pos1 ());
    _add_string_feature (m.pos2 ());
    _add_string_feature (m.infl ());
    if (_opt.clen) _add_cluster_feature (m);
  }
  inline void parser::_add_morpheme_feature (const ny::uint idx, const int offset, const bool flag) {
    if (flag) _add_morpheme_feature (idx, offset);
    else      _fi += (4 + _opt.clen) * _dict->num_lexical_features ();
  }
  void parser::_add_lexical_feature (const Bunsetsu &b) {
    const Morph *  ms  = &_s->morph[b.mpos];
    const ny::fv_it::difference_type pos = _fv.end () - _fv.begin ();
    for (ny::uint i = 0; i < b.mlen; ++i)
      _fv.push_back (_fi + ms[i].surf ());
    std::sort (_fv.begin () + pos, _fv.end ());
    _fv.erase (std::unique (_fv.begin () + pos, _fv.end ()), _fv.end ());
    _fi += _dict->num_lexical_features ();
    if (_opt.clen) _add_cluster_feature (b, ms);
  }
  inline void parser::_add_lexical_feature (const Bunsetsu &b, bool flag) {
    if (flag) _add_lexical_feature (b);
    else      _fi += (1 + _opt.clen) * _dict->num_lexical_features ();
  }
  void parser::_event_gen_from_tuple (const ny::uint i) { // for chunking
    _fi = 1;
    _fv.clear ();
    if (i >= 2) {
      const Morph &m = _s->morph[i - 2];
      _add_string_feature (m.surf ());
      _add_string_feature (m.pos2 ());
      _add_string_feature (m.infl ());
    } else
      _fi += 3 * _dict->num_lexical_features ();
    _add_morpheme_feature (0, static_cast <int> (i) - 1);
    _add_morpheme_feature (0, static_cast <int> (i));
    if (i <= _s->mlen - 2) _add_string_feature (_s->morph[i + 1].surf ());
    else                   _fi += _dict->num_lexical_features ();
    if (i <= _s->mlen - 3) _add_string_feature (_s->morph[i + 2].surf ());
    else                   _fi += _dict->num_lexical_features ();
  }
  void parser::_event_gen_from_tuple (const ny::uint i, const ny::uint j) {
    _fi = 1;
    _fv.clear ();
    const Bunsetsu &bi (_s->bun[i]), &bj (_s->bun[j]);
    if (bi.head >= 0) {
      const Morph &m = _s->morph[bi.mpos + static_cast <ny::uint> (bi.head)];
      _add_string_feature (m.surf ());
      _add_string_feature (m.pos2 ());
    } else
      _fi += 2 * _dict->num_lexical_features ();
    _add_morpheme_feature (bi.mpos, bi.tail,  bi.tail >= 0);
    {
      const Morph &m = _s->morph[bj.mpos];
      _add_string_feature (m.surf ());
      _add_string_feature (m.pos2 ());
    }
    _add_morpheme_feature (bj.mpos, bj.head, bj.head >= 0);
    _add_morpheme_feature (bj.mpos, bj.tail, bj.tail >= 0); // * maybe useless
    if (j < _s->blen - 1) {
      const Bunsetsu &bk = _s->bun[j + 1];
      const Morph * m = &_s->morph[bk.mpos];
      _add_string_feature (m[0].surf ());
      if (bk.head >= 0) _add_string_feature (m[bk.head].surf ());
      _add_string_feature (m[bk.mlen - 1].surf ());
    } else
      _fi += 3 * _dict->num_lexical_features ();
    if (i > 0) {
      const Bunsetsu &bk = _s->bun[i - 1];
      const Morph * m = &_s->morph[bk.mpos];
      _add_string_feature (m[bk.mlen - 1].surf ());
    } else
      _fi += 1 * _dict->num_lexical_features ();
    _add_local_feature (i, bi, 0);
    _add_local_feature (j, bj, _s->blen - 1);
    _add_global_feature (i, j, bi, bj);
#ifdef USE_STACKING
#endif
    // add new features (>= _fi)  here
  }
  void parser::_event_gen_from_tuple (const ny::uint i, const ny::uint j, const ny::uint k) {
    _fi = 1;
    _fv.clear ();
    const Bunsetsu &bi (_s->bun[i]), &bj (_s->bun[j]), &bk (_s->bun[k]);
    if (bi.head >= 0) {
      const Morph &m = _s->morph[bi.mpos + static_cast <ny::uint> (bi.head)];
      _add_string_feature (m.surf ());
      _add_string_feature (m.pos2 ());
    } else
      _fi += 2 * _dict->num_lexical_features ();
    _add_morpheme_feature (bi.mpos, bi.tail,  bi.tail >= 0);
    {
      const Morph &m = _s->morph[bj.mpos];
      _add_string_feature (m.surf ());
      _add_string_feature (m.pos2 ());
    }
    _add_morpheme_feature (bj.mpos, bj.head, bj.head >= 0);
    _add_morpheme_feature (bj.mpos, bj.tail, bj.tail >= 0); // * maybe useless
    {
      const Morph &m = _s->morph[bk.mpos];
      _add_string_feature (m.surf ());
      _add_string_feature (m.pos2 ());
    }
    _add_morpheme_feature (bk.mpos, bk.head, bk.head >= 0);
    _add_morpheme_feature (bk.mpos, bk.tail, bk.tail >= 0); // * maybe useless
    if (j < _s->blen - 1) {
      const Bunsetsu &bl = _s->bun[j + 1];
      const Morph * m = &_s->morph[bl.mpos];
      _add_string_feature (m[0].surf ());
      if (bl.head >= 0) _add_string_feature (m[bl.head].surf ());
      _add_string_feature (m[bl.mlen - 1].surf ());
    } else
      _fi += 3 * _dict->num_lexical_features ();
    if (k < _s->blen - 1) {
      const Bunsetsu &bl = _s->bun[k + 1];
      const Morph * m = &_s->morph[bl.mpos];
      _add_string_feature (m[0].surf ());
      if (bl.head >= 0) _add_string_feature (m[bl.head].surf ());
      _add_string_feature (m[bl.mlen - 1].surf ());
    } else
      _fi += 3 * _dict->num_lexical_features ();
    if (i > 0) {
      const Bunsetsu &bl = _s->bun[i - 1];
      const Morph *m = &_s->morph[bl.mpos];
      _add_string_feature (m[bl.mlen - 1].surf ());
    } else
      _fi += 1 * _dict->num_lexical_features ();
    //
    _add_local_feature (i, bi, 0);
    _add_local_feature (j, bj, _s->blen - 1);
    _add_local_feature (k, bk, _s->blen - 1);
    _add_global_feature (i, j, bi, bj);
    _add_global_feature (i, k, bi, bk);
#ifdef USE_STACKING
    _add_lexical_feature (_s->bun[_s->bun[i].cand], j != _s->bun[i].cand);
#endif
    // add new features (>= _fi)  here
  }

#if defined (USE_OPAL) || defined (USE_MAXENT)
  inline void parser::_processSample (const bool flag) {
    if (_opt.learner == OPAL) {
#ifdef USE_OPAL
      opal::ex_t x;
      _opal->set_ex (x, flag ? +1 : -1, _fv, true, _opal_opt.kernel == opal::POLY);
      _ex_pool.put (x);
#endif
    } else if (_opt.learner == MAXENT) {
#ifdef USE_MAXENT
      static ME_Sample ms;
      ms.label = flag ? "+1" : "-1";
      ms.features.clear ();
      _project (ms);
      _libme->add_training_sample (ms);
#endif
    }
  }
#endif
  template <const process_t MODE>
  void parser::_parseLinear () {
    static std::stack <ny::uint> stack;
    const ny::uint len = _s->blen;
    for (ny::uint j = 1; j < len; ++j) {
      stack.push (j - 1);
      while (! stack.empty ()) {
        const ny::uint i = stack.top ();
        Bunsetsu &b = _s->bun[i];
        b.depnd_prob = 0.0;
        if (j != len - 1) {
          _event_gen_from_tuple (i, j);
          bool flag = (j == static_cast <ny::uint> (b.did));
          // output example for training / fstrie construction
          if (MODE != PARSE) _print_ex (flag);
          if (MODE != LEARN)
            flag = (b.depnd_prob = _pecco->getProbability (_fv)) > 0.5;
#if defined (USE_OPAL) || defined (USE_MAXENT)
          else
            _processSample (flag); // learn
#endif
          if (! flag) break;
        }
        b.est_did = static_cast <int> (j);
        stack.pop ();
      }
    }
  }
  template <const process_t MODE>
  void parser::_parseChunking () {
    static std::list <chunk_info> cinfo;
    const ny::uint len = _s->blen;
    for (ny::uint i = 0; i < len - 1; ++i)
      cinfo.push_back (chunk_info (i));
    // cinfo.push_back (chunk_info (i, false)); // deterministic
    while (! cinfo.empty ()) {
      typename std::list <chunk_info>::reverse_iterator it = cinfo.rbegin ();
      bool next = true; // D
      ny::uint i = 0; // bug fix
      ny::uint j = it->id;
      _s->bun[j].depnd_prob = 0.0;
      _s->bun[j].est_did = static_cast <int> (len - 1);
      for (++it; it != cinfo.rend (); ++it) {
        i = it->id;
        // if (! it->done) { // deterministic (assuming no dynamic features)
        _event_gen_from_tuple (i, j);
        _s->bun[i].depnd_prob = 0.0;
        bool flag = (j == static_cast <ny::uint> (_s->bun[i].did));
        if (MODE != PARSE) _print_ex (flag);
        if (MODE != LEARN)
          flag = (_s->bun[i].depnd_prob = _pecco->getProbability (_fv)) > 0.5;
#if defined (USE_OPAL) || defined (USE_MAXENT)
        else               _processSample (flag);
#endif
        if (flag) _s->bun[i].est_did = static_cast <int> (j);
        // it->done = true; // deterministic
        // }
        if (_s->bun[i].est_did == -1 && next) { // O & D
          it = std::list <chunk_info>::reverse_iterator (cinfo.erase (it.base ()));
          // it->done = false; // deterministic
        }
        next = _s->bun[i].est_did != -1;
        j = i;
      }
      if (_s->bun[i].est_did != -1) cinfo.erase (it.base ());
    }
  }
  template <const process_t MODE>
  void parser::_parseBackward () {
    const ny::uint len = _s->blen;
    if (len < 2) return;
    for (ny::uint i = len - 2;; --i) {
      _s->bun[i].depnd_prob = 0.0; // non-deterministic
      for (ny::uint j = i + 1;; j = static_cast <ny::uint> (_s->bun[j].est_did)) { // non-deterministic
        // for (ny::uint j = i + 1; j != len - 1; j = static_cast <ny::uint> (_s->bun[j].est_did)) { // deterministic
      
        double prob = 0;
        _event_gen_from_tuple (i, j);
        bool flag = (j == static_cast <ny::uint> (_s->bun[i].did));
        // output example for training / fstrie construction
        if (MODE != PARSE) _print_ex (flag);
        if (MODE == LEARN) { // learn
#if defined (USE_OPAL) || defined (USE_MAXENT)
          _processSample (flag);
#endif
          if (flag) // non-deterministic
            { _s->bun[i].est_did = static_cast <int> (j); prob = 1.0; }
          // if (flag) // deterministic
          //   { _s->bun[i].est_did = static_cast <int> (j); prob = 1.0; break; }
        } else {
          prob = _pecco->getProbability (_fv);
        }
        if (prob > _s->bun[i].depnd_prob) {// non-deterministic
          _s->bun[i].est_did = static_cast <int> (j);
          _s->bun[i].depnd_prob = prob;
        }
        // if (prob >= 0.5) { // deterministic
        //   _s->bun[i].est_did = static_cast <int> (j);
        //   _s->bun[i].depnd_prob = prob;
        //   break;
        // } 
        if (_s->bun[j].est_did == -1) break; // non-deterministic
      }
      // if (_s->bun[i].est_did == -1) // deterministic
      //   _s->bun[i].est_did = static_cast <int> (len) - 1; 
      if (i == 0) break;
    }
  }
  template <const process_t MODE>
  void parser::_parseTournament () {
    const ny::uint len = _s->blen;
    if (len < 2) return;
    if (MODE == LEARN) { // I don't merge the two outermost loops
      for (ny::uint i = 0; i < len - 2; ++i) {
        const ny::uint h = static_cast <ny::uint> (_s->bun[i].did); // head
        for (ny::uint j = i + 1; j <= len - 1; ++j) {
          bool flag = true;
          if      (j < h) { flag = true;  _event_gen_from_tuple (i, j, h); }
          else if (j > h) { flag = false; _event_gen_from_tuple (i, h, j); }
          else continue;
#if defined (USE_OPAL) || defined (USE_MAXENT)
          _processSample (flag);
#endif
          _print_ex (flag);
        }
      }
    } else {
      for (ny::uint i = len - 2;; --i) {
        ny::uint did = i + 1; // head
        ny::uint j   = did;
        while (_s->bun[j].est_did != -1) { // head of head
          j = static_cast <ny::uint> (_s->bun[j].est_did);
          _event_gen_from_tuple (i, did, j);
          if (MODE == CACHE)
            _print_ex (did < static_cast <ny::uint> (_s->bun[i].did));
          bool flag = (_s->bun[i].depnd_prob = _pecco->getProbability (_fv)) > 0.5;
          if (flag) did = j; // RIGHT
        }
        _s->bun[i].est_did = static_cast <int> (did);
        if (i == 0) break;
      }
    }
  }
  template <const process_t MODE>
  void parser::_parse () {
    TIMER (if (MODE != LEARN) _depnd_t->startTimer ());
    if (MODE != LEARN) _switch_classifier (DEPND);
    switch (_opt.parser) {
      case LINEAR:     _parseLinear     <MODE> (); break;
      case CHUNKING:   _parseChunking   <MODE> (); break;
      case TOURNAMENT: _parseTournament <MODE> (); break;
      case BACKWARD:   _parseBackward   <MODE> (); break;
    }
    TIMER (if (MODE != LEARN) _depnd_t->stopTimer ());
    if (MODE == PARSE && _opt.input != RAW) _collectStat <DEPND> ();
  }
  template <const process_t MODE>
  inline void parser::_chunk () {
    TIMER (if (MODE != LEARN) _chunk_t->startTimer ());
    if (MODE != LEARN) _switch_classifier (CHUNK);
    _s->addBunsetsu (0);
    _s->morph[0].chunk = true;
    const ny::uint mlen = _s->mlen;
    for (ny::uint i = 1; i < mlen; ++i) {
      Morph &m = _s->morph[i];
      m.chunk = m.chunk_ref;
      _event_gen_from_tuple (i);
      // output example for training / fstrie construction
      if (MODE != PARSE) _print_ex (m.chunk_ref);
      if (MODE != LEARN)
        m.chunk = (m.chunk_prob  = _pecco->getProbability (_fv)) > 0.5;
#if defined (USE_OPAL) || defined (USE_MAXENT)
      else
        _processSample (m.chunk_ref);
#endif
      if (m.chunk) _s->addBunsetsu (i);
    }
    TIMER (if (MODE != LEARN) _chunk_t->stopTimer ());
    if (MODE == PARSE && _opt.input != RAW) _collectStat <CHUNK> ();
  }
  template <>
  void parser::_collectStat <CHUNK> () {
    ++_chunk_stat.snum;
    bool flag (true), prev (true);
    for (ny::uint i = 1; i < _s->mlen; ++i) {
      const Morph &m = _s->morph[i];
      if (m.chunk & m.chunk_ref) { // 'B' & 'B'
        if (prev) ++_chunk_stat.pp; else ++_chunk_stat.np, ++_chunk_stat.pn;
        prev = true;
      } else if (m.chunk | m.chunk_ref) { // 'B' & 'I' or 'I' & 'B'
        if (m.chunk) ++_chunk_stat.np; else ++_chunk_stat.pn;
        flag = false;
        prev = false;
      }
    }
    if (prev) ++_chunk_stat.pp; else  ++_chunk_stat.np, ++_chunk_stat.pn;
    if (flag) ++_chunk_stat.scorr;
  }
  template <> void parser::_collectStat <DEPND> () {
    if (_s->blen >= 1) { // set >= 2 to ignore sentene w/ a single bunsetsu
      ++_depnd_stat.snum;
      _depnd_stat.bnum += _s->blen - 1;
      ny::uint bcorr = 0;
      const ny::uint len = _s->blen;
      for (ny::uint i = 0; i < len - 1; ++i)
        if (_s->bun[i].est_did == _s->bun[i].did) ++bcorr;
      // std::printf ("%d\n", b.est_did == b.did ? 1 : 0); // for McNemar
      _depnd_stat.bcorr += bcorr;
      if (bcorr == len - 1) ++_depnd_stat.scorr;
      // std::printf ("%d\n", _flag & 0x2 ? 1 : 0); // for McNemar
    }
  }
  template <const process_t MODE>
  void parser::_analyze () {
    setSentence ();
    size_t n = 0;
    if (MODE == PARSE)
      std::fprintf (stderr, "(input: STDIN [-I %d])\n", _opt.input);
#ifdef USE_AS_STANDALONE
    if (_opt.input == RAW) {
      std::string mecab_opt ("$0");
      if (_opt.mecab_dic) mecab_opt += std::string (" -d ") + _opt.mecab_dic;
      MeCab::Tagger * tagger = MeCab::createTagger (mecab_opt.c_str ());
      if (! tagger)
        ny::print_err (HERE "fail to invoke MeCab; you may want to set [-d].\n");
      char header[1024];
      char * line = 0;
      size_t read = 0;
      char buf[IOBUF_SIZE]; std::setvbuf (stdin, &buf[0], _IOFBF, IOBUF_SIZE);
      while (1) {
        TIMER (_io_t->startTimer ());
        const bool input_ok = ny::getLine (stdin, line, read);
        TIMER (_io_t->stopTimer ());
        if (! input_ok) break;
        TIMER (_preproc_t->startTimer ());
        size_t header_len = std::sprintf (header, "# S-ID: %ld; J.DepP\n", ++n);
        _s->setHeader (header, header_len);
        for (MeCab::Node * m = tagger->parseToNode (line, read - 1)->next;
             m->stat != MECAB_EOS_NODE; m = m->next) // set_feature
          _s->addMorph (m->length, m->surface, m->feature, _dict);
        TIMER (_preproc_t->stopTimer ());
        _chunk <PARSE> ();
        _s->setup (_dict); // bug?
        _parse <PARSE> ();
        TIMER (_io_t->startTimer ());
        _s->print (_opt.input);
        TIMER (_io_t->stopTimer  ());
        _s->clear ();
      }
      std::fclose (stdin);
      delete tagger;
    } else {
#endif
      static const char * mode[] = { "learn", "parse", "both", "cache" };
      int fd = 0; // stdin
      if (MODE != PARSE) {
        if ((fd = open (_opt.train, O_RDONLY)) == -1)
          ny::print_err (HERE "no such file: %s\n", _opt.train);
        std::string model (_opt.model_dir + "/" + input0[_opt.input]);
        if (_opt.input == DEPND) {
          char sigparse[16]; std::sprintf (sigparse, ".p%d", _opt.parser);
          model += sigparse;
        }
        if (MODE == LEARN) {
          const std::string train (model + ".train");
          _writer = std::fopen (train.c_str (), "wb");
        } else {
          std::string event (model + ".event");
#ifdef USE_MODEL_SUFFIX
          if (_opt.learner == SVM)
            event += std::string (".s") + _pecco_opt.sigma;
#endif
          _writer = std::fopen (event.c_str (), "wb");
        }
      }
      const bool output = _opt.input == RAW ||
                          (MODE == PARSE && _opt.verbose == -1);
      char buf[IOBUF_SIZE];
      char * q (&buf[0]), * q_end (&buf[0] + IOBUF_SIZE);
      ssize_t avail = 0;
      std::vector <char *> pos;
      while (1) {
        TIMER (if (MODE != LEARN) _io_t->startTimer ());
        const bool input_ok
          = (avail = read (fd, q, static_cast <size_t> (q_end - q))) > 0;
        TIMER (if (MODE != LEARN) _io_t->stopTimer ());
        if (! input_ok) break;
        q_end = q + avail;
        q = &buf[0];
        while (1) {
          TIMER (if (MODE != LEARN) _io_t->startTimer ());
          // read input
          pos.clear ();
          for (char * r = q; r + 3 < q_end; ++r) {
            pos.push_back (r);
            if (*r == 'E') { q = r + 4; break; }  // find EOS\n
            while (r != q_end && *r != '\n') ++r; // next line
          }
          if (pos.empty () || *pos.back () != 'E') { // premature input
            std::memmove (&buf[0], q, static_cast <size_t> (q_end - q));
            q     = &buf[0] + (q_end - q);
            q_end = &buf[0] + IOBUF_SIZE;
            break;
          }
          TIMER (if (MODE != LEARN) _io_t->stopTimer ());
          // process
          TIMER (if (MODE != LEARN) _preproc_t->startTimer ());
          _s->clear ();
          ++n;
          if (_opt.input == RAW) {
            char header[1024];
            size_t header_len
              = static_cast <size_t> (std::sprintf (header, "# S-ID: %ld; J.DepP\n", n));
            _s->setHeader (header, header_len);
          }
          bool flag = false; // reference chunk annotation
          for (size_t i = 0; i < pos.size () - 1; ++i) {
            char * line = pos[i];
            const size_t read = static_cast <size_t> (pos[i + 1] - pos[i]);
            if (_opt.input == RAW) _s->addMorph (line, read - 1, _dict);
            else
              switch (*line) {
                case '#': if (output) _s->setHeader (line, read); break;
                case '*':
                  if (_opt.input == DEPND)
                    _s->addBunsetsu (line, read - 1, _s->mlen);
                  else
                    flag = true;
                  break;
                default:
                  _s->addMorph (line, read - 1, _dict, flag);
                  flag = false;
              }
          }
          TIMER (if (MODE != LEARN) _preproc_t->stopTimer ());
          if (_opt.input != DEPND) _chunk <MODE> ();
          _s->setup (_dict); // bug?
          if (_opt.input != CHUNK) _parse <MODE> ();
          TIMER (if (MODE != LEARN) _io_t->startTimer ());
          if (output) _s->print (_opt.input);
          TIMER (if (MODE != LEARN) _io_t->stopTimer ());
          if (MODE != PARSE && n % 1000 == 0)
            std::fprintf (stderr, "\r%s: %ld sent. processed", mode[MODE], n);
          if (MODE != PARSE && _opt.max_sent && n >= _opt.max_sent)
            { lseek (fd, 0, SEEK_END); break; } // go to the end
        }
      }
      close (fd);
      if (avail == 0 && q == q_end)
        ny::print_err ("set a larger value to IOBUF_SIZE.\n");
      if (MODE != PARSE) {
        std::fprintf (stderr, "\r%s: %ld sent. processed.", mode[MODE], n);
        std::fclose (_writer);
      }
      std::fprintf (stderr, "\n");
#ifdef USE_AS_STANDALONE
    }
#endif
    delete _s;
  }
  void parser::_learn () {
    std::string model (_opt.model_dir + "/" + input0[_opt.input]);
    if (_opt.input == DEPND) {
      char sigparse[16]; std::sprintf (sigparse, ".p%d", _opt.parser);
      model += sigparse;
    }
    switch (_opt.learner) { // they should be newed earlier
#ifdef USE_OPAL
      case OPAL:
        _opal->train (_ex_pool, _opal_opt.iter);
        _opal->save (model.c_str ());
        break;
#endif
#ifdef USE_SVM
      case SVM:
        {
          const std::string train (model + ".train");
          TinySVM::Example * ex = new TinySVM::Example;
          ex->read (train.c_str ());
          _tinysvm = ex->learn (_tiny_param);
          _tinysvm->write (model.c_str ());
          delete ex;
          break;
        }
#endif
#ifdef USE_MAXENT
      case MAXENT:
        switch (_maxent_opt.algo) {
          case SGD:   _libme->use_SGD ();
          case OWLQN: _libme->use_l1_regularizer (_maxent_opt.reg_cost); break;
          case LBFGS: _libme->use_l2_regularizer (_maxent_opt.reg_cost); break;
          default: ny::print_err ("MAXENT optimizer disabled.\n");
        }
        _libme->train ();
        _libme->save_to_file (model);
        break;
#endif
      default: ny::print_err ("learner disabled.\n");
    }
  }
  void parser::_registerMorph (char * cs, const size_t &len, sbag_t &sbag,
                               std::set <ny::uint> &context_feature_ids) {
    ny::uint surf (0), i (0);
    const char * post_particle = _opt.utf8 ? UTF8_POST_PARTICLE : EUC_POST_PARTICLE;
    for (char * p (cs), * const p_end (p + len);
         cs < p_end && i < NUM_FIELD; cs = ++p, ++i) {
      if (i == SURF) while (p != p_end && *p != SURFACE_END) ++p;
      else           while (p != p_end && *p != FEATURE_SEP) ++p;
      *p = '\0';
      // if (i <= POS2 || i == INFL) {
      if (i == SURF || i == POS1 || i == POS2 || i == INFL) {
        sbag_t::iterator it = sbag.find (cs);
        if (it == sbag.end ()) { // new entry
          char * copy = new char[p - cs + 1];
          std::strcpy (copy, cs);
          it = sbag.insert (sbag_t::value_type (copy, static_cast <ny::uint> (sbag.size ()))).first;
        }
        const ny::uint id = it->second;
        switch (i) {
          case SURF: surf = id; break;
          case POS1:
            if (std::strcmp (cs, post_particle) == 0)
              context_feature_ids.insert (surf);
            break;
          default: break;
        }
      }
    }
#ifdef USE_JUMAN_POS
    if (i < NUM_FIELD)
      ny::print_err (HERE "# fields is less than %d.\n", NUM_FIELD);
#endif
  }
  // morphological dictionary are extracted from the training data
  void parser::_setMorphDic () {
    if (_opt.verbose > 0) std::fprintf (stderr, "Loading dict..");
    const std::string dict (_opt.model_dir + "/dic" + (_opt.utf8 ? ".utf8" : ".euc"));
    struct stat st;
    if (stat (dict.c_str (), &st) != 0) {
      if (_opt.verbose > 0)
        std::fprintf (stderr, "not found; reading %s..", _opt.train);
      sbag_t              sbag;
      std::set <ny::uint> context_feature_ids;
      FILE * reader = std::fopen (_opt.train, "r");
      if (! reader) ny::print_err (HERE "no such file: %s\n", _opt.train);
      char * line = 0;
      size_t read = 0;
      while (ny::getLine (reader, line, read))
        if (*line != '*' && *line != '#' && *line != 'E')
          _registerMorph (line, read - 1, sbag, context_feature_ids);
      std::fclose (reader);
      // add some mandatory features
      sbag.insert (sbag_t::value_type (_opt.utf8 ? UTF8_COMMA : EUC_COMMA,
                                       static_cast <ny::uint> (sbag.size ())));
      sbag.insert (sbag_t::value_type (_opt.utf8 ? UTF8_PERIOD : EUC_PERIOD,
                                       static_cast <ny::uint> (sbag.size ())));
      sbag.insert (sbag_t::value_type (_opt.utf8 ? UTF8_POST_PARTICLE : EUC_POST_PARTICLE,
                                       static_cast <ny::uint> (sbag.size ())));
      sbag.insert (sbag_t::value_type (_opt.utf8 ? UTF8_BRACKET_START : EUC_BRACKET_START,
                                       static_cast <ny::uint> (sbag.size ())));
      sbag.insert (sbag_t::value_type (_opt.utf8 ? UTF8_BRACKET_END : EUC_BRACKET_END,
                                       static_cast <ny::uint> (sbag.size ())));
      sbag.insert (sbag_t::value_type (_opt.utf8 ? UTF8_SPECIAL : EUC_SPECIAL,
                                       static_cast <ny::uint> (sbag.size ())));
#ifdef USE_JUMAN_POS
      sbag.insert (sbag_t::value_type (_opt.utf8 ? UTF8_SUFFIX : EUC_SUFFIX,
                                       static_cast <ny::uint> (sbag.size ())));
#endif
      if (_opt.verbose > 0)
        std::fprintf (stderr, "done.\n");
      const ny::uint num_lexical_features = static_cast <ny::uint> (sbag.size ());
      const ny::uint num_context_features = static_cast <ny::uint> (context_feature_ids.size ());
      if (num_context_features == 0)
        ny::print_err (HERE "no particles found in %s\n"
                       "\tthe charset / posset may mismatch with jdepp_conf.h\n",
                       _opt.train);
      FILE * writer = std::fopen (dict.c_str (), "wb");
      std::fwrite (&num_lexical_features, sizeof (ny::uint), 1, writer);
      std::fwrite (&num_context_features, sizeof (ny::uint), 1, writer);
      std::fclose (writer);
      //
      std::vector <const char *> str;
      std::vector <size_t>       len;
      std::vector <int>          val;
      for (sbag_t::const_iterator it = sbag.begin (); it != sbag.end (); ++it) {
        str.push_back (it->first);
        len.push_back (std::strlen (it->first));
        // move (gap) lexical featuress forward and assign smaller IDs
        const ny::uint id = it->second;
        std::set <ny::uint>::iterator jt = context_feature_ids.lower_bound (id);
        const int index
          = static_cast <int> (std::distance (context_feature_ids.begin (), jt));
        val.push_back (*jt == id ? index: static_cast <int> (id + num_context_features) - index);
      }
      ny::trie t;
      pecco::build_trie (&t, "dict trie", dict, str, len, val, _opt.verbose > 0, "ab");
      for (sbag_t::iterator it = sbag.begin (); it != sbag.end (); ++it)
        delete [] it->first;
    }
    _dict = new dict_t (dict.c_str (), _opt.utf8);
    if (_opt.verbose > 0)
      std::fprintf (stderr, "done. (# strings + 1 (unseen) = %d).\n",
                    _dict->num_lexical_features ());
    _context_feature_bits.resize (_dict->context_feature_bit_len (), 0);
  }
  void parser::_setup_learner () {
    switch (_opt.learner) { // they should be newed earlier
#ifdef USE_OPAL
      case OPAL:
        _opal_opt.set (_opt.learner_argc, _opt.learner_argv);
        _opal = new opal::Model (_opal_opt);
        break;
#endif
#ifdef USE_SVM
      case SVM:
        _tiny_param.set (_opt.learner_argc, _opt.learner_argv);
        if (_tiny_param.kernel_type != TinySVM::POLY || _tiny_param.degree > 4)
          ny::print_err ("only polynomial kernel [-t 1] of [-d] <= 4 is supported in SVM.\n");
        break;
#endif
#ifdef USE_MAXENT
      case MAXENT:
        _maxent_opt.set (_opt.learner_argc, _opt.learner_argv);
        _libme = new ME_Model;
        break;
#endif
      default: ny::print_err ("learner disabled.\n");
        break;
    }
  }
  void parser::_cleanup_learner () {
    switch (_opt.learner) {
#ifdef USE_OPAL
      case OPAL:
        // _opal->printStat ();
        delete _opal;    break;
#endif
#ifdef USE_SVM
      case SVM:    delete _tinysvm; break;
#endif
#ifdef USE_MAXENT
      case MAXENT: delete _libme;   break;
#endif
      default: ny::print_err ("learner disabled.\n");
    }
  }
  void parser::_switch_classifier (const input_t in) // exchange classifiers
  { _pecco = in == CHUNK ? _pecco_chunk : _pecco_depnd; }
  void parser::_setup_classifier (const input_t in, int argc, char ** argv) {
    std::string model (_opt.model_dir + "/" + input0[in]);
    if (in == DEPND) {
      char sigparse[16]; std::sprintf (sigparse, ".p%d", _opt.parser);
      model += sigparse;
    }
    if (_opt.mode == BOTH) { // induce an appropriate classifier from input model
#ifdef USE_OPAL
      if (_opt.learner == OPAL && _opal_opt.kernel == opal::POLY)
        _opt.learner = SVM; // fake
#endif
    } else {
      FILE * fp = std::fopen (model.c_str (), "r");
      if (! fp || std::feof (fp))
        ny::print_err (HERE "no model found: %s; train a model first [-t 0]\n",
                       model.c_str ());
      switch (std::fgetc (fp)) {
        case  0 :
        case '#': _opt.learner = OPAL;   break;
        case 'o': // delegate
        case 'T': _opt.learner = SVM;    break;
        case '-':
        case '+': _opt.learner = MAXENT; break;
        default:  ny::print_err (HERE "unknown model type found: %s;");
      }
#ifndef USE_OPAL
      if (_opt.learner == OPAL)
        ny::print_err (HERE "unsupported model found, try jdepp_pa\n");
#endif
      std::fclose (fp);
    }
    if  (_opt.learner == OPAL) { // opal as a classifier
#ifdef USE_OPAL
      opal::option opal_opt (argc, argv);
      opal_opt.model = model.c_str ();
      _pecco = new pecco::pecco (opal_opt, reinterpret_cast <opal::Model *> (0));
      _pecco->load (opal_opt.model);
#endif
    } else {
      const std::string train (model + ".train");
      const std::string event (model + ".event");
      _pecco_opt.set (argc, argv);
      _pecco_opt.model = model.c_str ();
      _pecco_opt.train = train.c_str ();
      _pecco_opt.event = event;
      if (_opt.learner == SVM) {
#if defined (USE_SVM) || defined (USE_OPAL)
#ifdef USE_MODEL_SUFFIX
        _pecco_opt.event += std::string (".s") + _pecco_opt.sigma; // approx
#endif
        _pecco = new pecco::pecco (_pecco_opt,
                                   static_cast <pecco::kernel_model *> (0));
#endif
      } else {
#ifdef USE_MAXENT
        _pecco = new pecco::pecco (_pecco_opt,
                                   static_cast <pecco::linear_model *> (0));
#endif
      }
      _pecco->load (_pecco_opt.model);
    }
    if (in == CHUNK)  _pecco_chunk = _pecco; else _pecco_depnd = _pecco;
  }
  void parser::_cleanup_classifier (const input_t in) {
    _switch_classifier (in);
    // _pecco->printStat ();
    delete _pecco;
  }
#ifdef USE_MAXENT
  void parser::_project (ME_Sample &ms) const {
    char key[64];
    for (ny::fv_it beg (_fv.begin ()), end (_fv.end ()), it (beg);
         it != end; ++it) {
      const int pos_i = std::sprintf (key, "%d", *it);
      ms.add_feature (key);
      if (_maxent_opt.degree != 1)
        for (ny::fv_it jt = beg; jt != it; ++jt) {
          const int pos_j = pos_i + std::sprintf (key + pos_i, ":%d", *jt);
          ms.add_feature (key);
          if (_maxent_opt.degree != 2)
            for (ny::fv_it kt = beg; kt != jt; ++kt) {
              std::sprintf (key + pos_j, ":%d", *kt);
              ms.add_feature (key);
            }
        }
    }
  }
#endif
  void parser::run () {
#ifdef USE_TIMER
    if (_opt.verbose > 0)
      std::fprintf (stderr, "Processor Speed: %.3LfGHz.\n",
                    ny::Timer::clock / 1000);
#endif
    if (_opt.input == RAW && _opt.mode != PARSE)
#ifdef USE_AS_STANDALONE
      ny::print_err ("You can input RAW sentences [-I 0] only for parsing [-t 1].\n");
#else
    ny::print_err ("You can input POS-tagged sentences [-I 0] only for parsing [-t 1].\n");
#endif

    TIMER (_dict_t->startTimer ());
    _setMorphDic ();
    TIMER (_dict_t->stopTimer ());
    // learn
    if  (_opt.mode == LEARN || _opt.mode == BOTH) { // learn or both
      _setup_learner ();
      _analyze <LEARN> ();
      _learn ();
      _cleanup_learner ();
    }
    // parse, cache
    if (_opt.mode != LEARN) {
      if (_opt.input != DEPND)
        _setup_classifier (CHUNK, _opt.chunk_argc, _opt.chunk_argv);
      if (_opt.input != CHUNK)
        _setup_classifier (DEPND, _opt.depnd_argc, _opt.depnd_argv);
      if (_opt.mode == CACHE) {
        if (_opt.learner == OPAL)
          ny::print_err ("needless to cache in opal classifier [-t 0].\n");
        _analyze <CACHE> ();
      } else {
        _analyze <PARSE> ();
        if (_opt.input == CHUNK) _chunk_stat.print ();
        if (_opt.input == DEPND) _depnd_stat.print ();
      }
      // cleanup
      TIMER (_timer_pool.print ());
      if (_opt.input != DEPND) _cleanup_classifier (CHUNK);
      if (_opt.input != CHUNK) _cleanup_classifier (DEPND);
    }
  }
}
// ToDo
// - output model info, or at least add signature
// - implement POS tagger / word segmenter faster than MeCab
