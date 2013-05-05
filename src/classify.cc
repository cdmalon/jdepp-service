// pecco -- please enjoy classification with conjunctive features
//  $Id: classify.cc 827 2012-05-10 14:12:37Z ynaga $
// Copyright (c) 2008-2012 Naoki Yoshinaga <ynaga@tkl.iis.u-tokyo.ac.jp>
#include "classify.h"
#ifdef USE_KERNEL
#include "kernel.h"
#endif
#ifdef USE_LINEAR
#include "linear.h"
#endif

#define SHIFT_LEFT_UNSIGNED_USES_SHL (((unsigned int)0xffffffff >> 1) == 0x7fffffff)

namespace pecco {
  
  // radix sort for very long keys
  template <typename T>
  template <typename value_type, typename SrcIterator, typename OutIterator>
  inline void ClassifierBase <T>::_radix (SrcIterator first, SrcIterator last,
                                          OutIterator dest, size_t shift) {
    value_type count[BIN_ELM];
    std::memset (count, 0, sizeof (count));
    for (SrcIterator it = first; it != last; ++it)
      ++count[((*it) >> shift) & BIN_MAX];
    value_type offset[BIN_ELM];
    offset[0] = 0;
    std::partial_sum (&count[0], &count[(BIN_ELM)-1], &offset[1]);
  
    for (SrcIterator it = first; it != last; ++it)
      *(dest + offset[((*it)>>shift) & BIN_MAX]++) = *it;
  }

  template <typename T>
  template <size_t SIZE, typename Iterator>
  inline void ClassifierBase <T>::_radix_sort (const Iterator &first, const Iterator &last)
  {
    typedef typename std::iterator_traits<Iterator>::value_type value_type;
    if (first == last) return;
    const size_t n = static_cast <size_t> (std::distance (first, last));
    static std::vector <value_type> temp;
    if (temp.size () < n) temp.resize (n, value_type (0));
    for (size_t i = 0; i < SIZE / BIN_BITS; ++i) {
      _radix <value_type> (first, last, &temp[0], i * BIN_BITS);
      _radix <value_type> (&temp[0], &temp[n], first, ++i * BIN_BITS);
    }
  }
  // for short feature vector < 50
  template <typename T>
  template <typename Iterator>
  inline void ClassifierBase <T>::_insertion_sort (const Iterator &first,
                                                   const Iterator &last) {
    typedef typename std::iterator_traits <Iterator>::value_type value_type;
    std::less <value_type> cmp;
    if (first == last) return;
    Iterator sorted = first;
    for (++sorted; sorted != last; ++sorted) {
      value_type temp = *sorted;
      Iterator prev, curr;
      prev = curr = sorted;
      for (--prev; curr != first && cmp (temp, *prev); --prev, --curr)
        *curr = *prev;
      *curr = temp;
    }
  }
  // Zipf-aware sorting funciton
  template <typename T>
  template <typename Iterator>
  inline void ClassifierBase <T>::_bucket_sort (const Iterator &first, const Iterator &last)
  {
    typedef typename std::iterator_traits <Iterator>::value_type value_type;
    if (first == last) return;
    ny::uint bucket = 0;
    // assuming most feature IDs are less than sizeof (U) * 8
    Iterator it, jt;
    for (it = jt = last; it > first;) {
      const value_type &fi = *--it;
      if (fi < sizeof (ny::uint) * 8)
        bucket |= (1 << fi);
      else // leave unsorted (large) feature IDs on tail
        *--jt = fi;
    }
    // do bucket sort
    it = first;
    // pick input numbers by twiddling bits
    //   cf. http://www-cs-faculty.stanford.edu/~uno/fasc1a.ps.gz (p. 8)
    while (bucket) {
      pecco::byte_4 b (static_cast <float> (bucket & - bucket)); // pull rightmost 1
      *it = (b.u >> 23) - 127; ++it; // pick it via floating number
      bucket &= (bucket - 1); // unset rightmost 1
    }
    // do insertion sort for the tail
    _insertion_sort (it, last);
  }
  // read a vector of feature indices
  template <typename T>
  void ClassifierBase <T>::_convertFvstr2Fv (char * &p, const char * const p_end,
                                             ny::fv_t &fv, const bool flag) {
    // convert an example in SVM-light format into a feature vector
    if (! flag) fv.reserve (_maf);
    while (p != p_end) {
      ++p;
      ny::uint fn = ny::strton <ny::uint> (p, &p);
      if (flag) {
        if (_fncnt.find (fn) == _fncnt.end ())
          _fncnt.insert (counter_t::value_type (fn, 0));
        fv.push_back (fn);
      } else {
        if (ny::uint fi = _fn2fi[fn < _fn2fi.size () ? fn : 0])
          fv.push_back (fi);
      }
      while (p != p_end && ! std::isspace (*p)) ++p;
    }
    if (! flag) std::sort (fv.begin (), fv.end ());
  }
  template <typename T>
  void ClassifierBase <T>::_convertFv2Fv (ny::fv_t &fnv) const {
    // convert feature numbers to feature indices
    const ny::fv_t::iterator it_end = fnv.end ();
    ny::fv_t::iterator jt = fnv.begin ();
    for (ny::fv_it it = jt; it != it_end; ++it)
      if (ny::uint fi = _fn2fi[*it < _fn2fi.size () ? *it : 0])
        { *jt = fi; ++jt; }
    fnv.erase (jt, it_end);
  }
  //
  template <typename T>
  void ClassifierBase <T>::_sortFv (ny::fv_t &fnv) const {
#ifdef SHIFT_LEFT_UNSIGNED_USES_SHL
    // for feature vectors following Zipf
    _bucket_sort (fnv.begin (), fnv.end ());
#else
    // for general feature vectors
    std::sort (fnv.begin (), fnv.end ());
    // for general short feature vectors
    // _insertion_sort (fnv.begin (), fnv.end ());
    // for general long feacture vectors
    // _radix_sort <SORT_MAX> (fnv.begin (), fnv.end ());
#endif
  }
  // * a function that estimates the additional cost
  //   when eliminating nodes
  template <typename T>
  inline size_t ClassifierBase <T>::_cost_fun (size_t m, size_t n) {
    // by frequency * # of features
    switch (_d) {
      case 1: return 0;
      case 2: return ((n*n - n) - (m*m - m)) / 2;
      case 3: return ((n*n*n - n) - (m*m*m - m)) / 6;
      case 4: return ((n*n*n*n - 2*n*n*n + 11*n*n - 10*n)
                      -(m*m*m*m - 2*m*m*m + 11*m*m - 10*m)) / 24;
      default: ny::print_err ("please add case statement."); return 0; // dummy;
    }
  }
  // assign younger feature indices to important feature numbers
  template <typename T>
  bool ClassifierBase <T>::_packingFeatures (std::vector <ny::fv_t> &fvv) {
    if (_opt.verbose > 0)
      std::fprintf (stderr, "packing feature id..");
    // reorder features according to their frequency in training data
    if (_opt.train) {
      FILE* reader =  std::fopen (_opt.train, "r");
      if (! reader) ny::print_err (HERE "no such file: %s\n", _opt.train);
      _nt         = 0;
      char * line = 0;
      size_t read = 0;
      while (ny::getLine (reader, line, read)) {
        if (*line != '\n') {
          // assume primitive feature vectors
          char * p (line), * const p_end (line + read - 1);
          while (p != p_end && ! std::isspace (*p)) ++p;
          while (p != p_end) {
            ++p;
            const ny::uint fn = ny::strton <ny::uint> (p, &p); // *p = ':'
            counter_t::iterator it = _fncnt.find (fn);
            if (it != _fncnt.end ()) ++it->second;
            while (p != p_end && ! std::isspace (*p)) ++p;
          }
          ++_nt;
        }
      }
      std::fclose (reader);
    } else {
      // reorder features according to their frequency in support vectors / conjunctive features
      for (std::vector <ny::fv_t>::const_iterator it = fvv.begin ();
           it != fvv.end (); ++it) // for each unit
        for (ny::fv_it fit = it->begin (); fit != it->end (); ++fit) {
          counter_t::iterator jt = _fncnt.find (*fit); // for each feature
          if (jt != _fncnt.end ()) ++jt->second;
        }
    }
    // building feature mapping
    ny::counter <ny::uint>::type fn_counter;
    fn_counter.reserve (_fncnt.size ());
    ny::uint fn_max = 0;
    for (counter_t::const_iterator it = _fncnt.begin ();
         it != _fncnt.end (); ++it) {
      const ny::uint fn  = it->first;
      const ny::uint cnt = it->second;
      fn_max = std::max (fn, fn_max);
      fn_counter.push_back (ny::counter <ny::uint>::type::value_type (cnt, fn));
    }
    std::sort (fn_counter.begin (), fn_counter.end ());
    _fn2fi.resize (fn_max + 1, 0);
    _fi2fn.resize (_nf + 1, 0); // _nf
    ny::uint fi = 1;
    for (ny::counter <ny::uint>::type::reverse_iterator it = fn_counter.rbegin ();
         it != fn_counter.rend (); ++it) {
      const ny::uint fn = it->second;
      _fi2fn[fi] = fn;
      _fn2fi[fn] = fi;
      ++fi;
    }
    // rename features in conjunctive features / support vectors
    for (std::vector <ny::fv_t>::iterator it = fvv.begin ();
         it != fvv.end (); ++it) { // for each unit
      _convertFv2Fv (*it);
      _sortFv (*it); // for pkb; bug
    }
    if (_opt.verbose > 0)
      std::fprintf (stderr, "done.\n");
    return true;
  }
  template <typename T>
  template <int D, binary_t FLAG>
  inline void ClassifierBase <T>::_pkePseudoInnerLoop (double * score, ny::fv_it it, const ny::fv_it &beg, const ny::fv_it &end, const ny::uint pos) {
    for (; it != end; ++it) {
      ny::uint pos_ = pos;
      const ny::uint j = *it - 1;
      switch (D) {
        case 4: pos_ += j * (j - 1) * (j - 2) * (j - 3) / 24;
        case 3: pos_ += j * (j - 1) * (j - 2) / 6;
        case 2: pos_ += j * (j - 1) / 2;
        case 1: pos_ += j;
      }
      PROFILE (++_traverse);
      PROFILE (for (ny::uint i (0), n (pos_ * _nl); i < _nl; ++i)
                 if (_fw[n + i]) { ++_lookup; break; });
      addScore <FLAG> (score, pos_);
      _pkePseudoInnerLoop <D - 1, FLAG> (score, beg, beg, it, pos_ + 1);
     }
  }
  // explicit specializations
#ifdef USE_KERNEL
  template <> template <>
  inline void ClassifierBase <kernel_model>::_pkePseudoInnerLoop <0, BINARY> (double *, ny::fv_it, const ny::fv_it &, const ny::fv_it &, const ny::uint) {}
  template <> template <>
  inline void ClassifierBase <kernel_model>::_pkePseudoInnerLoop <0, MULTI>  (double *, ny::fv_it, const ny::fv_it &, const ny::fv_it &, const ny::uint) {}
#endif
#ifdef USE_LINEAR
  template <> template <>
  inline void ClassifierBase <linear_model>::_pkePseudoInnerLoop <0, BINARY> (double *, ny::fv_it, const ny::fv_it &, const ny::fv_it &, const ny::uint) {}
  template <> template <>
  inline void ClassifierBase <linear_model>::_pkePseudoInnerLoop <0, MULTI>  (double *, ny::fv_it, const ny::fv_it &, const ny::fv_it &, const ny::uint) {}
#endif
  
  template <typename T>
  template <int D, binary_t FLAG>
  inline void ClassifierBase <T>::_pkeInnerLoop (double * score, ny::fv_it it, const ny::fv_it &beg, const ny::fv_it &end, const size_t pos) {
    byte_encoder encoder;
    for (; it != end; ++it) {
      size_t p (0), pos_ (pos);
      encoder.encode (*it);
      const int n = _ftrie.traverse (encoder.key (), pos_, p, encoder.len ());
      PROFILE (++_traverse);
      PROFILE (if (n >= 0) ++_lookup);
      // pke       : 0.0114 ms./classify (694.25668038/60766)
      if (n == -2) continue;
      if (n >=  0) addScore <FLAG> (score, n, _fw);
      _pkeInnerLoop <D - 1, FLAG> (score, beg, beg, it, pos_);
    }
  }
  // explicit specializations
#ifdef USE_KERNEL
  template <> template <>
  inline void ClassifierBase <kernel_model>::_pkeInnerLoop <0, BINARY> (double *, ny::fv_it, const ny::fv_it &, const ny::fv_it &, const size_t) {}
  template <> template <>
  inline void ClassifierBase <kernel_model>::_pkeInnerLoop <0, MULTI>  (double *, ny::fv_it, const ny::fv_it &, const ny::fv_it &, const size_t) {}
#endif
#ifdef USE_LINEAR
  template <> template <>
  inline void ClassifierBase <linear_model>::_pkeInnerLoop <0, BINARY> (double *, ny::fv_it, const ny::fv_it &, const ny::fv_it &, const size_t) {}
  template <> template <>
  inline void ClassifierBase <linear_model>::_pkeInnerLoop <0, MULTI>  (double *, ny::fv_it, const ny::fv_it &, const ny::fv_it &, const size_t) {}
#endif
  
  template <typename T>
  template <int D>
  void ClassifierBase <T>::_pkeClassify (double * score, ny::fv_it it, const ny::fv_it &beg, const ny::fv_it &end) {
    const ny::fv_it pend = std::lower_bound (it, end, 1 << PSEUDO_TRIE_N[_d]);
    PROFILE (_flen += end - beg);
    if (is_binary_classification ()) {
      _pkePseudoInnerLoop <D, BINARY> (score, it,   beg, pend, 0);
      _pkeInnerLoop       <D, BINARY> (score, pend, beg, end,  0);
    } else {
      _pkePseudoInnerLoop <D, MULTI>  (score, it,   beg, pend, 0);
      _pkeInnerLoop       <D, MULTI>  (score, pend, beg, end,  0);
    }
  }
  template <typename T>
  void ClassifierBase <T>::_pkeClassify (const ny::fv_t &fv, double * score, ny::fv_it it, const ny::fv_it &end) {
    const ny::fv_it beg  = fv.begin ();
    switch (_d) {
      case 1: _pkeClassify <1> (score, it, beg, end); break;
      case 2: _pkeClassify <2> (score, it, beg, end); break;
      case 3: _pkeClassify <3> (score, it, beg, end); break;
      case 4: _pkeClassify <4> (score, it, beg, end); break;
      default: ny::print_err ("please add case statement.");
    }
  }
  template <typename T>
  template <binary_t FLAG>
  void ClassifierBase <T>::_fstClassify (double * score, ny::fv_it &cit, const ny::fv_it &end) {
    PROFILE (const ny::fv_it beg = cit);
    size_t   pos  = 0;
    ny::uint prev = 0;
    byte_encoder encoder;
    for (; cit != end; ++cit) {
      size_t p = 0;
      encoder.encode (*cit - prev);
      PROFILE (++_traverse_sp);
      const int n = _fstrie.traverse (encoder.key (), pos, p, encoder.len ());
      prev = *cit;
      if (n < 0) break;
      addScore <FLAG> (score, n, _fsw);
    }
    PROFILE (_all +=
             _cost_fun (0, static_cast <size_t> (std::distance (beg, end))));
    PROFILE (_hit +=
             _cost_fun (static_cast <size_t> (std::distance (beg, cit)),
                        static_cast <size_t> (std::distance (beg, end))));
    PROFILE (++_lookup_sp);
  }
  template <typename T>
  void ClassifierBase <T>::_fstClassify (const ny::fv_t &fv, double * score) {
    ny::fv_it cit = fv.begin ();
    if (is_binary_classification ())
      _fstClassify <BINARY> (score, cit, fv.end ());
    else
      _fstClassify <MULTI>  (score, cit, fv.end ());
    _baseClassify (fv, score, cit, fv.end ()); // splitSVM + FST
    // return _pkeClassify (fv, score, cit, end); // used in EMNLP
  }
  // fstrie construction
  template <typename T>
  bool ClassifierBase <T>::_setFStrie () {
    const bool abuse = abuse_trie ();
    std::string fstrie_fn_ (_opt.event + (abuse ? "." : ".n") + ny::TRIE_SUFFIX);
    std::string fsw_fn_    (_opt.event + ".weight");
    std::ostringstream ss;
#ifdef USE_MODEL_SUFFIX
    ss << ".-" << _opt.fst_factor;
#endif
    std::string fstrie_fn  (fstrie_fn_ + ss.str ());
    std::string fsw_fn     (fsw_fn_    + ss.str ());
    // check the update time of events
    if (_opt.verbose > 0)
      std::fprintf (stderr, "loading fstrie..");
    if (! _opt.force && newer (fstrie_fn.c_str (), _opt.event.c_str ()) &&
        _fstrie.open (fstrie_fn.c_str ()) == 0) {
      if (! abuse) { // load the pre-computed weights
        FILE* reader = std::fopen (fsw_fn.c_str (), "rb");
        if (! reader) ny::print_err (HERE "no such file: %s\n", fsw_fn.c_str ());
        if (std::fseek (reader, 0, SEEK_END) != 0) return -1;
        const size_t nfs = static_cast <size_t> (std::ftell (reader)) / (_nl * sizeof (ny::fl_t));
        if (std::fseek (reader, 0, SEEK_SET) != 0) return -1;
        _fsw = new ny::fl_t [_nl * nfs];
        std::fread (&_fsw[0], sizeof (ny::fl_t), _nl * nfs, reader);
        std::fclose (reader);
      }
      if (_opt.verbose > 0) std::fprintf (stderr, "done.\n");
    } else {
      if (_opt.verbose > 0) std::fprintf (stderr, "not found.\n");
      FILE* reader = std::fopen (_opt.event.c_str (), "r");
      if (! reader)
        ny::print_err (HERE "no such event file: %s\n", _opt.event.c_str ());
      if (_opt.verbose > 0)
        std::fprintf (stderr, "building fstrie from %s..", _opt.event.c_str ());
      size_t nt (0), nkey (0);
      std::set <FstKey*, FeatKeypLess> fst_key;
      FstKey q;
      std::vector <ny::uchar> s;
      std::vector <double> w (_nl, 0);
      char * line = 0;
      size_t read = 0;
      ny::fv_t fv;
      while (ny::getLine (reader, line, read)) {
        if (*line != '\0') {
          // skip label
          char* p (line), * p_end (line + read - 1);
          while (p != p_end && ! std::isspace (*p)) ++p;
          fv.clear ();
          _convertFvstr2Fv (p, p_end, fv, false);
          // precompute score
          if (KEY_SIZE * fv.size () >= s.size ())
            s.resize (KEY_SIZE * fv.size () + 1);
          ny::uint len (0), prev (0);
          byte_encoder encoder;
          for (ny::fv_t::iterator jt = fv.begin (); jt != fv.end (); ) {
            len += encoder.encode (*jt - prev, &s[len]);
            prev = *jt;
            s[len] = '\0'; q.key = &s[0]; q.len = len;
            std::set <FstKey*, FeatKeypLess>::iterator it
              = fst_key.lower_bound (&q);
            ++jt;
            if (it == fst_key.end () || FeatKeypLess () (&q, *it)) {
              ++nkey;
              std::fill (&w[0], &w[_nl], 0);
              _baseClassify (fv, &w[0], jt - 1, jt);
              it = fst_key.insert (it, new FstKey (&s[0], 0, len, _nl));
              std::copy (w.begin (), w.end (), (*it)->cont);
              const size_t j
                = static_cast <size_t> (std::distance (fv.begin (), jt));
              (*it)->weight = _cost_fun (j - 1, j);
              if (jt == fv.end ()) (*it)->leaf = true;  // leaf (tentative)
            } else {
              if (jt != fv.end ()) (*it)->leaf = false; // this is not leaf
            }
            (*it)->count += 1;
          }
        }
        if (++nt % 1000 == 0 && _opt.verbose > 0)
          std::fprintf (stderr, "\r processing %ld events => %ld feature sequences",
                        nt, nkey);
      }
      q.key = 0; q.len = 0;
      std::fclose (reader);
      if (_opt.verbose > 0)
        std::fprintf (stderr, "\r processing %ld events => %ld feature sequences\n",
                      nt, nkey);
      std::set <FstKey*, FstKeypLess> fst_leaf;
      for (std::set <FstKey*, FstKeypLess>::const_iterator it = fst_key.begin ();
           it != fst_key.end (); ++it)
        if ((*it)->leaf) fst_leaf.insert (*it);
      if (_opt.verbose > 0)
        std::fprintf (stderr, " # leaf: %ld\n", fst_leaf.size ());
      std::vector <const char*> str;
      std::vector <size_t>      len;
      std::vector <int>         val;
      str.reserve (nkey);
      len.reserve (nkey);
      val.reserve (nkey);
      if (! abuse) _fsw = new ny::fl_t [_nl * nkey];
      size_t nkey_small = nkey;
      for (size_t j = 0; nkey_small > 0; ++j) {
        // sort dictionary order of key strings
        while (fst_key.size () > nkey_small) {
          // should memory release
          std::set <FstKey*, FstKeypLess>::iterator jt = fst_leaf.begin ();
          FstKey * const p = *jt;
          std::set <FstKey*, FeatKeypLess>::iterator it = fst_key.find (p);
          // WARNING: a leaf node is sorted by its impact; look left/right
          // add its longest prefix if there are no sibling
          // 1 2 3     l
          // 1 2 3 4   p
          // 1 2 3 5   r
          // 1 2 3 5 6
          if (it != fst_key.begin ()) { // not a root node
            FstKey * const l = *(--it); ++it;
            if (l->is_prefix (p)) { // prefix (next leaf)
              if (++it == fst_key.end () || ! l->is_prefix (*it)) // no sibling
                { fst_leaf.insert (jt, l); }
              --it;
            }
          }
          fst_leaf.erase (jt);
          fst_key.erase (it);
          delete p;
        }
        ny::uint i = 0;
        std::set <FstKey*, FeatKeypLess>::const_iterator it = fst_key.begin ();
        typedef std::map <std::vector <double>, size_t> w2id_t;
        w2id_t w2id;
        size_t uniq = 0;
        for (; it != fst_key.end (); ++it) {
          FstKey* p = *it;
          str.push_back (reinterpret_cast <char*> (p->key));
          len.push_back (p->len);
          ny::fl_t * wv = p->cont;
          if (abuse) {
            union byte_4 b4 (static_cast <float> (*wv));
            b4.u >>= 1;
            val.push_back (b4.i);
          } else {
            w.assign (&wv[0], &wv[_nl]);
            std::pair <w2id_t::iterator, bool> itb
              = w2id.insert (w2id_t::value_type (w, uniq * _nl));
            if (itb.second) {
              for (ny::uint li = 0; li < _nl; ++li)
                _fsw[uniq * _nl + li] = wv[li]; // 0:(li=0, i=0) 1:(li=1, i=0)
              ++uniq;
            }
            val.push_back (static_cast <int> (itb.first->second)); // i*_nl
          }
          ++i;
        }
        if (j == _opt.fst_factor || _opt.fst_verbose) {
          ss.str ("");
#ifdef USE_MODEL_SUFFIX
          ss << ".-" << j;
#endif
          if (! abuse) {
            std::string fsw_fn_j  (fsw_fn_ + ss.str ());
            FILE* writer = std::fopen (fsw_fn_j.c_str (), "wb");
            std::fwrite (&_fsw[0], sizeof (ny::fl_t), uniq *_nl, writer);
            std::fclose (writer);
          }
          std::string fstrie_fn_j (fstrie_fn_ + ss.str ());
          ss.str ("");
          ss << "feature sequence trie (with 2^-" << j << " feature sequences)";
          build_trie (&_fstrie, ss.str (), fstrie_fn_j, str, len, val, _opt.verbose > 0);
          _fstrie.clear (); // if (j == _opt.fs_factor)
        }
        str.clear  ();
        len.clear  ();
        val.clear  ();
        nkey_small >>= 1;
        if (_opt.fst_factor == j && ! _opt.fst_verbose) break;
      }
      // try reload
      if (_fstrie.open (fstrie_fn.c_str ()) != 0)
        ny::print_err (HERE "no such double array: %s\n", fstrie_fn.c_str ());
      if (!abuse) {
        delete [] _fsw;
        // load computed score
        reader = std::fopen (fsw_fn.c_str (), "rb");
        if (std::fseek (reader, 0, SEEK_END) != 0) return -1;
        const size_t nfs
          = static_cast <size_t> (std::ftell (reader)) / (_nl * sizeof (ny::fl_t));
        if (std::fseek (reader, 0, SEEK_SET) != 0) return -1;
        _fsw = new ny::fl_t [_nl * nfs];
        std::fread (&_fsw[0], sizeof (ny::fl_t), _nl * nfs, reader);
        std::fclose (reader);
      }
      for (std::set <FstKey*>::const_iterator it = fst_key.begin ();
           it != fst_key.end (); ++it)
        delete *it;
      if (_opt.verbose > 0) std::fprintf (stderr, "done.\n");
    }
    return true;
  }
  template <typename T>
  void ClassifierBase <T>::classify (char * p, double * score) {
    if (! p) { std::fill (&score[0], &score[_nl], 0); return; }
    _fv.clear ();
    while (1) {
      _fv.push_back (ny::strton <ny::uint> (p, &p)); // *p = ':'
      while (*p != '\0' && *p != ' ' && *p != '\t') ++p;
      if (*p) ++p; else break;
    }
    classify (_fv, &score[0]);
  }
  template <typename T>
  void ClassifierBase <T>::batch () { // batch classification
    if (_opt.verbose > 0) std::fprintf (stderr, "processing examples..");
    FILE * reader = _opt.test ? std::fopen (_opt.test, "r") : stdin;
    if (! reader) ny::print_err (HERE "no such file: %s\n", _opt.test);
    if (reader == stdin)
      std::fprintf (stderr, "(input: STDIN)\n");

    char * line = 0;
    size_t read = 0;
    double * score = new double[_nl];
    ny::uint pp (0), np (0), pn (0), nn (0);
    while (ny::getLine (reader, line, read)) {
      if (*line != '\0') {
        char * p (line), * const p_end (line + read - 1), * label (p);
        while (p != p_end && *p != ' ' && *p != '\t') ++p; *p = '\0';
        classify (++p, &score[0]);
        const ny::uint li = getLabel (score);
        if (_opt.verbose > 1) printScore (li, score);
        if (std::strcmp (label, _li2l[li]) == 0)
          if (li == _tli) ++pp; else ++nn;
        else
          if (li == _tli) ++np; else ++pn;
      }
    }
    delete [] score;
    if (reader != stdin) std::fclose (reader);
    if (_opt.verbose > 0) std::fprintf (stderr, "done.\n");
    std::fprintf (stderr, "acc. %.4f (pp %d) (pn %d) (np %d) (nn %d)\n",
                  (pp + nn) * 1.0 / (pp + pn + np + nn), pp, pn, np, nn);
    printStat ();
  }
  template <typename T>
  void ClassifierBase <T>::printStat () { // print classifier statistics
#ifdef USE_TIMER
      _timer_pool.print ();
#endif
#ifdef USE_PROFILING
      std::fprintf (stderr, "# active primitive features: %ld\n",_flen);
      if (_opt.algo == PKE || _opt.algo == FST) {
        std::fprintf (stderr, "weight lookup (succeeded): %ld\n",  _lookup);
        std::fprintf (stderr, "weight lookup: %ld\n", _traverse);
        if (_opt.algo == FST) {
          std::fprintf (stderr, "weight lookup (split): %ld\n", _lookup_sp);
          std::fprintf (stderr, "weight traverse (split): %ld\n", _traverse_sp);
          std::fprintf (stderr, "node reduction: %f (%ld / %ld)\n",
                        static_cast <ny::fl_t> (_hit) / static_cast <ny::fl_t> (_all), _hit, _all);
        }
      }
#endif
  }
  // explicit specialization
#ifdef USE_KERNEL
  template class ClassifierBase <kernel_model>;
#endif
#ifdef USE_LINEAR
  template class ClassifierBase <linear_model>;
#endif
}
// use pseudo array
