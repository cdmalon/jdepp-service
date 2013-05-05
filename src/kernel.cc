// pecco -- please enjoy classification with conjunctive features
//  $Id: kernel.cc 827 2012-05-10 14:12:37Z ynaga $
// Copyright (c) 2008-2012 Naoki Yoshinaga <ynaga@tkl.iis.u-tokyo.ac.jp>
#include "kernel.h"

namespace pecco {
  
  bool kernel_model::_splitInit () {
    if (_d == 1 && _opt.algo == FST)
      std::fprintf (stderr, "NOTE: FST [-t 1] is useless for [-d 1].\n");
    double ec (0), ec_c (0), ec_r (0);
    for (ny::uint fi = 1; fi <= _nf; ++fi) {
      const ny::uint fn    = _fi2fn[fi];
      const double &size = static_cast <double> (_fncnt[fn]);
      double ratio = size / (static_cast <double> (_opt.train ? _nt :_nunit));
      ratio *= static_cast <double> (100.0);
      if (ratio >= _fratio) ec_c += size; else ec_r += size;
      ec += size;
    }
    return true;
  }
  bool kernel_model::_precomputeKernel () {
    if (_opt.verbose > 0) std::fprintf (stderr, "precomputing kernel..");
    _polyk = new double[_maf + 1];
    if (_fratio)
      _spolyk = new double[_maf + 1];
    for (ny::uint i = 0; i <= _maf; ++i) { // bug; i -> dot product rare
      _polyk[i] = static_cast <double>
                  (std::pow (_s * i + _r, static_cast <double> (_d)));
      if (_fratio)
        _spolyk[i] = static_cast <double>
                     (std::pow (_s * (i + 1) + _r, static_cast <double> (_d))
                      - std::pow (_s * i + _r, static_cast <double> (_d)));
    }
    if (_opt.verbose > 0) std::fprintf (stderr, "done.\n");
    return true;
  }
  void kernel_model::_setup_binary_labels () {
    _nl = 1;
    char * label_pos = new char[3]; std::strcpy (label_pos, "+1");
    char * label_neg = new char[3]; std::strcpy (label_neg, "-1");
    _li2l.push_back (label_pos); _li2l.push_back (label_neg);
    _l2li.insert (lmap::value_type (label_pos, 0));
    _l2li.insert (lmap::value_type (label_neg, 1));
  }
  bool kernel_model::load (const char * model) {
    TIMER (_model_t->startTimer ());
    const std::string model_bin (std::string (model) + ".bin");

    if (! _opt.force && newer (model_bin.c_str (), model)) {
      if (_opt.verbose > 0)
        std::fprintf (stderr, "loading compiled model parameters..");
      FILE* reader = std::fopen (model_bin.c_str (), "rb");
      // model parameters
      std::fread (&_d,     sizeof (ny::uint),   1, reader);
      std::fread (&_nl,    sizeof (ny::uint),   1, reader);
      std::fread (&_nf,    sizeof (ny::uint),   1, reader);
      std::fread (&_nunit, sizeof (ny::uint),   1, reader);
      std::fread (&_maf,   sizeof (ny::uint),   1, reader);
      std::fread (&_r,     sizeof (double),     1, reader);
      std::fread (&_s,     sizeof (double),     1, reader);
      _b  = new double[_nl];
      _m0 = new double[_nl];
      std::fread (_b,      sizeof (double),   _nl, reader);
      std::fread (_m0,     sizeof (double),   _nl, reader);
      // label map
      if (_nl >= 2)
        for (ny::uint li = 0; li < _nl; ++li) {
          ny::uint len = 0;
          std::fread (&len, sizeof (ny::uint), 1, reader);
          char * p = new char[len + 1]; p[len] = '\0';
          std::fread (p, sizeof (char), len, reader);
          _li2l.push_back (p);
          _l2li.insert (lmap::value_type (p, li));
          if (std::strcmp (p, "+1") == 0) _tli = li;
        }
      else
        _setup_binary_labels ();
      // feature map
      _fi2fn.reserve (_nf + 1);
      _fi2fn.push_back (0);
      double ratio = 0;
      for (ny::uint fn = 0;
           _fi2fn.size () <= _nf &&
           std::fread (&fn, sizeof (ny::uint), 1, reader) &&
           std::fread (&ratio, sizeof (double), 1, reader); ) {
        if (fn >= _fn2fi.size ()) _fn2fi.resize (fn + 1, 0);
        _fn2fi[fn] = static_cast <ny::uint> (_fi2fn.size ());
        _fi2fn.push_back (fn);
        if (ratio >= _fratio) _f_r = static_cast <ny::uint> (_fi2fn.size ());
      }
      if (_opt.algo != PKI) _sv.reserve (_nunit);
      // if (_opt.algo == PKI)
      _alph.reserve (_nunit * _nl);
      if (_opt.algo == PKI || _fratio)
        _f2ss.resize (_nf + 1); // _nf
      std::vector <ny::fl_t> alpha (_nl, 0);
      ny::uint len   = 0;
      ny::uint i     = 0;
      std::vector <ny::uchar> s;
      ny::fv_t fv;
      std::ostringstream ss;
#ifdef USE_MODEL_SUFFIX
      ss << model << ".r" << _opt.fratio << ".s" << _opt.sigma << ".af";
#else
      ss << model << ".af";
#endif
      // we should read data when PKE with unknown sigma
      if (_opt.algo == PKI || _fratio || _opt.force || ! newer (ss.str ().c_str (), model))
        while (std::fread (&alpha[0], sizeof (ny::fl_t), _nl, reader) &&
               std::fread (&len,      sizeof (ny::uint), 1,   reader)) {
          // read
          if (s.size () < len) s.resize (len);
          std::fread (&s[0], sizeof (ny::uchar), len, reader);
          fv.clear ();
          ny::uint r (0), b (0), prev (0);
          for (ny::uint j = 0; j < len; ++j) {
            r += ((static_cast <ny::uint> (s[j]) & 0x7f) << b);
            if (s[j] & 0x80) b += 7; else { fv.push_back (prev += r); r = b = 0;}
          }
          for (ny::uint li = 0; li < _nl; ++li) _alph.push_back (alpha[li]);
          if (_opt.algo == PKI) {
            for (ny::fv_it it = fv.begin (); it != fv.end (); ++it)
              _f2ss[*it].push_back (i);
          } else {
            _sv.push_back (fv);
            for (ny::fv_it it = fv.begin (); it != fv.end (); ++it)
              if (*it >= _f_r) _f2ss[*it].push_back (i);
          }
          if (++i >= _nunit) break;
        }
      std::fclose (reader);
      if (_opt.verbose > 0) std::fprintf (stderr, "done.\n");
    } else {
      if (_opt.verbose > 0)
        std::fprintf (stderr, "loading/compiling model parameters..");
      FILE * reader = std::fopen (model, "r");
      if (! reader)
        ny::print_err (HERE "no such file: %s\n", model);
      bool   header = true;
      char * line   = 0;
      size_t read   = 0;
      ny::fv_t fnv;
      while (ny::getLine (reader, line, read)) { // parse
        if (*line != '\0') {
          if (! header) {
            char * p (line);
            for (ny::uint li = 0; li < _nl; ++li) {
              const ny::fl_t alpha = ny::strton <ny::fl_t> (p, &p); ++p;
              _m0[li] += alpha;
              _alph.push_back (alpha);
            }
            --p;
            fnv.clear ();
            _convertFvstr2Fv (p, line + read - 1, fnv, true);
            _maf = std::max (_maf, static_cast <ny::uint> (fnv.size ()));
            _sv.push_back (fnv);
            ++_nunit;
          } else if (std::strstr (line, " # kernel type") != NULL) {
            if (*line != '1') ny::print_err (HERE "unsupported kernel used.\n");
          } else if (std::strstr (line, " # kernel parameter") != NULL) {
            switch (line[read - 2]) {
              case 'r': _r = ny::strton <double>   (line, NULL); break;
              case 's': _s = ny::strton <double>   (line, NULL); break;
              case 'd': _d = ny::strton <ny::uint> (line, NULL); break;
            }
          } else if (std::strstr (line, " # labels: ") != NULL) { // multiclass
            char * p (line), * const p_end (line + read - 1);
            _nl = ny::strton <ny::uint> (p, &p); p += 10;
            while (++p) {
              char * ys = p; while (p != p_end && *p != ' ') ++p; *p = '\0';
              char * copy = new char[p - ys + 1];
              std::strcpy (copy, ys);
              _li2l.push_back (copy);
              _l2li.insert (lmap::value_type (copy, _nl));
              if (p == p_end) break;
            }
          } else if (std::strstr (line, " # threshold b") != NULL) {
            if (_l2li.empty ()) _setup_binary_labels ();
            _b  = new double [_nl] ();
            _m0 = new double [_nl] ();
            char * p = line;
            for (ny::uint li = 0; li < _nl; ++li)
              _b[li] = ny::strton <double> (p, &p), ++p;
            header = false;
          }
        }
      }
      std::fclose (reader);
      if (_opt.verbose > 0) std::fprintf (stderr, "done.\n");
      // Zipf-aware feature indexing
      _nf = static_cast <ny::uint> (_fncnt.size ());
      _packingFeatures (_sv);
      // report splitSVM status
      if (_fratio) _splitInit ();
      // output map
      FILE* writer = std::fopen (model_bin.c_str (), "wb");
      std::fwrite (&_d,     sizeof (ny::uint), 1, writer);
      std::fwrite (&_nl,    sizeof (ny::uint), 1, writer);
      std::fwrite (&_nf,    sizeof (ny::uint), 1, writer);
      std::fwrite (&_nunit, sizeof (ny::uint), 1, writer);
      std::fwrite (&_maf,   sizeof (ny::uint), 1, writer);
      std::fwrite (&_r,     sizeof (double),   1, writer);
      std::fwrite (&_s,     sizeof (double),   1, writer);
      std::fwrite (_b,      sizeof (double), _nl, writer);
      std::fwrite (_m0,     sizeof (double), _nl, writer);
      if (_nl >= 2)
        for (ny::uint li = 0; li < _nl; ++li) {
          const char * p = _li2l[li];
          ny::uint len = static_cast <ny::uint> (std::strlen (p));
          std::fwrite (&len, sizeof (ny::uint), 1,    writer);
          std::fwrite (p,    sizeof (char),     len, writer);
        }
      // common-to-rare ranking of features
      for (ny::uint fi = 1; fi <= _nf; ++fi) { // fi = 0; fi < _nf
        ny::uint fn = _fi2fn[fi];
        double ratio = static_cast <double> (_fncnt[fn]);
        ratio /= static_cast <double> (_opt.train ? _nt : _nunit);
        ratio *= 100.0;
        if (ratio >= _fratio) _f_r = fi + 1;
        std::fwrite (&fn,    sizeof (ny::uint), 1, writer);
        std::fwrite (&ratio, sizeof (double),   1, writer);
      }
      // support vectors
      _f2ss.resize (_nf + 1, ss_t ());
      std::vector <ny::uchar> s_;
      for (ny::uint i = 0; i < _nunit; ++i) {
        const ny::fv_t &s = _sv[i];
        // compile
        if (s_.size () < s.size () * KEY_SIZE)
          s_.resize (s.size () * KEY_SIZE);
        ny::uint len (0), prev (0);
        byte_encoder encoder;
        for (ny::fv_it it = s.begin (); it != s.end (); prev = *it, ++it)
          len += encoder.encode (*it - prev, &s_[len]);
        std::fwrite (&_alph[i * _nl], sizeof (ny::fl_t),  _nl, writer);
        std::fwrite (&len,            sizeof (ny::uint),    1, writer);
        std::fwrite (&s_[0],          sizeof (ny::uchar), len, writer);
        if (_opt.algo == PKI || _fratio)
          for (ny::fv_it it = _sv[i].begin (); it != _sv[i].end (); ++it)
            if (_opt.algo == PKI || *it >= _f_r)
              _f2ss[*it].push_back (i);
      }
      if (_opt.algo == PKI)
        std::vector <ny::fv_t> ().swap (_sv);
      _fncnt.clear ();
      std::fclose (writer);
    }
#ifdef ABUSE_TRIE
    if (! is_binary_classification ())
      ny::print_err (HERE "ABUSE_TRIE assumes a binary label.\n");
#endif
    _fv.reserve (_maf);
    // calculate dot product
    if (_opt.algo == PKI || _fratio)
      _precomputeKernel ();
    // used to indicate active features
    if (_opt.algo == PKI)
      _dot = new ny::uint[_nunit];
    else
      _fbit.resize (_nf + 1, false);
    // kernel expansion
    _nf_cut = _nf;
    if (_opt.algo == PKE || _opt.algo == FST) {
      _setPKEcoeff ();
      _setFtrie ();
    }
    if (_sigma || _f_r - 1 < _nf) // _f_r != _nf
      if (_opt.algo != PKI && _sigma && _opt.verbose > 1) {
        std::fprintf (stderr, "NOTE (approximated computation);");
        std::fprintf (stderr, " %d/%d features used, PKE sigma=%s\n",
                      _nf_cut, _nf, _opt.sigma); // _f_r - 1
      }
    // fstrie construction
    if (_opt.algo == FST) _setFStrie ();
    TIMER (_model_t->stopTimer ());
    if (_opt.verbose > 0) printParam ();
    if (_fratio)
      std::fill (_fbit.begin () + 1, _fbit.end (), false); // reuse it
    return true;
  }
  bool kernel_model::_setPKEcoeff () {
    // calculate weight coefficients for conjunctive features
    switch (_d) {
      case 1:
        _coeff[0] = _r;
        _coeff[1] = _s; break;
      case 2:
        _coeff[0] = _r*_r;
        _coeff[1] = _s*_s+2*_r*_s;
        _coeff[2] =  2*_s*_s; break;
      case 3:
        _coeff[0] = _r*_r*_r;
        _coeff[1] = _s*_s*_s + 3*_r*_s*_s + 3*_r*_r*_s;
        _coeff[2] =  6*_s*_s*_s + 6*_r*_s*_s;
        _coeff[3] =  6*_s*_s*_s; break;
      case 4:
        _coeff[0] = _r*_r*_r*_r;
        _coeff[1] = _s*_s*_s*_s + 4*_r*_s*_s*_s + 6*_r*_r*_s*_s + 4*_r*_r*_r*_s;
        _coeff[2] = 14*_s*_s*_s*_s + 24*_r*_s*_s*_s  + 12*_r*_r*_s*_s;
        _coeff[3] = 36*_s*_s*_s*_s + 24*_r*_s*_s*_s;
        _coeff[4] = 24*_s*_s*_s*_s; break;
      default: ny::print_err ("please add case statement.");
    }
    _max_coeff = *std::max_element (&_coeff[0], &_coeff[_d+1]);
    return true;
  }
  // implementation of kernel expansion (Kudo & Matsumoto, 2003) for multi-class
  void kernel_model::_pkePrefixSpan (ny::fv_t &fc, std::vector <ny::fl_t> &fw,
                                     const std::vector <std::pair <ny::uint, int> > &proj,
                                     std::vector <FeatKey*> &pke_key) {
    // pseudo projection to mine conjunctive features
    typedef std::vector <std::pair <ny::uint, int> > proj_t;
    typedef ny::map <ny::uint, proj_t>::type proj_map_t;
    proj_map_t proj_map;
    for (proj_t::const_iterator pit = proj.begin (); pit != proj.end (); ++pit) {
      const ny::uint i   = pit->first;
      const ny::fv_t &fv = _sv[i];
      for (int j = pit->second; j >= 0; --j)
        proj_map[fv[static_cast <size_t> (j)]].push_back (proj_t::value_type (i, j - 1));
    }
    ny::uchar s[KEY_SIZE * MAX_KERNEL_DEGREE]; // need to new if copied?
    std::vector <double> mu_pos (_nl), mu_neg (_nl);
    std::vector <ny::fl_t> w (_nl);
    for (proj_map_t::const_iterator pmit = proj_map.begin ();
         pmit != proj_map.end (); ++pmit) {
      // extend span
      const ny::uint fi         = pmit->first;
      const proj_t  &proj_child = pmit->second;
      std::fill (&mu_pos[0], &mu_pos[_nl], 0);
      std::fill (&mu_neg[0], &mu_neg[_nl], 0);
      std::fill (&w[0],      &w[_nl],      0);
      fc.push_back (fi);
      bool terminate = (fc.size () >= _d || proj_child.empty ());
      for (proj_t::const_iterator pit = proj_child.begin ();
           pit != proj_child.end (); ++pit) {
        const ny::fl_t * alpha = &_alph[pit->first * _nl];
        for (ny::uint li = 0; li < _nl; ++li)
          w[li] += static_cast <ny::fl_t> (alpha[li] * _coeff[fc.size ()]);
        if (! terminate) {
          for (ny::uint li = 0; li < _nl; ++li)
            if (alpha[li] > 0)
              mu_pos[li] += alpha[li] * _max_coeff;
            else
              mu_neg[li] += alpha[li] * _max_coeff;
        }
      }
      bool terminate_ = true;
      bool exceed     = false;
      for (ny::uint li = 0; li < _nl; ++li) { // pruning
        terminate_ &= (mu_neg[li] > _sigma_neg[li] && mu_pos[li] < _sigma_pos[li]);
        exceed    |= w[li] <= _sigma_neg[li] || w[li] >= _sigma_pos[li];
      }
      terminate |= terminate_; // the extending features will be useless
      if (exceed) { // approximation
        ++_exnunit;
        if (fc.front () >> PSEUDO_TRIE_N[_d]) {
          ny::uint len = 0;
          byte_encoder encoder;
          for (ny::fv_it it = fc.begin (); it != fc.end (); ++it) {
            _fbit[*it] = true;
            len += encoder.encode (*it, s + len);
          }
          // we remain the weights when a weight for one label exceeds threshold
          // because it does not take extra memory in current implementation
          pke_key.push_back (new FeatKey (s, &w[0], len, _nl));
        } else {
          ny::uint pos = 0;
          for (ny::uint i = 0; i < fc.size (); ++i) {
            _fbit[fc[i]] = true;
            const ny::uint j = fc[i] - 1;
            switch (_d - i) {
              case 4: pos += j * (j - 1) * (j - 2) * (j - 3) / 24;
              case 3: pos += j * (j - 1) * (j - 2) / 6;
              case 2: pos += j * (j - 1) / 2;
              case 1: pos += j;
            }
            if (i > 0) ++pos;
          }
          for (ny::uint li = 0; li < _nl; ++li)
            fw[pos * _nl + li] = w[li];
        }
      }
      if (! terminate) _pkePrefixSpan (fc, fw, proj_child, pke_key);
      if (fc.size () == 1)
        if (_opt.verbose > 0)
          std::fprintf (stderr,
                        "\r processed %d features => %d conjunctive features",
                        ++_processed, _exnunit);
      fc.pop_back ();
    }
  }
  // polynomial kernel expanded
  bool kernel_model::_setFtrie () {
    std::ostringstream ss;
    ss << _opt.model;
#ifdef USE_MODEL_SUFFIX
    ss << ".r" << _opt.fratio << ".s" << _opt.sigma;
#endif
#ifdef ABUSE_TRIE
    const std::string ftrie_fn (ss.str () + ".e"  + ny::TRIE_SUFFIX);
#else
    const std::string ftrie_fn (ss.str () + ".ne" + ny::TRIE_SUFFIX); // w/o alpha
#endif
    const std::string fw_fn (ss.str () + ".weight");
    const std::string af_fn (ss.str () + ".af");
    // set PKE coeff
    if (_opt.verbose > 0) std::fprintf (stderr, "loading feature weight trie..");
    if (! _opt.force && newer (af_fn.c_str (), _opt.model)) {
      // read pre-computed weights
      _ftrie.open (ftrie_fn.c_str ()); // may fail
      FILE* reader = std::fopen (af_fn.c_str (), "rb");
      if (! reader) ny::print_err (HERE "no such file: %s\n", af_fn.c_str ());
      std::fread (&_exnunit, sizeof (ny::uint), 1, reader);
      bool * fbit = new bool[_nf + 1];
      std::fread (&fbit[0], sizeof (bool), _nf + 1, reader); // _nf
      for (size_t i = 0; i <= _nf; ++i) _fbit[i] = fbit[i];
      delete [] fbit;
      std::fclose (reader);
      reader = std::fopen (fw_fn.c_str (), "rb");
      if (! reader) ny::print_err (HERE "no such file: %s\n", fw_fn.c_str ());
      if (std::fseek (reader, 0, SEEK_END) != 0) return -1;
      const size_t uniq
        = static_cast <size_t> (std::ftell (reader)) / (_nl * sizeof (ny::fl_t));
      if (std::fseek (reader, 0, SEEK_SET) != 0) return -1;
      _fw = new ny::fl_t [uniq * _nl];
      std::fread (&_fw[0], sizeof (ny::fl_t), uniq * _nl, reader);
      std::fclose (reader);
      if (_opt.verbose > 0) std::fprintf (stderr, "done.\n");
      _exnunit_cut = _exnunit;
    } else {
      if (_opt.verbose > 0) {
        std::fprintf (stderr, "not found.\n");
        std::fprintf (stderr, "calculating conjunctive feature weight..");
      }
      if (_f_r == 0) { // no common feature
        if (_opt.verbose > 0) std::fprintf (stderr, "skipped.\n");
        return false;
      }
      if (_opt.verbose > 0) std::fprintf (stderr, "\n");
      size_t uniq = 0;
      size_t N = std::min ((1U << PSEUDO_TRIE_N[_d]) - 1, _f_r);
      switch (_d) {
        case 4: uniq += N * (N - 1) * (N - 2) * (N - 3) / 24;
        case 3: uniq += N * (N - 1) * (N - 2) / 6;
        case 2: uniq += N * (N - 1) / 2;
        case 1: uniq += N;
      }
      std::vector <ny::fl_t> fw (uniq * _nl, 0); // initialized to 0;
      std::vector <size_t> pos_num (_nl, 0);
      std::vector <size_t> neg_num (_nl, 0);
      // pseudo projection (PrefixSpan)
      std::vector <std::pair <ny::uint, int> > proj;
      for (ny::uint i = 0; i < _nunit; ++i) {
        const int tail =
          static_cast <int> (std::lower_bound (_sv[i].begin (), _sv[i].end (), _f_r) - _sv[i].begin () - 1);
        proj.push_back (std::make_pair (i, tail));
        for (ny::uint li = 0; li < _nl; ++li)
          if (_alph[i * _nl + li] > 0) ++pos_num[li]; else ++neg_num[li];
      }
      _sigma_pos.resize (_nl, 0);
      _sigma_neg.resize (_nl, 0);
      for (ny::uint li = 0; li < _nl; ++li)  {
        _sigma_pos[li] =  _sigma * static_cast <double> (pos_num[li]) / _nunit;
        _sigma_neg[li] = -_sigma * static_cast <double> (neg_num[li]) / _nunit;
      }
      if (_opt.verbose > 0) {
        std::fprintf (stderr, " %d / %d features,\n", _f_r - 1, _nf);
        for (ny::uint li = 0; li < _nl; ++li)
          std::fprintf (stderr, "  _sigma_pos%d=%f, _sigma_neg%d=%f\n",
                        li, _sigma_pos[li], li, _sigma_neg[li]);
      }
      _processed = 0;
      ny::fv_t fc;
      std::vector <FeatKey *> pke_key;
      _pkePrefixSpan (fc, fw, proj, pke_key);
      _exnunit_cut = _exnunit;
      _nf_cut = 0;
#ifndef USE_ABUSE
      fw.resize ((uniq + _exnunit) * _nl, 0);
#endif
      for (ny::uint fi = 1; fi <= _nf; ++fi) if (_fbit[fi]) ++_nf_cut;
      if (_opt.verbose > 0) {
        std::fprintf (stderr, "\r processed %d features => %d conjunctive features\n",
                      _f_r - 1, _exnunit);
        std::fprintf (stderr, " # active features after thresholding %d => %d\n",
                      _nf, _nf_cut);
      }
      if (_opt.verbose > 0) std::fprintf (stderr, "done.\n");
      typedef std::map <std::vector <ny::fl_t>, size_t> w2id_t;
      w2id_t w2id;
      std::vector <ny::fl_t> w (_nl);
      for (size_t i = 0; i < uniq; ++i) {
        w.assign (&fw[i * _nl], &fw[(i + 1) * _nl]);
        w2id.insert (w2id_t::value_type (w, i * _nl));
      }
      // build feature weight trie
      if (! pke_key.empty ()) {
        std::vector <const char*> str; str.reserve (pke_key.size ());
        std::vector <size_t>      len; len.reserve (pke_key.size ());
        std::vector <int>         val; val.reserve (pke_key.size ());
        FeatKeypLess  featkeypless;
        std::sort (pke_key.begin (), pke_key.end (), featkeypless);
        std::vector <FeatKey*>::const_iterator dit = pke_key.begin ();
        for (; dit != pke_key.end (); ++dit) {
          str.push_back (reinterpret_cast <char*> ((*dit)->key));
          len.push_back ((*dit)->len);
          ny::fl_t * wv = (*dit)->cont;
#ifdef ABUSE_TRIE
          union byte_4 b4 (static_cast <float> (*wv)); b4.u >>= 1;
          val.push_back (b4.i);
#else
          // uniq identical weights
          // (minimize the double array with a slight cache-loss increase)
          w.assign (&wv[0], &wv[_nl]);
          std::pair <w2id_t::iterator, bool> itb
            = w2id.insert (w2id_t::value_type (w, uniq * _nl));
          if (itb.second) {
            for (ny::uint li = 0; li < _nl; ++li)
              fw[uniq * _nl + li] = wv[li]; // put weights for all labels
            ++uniq;
          }
          val.push_back (static_cast <int> (itb.first->second));
#endif
        }
        build_trie (&_ftrie, "feature weight trie", ftrie_fn, str, len, val,
                    _opt.verbose > 0);
        for (dit = pke_key.begin (); dit != pke_key.end (); ++dit) delete *dit;
      }
      FILE * writer = std::fopen (fw_fn.c_str (), "wb");
      std::fwrite (&fw[0], sizeof (ny::fl_t), uniq * _nl, writer);
      std::fclose (writer);
      _fw = new ny::fl_t [uniq * _nl];
      std::copy (&fw[0], &fw[uniq * _nl], _fw);
      //
      writer = std::fopen (af_fn.c_str (), "wb");
      std::fwrite (&_exnunit, sizeof (ny::uint), 1, writer); // omit it
      bool * fbit = new bool[_nf + 1];
      for (size_t i = 0; i <= _nf; ++i) fbit[i] = _fbit[i];
      std::fwrite (&fbit[0], sizeof (bool), _nf + 1, writer);
      delete [] fbit;
      std::fclose (writer);
      //
      if (! _fratio)
        std::vector <ny::fl_t> ().swap (_alph); // clear
    }
    _nf_cut = 0;
    for (ny::uint fi = 1; fi <= _nf; ++fi) if (_fbit[fi]) ++_nf_cut;
    for (ny::uint li = 0; li < _nl;  ++li) _m0[li] *= _coeff[0]; // bug fix
    return true;
  }
  template <binary_t FLAG>
  void kernel_model::_pkiClassify (double * score, const ny::fv_it &beg, const ny::fv_it &end) const {
    for (ny::fv_it it = beg; it != end; ++it) {
      const ss_t &ss = _f2ss[*it];
      for (ss_it st = ss.begin (); st != ss.end (); ++st) ++_dot[*st];
    }
    for (ny::uint i = 0; i < _nunit; ++i) // can be faster if alpha in sv
      addScore <FLAG> (score, i, _polyk[_dot[i]]);
  }
  void kernel_model::_pkiClassify (const ny::fv_t &fv, double * score) const {
    std::fill (&_dot[0], &_dot[_nunit], 0);
    if (is_binary_classification ())
      _pkiClassify <BINARY> (score, fv.begin (), fv.end ());
    else
      _pkiClassify <MULTI>  (score, fv.begin (), fv.end ());
  }
  // implementation of splitSVM (Goldberg & Elhadad, 2008)
  template <binary_t FLAG>
  void kernel_model::_splitClassify (double * score, ny::fv_it cit, const ny::fv_it &beg, const ny::fv_it &end) {
    // rare part by incremental split polynomial kernel
    for (ny::fv_it it = beg; it != cit; ++it) _fbit[*it] = true;
    for (ny::fv_it it = cit; it != end; ++it) {
      const ss_t &ss = _f2ss[*it];
      for (ss_it st = ss.begin (); st != ss.end (); ++st) {
        const ny::fv_t &s = _sv[*st];
        ny::uint dot_c = 0;
        for (ny::fv_it sit = s.begin (); sit != s.end (); ++sit)
          dot_c += _fbit[*sit];
        addScore <FLAG> (score, *st, _spolyk[dot_c]);
      }
      _fbit[*it] = true;
    }
    for (ny::fv_it it = beg; it != end; ++it) _fbit[*it] = false;
  }
  void kernel_model::_splitClassify (const ny::fv_t &fv, double * score, ny::fv_it cit, const ny::fv_it &end) {
    if (_f_r > _nf) {
      _pkeClassify (fv, score, cit, end); return;
    } else {
      const ny::fv_it rit = std::lower_bound (fv.begin (), fv.end (), _f_r);
      if (rit >= end)     { _pkeClassify (fv, score, cit, end); return;    }
      else if (cit < rit) { _pkeClassify (fv, score, cit, rit); cit = rit; }
    }
    PROFILE (_hit += cit - fv.begin ());
    if (is_binary_classification ())
      _splitClassify <BINARY> (score, cit, fv.begin (), end);
    else
      _splitClassify <MULTI>  (score, cit, fv.begin (), end);
  }
}
// score to weight
