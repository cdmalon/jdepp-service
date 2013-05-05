// pecco -- please enjoy classification with conjunctive features
//  $Id: linear.cc 827 2012-05-10 14:12:37Z ynaga $
// Copyright (c) 2008-2012 Naoki Yoshinaga <ynaga@tkl.iis.u-tokyo.ac.jp>
#include "linear.h"

namespace pecco {
  
  void linear_model::_convertCfstr2Cf (char * &p, ny::fv_t &fv, const bool flag) {
    // convert a feature number string to feature indices
    if (! flag) fv.reserve (_maf);
    while (*p != '\t') {
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
      while (*p != ':' && *p != '\t') ++p;
    }
    if (! flag) std::sort (fv.begin (), fv.end ());
  }
  bool linear_model::_setFtrie () {
    const std::string ftrie_fn (std::string (_opt.model) + ".ne" + ny::TRIE_SUFFIX);
    const std::string fw_fn    (std::string (_opt.model) + ".weight");

    if (_opt.verbose > 0)
      std::fprintf (stderr, "loading feature weight trie..");
    if (! _opt.force && newer (fw_fn.c_str (), _opt.model)) {
      // read feature weight trie
      _ftrie.open (ftrie_fn.c_str ()); // may fail
      FILE * reader = std::fopen (fw_fn.c_str (), "rb");
      if (! reader) ny::print_err (HERE "no such file: %s\n", fw_fn.c_str ());
      if (std::fseek (reader, 0, SEEK_END) != 0) return -1;
      const size_t uniq
        = static_cast <size_t> (std::ftell (reader)) / (_nl * sizeof (ny::fl_t));
      if (std::fseek (reader, 0, SEEK_SET) != 0) return -1;
      _fw = new ny::fl_t [_nl * uniq]; // one-dimentional array is faster
      std::fread (&_fw[0], sizeof (ny::fl_t), uniq * _nl, reader);
      std::fclose (reader);
      if (_opt.verbose > 0) std::fprintf (stderr, "done.\n");
    } else {
      if (_opt.verbose > 0) std::fprintf (stderr, "not found.\n");
      std::vector <FeatKey *> pke_key;
      size_t uniq = 0;
      size_t N = (1U << PSEUDO_TRIE_N[_d]) - 1;
      switch (_d) {
        case 4: uniq += N * (N - 1) * (N - 2) * (N - 3) / 24;
        case 3: uniq += N * (N - 1) * (N - 2) / 6;
        case 2: uniq += N * (N - 1) / 2;
        case 1: uniq += N;
      }
      std::vector <ny::fl_t> fw ((uniq + _nunit) * _nl, 0);
      ny::uchar s[KEY_SIZE * MAX_KERNEL_DEGREE];
      // should remember feature combination including rare feature
      for (ny::uint fi = 0; fi < _nunit; ++fi) {
        ny::fv_t &fc = _fcv[fi];
        std::sort (fc.rbegin (), fc.rend ()); // don't comment out
        if (fc.front () >> PSEUDO_TRIE_N[_d]) {
          ny::uint len = 0;
          byte_encoder encoder;
          for (ny::fv_it it = fc.begin (); it != fc.end (); ++it)
            len += encoder.encode (*it, s + len);
          pke_key.push_back (new FeatKey (s, &_fw_tmp[fi][0], len, _nl));
        } else {
          ny::uint pos = 0;
          for (ny::uint i = 0; i < fc.size (); ++i) {
            const ny::uint j = fc[i] - 1;
            switch (_d - i) {
              case 3: pos += j * (j - 1) * (j - 2) / 6;
              case 2: pos += j * (j - 1) / 2;
              case 1: pos += j;
            }
            if (i > 0) ++pos;
          }
          for (ny::uint li = 0; li < _nl; ++li)
            fw[pos * _nl + li] = _fw_tmp[fi][li];
        }
      }
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
          // uniq identical weights
          // (minimize the double array with a slight cahce-miss increase)
          w.assign (&wv[0], &wv[_nl]);
          std::pair <w2id_t::iterator, bool> itb
            = w2id.insert (w2id_t::value_type (w, uniq * _nl));
          if (itb.second) {
            for (ny::uint li = 0; li < _nl; ++li)
              fw[uniq * _nl + li] = wv[li]; // put weights for all labels
            ++uniq;
          }
          val.push_back (static_cast <int> (itb.first->second));
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
    }
    return true;
  }
  bool linear_model::load (const char * model) {
    TIMER (_model_t->startTimer ());
    const std::string model_bin (std::string (model) + ".bin"); // compiled model
    if (! _opt.force && newer (model_bin.c_str (), model)) {
      if (_opt.verbose > 0)
        std::fprintf (stderr, "loading compiled model parameters..");
      FILE* reader = std::fopen (model_bin.c_str (), "rb");
      std::fread (&_d,     sizeof (ny::uint), 1, reader);
      std::fread (&_nl,    sizeof (ny::uint), 1, reader);
      std::fread (&_nf,    sizeof (ny::uint), 1, reader);
      std::fread (&_nunit, sizeof (ny::uint), 1, reader);
      // label map
      for (ny::uint li = 0; li < _nl; ++li) {
        ny::uint len = 0;
        std::fread (&len, sizeof (ny::uint), 1, reader);
        char * p = new char[len + 1];
        std::fread (p, sizeof (char), len, reader);
        p[len] = '\0';
        _li2l.push_back (p);
        _l2li.insert (lmap::value_type (p, li));
        if (std::strcmp (p, "+1") == 0) _tli = li;
      }
      // feature map
      _fi2fn.reserve (_nf + 1);
      _fi2fn.push_back (0);
      for (ny::uint fn (0); std::fread (&fn, sizeof (ny::uint), 1, reader);) {
        if (fn >= _fn2fi.size ()) _fn2fi.resize (fn + 1, 0);
        _fn2fi[fn] = static_cast <ny::uint> (_fi2fn.size ());
        _fi2fn.push_back (fn);
      }
      std::fclose (reader);
      if (_opt.verbose > 0) std::fprintf (stderr, "done.\n");
    } else {
      FILE *reader = std::fopen (model, "r");
      if (_opt.verbose > 0)
        std::fprintf (stderr, "loading/compiling model parameters..");
      if (! reader)
        ny::print_err (HERE "no such file: %s\n", model);
      std::map <ny::fv_t, ny::uint> fc2fci; // conjunctive features
      char*  line = 0;
      size_t read = 0;
      while (ny::getLine (reader, line, read)) {
        if (*line != '\0') {
          char * p (line), * const p_end (p + read - 1);
          while (p != p_end && *p != '\t') ++p; *p = '\0';
          lmap::iterator it = _l2li.find (line);
          if (it == _l2li.end ()) { // identify label
            char * label = new char[read];
            std::strcpy (label, line);
            if (std::strcmp (label, "+1") == 0) _tli = _nl;
            it = _l2li.insert (std::make_pair (label, _nl++)).first;
            _li2l.push_back (label);
          }
          // read conjunctive feature
          const ny::uint li = it->second; // bug if all alpha null
          ny::fv_t fnv;
          _convertCfstr2Cf (p, fnv, true);
          _d = std::max (_d, static_cast <ny::uint> (fnv.size ())); // degree
          std::map <ny::fv_t, ny::uint>::iterator jt = fc2fci.lower_bound (fnv);
          if (jt == fc2fci.end () || jt->first > fnv) {
            jt = fc2fci.insert (jt, std::make_pair (fnv, _nunit));
            _fw_tmp.push_back (std::vector <ny::fl_t> ());
            _fcv.push_back (fnv);
            ++_nunit;
          }
          const ny::uint uid = jt->second;
          const ny::fl_t w = ny::strton <ny::fl_t> (++p);
          if (w == 0) continue;
          _fw_tmp[uid].resize (_nl, 0);
          _fw_tmp[uid][li] = w;
        }
      }
      std::fclose (reader);
      for (ny::uint i = 0; i < _nunit; ++i)
        _fw_tmp[i].resize (_nl, 0);
      if (_opt.verbose > 0) std::fprintf (stderr, "done.\n");
    
      // construct dense feature;
      _nf = static_cast <ny::uint> (_fncnt.size ());
      _packingFeatures (_fcv); // need afn

      // output map
      FILE* writer = std::fopen (model_bin.c_str (), "wb");
      std::fwrite (&_d,     sizeof (ny::uint), 1, writer);
      std::fwrite (&_nl,    sizeof (ny::uint), 1, writer);
      std::fwrite (&_nf,    sizeof (ny::uint), 1, writer);
      std::fwrite (&_nunit, sizeof (ny::uint), 1, writer);
      for (ny::uint li = 0; li < _nl; ++li) {
        const char * p = _li2l[li];
        ny::uint len = static_cast <ny::uint> (std::strlen (p));
        std::fwrite (&len, sizeof (ny::uint), 1,    writer);
        std::fwrite (p,    sizeof (char),     len, writer);
      }
      for (ny::uint fi = 1; fi <= _nf; ++fi)
        std::fwrite (&_fi2fn[fi], sizeof (ny::uint), 1, writer);
      std::fclose (writer);
      _fncnt.clear ();
    }
    _setFtrie ();
    // fstrie construction
    if (_opt.algo == FST) _setFStrie ();
    TIMER (_model_t->stopTimer ());
    if (_opt.verbose > 0) printParam ();
    return true;
  }
}
