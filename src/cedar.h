// CEDAR -- C++-implemented Efficient Double ARray library
//  $Id: cedar.h 808 2012-04-16 18:46:44Z ynaga $
// Copyright (c) 2009-2012 Naoki Yoshinaga <ynaga@tkl.iis.u-tokyo.ac.jp>
#ifndef CEDAR_H
#define CEDAR_H

#include <stdint.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <algorithm>

namespace cedar {
  // type alias
  typedef unsigned int   uint;
  typedef unsigned char  uchar;
  // dynamic double array
  template <typename value_type, const int NAN1, const int NAN2,
            const size_t MAX_TRIAL = 8, const size_t SIZE = 0>
  class da {
  public:
    typedef value_type result_type;
    struct result_pair_t {
      value_type value;
      size_t     length;
    };
    struct node {
      union { int base; value_type value; }; // neg to store prev empty index
      int check;                             // neg to store next empty index
      node () : base (0), check (0) {};
      node (const int base_, const int check_) : base (base_), check (check_) {};
    };
    struct ninfo {   // x1.5 update speed; +.25 % memory (8n -> 10n)
      uchar child;   // first child
      uchar sibling; // right sibling (= 0 if not exist)
      ninfo () : child (0), sibling (0) {};
      ninfo (const uchar child_, const uchar sibling_) : child (child_), sibling (sibling_) {};
    };
    struct block {  // a block w/ 256 element
      uint16_t num;   // # empty element; 0 - 256
      uint8_t  trial; // # trial
      uint8_t  ehead; // offset to first empty item
      int      prev;  // prev block; 3bit
      int      next;  // next block; 3bit
      block () : num (256), trial (0), ehead (0), prev (0), next (0) {};
    };
    const int NO_VALUE;
    const int NO_PATH;
    size_t pid[SIZE + 1]; // for traverse
    da () : NO_VALUE (NAN1), NO_PATH (NAN2), _array (0), _ninfo (0), _block (0), _bheadF (0), _bheadC (0), _bheadO (0), _size (0), _capacity (0), _no_delete (false)
    { _initialize (); }
    ~da () { clear (false); }
    size_t capacity () const { return static_cast <size_t> (_capacity); }
    size_t size     () const { return static_cast <size_t> (_size); }
    size_t nonzero_size () const {
      uint n = _size - _block[0].num;
      if (uint bi = _bheadC)
        do n -= _block[bi].num; while ((bi = _block[bi].next) != _bheadC);
      if (uint bi = _bheadO)
        do n -= _block[bi].num; while ((bi = _block[bi].next) != _bheadO);
      return n;
    }
    // interfance
    template <class T>
    T exactMatchSearch (const char* key, size_t len = 0, size_t from = 0) const {
      union { int i; value_type x; } b;
      size_t pos = 0;
      T result;
      if ((b.i = _find (key, from, pos, len)) == NO_PATH) b.i = NO_VALUE;
      _set_result (&result, b.x, len);
      return result;
    }
    template <typename U>
    size_t commonPrefixSearch (const char* key, U* result, size_t result_len, size_t len = 0, size_t from = 0) const {
      size_t num (0), pos (0);
      while (pos < len) {
        union { int i; value_type x; } b;
        b.i = _find (key, from, pos, pos + 1);
        if (b.i == NO_VALUE) continue;
        if (b.i == NO_PATH)  return num;
        if (num < result_len) { set_result (&result[num], b.x, len); ++result; }
        ++num;
      }
      return num;
    }
    value_type traverse (const char* key, size_t &from, size_t &pos, size_t len = 0) const {
      union { int i; value_type x; } b;
      b.i = _find (key, from, pos, len);
      return b.x;
    }
    value_type& update (const char *key, size_t len = 0, value_type val = value_type (0)) {
      size_t from = 0;
      return update (key, from, len, val);
    }
    value_type& update (const char *key, size_t &from, size_t len, value_type val) {
      if (! len) len = std::strlen (key);
      // follow link
      for (size_t pos = 0;; ++pos) {
        const int   base = _array[from].base;
        const uchar c    = static_cast <uchar> (pos == len ? 0 : key[pos]);
        int         to   = base ^ c;
        if (base < 0 || _array[to].check < 0) {
          to = _pop_enode (base, c, static_cast <int> (from));
          _push_sibling (from, base, to, c);
        } else if (_array[to].check != static_cast <int> (from)) {
          to = _resolve (from, base, c);
        }
        if (pos == len) return _array[to].value += val;
        from = static_cast <size_t> (to);
      }
    }
    // recover prefix from the given node
    void reverse (size_t to, std::string &val) {
      for (; to; to = _array[to].check)
        val += _array[_array[to].check].base ^ to;
      std::reverse (val.begin (), val.end ());
    }
    // # keys inserted
    size_t nkeys () {
      size_t m = 0;
      for (int i = 1; i < static_cast <int> (_size); ++i) // check bugs
        if (_array[i].check >= 0 && _array[_array[i].check].base == i) ++m;
      return m;
    }
    int save (const char* fn, const char* mode = "wb") const {
      FILE * writer = std::fopen (fn, mode);
      if (! writer) return -1;
      std::fwrite (_array, sizeof (node), static_cast <size_t> (_size), writer);
      std::fclose (writer);
      return 0;
    }
    int open (const char * fn, const char * mode = "rb",
              const size_t offset = 0, size_t size = 0) {
      FILE * reader = std::fopen (fn, mode);
      if (! reader) return -1;
      // get size
      if (! size) {
        if (std::fseek (reader, 0, SEEK_END) != 0) return -1;
        size = static_cast <size_t> (std::ftell (reader)) / sizeof (node);
        if (std::fseek (reader, 0, SEEK_SET) != 0) return -1;
        if (size <= offset) return -1;
        size -= offset;
      }
      // set array
      if (_array) delete [] _array;
      _array = new node[size];
      if (std::fseek (reader, static_cast <long> (offset), SEEK_SET) != 0)
        return -1;
      if (size != std::fread (reinterpret_cast <node*> (_array),
                               sizeof (node), size, reader)) return -1;
      std::fclose (reader);
      _size = static_cast <int> (size);
      return 0;
    }
    int build (size_t num, const char** key, const size_t* len = 0, const value_type* val = 0) {
      for (size_t i = 0; i < num; ++i)
        update (key[i], len ? len[i] : 0, val ? val[i] : 0);
      return 0;
    }
    size_t unit_size () const { return sizeof (node); }
    void set_array (void *p, size_t size = 0) { // TODO: immutable
      clear (false);
      _no_delete = true;
      _array = reinterpret_cast <node*> (p);
      _size  = static_cast <int> (size);
    }
    void clear (bool reuse = true) {
      if (_array && ! _no_delete) delete [] _array; _array = 0;
      if (_ninfo && ! _no_delete) delete [] _ninfo; _ninfo = 0;
      if (_block && ! _no_delete) delete [] _block; _block = 0;
      _bheadF = _bheadC = _bheadO = _size = _capacity = 0; // *
      if (reuse)_initialize ();
      _no_delete = false;
    }
  private:
    // currently disabled; implement these if you need
    da (const da&);
    da& operator= (const da&);
    node*  _array;
    ninfo* _ninfo;
    block* _block;
    int    _bheadF; // first block of Full;   0
    int    _bheadC; // first block of Closed; 0 if no Closed
    int    _bheadO; // first block of Open;   0 if no Open
    int    _size;
    int    _capacity;
    bool   _no_delete;

    template <typename U>
    void _realloc_array (U* &array, const int size_n, const int size_p) {
      U* tmp = new U[size_n];
      if (array) std::copy (array, array + size_p, tmp);
      std::swap (tmp, array);
      if (tmp) delete [] tmp;
    }
    void _initialize () { // initilize the first special block
      _realloc_array (_array, 256, _size);
      _realloc_array (_ninfo, 256, _size);
      _realloc_array (_block, 1,   _size);
      for (int i = 1; i < 256; ++i)
        _array[i] = node (i == 1 ? -255 : -(i-1), i == 255 ? -1 : -(i+1));
      _size = _capacity = 256;
      std::memset (&pid[0], 0, sizeof (size_t) * (SIZE + 1));
    }
    void _pop_block (const bool last, int &head_in, const int bi) {
      const block &b = _block[bi];
      if (last) { // last one poped; Closed or Open
        head_in = 0;
      } else {
        _block[b.prev].next = b.next;
        _block[b.next].prev = b.prev;
        if (bi == head_in) head_in = b.next; // head_in's next
      }
    }
    void _push_block (const bool empty, int &head_out, const int bi) {
      block &b = _block[bi];
      if (empty) { // the destination is empty
        head_out = b.prev = b.next = bi;
      } else {
        int &tail_out = _block[head_out].prev;
        b.prev = tail_out;
        b.next = head_out;
        tail_out = _block[tail_out].next = bi;
      }
    }
    int _add_block () {
      const int n = _size + 256;
      if (n > _capacity) { // allocate memory if needed
        while (n > _capacity) _capacity *= 2;
        _realloc_array (_array, _capacity, _size);           // copy _array
        _realloc_array (_ninfo, _capacity, _size);           // copy _ninfo
        _realloc_array (_block, _capacity >> 8, _size >> 8); // copy _block
      }
      for (int i = _size; i < n; ++i)
        _array[i] = node ((i & 0xff) == 0    ? -(i+255) : -(i-1),
                          (i & 0xff) == 0xff ? -(i-255) : -(i+1));
      _push_block (!_bheadO, _bheadO, _size >> 8); // append to block Open
      _size += 256;
      return (_size >> 8) - 1;
    }
    // find key from double _array
    int _find (const char* key, size_t &from, size_t &pos, size_t len) const {
      if (! len) len = std::strlen (key);
      // follow edge to get value of key
      for (; pos < len; ++pos) {
        const int to = _array[from].base ^ static_cast <uchar> (key[pos]);
        if (_array[to].check != static_cast <int> (from)) return NO_PATH;
        from = static_cast <size_t> (to);
      }
      // get value
      const node &n = _array[_array[from].base ^ 0];
      return (n.check == static_cast <int> (from) ? n.base : NO_VALUE);
    };
    void _set_result (result_type *x, value_type r, size_t) const
    { *x = r; }
    void _set_result (result_pair_t *x, value_type r, size_t l) const
    { x->value = r; x->length = l; }
    // pop empty node from block; never transfer the special block (bi = 0)
    int _pop_enode (int base, const uchar label, const int from) {
      int bi = 0;
      if (base < 0) { // choose a block to put a single element
        bi = (_bheadC ? _bheadC : (_bheadO ? _bheadO : _add_block ()));
        base = _array[from].base = ((bi << 8) ^ _block[bi].ehead) ^ label;
      } else {
        bi = base >> 8;
      }
      const int e = base ^ label;
      node  &n = _array[e];
      block &b = _block[bi];
      if (--b.num == 0) {
        if (bi) _transfer_block (bi, _bheadC, _bheadF); // Closed to Full
      } else { // release empty node from empty ring
        _array[-n.base].check = n.check;
        _array[-n.check].base = n.base;
        if ((e & 0xff) == b.ehead) b.ehead = static_cast <uchar> (-n.check & 0xff); // set ehead
        if (bi && b.num == 1 && b.trial != MAX_TRIAL) // Open to Closed
          _transfer_block (bi, _bheadO, _bheadC);
      }
      n = node (label ? NO_PATH : 0, static_cast <int> (from)); // initialize the released node
      return e;
    }
    // push empty node into empty ring
    void _push_enode (const int e) {
      const int bi = e >> 8;
      block &b = _block[bi];
      if (++b.num == 1) { // Full to Closed
        b.ehead = static_cast <uint8_t> (e & 0xff);
        _array[e] = node (-e, -e);
        _transfer_block (bi, _bheadF, _bheadC); // Full to Closed
      } else {
        const int next = (bi << 8) ^ b.ehead;
        const int prev = -_array[next].base;
        _array[prev].check = _array[next].base = -e;
        _array[e] = node (-prev, -next);
        if (b.num == 2 || b.trial == MAX_TRIAL)  // Closed to Open
          _transfer_block (bi, _bheadC, _bheadO);
        b.trial = 0;
      }
      _ninfo[e] = ninfo (0, 0); // reset ninfo; no child, no brother
    }
    // transfer block from one start w/ head_in to one start w/ head_out
    void _transfer_block (const int bi, int &head_in, int &head_out) {
      _pop_block  (bi == _block[bi].next, head_in, bi);
      _push_block (!head_out && _block[bi].num, head_out, bi);
    }
    // push label to from's child
    void _push_sibling (const size_t from, const int base, const int to, const uchar label) {
      uchar &child = _ninfo[from].child;
      if (base < 0 || child) { // first child or non-zero child
        _ninfo[to] = ninfo (0, child);
        child = label;
      } else {
        _ninfo[to] = ninfo (0, _ninfo[base ^ child].sibling);
        _ninfo[base ^ child].sibling = label;
      }
    }
    // enumerate (more than one) children
    int _get_child (uchar *first, const int base, uchar c, uint l = 256) {
      uchar *last = first;
      if (!c)   { *last++ = c; c = _ninfo[base ^ c].sibling; } // 0 (terminal)
      if (l < 256) *last++ = static_cast <uchar> (l);
      while (c) { *last++ = c; c = _ninfo[base ^ c].sibling; }
      return static_cast <int> (last - first);
    }
    // explore new block to settle down
    int _find_place (const int nc, uchar* const &firstborn, uchar* const &lastborn) {
      if (_bheadC && nc == 1)
        return (_bheadC << 8) ^ _block[_bheadC].ehead ^ *firstborn;
      if (int bi = _bheadO) {
        if (nc == 1) return (bi << 8) ^ _block[bi].ehead ^ *firstborn;
        do { // set candidate block
          block &b = _block[bi];
          if (b.num >= nc) { // explore configuration
            int e = (bi << 8) ^ b.ehead;
            do {
              const int base = e ^ *firstborn;
              for (const uchar *p = firstborn + 1; _array[base ^ *p].check < 0; )
                if (++p == lastborn) return base; // no conflict
            } while (((e = -_array[e].check) & 0xff) != b.ehead);
          }
          const int bnext = b.next;
          if (++b.trial == MAX_TRIAL) _transfer_block (bi, _bheadO, _bheadC);
          bi = bnext;
        } while (_bheadO && bi != _bheadO);
      }
      return _size ^ *firstborn;
    }
    // resolve conflict on base_n ^ label_n = base_p ^ label_p
    int _resolve (size_t &from_n, const int base_n, const uchar label_n) {
      // examine sibling of conflicted nodes
      const int to_pn  = base_n ^ label_n;
      const int from_p = _array[to_pn].check;
      const int base_p = _array[from_p].base;
      // siblings of current occupier
      uchar child_p[256]; // <= 255
      const int nc_p = _get_child (child_p, base_p, _ninfo[from_p].child);
      // siblings of incoming occupier
      uchar child_n[256]; // <= 255
      // const uint nc_n = _get_child (child_n, base_n, _ninfo[from_n].child, label_n);
      const int nc_n = (nc_p == 1 ? 256  // always replace single state
                        : _get_child (child_n, base_n, _ninfo[from_n].child, label_n));
      // search new address
      const bool flag = (nc_n <= nc_p); // replace new
      uchar* const firstborn = flag ? &child_n[0] : &child_p[0];
      uchar* const lastborn  = firstborn + (flag ? nc_n : nc_p);
      const int base = _find_place (flag ? nc_n : nc_p, firstborn, lastborn);
      const int bi   = base >> 8;
      if (base >= _size) _add_block ();
      if (_block[bi].num > 1 && _block[bi].trial == MAX_TRIAL) // bug fix
        _transfer_block (bi, _bheadC, _bheadO); // Closed to Open
      _block[bi].trial = 0; // reset
      // replace & modify empty list
      const int from  = flag ? static_cast <int> (from_n) : from_p;
      const int base_ = flag ? base_n : base_p;
      if (flag && *firstborn == label_n) _ninfo[from].child = label_n;
      _array[from].base = base;
      for (const uchar* p = firstborn; p != lastborn; ++p) { // to_ => to
        const int to  = _pop_enode (base, *p, from);
        const int to_ = base_ ^ *p;
        node  &n = _array[to];
        ninfo &i = _ninfo[to];
        i.sibling = (p+1 == lastborn ? static_cast <uchar> (0) : *(p+1));
        if (flag && to_ == to_pn) continue; // skip newcomer
        for (uint j = 0; pid[j] != 0; ++j) // keep updated for traverse ()
          if (pid[j] == static_cast <size_t> (to_))
            { pid[j] = static_cast <size_t> (to); break; }
        if ((n.base = _array[to_].base) > 0 && *p) { // copy base; bugfix
          uchar c = i.child = _ninfo[to_].child;
          do _array[n.base ^ c].check = to; // adjust grand son's check
          while ((c = _ninfo[n.base ^ c].sibling));
        }
        if (! flag && to_ == static_cast <int> (from_n))
          from_n = static_cast <size_t> (to); // bug fix
        if (! flag && to_ == to_pn) { // the address is immediately used
          _array[to_] = node (label_n ? NO_PATH : 0, static_cast <int> (from_n));
          _push_sibling (from_n, base_n, to_pn, label_n);
        } else {
          _push_enode (to_);
        }
      }
      return flag ? base ^ label_n : to_pn;
    }
  };
}
#endif
