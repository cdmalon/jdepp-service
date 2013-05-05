// pecco -- please enjoy classification with conjunctive features
//  $Id: typedef.cc 827 2012-05-10 14:12:37Z ynaga $
// Copyright (c) 2008-2012 Naoki Yoshinaga <ynaga@tkl.iis.u-tokyo.ac.jp>
#include "typedef.h"

namespace ny {
  // (a slightly) secure conversion from string to numerical
  template <> float  strton <float>  (const char * s, char ** error)
  { return std::strtof (s, error); }
  template <> double strton <double> (const char * s, char ** error)
  { return std::strtod (s, error); }
  template <> uint   strton <uint>   (const char * s, char ** error) {
    uint64_t ret = std::strtoul (s, error, 10);
    if (ret > UINT_MAX)
      print_err (HERE "overflow: %s\n", s);
    return static_cast <uint> (ret);
  }
  template <> int    strton <int>    (const char * s, char ** error) {
    int64_t ret = std::strtol (s, error, 10);
    if (ret > INT_MAX)      print_err (HERE "overflow: %s\n", s);
    else if (ret < INT_MIN) print_err (HERE "underflow: %s\n", s);
    return static_cast <int> (ret);
  }
  template <> long   strton <long>   (const char * s, char ** error)
  { return std::strtol (s, error, 10); }
  template <typename T> T strton (const char * s) {
    char* err;
    T n = strton <T> (s, &err);
    if (*err != '\0') print_err (HERE "invalid conversion: %s\n", s);
    return n;
  }
  // explict specialization
  template float  strton (const char *);
  template double strton (const char *);
  template uint   strton (const char *);
  template int    strton (const char *);
}
