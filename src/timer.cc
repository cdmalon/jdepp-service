// pecco -- please enjoy classification with conjunctive features
//  $Id: timer.cc 827 2012-05-10 14:12:37Z ynaga $
// Copyright (c) 2008-2012 Naoki Yoshinaga <ynaga@tkl.iis.u-tokyo.ac.jp>
#include "timer.h"

namespace ny {
  
  const long double Timer::clock = Timer::getCPUClock ();

  long double Timer::getCPUClock () {
    register uint64_t start, end;
    timeval _s_t, _e_t;
    register int min_interval = 5000;
    uint64_t interval = 0;
    gettimeofday (&_s_t, NULL);
    start = rdtsc ();
    gettimeofday (&_e_t, NULL);
    while ((_e_t.tv_sec - _s_t.tv_sec) * 1000000 +
           _e_t.tv_usec - _s_t.tv_usec < min_interval)
      gettimeofday (&_e_t, NULL);
    end = rdtsc ();
    interval = static_cast <uint64_t> ((_e_t.tv_sec - _s_t.tv_sec) * 1000000
                                       + _e_t.tv_usec - _s_t.tv_usec);
    return (end - start) / interval;
  }
  void Timer::printElapsed () const {
    if (_trial > 0) {
      if (_trial == 1)
        std::fprintf (stderr, "%-10s: %.4Lf ms.\n",
                      _label.c_str (), ( _elapsed / clock) / 1000);
      else
        std::fprintf (stderr, "%-10s: %.4Lf ms./%s (%.8Lf/%d)\n",
                      _label.c_str (), (_elapsed / clock) / 1000 / _trial,
                      _unit.c_str (),  (_elapsed / clock) / 1000, _trial);
    }
  }
}
