// J.DepP -- Japanese Dependency Parsers
//  $Id: jdepp.cc 788 2012-04-06 14:31:21Z ynaga $
// Copyright (c) 2008-2012 Naoki Yoshinaga <ynaga@tkl.iis.u-tokyo.ac.jp>
#include "pdep.h"

int main (int argc, char* argv[]) {

  pdep::option opt (argc, argv);
  pdep::parser parser (opt);

  parser.run ();
  return 0;
}
