# AM_CXXFLAGS = -Wall -Wextra -Wformat=2 -Wcast-qual -Wcast-align -Wwrite-strings -Wconversion -Wpointer-arith -pedantic
AM_CXXFLAGS = -Wall
bin_PROGRAMS = jdepp

jdepp_SOURCES = typedef.cc typedef.h \
                classify.cc  classify.h \
                pecco.h \
                pdep.cc pdep.h \
                jdepp.cc

EXTRA_jdepp_SOURCES = timer.cc timer.h \
                      kernel.cc kernel.h \
		      linear.cc linear.h \
		      pa.h \
                      cedar.h \
		      lbfgs.cpp lbfgs.h mathvec.h sgd.cpp owlqn.cpp \
		      maxent.cpp maxent.h

jdepp_CPPFLAGS = @JDEPP_CPPFLAGS@
jdepp_LDFLAGS = @JDEPP_LDFLAGS@
jdepp_DEPENDENCIES = @JDEPP_OBJ@
jdepp_LDADD = @JDEPP_OBJ@
