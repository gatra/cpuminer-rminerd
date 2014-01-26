
if WANT_JANSSON
JANSSON_INCLUDES= -I$(top_srcdir)/compat/jansson
else
JANSSON_INCLUDES=
endif

EXTRA_DIST	= example-cfg.json nomacro.pl

SUBDIRS		= compat

INCLUDES	= $(PTHREAD_FLAGS) -fno-strict-aliasing $(JANSSON_INCLUDES)

bin_PROGRAMS	= rminerd

rminerd_SOURCES	= elist.h miner.h compat.h \
		  cpu-miner.c util.c \
		  scanriecoin.c \
		  sha2.c scrypt.c
if ARCH_x86
rminerd_SOURCES += sha2-x86.S scrypt-x86.S
endif
if ARCH_x86_64
rminerd_SOURCES += sha2-x64.S scrypt-x64.S
endif
if ARCH_ARM
rminerd_SOURCES += sha2-arm.S scrypt-arm.S
endif
rminerd_LDFLAGS	= $(PTHREAD_FLAGS)
rminerd_LDADD	= @LIBCURL@ @JANSSON_LIBS@ @PTHREAD_LIBS@ @WS2_LIBS@ -lgmp -lm
rminerd_CPPFLAGS = @LIBCURL_CPPFLAGS@

