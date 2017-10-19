/****************************************************************************
 *
 * MODULE:	Compatibility with GRASS GeoPAT
 * AUTHOR(S):	Pawel Netzel
 * COPYRIGHT:	(C) Space Informatics Lab, University of Cincinnati
 *
 *		This program is free software under the GNU General Public
 *		License (>=v2). Read the file COPYING that comes with GRASS
 *		for details.
 *
 *****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <time.h>
#include <math.h>

#include <sml.h>
#include <ezgdal.h>

#include "geopat_compatibility.h"

#define WARN 0
#define ERR  1
#define MSG  2

static void print_msg(int type, const char *template, va_list ap) {
  char buffer[2000];

  if(type==ERR)
    sprintf(buffer, "ERROR: %s\n", template);
  else if(type==WARN)
    sprintf(buffer, "WARNING: %s\n", template);
  else
    sprintf(buffer, "%s\n", template);

  vprintf(buffer, ap);
}

void G_message(const char *msg, ...) {
  va_list ap;

  va_start(ap, msg);
  print_msg(MSG, msg, ap);
  va_end(ap);

}

void G_fatal_error(const char *msg, ...) {
  static int busy;
  va_list ap;

  if (busy)
    exit(1);
  busy = 1;

  va_start(ap, msg);
  print_msg(ERR, msg, ap);
  va_end(ap);

  exit(1);
}

void G_warning(const char *msg, ...) {
  va_list ap;

  va_start(ap, msg);
  print_msg(WARN, msg, ap);
  va_end(ap);
}

