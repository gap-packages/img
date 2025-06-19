/****************************************************************************
 *
 * img.h                                                    Laurent Bartholdi
 *
 * Copyright (C) 2010-2013, Laurent Bartholdi
 *
 ****************************************************************************
 *
 * header / type declarations for IMG DLL add-on
 *
 ****************************************************************************/

#undef VERY_LONG_DOUBLES

#define _GNU_SOURCE
#include "gap_all.h"

#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "poly.h"

void InitP1Kernel(void);
void InitP1Library(void);

static inline Obj ALLOC_PLIST (UInt len)
{
  Obj f = NEW_PLIST(T_PLIST, len);
  SET_LEN_PLIST(f, len);
  return f;
}

static void set_elm_plist(Obj list, UInt pos, Obj obj) /* safe to nest */
{
  SET_ELM_PLIST(list,pos,obj);
  CHANGED_BAG(list);
}

/* img.h . . . . . . . . . . . . . . . . . . . . . . . . . . . ends here */
