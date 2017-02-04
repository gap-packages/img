/****************************************************************************
 *
 * fr_dll.h                                                 Laurent Bartholdi
 *
 * Copyright (C) 2010-2013, Laurent Bartholdi
 *
 ****************************************************************************
 *
 * header / type declarations for FR DLL add-on
 *
 ****************************************************************************/

#undef VERY_LONG_DOUBLES

/* #define sigjmp_buf int */
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <setjmp.h>
#include "src/compiled.h"
#include "src/macfloat.h"
#include "poly.h"

#ifdef MALLOC_HACK
#include <malloc.h>
#endif

void InitP1Kernel(void);
void InitP1Library(void);

/****************************************************************************
 * stolen from src/float.c
 ****************************************************************************/
#define VAL_FLOAT(obj) (*(Double *)ADDR_OBJ(obj))
#define SIZE_FLOAT   sizeof(Double)
#ifndef T_FLOAT
#define T_FLOAT T_MACFLOAT
#endif
static inline Obj NEW_FLOAT (Double val)
{
  Obj f = NewBag(T_FLOAT, SIZE_FLOAT);
  *(Double *)ADDR_OBJ(f) = val;
  return f;
}

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

/* fr_dll.h . . . . . . . . . . . . . . . . . . . . . . . . . . . ends here */
