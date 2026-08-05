//
// (C) Copyright 2011-2018 Sergey A. Babkin.
// This file is a part of Triceps.
// See the file COPYRIGHT for the copyright notice and license information
//
//
// Assorted functions working on strings.

#ifndef __Triceps_StringTools_h__
#define __Triceps_StringTools_h__

#include <stdarg.h>
#include <string>

// like sprintf() but returns a C++ string
std::string strprintf(const char *fmt, ...)
	__attribute__((format(printf, 1, 2)));

// like vsprintf() but returns a C++ string
std::string vstrprintf(const char *fmt, va_list ap)
	__attribute__((format(printf, 1, 0)));

#endif // __Triceps_StringTools_h__
