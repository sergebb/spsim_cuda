/* A Bison parser, made by GNU Bison 2.1.  */

/* Skeleton parser for Yacc-like parsing with Bison,
   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005 Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor,
   Boston, MA 02110-1301, USA.  */

/* As a special exception, when this file is copied by Bison into a
   Bison output file, you may use that output file without restriction.
   This special exception was added by the Free Software Foundation
   in version 1.24 of Bison.  */

/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     TOK_BOOLEAN = 258,
     TOK_INTEGER = 259,
     TOK_HEX = 260,
     TOK_FLOAT = 261,
     TOK_STRING = 262,
     TOK_NAME = 263,
     TOK_EQUALS = 264,
     TOK_NEWLINE = 265,
     TOK_ARRAY_START = 266,
     TOK_ARRAY_END = 267,
     TOK_LIST_START = 268,
     TOK_LIST_END = 269,
     TOK_COMMA = 270,
     TOK_GROUP_START = 271,
     TOK_GROUP_END = 272,
     TOK_END = 273,
     TOK_GARBAGE = 274
   };
#endif
/* Tokens.  */
#define TOK_BOOLEAN 258
#define TOK_INTEGER 259
#define TOK_HEX 260
#define TOK_FLOAT 261
#define TOK_STRING 262
#define TOK_NAME 263
#define TOK_EQUALS 264
#define TOK_NEWLINE 265
#define TOK_ARRAY_START 266
#define TOK_ARRAY_END 267
#define TOK_LIST_START 268
#define TOK_LIST_END 269
#define TOK_COMMA 270
#define TOK_GROUP_START 271
#define TOK_GROUP_END 272
#define TOK_END 273
#define TOK_GARBAGE 274




#if ! defined (YYSTYPE) && ! defined (YYSTYPE_IS_DECLARED)
#line 59 "grammar.y"
typedef union YYSTYPE {
  long ival;
  double fval;
  char *sval;
} YYSTYPE;
/* Line 1447 of yacc.c.  */
#line 82 "grammar.h"
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif





