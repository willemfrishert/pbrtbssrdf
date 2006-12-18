/* A Bison parser, made by GNU Bison 1.875b.  */

/* Skeleton parser for Yacc-like parsing with Bison,
   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003 Free Software Foundation, Inc.

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
   Foundation, Inc., 59 Temple Place - Suite 330,
   Boston, MA 02111-1307, USA.  */

/* As a special exception, when this file is copied by Bison into a
   Bison output file, you may use that output file without restriction.
   This special exception was added by the Free Software Foundation
   in version 1.24 of Bison.  */

/* Written by Richard Stallman by simplifying the original so called
   ``semantic'' parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Using locations.  */
#define YYLSP_NEEDED 0



/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     STRING = 258,
     ID = 259,
     NUM = 260,
     LBRACK = 261,
     RBRACK = 262,
     ACCELERATOR = 263,
     AREALIGHTSOURCE = 264,
     ATTRIBUTEBEGIN = 265,
     ATTRIBUTEEND = 266,
     CAMERA = 267,
     CONCATTRANSFORM = 268,
     COORDINATESYSTEM = 269,
     COORDSYSTRANSFORM = 270,
     FILM = 271,
     IDENTITY = 272,
     LIGHTSOURCE = 273,
     LOOKAT = 274,
     MATERIAL = 275,
     OBJECTBEGIN = 276,
     OBJECTEND = 277,
     OBJECTINSTANCE = 278,
     PIXELFILTER = 279,
     REVERSEORIENTATION = 280,
     ROTATE = 281,
     SAMPLER = 282,
     SCALE = 283,
     SEARCHPATH = 284,
     SHAPE = 285,
     SURFACEINTEGRATOR = 286,
     TEXTURE = 287,
     TRANSFORMBEGIN = 288,
     TRANSFORMEND = 289,
     TRANSFORM = 290,
     TRANSLATE = 291,
     VOLUME = 292,
     VOLUMEINTEGRATOR = 293,
     WORLDBEGIN = 294,
     WORLDEND = 295,
     HIGH_PRECEDENCE = 296
   };
#endif
#define STRING 258
#define ID 259
#define NUM 260
#define LBRACK 261
#define RBRACK 262
#define ACCELERATOR 263
#define AREALIGHTSOURCE 264
#define ATTRIBUTEBEGIN 265
#define ATTRIBUTEEND 266
#define CAMERA 267
#define CONCATTRANSFORM 268
#define COORDINATESYSTEM 269
#define COORDSYSTRANSFORM 270
#define FILM 271
#define IDENTITY 272
#define LIGHTSOURCE 273
#define LOOKAT 274
#define MATERIAL 275
#define OBJECTBEGIN 276
#define OBJECTEND 277
#define OBJECTINSTANCE 278
#define PIXELFILTER 279
#define REVERSEORIENTATION 280
#define ROTATE 281
#define SAMPLER 282
#define SCALE 283
#define SEARCHPATH 284
#define SHAPE 285
#define SURFACEINTEGRATOR 286
#define TEXTURE 287
#define TRANSFORMBEGIN 288
#define TRANSFORMEND 289
#define TRANSFORM 290
#define TRANSLATE 291
#define VOLUME 292
#define VOLUMEINTEGRATOR 293
#define WORLDBEGIN 294
#define WORLDEND 295
#define HIGH_PRECEDENCE 296




/* Copy the first part of user declarations.  */
#line 11 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"

#include "api.h"
#include "pbrt.h"
#include "paramset.h"
#include <stdarg.h>

extern int yylex( void );
int line_num = 0;
string current_file;

#define YYMAXDEPTH 100000000

void yyerror( char *str ) {
	Severe( "Parsing error: %s", str);
}

void ParseError( const char *format, ... ) PRINTF_FUNC;

void ParseError( const char *format, ... ) {
	char error[4096];
	va_list args;
	va_start( args, format );
	vsprintf( error, format, args );
	yyerror(error);
	va_end( args );
}

int cur_paramlist_allocated = 0;
int cur_paramlist_size = 0;
const char **cur_paramlist_tokens = NULL;
void **cur_paramlist_args = NULL;
int *cur_paramlist_sizes = NULL;
bool *cur_paramlist_texture_helper = NULL;

#define CPS cur_paramlist_size
#define CPT cur_paramlist_tokens
#define CPA cur_paramlist_args
#define CPTH cur_paramlist_texture_helper
#define CPSZ cur_paramlist_sizes

typedef struct ParamArray {
	int element_size;
	int allocated;
	int nelems;
	void *array;
} ParamArray;

ParamArray *cur_array = NULL;
bool array_is_single_string = false;

#define NA(r) ((float *) r->array)
#define SA(r) ((const char **) r->array)

void AddArrayElement( void *elem ) {
	if (cur_array->nelems >= cur_array->allocated) {
		cur_array->allocated = 2*cur_array->allocated + 1;
		cur_array->array = realloc( cur_array->array,
			cur_array->allocated*cur_array->element_size );
	}
	char *next = ((char *)cur_array->array) + cur_array->nelems *
		cur_array->element_size;
	memcpy( next, elem, cur_array->element_size );
	cur_array->nelems++;
}

ParamArray *ArrayDup( ParamArray *ra )
{
	ParamArray *ret = new ParamArray;
	ret->element_size = ra->element_size;
	ret->allocated = ra->allocated;
	ret->nelems = ra->nelems;
	ret->array = malloc(ra->nelems * ra->element_size);
	memcpy( ret->array, ra->array, ra->nelems * ra->element_size );
	return ret;
}

void ArrayFree( ParamArray *ra )
{
	free(ra->array);
	delete ra;
}

void FreeArgs()
{
	for (int i = 0; i < cur_paramlist_size; ++i)
		delete[] ((char *)cur_paramlist_args[i]);
}

static bool VerifyArrayLength( ParamArray *arr, int required,
	const char *command ) {
	if (arr->nelems != required) {
		ParseError( "%s requires a(n) %d element array!", command, required);
		return false;
	}
	return true;
}
enum { PARAM_TYPE_INT, PARAM_TYPE_BOOL, PARAM_TYPE_FLOAT, PARAM_TYPE_POINT,
	PARAM_TYPE_VECTOR, PARAM_TYPE_NORMAL, PARAM_TYPE_COLOR,
	PARAM_TYPE_STRING, PARAM_TYPE_TEXTURE };
static void InitParamSet(ParamSet &ps, int count, const char **tokens,
	void **args, int *sizes, bool *texture_helper);
static bool lookupType(const char *token, int *type, string &name);
#define YYPRINT(file, type, value)  \
{ \
	if ((type) == ID || (type) == STRING) \
		fprintf ((file), " %s", (value).string); \
	else if ((type) == NUM) \
		fprintf ((file), " %f", (value).num); \
}


/* Enabling traces.  */
#ifndef YYDEBUG
# define YYDEBUG 1
#endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

#if ! defined (YYSTYPE) && ! defined (YYSTYPE_IS_DECLARED)
#line 122 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"
typedef union YYSTYPE {
char string[1024];
float num;
ParamArray *ribarray;
} YYSTYPE;
/* Line 191 of yacc.c.  */
#line 275 "c:\\Downloads\\pbrt-src-1.02\\core\\/pbrtparse.cpp"
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif



/* Copy the second part of user declarations.  */


/* Line 214 of yacc.c.  */
#line 287 "c:\\Downloads\\pbrt-src-1.02\\core\\/pbrtparse.cpp"

#if ! defined (yyoverflow) || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# if YYSTACK_USE_ALLOCA
#  define YYSTACK_ALLOC alloca
# else
#  ifndef YYSTACK_USE_ALLOCA
#   if defined (alloca) || defined (_ALLOCA_H)
#    define YYSTACK_ALLOC alloca
#   else
#    ifdef __GNUC__
#     define YYSTACK_ALLOC __builtin_alloca
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's `empty if-body' warning. */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (0)
# else
#  if defined (__STDC__) || defined (__cplusplus)
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   define YYSIZE_T size_t
#  endif
#  define YYSTACK_ALLOC malloc
#  define YYSTACK_FREE free
# endif
#endif /* ! defined (yyoverflow) || YYERROR_VERBOSE */


#if (! defined (yyoverflow) \
     && (! defined (__cplusplus) \
	 || (YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  short yyss;
  YYSTYPE yyvs;
  };

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (short) + sizeof (YYSTYPE))				\
      + YYSTACK_GAP_MAXIMUM)

/* Copy COUNT objects from FROM to TO.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if 1 < __GNUC__
#   define YYCOPY(To, From, Count) \
      __builtin_memcpy (To, From, (Count) * sizeof (*(From)))
#  else
#   define YYCOPY(To, From, Count)		\
      do					\
	{					\
	  register YYSIZE_T yyi;		\
	  for (yyi = 0; yyi < (Count); yyi++)	\
	    (To)[yyi] = (From)[yyi];		\
	}					\
      while (0)
#  endif
# endif

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack)					\
    do									\
      {									\
	YYSIZE_T yynewbytes;						\
	YYCOPY (&yyptr->Stack, Stack, yysize);				\
	Stack = &yyptr->Stack;						\
	yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
	yyptr += yynewbytes / sizeof (*yyptr);				\
      }									\
    while (0)

#endif

#if defined (__STDC__) || defined (__cplusplus)
   typedef signed char yysigned_char;
#else
   typedef short yysigned_char;
#endif

/* YYFINAL -- State number of the termination state. */
#define YYFINAL  65
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   104

/* YYNTOKENS -- Number of terminals. */
#define YYNTOKENS  42
/* YYNNTS -- Number of nonterminals. */
#define YYNNTS  22
/* YYNRULES -- Number of rules. */
#define YYNRULES  61
/* YYNRULES -- Number of states. */
#define YYNSTATES  124

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   296

#define YYTRANSLATE(YYX) 						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const unsigned char yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,    36,    37,    38,    39,    40,    41
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const unsigned char yyprhs[] =
{
       0,     0,     3,     5,     6,     7,     8,    10,    12,    14,
      16,    21,    24,    27,    29,    32,    34,    36,    41,    44,
      47,    49,    52,    55,    56,    59,    60,    63,    66,    68,
      72,    76,    78,    80,    84,    87,    90,    93,    97,    99,
     103,   114,   118,   121,   123,   126,   130,   132,   138,   142,
     147,   150,   154,   158,   164,   166,   168,   171,   176,   180,
     184,   186
};

/* YYRHS -- A `-1'-separated list of the rules' RHS. */
static const yysigned_char yyrhs[] =
{
      43,     0,    -1,    62,    -1,    -1,    -1,    -1,    48,    -1,
      53,    -1,    49,    -1,    50,    -1,    44,     6,    51,     7,
      -1,    44,    52,    -1,    51,    52,    -1,    52,    -1,    45,
       3,    -1,    54,    -1,    55,    -1,    44,     6,    56,     7,
      -1,    44,    57,    -1,    56,    57,    -1,    57,    -1,    46,
       5,    -1,    59,    60,    -1,    -1,    61,    60,    -1,    -1,
       3,    47,    -1,    62,    63,    -1,    63,    -1,     8,     3,
      58,    -1,     9,     3,    58,    -1,    10,    -1,    11,    -1,
      12,     3,    58,    -1,    13,    53,    -1,    14,     3,    -1,
      15,     3,    -1,    16,     3,    58,    -1,    17,    -1,    18,
       3,    58,    -1,    19,     5,     5,     5,     5,     5,     5,
       5,     5,     5,    -1,    20,     3,    58,    -1,    21,     3,
      -1,    22,    -1,    23,     3,    -1,    24,     3,    58,    -1,
      25,    -1,    26,     5,     5,     5,     5,    -1,    27,     3,
      58,    -1,    28,     5,     5,     5,    -1,    29,     3,    -1,
      30,     3,    58,    -1,    31,     3,    58,    -1,    32,     3,
       3,     3,    58,    -1,    33,    -1,    34,    -1,    35,    54,
      -1,    36,     5,     5,     5,    -1,    38,     3,    58,    -1,
      37,     3,    58,    -1,    39,    -1,    40,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const unsigned short yyrline[] =
{
       0,   144,   144,   148,   158,   163,   168,   172,   177,   181,
     187,   192,   196,   199,   203,   209,   213,   218,   223,   227,
     230,   234,   240,   244,   249,   253,   256,   274,   277,   281,
     288,   295,   299,   303,   310,   316,   320,   324,   331,   335,
     342,   346,   353,   357,   361,   365,   372,   376,   380,   387,
     391,   395,   402,   409,   416,   420,   424,   430,   434,   441,
     448,   452
};
#endif

#if YYDEBUG || YYERROR_VERBOSE
/* YYTNME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals. */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "STRING", "ID", "NUM", "LBRACK", "RBRACK", 
  "ACCELERATOR", "AREALIGHTSOURCE", "ATTRIBUTEBEGIN", "ATTRIBUTEEND", 
  "CAMERA", "CONCATTRANSFORM", "COORDINATESYSTEM", "COORDSYSTRANSFORM", 
  "FILM", "IDENTITY", "LIGHTSOURCE", "LOOKAT", "MATERIAL", "OBJECTBEGIN", 
  "OBJECTEND", "OBJECTINSTANCE", "PIXELFILTER", "REVERSEORIENTATION", 
  "ROTATE", "SAMPLER", "SCALE", "SEARCHPATH", "SHAPE", 
  "SURFACEINTEGRATOR", "TEXTURE", "TRANSFORMBEGIN", "TRANSFORMEND", 
  "TRANSFORM", "TRANSLATE", "VOLUME", "VOLUMEINTEGRATOR", "WORLDBEGIN", 
  "WORLDEND", "HIGH_PRECEDENCE", "$accept", "start", "array_init", 
  "string_array_init", "num_array_init", "array", "string_array", 
  "real_string_array", "single_element_string_array", "string_list", 
  "string_list_entry", "num_array", "real_num_array", 
  "single_element_num_array", "num_list", "num_list_entry", "paramlist", 
  "paramlist_init", "paramlist_contents", "paramlist_entry", 
  "ri_stmt_list", "ri_stmt", 0
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const unsigned short yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,   287,   288,   289,   290,   291,   292,   293,   294,
     295,   296
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const unsigned char yyr1[] =
{
       0,    42,    43,    44,    45,    46,    47,    47,    48,    48,
      49,    50,    51,    51,    52,    53,    53,    54,    55,    56,
      56,    57,    58,    59,    60,    60,    61,    62,    62,    63,
      63,    63,    63,    63,    63,    63,    63,    63,    63,    63,
      63,    63,    63,    63,    63,    63,    63,    63,    63,    63,
      63,    63,    63,    63,    63,    63,    63,    63,    63,    63,
      63,    63
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const unsigned char yyr2[] =
{
       0,     2,     1,     0,     0,     0,     1,     1,     1,     1,
       4,     2,     2,     1,     2,     1,     1,     4,     2,     2,
       1,     2,     2,     0,     2,     0,     2,     2,     1,     3,
       3,     1,     1,     3,     2,     2,     2,     3,     1,     3,
      10,     3,     2,     1,     2,     3,     1,     5,     3,     4,
       2,     3,     3,     5,     1,     1,     2,     4,     3,     3,
       1,     1
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const unsigned char yydefact[] =
{
       0,     0,     0,    31,    32,     0,     3,     0,     0,     0,
      38,     0,     0,     0,     0,    43,     0,     0,    46,     0,
       0,     0,     0,     0,     0,     0,    54,    55,     3,     0,
       0,     0,    60,    61,     0,     2,    28,    23,    23,    23,
       5,    34,    15,    16,    35,    36,    23,    23,     0,    23,
      42,    44,    23,     0,    23,     0,    50,    23,    23,     0,
       0,    56,     0,    23,    23,     1,    27,    29,    25,    30,
      33,     5,     0,    18,    37,    39,     0,    41,    45,     0,
      48,     0,    51,    52,     0,     0,    59,    58,     3,    22,
      25,     5,    20,    21,     0,     0,    49,    23,    57,     4,
      26,     6,     8,     9,     7,    24,    17,    19,     0,    47,
      53,     4,     0,    11,     0,     4,    13,    14,     0,    10,
      12,     0,     0,    40
};

/* YYDEFGOTO[NTERM-NUM]. */
static const yysigned_char yydefgoto[] =
{
      -1,    34,    40,   112,    72,   100,   101,   102,   103,   115,
     113,    41,    42,    43,    91,    73,    67,    68,    89,    90,
      35,    36
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -111
static const yysigned_char yypact[] =
{
      54,     5,     6,  -111,  -111,     9,  -111,    11,    12,    14,
    -111,    16,    15,    20,    22,  -111,    23,    26,  -111,    25,
      28,    27,    30,    31,    32,    33,  -111,  -111,  -111,    34,
      35,    37,  -111,  -111,    41,    54,  -111,  -111,  -111,  -111,
      36,  -111,  -111,  -111,  -111,  -111,  -111,  -111,    38,  -111,
    -111,  -111,  -111,    40,  -111,    42,  -111,  -111,  -111,    43,
      36,  -111,    44,  -111,  -111,  -111,  -111,  -111,    45,  -111,
    -111,  -111,    46,  -111,  -111,  -111,    47,  -111,  -111,    48,
    -111,    49,  -111,  -111,    52,    51,  -111,  -111,  -111,  -111,
      45,    50,  -111,  -111,    53,    90,  -111,  -111,  -111,     1,
    -111,  -111,  -111,  -111,  -111,  -111,  -111,  -111,    91,  -111,
    -111,    92,    56,  -111,    93,    94,  -111,  -111,    95,  -111,
    -111,    97,    98,  -111
};

/* YYPGOTO[NTERM-NUM].  */
static const yysigned_char yypgoto[] =
{
    -111,  -111,   -28,  -111,  -111,  -111,  -111,  -111,  -111,  -111,
    -110,   -51,    71,  -111,  -111,   -67,   -36,  -111,   -40,  -111,
    -111,    69
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -6
static const yysigned_char yytable[] =
{
      60,   116,    69,    70,    92,   120,    -5,   111,    37,    38,
      74,    75,    39,    77,    44,    45,    78,    46,    80,    47,
      48,    82,    83,    49,   107,    50,    51,    86,    87,    52,
      53,    54,    55,    56,    57,    58,    59,   104,    63,    62,
      64,    65,    71,    76,    92,    79,    84,    81,    88,    85,
     105,    93,    94,    95,    96,    97,    98,   106,   108,   117,
      99,   110,     1,     2,     3,     4,     5,     6,     7,     8,
       9,    10,    11,    12,    13,    14,    15,    16,    17,    18,
      19,    20,    21,    22,    23,    24,    25,    26,    27,    28,
      29,    30,    31,    32,    33,   109,   114,    -5,   118,    61,
     121,   119,   122,   123,    66
};

static const unsigned char yycheck[] =
{
      28,   111,    38,    39,    71,   115,     5,     6,     3,     3,
      46,    47,     3,    49,     3,     3,    52,     3,    54,     3,
       5,    57,    58,     3,    91,     3,     3,    63,    64,     3,
       5,     3,     5,     3,     3,     3,     3,    88,     3,     5,
       3,     0,     6,     5,   111,     5,     3,     5,     3,     5,
      90,     5,     5,     5,     5,     3,     5,     7,     5,     3,
      88,    97,     8,     9,    10,    11,    12,    13,    14,    15,
      16,    17,    18,    19,    20,    21,    22,    23,    24,    25,
      26,    27,    28,    29,    30,    31,    32,    33,    34,    35,
      36,    37,    38,    39,    40,     5,     5,     5,     5,    28,
       5,     7,     5,     5,    35
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const unsigned char yystos[] =
{
       0,     8,     9,    10,    11,    12,    13,    14,    15,    16,
      17,    18,    19,    20,    21,    22,    23,    24,    25,    26,
      27,    28,    29,    30,    31,    32,    33,    34,    35,    36,
      37,    38,    39,    40,    43,    62,    63,     3,     3,     3,
      44,    53,    54,    55,     3,     3,     3,     3,     5,     3,
       3,     3,     3,     5,     3,     5,     3,     3,     3,     3,
      44,    54,     5,     3,     3,     0,    63,    58,    59,    58,
      58,     6,    46,    57,    58,    58,     5,    58,    58,     5,
      58,     5,    58,    58,     3,     5,    58,    58,     3,    60,
      61,    56,    57,     5,     5,     5,     5,     3,     5,    44,
      47,    48,    49,    50,    53,    60,     7,    57,     5,     5,
      58,     6,    45,    52,     5,    51,    52,     3,     5,     7,
      52,     5,     5,     5
};

#if ! defined (YYSIZE_T) && defined (__SIZE_TYPE__)
# define YYSIZE_T __SIZE_TYPE__
#endif
#if ! defined (YYSIZE_T) && defined (size_t)
# define YYSIZE_T size_t
#endif
#if ! defined (YYSIZE_T)
# if defined (__STDC__) || defined (__cplusplus)
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# endif
#endif
#if ! defined (YYSIZE_T)
# define YYSIZE_T unsigned int
#endif

#define yyerrok		(yyerrstatus = 0)
#define yyclearin	(yychar = YYEMPTY)
#define YYEMPTY		(-2)
#define YYEOF		0

#define YYACCEPT	goto yyacceptlab
#define YYABORT		goto yyabortlab
#define YYERROR		goto yyerrlab1


/* Like YYERROR except do call yyerror.  This remains here temporarily
   to ease the transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  */

#define YYFAIL		goto yyerrlab

#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)					\
do								\
  if (yychar == YYEMPTY && yylen == 1)				\
    {								\
      yychar = (Token);						\
      yylval = (Value);						\
      yytoken = YYTRANSLATE (yychar);				\
      YYPOPSTACK;						\
      goto yybackup;						\
    }								\
  else								\
    { 								\
      yyerror ("syntax error: cannot back up");\
      YYERROR;							\
    }								\
while (0)

#define YYTERROR	1
#define YYERRCODE	256

/* YYLLOC_DEFAULT -- Compute the default location (before the actions
   are run).  */

#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)         \
  Current.first_line   = Rhs[1].first_line;      \
  Current.first_column = Rhs[1].first_column;    \
  Current.last_line    = Rhs[N].last_line;       \
  Current.last_column  = Rhs[N].last_column;
#endif

/* YYLEX -- calling `yylex' with the right arguments.  */

#ifdef YYLEX_PARAM
# define YYLEX yylex (YYLEX_PARAM)
#else
# define YYLEX yylex ()
#endif

/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)			\
do {						\
  if (yydebug)					\
    YYFPRINTF Args;				\
} while (0)

# define YYDSYMPRINT(Args)			\
do {						\
  if (yydebug)					\
    yysymprint Args;				\
} while (0)

# define YYDSYMPRINTF(Title, Token, Value, Location)		\
do {								\
  if (yydebug)							\
    {								\
      YYFPRINTF (stderr, "%s ", Title);				\
      yysymprint (stderr, 					\
                  Token, Value);	\
      YYFPRINTF (stderr, "\n");					\
    }								\
} while (0)

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (cinluded).                                                   |
`------------------------------------------------------------------*/

#if defined (__STDC__) || defined (__cplusplus)
static void
yy_stack_print (short *bottom, short *top)
#else
static void
yy_stack_print (bottom, top)
    short *bottom;
    short *top;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (/* Nothing. */; bottom <= top; ++bottom)
    YYFPRINTF (stderr, " %d", *bottom);
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)				\
do {								\
  if (yydebug)							\
    yy_stack_print ((Bottom), (Top));				\
} while (0)


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

#if defined (__STDC__) || defined (__cplusplus)
static void
yy_reduce_print (int yyrule)
#else
static void
yy_reduce_print (yyrule)
    int yyrule;
#endif
{
  int yyi;
  unsigned int yylno = yyrline[yyrule];
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %u), ",
             yyrule - 1, yylno);
  /* Print the symbols being reduced, and their result.  */
  for (yyi = yyprhs[yyrule]; 0 <= yyrhs[yyi]; yyi++)
    YYFPRINTF (stderr, "%s ", yytname [yyrhs[yyi]]);
  YYFPRINTF (stderr, "-> %s\n", yytname [yyr1[yyrule]]);
}

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (yydebug)				\
    yy_reduce_print (Rule);		\
} while (0)

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YYDSYMPRINT(Args)
# define YYDSYMPRINTF(Title, Token, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef	YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   SIZE_MAX < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#if YYMAXDEPTH == 0
# undef YYMAXDEPTH
#endif

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif



#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined (__GLIBC__) && defined (_STRING_H)
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
static YYSIZE_T
#   if defined (__STDC__) || defined (__cplusplus)
yystrlen (const char *yystr)
#   else
yystrlen (yystr)
     const char *yystr;
#   endif
{
  register const char *yys = yystr;

  while (*yys++ != '\0')
    continue;

  return yys - yystr - 1;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined (__GLIBC__) && defined (_STRING_H) && defined (_GNU_SOURCE)
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
static char *
#   if defined (__STDC__) || defined (__cplusplus)
yystpcpy (char *yydest, const char *yysrc)
#   else
yystpcpy (yydest, yysrc)
     char *yydest;
     const char *yysrc;
#   endif
{
  register char *yyd = yydest;
  register const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

#endif /* !YYERROR_VERBOSE */



#if YYDEBUG
/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

#if defined (__STDC__) || defined (__cplusplus)
static void
yysymprint (FILE *yyoutput, int yytype, YYSTYPE *yyvaluep)
#else
static void
yysymprint (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE *yyvaluep;
#endif
{
  /* Pacify ``unused variable'' warnings.  */
  (void) yyvaluep;

  if (yytype < YYNTOKENS)
    {
      YYFPRINTF (yyoutput, "token %s (", yytname[yytype]);
# ifdef YYPRINT
      YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# endif
    }
  else
    YYFPRINTF (yyoutput, "nterm %s (", yytname[yytype]);

  switch (yytype)
    {
      default:
        break;
    }
  YYFPRINTF (yyoutput, ")");
}

#endif /* ! YYDEBUG */
/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

#if defined (__STDC__) || defined (__cplusplus)
static void
yydestruct (int yytype, YYSTYPE *yyvaluep)
#else
static void
yydestruct (yytype, yyvaluep)
    int yytype;
    YYSTYPE *yyvaluep;
#endif
{
  /* Pacify ``unused variable'' warnings.  */
  (void) yyvaluep;

  switch (yytype)
    {

      default:
        break;
    }
}


/* Prevent warnings from -Wmissing-prototypes.  */

#ifdef YYPARSE_PARAM
# if defined (__STDC__) || defined (__cplusplus)
int yyparse (void *YYPARSE_PARAM);
# else
int yyparse ();
# endif
#else /* ! YYPARSE_PARAM */
#if defined (__STDC__) || defined (__cplusplus)
int yyparse (void);
#else
int yyparse ();
#endif
#endif /* ! YYPARSE_PARAM */



/* The lookahead symbol.  */
int yychar;

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval;

/* Number of syntax errors so far.  */
int yynerrs;



/*----------.
| yyparse.  |
`----------*/

#ifdef YYPARSE_PARAM
# if defined (__STDC__) || defined (__cplusplus)
int yyparse (void *YYPARSE_PARAM)
# else
int yyparse (YYPARSE_PARAM)
  void *YYPARSE_PARAM;
# endif
#else /* ! YYPARSE_PARAM */
#if defined (__STDC__) || defined (__cplusplus)
int
yyparse (void)
#else
int
yyparse ()

#endif
#endif
{
  
  register int yystate;
  register int yyn;
  int yyresult;
  /* Number of tokens to shift before error messages enabled.  */
  int yyerrstatus;
  /* Lookahead token as an internal (translated) token number.  */
  int yytoken = 0;

  /* Three stacks and their tools:
     `yyss': related to states,
     `yyvs': related to semantic values,
     `yyls': related to locations.

     Refer to the stacks thru separate pointers, to allow yyoverflow
     to reallocate them elsewhere.  */

  /* The state stack.  */
  short	yyssa[YYINITDEPTH];
  short *yyss = yyssa;
  register short *yyssp;

  /* The semantic value stack.  */
  YYSTYPE yyvsa[YYINITDEPTH];
  YYSTYPE *yyvs = yyvsa;
  register YYSTYPE *yyvsp;



#define YYPOPSTACK   (yyvsp--, yyssp--)

  YYSIZE_T yystacksize = YYINITDEPTH;

  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;


  /* When reducing, the number of symbols on the RHS of the reduced
     rule.  */
  int yylen;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY;		/* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */

  yyssp = yyss;
  yyvsp = yyvs;

  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed. so pushing a state here evens the stacks.
     */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
	/* Give user a chance to reallocate the stack. Use copies of
	   these so that the &'s don't force the real ones into
	   memory.  */
	YYSTYPE *yyvs1 = yyvs;
	short *yyss1 = yyss;


	/* Each stack pointer address is followed by the size of the
	   data in use in that stack, in bytes.  This used to be a
	   conditional around just the two extra args, but that might
	   be undefined if yyoverflow is a macro.  */
	yyoverflow ("parser stack overflow",
		    &yyss1, yysize * sizeof (*yyssp),
		    &yyvs1, yysize * sizeof (*yyvsp),

		    &yystacksize);

	yyss = yyss1;
	yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyoverflowlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
	goto yyoverflowlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
	yystacksize = YYMAXDEPTH;

      {
	short *yyss1 = yyss;
	union yyalloc *yyptr =
	  (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
	if (! yyptr)
	  goto yyoverflowlab;
	YYSTACK_RELOCATE (yyss);
	YYSTACK_RELOCATE (yyvs);

#  undef YYSTACK_RELOCATE
	if (yyss1 != yyssa)
	  YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;


      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
		  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
	YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

/* Do appropriate processing given the current state.  */
/* Read a lookahead token if we need one and don't already have one.  */
/* yyresume: */

  /* First try to decide what to do without reference to lookahead token.  */

  yyn = yypact[yystate];
  if (yyn == YYPACT_NINF)
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = YYLEX;
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YYDSYMPRINTF ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yyn == 0 || yyn == YYTABLE_NINF)
	goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  /* Shift the lookahead token.  */
  YYDPRINTF ((stderr, "Shifting token %s, ", yytname[yytoken]));

  /* Discard the token being shifted unless it is eof.  */
  if (yychar != YYEOF)
    yychar = YYEMPTY;

  *++yyvsp = yylval;


  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  yystate = yyn;
  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     `$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 2:
#line 145 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"
    {
;}
    break;

  case 3:
#line 149 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"
    {
	if (cur_array) ArrayFree( cur_array );
	cur_array = new ParamArray;
	cur_array->allocated = 0;
	cur_array->nelems = 0;
	cur_array->array = NULL;
	array_is_single_string = false;
;}
    break;

  case 4:
#line 159 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"
    {
	cur_array->element_size = sizeof( const char * );
;}
    break;

  case 5:
#line 164 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"
    {
	cur_array->element_size = sizeof( float );
;}
    break;

  case 6:
#line 169 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"
    {
	yyval.ribarray = yyvsp[0].ribarray;
;}
    break;

  case 7:
#line 173 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"
    {
	yyval.ribarray = yyvsp[0].ribarray;
;}
    break;

  case 8:
#line 178 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"
    {
	yyval.ribarray = yyvsp[0].ribarray;
;}
    break;

  case 9:
#line 182 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"
    {
	yyval.ribarray = ArrayDup(cur_array);
	array_is_single_string = true;
;}
    break;

  case 10:
#line 188 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"
    {
	yyval.ribarray = ArrayDup(cur_array);
;}
    break;

  case 11:
#line 193 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"
    {
;}
    break;

  case 12:
#line 197 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"
    {
;}
    break;

  case 13:
#line 200 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"
    {
;}
    break;

  case 14:
#line 204 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"
    {
	char *to_add = strdup(yyvsp[0].string);
	AddArrayElement( &to_add );
;}
    break;

  case 15:
#line 210 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"
    {
	yyval.ribarray = yyvsp[0].ribarray;
;}
    break;

  case 16:
#line 214 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"
    {
	yyval.ribarray = ArrayDup(cur_array);
;}
    break;

  case 17:
#line 219 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"
    {
	yyval.ribarray = ArrayDup(cur_array);
;}
    break;

  case 18:
#line 224 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"
    {
;}
    break;

  case 19:
#line 228 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"
    {
;}
    break;

  case 20:
#line 231 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"
    {
;}
    break;

  case 21:
#line 235 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"
    {
	float to_add = yyvsp[0].num;
	AddArrayElement( &to_add );
;}
    break;

  case 22:
#line 241 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"
    {
;}
    break;

  case 23:
#line 245 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"
    {
	cur_paramlist_size = 0;
;}
    break;

  case 24:
#line 250 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"
    {
;}
    break;

  case 25:
#line 253 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"
    {
;}
    break;

  case 26:
#line 257 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"
    {
	void *arg = new char[ yyvsp[0].ribarray->nelems * yyvsp[0].ribarray->element_size ];
	memcpy(arg, yyvsp[0].ribarray->array, yyvsp[0].ribarray->nelems * yyvsp[0].ribarray->element_size);
	if (cur_paramlist_size >= cur_paramlist_allocated) {
		cur_paramlist_allocated = 2*cur_paramlist_allocated + 1;
		cur_paramlist_tokens = (const char **) realloc(cur_paramlist_tokens, cur_paramlist_allocated*sizeof(const char *) );
		cur_paramlist_args = (void * *) realloc( cur_paramlist_args, cur_paramlist_allocated*sizeof(void *) );
		cur_paramlist_sizes = (int *) realloc( cur_paramlist_sizes, cur_paramlist_allocated*sizeof(int) );
		cur_paramlist_texture_helper = (bool *) realloc( cur_paramlist_texture_helper, cur_paramlist_allocated*sizeof(bool) );
	}
	cur_paramlist_tokens[cur_paramlist_size] = yyvsp[-1].string;
	cur_paramlist_sizes[cur_paramlist_size] = yyvsp[0].ribarray->nelems;
	cur_paramlist_texture_helper[cur_paramlist_size] = array_is_single_string;
	cur_paramlist_args[cur_paramlist_size++] = arg;
	ArrayFree( yyvsp[0].ribarray );
;}
    break;

  case 27:
#line 275 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"
    {
;}
    break;

  case 28:
#line 278 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"
    {
;}
    break;

  case 29:
#line 282 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"
    {
	ParamSet params;
	InitParamSet(params, CPS, CPT, CPA, CPSZ, CPTH);
	pbrtAccelerator(yyvsp[-1].string, params);
	FreeArgs();
;}
    break;

  case 30:
#line 289 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"
    {
	ParamSet params;
	InitParamSet(params, CPS, CPT, CPA, CPSZ, CPTH);
	pbrtAreaLightSource(yyvsp[-1].string, params);
	FreeArgs();
;}
    break;

  case 31:
#line 296 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"
    {
	pbrtAttributeBegin();
;}
    break;

  case 32:
#line 300 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"
    {
	pbrtAttributeEnd();
;}
    break;

  case 33:
#line 304 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"
    {
	ParamSet params;
	InitParamSet(params, CPS, CPT, CPA, CPSZ, CPTH);
	pbrtCamera(yyvsp[-1].string, params);
	FreeArgs();
;}
    break;

  case 34:
#line 311 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"
    {
	if (VerifyArrayLength( yyvsp[0].ribarray, 16, "ConcatTransform" ))
		pbrtConcatTransform( (float *) yyvsp[0].ribarray->array );
	ArrayFree( yyvsp[0].ribarray );
;}
    break;

  case 35:
#line 317 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"
    {
	pbrtCoordinateSystem( yyvsp[0].string );
;}
    break;

  case 36:
#line 321 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"
    {
	pbrtCoordSysTransform( yyvsp[0].string );
;}
    break;

  case 37:
#line 325 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"
    {
	ParamSet params;
	InitParamSet(params, CPS, CPT, CPA, CPSZ, CPTH);
	pbrtFilm(yyvsp[-1].string, params);
	FreeArgs();
;}
    break;

  case 38:
#line 332 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"
    {
	pbrtIdentity();
;}
    break;

  case 39:
#line 336 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"
    {
	ParamSet params;
	InitParamSet(params, CPS, CPT, CPA, CPSZ, CPTH);
	pbrtLightSource(yyvsp[-1].string, params);
	FreeArgs();
;}
    break;

  case 40:
#line 343 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"
    {
	pbrtLookAt(yyvsp[-8].num, yyvsp[-7].num, yyvsp[-6].num, yyvsp[-5].num, yyvsp[-4].num, yyvsp[-3].num, yyvsp[-2].num, yyvsp[-1].num, yyvsp[0].num);
;}
    break;

  case 41:
#line 347 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"
    {
	ParamSet params;
	InitParamSet(params, CPS, CPT, CPA, CPSZ, CPTH);
	pbrtMaterial(yyvsp[-1].string, params);
	FreeArgs();
;}
    break;

  case 42:
#line 354 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"
    {
	pbrtObjectBegin(yyvsp[0].string);
;}
    break;

  case 43:
#line 358 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"
    {
	pbrtObjectEnd();
;}
    break;

  case 44:
#line 362 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"
    {
	pbrtObjectInstance(yyvsp[0].string);
;}
    break;

  case 45:
#line 366 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"
    {
	ParamSet params;
	InitParamSet(params, CPS, CPT, CPA, CPSZ, CPTH);
	pbrtPixelFilter(yyvsp[-1].string, params);
	FreeArgs();
;}
    break;

  case 46:
#line 373 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"
    {
	pbrtReverseOrientation();
;}
    break;

  case 47:
#line 377 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"
    {
	pbrtRotate(yyvsp[-3].num, yyvsp[-2].num, yyvsp[-1].num, yyvsp[0].num);
;}
    break;

  case 48:
#line 381 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"
    {
	ParamSet params;
	InitParamSet(params, CPS, CPT, CPA, CPSZ, CPTH);
	pbrtSampler(yyvsp[-1].string, params);
	FreeArgs();
;}
    break;

  case 49:
#line 388 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"
    {
	pbrtScale(yyvsp[-2].num, yyvsp[-1].num, yyvsp[0].num);
;}
    break;

  case 50:
#line 392 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"
    {
	pbrtSearchPath(yyvsp[0].string);
;}
    break;

  case 51:
#line 396 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"
    {
	ParamSet params;
	InitParamSet(params, CPS, CPT, CPA, CPSZ, CPTH);
	pbrtShape(yyvsp[-1].string, params);
	FreeArgs();
;}
    break;

  case 52:
#line 403 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"
    {
	ParamSet params;
	InitParamSet(params, CPS, CPT, CPA, CPSZ, CPTH);
	pbrtSurfaceIntegrator(yyvsp[-1].string, params);
	FreeArgs();
;}
    break;

  case 53:
#line 410 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"
    {
	ParamSet params;
	InitParamSet(params, CPS, CPT, CPA, CPSZ, CPTH);
	pbrtTexture(yyvsp[-3].string, yyvsp[-2].string, yyvsp[-1].string, params);
	FreeArgs();
;}
    break;

  case 54:
#line 417 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"
    {
	pbrtTransformBegin();
;}
    break;

  case 55:
#line 421 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"
    {
	pbrtTransformEnd();
;}
    break;

  case 56:
#line 425 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"
    {
	if (VerifyArrayLength( yyvsp[0].ribarray, 16, "Transform" ))
		pbrtTransform( (float *) yyvsp[0].ribarray->array );
	ArrayFree( yyvsp[0].ribarray );
;}
    break;

  case 57:
#line 431 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"
    {
	pbrtTranslate(yyvsp[-2].num, yyvsp[-1].num, yyvsp[0].num);
;}
    break;

  case 58:
#line 435 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"
    {
	ParamSet params;
	InitParamSet(params, CPS, CPT, CPA, CPSZ, CPTH);
	pbrtVolumeIntegrator(yyvsp[-1].string, params);
	FreeArgs();
;}
    break;

  case 59:
#line 442 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"
    {
	ParamSet params;
	InitParamSet(params, CPS, CPT, CPA, CPSZ, CPTH);
	pbrtVolume(yyvsp[-1].string, params);
	FreeArgs();
;}
    break;

  case 60:
#line 449 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"
    {
	pbrtWorldBegin();
;}
    break;

  case 61:
#line 453 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"
    {
	pbrtWorldEnd();
;}
    break;


    }

/* Line 999 of yacc.c.  */
#line 1743 "c:\\Downloads\\pbrt-src-1.02\\core\\/pbrtparse.cpp"

  yyvsp -= yylen;
  yyssp -= yylen;


  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;


  /* Now `shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*------------------------------------.
| yyerrlab -- here on detecting error |
`------------------------------------*/
yyerrlab:
  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if YYERROR_VERBOSE
      yyn = yypact[yystate];

      if (YYPACT_NINF < yyn && yyn < YYLAST)
	{
	  YYSIZE_T yysize = 0;
	  int yytype = YYTRANSLATE (yychar);
	  const char* yyprefix;
	  char *yymsg;
	  int yyx;

	  /* Start YYX at -YYN if negative to avoid negative indexes in
	     YYCHECK.  */
	  int yyxbegin = yyn < 0 ? -yyn : 0;

	  /* Stay within bounds of both yycheck and yytname.  */
	  int yychecklim = YYLAST - yyn;
	  int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
	  int yycount = 0;

	  yyprefix = ", expecting ";
	  for (yyx = yyxbegin; yyx < yyxend; ++yyx)
	    if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
	      {
		yysize += yystrlen (yyprefix) + yystrlen (yytname [yyx]);
		yycount += 1;
		if (yycount == 5)
		  {
		    yysize = 0;
		    break;
		  }
	      }
	  yysize += (sizeof ("syntax error, unexpected ")
		     + yystrlen (yytname[yytype]));
	  yymsg = (char *) YYSTACK_ALLOC (yysize);
	  if (yymsg != 0)
	    {
	      char *yyp = yystpcpy (yymsg, "syntax error, unexpected ");
	      yyp = yystpcpy (yyp, yytname[yytype]);

	      if (yycount < 5)
		{
		  yyprefix = ", expecting ";
		  for (yyx = yyxbegin; yyx < yyxend; ++yyx)
		    if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
		      {
			yyp = yystpcpy (yyp, yyprefix);
			yyp = yystpcpy (yyp, yytname[yyx]);
			yyprefix = " or ";
		      }
		}
	      yyerror (yymsg);
	      YYSTACK_FREE (yymsg);
	    }
	  else
	    yyerror ("syntax error; also virtual memory exhausted");
	}
      else
#endif /* YYERROR_VERBOSE */
	yyerror ("syntax error");
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
	 error, discard it.  */

      /* Return failure if at end of input.  */
      if (yychar == YYEOF)
        {
	  /* Pop the error token.  */
          YYPOPSTACK;
	  /* Pop the rest of the stack.  */
	  while (yyss < yyssp)
	    {
	      YYDSYMPRINTF ("Error: popping", yystos[*yyssp], yyvsp, yylsp);
	      yydestruct (yystos[*yyssp], yyvsp);
	      YYPOPSTACK;
	    }
	  YYABORT;
        }

      YYDSYMPRINTF ("Error: discarding", yytoken, &yylval, &yylloc);
      yydestruct (yytoken, &yylval);
      yychar = YYEMPTY;

    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*----------------------------------------------------.
| yyerrlab1 -- error raised explicitly by an action.  |
`----------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;	/* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (yyn != YYPACT_NINF)
	{
	  yyn += YYTERROR;
	  if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
	    {
	      yyn = yytable[yyn];
	      if (0 < yyn)
		break;
	    }
	}

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
	YYABORT;

      YYDSYMPRINTF ("Error: popping", yystos[*yyssp], yyvsp, yylsp);
      yydestruct (yystos[yystate], yyvsp);
      yyvsp--;
      yystate = *--yyssp;

      YY_STACK_PRINT (yyss, yyssp);
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  YYDPRINTF ((stderr, "Shifting error token, "));

  *++yyvsp = yylval;


  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

#ifndef yyoverflow
/*----------------------------------------------.
| yyoverflowlab -- parser overflow comes here.  |
`----------------------------------------------*/
yyoverflowlab:
  yyerror ("parser stack overflow");
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
  return yyresult;
}


#line 456 "c:\\Downloads\\pbrt-src-1.02\\core\\pbrtparse.y"

static void InitParamSet(ParamSet &ps, int count, const char **tokens,
		void **args, int *sizes, bool *texture_helper) {
	ps.Clear();
	for (int i = 0; i < count; ++i) {
		int type;
		string name;
		if (lookupType(tokens[i], &type, name)) {
			if (texture_helper && texture_helper[i] && type != PARAM_TYPE_TEXTURE && type != PARAM_TYPE_STRING)
			{
				Warning( "Bad type for %s. Changing it to a texture.", name.c_str());
				type = PARAM_TYPE_TEXTURE;
			}
			void *data = args[i];
			int nItems = sizes[i];
			if (type == PARAM_TYPE_INT) {
				// parser doesn't handle ints, so convert from floats here....
				int nAlloc = sizes[i];
				int *idata = new int[nAlloc];
				float *fdata = (float *)data;
				for (int j = 0; j < nAlloc; ++j)
					idata[j] = int(fdata[j]);
				ps.AddInt(name, idata, nItems);
				delete[] idata;
			}
			else if (type == PARAM_TYPE_BOOL) {
				// strings -> bools
				int nAlloc = sizes[i];
				bool *bdata = new bool[nAlloc];
				for (int j = 0; j < nAlloc; ++j) {
					string s(*((const char **)data));
					if (s == "true") bdata[j] = true;
					else if (s == "false") bdata[j] = false;
					else {
						Warning("Value \"%s\" unknown for boolean parameter \"%s\"."
							"Using \"false\".", s.c_str(), tokens[i]);
						bdata[j] = false;
					}
				}
				ps.AddBool(name, bdata, nItems);
				delete[] bdata;
			}
			else if (type == PARAM_TYPE_FLOAT) {
				ps.AddFloat(name, (float *)data, nItems);
			} else if (type == PARAM_TYPE_POINT) {
				ps.AddPoint(name, (Point *)data, nItems / 3);
			} else if (type == PARAM_TYPE_VECTOR) {
				ps.AddVector(name, (Vector *)data, nItems / 3);
			} else if (type == PARAM_TYPE_NORMAL) {
				ps.AddNormal(name, (Normal *)data, nItems / 3);
			} else if (type == PARAM_TYPE_COLOR) {
				ps.AddSpectrum(name, (Spectrum *)data, nItems / COLOR_SAMPLES);
			} else if (type == PARAM_TYPE_STRING) {
				string *strings = new string[nItems];
				for (int j = 0; j < nItems; ++j)
					strings[j] = string(*((const char **)data+j));
				ps.AddString(name, strings, nItems);
				delete[] strings;
			}
			else if (type == PARAM_TYPE_TEXTURE) {
				if (nItems == 1) {
					string val(*((const char **)data));
					ps.AddTexture(name, val);
				}
				else
					Error("Only one string allowed for \"texture\" parameter \"%s\"",
						name.c_str());
			}
		}
		else
			Warning("Type of parameter \"%s\" is unknown",
				tokens[i]);
	}
}
static bool lookupType(const char *token, int *type, string &name) {
	Assert(token != NULL);
	*type = 0;
	const char *strp = token;
	while (*strp && isspace(*strp))
		++strp;
	if (!*strp) {
		Error("Parameter \"%s\" doesn't have a type declaration?!", token);
		return false;
	}
	#define TRY_DECODING_TYPE(name, mask) \
		if (strncmp(name, strp, strlen(name)) == 0) { \
			*type = mask; strp += strlen(name); \
		}
	     TRY_DECODING_TYPE("float",    PARAM_TYPE_FLOAT)
	else TRY_DECODING_TYPE("integer",  PARAM_TYPE_INT)
	else TRY_DECODING_TYPE("bool",     PARAM_TYPE_BOOL)
	else TRY_DECODING_TYPE("point",    PARAM_TYPE_POINT)
	else TRY_DECODING_TYPE("vector",   PARAM_TYPE_VECTOR)
	else TRY_DECODING_TYPE("normal",   PARAM_TYPE_NORMAL)
	else TRY_DECODING_TYPE("string",   PARAM_TYPE_STRING)
	else TRY_DECODING_TYPE("texture",  PARAM_TYPE_TEXTURE)
	else TRY_DECODING_TYPE("color",    PARAM_TYPE_COLOR)
	else {
		Error("Unable to decode type for token \"%s\"", token);
		return false;
	}
	while (*strp && isspace(*strp))
		++strp;
	name = string(strp);
	return true;
}

