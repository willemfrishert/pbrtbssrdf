typedef union {
char string[1024];
float num;
ParamArray *ribarray;
} YYSTYPE;
#define	STRING	258
#define	ID	259
#define	NUM	260
#define	LBRACK	261
#define	RBRACK	262
#define	ACCELERATOR	263
#define	AREALIGHTSOURCE	264
#define	ATTRIBUTEBEGIN	265
#define	ATTRIBUTEEND	266
#define	CAMERA	267
#define	CONCATTRANSFORM	268
#define	COORDINATESYSTEM	269
#define	COORDSYSTRANSFORM	270
#define	FILM	271
#define	IDENTITY	272
#define	LIGHTSOURCE	273
#define	LOOKAT	274
#define	MATERIAL	275
#define	OBJECTBEGIN	276
#define	OBJECTEND	277
#define	OBJECTINSTANCE	278
#define	PIXELFILTER	279
#define	REVERSEORIENTATION	280
#define	ROTATE	281
#define	SAMPLER	282
#define	SCALE	283
#define	SEARCHPATH	284
#define	SHAPE	285
#define	SURFACEINTEGRATOR	286
#define	TEXTURE	287
#define	TRANSFORMBEGIN	288
#define	TRANSFORMEND	289
#define	TRANSFORM	290
#define	TRANSLATE	291
#define	VOLUME	292
#define	VOLUMEINTEGRATOR	293
#define	WORLDBEGIN	294
#define	WORLDEND	295
#define	HIGH_PRECEDENCE	296


extern YYSTYPE yylval;
