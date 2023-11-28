#ifdef _SGI
#define dznrm2 dznrm2_
#define zdotc zdotc_
#define zaxpy zaxpy_
#define zreadmtc zreadmtc_
#define zcsrcsc zcsrcsc_
#define zroscal zroscal_
#define zcoscal zcoscal_
#define qsplit qsplit_
#else
#ifdef _LINUX
#define dznrm2 dznrm2_
#define zdotc zdotc_
#define zaxpy zaxpy_
#define zreadmtc zreadmtc_
#define zcsrcsc zcsrcsc_
#define zroscal zroscal_
#define zcoscal zcoscal_
#define qsplit qsplit_
#else
#ifdef _IBM 
#define dznrm2 dznrm2
#define zdotc zdotc
#define zaxpy zaxpy
#define zreadmtc zreadmtc
#define zcsrcsc zcsrcsc
#define zroscal zroscal
#define zcoscal zcoscal
#define qsplit qsplit
#else
#define dznrm2 dznrm2_
#define zdotc zdotc_
#define zaxpy zaxpy_
#define zreadmtc zreadmtc_
#define zcsrcsc zcsrcsc_
#define zroscal zroscal_
#define zcoscal zcoscal_
#define qsplit qsplit_
#endif
#endif
#endif

#ifndef min
#define min(a,b) (((a)>(b))?(b):(a))
#endif
#ifndef max
#define max(a,b) (((a)>(b))?(a):(b))
#endif

#define MAX_LINE        256
#define MAX_HBNAME      64
#define MAX_MAT	   100


