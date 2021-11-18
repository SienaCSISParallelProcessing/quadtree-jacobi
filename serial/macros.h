/*
   macros.h

   Based on util_macros originally written by

   Louis Ziantz
   Jim Teresco

   Rensselaer Polytechnic Institute
   Department of Computer Science
   Scientific Computation Research Center

   Macros which are generally useful.

   Created:
   Wed Nov 29 17:05:11 EST 1995

   Last change:
   Mon Jan 22 13:02:18 EST 1996

   $Id$

   Modification History
   07/16/1998  Added in FIELD_SET and OUTPUT_IN_ORDER from pmdb
*/

#ifndef __H_MACROS
#define __H_MACROS

#include <stdlib.h>
#include <stdio.h>
#include <errno.h>

#ifndef FALSE
#define FALSE 0
#endif
#ifndef TRUE
#define TRUE  1
#endif

extern void com_abort(char *, char *);

#define SAFE_MALLOC(v,type,size) \
 {  v = (type) malloc(size) ; \
    if ( v == NULL) { \
       fflush(stdout); \
       fprintf(stderr,"in file %s, line %d, failed to allocate %ld bytes",\
               __FILE__,__LINE__,size); \
       com_abort(NULL,NULL); \
    } \
 }

#define SAFE_REALLOC(v,type,size) \
 {  v = (type) realloc(v,size) ; \
    if ( v == NULL) { \
       fflush(stdout); \
       fprintf(stderr,"in file %s, line %d, failed to reallocate %ld bytes",\
               __FILE__,__LINE__,size); \
       com_abort(NULL,NULL); \
    } \
 }

#define GRACEFULLY_OPEN(fp, filename, mode) {                      \
  fp = fopen(filename,mode) ;                                      \
  if ( fp == NULL ) {                                              \
    char msg[FILENAME_MAX+25];                                     \
    sprintf(msg, "Could not open file %s.\n", filename);           \
    perror(" fopen");                                              \
    com_abort("GRACEFULLY_OPEN", msg);                             \
  }                                                                \
}

#define FILE_EXISTS(fp,filename)                                   \
        ( ((fp = fopen(filename,"r")) != NULL)                     \
        ? (! fclose(fp)) /* TRUE */ : FALSE) 

/*
 * Macro for checking return from fprintf
 *
 * Should only be used if routine doesn't need to clean up before
 * returning.  Calling routine should decare "int ret=0;" and do a
 * "return(ret)" at the end.  This should not be used in main, as it
 * will result in the program terminating and returning 0 on failure.
 *
 */

#define SAFE_FPRINTF  if (ret==EOF) return(0); else ret=fprintf

/*
 * Macro for checking return from fclose
 *
 * Should only be used if routine doesn't need to clean up before
 * returning.  Calling routine must return an int, 0 for failure.
 * It is the caller's repsonsibility to abort if necessary.
 *
 */


#define SAFE_FCLOSE(fp, fname) \
  if (fclose(fp)) { \
      fprintf(stderr,"error closing file %s\n",fname); \
      perror(" fclose"); \
      return(0); \
    }

#define GRACEFULLY_CLOSE(fp, fname) \
  if (fclose(fp)) { \
      char msg[FILENAME_MAX+25]; \
      sprintf(msg,"Error closing file %s\n",fname); \
      perror(" fclose"); \
      com_abort("GRACEFULLY_CLOSE", msg); \
    }

#define ASSERT(cond) \
  if (!(cond)) { \
      fflush(stdout); \
      fprintf(stderr,"ASSERT failed in file %s at line %d\n", \
         __FILE__,__LINE__); \
      com_abort(NULL,NULL); \
  }

#define FAIL \
      fflush(stdout); \
      fprintf(stderr,"FAIL in file %s at line %d\n", \
         __FILE__,__LINE__); \
      com_abort(NULL,NULL)

/* do not add anything below this line! */
#endif /* __H_MACROS */
