/* 
   com_abort.c, based on version from

   Parallel Communications Library

   Jim Teresco

   Rensselaer Polytechnic Institute
   Department of Computer Science
   Scientific Computation Research Center

   Adapted from PMDB, JDT, Wed Nov 22 15:36:53 EST 1995

   Last change: 
   Mon Nov 27 15:08:24 EST 1995

   $Id$
   $Log: com_abort.c,v $
   Revision 1.1  2006/03/06 05:28:23  terescoj
   First prelim version

   Revision 1.1.1.1  2002/06/15 01:04:29  terescoj
   CVSing scorec util library

   Revision 1.2  2000/08/14 19:02:23  oklaas
   update to latest version of com

   Revision 1.1.1.1  1996/12/18 14:53:28  jteresco
   initial sources

*/

#include <stdio.h>
#include <stdlib.h>

void com_abort(char *function, char *msg) {

  if (function || msg) {
    fflush(stdout);
    fprintf(stderr,"ABORT in %s!\n",function);
    if (msg)
      fprintf(stderr,"  ErrMsg: %s\n",msg);
  }
  abort();
}
