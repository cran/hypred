/* hypredRoutines.c --- 
 * 
 * Filename: hypredRoutines.c
 * Description: C routines to be used by functions in the package hypred
 * Author: Frank Technow
 * Maintainer: Frank Technow
 * Created: Fr Sep 17 13:24:30 2010 (+0200)
 * Version: 0.4
 * Last-Updated: Mo Nov 20 10:59:09 2013 (+0100)
 *           By: Technow
 *     Update #: 19
 * URL: 
 * Keywords: 
 * Compatibility: 
 * 
 */

/* Commentary: 
 * 
 *  Currently this file contains two routines, hypredCode_FUN_dom
 *  which codes a design matrix and meiosisFUNallChr which simulates
 *  the meiosis. Additionally this file contains some registration
 *  code.
 * 
 */

/* Change Log:
 * 
 *  V0.3 - array lCx was initiated with zero length in some
 *  cases. This is fixed, the minimum length is now one.
 *
 *  V0.4 - (*(SNPaftCx + numNextCx - 1) > *numLoci) was evaluated even
 *  if numNextCx was larger than nCx, fixed now
 *
 */

/* This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 3, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street, Fifth
 * Floor, Boston, MA 02110-1301, USA.
 */

/* Code: */


#include <stdio.h>
#include <string.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Applic.h> 
#include <R_ext/Rdynload.h> 






/***************************************************
 *
 * meiosisFUNallChr
 *
 ****************************************************/

/*
 *
 * Commentary: 
 * 
 *  This function produces 1 haploid progeny chromosome set from two homologous
 *  parental chromosomes sets. The function is basically identical to
 *  the function meiosisFUNCL, just that it now processed multiple chromosomes
 *
 *  The arguments are:
 * 
 *   - homologe1: an integer vector with one of the homologous
 *     parental chromosome sets
 *
 *   - homologe2: an integer vector with the other of the homologous
 *     parental chromosome sets
 *   
 *   - newHaplotype: an integer vector with the newly formed progeny
 *     haploid chromosome set
 *   
 *   - mapPos: the positions of the SNP on the genetic map in M (double
 *     vector) 
 *
 *   - numLoci: the number of loci on one chromosome(integer)
 *
 *   - numChr: the number of chromosomes 
 *
 *   - lenChr: the length of the chromosomes in M (double)
 *
 *
 *  Abbrevations in comments:
 *
 *  - Cx = crossing over
 *  - CL = count-location model
 *  - M  = Morgan
 *  - L  = length of the chromosome in M
 *
 *
 * Notes on the code:
 * 
 * 1. program structure (repeated for each chromosome):
 *
 * Part1: sample the recombination parameters
 *
 * Part2: find the SNP right after a Cx
 *
 * Part3: create the new haplotype by copying from the parental
 * homologes 
 *  	    
 * 1.  if multiple Cx occurred between two SNP X and Y the following
 *    will happen:
 *
 *    - code in Part2 will detect when mapPos Y is larger than the
 *      current and the next Cx. In this situation, countSNP will
 *      stay the same for one more loop (by decrementing it
 *      before the obligatory incrementing). Then everything will run
 *      as normal, just as if there were no multiple Cx. 
 *
 *    - in Part3: SNPaftCx will contain multiple entries of Y, as the SNP right
 *      behind the Cx (see above).
 *
 *    - when countSNP is < Y, everything will go on as normal 
 *      and the while loop below will run through until 
 *      countSNP + 1 == Y. 
 
 *    - at this point the Cx is performed by changing the 
 *      pointer to the other homologe and the location and id 
 *      of the next Cx is updated, all as normal 

 *    - then countSNP == Y and because Y >= Y (which gives the id of
 *      the "next Cx", the while loop can't be initiated and no copying
 *      of SNPs is performed and countSNP IS NOT INCREMENTED. But the
 *      pointer is still changed to the other homologe, which is the
 *      original homologe. Hence a double Cx occurred between X and Y
 *      and hence no recombination.
 
 *    - the id of the next Cx (numNextCx) is changed as well.
 
 *    - since countSNP was not incremented, still countSNP == 
 *      Y. So in the next iteration the same while condition is 
 *      evaluated again. If the next Cx is behind Y, the while 
 *      loop can be executed, and the program continues from Y, 
 *      just as if no Cx occurred. If, however, there were 3 Cx 
 *      between X and Y, still Y >= Y and the while loop can 
 *      not be initiated. Then however, the pointer is changed 
 *      again to the other homologe, so a triple Cx occurred 
 *      between X and Y and hence a recombination. 
 
 *    - Then in the next iteration, the program will continue 
 *      from Y, just as if only one Cx occurred. 
 
 *    - In short: if there was an even number of Cx between X 
 *      and Y, no recombination will occur and if there was 
 *      an odd number one recombination. ALL OF THIS IS JUST AS 
 *      IT SHOULD BE! 
 *
 */



/* #define DEBUG 			/\* turn debugging on *\/ */
#undef  DEBUG			/* turn debugging of */


void meiosisFUNallChr (int *homologe1,	            /* one parental homologe
						       chromosome sets */
		       int *homologe2,	            /* the other parental
						       homologe chromosome sets */
		       int *newHaplotype,           /* the newly formed progeny
						       haplotype */
		       double *mapPos,              /* positions of SNP on
						       chromosome in M */
		       int *numChr,	            /* the number of chromosomes  */
		       double *lenChr,              /* the length of the
						       chromosome in M  */
		       int *numLoci	            /* the number of loci */
		   )

{

  /* initiate some variables */
  
  double u;			/* random uniform[0,1] derivate  */

  int pointingAt;		/* variable indicating which homologe
				   is pointed at */
  int *currentHomologe_ptr;	/* pointer to the current homologe */

  int countSNP;	        	/* while loop counter for SNP */

  int nCx;			/* number of Cx */

  int countChr;			/* counter for the chromosomes */

  for(countChr = 0; countChr < *numChr; ++countChr)
    {
      
      /***************************************************
       *
       * Part 1
       *
       ****************************************************/


      /* sample the recombination parameters */
  
      /* number of Cx (nCx) on random meiotic product */
      GetRNGstate();
      nCx = rpois(*(lenChr + countChr));
      PutRNGstate();

  
      /* sample the locations of the Cx (and sort them) */
      
      int position_dummy;
      if(nCx > 0)
	position_dummy = 0;
      else
	position_dummy = 1;
      
      // vector of Cx positions (positition_dummy=1 is added so that
      // the array has length > 0 even if nCx == 0 (the element in the
      // last position is ignored)
      double lCx[nCx + position_dummy];     
      	

      if(nCx > 0)
	{

	  int countCx;
	  for(countCx = 0; countCx < nCx; ++countCx)
	    {
	      GetRNGstate();
	      lCx[countCx] = runif(0,*(lenChr + countChr));
	      PutRNGstate();
	    } /* end for */
	  /* sort */
	  R_rsort(lCx, nCx);
	}	  /* end if */


      /*    randomly chose the homologe that will be followed through the */
      /*    recombination process  */
      GetRNGstate();
      u = runif(0,1);
      PutRNGstate();
      /*  if u < 0.5 point to homologe1 else to homologe2 */
      if(u < 0.5)
	{
	  currentHomologe_ptr = homologe1;
	  pointingAt = 1;
	}
      else
	{
	  currentHomologe_ptr = homologe2;
	  pointingAt = 2;
	}

      /***************************************************
       *
       * Part 2
       *
       ****************************************************/



      if(nCx > 0)			/* only if there were Cx */
	{
  
	  int SNPaftCx[nCx];		/* stores the index of the SNP
					   right after the Cx */

  
	  double locNextCx;    /* the location of the next Cx */

	  locNextCx = lCx[0];
  
	  /* number of the Cx ahead (at start, this is */
	  /* 1, because we are at the beginning of */
	  /* the chromosome, if Cx 1 is passed */
	  /* countCx will change to 2, because Cx */
	  /* number 2 is now ahead) */
  
	  int numNextCx = 1;

	  /* loop over the whole chromosome, forever or until a break is
	     encountered */

	  countSNP = 0;

	  while(1)
	    {
	      /* if the ID of the Cx ahead is larger than nCx, break, */
	      /* because the last Cx was already passed */
	      if(numNextCx > nCx)
		{
		  break;
		}
	
	      /* if countSNP > (numLoci - 1), break. This situation occurs when */
	      /* there are Cx behind the last SNP, in this situation the if */
	      /* statement above (numNextCx > nCx) == FALSE since there are */
	      /* still nCx left, not including the following statement would */
	      /* result in an error further down, since countSNP is raised */
	      /* above numLoci, which is the length of mapPos. The SNP */
	      /* behind this Cx will be at a position that is non existing. */
	      /* The code in the copy part will detect this and act */
	      /* accordingly in such a way that the Cx behind the last SNP */
	      /* are ignored, as they should be, since they could not be */
	      /* observed. */
	      if(countSNP > (*numLoci - 1))
		{
		  *(SNPaftCx + numNextCx-1) = countSNP;
		  break;
		}
      
	      /* if the map position of the current SNP is behind the */
	      /* location of the next Cx */
      
	      if(*(mapPos + countSNP) > locNextCx)
		{
		  /* store the position of this SNP */
		  *(SNPaftCx + numNextCx-1) = countSNP;
		  	      
		  /* change the ID of the next Cx that is ahead */
		  ++numNextCx;

		  /* break if we are beyond the last Cx*/
		  if(numNextCx > nCx)
		    {
		      break;
		    }
            
		  /* change to the location of the next Cx */
		  locNextCx = lCx[numNextCx - 1];


		  /* if the mapPos of the SNP is also larger than the next Cx */
		  /* location then decrement countSNP. Outside this if block */
		  /* countSNP will incremented again, so that countSNP stays the */
		  /* same, until the mapPos of countSNP is not larger than the */
		  /* next Cx anymore. This assures that multiple Cx between two */
		  /* SNP are detected. */
		  if(countSNP <= (*numLoci - 1)){
		    if(*(mapPos + countSNP) > locNextCx)
		      --countSNP;
		  }
		} /* end if */
  
	      /* go to the next SNP */
	      ++countSNP;
	    } /* end while */

#ifdef DEBUG
	  int i;
	  for(i = 0; i < nCx; ++i)
	    {
	      Rprintf("SNPaftCx %d\n", *(SNPaftCx + i));
	    }
#endif	/* DEBUG */


	  /***************************************************
	   *
	   * Part 3
	   *
	   ****************************************************/


	  /*    start a while loop over the chromosome */

	  numNextCx = 1;

	  countSNP = 0;

	  while(TRUE)
	    {
	      /* copy all SNP between the last and the next Cx from the */
	      /* currentHomologe */
        
	      /* if the last Cx was passed copy until the end of the */
	      /* chromosome OR */
  
	      /* if the SNP behind the next Cx is a SNP that is out of the */
	      /* range of existing SNP (this happens when there are Cx after */
	      /* the last existing SNP) then ignore the Cx and copy until */
	      /* the end of the chromosome as well */
        
	      if((numNextCx > nCx) || (*(SNPaftCx + numNextCx - 1) > *numLoci))
		{
		  while(countSNP <= (*numLoci - 1)  )
		    {
		      *(newHaplotype + countSNP) = *(currentHomologe_ptr + countSNP); 
		      ++countSNP;
		    }
		  /* and break */
		  break;
		}
	      else
		/* else  copy until the next Cx; */  
		{
		  while(countSNP <  *(SNPaftCx + numNextCx - 1))
		    {
		      *(newHaplotype + countSNP) = *(currentHomologe_ptr + countSNP); 
		      ++countSNP;
		    }
		} /* end else */

	      /* after the next Cx, point to the other homologe */
	      if(pointingAt == 1)
		{
		  currentHomologe_ptr = homologe2;
		  pointingAt = 2;
		}
	      else
		{
		  currentHomologe_ptr = homologe1;
		  pointingAt = 1;
		}
	  
	      /* change the id of the next Cx */
	      ++numNextCx;


	    } /* end while */
	} /* end if nCx > 0 */
      /* if no Cx occured at all, the newHaplotype is just one of the
	 parental homologes*/
      if(nCx == 0)
	{
	  int countSNP;
	  for(countSNP = 0; countSNP < (*numLoci); ++countSNP)
	    {
	      *(newHaplotype + countSNP) = *(currentHomologe_ptr + countSNP);
	    } /* end for */
      
	}	  /* end if */

#ifdef DEBUG
      Rprintf("nCx %d\n", nCx);

      int i;

      for(i = 0; i < nCx; ++i)
	{
	  Rprintf("lCx %f\n", lCx[i]);
	}
#endif	/* DEBUG */

            
      /* move the pointers to next chromosomes */

      if(countChr != (*numChr - 1)) /* don't in the last loop */
	{
	  homologe1 += *numLoci;
	  homologe2 += *numLoci;
	  newHaplotype += *numLoci;
	  mapPos += *numLoci;
	}

    }       /* end for countChr */
  
}	  /* end meiosisFUNallChr */


/***************************************************
 *
 * hypredCode_FUN_dom
 *
 ****************************************************/


/* Commentary: 
 * 
 * This function is intended to be called from R. 
 *
 * The arguments are: 
 *
 * - genotypes: integer matrix with the SNP genotypes (individuals in
     rows, SNPs in columns)
 *
 * - designX: the design matrix for add. effects (double)
 *
 * - designW: the design matrix for dom. effects (double)
 *
 * - nSNP: number of markers (integer)
 *
 * - N: number of records (integer)
 *
 *
 * The coding is as follows (see Xu, 2003):
 *
 *    x = sqrt(2), w = -1 for 11;
 *    x = 0, w = 1 for 10 and x = -sqrt(2), w = -1 for 00
 *
 *
 * code design: 
 *
 * - start at the first SNP, loop over each individual for this SNP
 *
 * - look for its state
 *
 * - code the X and W accordingly
 *
 * - go to the next SNP
 *
 * The variable countSNP will code for the SNP by individual
 * combination and index the matrices designX and designW, the
 * variables first_allele and second_allele will index the first and
 * second allele of a given SNP and individual in the matrix genotypes
 */


void hypredCode_FUN_dom (int *genotypes,      /* the raw SNP genotypes */
			 double *designX,     /* add. design
						 matrix */
			 double *designW,     /* dom design matrix */
			 int *nSNP, 	      /* number of markers */
			 int *N	              /* number of records */
			   )
{
  int count_snp;		/* variable for individual by SNP combination
				   of genotypes */
  int first_allele = 0;	        /* subscript for first allele of a given SNP
			           in a given individual */
  int second_allele = 1;	/* subscript for second allele of a given SNP
			           in a given individual */
  
  /* loop over genotypes (column major wise) */
  for(count_snp = 0; count_snp < *nSNP * *N; ++count_snp)
    {
      /* if 11 */
      if((*(genotypes + first_allele) == 1) & (*(genotypes + second_allele) == 1))
	{
	  *(designX + count_snp) = M_SQRT2;
	  *(designW + count_snp) =  -1.0;
	} /* 00 */
      else if((*(genotypes + first_allele) == 0 )& (*(genotypes + second_allele) == 0))
	{
	  *(designX + count_snp) = -M_SQRT2;
	  *(designW + count_snp) = -1.0;
	}
      else 			/* 01 or 10 */
	{
	  *(designX + count_snp) = 0;
	  *(designW + count_snp) = 1.0;
	}
	
      first_allele += 2; 
      second_allele += 2;
    } /* end count_snp */

} /* end hypredCode_FUN_dom */





/***************************************************
 *
 * registration part
 *
 ****************************************************/

static R_NativePrimitiveArgType
meiosisFUNallChr_type[7] = {INTSXP, INTSXP, INTSXP, REALSXP, INTSXP, REALSXP, INTSXP };

static R_NativePrimitiveArgType
hypredCode_FUN_dom_type[5] = {INTSXP, REALSXP, REALSXP, INTSXP, INTSXP };
      

static const
R_CMethodDef cMethods[] = {
  {"meiosisFUNallChr", (DL_FUNC) &meiosisFUNallChr, 7, meiosisFUNallChr_type},
  {"hypredCode_FUN_dom", (DL_FUNC) &hypredCode_FUN_dom, 5, hypredCode_FUN_dom_type},
  {NULL, NULL, 0}
};

void R_init_hypred(DllInfo *info)
{
  /* register the .C routines, no .Call, .Fortran of .External
     routines, so pass those arrays as NULL */
  R_registerRoutines(info,
		     cMethods,
		     NULL, NULL, NULL);
}


void
R_unload_hypred(DllInfo *info)
{
  /* Release resources. */
}


/***************************************************
 * end registration part
 ****************************************************/



/* hypredRoutines.c ends here */
