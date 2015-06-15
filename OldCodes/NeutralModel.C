// NEUTRAL MODEL
// PRINTS THEORETICAL VS. COMPUTED COVER

/**** headers ****/
#include <iostream>
//#include "../../util/ran2_NR.h"
// check the random generator
using namespace std;

/**** globals ****/
const int   EMP       = 0;
const int   VEG       = 1;
const int   W         = 100;	// grid size
const int   NP        = W*W;
const int   MXT       = 1000;	// time length
const int   RMX       = 100;	// independent runs
const float PF        = 0.1;	// initial (random) cover
const float M         = 0.25;	// m parameter

/**** main routine ****/
main(){
   int i, j, n, t, l, **grid;
   float r, p, ntot, cover;
   seed();			// random number generator seeding
   grid = new int* [W];
   for(i=0;i<W;i++) grid[i] = new int [W];

   // initialize grid
   for(i=0;i<W;i++)
      for(j=0;j<W;j++)
         if(randm()<PF) grid[i][j] = VEG;
         else 		grid[i][j] = EMP;

   // begin dynamics
   for(p=0.4;p<0.801;r+=0.05){	// r parameter loop
   
      r = p*M/(1.0 - p);	// computing r from known cover 
      cover = 0.0;
      for(n=1;n<=RMX;n++){	// independent run loop
   	 for(t=1;t<=MXT;t++){	// time loop
      	    for(l=0;l<NP;l++){       // asynchronous updating begins
	       i = int(W*randm());
	       j = int(W*randm());
	       if(grid[i][j]==EMP){
	          if(randm()<r) grid[i][j] = VEG;	// EMP --> VEG
               }
               else if(randm()<M) grid[i][j] = EMP;	// VEG --> EMP
            }
         }

         // compute cover
	 ntot = 0.0;
	 for(i=0;i<W;i++)
	    for(j=0;j<W;j++) if(grid[i][j]==VEG) ntot++;
         ntot /= NP;
	 cover += ntot;
      }
      cout << "cover(p) = " << p << ", cover(obs) = " << cover/RMX << "\n";
   }

   return 0;
}

