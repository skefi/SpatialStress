// FACILITATION MODEL
// PRINTS B VS. COVER

/**** headers ****/
#include <iostream>
#include "../../util/ran2_NR.h"
using namespace std;

/**** globals ****/
const int   DES       = 0;
const int   EMP       = 1;
const int   VEG       = 2;
const int   NEIGHBORS = 4;
const int   W         = 100;
const int   NP        = W*W;
const int   MXT       = 1000;
const int   RMX       = 100;
const float PF        = 0.2;
const float PS        = 0.5;
const int   X[]       = {1,  0, -1,  0};
const int   Y[]       = {0,  1,  0, -1};
const float DEL       = 0.1;
const float DEL2      = 1.0 - DEL;
const float M         = 0.1;
const float D         = 0.2;
const float R         = 0.0001;
const float F         = 0.9;
const float C         = 0.3;

/**** main routine ****/
main(){
   char mf='n';		// mean-field(='y'), spatial(='n')
   int i, j, n, t, i1, j1, k, l, t2, **grid, count;
   float B, b, m, d, r, f, rn, rt, ntot, max_rate, cover;
   seed();
   grid = new int* [W];
   for(i=0;i<W;i++) grid[i] = new int [W];

   for(B=0.4;B<=0.71;B+=0.05){
      // normalize by the largest process rate
      max_rate = R + F;
      if(max_rate<(D + B)){
	 max_rate = D + B;
	 if(max_rate<M) max_rate = M;
      }
      else if(max_rate<M) max_rate = M;
      m = M/max_rate;
      b = B/max_rate;
      d = D/max_rate;
      r = R/max_rate;
      f = F/max_rate;

      // initialize grid
      for(i=0;i<W;i++)
	 for(j=0;j<W;j++){
         rn = randm();
         if(rn<PF)          grid[i][j] = EMP;
         else if(rn>=1.-PS) grid[i][j] = VEG;
         else               grid[i][j] = DES;
      }

      // begin dynamics
      cover = 0.0;
      for(n=1;n<=RMX;n++){
      
         for(t=1;t<=MXT;t++)
	    for(l=0;l<NP;l++){       // asynchronous updating begins
	       i = int(W*randm());
	       j = int(W*randm());

	       /**** DES --> EMP transition (facilitation)  ****/
       	       if(grid[i][j]==DES){
	          count = 0;
	          for(k=0;k<NEIGHBORS;k++){
	       	     if(mf=='y'){		// if mean-field
                        i1 = int(W*randm());
		        j1 = int(W*randm());
	             }
	             else{			// if spatial
                        i1 = (i + X[k] + W) % W;
                        j1 = (j + Y[k] + W) % W;
	             }
                     if(grid[i1][j1]==VEG) count++;
	          }
                  if(randm()<(r + (f*count)/NEIGHBORS)) grid[i][j] = EMP;
               }

               /**** EMP --> DES/VEG transition ****/
               else if(grid[i][j]==EMP){
	          rn = randm();
	          if(rn<d) grid[i][j] = DES;		// degeneration
	          else{					// colonization
                     count = 0;
                     for(k=0;k<NEIGHBORS;k++){
		        if(mf=='y'){			// if mean-field
                           i1 = int(W*randm());
                           j1 = int(W*randm());
		        }
		        else{				// if neighborhood
                           i1 = (i + X[k] + W) % W;
                           j1 = (j + Y[k] + W) % W;
		        }
                        if(grid[i1][j1]==VEG) count++;
	             }
	             rt = d+(DEL*ntot+(DEL2*count)/NEIGHBORS)*(b-C*ntot); 
	             if(rn<rt) grid[i][j] = VEG;
                  }
               }

	       /**** VEG --> EMP transition (death) ***/
               else if(randm()<m) grid[i][j] = EMP;
            }

         // compute vegetation cover
         ntot = 0.0;
         for(i=0;i<W;i++)
  	    for(j=0;j<W;j++) if(grid[i][j]==VEG) ntot++;
         cover += ntot/NP;
      }
      
      cout << "b = " << B << ",  cover = " << cover/RMX << "\n";
   }

   return 0;
}

