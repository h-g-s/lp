#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <float.h>
#include <string.h>
#include <math.h>
extern "C"
{
#include "../lp.h"
}

#define EPS 1e-5

int main(int argc, const char *argv[])
{
   LinearProgram *lp = lp_create();

   double obj[] = { -1.0, 1.0 };
   double lb[]  = {  0.0, 0.0 };
   double ub[]  = {  DBL_MAX, DBL_MAX };
   char isInteger[] = { 1, 1 };

   char vNames[][64] = { "v1", "v2" };
   char *pstr[2];
   pstr[0] = vNames[0];
   pstr[1] = vNames[1];
   lp_add_cols( lp, 2, obj, lb, ub, isInteger, pstr );

   double A[8][2] = { { 1.0,  1.0 },
                      {-1.0,  5.0 },
                      { 2.5,  1.0 },
                      { 1.0, -5.5 },
                      {-1.0,  1.0 },
                      { 8.5,  1.0 },
                      { 2.5,  5.5 },
                      {-1.5,  5.0 } };

   char sense[8] = { 'G', 'L', 'L', 'L', 'L', 'L', 'G', 'L' };
   double rhs[8] = { 3.5, 22.5, 19.5, 2.5, 2.5, 58.5, 12, 17.75 };
   int idx[] = { 0, 1 };

   for (int i = 0; (i<8) ; i++)
   {
      char rName[64];
      sprintf( rName, "rest%d", i+1 );

      lp_add_row( lp, 2, idx, A[i], rName, sense[i], rhs[i] );
   }

   lp_write_lp( lp, "out1.lp" );
   
   const double *obj1 = lp_obj_coef( lp );
   assert( obj1[0] == -1 );
   assert( obj1[1] == 1 );
    
   /* checking query */
   {
       int cidx[8]; double ccoef[8];
       int cnz = lp_col( lp, 0, cidx, ccoef );
       assert( cnz == 8 );
       for ( int i=0; (i<8) ; i++ ) 
           assert( A[cidx[i]][0] == ccoef[i] );

       int nz = lp_row( lp, 3, cidx, ccoef);
       assert( nz == 2 );
       assert( cidx[0]  == 0 || cidx[0] == 1);
       assert( cidx[1]  == 0 || cidx[1] == 1);
       assert( cidx[0]  != cidx[1] );
       assert( ccoef[0] == 1.0 || ccoef[0] == -5.5 );
       assert( ccoef[1] == 1.0 || ccoef[1] == -5.5 );
       assert( ccoef[0] != ccoef[1] );
       assert( lp_rhs(lp,3)==2.5 );
       assert( lp_sense(lp,3)=='L' );
   }
    
   int status = lp_optimize( lp );

   switch (status) 
   {
      case LP_OPTIMAL:
         break;
      default:
         fprintf( stderr, "status should be optimal. exiting\n");
         abort();
   }


   printf("Obj value: %g\n", lp_obj_value(lp) );

   assert( fabs(-5-lp_obj_value(lp))<EPS );

   double expectedX[] = { 6, 1 };

   for ( int i=0 ; (i<2) ; i++ )
   {
      char name[256];
      printf( "%s : %g\n", lp_col_name(lp,i,name), lp_x(lp)[i] );
      assert( fabs(expectedX[i] - lp_x(lp)[i]) < EPS );
   }

   status = lp_optimize_as_continuous( lp );
   
   switch (status) 
   {
      case LP_OPTIMAL:
         break;
      default:
         fprintf( stderr, "status should be optimal. exiting\n");
         abort();
   }

   for ( int i=0 ; (i<8) ; i++ )
   {
      char name[256];
      printf( "%s : %g\n", lp_row_name(lp,i,name), lp_row_price(lp)[i] );
   }

   lp_write_lp( lp, "a.lp" );

   printf("Obj value: %g\n", lp_obj_value(lp) );
   assert( fabs(-6.01047-lp_obj_value(lp))<EPS );

   status = lp_optimize( lp );
   switch (status) 
   {
      case LP_OPTIMAL:
         break;
      default:
         fprintf( stderr, "status should be optimal. exiting\n");
         abort();
   }
   printf("Obj value: %g\n", lp_obj_value(lp) );
   assert( fabs(-5.0-lp_obj_value(lp))<EPS );

   lp_set_col_bounds( lp, 1, 2.0, 3.0 );

   lp_write_lp( lp, "b.lp" );
   status = lp_optimize( lp );
   switch (status) 
   {
      case LP_OPTIMAL:
         break;
      default:
         fprintf( stderr, "status should be optimal. exiting\n");
         abort();
   }

   printf("Obj value: %g\n", lp_obj_value(lp) );
   assert( fabs(-4.0-lp_obj_value(lp))<EPS );

   assert( fabs( lp_x(lp)[0]-6.0 ) < EPS );
   printf( "x[1]: %.8f\n", lp_x(lp)[1] );
   assert( lp_x(lp)[1]>=1.999 && lp_x(lp)[1]<=3.0001);

   assert( fabs(lp_row_slack(lp)[0]-4.5)<=1e-8 );
   assert( fabs(lp_row_slack(lp)[2]-2.5)<=1e-8 );

   // getting obj
   double *objv = new double[ lp_cols(lp) ];
   memcpy( objv, lp_obj_coef(lp), sizeof(double)*lp_cols(lp) );

   // adding obj as constraint
   {
       int *oidx = new int[lp_cols(lp)];
       double *ocoef = new double[lp_cols(lp)];
       int nz = 0;
       for ( int i=0 ; (i<lp_cols(lp)) ; ++i )
           if (fabs(objv[i])>1e-10)
           {
               oidx[nz] = i;
               ocoef[nz] = objv[i];
               ++nz;
           }

       lp_add_row( lp, nz, oidx, ocoef, "objcut", 'G', -2.0 );

       status  = lp_optimize( lp  );

       assert( status == LP_OPTIMAL );

       assert( lp_obj_value(lp) >= -2.0 -1e-10 );

       lp_write_lp( lp, "c.lp" );

       int nridx[] =     {0, 1, 0 };
       double nrcoef[] = {7, 9, 1 };
       int starts[] = { 0, 2, 3 };
       double nrrhs[] = { 61, 3 };
       char sense[] = { 'L', 'L' };

       char n1[] = "newr1";
       char n2[] = "newr2";

        char *nrnames[2] = { n1, n2 };

       lp_add_rows( lp, 2, starts, nridx, nrcoef, sense, nrrhs, (const char **) nrnames );

       lp_write_lp( lp, "d.lp" );

       status = lp_optimize( lp );

       assert( fabs(lp_obj_value(lp) + 1.0)<1e-10 );

       // removing rows
       int nrr = 6;
       int rr[] = { 0, 2, 4, 6, 8, 10 };

       lp_remove_rows( lp, nrr, rr );

       status = lp_optimize( lp );

       assert( fabs(lp_obj_value(lp) + 4.0)<1e-10 );

       lp_write_lp( lp, "e.lp" );

       delete[] oidx;
       delete[] ocoef;
   }

   delete[] objv;

   lp_free( &lp );

   lp = lp_create();
   lp_read( lp, "abz5.lp");
   lp_set_max_nodes( lp, 500 );
   lp_set_max_seconds( lp, 30 );
   lp_set_cuts( lp, 1 );
   lp_set_heur_fp_passes( lp, 100 );
   lp_set_heur_proximity( lp, 0 );
   lp_set_max_saved_sols( lp, 50 );
   lp_set_print_messages( lp, 1 );
   lp_optimize( lp );
   /*
   printf("%d solutions were found.\n", lp_num_saved_sols(lp) );
   for (int i = 0; i<lp_num_saved_sols(lp) ; i++) 
   {
      printf("checking validity of solution %d\n", i);
      LinearProgram *lpfix = lp_clone( lp );
      for ( int j=0 ; (j<lp_cols(lpfix)) ; ++j )
      {
         if ( lp_is_integer( lpfix, j ) )
         {
            double val = lp_saved_sol_x( lp, i )[j];
            lp_set_col_bounds( lpfix, j, val, val );
         }
      }
      lp_free( &lpfix );
   }
   */
   lp_free( &lp );

   /*
   lp = lp_create();
   lp_read( lp, "abz5.lp");
   lp_fix_col( lp, 0, lp_col_ub( lp, 0 ) );
 //  lp_fix_col( lp, 1, lp_col_ub( lp, 1 ) );
   LinearProgram *lpp = lp_pre_process( lp );

   printf( "After pre-processing columns went from %d to %d\n", lp_cols( lp ), lp_cols( lpp ) );

   lp_free( &lp );
   lp_free( &lpp ); */

   printf("test completed with success.\n");
   
   return 0;
}

