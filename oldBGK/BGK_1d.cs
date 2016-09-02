using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;


namespace BGK.Models
{
    /// <summary>
    /// BGK model solver class
    /// </summary>
    class BGK_1d
    {
         //---distribution function;
         //--first dimension - time; second - x coordinate; third - velocity
         public float[,,] F;
         //--sound velocity squared
         // private float cs = 1f;

         private float vw; // velocity weight
         //--invtau
         private float invtau;
         private const int tsteps = 2; // number of time steps, technically=2 for 1 step solver
         private int nspace; // number of cells
         private int nvel; // number of velocities 
         private float dx; // space step;
         private float dt; // time step;

         private float[,] Fnew; // time layer for one step solution

         //--constructor: define initial conditions
         public BGK_1d(float[,] Finit, float Kn, float vlcweight) {

             invtau = 1.0f / Kn; // Knudsen number
             vw = vlcweight; // dist. between velocities
             nspace = Finit.GetLength(0);
             nvel = Finit.GetLength(1);

             dx = 1.0f / (nspace - 1.0f); // space step
             dt =  dx / (vw * (nvel - 1.0f) * 0.5f) ; // time step (CFL condition)


             Fnew = new float[nspace, nvel];  // time layer for one step solution
             F = new float[tsteps, nspace, nvel]; // tsteps = 2

             //---set the initial distribution
             for (int j = 0; j < Finit.GetLength(0); j++) 
             { 
                 for(int k = 0; k < nvel; k++)
                 {
                     F[0, j, k] = Finit[j, k];
                 }
             }
         
         }

         private float _dens(int t, int j) // equilibrium density evaluation
         {
             float result_dens = 0;
             for (int k = 0; k < nvel; k++)
             {
                 result_dens += F[t, j, k];
             }
             return result_dens;
         }

         private float _vlc(int t, int j, float rho) // equilibrium velocity evaluation
         {
             if (rho == 0) //--!warning comparison with float zero
             {
                 return 0f;
             }
             else
             {
                 float result_vlc = 0;
                 for (int k = 0; k < nvel; k++)
                 {
                     result_vlc += F[t, j, k] * (k - (nvel - 1.0f) / 2.0f) * vw;
                 }
                 return result_vlc / rho;
             }
         }

         private float _temp(int t, int j, float rho, float u)
         { // equilibrium temperature evaluation
            float result_temp = 0;
            for (int k = 0; k < nvel; k++) {
                result_temp += ((k - (nvel - 1.0f) / 2.0f) * vw - u) *
                               ((k - (nvel - 1.0f) / 2.0f) * vw - u) * 
                               F[t, j, k];
            }
            if (rho == 0.0)
                return 0.0f;
            else
                return result_temp / rho;
        }

        /* normal (Maxwell) distribution function (1d) */
        /* warning, this function does not conserve discrete density,
         * so before using we shoud check density using function below !! */
         static public float _n_eq(float T, float rho, float u, float vw, int nvel, int k)
         {
             if (rho == 0.0f)
                 return 0.0f;
             double _tmp = rho / Math.Sqrt(2.0 * Math.PI * T) * Math.Exp(-
                 ((k - (nvel - 1.0) / 2.0) * vw - u) *
                 ((k - (nvel - 1.0) / 2.0) * vw - u) / (2.0 * T));
             return (float)_tmp;
         }

         // checks actual density on defined grid
         static public float _n_eq_check_dens(float T, float u, float vw, int nvel)
         {
             double _tmp = 0;
             for (int k = 0; k < nvel; ++k)
             {
                 _tmp += 1.0 / Math.Sqrt(2.0 * Math.PI * T) * Math.Exp(-
                     ((k - (nvel - 1.0) / 2.0) * vw - u) *
                     ((k - (nvel - 1.0) / 2.0) * vw - u) / (2.0 * T));
             }
             return (float)_tmp;
         }

         //-- time evolution for 1 step
         public float[,] Solve_one_step() 
         {

                 const int t = 1; // only one step considered

                 // copy initial values to Fnew
                 for (int k = 0; k < nvel; ++k) {
                     Fnew[0, k] = F[1, 0, k] = F[0, 0, k];
                     Fnew[nspace - 1, k] = F[1, nspace - 1, k] = F[0, nspace - 1, k];
                 }

                 // 1. Evaluation with no collision integral
                 for (int j = 1; j + 1 < nspace; ++j) {//--all cells (0, 1) 
                     // copy all values to new time layer 
                     for (int k = 0; k < nvel; ++k) {
                         F[1, j, k] = F[0, j, k];
                     }
                     
                     // positive velocities
                     for (int k = nvel / 2; k < nvel; ++k)
                         F[t, j, k] += (vw * (k - (nvel - 1.0f) / 2.0f)) * dt / dx 
                                     * (F[t - 1, j - 1, k] - F[t - 1, j, k]);
                     // negative velocities
                     for (int k = 0; k < nvel / 2; ++k)
                         F[t, j, k] -= (vw * (k - (nvel - 1.0f) / 2.0f)) * dt / dx 
                                     * (F[t - 1, j + 1, k] - F[t - 1, j, k]);

                 }
                 /* part of evaluation */
                 // positive velocities for x = 1
                 for (int j = nspace - 1, k = nvel / 2; k < nvel; ++k)
                     F[t, j, k] += (vw * (k - (nvel - 1.0f) / 2.0f)) * dt / dx
                                 * (F[t - 1, j - 1, k] - F[t - 1, j, k]);
                 // negative velocities for x = 0
                 for (int j = 0,          k = 0; k < nvel / 2; ++k)
                     F[t, j, k] -= (vw * (k - (nvel - 1.0f) / 2.0f)) * dt / dx
                                 * (F[t - 1, j + 1, k] - F[t - 1, j, k]);


                 // implicit (1, j), t or explicit (0, j), t-1 collision integral evaluation
                 for (int j = 0; j < nspace; ++j) // x in [0, 1]
                 {
                     float rho = _dens(1, j);
                     float u = _vlc(1, j, rho);
                     float T = _temp(1, j, rho, u);
                     float check_dens = _n_eq_check_dens(T, u, vw, nvel);

                     if (check_dens > 0.0f)
                     {
                         if (j > 0)
                             for (int k = nvel / 2; k < nvel; ++k)
                                 F[t, j, k] += invtau * (_n_eq(T, rho / check_dens, u, vw, nvel, k) - F[t, j, k]) * dt;
                         if (j + 1 < nspace)
                             for (int k = 0; k < nvel / 2; ++k)
                                 F[t, j, k] += invtau * (_n_eq(T, rho / check_dens, u, vw, nvel, k) - F[t, j, k]) * dt;
                     }

                     for (int k = 0; k < nvel; ++k)
                         Fnew[j, k] = F[t, j, k];
                 }

                 return Fnew;
         }//--end of the method


         //--evaluate density (x)
         public float[] P()
         {
             float[] dens = new float[nspace];
    
                 for (int j = 0; j < nspace; j++)
                     dens[j] = _dens(1, j);

             return dens;
         }//--end of the method

         //--evaluate bulk velocity (x)
         public float[] U()
         {
             float[] u = new float[nspace];

             for (int j = 0; j < nspace; j++)
                 u[j] = _vlc(1, j, _dens(1, j));

             return u;
         }//--end of the method


         //--evaluate temperature (x)
         public float[] Temperatures()
         {
             float[] T = new float[nspace];

             for (int j = 0; j < nspace; j++)
                 T[j] = _temp(1, j, _dens(1, j), _vlc(1, j, _dens(1, j)));

             return T;
         }//--end of the method

     }

}
