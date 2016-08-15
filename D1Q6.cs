using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;


namespace Latticeb.Models
{
    /// <summary>
    /// D1Q6 model solver class
    /// </summary>
    class D1Q6
    {
         //---distribution function;
         //--first dimension - time; second - x coordinate; third - velocity
         //-- delta t == 1.0
         //-- delta x == 0.5
         public float[,,] F;
         public float[,] Fnew; //--to return result
         //--sound velocity squared
         // private float cs = 1f;
         //--invtau
         private float invtau;
         private int tsteps;
         private int N;

         //--constructor: define initial conditions
         public D1Q6(float[,] Finit, float ttau) {
             invtau = 1f / ttau;
             tsteps = 1;
             N = Finit.GetLength(0);
             F = new float[tsteps + 1, N, 6];
             Fnew = new float[N, 6];

             //---set the initial distribution
             for (int j = 0; j < N; j++) 
             { 
                 for(int k = 0; k < 6; k++)
                 {
                     F[0, j, k] = Finit[j, k];
                 }
             }
         
         }

         private float _dens(int t, int j)
         {
             return F[t, j, 0] + F[t, j, 1] +
                    F[t, j, 2] + F[t, j, 3] +
                    F[t, j, 4] + F[t, j, 5];
         }

         private float _vlc(int t, int j, float rho)
         {
             if (rho == 0) //--!warning comparison with float zero
             {
                 return 0f;
             }
             else
             {
                 return (F[t, j, 4] - F[t, j, 1]
                 + 0.5f * (F[t, j, 3] - F[t, j, 2])
                 + 3.0f * (F[t, j, 5] - F[t, j, 0])) / rho;
             }
         }

         static public float _n_eq_k(float rho, float u, int k)
         {
             float c_i, w_i;
             const float cs = 1.0f;
             switch (k) {
                 case 0: c_i = -3.0f; w_i = 0.0303f; break;
                 case 1: c_i = -1.0f; w_i = 0.1463f; break;
                 case 2: c_i = -0.5f; w_i = 0.3234f; break;
                 case 3: c_i = +0.5f; w_i = 0.3234f; break;
                 case 4: c_i = +1.0f; w_i = 0.1463f; break;
                 case 5: c_i = +3.0f; w_i = 0.0303f; break;
                 default: c_i = 0.0f; w_i = 0.0f; break;
             }
             return rho * w_i * (1 + (u * c_i / cs) + 0.5f * (c_i * c_i - cs) * u * u / (cs * cs));
         }

         private float _n_eq_0(float rho, float u) { return _n_eq_k(rho, u, 0); }
         private float _n_eq_1(float rho, float u) { return _n_eq_k(rho, u, 1); }
         private float _n_eq_2(float rho, float u) { return _n_eq_k(rho, u, 2); }
         private float _n_eq_3(float rho, float u) { return _n_eq_k(rho, u, 3); }
         private float _n_eq_4(float rho, float u) { return _n_eq_k(rho, u, 4); }
         private float _n_eq_5(float rho, float u) { return _n_eq_k(rho, u, 5); }

         //-- time evolution of distribution function
         public float[,] Solve() 
         {
             //--left border
             for (int t = 0; t < tsteps; t++) //--time steps
             {
                 float rho, u; //--density and velocity in current cell

                 F[t + 1, 0, 5] = F[t, 0, 5];
                 F[t + 1, 1, 5] = F[t, 0, 5];
                 F[t + 1, 2, 5] = F[t, 0, 5];
                 F[t + 1, 3, 5] = F[t, 0, 5];
                 F[t + 1, 4, 5] = F[t, 0, 5];
                 F[t + 1, 5, 5] = F[t, 0, 5];

                 F[t + 1, 0, 4] = F[t, 0, 4];
                 F[t + 1, 1, 4] = F[t, 0, 4];
               
                 F[t + 1, 0, 3] = F[t, 0, 3];

                 F[t + 1, N - 1, 0] = F[t, N - 1, 0];   
                 F[t + 1, N - 2, 0] = F[t, N - 1, 0];
                 F[t + 1, N - 3, 0] = F[t, N - 1, 0];
                 F[t + 1, N - 4, 0] = F[t, N - 1, 0];
                 F[t + 1, N - 5, 0] = F[t, N - 1, 0];
                 F[t + 1, N - 6, 0] = F[t, N - 1, 0];

                 F[t + 1, N - 1, 0] = F[t, N - 1, 1];
                 F[t + 1, N - 2, 0] = F[t, N - 1, 1];

                 F[t + 1, N - 1, 0] = F[t, N - 1, 2];


                 for (int j = 0; j < 6; ++j) //--cells on the left side
                 {
                     //--velocities 0 - 5 for j in [3, n - 2)
                     rho = -_dens(t, j);
                     u = _vlc(t, j, rho);
                     F[t + 1, j + 6, 5] = invtau * (_n_eq_5(rho, u) - F[t, j, 5]) + F[t, j, 5];
                     F[t + 1, j + 2, 4] = invtau * (_n_eq_4(rho, u) - F[t, j, 4]) + F[t, j, 4];
                     F[t + 1, j + 1, 3] = invtau * (_n_eq_3(rho, u) - F[t, j, 3]) + F[t, j, 3];
                     if (j >= 1) F[t + 1, j - 1, 2] = invtau * (_n_eq_2(rho, u) - F[t, j, 2]) + F[t, j, 2];
                     if (j >= 2) F[t + 1, j - 2, 1] = invtau * (_n_eq_1(rho, u) - F[t, j, 1]) + F[t, j, 1];
                     // if (j >= 6) F[t + 1, j - 6, 0] = invtau * (_n_eq_0(rho, u) - F[t, j, 0]) + F[t, j, 0];
                 }

                 for (int j = 6; j < N - 6; ++j) //--all internal cells
                 {
                     //--velocities 0 - 5 for j in [3, n - 2)
                     rho = -_dens(t, j);
                     u = _vlc(t, j, rho);
                     F[t + 1, j - 6, 0] = invtau * (_n_eq_0(rho, u) - F[t, j, 0]) + F[t, j, 0];
                     F[t + 1, j - 2, 1] = invtau * (_n_eq_1(rho, u) - F[t, j, 1]) + F[t, j, 1];
                     F[t + 1, j - 1, 2] = invtau * (_n_eq_2(rho, u) - F[t, j, 2]) + F[t, j, 2];
                     F[t + 1, j + 1, 3] = invtau * (_n_eq_3(rho, u) - F[t, j, 3]) + F[t, j, 3];
                     F[t + 1, j + 2, 4] = invtau * (_n_eq_4(rho, u) - F[t, j, 4]) + F[t, j, 4];
                     F[t + 1, j + 6, 5] = invtau * (_n_eq_5(rho, u) - F[t, j, 5]) + F[t, j, 5];
                 }

                 for (int j = N - 6; j < N; ++j) //--cells on the right side
                 {
                     //--velocities 0 - 5 for j in [3, n - 2)
                     rho = -_dens(t, j);
                     u = _vlc(t, j, rho);
                     F[t + 1, j - 6, 0] = invtau * (_n_eq_0(rho, u) - F[t, j, 0]) + F[t, j, 0];
                     F[t + 1, j - 2, 1] = invtau * (_n_eq_1(rho, u) - F[t, j, 1]) + F[t, j, 1];
                     F[t + 1, j - 1, 2] = invtau * (_n_eq_2(rho, u) - F[t, j, 2]) + F[t, j, 2];
                     if (j + 1 < N) F[t + 1, j + 1, 3] = invtau * (_n_eq_3(rho, u) - F[t, j, 3]) + F[t, j, 3];
                     if (j + 2 < N) F[t + 1, j + 2, 4] = invtau * (_n_eq_4(rho, u) - F[t, j, 4]) + F[t, j, 4];
                     // if (j + 6 < N) F[t + 1, j + 6, 5] = invtau * (_n_eq_5(rho, u) - F[t, j, 5]) + F[t, j, 5];
                 }
             }
             //--there is no boundary conditions

             for (int j = 0; j < N; j++)
             {
                 for (int k = 0; k < 6; k++)
                 {
                     Fnew[j, k] = F[1, j, k];
                 }
             }
             return Fnew;
         } //--end of the method

         //--evaluate density (t,x)
         public float[,] P()
         {
             float[,] dens = new  float[F.GetLength(0) , F.GetLength(1)];
          
             for (int t = 0; t < F.GetLength(0); t++)
             {
                 for (int j = 0; j < F.GetLength(1); j++)
                 {

                     dens[t, j] = _dens(t, j);
                 }
             }
             return dens;
         }//--end of the method

         //--evaluate bulk velocity (t,x)
         public float[,] U()
         {
             float[,] u = new float[F.GetLength(0), F.GetLength(1)];

             for (int t = 0; t < F.GetLength(0); t++)
             {
                 for (int j = 0; j < F.GetLength(1); j++)
                 {
                     u[t, j] = _vlc(t, j, _dens(t, j));
                 }
             }
             return u;
         }//--end of the method


     }

}
