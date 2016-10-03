using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;


namespace latticeb
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
        public float[, ,] F;
        public float[,] Fnew; //-- array to be returned as result
        //--sound velocity squared
        // private float cs = 1f;
        //--invtau
        private float invtau;
        private int tsteps;
        private int N;

        //--constructor: define initial conditions
        public D1Q6(float[,] Finit, float ttau)
        {
            invtau = 1f / ttau;
            tsteps = 1;
            N = Finit.GetLength(0);
            F = new float[tsteps + 1, N, 6];
            Fnew = new float[N, 6];

            //---set the initial distribution
            for (int j = 0; j < N; j++)
            {
                for (int k = 0; k < 6; k++)
                {
                    F[0, j, k] = Finit[j, k];
                }
            }

        }

        // Model parameter constants
        const float w0 = 0.3234f;
        const float w1 = 0.1463f;
        const float w2 = 0.0303f;

        /* inverse matrix constants
        {{0.3234 * 1,    0.1463 * 1,  0.0303 * 1},
         {0.3234 * 0.5,  0.1463 * 1,  0.0303 * 3},
         {0.3234 * 0.25, 0.1463 * 1,  0.0303 * 9}}^-1 = 
         * 
         {{7.42115, -9.89487, 2.47372},
         {-10.2529, 23.9234, -6.83527},
         {3.30033, -9.90099, 6.60066}}
          wolframalpha.com 
         */

      /*const float a11 = 7.42115f; const float a12 = -9.89487f; const float a13 = 2.47372f;
        const float a21 = -10.2529f; const float a22 = 23.9234f; const float a23 = -6.83527f;
        const float a31 = 3.30033f; const float a32 = -9.90099f; const float a33 = 6.60066f;
        */

        /*Oleg's correction  -  one should use the Vandermonde matrix  (NO weights w1 w2 w3 in the matrix), see below

         M=[
            1      1    1;
            0.5    1    3;
            0.5^2  1^2  3^2;
           
         inv(M)=[
             2.4    -3.2   0.8;   
            -1.5     3.5    -1;
             0.1    -0.3   0.2
            ]
         */

        const float a11 = 2.4f; const float a12 =-3.2f; const float a13 = 0.8f;
        const float a21 =-1.5f; const float a22 = 3.5f; const float a23 = -1.0f;
        const float a31 = 0.1f; const float a32 = -0.3f; const float a33 = 0.2f;

        // distribution and moment evaluation functions
        static public float distr0(float m0, float m1, float m2)
        {
            return a11 * m0 + a12 * m1 + a13 * m2;
        }

        static public float distr1(float m0, float m1, float m2)
        {
            return a21 * m0 + a22 * m1 + a23 * m2;
        }

        static public float distr2(float m0, float m1, float m2)
        {
            return a31 * m0 + a32 * m1 + a33 * m2;
        }


        static public float moment0(float v0, float v1, float v2)
        {
            return w0 * 1.0f * v0 + w1 * v1 + w2 * 1.0f * v2;
        }

        static public float moment1(float v0, float v1, float v2)
        {
            return w0 * 0.5f * v0 + w1 * v1 + w2 * 3.0f * v2;
        }

        static public float moment2(float v0, float v1, float v2)
        {
            return w0 * 0.25f * v0 + w1 * v1 + w2 * 9.0f * v2;
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


        private float _tmpr(int t, int j, float rho)
        {
            if (rho == 0) //--!warning comparison with float zero
            {
                return 0f;
            }
            else
            {
                return (F[t, j, 4] + F[t, j, 1]
                + 0.25f * (F[t, j, 3] + F[t, j, 2])
                + 9.0f * (F[t, j, 5] + F[t, j, 0])) / rho - _vlc(t, j, rho) * _vlc(t, j, rho);
            }
        }


        static public float _n_eq_k(float rho, float u, int k) // equilibrium distribution
        {
            float c_i, w_i;
            const float cs = 1.0f;
            switch (k)
            {
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

        // equilibrium distibution function for each velocity
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

                // this part of code might be reviewed, better variant of with compuatation may be applied 
                rho = _dens(t, 0);
                u = _vlc(t, 0, rho);

                const float flag_a = 0.0f; // near border correction flag (some affect found)
                const float flag_b = 0.0f; // border correction flag (no high affect found)

                F[t + 1, 0, 5] = F[t, 0, 5] + flag_b * invtau * (_n_eq_5(rho, u) - F[t, 0, 5]);
                F[t + 1, 1, 5] = F[t, 0, 5] + flag_a * invtau * (_n_eq_5(rho, u) - F[t, 0, 5]);
                F[t + 1, 2, 5] = F[t, 0, 5] + flag_a * invtau * (_n_eq_5(rho, u) - F[t, 0, 5]);
                F[t + 1, 3, 5] = F[t, 0, 5] + flag_a * invtau * (_n_eq_5(rho, u) - F[t, 0, 5]);
                F[t + 1, 4, 5] = F[t, 0, 5] + flag_a * invtau * (_n_eq_5(rho, u) - F[t, 0, 5]);
                F[t + 1, 5, 5] = F[t, 0, 5] + flag_a * invtau * (_n_eq_5(rho, u) - F[t, 0, 5]);

                F[t + 1, 0, 4] = F[t, 0, 4] + flag_b * invtau * (_n_eq_4(rho, u) - F[t, 0, 4]);
                F[t + 1, 1, 4] = F[t, 0, 4] + flag_a * invtau * (_n_eq_4(rho, u) - F[t, 0, 4]);

                F[t + 1, 0, 3] = F[t, 0, 3] + flag_b * invtau * (_n_eq_3(rho, u) - F[t, 0, 3]);

                rho = _dens(t, N - 1);
                u = _vlc(t, N - 1, rho);

                F[t + 1, N - 1, 0] = F[t, N - 1, 0] + flag_b * invtau * (_n_eq_0(rho, u) - F[t, N - 1, 0]);
                F[t + 1, N - 2, 0] = F[t, N - 1, 0] + flag_a * invtau * (_n_eq_0(rho, u) - F[t, N - 1, 0]);
                F[t + 1, N - 3, 0] = F[t, N - 1, 0] + flag_a * invtau * (_n_eq_0(rho, u) - F[t, N - 1, 0]);
                F[t + 1, N - 4, 0] = F[t, N - 1, 0] + flag_a * invtau * (_n_eq_0(rho, u) - F[t, N - 1, 0]);
                F[t + 1, N - 5, 0] = F[t, N - 1, 0] + flag_a * invtau * (_n_eq_0(rho, u) - F[t, N - 1, 0]);
                F[t + 1, N - 6, 0] = F[t, N - 1, 0] + flag_a * invtau * (_n_eq_0(rho, u) - F[t, N - 1, 0]);

                F[t + 1, N - 1, 1] = F[t, N - 1, 1] + flag_b * invtau * (_n_eq_1(rho, u) - F[t, N - 1, 1]);
                F[t + 1, N - 2, 1] = F[t, N - 1, 1] + flag_a * invtau * (_n_eq_1(rho, u) - F[t, N - 1, 1]);

                F[t + 1, N - 1, 2] = F[t, N - 1, 2] + flag_b * invtau * (_n_eq_2(rho, u) - F[t, N - 1, 2]);


                for (int j = 0; j < 6; ++j) //--cells on the left side
                {
                    //--velocities 0 - 5 for j in [3, n - 2)
                    rho = _dens(t, j);
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
                    rho = _dens(t, j);
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
                    rho = _dens(t, j);
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
            float[,] dens = new float[F.GetLength(0), F.GetLength(1)];

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


        //--evaluate  temperature (t,x)
        public float[,] T()
        {
            float[,] temper = new float[F.GetLength(0), F.GetLength(1)];

            for (int t = 0; t < F.GetLength(0); t++)
            {
                for (int j = 0; j < F.GetLength(1); j++)
                {
                    temper[t, j] = _tmpr(t, j, _dens(t, j));
                }
            }
            return temper;
        }//--end of the method


    }

}


/*

M=...
    [
    0.3559                 0.1298                0.0143;
    0.3559*0.4249          0.1298*1.6313         0.0143*2.5106;
    0.3559*(0.4249^2)      0.1298*(1.6313^2)     0.0143*(2.5106^2);
    ]
 * 
 inv(M) =

    4.5734   -4.6252    1.1167
   -7.7475   21.3196   -7.2627
   26.4299  -78.4044   38.1307

*/