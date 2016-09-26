using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;


namespace latticeb
{
    /// <summary>
    /// BGK model solver class
    /// </summary>
    class BGK_1d
    {
        //---distribution function;
        //--first dimension - time; second - x coordinate; third - velocity
        public float[, ,] F;
        //--sound velocity squared
        // private float cs = 1f;

        private float vw; // velocity weight
        //--invtau
        private float invtau;
        private int tsteps; // number of time steps, technically=2 for 1 step solver
        private int nspace; // number of cells
        private int nvel; // number of velocities 
        private float dx; // space step;
        private float dt; // time step;


        private float[,] Faux;

        public float dtime()
        {
            return dt;
        }

        //--constructor: define initial conditions
        public BGK_1d(float[,] Finit, float Kn, float vlcweight, float physlength, int tm)
        {

            invtau = 1.0f / Kn; // Knudsen number
            vw = vlcweight; // dist. between velocities
            tsteps = tm;
            nspace = Finit.GetLength(0);
            nvel = Finit.GetLength(1);

            dx = physlength / (nspace - 1.0f); // space step
            dt = 0.85f * dx / (vw * (nvel - 1.0f) * 0.5f); // time step (CFL condition)


            F = new float[tsteps, nspace, nvel];
            Faux = new float[nspace, nvel];

            //---set the initial distribution
            for (int j = 0; j < Finit.GetLength(0); j++)
            {
                for (int k = 0; k < nvel; k++)
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
                result_dens += F[t, j, k] * vw;
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
                    result_vlc += F[t, j, k] * (k - (nvel - 1.0f) / 2.0f) * vw * vw;
                }
                return result_vlc / rho;
            }
        }

        private float _temp(int t, int j, float rho, float u)
        { // equilibrium temperature evaluation
            float result_temp = 0;
            for (int k = 0; k < nvel; k++)
            {
                result_temp += ((k - (nvel - 1.0f) / 2.0f) * vw - u) *
                               ((k - (nvel - 1.0f) / 2.0f) * vw - u) *
                               F[t, j, k] * vw;
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



        //-- time evolution 

        public void Solve(float time)
        {
            float c;

            //-- start the time loop
            for (int i = 0; i < time - 1; ++i)
            {
                //--step 0: set the left and right boundary vals
                for (int k = 0; k < nvel; ++k)
                {
                    if (k < nvel / 2.0) { F[i + 1, nspace - 1, k] = F[i, nspace - 1, k]; }
                    if (k >= nvel / 2.0) { F[i + 1, 0, k] = F[i, 0, k]; }
                }

                //---step 1:  advection

                for (int j = 0; j < nspace; j++)
                {
                    for (int k = 0; k < nvel; k++)
                    {
                        c = (vw * (k - (nvel - 1) / 2.0f)) * dt / dx;

                        if (k < nvel / 2.0 & j < nspace - 1) { Faux[j, k] = (1 + c) * F[i, j, k] - c * F[i, j + 1, k]; }
                        if (k >= nvel / 2.0 & j > 0) { Faux[j, k] = (1 - c) * F[i, j, k] + c * F[i, j - 1, k]; }
                    }

                }

                //---step2: collision
                for (int j = 0; j < nspace; j++)
                {
                    float rho = _dens(i, j);
                    float u = _vlc(i, j, rho);
                    float T = _temp(i, j, rho, u);
                    float rhs = 0;

                    for (int k = 0; k < nvel; ++k)
                    {
                        rhs = invtau * (_n_eq(T, rho, u, vw, nvel, k) - F[i, j, k]) * dt;
                        //Console.WriteLine("i=" + i + " j=" + j + " k=" + k + " rhs=" + rhs);
                        if (j == 0 & k < nvel / 2.0) { F[i + 1, j, k] = rhs + Faux[j, k]; }
                        if (j == nspace - 1 & k >= nvel / 2.0) { F[i + 1, j, k] = rhs + Faux[j, k]; }
                        if (j > 0 & j < nspace - 1) { F[i + 1, j, k] = rhs + Faux[j, k]; }

                    }

                }

            }

        }

        //--evaluate density (x)
        public float[,] P()
        {
            float[,] dens = new float[tsteps, nspace];

            for (int t = 0; t < tsteps; t++)
            {
                for (int j = 0; j < nspace; j++) dens[t, j] = _dens(t, j);
            }

            return dens;
        }//--end of the method



        //--evaluate bulk velocity (x)
        public float[,] U()
        {
            float[,] u = new float[tsteps, nspace];
            for (int t = 0; t < tsteps; t++)
            {
                for (int j = 0; j < nspace; j++) u[t, j] = _vlc(t, j, _dens(t, j));
            }
            return u;
        }//--end of the method


        //--evaluate temperature (x)
        public float[,] Temperatures()
        {
            float[,] T = new float[tsteps, nspace];
            for (int t = 0; t < tsteps; t++)
            {
                for (int j = 0; j < nspace; j++) T[t, j] = _temp(t, j, _dens(t, j), _vlc(t, j, _dens(t, j)));
            }
            return T;
        }//--end of the method

    }

}
