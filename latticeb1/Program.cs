using System;
using System.Collections.Generic;
using System.Linq;
using System.Windows.Forms;

namespace latticeb
{
    static class Program
    {
        [STAThread]
        static void Main()
        {
            const int tsteps = 2000;  // time steps
            const float physlength = 10.0f; // phys length
            const int xsteps = 20; // space nodes quantity for BGK
            const int nvlc = 60; // velocity nodes quantity for BGK

            //--BGK segment
            float[,] f = new float[xsteps, nvlc];
            const float Kn = 1.0f;
            const float vw = 0.2f; // distance between velocities
            float totmass = 0f; // total mass

            // initial boundary conditions
            float u1 = 0.0f;
            float T1 = 1.0f;
            float rho1 = 1.0f;

            float u2 = 0.0f;
            float T2 = 1.0f;
            float rho2 = 1.0f;

            //--BGK boundaries, time = 0;
            for (int k = 0; k < nvlc; ++k)
            {
                //--hybrid
                if (k <= (nvlc - 1) / 2)
                {
                    f[xsteps - 1, k] = BGK_1d._n_eq(T2, rho2, u2, vw, nvlc, k);
                    f[0, k] = 0;
                }
                else { 
                    f[xsteps - 1, k] = 0;
                    f[0, k] = BGK_1d._n_eq(T1, rho1, u1, vw, nvlc, k);
                }
            }


            BGK_1d bgk = new BGK_1d(f, Kn, vw, physlength, tsteps);
            bgk.Solve(tsteps);

            //--results
            float[,] P = bgk.P();
           // float[,] Ps = bgk.Ps();
            float[,] U = bgk.U();
            float[,] T = bgk.Temperatures();

            totmass = 0; //--diagnostics of the left BGK area
            for (int i = 0; i < xsteps; i++)
            {
                Console.WriteLine(i / (xsteps - 1.0f) + " " + P[tsteps - 1, i] + " " + U[tsteps - 1, i] + " " + T[tsteps - 1, i] + " " + P[tsteps - 1, i] * U[tsteps-1,i]);
                totmass = totmass + P[tsteps-1,i];
            }
           
            
            Console.ReadLine();
        }
    }
}
