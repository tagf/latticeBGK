using System;
using System.Collections.Generic;
using System.Linq;
using System.Windows.Forms;

using BGK.Models;
using Latticeb.Models;


namespace Latticeb
{
    static class Program
    {
        [STAThread]
        static void Main()
        {
            const int tsteps = 1;  // time steps
            const float physlength = 100.0f; // phys length
            const int xsteps = 101; // space nodes quantity for BGK
            const int nvlc = 100; // velocity nodes quantity for BGK

            //--BGK segment
            float[,] f = new float[xsteps, nvlc];
            const float Kn = 1.0f;
            const float vw = 0.05f; // distance between velocities
            float totmass = 0f; // total mass

            // initial boundary conditions
            float u1 = 0.8f;
            float T1 = 1.0f;
            float rho1 = 1.0f / BGK_1d._n_eq_check_dens(T1, u1, vw, nvlc);

            float u2 = 0.8f;
            float T2 = 1.0f;
            float rho2 = 1.0f / BGK_1d._n_eq_check_dens(T2, u2, vw, nvlc);

            //--BGK boundaries, time = 0;
            for (int k = 0; k < nvlc; ++k)
            {
                //--hybrid
                f[0, k] = BGK_1d._n_eq(T1, rho1, u1, vw, nvlc, k);
                f[xsteps - 1, k] = BGK_1d._n_eq(T2, rho2, u2, vw, nvlc, k);
            }

            BGK_1d bgk = new BGK_1d(f, Kn, vw, physlength);
            f = bgk.Solve(1000.0f);

            //--results
            float[] P = bgk.P();
            float[] U = bgk.U();
            float[] T = bgk.Temperatures();

            totmass = 0; //--diagnostics of the left BGK area
            for (int i = 0; i < xsteps; i++)
            {
                Console.WriteLine(i / (xsteps - 1.0f) + " " + P[i] + " " + U[i] + " " + T[i] + " " + P[i] * U[i]);
                totmass = totmass + P[i];
            }


            Console.WriteLine("mass=" + totmass);
            Console.ReadLine();

        }
    }
}
