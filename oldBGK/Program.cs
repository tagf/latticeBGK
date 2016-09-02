using System;
using System.Collections.Generic;
using System.Linq;
using System.Windows.Forms;
using Latticeb.Models;
using BGK.Models;


namespace Latticeb
{
    static class Program
    {
        [STAThread]
        static void Main()
        {
            const int tsteps = 100000, xsteps=51; // time steps, space nodes quantity
            const int nvlc = 40; // velocity nodes quantity

            float[,] finit=new float[xsteps, nvlc];
            
            const float Kn = 1.0f;
            const float vw = 0.2f; // distance between velocities

            // initial boundary conditions
            float u1 = 0.5f;
            float T1 = 1.0f;
            float rho1 = 1.0f / BGK_1d._n_eq_check_dens(T1, u1, vw, nvlc);

            float u2 = -0.3f;
            float T2 = 1.0f;
            float rho2 = 1.0f / BGK_1d._n_eq_check_dens(T2, u2, vw, nvlc);

            //--initial data and boundaries, time = 0;
            for(int k = 0; k < nvlc ; ++k)
            {
                finit[0, k] = BGK_1d._n_eq(T1, rho1, u1, vw, nvlc, k);
                finit[5, k] = BGK_1d._n_eq(0.8f, 1.2f / BGK_1d._n_eq_check_dens(0.8f, 0.3f, vw, nvlc), 0.3f, vw, nvlc, k); // test
                finit[xsteps - 1, k] = BGK_1d._n_eq(T2, rho2, u2, vw, nvlc, k);
            }

            //-- solving !
            //-- in integrated version with lattice model only one_step solver would be actual
            BGK_1d flow = new BGK_1d(finit, Kn, vw); // solver
            float[,] Fnew = flow.Solve_one_step(); // time layer


            for (int time = 0; time < tsteps / 1000; ++time)
            {
                for (int tmp = 0; tmp < 1000; ++tmp) 
                {
                  flow = new BGK_1d(Fnew, Kn, vw);
                  Fnew = flow.Solve_one_step();
                }
                Console.WriteLine(time);
            }

            //--results

            float[] P = flow.P();
            float[] U = flow.U();
            float[] T = flow.Temperatures();
            
            float totmass = 0f;
            for (int i = 0; i < xsteps; i++)
            {
                Console.WriteLine(i / (xsteps - 1.0f) + " " + P[i] + " " + U[i] + " " + T[i] + " " + P[i]*U[i]);
                totmass = totmass + P[i];
            }

            Console.WriteLine("mass=" + totmass);
            Console.ReadLine();

         //   Application.EnableVisualStyles();
         //   Application.SetCompatibleTextRenderingDefault(false);
         //   Application.Run(new Form1());
        }
    }
}
