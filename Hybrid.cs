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
            const float physlength = 10.0f; // phys length
            const int xsteps = 40; // space nodes quantity for BGK
            const int nvlc = 60; // velocity nodes quantity for BGK

            //--hybrid Left and Right BGK segments
            float[,] fleft = new float[xsteps, nvlc];
            float[,] fright = new float[xsteps, nvlc];
            const float Kn = 1.0f;
            const float vw = 0.024f; // distance between velocities
            float totmass = 0f; // total mass

            // initial boundary conditions
            float u1 = 0.08f;
            float T1 = 0.64f;
            float rho1 = 1.0f / BGK_1d._n_eq_check_dens(T1, u1, vw, nvlc);

            float u2 = -0.055f;
            float T2 = 0.64f;
            float rho2 = 1.0f / BGK_1d._n_eq_check_dens(T2, u2, vw, nvlc);

            //--BGK boundaries, time = 0;
            for (int k = 0; k < nvlc; ++k)
            {
                //--hybrid
                fleft[0, k] = BGK_1d._n_eq(T1, rho1, u1, vw, nvlc, k);
                fright[xsteps - 1, k] = BGK_1d._n_eq(T2, rho2, u2, vw, nvlc, k);
            }


            const int nlattice = 161; // space nodes quantity for lattice
            float[,] flattice = new float[nlattice, 6];

            //--hybrid declaration
            BGK_1d bgk_left, bgk_right;
            D1Q6 lattice;

            float m0 = 0, m1 = 0, m2 = 0;

            const float w0 = 0.3234f;
            const float w1 = 0.1463f;
            const float w2 = 0.0303f;

            /*
            {{0.3234 * 1,0.1463 * 1,0.0303 * 1},
             {0.3234 * 0.5,0.1463 * 1,0.0303 * 3},
             {0.3234 * 0.25,0.1463 * 1,0.0303 * 9}}^-1 = 
             {{7.42115, -9.89487, 2.47372},
             {-10.2529, 23.9234, -6.83527},
             {3.30033, -9.90099, 6.60066}}
              wolframalpha.com */

            const float a11 = 7.42115f;
            const float a12 = -9.89487f;
            const float a13 = 2.47372f;
            const float a21 = -10.2529f;
            const float a22 = 23.9234f;
            const float a23 = -6.83527f;
            const float a31 = 3.30033f;
            const float a32 = -9.90099f;
            const float a33 = 6.60066f;

            for (int i = 1; i < nlattice - 1; ++i)
            {
                flattice[i, 0] = D1Q6._n_eq_k(1.0f, u1, 0);
                flattice[i, 1] = D1Q6._n_eq_k(1.0f, u1, 1);
                flattice[i, 2] = D1Q6._n_eq_k(1.0f, u1, 2);
                flattice[i, 3] = D1Q6._n_eq_k(1.0f, u1, 3);
                flattice[i, 4] = D1Q6._n_eq_k(1.0f, u1, 4);
                flattice[i, 5] = D1Q6._n_eq_k(1.0f, u1, 5);
            }

            for (int time = 0; time < tsteps; ++time)
            {
                bgk_left = new BGK_1d(fleft, Kn, vw, physlength);
                bgk_right = new BGK_1d(fright, Kn, vw, physlength);

                //--hybrid, solve BGK for t = 1.0
                fleft = bgk_left.Solve(1.0f);
                fright = bgk_right.Solve(1.0f);

                m0 = m1 = m2 = 0;
                for (int i = nvlc / 2; i < nvlc; ++i)
                {
                    m0 += fleft[xsteps - 1, i];
                    m1 += fleft[xsteps - 1, i] * (vw * (i - (nvlc - 1.0f) / 2.0f));
                    m2 += fleft[xsteps - 1, i] * (vw * (i - (nvlc - 1.0f) / 2.0f))
                                      * (vw * (i - (nvlc - 1.0f) / 2.0f));
                }
                // m0 *= 1.0f; m1 *= vw; m2 *= vw;

                //--TODO: set LBE condition (flattice[,] = F(flattice, fleft, fright) )
                flattice[0, 3] = a11 * m0 + a12 * m1 + a13 * m2;
                flattice[0, 4] = a21 * m0 + a22 * m1 + a23 * m2;
                flattice[0, 5] = a31 * m0 + a32 * m1 + a33 * m2;

                m0 = m1 = m2 = 0;
                for (int i = 0; i < nvlc / 2; ++i)
                {
                    m0 += fright[0, i];
                    m1 -= fright[0, i] * (vw * (i - (nvlc - 1.0f) / 2.0f));
                    m2 += fright[0, i] * (vw * (i - (nvlc - 1.0f) / 2.0f))
                                               * (vw * (i - (nvlc - 1.0f) / 2.0f));
                }
                // m0 *= 1.0; m1 *= vw; m2 *= vw;

                //--TODO: set LBE condition (flattice[,] = F(flattice, fleft, fright) )
                flattice[nlattice - 1, 2] = a11 * m0 + a12 * m1 + a13 * m2;
                flattice[nlattice - 1, 1] = a21 * m0 + a22 * m1 + a23 * m2;
                flattice[nlattice - 1, 0] = a31 * m0 + a32 * m1 + a33 * m2;

                lattice = new D1Q6(flattice, 1.0f);
                flattice = lattice.Solve();

                //--TODO: set new BGK internal boundary data fleft, fright = g(fleft(right), flattice)

                //--левая точка сращивания
                m0 = w0 * 1.0f * flattice[0, 2] + w1 * flattice[0, 1] + w2 * 1.0f * flattice[0, 0];
                m1 = -w0 * 0.5f * flattice[0, 2] - w1 * flattice[0, 1] - w2 * 3.0f * flattice[0, 0];
                m2 = w0 * 0.25f * flattice[0, 2] + w1 * flattice[0, 1] + w2 * 9.0f * flattice[0, 0];
                for (int i = 0; i < nvlc / 2; ++i)
                {
                    fleft[xsteps - 1, i] = vw * 1.0f / (float)Math.Sqrt(2.0 * Math.PI * 1.0f)
                        * (float)Math.Exp(-
                 ((i - (nvlc - 1.0) / 2.0) * vw) *
                 ((i - (nvlc - 1.0) / 2.0) * vw) / (2.0))
                 * (float)(m0 + m1 * (i - (nvlc - 1.0) / 2.0) * vw +
                 (m1 - m0) *
                 (((i - (nvlc - 1.0) / 2.0) * vw) *
                 ((i - (nvlc - 1.0) / 2.0) * vw) - 1.0)
                 );
                }

                //--правая точка сращивания
                m0 = w0 * 1.0f * flattice[nlattice - 1, 3] + w1 * flattice[nlattice - 1, 4] + w2 * 1.0f * flattice[nlattice - 1, 5];
                m1 = w0 * 0.5f * flattice[nlattice - 1, 3] + w1 * flattice[nlattice - 1, 4] + w2 * 3.0f * flattice[nlattice - 1, 5];
                m2 = w0 * 0.25f * flattice[nlattice - 1, 3] + w1 * flattice[nlattice - 1, 4] + w2 * 9.0f * flattice[nlattice - 1, 5];
                for (int i = nvlc / 2; i < nvlc; ++i)
                {
                    fright[0, i] = vw * 1.0f / (float)Math.Sqrt(2.0 * Math.PI * 1.0f)
                        * (float)Math.Exp(-
                 ((i - (nvlc - 1.0) / 2.0) * vw) *
                 ((i - (nvlc - 1.0) / 2.0) * vw) / (2.0))
                 * (float)(m0 + m1 * (i - (nvlc - 1.0) / 2.0) * vw +
                 (m1 - m0) *
                 (((i - (nvlc - 1.0) / 2.0) * vw) *
                 ((i - (nvlc - 1.0) / 2.0) * vw) - 1.0)
                 );
                }
            }

            //--results
            bgk_left = new BGK_1d(fleft, Kn, vw, physlength);
            bgk_right = new BGK_1d(fright, Kn, vw, physlength);
            bgk_left.Solve_one_step();
            bgk_right.Solve_one_step();
            float[] P = bgk_left.P();
            float[] U = bgk_left.U();
            float[] T = bgk_left.Temperatures();

            totmass = 0; //--diagnostics of the left BGK area
            for (int i = 0; i < xsteps; i++)
            {
                Console.WriteLine(i / (xsteps - 1.0f) + " " + P[i] + " " + U[i] + " " + T[i] + " " + P[i] * U[i]);
                totmass = totmass + P[i];
            }

            //--diagnostics lattice
            for (int i = 0; i < nlattice; i++)
            {
                Console.WriteLine(flattice[i, 0] + " " + flattice[i, 1] + " " + flattice[i, 2] +
                            " " + flattice[i, 3] + " " + flattice[i, 4] + " " + flattice[i, 5] + " x" + i * 0.5);
            }

            Console.WriteLine("mass=" + totmass);
            Console.ReadLine();

            //   Application.EnableVisualStyles();
            //   Application.SetCompatibleTextRenderingDefault(false);
            //   Application.Run(new Form1());
        }
    }
}
