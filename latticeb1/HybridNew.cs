using System;
using System.Collections.Generic;
using System.Linq;
using System.Windows.Forms;

// using BGK.Models;
using Latticeb.Models;

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

            //--hybrid Left and Right BGK segments
            float[,] fleft = new float[xsteps, nvlc];
            float[,] fright = new float[xsteps, nvlc];
            const float Kn = 1.0f;
            const float vw = 0.02f; // distance between velocities
            float totmass = 0f; // total mass

            // initial boundary conditions
            float u1 = 0.0f;
            float T1 = 1.0f;
            float rho1 = 1.0f; // нормировка ? / BGK_1d._n_eq_check_dens(T1, u1, vw, nvlc);

            float u2 = 0.0f;
            float T2 = 1.0f;
            float rho2 = 1.0f; //? / BGK_1d._n_eq_check_dens(T2, u2, vw, nvlc);

            //--BGK boundaries, time = 0;
            for (int k = 0; k < nvlc; ++k)
            {
                //--hybrid
                if (k <= (nvlc - 1) / 2) // условие может быть не актуальным
                {
                    fright[xsteps - 1, k] = BGK_1d._n_eq(T2, rho2, u2, vw, nvlc, k);
                    // f[0, k] = 0; не актуально
                }
                else
                {
                    // f[xsteps - 1, k] = 0; не актуально
                    fleft[0, k] = BGK_1d._n_eq(T1, rho1, u1, vw, nvlc, k);
                }
            }

            const int nlattice = 161; // space nodes quantity for lattice
            float[,] flattice = new float[nlattice, 6];

            //--hybrid declaration
            BGK_1d bgk_left, bgk_right;
            D1Q6E lattice;

            // stabilising values for lattice
            for (int i = 1; i < nlattice - 1; ++i)
            {
                flattice[i, 0] = D1Q6E._n_eq_k(1.0f, u1, 0);
                flattice[i, 1] = D1Q6E._n_eq_k(1.0f, u1, 1);
                flattice[i, 2] = D1Q6E._n_eq_k(1.0f, u1, 2);
                flattice[i, 3] = D1Q6E._n_eq_k(1.0f, u1, 3);
                flattice[i, 4] = D1Q6E._n_eq_k(1.0f, u1, 4);
                flattice[i, 5] = D1Q6E._n_eq_k(1.0f, u1, 5);
            }

            BGK_1d testBGK = new BGK_1d(fleft, Kn, vw, physlength, 1);
            int timeindex = (int)(1.0f / testBGK.dtime());
            // I recommend to check this value

            for (int time = 0; time < tsteps; ++time)
            {
                // константа 2000 взята наобум, т.к. заранее это значение неизвестно.
                bgk_left = new BGK_1d(fleft, Kn, vw, physlength, timeindex + 1);
                bgk_right = new BGK_1d(fright, Kn, vw, physlength, timeindex + 1);

                //--hybrid, solve BGK for t = 1.0
                bgk_left.Solve(timeindex);
                bgk_right.Solve(timeindex);
                for (int i = 0; i < xsteps; ++i)
                {
                    for (int j = 0; j < nvlc; ++j)
                    {
                        fleft[i, j] = bgk_left.F[timeindex, i, j];
                        fright[i, j] = bgk_right.F[timeindex, i, j];
                    }
                }

                float m0, m1, m2;
                m0 = m1 = m2 = 0;
                // этот код отвечает за подсчет моментов при склейке, его стоит проверить
                for (int i = nvlc / 2; i < nvlc; ++i)
                {
                    m0 += fleft[xsteps - 1, i];
                    m1 += fleft[xsteps - 1, i] * (vw * (i - (nvlc - 1.0f) / 2.0f));
                    m2 += fleft[xsteps - 1, i] * (vw * (i - (nvlc - 1.0f) / 2.0f))
                                      * (vw * (i - (nvlc - 1.0f) / 2.0f));
                }
                // m0 *= 1.0f; m1 *= vw; m2 *= vw; проверить

                //--setting LBE condition (flattice[,] = F(flattice, fleft, fright) )
                flattice[0, 3] = D1Q6E.distr0(m0, m1, m2);
                flattice[0, 4] = D1Q6E.distr1(m0, m1, m2);
                flattice[0, 5] = D1Q6E.distr2(m0, m1, m2);

                m0 = m1 = m2 = 0;
                for (int i = 0; i < nvlc / 2; ++i) // проверить !
                {
                    m0 += fright[0, i];
                    m1 -= fright[0, i] * (vw * (i - (nvlc - 1.0f) / 2.0f));
                    m2 += fright[0, i] * (vw * (i - (nvlc - 1.0f) / 2.0f))
                                               * (vw * (i - (nvlc - 1.0f) / 2.0f));
                }
                // m0 *= 1.0; m1 *= vw; m2 *= vw;

                flattice[nlattice - 1, 2] = D1Q6E.distr0(m0, m1, m2);
                flattice[nlattice - 1, 1] = D1Q6E.distr1(m0, m1, m2);
                flattice[nlattice - 1, 0] = D1Q6E.distr2(m0, m1, m2);

                lattice = new D1Q6E(flattice, 1.0f);
                flattice = lattice.Solve(); // считает единицу во времени


                //--левая точка сращивания
                m0 = D1Q6E.moment0(flattice[0, 2], flattice[0, 1], flattice[0, 0]);
                m1 = -D1Q6E.moment1(flattice[0, 2], flattice[0, 1], flattice[0, 0]);
                m2 = D1Q6E.moment2(flattice[0, 2], flattice[0, 1], flattice[0, 0]);
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
                m0 = D1Q6E.moment0(flattice[nlattice - 1, 3], flattice[nlattice - 1, 4], flattice[nlattice - 1, 5]);
                m1 = D1Q6E.moment1(flattice[nlattice - 1, 3], flattice[nlattice - 1, 4], flattice[nlattice - 1, 5]);
                m2 = D1Q6E.moment2(flattice[nlattice - 1, 3], flattice[nlattice - 1, 4], flattice[nlattice - 1, 5]);

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
            bgk_left = new BGK_1d(fleft, Kn, vw, physlength, 1);
            bgk_right = new BGK_1d(fright, Kn, vw, physlength, 1);
            float[,] P = bgk_left.P();
            float[,] U = bgk_left.U();
            float[,] T = bgk_left.Temperatures();

            totmass = 0; //--diagnostics of the left BGK area
            for (int i = 0; i < xsteps; i++)
            {
                Console.WriteLine(i / (xsteps - 1.0f) + " "
                    + P[0, i] + " " 
                    + U[0, i] + " " 
                    + T[0, i] + " " 
                    + P[0, i] * U[0, i]);
                totmass = totmass + P[0, i];
            }

            //--diagnostics lattice
            for (int i = 0; i < nlattice; i++)
            {
                Console.WriteLine(flattice[i, 0] + " " + flattice[i, 1] + " " + flattice[i, 2] +
                            " " + flattice[i, 3] + " " + flattice[i, 4] + " " + flattice[i, 5] + " x" + i * 0.5);
            }

            Console.WriteLine("mass=" + totmass);
            Console.ReadLine();

        }
    }
}
