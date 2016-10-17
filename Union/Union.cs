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
            // General
            const int tsteps = 10000;  // time steps for initial run of BGK / lattice
            const float physlength = 1.0f; // phys length
            const float Kn = 1.0f; // Knudsen number

            //--BGK segment
            const int Init_BGK_Steps = 51; // space nodes quantity for initial BGK run (50 + 1)
 
            const float vw = 0.2f; // distance between velocities
            const int nvlc = 60; // velocity nodes quantity for BGK
            float totmass = 0f; // total mass

            // boundary conditions for problem
            float u1 = 0.5f;
            float T1 = 0.7f;
            float rho1 = 1.0f;

            float u2 = -0.3f;
            float T2 = 0.7f;
            float rho2 = 1.0f;

            //-- half moments of BGK
            float lm0, lm1, lm2; 
            float rm0, rm1, rm2;

            // HYBRIDISATION PARAMETERS / VARIABLES
            //--hybrid Left and Right BGK segments
            const int Xsteps = Init_BGK_Steps / 5 + 1; // space nodes quantity for BGK left & right (0.2 * 50 + 1)
            const int Nlattice = 121; // space nodes quantity for lattice
            // data transmission variables
            float[,] fleft = new float[Xsteps, nvlc];
            float[,] fright = new float[Xsteps, nvlc];
            float[,] flattice = new float[Nlattice, 6];
            
            //-- initial 0.2 and 0.8 half moments computation
            {// initial full domain BGK run, incapsulation in scope
                float[,] f = new float[Init_BGK_Steps, nvlc];

                //--BGK boundaries, time = 0;
                for (int k = 0; k < nvlc; ++k)
                {
                    //--non-hybrid
                    if (k <= (nvlc - 1) / 2)
                    {
                        f[Init_BGK_Steps - 1, k] = BGK_1d._n_eq(T2, rho2, u2, vw, nvlc, k);
                        f[0, k] = 0;
                    }
                    else
                    {
                        f[Init_BGK_Steps - 1, k] = 0;
                        f[0, k] = BGK_1d._n_eq(T1, rho1, u1, vw, nvlc, k);
                    }
                }

                BGK_1d bgk = new BGK_1d(f, Kn, vw, physlength, tsteps);
                bgk.Solve(tsteps);

                lm0 = lm1 = lm2 = 0;
                int tstep = 10000 - 1;
                int xstep = 10;

                // этот код отвечает за подсчет моментов при склейке
                for (int i = nvlc / 2; i < nvlc; ++i)
                {
                    lm0 += bgk.F[tstep, xstep, i];
                    lm1 += bgk.F[tstep, xstep, i] * (vw * (i - (nvlc - 1.0f) / 2.0f));
                    lm2 += bgk.F[tstep, xstep, i] * (vw * (i - (nvlc - 1.0f) / 2.0f))
                                      * (vw * (i - (nvlc - 1.0f) / 2.0f));
                }
                lm0 *= vw; lm1 *= vw; lm2 *= vw;

                rm0 = rm1 = rm2 = 0;
                xstep = 40;

                // этот код отвечает за подсчет моментов при склейке
                for (int i = 0; i < nvlc / 2; ++i)
                {
                    rm0 += bgk.F[tstep, xstep, i];
                    rm1 += bgk.F[tstep, xstep, i] * (vw * (i - (nvlc - 1.0f) / 2.0f));
                    rm2 += bgk.F[tstep, xstep, i] * (vw * (i - (nvlc - 1.0f) / 2.0f))
                                      * (vw * (i - (nvlc - 1.0f) / 2.0f));
                }
                rm0 *= vw; rm1 *= vw; rm2 *= vw;

                // provide initial values for left and right BGK domains for iterations
                for (int i = 0; i < Xsteps; ++i)
                {
                    for (int j = 0; j < nvlc; ++j)
                    {
                        fleft[i, j] = f[i, j];
                        fright[i, j] = f[Init_BGK_Steps - Xsteps + i, j];
                    }
                }
            }

            // what we have computed now are initialising moments for lattice
            Console.WriteLine("m0=" + lm0 + " m1=" + lm1 + " m2=" + lm2);
            Console.WriteLine("m0=" + rm0 + " m1=" + rm1 + " m2=" + rm2);

            /* initial boundary data for lattice block */
            {
                float m0, m1, m2;

                // m0 = 0.7572213f; m1 = 0.6452931f; m2 = 0.8266699f; should be
                m0 = lm0;
                m1 = lm1;
                m2 = lm2;

                //--setting LBE condition (flattice[,] = F(flattice, fleft, fright) )
                flattice[0, 3] = D1Q6.distr0(m0, m1, m2);
                flattice[0, 4] = D1Q6.distr1(m0, m1, m2);
                flattice[0, 5] = D1Q6.distr2(m0, m1, m2);

                // m0 = 0.6613145f; m1 = 0.5069057f; m2 = 0.5982312f; should be
                m0 = rm0;
                m1 = - /* minus required ! */ rm1;
                m2 = rm2;

                flattice[Nlattice - 1, 2] = D1Q6.distr0(m0, m1, m2);
                flattice[Nlattice - 1, 1] = D1Q6.distr1(m0, m1, m2);
                flattice[Nlattice - 1, 0] = D1Q6.distr2(m0, m1, m2);
            }

            // Lattice domain initial computation
            for (int time = 0; time < tsteps; ++time)
            {
                D1Q6 lattice = new D1Q6(flattice, 100.0f /* WARNING, this number should depend on KN and other params! */);
                flattice = lattice.Solve(); // считает единицу во времени
            }

            // HYBRIDISATION ITERATIONS
            // check for consistency
            const int NumIterations = 1;
            for (int i = 0; i < NumIterations; ++i)
            {

                // 1. solving prepared BGK subdomains
                BGK_1d lbgk = new BGK_1d(fleft, Kn, vw, physlength, tsteps);
                BGK_1d rbgk = new BGK_1d(fright, Kn, vw, physlength, tsteps);
                lbgk.Solve(tsteps);
                rbgk.Solve(tsteps);

                // extraction
                for (int k = 0; k < Xsteps; ++k)
                {
                    for (int j = 0; j < nvlc; ++j)
                    {
                        fleft[k, j] = lbgk.F[tsteps - 1, k, j];
                        fright[k, j] = rbgk.F[tsteps - 1, k, j];
                    }
                }

                // 2. computation of half-moments from BGKs 

                lm0 = lm1 = lm2 = 0;
                int xstep;
                int tstep = tsteps - 1; 
                xstep = Xsteps - 1;
                // этот код отвечает за подсчет моментов при склейке
                for (int j = nvlc / 2; i < nvlc; ++i)
                {
                    lm0 += lbgk.F[tstep, xstep, i];
                    lm1 += lbgk.F[tstep, xstep, i] * (vw * (i - (nvlc - 1.0f) / 2.0f));
                    lm2 += lbgk.F[tstep, xstep, i] * (vw * (i - (nvlc - 1.0f) / 2.0f))
                                      * (vw * (i - (nvlc - 1.0f) / 2.0f));
                }
                lm0 *= vw; lm1 *= vw; lm2 *= vw;

                //--setting LBE condition (flattice[,] = F(flattice, fleft, fright) )
                //Console.WriteLine("m0=" + lm0 + " m1=" + lm1 + " m2=" + lm2);

                rm0 = rm1 = rm2 = 0;
                xstep = 0;

                // этот код отвечает за подсчет моментов при склейке
                for (int j = 0; i < nvlc / 2; ++i)
                {
                    rm0 += rbgk.F[tstep, xstep, i];
                    rm1 += rbgk.F[tstep, xstep, i] * (vw * (i - (nvlc - 1.0f) / 2.0f));
                    rm2 += rbgk.F[tstep, xstep, i] * (vw * (i - (nvlc - 1.0f) / 2.0f))
                                      * (vw * (i - (nvlc - 1.0f) / 2.0f));
                }
                rm0 *= vw; rm1 *= vw; rm2 *= vw;

                /* boundary data for lattice block */
                {
                    float m0, m1, m2;

                    // m0 = 0.7572213f; m1 = 0.6452931f; m2 = 0.8266699f; should be
                    m0 = lm0;
                    m1 = lm1;
                    m2 = lm2;

                    //--setting LBE condition (flattice[,] = F(flattice, fleft, fright) )
                    flattice[0, 3] = D1Q6.distr0(m0, m1, m2);
                    flattice[0, 4] = D1Q6.distr1(m0, m1, m2);
                    flattice[0, 5] = D1Q6.distr2(m0, m1, m2);

                    // m0 = 0.6613145f; m1 = 0.5069057f; m2 = 0.5982312f; should be
                    m0 = rm0;
                    m1 = - /* minus required ! */ rm1;
                    m2 = rm2;

                    flattice[Nlattice - 1, 2] = D1Q6.distr0(m0, m1, m2);
                    flattice[Nlattice - 1, 1] = D1Q6.distr1(m0, m1, m2);
                    flattice[Nlattice - 1, 0] = D1Q6.distr2(m0, m1, m2);
                }

                // 3. lattice step
                D1Q6 lattice = new D1Q6(flattice, 100.0f);
                flattice = lattice.Solve();

                // 4. back from lattice to updating BGK subdomain border data
                // to the left
                lm0 = D1Q6.moment0(flattice[0, 2], flattice[0, 1], flattice[0, 0]);
                lm1 = D1Q6.moment1(flattice[0, 2], flattice[0, 1], flattice[0, 0]);
                lm2 = D1Q6.moment2(flattice[0, 2], flattice[0, 1], flattice[0, 0]);
                // to the right
                rm0 = D1Q6.moment0(flattice[Nlattice - 1, 3], flattice[Nlattice - 1, 4], flattice[Nlattice - 1, 5]);
                rm1 = D1Q6.moment1(flattice[Nlattice - 1, 3], flattice[Nlattice - 1, 4], flattice[Nlattice - 1, 5]);
                rm2 = D1Q6.moment2(flattice[Nlattice - 1, 3], flattice[Nlattice - 1, 4], flattice[Nlattice - 1, 5]);

                // TODO: set BGK half-distributions using m0-m2 for left and right subdomains respectively

            }

            // RESULTS
            //--check whether obtaining of results is correct 
            BGK_1d bgk_left = new BGK_1d(fleft, Kn, vw, physlength, 1);
            BGK_1d bgk_right = new BGK_1d(fright, Kn, vw, physlength, 1);
            float[,] P = bgk_left.P();
            float[,] U = bgk_left.U();
            float[,] T = bgk_left.Temperatures();

            totmass = 0; //--diagnostics of the left BGK area
            Console.WriteLine("FINAL DIAGNOSTICS, lbgk");
            for (int i = 0; i < Xsteps; i++)
            {
                Console.WriteLine(i / (Xsteps - 1.0f) + " "
                    + P[0, i] + " " 
                    + U[0, i] + " " 
                    + T[0, i] + " " 
                    + P[0, i] * U[0, i]);
                totmass = totmass + P[0, i];
            }

            //--diagnostics lattice
            D1Q6 Lattice = new D1Q6(flattice, 100f /* 100f = f(Kn, ...) !! */);
            Lattice.Solve();
            float[,] Temp = Lattice.T();

            Console.WriteLine("FINAL DIAGNOSTICS, lattice");
            for (int i = 0; i < Nlattice; i++)
            {
                float dens = flattice[i, 0] + flattice[i, 1] +
                            + flattice[i, 2] + flattice[i, 3]
                            + flattice[i, 4] + flattice[i, 5];
                float velc = flattice[i, 4] - flattice[i, 1] +
                            +0.5f * (flattice[i, 3] - flattice[i, 2])
                            + 3.0f * (flattice[i, 5] - flattice[i, 0]);
                velc /= dens;
                Console.WriteLine(flattice[i, 0] + " " + flattice[i, 1] + " " + flattice[i, 2] +
                            " " + flattice[i, 3] + " " + flattice[i, 4] + " " + flattice[i, 5] + " x " + (0.2f + i * 0.6f / 120f)
                            + " d= " + dens + " v= " + velc + " T= " + Temp[0, i]);
            }

            // Console.WriteLine("mass=" + totmass);
            Console.ReadLine();

        }
    }
}
