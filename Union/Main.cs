//Rextester.Program.Main is the entry point for your code. Don't change it.
//Compiler version 4.0.30319.17929 for Microsoft (R) .NET Framework 4.5

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.RegularExpressions;


namespace Rextester
{
    public class Program
    {
        [STAThread]
        public static void Main(string[] args)
        {
            // General
            const int tsteps = 1000;  // time steps for initial run of BGK / lattice
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
                int tstep = tsteps - 1;
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
            for (int iteration = 0; iteration < NumIterations; ++iteration)
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
                for (int i = nvlc / 2; i < nvlc; ++i)
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
                for (int i = 0; i < nvlc / 2; ++i)
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
                lm1 = - /* minus required */ D1Q6.moment1(flattice[0, 2], flattice[0, 1], flattice[0, 0]);
                lm2 = D1Q6.moment2(flattice[0, 2], flattice[0, 1], flattice[0, 0]);
                // to the right
                rm0 = D1Q6.moment0(flattice[Nlattice - 1, 3], flattice[Nlattice - 1, 4], flattice[Nlattice - 1, 5]);
                rm1 = D1Q6.moment1(flattice[Nlattice - 1, 3], flattice[Nlattice - 1, 4], flattice[Nlattice - 1, 5]);
                rm2 = D1Q6.moment2(flattice[Nlattice - 1, 3], flattice[Nlattice - 1, 4], flattice[Nlattice - 1, 5]);

                // Setting BGK half-distributions using m0-m2 for left and right subdomains respectively
                // TODO: check the code below, should lm2, rm2 not be used ?
                //--левая точка сращивания
                for (int i = 0; i < nvlc / 2; ++i)
                {
                    double phys_veloc = ((i - (nvlc - 1.0) / 2.0) * vw);
                    fleft[Xsteps - 1, i] = vw * 1.0f / (float) Math.Sqrt(2.0 * Math.PI * 1.0f)
                        * (float) Math.Exp( - phys_veloc * phys_veloc / (2.0))
                        * (float) (lm0 + lm1 * phys_veloc + (lm1 - lm0) * (phys_veloc * phys_veloc - 1.0));
                }

                //--правая точка сращивания
                for (int i = nvlc / 2; i < nvlc; ++i)
                {
                    double phys_veloc = ((i - (nvlc - 1.0) / 2.0) * vw);
                    fright[0, i] = vw * 1.0f / (float) Math.Sqrt(2.0 * Math.PI * 1.0f)
                        * (float) Math.Exp( - phys_veloc * phys_veloc / (2.0))
                        * (float) (rm0 + rm1 * phys_veloc + (rm1 - rm0) * (phys_veloc * phys_veloc - 1.0));
                }
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
