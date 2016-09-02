using System;
using System.Collections.Generic;
using System.Linq;
using System.Windows.Forms;

using Latticeb.Models;


namespace Latticeb
{
    static class Program
    {
        [STAThread]
        static void Main()
        {
            const int tsteps = 10000;  // time steps

            // internal condition
            float u0 = 0.48f;

            // initial boundary conditions, Temperature = 1.
            float u1 = 0.35f;
            float rho1 = 1.0f;

            float u2 = -0.28f;
            float rho2 = 1.0f;

            const int nlattice = 161; // space nodes quantity for lattice
            float[,] flattice = new float[nlattice, 6];

            //--hybrid declaration
            D1Q6 lattice;

            for (int i = 1; i < nlattice - 1; ++i)
            {
                flattice[i, 0] = D1Q6._n_eq_k(1.0f, u0, 0);
                flattice[i, 1] = D1Q6._n_eq_k(1.0f, u0, 1);
                flattice[i, 2] = D1Q6._n_eq_k(1.0f, u0, 2);
                flattice[i, 3] = D1Q6._n_eq_k(1.0f, u0, 3);
                flattice[i, 4] = D1Q6._n_eq_k(1.0f, u0, 4);
                flattice[i, 5] = D1Q6._n_eq_k(1.0f, u0, 5);
            }

            for (int time = 0; time < tsteps; ++time)
            {

                //--TODO: set LBE condition (flattice[,] = F(flattice, fleft, fright) )
                flattice[0, 3] = D1Q6._n_eq_k(rho1, u1, 3);
                flattice[0, 4] = D1Q6._n_eq_k(rho1, u1, 4);
                flattice[0, 5] = D1Q6._n_eq_k(rho1, u1, 5);

                flattice[nlattice - 1, 2] = D1Q6._n_eq_k(rho2, u2, 2);
                flattice[nlattice - 1, 1] = D1Q6._n_eq_k(rho2, u2, 1);
                flattice[nlattice - 1, 0] = D1Q6._n_eq_k(rho2, u2, 0);

                lattice = new D1Q6(flattice, 1.0f);
                flattice = lattice.Solve();
            }

            //--results

            //--diagnostics lattice
            for (int i = 0; i < nlattice; i++)
            {
                float vrho = flattice[i, 4] - flattice[i, 1]
                    + 0.5f * (flattice[i, 3] - flattice[i, 2])
                    + 3.0f * (flattice[i, 5] - flattice[i, 0]);

                float _rho = flattice[i, 0] + flattice[i, 1] +
                    flattice[i, 2] + flattice[i, 3] +
                    flattice[i, 4] + flattice[i, 5];

                Console.WriteLine(flattice[i, 0] + " " + flattice[i, 1] + " " + flattice[i, 2] +
                            " " + flattice[i, 3] + " " + flattice[i, 4] + " " + flattice[i, 5] + " x" + i * 0.5
                            + " v =" + vrho / _rho);

            }

            Console.ReadLine();
        }
    }
}
