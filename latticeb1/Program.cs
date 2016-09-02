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
            int tsteps = 100, xsteps=10;

            float[,] finit=new float[xsteps,3];
            float[,] P = new  float[tsteps, xsteps];

            //--initial data - step-like;
            for(int i=0;i<finit.GetLength(0);i++)
            {
               if(i<5){
                  finit[i,0]=1/6f;
                  finit[i,1]=1*4/6f;
                  finit[i,2]=1*1/6f;
               }
               else{
                  finit[i,0]=2*1/6f;
                  finit[i,1]=2*4/6f;
                  finit[i,2]=2*1/6f;
               
               }
   
            }

            //-- solving
            D1Q3 flow = new D1Q3(finit, 1.0f, tsteps);
            flow.Solve();
            P = flow.P();

            //--results
            float totmass = 0f;
            for (int i = 0; i < xsteps; i++)
            {
                Console.WriteLine("density=" + P[tsteps-1, i]);
                totmass = totmass + P[tsteps - 1, i];
            }

            Console.WriteLine("mass=" + totmass);
            Console.ReadLine();

         //   Application.EnableVisualStyles();
         //   Application.SetCompatibleTextRenderingDefault(false);
         //   Application.Run(new Form1());
        }
    }
}
