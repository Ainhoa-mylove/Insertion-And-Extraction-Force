using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Analysis of the IAEF
{
    public class calculate
    {
//==============================================================================================================
public static double NQ(double Φ)
        {
            double Nq;
            Nq = Math.Exp(Math.PI * Math.Tan(Φ*(Math.PI/180))) *Math.Pow(Math.Tan((45+ Φ/2) *(Math.PI/180)),2);
            return Nq;

        }
//attention:Tan(Φ/180*Math.pi)=Tan(Φ*(Math.PI/180)) It only can be used as Tan(Φ*(Math.PI/180)) cause C#/.Net grammar
        public static double NC(double Φ)
        {
            double Nc;
            if (Φ == 0)
            {
                Nc = 5.14;
            }
            else
            {
                Nc = (calculate.NQ(Φ) - 1.0) / Math.Tan(Φ*(Math.PI/180));
            }
                
            return Nc;
        }

        public static double NR(double Φ)
        {
            double Nr;
            Nr = (calculate.NQ(Φ) - 1) * Math.Tan(1.4 * Φ *(Math.PI/180));
            return Nr;
        }
//NQ,NR,NC coefficient of foundation bearing capacity,they connect with Φ（Φ is Angle of soil friction）
//==============================================================================================================
        public static double SC(double  B, double  L, double  Φ)
        {
            double Sc;
            Sc = 1.0 + 0.2 * (B / L) * Math.Pow(Math.Tan((45 + Φ / 2) * (Math.PI/180)),2);
            return Sc;
        }

        public static double SQ(double  B, double  L, double  Φ)
        {
            double Sq;
            if (Φ == 0)
            {
                Sq = 1.0;

            }
            else
            {
                Sq = 1 + 0.1 * (B / L) * Math.Pow(Math.Tan((45+ Φ/2)*(Math.PI/180)),2);

            }
            return Sq;
        }

        public static double SR(double  B, double  L, double  Φ)
        {
            double Sr;
            Sr = calculate.SQ(B, L, Φ);
            return Sr;
        }
//SC,SQ,SR are Base shape correction factor.they connect with B,L,Φ(B，L are shipshoe-size.Φ is Angle of soil friction)
//==============================================================================================================
        public static double DC(double D, double  B,double  Φ)
        {
            double Dc;
            Dc = 1.0 + 0.2 * (D / B) * Math.Tan((45+ Φ/2)*(Math.PI/180));
            return Dc;
        }

        public static double DQ(double D, double  B, double  Φ)
        {
            double Dq;
            if (Φ == 0)
            {
                Dq = 1.0;
            }
            else
            {
                Dq=1.0+0.1* (D / B) * Math.Tan((45 + Φ / 2) * (Math.PI/180));
            }
            return Dq;
        }
//DC,DQ are Foundation depth correction factor.there is an other factor named DR which is not calculating.cause it has the same
  values with DQ(parameter D is insert deeping.)
//==============================================================================================================       
        public static double Q(double D, double γ)
        {
            double q;
            q = D * γ;
            return q;

        }
//==============================================================================================================
        public static List<double > QU(double c,double[] q,double  B,double γ,double nq, double nc, double nr, double sc, double sq, double sr,double[] dc, double[] dq, double[] dr)
        {
            List<double> Qu = new List<double>();
            for(int i=0;i<dc.Length;)
            {
                Qu.Add(c * nc * sc * dc[i] + q[i] * nq * sq * dq[i] + 0.5 * γ * B * nr * sr * dr[i]);
                i += 1;
            }
            return Qu;

        }
        public static List<double> Final_results (double B, double L, double[] qu)
        {
            List<double> Fianl = new List<double>();
            for (int i = 0; i < qu.Length;)
            {
                Fianl.Add(B*L*qu[i]/1000);
                i += 1;
            }
            return Fianl;
         }

    }
}
//=======================================================================================================================
