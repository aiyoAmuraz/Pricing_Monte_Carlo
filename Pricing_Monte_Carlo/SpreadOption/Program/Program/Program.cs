using System;
using MathNet.Numerics.Distributions;


namespace Program
{
    public static class CF_Spread_Option
    {
        public static double CF_Kirk(double S1_0, double S2_0,
                                     double K, double rho, double T, double r,
                                     double sigma1, double sigma2)
        {
            //Here we calculate our volatility
            double a_kirk = Math.Sqrt(Math.Pow(sigma1, 2) - 2 * rho * sigma1 * sigma2 * (S2_0 / (S2_0 + K)) + Math.Pow(sigma2 * (S2_0 / (S2_0 + K)), 2));

            //Here we calculate our d1 and d2
            double S = (S1_0 / (S2_0 + K));
            double d1 = (Math.Log(S) + 0.5 * Math.Pow(a_kirk,2) * T) / (a_kirk * Math.Sqrt(T));
            double d2 = d1 - a_kirk * Math.Sqrt(T);

            //Fit a normal distribution
            var normal = new Normal(0, 1);

            //#Here we use the above calculations to approximate the call spread using Kirk's formula
            double C_kirk = S1_0 * normal.CumulativeDistribution(d1) - (S2_0 + K) * normal.CumulativeDistribution(d2);
            C_kirk *= Math.Exp(-r * T);

            return C_kirk;
        }

        public static double CF_Modified_Kirk(double S1_0, double S2_0,
                                              double K, double rho, double T, double r,
                                              double sigma1, double sigma2)
        {
            //Here we calculate our volatility
            double a_kirk = Math.Sqrt(Math.Pow(sigma1, 2) - 2 * rho * sigma1 * sigma2 * (S2_0 / (S2_0 + K)) + Math.Pow(sigma2 * (S2_0 / (S2_0 + K)), 2));
            double X_t = Math.Log(S1_0);
            double x_aster = Math.Log(S2_0 + K);
            double I_t = Math.Sqrt(Math.Pow(a_kirk, 2)) + 0.5 * (Math.Pow((sigma2 * S2_0 / (S2_0 + K)) - rho * sigma1, 2)) * (1 / (Math.Pow(Math.Sqrt(Math.Pow(a_kirk, 2)), 3))) * Math.Pow(sigma2, 2) * ((S2_0 * K) / (Math.Pow(S2_0 + K,2))) * (X_t - x_aster);

            //Here we calculate our d1 and d2
            double S = S1_0 / (S2_0 + K);
            double d1 = (Math.Log(S) + 0.5 * Math.Pow(I_t,2) * T) / (I_t * Math.Sqrt(T));
            double d2 = d1 - I_t * Math.Sqrt(T);

            //Fit a normal distribution
            var normal = new Normal(0, 1);

            //#Here we use the above calculations to approximate the call spread using Kirk's formula
            double C_modified_kirk = S1_0 * normal.CumulativeDistribution(d1) - (S2_0 + K) * normal.CumulativeDistribution(d2);
            C_modified_kirk *= Math.Exp(-r * T);

            return C_modified_kirk;
        }

        static void Main(string[] args)
        {
            double s1 = 100.0;
            double s2 = 100.0;
            double T = 0.5;
            double sigma1 = 0.3;
            double sigma2 = 0.2;
            double K = 5.0;
            double r = 0.0;
            double rho = 0.999;

            Console.WriteLine(" spread_option_price_CF = " + CF_Spread_Option.CF_Kirk(s1, s2, K, rho, T, r, sigma1, sigma2));
            Console.WriteLine(" spread_option_price_CF = " + CF_Spread_Option.CF_Modified_Kirk(s1, s2, K, rho, T, r, sigma1, sigma2));

        }
    }
}