using System;
using MathNet.Numerics.Distributions;

namespace SpreadOption
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
            double d1 = (Math.Log(S) + 0.5 * Math.Pow(a_kirk, 2) * T) / (a_kirk * Math.Sqrt(T));
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
            double I_t = Math.Sqrt(Math.Pow(a_kirk, 2)) + 0.5 * (Math.Pow((sigma2 * S2_0 / (S2_0 + K)) - rho * sigma1, 2)) * (1 / (Math.Pow(Math.Sqrt(Math.Pow(a_kirk, 2)), 3))) * Math.Pow(sigma2, 2) * ((S2_0 * K) / (Math.Pow(S2_0 + K, 2))) * (X_t - x_aster);

            //Here we calculate our d1 and d2
            double S = S1_0 / (S2_0 + K);
            double d1 = (Math.Log(S) + 0.5 * Math.Pow(I_t, 2) * T) / (I_t * Math.Sqrt(T));
            double d2 = d1 - I_t * Math.Sqrt(T);

            //Fit a normal distribution
            var normal = new Normal(0, 1);

            //#Here we use the above calculations to approximate the call spread using Kirk's formula
            double C_modified_kirk = S1_0 * normal.CumulativeDistribution(d1) - (S2_0 + K) * normal.CumulativeDistribution(d2);
            C_modified_kirk *= Math.Exp(-r * T);

            return C_modified_kirk;
        }

        public static double MC_Spread_Option(double S1_0, double S2_0,
                                              double K, double rho, double T, double r,
                                              double sigma1, double sigma2, double q,
                                              int numberofsim, int numberofstep)
        {
            double[,] S1 = new double[numberofsim, numberofstep];
            double[,] S2 = new double[numberofsim, numberofstep];
            double step = T / numberofstep;

            Random PRNG = new Random();
            double Gaussian_U;
            double Gaussian_V;

            double sum = 0.0;

            for (int j = 0; j < numberofsim; j++)
            {
                //Initialization
                S1[j, 0] = S1_0;
                S2[j, 0] = S2_0;
                double sum1_j = Math.Log(S1[j, 0]);
                double sum2_j = Math.Log(S2[j, 0]);

                for (int i = 1; i < numberofstep; i++)
                {
                    //Generate two independent standard Gaussian random variables
                    Gaussian_U = BoxMullerCart(PRNG);
                    Gaussian_V = BoxMullerCart(PRNG);

                    //Simulation of sample paths
                    S1[j, i] = S1[j, i - 1] * Math.Exp(sigma1 * Math.Sqrt(step) * Gaussian_U + sigma2 * Math.Sqrt(1 - Math.Pow(rho, 2)) * Math.Sqrt(step) * Gaussian_V + (r - q - 0.5 * Math.Pow(sigma1, 2) * step));
                    S2[j, i] = S2[j, i - 1] * Math.Exp(sigma2 * Math.Sqrt(step) * Gaussian_U + (r - q - 0.5 * Math.Pow(sigma2, 2) * step));

                    sum1_j += Math.Log(S1[j, i]);
                    sum2_j += Math.Log(S2[j, i]);

                }

                //Average on sample paths
                double S1_T = Math.Exp(sum1_j / numberofstep);
                double S2_T = Math.Exp(sum2_j / numberofstep);

                sum += Math.Exp(-r * T) * Math.Max(S1_T - S2_T - K, 0);
            }

            //Average on simulations
            return sum / (double)numberofsim;

        }

        public static double BoxMullerCart(Random PRNG, double mu = 0.0, double sig = 1.0)
        {

            double X = (PRNG.NextDouble() * 2 - 1);
            double Y = (PRNG.NextDouble() * 2 - 1);
            //double X = ((1 - PRNG.genrand_real2()) * 2 - 1);
            //double Y = ((1 - PRNG.genrand_real2()) * 2 - 1);

            double S = X * X + Y * Y;

            while (S >= 1.0 || S == 0.0)
            //while (S >= 1.0)
            {
                X = (PRNG.NextDouble() * 2 - 1);
                Y = (PRNG.NextDouble() * 2 - 1);
                S = X * X + Y * Y;
            }

            double Z1 = X * Math.Sqrt((-2 * Math.Log(S)) / S);
            double Z2 = Y * Math.Sqrt((-2 * Math.Log(S)) / S);

            double result1, result2;
            result1 = mu + sig * Z1;
            result2 = mu + sig * Z2;

            return result1;
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
            double q = 0.0;
            double rho = 0.999;
            int numberofsim = 10000;
            int numberofstep = 50;

            Console.WriteLine(" spread_option_price_CF = " + CF_Spread_Option.CF_Kirk(s1, s2, K, rho, T, r, sigma1, sigma2));
            Console.WriteLine(" spread_option_price_CF = " + CF_Spread_Option.CF_Modified_Kirk(s1, s2, K, rho, T, r, sigma1, sigma2));
            Console.WriteLine(" spread_option_price_MC = " + CF_Spread_Option.MC_Spread_Option(s1, s2, K, rho, T, r, sigma1, sigma2, q, numberofsim, numberofstep));

        }
    }
}