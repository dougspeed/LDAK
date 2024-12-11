/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Code I've found for evaluating the digamma, trigamma and incomplete gamma functions and sampling from gamma dist
//the latter is based on G. Marsaglia and W. Tsang. A simple method for generating gamma variables, 2000.

///////////////////////////

/* C implementations of useful functions.
 * Written by Tom Minka (unless otherwise noted).
 */

/* The digamma function is the derivative of gammaln.

   Reference:
    J Bernardo,
    Psi ( Digamma ) Function,
    Algorithm AS 103,
    Applied Statistics,
    Volume 25, Number 3, pages 315-317, 1976.

    From http://www.psc.edu/~burkardt/src/dirichlet/dirichlet.f
    (with modifications for negative numbers and extra precision)
*/

/* The trigamma function is the derivative of the digamma function.

   Reference:

    B Schneider,
    Trigamma Function,
    Algorithm AS 121,
    Applied Statistics, 
    Volume 27, Number 1, page 97-99, 1978.

    From http://www.psc.edu/~burkardt/src/dirichlet/dirichlet.f
    (with modification for negative arguments and extra precision)
*/

//At bottom are two implementations of incomplete gamma function
//http://people.sc.fsu.edu/~jburkardt/c_src/asa032/asa032.html
//http://people.sc.fsu.edu/~jburkardt/c_src/asa147/asa147.html

////////////////////////////

double gamain ( double x, double p, int *ifault );
void gamma_inc_values ( int *n_data, double *a, double *x, double *fx );
void timestamp ( );

////////////////////////////

double digamma(double x)
{
//double neginf = -INFINITY;
  static const double c = 12,
    digamma1 = -0.57721566490153286,
    trigamma1 = 1.6449340668482264365, /* pi^2/6 */
    s = 1e-6,
    s3 = 1./12,
    s4 = 1./120,
    s5 = 1./252,
    s6 = 1./240,
    s7 = 1./132;
//s8 = 691./32760, s9 = 1./12, s10 = 3617./8160;
  double result;

//DOUG has turned off some lines
  /* Illegal arguments  */
//  if((x == neginf) || isnan(x)) { return NAN; }
  /* Singularities */
//  if((x <= 0) && (floor(x) == x)) { return neginf; }

  /* Negative values */
  /* Use the reflection formula (Jeffrey 11.1.6):
   * digamma(-x) = digamma(x+1) + pi*cot(pi*x)
   *
   * This is related to the identity
   * digamma(-x) = digamma(x+1) - digamma(z) + digamma(1-z)
   * where z is the fractional part of x
   * For example:
   * digamma(-3.1) = 1/3.1 + 1/2.1 + 1/1.1 + 1/0.1 + digamma(1-0.1)
   *               = digamma(4.1) - digamma(0.1) + digamma(1-0.1)
   * Then we use
   * digamma(1-z) - digamma(z) = pi*cot(pi*z)
   */
  if(x < 0) {return digamma(1-x) + M_PI/tan(-M_PI*x);}

  /* Use Taylor series if argument <= S */
  if(x <= s) {return digamma1 - 1/x + trigamma1*x;}

  /* Reduce to digamma(X + N) where (X + N) >= C */
  result = 0;
  while(x < c) { result -= 1/x; x++; }

  /* Use de Moivre's expansion if argument >= C */
  /* This expansion can be computed in Maple via asympt(Psi(x),x) */
  if(x >= c) {
    double r = 1/x, t;
    result += log(x) - 0.5*r;
    r *= r;
#if 0
    result -= r * (s3 - r * (s4 - r * (s5 - r * (s6 - r * s7))));
#else
    /* this version for lame compilers */
    t = (s5 - r * (s6 - r * s7));
    result -= r * (s3 - r * (s4 - r * t));
#endif
  }
  return result;
}

////////////////////////////

double trigamma(double x)
{
//  double neginf = -INFINITY,
  double    small = 1e-4,
    large = 8,
    trigamma1 = 1.6449340668482264365, /* pi^2/6 = Zeta(2) */
    tetragamma1 = -2.404113806319188570799476,  /* -2 Zeta(3) */
    b2 =  1./6,  /* B_2 */
    b4 = -1./30, /* B_4 */
    b6 =  1./42, /* B_6 */
    b8 = -1./30, /* B_8 */
    b10 = 5./66; /* B_10 */
  double result;

//DOUG has turned off some lines
  /* Illegal arguments */
//  if((x == neginf) || isnan(x)) { return NAN; }
  /* Singularities */
//  if((x <= 0) && (floor(x) == x)) { return neginf; }

  /* Negative values */
  /* Use the derivative of the digamma reflection formula:
   * -trigamma(-x) = trigamma(x+1) - (pi*csc(pi*x))^2
   */
  if(x < 0) { result = M_PI/sin(-M_PI*x); return -trigamma(1-x) + result*result; }

  /* Use Taylor series if argument <= small */
  if(x <= small) { return 1/(x*x) + trigamma1 + tetragamma1*x; }

  /* Reduce to trigamma(x+n) where ( X + N ) >= B */
  result = 0;
  while(x < large) { result += 1/(x*x); x++; }

  /* Apply asymptotic formula when X >= B */
  /* This expansion can be computed in Maple via asympt(Psi(1,x),x) */
  if(x >= large) {
    double r = 1/(x*x), t;
#if 0
    result += 0.5*r + (1 + r*(b2 + r*(b4 + r*(b6 + r*(b8 + r*b10)))))/x;
#else
    t = (b4 + r*(b6 + r*(b8 + r*b10)));
    result += 0.5*r + (1 + r*(b2 + r*t))/x;
#endif
  }
  return result;
}

////////////////////////////

double gamain ( double x, double p, int *ifault )

/*
  Purpose:

    GAMAIN computes the incomplete gamma ratio.

  Discussion:

    A series expansion is used if P > X or X <= 1.  Otherwise, a
    continued fraction approximation is used.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    29 June 2014

  Author:

    Original FORTRAN77 version by G Bhattacharjee.
    C version by John Burkardt.

  Reference:

    G Bhattacharjee,
    Algorithm AS 32:
    The Incomplete Gamma Integral,
    Applied Statistics,
    Volume 19, Number 3, 1970, pages 285-287.

  Parameters:

    Input, double X, P, the parameters of the incomplete 
    gamma ratio.  0 <= X, and 0 < P.

    Output, int *IFAULT, error flag.
    0, no errors.
    1, P <= 0.
    2, X < 0.
    3, underflow.
    4, error return from the Log Gamma routine.

    Output, double GAMAIN, the value of the incomplete gamma ratio.
*/
{
  double a;
  double acu = 1.0E-08;
  double an;
  double arg;
  double b;
  double dif;
  double factor;
  double g;
  double gin;
  int i;
  double oflo = 1.0E+37;
  double pn[6];
  double rn;
  double term;
  double uflo = 1.0E-37;
  double value;

  *ifault = 0;
/*
  Check the input.
*/
  if ( p <= 0.0 )
  {
    *ifault = 1;
    value = 0.0;
    return value;
  }

  if ( x < 0.0 )
  {
    *ifault = 2;
    value = 0.0;
    return value;
  }

  if ( x == 0.0 )
  {
    *ifault = 0;
    value = 0.0;
    return value;
  }

  g = lgamma ( p );

  arg = p * log ( x ) - x - g;

  if ( arg < log ( uflo ) )
  {
    *ifault = 3;
    value = 0.0;
    return value;
  }

  *ifault = 0;
  factor = exp ( arg );
/*
  Calculation by series expansion.
*/
  if ( x <= 1.0 || x < p )
  {
    gin = 1.0;
    term = 1.0;
    rn = p;

    for ( ; ; )
    {
      rn = rn + 1.0;
      term = term * x / rn;
      gin = gin + term;

      if ( term <= acu )
      {
        break;
      }
    }

    value = gin * factor / p;
    return value;
  }
/*
  Calculation by continued fraction.
*/
  a = 1.0 - p;
  b = a + x + 1.0;
  term = 0.0;

  pn[0] = 1.0;
  pn[1] = x;
  pn[2] = x + 1.0;
  pn[3] = x * b;

  gin = pn[2] / pn[3];

  for ( ; ; )
  {
    a = a + 1.0;
    b = b + 2.0;
    term = term + 1.0;
    an = a * term;
    for ( i = 0; i <= 1; i++ )
    {
      pn[i+4] = b * pn[i+2] - an * pn[i];
    }

    if ( pn[5] != 0.0 )
    {
      rn = pn[4] / pn[5];
      dif = fabs ( gin - rn );
/*
  Absolute error tolerance satisfied?
*/
      if ( dif <= acu )
      {
/*
  Relative error tolerance satisfied?
*/
        if ( dif <= acu * rn )
        {
          value = 1.0 - factor * gin;
          break;
        }
      }
      gin = rn;
    }

    for ( i = 0; i < 4; i++ )
    {
      pn[i] = pn[i+2];
    }

    if ( oflo <= fabs ( pn[4] ) )
    {
      for ( i = 0; i < 4; i++ )
      {
        pn[i] = pn[i] / oflo;
      }
    }
  }

  return value;
}

////////////////////////////

double gammds ( double x, double p, int *ifault )	//worse than gamain

/*
  Purpose:

    GAMMDS computes the incomplete Gamma integral.

  Discussion:

    The parameters must be positive.  An infinite series is used.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 November 2010

  Author:

    Original FORTRAN77 version by Chi Leung Lau.
    C version by John Burkardt.

  Reference:

    Chi Leung Lau,
    Algorithm AS 147:
    A Simple Series for the Incomplete Gamma Integral,
    Applied Statistics,
    Volume 29, Number 1, 1980, pages 113-114.

  Parameters:

    Input, double X, P, the arguments of the incomplete
    Gamma integral.  X and P must be greater than 0.

    Output, int *IFAULT, error flag.
    0, no errors.
    1, X <= 0 or P <= 0.
    2, underflow during the computation.

    Output, double GAMMDS, the value of the incomplete
    Gamma integral.
*/
{
  double a;
  double arg;
  double c;
  double e = 1.0E-09;
  double f;
  double uflo = 1.0E-37;
  double value;
/*
  Check the input.
*/
  if ( x <= 0.0 )
  {
    *ifault = 1;
    value = 0.0;
    return value;
  }

  if ( p <= 0.0 )
  {
    *ifault = 1;
    value = 0.0;
    return value;
  }
/*
  LGAMMA is the natural logarithm of the gamma function.
*/
  arg = p * log ( x ) - lgamma ( p + 1.0 ) - x;

  if ( arg < log ( uflo ) )
  {
    value = 0.0;
    *ifault = 2;
    return value;
  }

  f = exp ( arg );

  if ( f == 0.0 )
  {
    value = 0.0;
    *ifault = 2;
    return value;
  }

  *ifault = 0;
/*
  Series begins.
*/
  c = 1.0;
  value = 1.0;
  a = p;

  for ( ; ; )
  {
    a = a + 1.0;
    c = c * x / a;
    value = value + c;

    if ( c <= e * value )
    {
      break;
    }
  }

  value = value * f;

  return value;
}

////////////////////////////

double rgamma(double a)
{
int count;
double d,c,x,v,u,v3,x2;

d=a-1.0/3;
c=pow(9*d,-.5);

count=1;
while(1)
{
v=-1;
while(v<=0)
{
x=rnorm_safe();
v=1+c*x;
}
v3=v*v*v;
//u=(double)rand()/RAND_MAX;
u=genrand_real2();
x2=x*x;
if (u<1-0.0331*x2*x2){return (d*v3);}
if(log(u)<0.5*x2+d*(1-v3+log(v3))){return(d*v3);}

count++;
if(count==1000){printf("Error sampling from gamma distribution with shape %f, please tell Doug\n\n", a);exit(1);}
}

}

////////////////////////////

