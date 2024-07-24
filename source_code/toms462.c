/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//code I've found for computing the cdf of a bivariate normal distribution
//taken from https://people.sc.fsu.edu/~jburkardt/c_src/toms462/toms462.html

///////////////////////////

/******************************************************************************/

double gauss ( double t )

/******************************************************************************/
/*
  Purpose:

    GAUSS returns the area of the lower tail of the normal curve.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 April 2012

  Author:

    John Burkardt

  Parameters:

    Input, double T, the evaluation point.

    Output, double GAUSS, the lower normal tail area.
*/
{
  double value;

  value = ( 1.0 + erf ( t / sqrt ( 2.0 ) ) ) / 2.0;

  return value;
}

////////

/******************************************************************************/

double bivnor ( double ah, double ak, double r )

/******************************************************************************/
/*
  Purpose:

    BIVNOR computes the bivariate normal CDF.

  Discussion:

    BIVNOR computes the probability for two normal variates X and Y
    whose correlation is R, that AH <= X and AK <= Y.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 April 2012

  Author:

    Original FORTRAN77 version by Thomas Donnelly.
    C version by John Burkardt.

  Reference:

    Thomas Donnelly,
    Algorithm 462: Bivariate Normal Distribution,
    Communications of the ACM,
    October 1973, Volume 16, Number 10, page 638.

  Parameters:

    Input, double AH, AK, the lower limits of integration.

    Input, double R, the correlation between X and Y.

    Output, double BIVNOR, the bivariate normal CDF.

  Local Parameters:

    Local, int IDIG, the number of significant digits
    to the right of the decimal point desired in the answer.
*/
{
  double a2;
  double ap;
  double b;
  double cn;
  double con;
  double conex;
  double ex;
  double g2;
  double gh;
  double gk;
  double gw;
  double h2;
  double h4;
  int i;
  static int idig = 15;
  int is;
  double rr;
  double s1;
  double s2;
  double sgn;
  double sn;
  double sp;
  double sqr;
  double t;
  static double twopi = 6.283185307179587;
  double w2;
  double wh;
  double wk;

  b = 0.0;

  gh = gauss ( - ah ) / 2.0;
  gk = gauss ( - ak ) / 2.0;

  if ( r == 0.0 )
  {
    b = 4.00 * gh * gk;
    b = fmax ( b, 0.0 );
    b = fmin ( b, 1.0 );
    return b;
  }

  rr = ( 1.0 + r ) * ( 1.0 - r );

  if ( rr < 0.0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "BIVNOR - Fatal error!\n" );
    fprintf ( stderr, "  1 < |R|.\n" );
    exit ( 0 );
  }

  if ( rr == 0.0 )
  {
    if ( r < 0.0 )
    {
      if ( ah + ak < 0.0 )
      {
        b = 2.0 * ( gh + gk ) - 1.0;
      }
    }
    else
    {
      if ( ah - ak < 0.0 )
      {
        b = 2.0 * gk;
      }
      else
      {
        b = 2.0 * gh;
      }
    }
    b = fmax ( b, 0.0 );
    b = fmin ( b, 1.0 );
    return b;
  }

  sqr = sqrt ( rr );

  if ( idig == 15 )
  {
    con = twopi * 1.0E-15 / 2.0;
  }
  else
  {
    con = twopi / 2.0;
    for ( i = 1; i <= idig; i++ )
    {
      con = con / 10.0;
    }
  }
/*
  (0,0)
*/
  if ( ah == 0.0 && ak == 0.0 )
  {
    b = 0.25 + asin ( r ) / twopi;
    b = fmax ( b, 0.0 );
    b = fmin ( b, 1.0 );
    return b;
  }
/*
  (0,nonzero)
*/
  if ( ah == 0.0 && ak != 0.0 )
  {
    b = gk;
    wh = -ak;
    wk = ( ah / ak - r ) / sqr;
    gw = 2.0 * gk;
    is = 1;
  }
/*
  (nonzero,0)
*/
  else if ( ah != 0.0 && ak == 0.0 )
  {
    b = gh;
    wh = -ah;
    wk = ( ak / ah - r ) / sqr;
    gw = 2.0 * gh;
    is = -1;
  }
/*
  (nonzero,nonzero)
*/
  else if ( ah != 0.0 && ak != 0.0 )
  {
    b = gh + gk;
    if ( ah * ak < 0.0 )
    {
      b = b - 0.5;
    }
    wh = - ah;
    wk = ( ak / ah - r ) / sqr;
    gw = 2.0 * gh;
    is = -1;
  }

  for ( ; ; )
  {
    sgn = -1.0;
    t = 0.0;

    if ( wk != 0.0 )
    {
      if ( fabs ( wk ) == 1.0 )
      {
        t = wk * gw * ( 1.0 - gw ) / 2.0;
        b = b + sgn * t;
      }
      else
      {
        if ( 1.0 < fabs ( wk ) )
        {
          sgn = -sgn;
          wh = wh * wk;
          g2 = gauss ( wh );
          wk = 1.0 / wk;

          if ( wk < 0.0 )
          {
            b = b + 0.5;
          }
          b = b - ( gw + g2 ) / 2.0 + gw * g2;
        }
        h2 = wh * wh;
        a2 = wk * wk;
        h4 = h2 / 2.0;
        ex = exp ( - h4 );
        w2 = h4 * ex;
        ap = 1.0;
        s2 = ap - ex;
        sp = ap;
        s1 = 0.0;
        sn = s1;
        conex = fabs ( con / wk );

        for ( ; ; )
        {
          cn = ap * s2 / ( sn + sp );
          s1 = s1 + cn;

          if ( fabs ( cn ) <= conex )
          {
            break;
          }
          sn = sp;
          sp = sp + 1.0;
          s2 = s2 - w2;
          w2 = w2 * h4 / sp;
          ap = - ap * a2;
        }
        t = ( atan ( wk ) - wk * s1 ) / twopi;
        b = b + sgn * t;
      }
    }
    if ( 0 <= is )
    {
      break;
    }
    if ( ak == 0.0 )
    {
      break;
    }
    wh = -ak;
    wk = ( ah / ak - r ) / sqr;
    gw = 2.0 * gk;
    is = 1;
  }

  b = fmax ( b, 0.0 );
  b = fmin ( b, 1.0 );

  return b;
}

///////////////////////////
