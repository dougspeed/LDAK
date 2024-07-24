/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Code I've found for sampling from a Normal (and exponential) Distribution
//I made my own version of RNOR that does not use global static variables
//at bottom are codes for normal cdfs (but now I instead use erfc)

///////////////////////////

/* The ziggurat method for RNOR and REXP
Combine the code below with the main program in which you want
normal or exponential variates.   Then use of RNOR in any expression
will provide a standard normal variate with mean zero, variance 1,
while use of REXP in any expression will provide an exponential variate
with density exp(-x),x>0.
Before using RNOR or REXP in your main, insert a command such as
zigset(86947731 );
with your own choice of seed value>0, rather than 86947731.
(If you do not invoke zigset(...) you will get all zeros for RNOR and REXP.)
For details of the method, see Marsaglia and Tsang, "The ziggurat method
for generating random variables", Journ. Statistical Software.
*/

/* Have cdf and inverse cdf using code from
http://www.wilmott.com/messageview.cfm?catid=44&threadid=95982
//  Copyright (c) 2012-2013 M.A. (Thijs) van den Berg, http://sitmo.com/
//
//  Use, modification and distribution are subject to the MIT Software License. 
//  
//  The MIT License (MIT)
//  Permission is hereby granted, free of charge, to any person obtaining a copy
//  of this software and associated documentation files (the "Software"), to deal
//  in the Software without restriction, including without limitation the rights
//  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//  copies of the Software, and to permit persons to whom the Software is
//  furnished to do so, subject to the following conditions:
//  
//  The above copyright notice and this permission notice shall be included in
//  all copies or substantial portions of the Software.
//  
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
//  THE SOFTWARE.

// This file contains two normal distribution samplers that use the inverse 
// cumulative function to generate normal distributed samples.
//
// The inversion routine is based on the work of Peter John Acklam
// http://home.online.no/~pjacklam/notes/invnorm
// 
// sitmo::normal_distribution_inv              Full machine precision
// sitmo::normal_distribution_inv_single       Single precision, Err < 1.15E-9
*/

/* Also not used but have another cdf function for standard normal
Code obtained from http://www.sitmo.com/doc/Calculating_the_Cumulative_Normal_Distribution
based on Abromowitz and Stegun approximation approximation explained
http://www.math.sfu.ca/~cbm/aands/page_932.htm*/

///////////////////////////

static unsigned int kn[128];
static float wn[128], fn[128];

void zigset_safe(unsigned long jsrseed)	//edited this simply to remove exponential function parts
{  const double m1 = 2147483648.0;
   double dn=3.442619855899,tn=dn,vn=9.91256303526217e-3, q;
   int i;

   q=vn/exp(-.5*dn*dn);
   kn[0]=(dn/q)*m1;
   kn[1]=0;

   wn[0]=q/m1;
   wn[127]=dn/m1;

   fn[0]=1.;
   fn[127]=exp(-.5*dn*dn);

    for(i=126;i>=1;i--)
    {dn=sqrt(-2.*log(vn/dn+exp(-.5*dn*dn)));
     kn[i+1]=(dn/tn)*m1;
     tn=dn;
     fn[i]=exp(-.5*dn*dn);
     wn[i]=dn/m1;
    }
}

double rnorm_safe()	//this is my own version of RNOR and nfix (below) that avoids using global static variables
{
const float r = 3.442620f;

int inum;
unsigned int unum;
float unifrand, unifrand2, x, y;


while(1)
{
//unifrand=(double)rand()/RAND_MAX;
unifrand=genrand_real1();
if(unifrand<0.5){inum=rand();}
else{inum=-rand();}

unum=inum&127;
if(fabs(inum)<kn[unum]){return(inum*wn[unum]);}

//first attempt did not return, so try to "fix"

if(unum==0)	//will always return
{
x=1;y=0;
while(y+y<x*x)
{
//unifrand=(double)rand()/RAND_MAX;
//unifrand2=(double)rand()/RAND_MAX;
unifrand=genrand_real1();
unifrand2=genrand_real1();
x=-log(1.0-unifrand)*0.2904764;
y=-log(1.0-unifrand2);
}

if(inum>0){return(r+x);}
else{return(-r-x);}
}

//so unum>0 - might be able to return, else start process again
x=inum*wn[unum];
//unifrand=(double)rand()/RAND_MAX;
unifrand=genrand_real1();
if(fn[unum]+unifrand*(fn[unum-1]-fn[unum])<exp(-.5*x*x)){return(x);}
}

}

///////////////////////////

//cdf of Normal code 1 - this one looks more impressive!
    
        double normal_cdf(double x)
        {
            double cdf;
            double poly;

            double xabs = fabs(x);
            if (xabs > 37.0) 
                cdf = 0.0;
            else {  
                double exponential = exp( -xabs*xabs / 2.0);
                if (xabs < 7.07106781186547) { 
                    poly = 3.52624965998911E-02 * xabs + 0.700383064443688;
                    poly = poly * xabs + 6.37396220353165;
                    poly = poly * xabs + 33.912866078383;
                    poly = poly * xabs + 112.079291497871;
                    poly = poly * xabs + 221.213596169931;
                    poly = poly * xabs + 220.206867912376;
                    cdf = exponential * poly;
                    poly = 8.83883476483184E-02 * xabs + 1.75566716318264;
                    poly = poly * xabs + 16.064177579207;
                    poly = poly * xabs + 86.7807322029461;
                    poly = poly * xabs + 296.564248779674;
                    poly = poly * xabs + 637.333633378831;
                    poly = poly * xabs + 793.826512519948;
                    poly = poly * xabs + 440.413735824752;
                    cdf = cdf / poly;
                } else {
                    poly = xabs + 0.65;
                    poly = xabs + 4 / poly;
                    poly = xabs + 3 / poly;
                    poly = xabs + 2 / poly;
                    poly = xabs + 1 / poly;
                    cdf = exponential / poly / 2.506628274631;
                }
            }
            if (x>0) cdf = 1.0 - cdf;
            return cdf;
        } // normal_cdf
        
        
        //-----------------------------------------
        // Standard Inverse Cumulative Normal Distribution
        // BASED ON http://home.online.no/~pjacklam/notes/invnorm/
        //-----------------------------------------

        double normal_inv_single(double p)
        {
            double x;
            double q, r;

            if ( (0.0 < p)  && (p < 0.02425) ) {
                q = sqrt( -2.0 * log(p) );
                x = (((((   -7.784894002430293e-03 *q 
                            -3.223964580411365e-01)*q
                            -2.400758277161838e+00)*q
                            -2.549732539343734e+00)*q
                            +4.374664141464968e+00)*q
                            +2.938163982698783e+00) 
                    /
                    ((((    +7.784695709041462e-03 *q
                            +3.224671290700398e-01)*q
                            +2.445134137142996e+00)*q
                            +3.754408661907416e+00)*q
                            +1.0
                    );
            } else {
                if ( (0.97575 < p) && (p < 1.0) ) {
                    q = sqrt( -2.0 * log(1.0-p) );
                    x = -(((((  -7.784894002430293e-03 *q 
                                -3.223964580411365e-01)*q
                                -2.400758277161838e+00)*q
                                -2.549732539343734e+00)*q
                                +4.374664141464968e+00)*q
                                +2.938163982698783e+00) 
                        /
                        ((((    +7.784695709041462e-03 *q
                                +3.224671290700398e-01)*q
                                +2.445134137142996e+00)*q
                                +3.754408661907416e+00)*q
                                +1.0
                        );
                } else {
                  if ((0.02425 <= p) && (p <= 0.97575)) {

                           q = p - 0.5;
                           r = q*q;
                           x = (((((    -3.969683028665376e+01*r
                                        +2.209460984245205e+02)*r
                                        -2.759285104469687e+02)*r
                                        +1.383577518672690e+02)*r
                                        -3.066479806614716e+01)*r
                                        +2.506628277459239e+00)*q 
                                /
                                (((((   -5.447609879822406e+01*r
                                        +1.615858368580409e+02)*r
                                        -1.556989798598866e+02)*r
                                        +6.680131188771972e+01)*r
                                        -1.328068155288572e+01)*r
                                        +1.0
                                );
                   
                    } else {
                    }
                }
            }
            return x;
        }
        
        double normal_inv(double p)
        {
            double x = normal_inv_single(p);
            //One iteration of Halleys rational method 
            //(third order) gives full machine precision.
            double e,u;
            if ( (0.0 < p) && (p < 1.0) ) {
                e = normal_cdf(x) - p;
                u = e * 2.50662827463100050242 * exp(x*x/2.0);
                x = x - u/(1.0 + x*u/2.0);
            }
            return x;
        }
        
///////////////////////////

//cdf of Normal code 2

double cdfN(const double x)
{
  const double b1 =  0.319381530;
  const double b2 = -0.356563782;
  const double b3 =  1.781477937;
  const double b4 = -1.821255978;
  const double b5 =  1.330274429;
  const double p  =  0.2316419;
  const double c  =  0.39894228;

  if(x >= 0.0) {
      double t = 1.0 / ( 1.0 + p * x );
      return (1.0 - c * exp( -x * x / 2.0 ) * t *
      ( t *( t * ( t * ( t * b5 + b4 ) + b3 ) + b2 ) + b1 ));
  }
  else {
      double t = 1.0 / ( 1.0 - p * x );
      return ( c * exp( -x * x / 2.0 ) * t *
      ( t *( t * ( t * ( t * b5 + b4 ) + b3 ) + b2 ) + b1 ));
    }
}

///////////////////////////

//old code for normal (and exponential) functions

/*
static unsigned int jz,jsr=123456789;

#define SHR3 (jz=jsr, jsr^=(jsr<<13), jsr^=(jsr>>17), jsr^=(jsr<<5),jz+jsr)
#define UNI (.5 + (signed) SHR3*.2328306e-9)
#define IUNI SHR3

static int hz;
//static unsigned int iz, kn[128], ke[256];
static unsigned int iz, ke[256];
//static float wn[128],fn[128],we[256],fe[256];
static float we[256],fe[256];

#define RNOR (hz=SHR3, iz=hz&127, (fabs(hz)<kn[iz])? hz*wn[iz] : nfix())
#define REXP (jz=SHR3, iz=jz&255, (    jz <ke[iz])? jz*we[iz] : efix())

//nfix() generates variates from the residue when rejection in RNOR occurs
float nfix(void)
{
const float r = 3.442620f;     //The start of the right tail
static float x, y;
 for(;;)
  {  x=hz*wn[iz];      //iz==0, handles the base strip
     if(iz==0)
       { do{ x=-log(UNI)*0.2904764; y=-log(UNI);}	//.2904764 is 1/r
        while(y+y<x*x);
        return (hz>0)? r+x : -r-x;
       }

// iz>0, handle the wedges of other strips
      if( fn[iz]+UNI*(fn[iz-1]-fn[iz]) < exp(-.5*x*x) ) return x;

//initiate, try to exit for(;;) for loop
      hz=SHR3;
      iz=hz&127;
      if(fabs(hz)<kn[iz]) return (hz*wn[iz]);
  }

}

//efix() generates variates from the residue when rejection in REXP occurs
float efix(void)
{ float x;
 for(;;)
  {  if(iz==0) return (7.69711-log(UNI));          //iz==0
     x=jz*we[iz]; if( fe[iz]+UNI*(fe[iz-1]-fe[iz]) < exp(-x) ) return (x);

//initiate, try to exit for(;;) loop
   jz=SHR3;
   iz=(jz&255);
   if(jz<ke[iz]) return (jz*we[iz]);
  }
}

//This procedure sets the seed and creates the tables
void zigset(unsigned long jsrseed)
{  const double m1 = 2147483648.0, m2 = 4294967296.;
   double dn=3.442619855899,tn=dn,vn=9.91256303526217e-3, q;
   double de=7.697117470131487, te=de, ve=3.949659822581572e-3;
   int i;
   jsr^=jsrseed;

//Set up tables for RNOR
   q=vn/exp(-.5*dn*dn);
   kn[0]=(dn/q)*m1;
   kn[1]=0;

   wn[0]=q/m1;
   wn[127]=dn/m1;

   fn[0]=1.;
   fn[127]=exp(-.5*dn*dn);

    for(i=126;i>=1;i--)
    {dn=sqrt(-2.*log(vn/dn+exp(-.5*dn*dn)));
     kn[i+1]=(dn/tn)*m1;
     tn=dn;
     fn[i]=exp(-.5*dn*dn);
     wn[i]=dn/m1;
    }

//Set up tables for REXP
    q = ve/exp(-de);
    ke[0]=(de/q)*m2;
    ke[1]=0;

    we[0]=q/m2;
    we[255]=de/m2;

    fe[0]=1.;
    fe[255]=exp(-de);

   for(i=254;i>=1;i--)
  {de=-log(ve/de+exp(-de));
   ke[i+1]= (de/te)*m2;
   te=de;
   fe[i]=exp(-de);
   we[i]=de/m2;
  }
}
*/

///////////////////////////

