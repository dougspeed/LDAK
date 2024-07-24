/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Code I've found for the Nelder-Mead algorithm (I use to get mle estimates for a gamma distribution)

///////////////////////////

//Mainly taken from https://people.sc.fsu.edu/~jburkardt/c_src/asa047/asa047.c

double gam_like(double x[2], double S, double T, double cut)
//returns minus log likelihood (divided by n) for truncated gamma
{
int info;
double alpha, beta, prob, value;

alpha=x[0];beta=x[1];
prob=(1-gamain(beta*cut,alpha,&info));
if(info!=0){printf("Error with gamain, please tell Doug; %d\n\n", info);exit(1);}
value=beta*S-(alpha-1)*T-alpha*log(beta)+log(prob)+lgamma(alpha);

return(value);
}	//end of gam_like

////////////////////////////

double gam_like2(double x[2], double *stats, int count2, int count3)
//returns sum of squared differences between log pvalues and log quantiles
//consider top count3 out of count2 positive stats
{
int j, info, count, wcount;
double alpha, beta, *probs, sum;

alpha=x[0];beta=x[1];
probs=malloc(sizeof(double)*count3);
count=0;sum=0;wcount=0;
for(j=0;j<count3;j++)
{
probs[j]=1-gamain(stats[j]*beta, alpha, &info);
if(probs[j]>0){probs[j]=log(probs[j]);}
else{probs[j]=-100;}
if(info>2)
{
if(wcount<5){printf("Warning with Gamma Function; stat %.2f, alpha %.2f, beta %.2f\n", stats[j], alpha, beta);}
wcount++;
}
else{sum+=pow(probs[j]-log(j+1)+log(count2+1),2);count++;}
}
if(wcount>5){printf("In total, problems with %d out of %d genes/chunks\n", wcount, count3);}
if(wcount>0){printf("\n");}

free(probs);

return(sum*count3/count);
}	//end of gam_like2

////////////////////////////

int nm_optim(double *alpha, double *beta, double *stats, int count2, int count3)
{
//will be using count3 out of count2 stats
double ccoeff = 0.5, del, dn, dnn, ecoeff = 2.0, eps = 0.001;
int i, ihi, ilo, j, jcount, l, nn;
double *p, *p2star, *pbar, *pstar, rcoeff = 1.0, rq, x, *y, y2star, ylo, ystar, z;

int n=2;	//number of variables
double start[2]={.5,.5};
double xmin[2];
double ynewlo[1];
double reqmin=1.0E-08;
double step[2]={1,1};

int konvge = 10;
int kcount = 1000;
int icount[1], numres[1], ifault[1];

/*
  Check the input parameters.
*/
  if ( reqmin <= 0.0 )
  {
    *ifault = 1;
    return(*ifault);
  }

  if ( n < 1 )
  {
    *ifault = 1;
    return(*ifault);
  }

  if ( konvge < 1 )
  {
    *ifault = 1;
    return(*ifault);
  }

//will need to make sure p, pstar, p2star and start stay positive

  p = ( double * ) malloc ( n * ( n + 1 ) * sizeof ( double ) );
  pstar = ( double * ) malloc ( n * sizeof ( double ) );
  p2star = ( double * ) malloc ( n * sizeof ( double ) );
  pbar = ( double * ) malloc ( n * sizeof ( double ) );
  y = ( double * ) malloc ( ( n + 1 ) * sizeof ( double ) );

  *icount = 0;
  *numres = 0;

  jcount = konvge; 
  dn = ( double ) ( n );
  nn = n + 1;
  dnn = ( double ) ( nn );
  del = 1.0;
  rq = reqmin * dn;
/*
  Initial or restarted loop.
*/
  for ( ; ; )
  {
    for ( i = 0; i < n; i++ )
    { 
      p[i+n*n] = start[i];
    }
//    y[n] = gam_like(start, S, T, cut);
    y[n] = gam_like2(start, stats, count2, count3);
    *icount = *icount + 1;

    for ( j = 0; j < n; j++ )
    {
      x = start[j];
      start[j] = start[j] + step[j] * del;
if(start[j]<=0){start[j]=0.0001;}
      for ( i = 0; i < n; i++ )
      {
        p[i+j*n] = start[i];
      }
//      y[j] = gam_like(start, S, T, cut);
      y[j] = gam_like2(start, stats, count2, count3);
      *icount = *icount + 1;
      start[j] = x;
    }
/*                 
  The simplex construction is complete.
                    
  Find highest and lowest Y values.  YNEWLO = Y(IHI) indicates
  the vertex of the simplex to be replaced.
*/                
    ylo = y[0];
    ilo = 0;

    for ( i = 1; i < nn; i++ )
    {
      if ( y[i] < ylo )
      {
        ylo = y[i];
        ilo = i;
      }
    }
/*
  Inner loop.
*/
    for ( ; ; )
    {
      if ( kcount <= *icount )
      {
        break;
      }
      *ynewlo = y[0];
      ihi = 0;

      for ( i = 1; i < nn; i++ )
      {
        if ( *ynewlo < y[i] )
        {
          *ynewlo = y[i];
          ihi = i;
        }
      }
/*
  Calculate PBAR, the centroid of the simplex vertices
  excepting the vertex with Y value YNEWLO.
*/
      for ( i = 0; i < n; i++ )
      {
        z = 0.0;
        for ( j = 0; j < nn; j++ )
        { 
          z = z + p[i+j*n];
        }
        z = z - p[i+ihi*n];  
        pbar[i] = z / dn;
      }
/*
  Reflection through the centroid.
*/
      for ( i = 0; i < n; i++ )
      {
        pstar[i] = pbar[i] + rcoeff * ( pbar[i] - p[i+ihi*n] );
if(pstar[i]<=0){pstar[i]=0.0001;}
      }
//      ystar = gam_like(pstar, S, T, cut);
      ystar = gam_like2(pstar, stats, count2, count3);
      *icount = *icount + 1;
/*
  Successful reflection, so extension.
*/
      if ( ystar < ylo )
      {
        for ( i = 0; i < n; i++ )
        {
          p2star[i] = pbar[i] + ecoeff * ( pstar[i] - pbar[i] );
if(p2star[i]<=0){p2star[i]=0.0001;}
        }
//        y2star = gam_like(p2star, S, T, cut);
      y2star = gam_like2(p2star, stats, count2, count3);
        *icount = *icount + 1;
/*
  Check extension.
*/
        if ( ystar < y2star )
        {
          for ( i = 0; i < n; i++ )
          {
            p[i+ihi*n] = pstar[i];
          }
          y[ihi] = ystar;
        }
/*
  Retain extension or contraction.
*/
        else
        {
          for ( i = 0; i < n; i++ )
          {
            p[i+ihi*n] = p2star[i];
          }
          y[ihi] = y2star;
        }
      }
/*
  No extension.
*/
      else
      {
        l = 0;
        for ( i = 0; i < nn; i++ )
        {
          if ( ystar < y[i] )
          {
            l = l + 1;
          }
        }

        if ( 1 < l )
        {
          for ( i = 0; i < n; i++ )
          {
            p[i+ihi*n] = pstar[i];
          }
          y[ihi] = ystar;
        }
/*
  Contraction on the Y(IHI) side of the centroid.
*/
        else if ( l == 0 )
        {
          for ( i = 0; i < n; i++ )
          {
            p2star[i] = pbar[i] + ccoeff * ( p[i+ihi*n] - pbar[i] );
if(p2star[i]<=0){p2star[i]=0.0001;}
          }
//          y2star = gam_like(p2star, S, T, cut);
      y2star = gam_like2(p2star, stats, count2, count3);
          *icount = *icount + 1;
/*
  Contract the whole simplex.
*/
          if ( y[ihi] < y2star )
          {
            for ( j = 0; j < nn; j++ )
            {
              for ( i = 0; i < n; i++ )
              {
                p[i+j*n] = ( p[i+j*n] + p[i+ilo*n] ) * 0.5;
                xmin[i] = p[i+j*n];
              }
//              y[j] = gam_like(xmin, S, T, cut);
      y[j] = gam_like2(xmin, stats, count2, count3);
              *icount = *icount + 1;
            }
            ylo = y[0];
            ilo = 0;

            for ( i = 1; i < nn; i++ )
            {
              if ( y[i] < ylo )
              {
                ylo = y[i];
                ilo = i;
              }
            }
            continue;
          }
/*
  Retain contraction.
*/
          else
          {
            for ( i = 0; i < n; i++ )
            {
              p[i+ihi*n] = p2star[i];
            }
            y[ihi] = y2star;
          }
        }
/*
  Contraction on the reflection side of the centroid.
*/
        else if ( l == 1 )
        {
          for ( i = 0; i < n; i++ )
          {
            p2star[i] = pbar[i] + ccoeff * ( pstar[i] - pbar[i] );
if(p2star[i]<=0){p2star[i]=0.0001;}
          }
//          y2star = gam_like(p2star, S, T, cut);
          y2star = gam_like2(p2star, stats, count2, count3);
          *icount = *icount + 1;
/*
  Retain reflection?
*/
          if ( y2star <= ystar )
          {
            for ( i = 0; i < n; i++ )
            {
              p[i+ihi*n] = p2star[i];
            }
            y[ihi] = y2star;
          }
          else
          {
            for ( i = 0; i < n; i++ )
            {
              p[i+ihi*n] = pstar[i];
            }
            y[ihi] = ystar;
          }
        }
      }
/*
  Check if YLO improved.
*/
      if ( y[ihi] < ylo )
      {
        ylo = y[ihi];
        ilo = ihi;
      }
      jcount = jcount - 1;

      if ( 0 < jcount )
      {
        continue;
      }
/*
  Check to see if minimum reached.
*/
      if ( *icount <= kcount )
      {
        jcount = konvge;

        z = 0.0;
        for ( i = 0; i < nn; i++ )
        {
          z = z + y[i];
        }
        x = z / dnn;

        z = 0.0;
        for ( i = 0; i < nn; i++ )
        {
          z = z + pow ( y[i] - x, 2 );
        }

        if ( z <= rq )
        {
          break;
        }
      }
    }
/*
  Factorial tests to check that YNEWLO is a local minimum.
*/
    for ( i = 0; i < n; i++ )
    {
      xmin[i] = p[i+ilo*n];
    }
    *ynewlo = y[ilo];

    if ( kcount < *icount )
    {
      *ifault = 2;
      break;
    }

    *ifault = 0;

    for ( i = 0; i < n; i++ )
    {
      del = step[i] * eps;
      xmin[i] = xmin[i] + del;
//      z = gam_like(xmin, S, T, cut);
      z = gam_like2(xmin, stats, count2, count3);
      *icount = *icount + 1;
      if ( z < *ynewlo )
      {
        *ifault = 2;
        break;
      }
      xmin[i] = xmin[i] - del - del;
//      z = gam_like(xmin, S, T, cut);
      z = gam_like2(xmin, stats, count2, count3);
      *icount = *icount + 1;
      if ( z < *ynewlo )
      {
        *ifault = 2;
        break;
      }
      xmin[i] = xmin[i] + del;
    }

    if ( *ifault == 0 )
    {
      break;
    }
/*
  Restart the procedure.
*/
    for ( i = 0; i < n; i++ )
    {
      start[i] = xmin[i];
    }
    del = eps;
    *numres = *numres + 1;
  }
  free ( p );
  free ( pstar );
  free ( p2star );
  free ( pbar );
  free ( y );

if(*ifault==0){*alpha=xmin[0];*beta=xmin[1];}

return(*ifault);
}	//end of nm_optim

///////////////////////////

