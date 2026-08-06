/*
Copyright 2026 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Top computes invV b using sparse CGD - either for REML or computing PRS / lambda
//dichot=0, linear; dichot=1, logistic; dichot=2, quasi-logistic (not currently used)

///////////////////////////

void sparse_cgd(int ns, int total, int num_rels, int *firsts, int *seconds, double *relats, double *cX, double *cR, int num_resps_use, double *pedhers, double *Yadj, int dichot, double *nullweights, int type, double tol, double her, double *polates, int nmcmc, double *gaussian, int ncal, double *cdata, double *cmults, double *pedgammas, double *pedsds, double *pedeffs, int nhers, double *tryhers, double *ratios, double *X, double *Y)
//type=0 - standard, type=1 - MC REML, type=2 - invV and lambda, type=3 - get MC REML denominator, type=4 - get MC REML numerator 
//type=5 - convergence based on test statistic (not used)
{
int i, j, g, m, p, p2, count, count2, cflag, cflag2, one=1;
double value, value2, value3, sum, sum2, sum3, sumsq, sumsq2, sumsq3, sumsq4, mean, var;

double *fractions, *cP, *cVP, *ctops, *cbots, *calphas, *cbetas, *ratios2, *stats;


fractions=malloc(sizeof(double)*total);
cP=malloc(sizeof(double)*ns*total);
cVP=malloc(sizeof(double)*ns*total);
ctops=malloc(sizeof(double)*total);
cbots=malloc(sizeof(double)*total);
calphas=malloc(sizeof(double)*total);
cbetas=malloc(sizeof(double)*total);
ratios2=malloc(sizeof(double)*total);
stats=malloc(sizeof(double)*total);

//fill fractions

if(type==0)	//just assume her constant for each phenotype
{
for(p=0;p<total;p++)
{
p2=p%(total/num_resps_use);
m=(p-p2)*num_resps_use/total;
fractions[p]=pedhers[m];
}
}

if(type==1)
{
if(dichot==0||dichot==2)	//use her for all models
{
for(p=0;p<total;p++){fractions[p]=her;}
}
else	//must work out phenotype-specific sigmas
{
for(m=0;m<num_resps_use;m++)
{
sum=0;for(i=0;i<ns;i++){sum+=pow(nullweights[i+m*ns],-1);}
mean=sum/ns;
for(g=0;g<nmcmc+1;g++){fractions[g+m*(nmcmc+1)]=mean*her/(1-her);}
}
}
}

if(type==2)
{
if(dichot==0||dichot==2)	//her is constant for each phenotype
{
for(p=0;p<total;p++)
{
p2=p%(2+ncal);
m=(p-p2)/(2+ncal);
fractions[p]=pedhers[m];
}
}
else	//sigma is constant for each phenotype
{
for(p=0;p<total;p++)
{
p2=p%(2+ncal);
m=(p-p2)/(2+ncal);

sum=0;for(i=0;i<ns;i++){sum+=pow(nullweights[i+m*ns],-1);}
mean=sum/ns;
fractions[p]=mean*pedhers[m]/(1-pedhers[m]);
}
}
}

if(type==3) //use her for all models
{
for(p=0;p<total;p++){fractions[p]=her;}
}

if(type==4) //load from tryhers (once each)
{
for(p=0;p<total;p++)
{
fractions[p]=tryhers[p];
}
}

if(type==5)	//load from pedhers
{
for(p=0;p<total;p++)
{
fractions[p]=pedhers[p];
}
}

//see if necessary to update residuals (only required if starting values in cX are not zero)
for(p=0;p<total;p++)
{
count=0;for(i=0;i<ns;i++){count+=(cX[(size_t)p*ns+i]!=0);}

if(count>0)	//cR becomes cR - V cX
{
if(dichot==0)
{
//subtract I (1-h2) cX
for(i=0;i<ns;i++){cR[(size_t)p*ns+i]-=cX[(size_t)p*ns+i]*(1-fractions[p]);}

//subtract K h2 cX (for off diagonals, must subtract twice)
for(i=0;i<num_rels;i++)
{
cR[(size_t)p*ns+firsts[i]]-=cX[(size_t)p*ns+seconds[i]]*relats[i]*fractions[p];
if(firsts[i]!=seconds[i]){cR[(size_t)p*ns+seconds[i]]-=cX[(size_t)p*ns+firsts[i]]*relats[i]*fractions[p];}
}
}
if(dichot==1)
{
p2=p%(total/num_resps_use);
m=(p-p2)*num_resps_use/total;

//subtract invW cX
for(i=0;i<ns;i++){cR[(size_t)p*ns+i]-=cX[(size_t)p*ns+i]/nullweights[i+m*ns];}

//subtract K sig cX (for off diagonals, must subtract twice)
for(i=0;i<num_rels;i++)
{
cR[(size_t)p*ns+firsts[i]]-=cX[(size_t)p*ns+seconds[i]]*relats[i]*fractions[p];
if(firsts[i]!=seconds[i]){cR[(size_t)p*ns+seconds[i]]-=cX[(size_t)p*ns+firsts[i]]*relats[i]*fractions[p];}
}
}
if(dichot==2)
{
//subtract (1-h2) invW cX
for(i=0;i<ns;i++){cR[(size_t)p*ns+i]-=cX[(size_t)p*ns+i]/nullweights[i+m*ns]*(1-fractions[p]);}

//subtract K h2 cX (for off diagonals, must subtract twice)
for(i=0;i<num_rels;i++)
{
cR[(size_t)p*ns+firsts[i]]-=cX[(size_t)p*ns+seconds[i]]*relats[i]*fractions[p];
if(firsts[i]!=seconds[i]){cR[(size_t)p*ns+seconds[i]]-=cX[(size_t)p*ns+firsts[i]]*relats[i]*fractions[p];}
}
}
}
}

////////

//set cP=cR
copy_matrix(ns, total, cR, cP, 0, NULL);

//tops are cR x cR
for(p=0;p<total;p++){ctops[p]=ddot_(&ns, cR+(size_t)p*ns, &one, cR+(size_t)p*ns, &one);}

count=0;
while(1)
{
if(count==100){printf("Warning, CGD did not converge within %d iterations\n\n", count);break;}

if(dichot==0)	//get V P, where V = K h2 + I (1-h2)
{
#pragma omp parallel for private(p,p2,m,i) schedule(static)
for(p=0;p<total;p++)
{
//set cVP = I (1-h2) cP
for(i=0;i<ns;i++){cVP[(size_t)p*ns+i]=cP[(size_t)p*ns+i]*(1-fractions[p]);}

if(fractions[p]!=0&&ctops[p]>0)	//add on K h2 cP (for off diagonals, must add twice)
{
for(i=0;i<num_rels;i++)
{
cVP[(size_t)p*ns+firsts[i]]+=cP[(size_t)p*ns+seconds[i]]*relats[i]*fractions[p];
if(firsts[i]!=seconds[i]){cVP[(size_t)p*ns+seconds[i]]+=cP[(size_t)p*ns+firsts[i]]*relats[i]*fractions[p];}
}
}
}
}
if(dichot==1)	//get V P, where V = K sig + invW
{
#pragma omp parallel for private(p,p2,m,i) schedule(static)
for(p=0;p<total;p++)
{
p2=p%(total/num_resps_use);
m=(p-p2)*num_resps_use/total;

//set cVP = invW cP
for(i=0;i<ns;i++){cVP[(size_t)p*ns+i]=cP[(size_t)p*ns+i]/nullweights[i+m*ns];}

if(fractions[p]!=0&&ctops[p]>0)	//add on K sig cP (for off diagonals, must add twice)
{
for(i=0;i<num_rels;i++)
{
cVP[(size_t)p*ns+firsts[i]]+=cP[(size_t)p*ns+seconds[i]]*relats[i]*fractions[p];
if(firsts[i]!=seconds[i]){cVP[(size_t)p*ns+seconds[i]]+=cP[(size_t)p*ns+firsts[i]]*relats[i]*fractions[p];}
}
}
}
}
if(dichot==2)	//get V P, where V = K h2 + invW (1-h2)
{
#pragma omp parallel for private(p,p2,m,i) schedule(static)
for(p=0;p<total;p++)
{
p2=p%(total/num_resps_use);
m=(p-p2)*num_resps_use/total;

//set cVP = (1-h2) invW cP
for(i=0;i<ns;i++){cVP[(size_t)p*ns+i]=cP[(size_t)p*ns+i]/nullweights[i+m*ns]*(1-fractions[p]);}

if(fractions[p]!=0&&ctops[p]>0)	//add on K h2 cP (for off diagonals, must add twice)
{
for(i=0;i<num_rels;i++)
{
cVP[(size_t)p*ns+firsts[i]]+=cP[(size_t)p*ns+seconds[i]]*relats[i]*fractions[p];
if(firsts[i]!=seconds[i]){cVP[(size_t)p*ns+seconds[i]]+=cP[(size_t)p*ns+firsts[i]]*relats[i]*fractions[p];}
}
}
}
}

//bottoms are cP x cVP
for(p=0;p<total;p++){cbots[p]=ddot_(&ns, cP+(size_t)p*ns, &one, cVP+(size_t)p*ns, &one);}

//alphas are tops/bottoms
for(p=0;p<total;p++)
{
if(cbots[p]>0){calphas[p]=ctops[p]/cbots[p];}
else{calphas[p]=0;}
}

//cX = cX + alpha cP
#pragma omp parallel for private(p,i) schedule(static)
for(p=0;p<total;p++)
{
for(i=0;i<ns;i++){cX[(size_t)p*ns+i]+=calphas[p]*cP[(size_t)p*ns+i];}
}

//cR = cR - alpha cVP
#pragma omp parallel for private(p,i) schedule(static)
for(p=0;p<total;p++)
{
for(i=0;i<ns;i++){cR[(size_t)p*ns+i]-=calphas[p]*cVP[(size_t)p*ns+i];}
}

//move tops into bottoms
for(p=0;p<total;p++){cbots[p]=ctops[p];}

//tops are cR x cR
for(p=0;p<total;p++){ctops[p]=ddot_(&ns, cR+(size_t)p*ns, &one, cR+(size_t)p*ns, &one);}

//betas are tops/bottoms
for(p=0;p<total;p++)
{
if(cbots[p]>0){cbetas[p]=ctops[p]/cbots[p];}
else{cbetas[p]=0;}
}

//cP = cR + beta cP
#pragma omp parallel for private(p,i) schedule(static)
for(p=0;p<total;p++)
{
for(i=0;i<ns;i++){cP[(size_t)p*ns+i]=cR[(size_t)p*ns+i]+cbetas[p]*cP[(size_t)p*ns+i];}
}

////////

if(type==0)	//check for convergence based on residuals
{
cflag=0;
for(p=0;p<total;p++){cflag+=(fabs(ctops[p]/ns)<tol);}

printf("Completed Scan %d; %d of the %d models have converged\n", count+1, cflag, total);

if(cflag==total){printf("\n");break;}
}

if(type==1)	//check for convergence based on delta ratios
{
cflag=0;
for(m=0;m<num_resps_use;m++)
{
if(dichot==0)	//numerator is invV Y (Y - invV Y (1-h2)) / invV Y invV Y for real phenotype
{
p=nmcmc+m;
sumsq=0;sumsq2=0;
for(i=0;i<ns;i++)
{
sumsq+=cX[(size_t)p*ns+i]*(Yadj[i+m*ns]-cX[(size_t)p*ns+i]*(1-her));
sumsq2+=pow(cX[(size_t)p*ns+i],2);
}
value=sumsq/sumsq2;
}
if(dichot==1)	//numerator is invV Y (Y - invW invV Y) for real phenotype
{
p=nmcmc+m*(nmcmc+1);
sumsq=0;
for(i=0;i<ns;i++)
{
sumsq+=cX[(size_t)p*ns+i]*(Yadj[i+m*ns]-cX[(size_t)p*ns+i]/nullweights[i+m*ns]);
}
value=sumsq;
}
if(dichot==2)	//numerator is invV Y (Y - invW invV Y (1-h2)) / invV Y invW invV Y for real phenotype
{
p=nmcmc+m;
sumsq=0;sumsq2=0;
for(i=0;i<ns;i++)
{
sumsq+=cX[(size_t)p*ns+i]*(Yadj[i+m*ns]-cX[(size_t)p*ns+i]/nullweights[i+m*ns]*(1-her));
sumsq2+=pow(cX[(size_t)p*ns+i],2)/nullweights[i+m*ns];
}
value=sumsq/sumsq2;
}

sum=0;
for(g=0;g<nmcmc;g++)
{
if(dichot==0)	//denominator is mean of r (r - invV Y (1-h2)) / r invV r for random vectors
{
p=g;
sumsq=0;sumsq2=0;
for(i=0;i<ns;i++)
{
sumsq+=gaussian[(size_t)g*ns+i]*(gaussian[(size_t)g*ns+i]-cX[(size_t)p*ns+i]*(1-her));
sumsq2+=gaussian[(size_t)g*ns+i]*cX[(size_t)p*ns+i];
}
sum+=sumsq/sumsq2;
}
if(dichot==1)	//denominator is mean of r (r - invW invV Y) for random vectors
{
p=g+m*(nmcmc+1);
sumsq=0;
for(i=0;i<ns;i++)
{
sumsq+=gaussian[(size_t)g*ns+i]*(gaussian[(size_t)g*ns+i]-cX[(size_t)p*ns+i]/nullweights[i+m*ns]);
}
sum+=sumsq;
}
if(dichot==2)	//denominator is mean of r (r - invW invV Y (1-h2)) / r invW invV r for random vectors
{
p=g;
sumsq=0;sumsq2=0;
for(i=0;i<ns;i++)
{
sumsq+=gaussian[(size_t)g*ns+i]*(gaussian[(size_t)g*ns+i]-cX[(size_t)p*ns+i]/nullweights[i+m*ns]*(1-her));
sumsq2+=gaussian[(size_t)g*ns+i]*cX[(size_t)p*ns+i]/nullweights[i+m*ns];
}
sum+=sumsq/sumsq2;
}
}
value2=sum/nmcmc;

value3=value/value2;

if(count>0)	//see if converged
{
if(fabs(value3-polates[m])<tol){cflag++;}
}
polates[m]=value3;
}

if(cflag==num_resps_use){break;}
}	//end of type=1

if(type==2)	//check for convergence based on residuals of invV Y and lambda
{
cflag=0;
cflag2=0;

for(m=0;m<num_resps_use;m++)
{
p=m*(2+ncal);
cflag+=(fabs(ctops[p]/ns)<tol);
p=1+m*(2+ncal);
cflag+=(fabs(ctops[p]/ns)<tol);

if(ncal>0)
{
//ridge scaling is GTG/(GT invV G) x var(invV Y) / t(Y) invV Y x n
//for dichot=1, it is simply GTWG/(GT invV G)
//for dichot=2, it is GTWG/(GT invV G) x var(1/sqrt(W) invV Y) / t(Y) invV Y x n
//effect calibration always first ratio

sum=0;sum2=0;value3=0;count2=0;
for(j=0;j<ncal;j++)
{
if(cmults[j]!=-9999)
{
p=m*(2+ncal);
p2=2+j+m*(2+ncal);

if(dichot==0)	//sumsq is t(G) G; sumsq2 is t(G) invVG; sum3 + sumsq3 give var(invV Y); sumsq4 is t(Y) invV Y
{
sum3=0;sumsq=0;sumsq2=0;sumsq3=0;sumsq4=0;
for(i=0;i<ns;i++)
{
sum3+=cX[(size_t)p*ns+i];
sumsq+=pow(cdata[i+(j+m*ncal)*ns],2);
sumsq2+=cdata[i+(j+m*ncal)*ns]*cX[(size_t)p2*ns+i];
sumsq3+=pow(cX[(size_t)p*ns+i],2);
sumsq4+=Yadj[i+m*ns]*cX[(size_t)p*ns+i];
}
value=sumsq/sumsq2*(sumsq3-sum3/ns*sum3)/sumsq4;
value2=sumsq/sumsq2;
}
if(dichot==1)	//sumsq is t(G) W G; sumsq2 is t(G) invVG
{
sumsq=0;sumsq2=0;
for(i=0;i<ns;i++)
{
sumsq+=pow(cdata[i+(j+m*ncal)*ns],2)*nullweights[i+m*ns];
sumsq2+=cdata[i+(j+m*ncal)*ns]*cX[(size_t)p2*ns+i];
}
value=sumsq/sumsq2;
value2=sumsq/sumsq2;
}
if(dichot==2)	//sumsq is t(G) W G; sumsq2 is t(G) invVG; sum3 + sumsq3 give var(1/sqrt(W) invV Y); sumsq4 is t(Y) invV Y
{
sum3=0;sumsq=0;sumsq2=0;sumsq3=0;sumsq4=0;
for(i=0;i<ns;i++)
{
sum3+=cX[(size_t)p*ns+i]*pow(nullweights[i+m*ns],-.5);
sumsq+=pow(cdata[i+(j+m*ncal)*ns],2)*nullweights[i+m*ns];
sumsq2+=cdata[i+(j+m*ncal)*ns]*cX[(size_t)p2*ns+i];
sumsq3+=pow(cX[(size_t)p*ns+i],2)/nullweights[i+m*ns];
sumsq4+=Yadj[i+m*ns]*cX[(size_t)p*ns+i];
}
value=sumsq/sumsq2*(sumsq3-sum3/ns*sum3)/sumsq4;
value2=sumsq/sumsq2;
}

sum+=value;
sum2+=pow(value,2);
value3+=value2;
count2++;
}
}

mean=sum/count2;
var=sum2/count2-pow(mean,2);
if(fabs(var)<1e-8){var=0;}

cflag2+=(fabs(mean-pedgammas[m])<tol);

pedgammas[m]=mean;
pedsds[m]=pow(var/count2,.5);
pedeffs[m]=value3/count2;
}
}	//end of m loop

if(ncal>0)
{
if(cflag+cflag2==3*num_resps_use){break;}
}
else
{
if(cflag==2*num_resps_use){break;}
}
}	//end of type=2

if(type==3)	//check for convergence based on bottom ratio (must have dichot=0)
{
if(count>0) //save current value
{value=ratios[0];}

//denominator is (truncated) mean of r (r - invV Y (1-h2)) / r invV r for random vectors
for(g=0;g<nmcmc;g++)
{
p=g;
sumsq=0;sumsq2=0;
for(i=0;i<ns;i++)
{
sumsq+=gaussian[(size_t)g*ns+i]*(gaussian[(size_t)g*ns+i]-cX[(size_t)p*ns+i]*(1-her));
sumsq2+=gaussian[(size_t)g*ns+i]*cX[(size_t)p*ns+i];
}
ratios2[g]=sumsq/sumsq2;
}

if(nmcmc<5) //get regular mean
{
sum=0;for(g=0;g<nmcmc;g++){sum+=ratios2[g];}
ratios[0]=sum/nmcmc;
}
else    //get mean with min and max excluded
{
qsort(ratios2, nmcmc, sizeof(double), compare_double);
sum=0;for(g=1;g<nmcmc-1;g++){sum+=ratios2[g];}
ratios[0]=sum/(nmcmc-2);
}

if(count>0)	//see if converged
{
if(fabs(ratios[0]-value)<tol){break;}
}
}	//end of type=3

if(type==4)	//check for convergence based on top ratios (must have dichot=0)
{
cflag=0;

for(j=0;j<nhers;j++)
{
//numerator is invV Y (Y - invV Y (1-h2)) / invV Y invV Y for real phenotype
p=j;
sumsq=0;sumsq2=0;
for(i=0;i<ns;i++)
{
sumsq+=cX[(size_t)p*ns+i]*(Yadj[i]-cX[(size_t)p*ns+i]*(1-tryhers[j]));
sumsq2+=pow(cX[(size_t)p*ns+i],2);
}
value=sumsq/sumsq2;

if(count>0)	//see if converged
{
if(fabs(value-ratios[j])<tol){cflag++;}
}
ratios[j]=value;
}

if(cflag==nhers){break;}
}	//end of type=4

if(type==5) //check for convergence based on (approximate) test statistic (might be trivial predictors)
{
cflag=0;

for(p=0;p<total;p++)
{
value=ddot_(&ns, cX+(size_t)p*ns, &one, Y+(size_t)p*ns, &one);
value2=ddot_(&ns, cX+(size_t)p*ns, &one, X+(size_t)p*ns, &one);
if(value2>0){value3=pow(value,2)/value2;}
else{value3=0;}

if(count>0)	//see if converged
{
if(fabs(value3-stats[p])<tol){cflag++;}
}
stats[p]=value3;
}

if(cflag==total){break;}
}   //end of type=5

count++;
}	//end of while loop

free(fractions);free(cP);free(cVP);free(ctops);free(cbots);free(calphas);free(cbetas);free(ratios2);free(stats);
}

///////////////////////////

