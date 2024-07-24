/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Computes invV b using sparse CGD - either for REML or computing PRS / lambda
//dichot=0, linear; dichot=1, logistic; dichot=2, quasi-logistic (not currently used)

///////////////////////////

void sparse_cgd(int ns, int total, int num_rels, int *firsts, int *seconds, double *kinships, double *cX, double *cR, int num_resps_use, double *pedhers, double *Yadj, int dichot, double *nullweights, int num_fixed, double *Z, int type, double tol, double her, double *polates, int nmcmc, double *gaussian, int ncal, double *cdata, double *cmults, double *pedgammas, double *pedsds)
//type=0 - standard, type=1 - MC REML, type=2 - invV and lambda
{
int i, j, g, m, p, p2, count, count2, cflag, cflag2;
double value, value2, value3, sum, sum2, sum3, sumsq, sumsq2, sumsq3, sumsq4, mean, var;

double *fractions, *cP, *cVP, *ctops, *cbots, *calphas, *cbetas;


fractions=malloc(sizeof(double)*total);
cP=malloc(sizeof(double)*ns*total);
cVP=malloc(sizeof(double)*ns*total);
ctops=malloc(sizeof(double)*total);
cbots=malloc(sizeof(double)*total);
calphas=malloc(sizeof(double)*total);
cbetas=malloc(sizeof(double)*total);

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

//see if necessary to update residuals
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
cR[(size_t)p*ns+firsts[i]]-=cX[(size_t)p*ns+seconds[i]]*kinships[i]*fractions[p];
if(firsts[i]!=seconds[i]){cR[(size_t)p*ns+seconds[i]]-=cX[(size_t)p*ns+firsts[i]]*kinships[i]*fractions[p];}
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
cR[(size_t)p*ns+firsts[i]]-=cX[(size_t)p*ns+seconds[i]]*kinships[i]*fractions[p];
if(firsts[i]!=seconds[i]){cR[(size_t)p*ns+seconds[i]]-=cX[(size_t)p*ns+firsts[i]]*kinships[i]*fractions[p];}
}
}
if(dichot==2)
{
//subtract (1-h2) invW cX
for(i=0;i<ns;i++){cR[(size_t)p*ns+i]-=cX[(size_t)p*ns+i]/nullweights[i+m*ns]*(1-fractions[p]);}

//subtract K h2 cX (for off diagonals, must subtract twice)
for(i=0;i<num_rels;i++)
{
cR[(size_t)p*ns+firsts[i]]-=cX[(size_t)p*ns+seconds[i]]*kinships[i]*fractions[p];
if(firsts[i]!=seconds[i]){cR[(size_t)p*ns+seconds[i]]-=cX[(size_t)p*ns+firsts[i]]*kinships[i]*fractions[p];}
}
}
}
}

////////

//set cP=cR
copy_matrix(ns, total, cR, cP, 0, NULL);

//tops are cR x cR
for(p=0;p<total;p++)
{
ctops[p]=0;
for(i=0;i<ns;i++){ctops[p]+=pow(cR[(size_t)p*ns+i],2);}
}

count=0;
while(1)
{
if(count==100){printf("Warning, CGD did not converge within %d iterations\n\n %f %f", count, ctops[0], ctops[1]);break;}

if(dichot==0)	//get V P, where V = K h2 + I (1-h2)
{
#pragma omp parallel for private(p,p2,m,i) schedule(static)
for(p=0;p<total;p++)
{
//set cVP = I (1-h2) cP
for(i=0;i<ns;i++){cVP[(size_t)p*ns+i]=cP[(size_t)p*ns+i]*(1-fractions[p]);}

//add on K h2 cP (for off diagonals, must add twice)
for(i=0;i<num_rels;i++)
{
cVP[(size_t)p*ns+firsts[i]]+=cP[(size_t)p*ns+seconds[i]]*kinships[i]*fractions[p];
if(firsts[i]!=seconds[i]){cVP[(size_t)p*ns+seconds[i]]+=cP[(size_t)p*ns+firsts[i]]*kinships[i]*fractions[p];}
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

//add on K sig cP (for off diagonals, must add twice)
for(i=0;i<num_rels;i++)
{
cVP[(size_t)p*ns+firsts[i]]+=cP[(size_t)p*ns+seconds[i]]*kinships[i]*fractions[p];
if(firsts[i]!=seconds[i]){cVP[(size_t)p*ns+seconds[i]]+=cP[(size_t)p*ns+firsts[i]]*kinships[i]*fractions[p];}
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

//add on K h2 cP (for off diagonals, must add twice)
for(i=0;i<num_rels;i++)
{
cVP[(size_t)p*ns+firsts[i]]+=cP[(size_t)p*ns+seconds[i]]*kinships[i]*fractions[p];
if(firsts[i]!=seconds[i]){cVP[(size_t)p*ns+seconds[i]]+=cP[(size_t)p*ns+firsts[i]]*kinships[i]*fractions[p];}
}
}
}

//bottoms are cP x cVP
for(p=0;p<total;p++)
{
cbots[p]=0;
for(i=0;i<ns;i++){cbots[p]+=cP[(size_t)p*ns+i]*cVP[(size_t)p*ns+i];}
}

//alphas are tops/bottoms
for(p=0;p<total;p++){calphas[p]=ctops[p]/cbots[p];}

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
for(p=0;p<total;p++)
{
ctops[p]=0;
for(i=0;i<ns;i++){ctops[p]+=pow(cR[(size_t)p*ns+i],2);}
}

//betas are tops/bottoms
for(p=0;p<total;p++){cbetas[p]=ctops[p]/cbots[p];}

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
for(p=0;p<total;p++){cflag+=(fabs(ctops[p])<tol);printf("Model %d has squared residuals %.4f\n", p+1, fabs(ctops[p]));}

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
cflag+=(fabs(ctops[p])<tol);
p=1+m*(2+ncal);
cflag+=(fabs(ctops[p])<tol);

if(ncal>0)
{
//ridge scaling is GTG/(GT invV G) x var(invV Y) / t(Y) invV Y x n
//for dichot=1, it is simply GTWG/(GT invV G)
//for dichot=2, it is GTWG/(GT invV G) x var(1/sqrt(W) invV Y) / t(Y) invV Y x n

sum=0;sum2=0;count2=0;
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
}

sum+=value;
sum2+=pow(value,2);
count2++;
}
}

mean=sum/count2;
var=sum2/count2-pow(mean,2);

cflag2+=(fabs(mean-pedgammas[m])<tol);
pedgammas[m]=mean;
pedsds[m]=pow(var/count2,.5);
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

count++;
}

free(fractions);free(cP);free(cVP);free(ctops);free(cbots);free(calphas);free(cbetas);
}

///////////////////////////

