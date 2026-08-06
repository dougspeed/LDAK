/*
Copyright 2026 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//allocate and set lambdas

///////////////////////////

lambdas=malloc(sizeof(double)*num_try);
lambdas2=malloc(sizeof(double)*num_try);
lambdas3=malloc(sizeof(double)*num_try);
lambdas4=malloc(sizeof(double)*num_try);

//set lambdas (values corresponding to exps=her=varphen=1 - will later scale)
for(p=0;p<num_try;p++)
{
//start at zero, then change when required
lambdas[p]=0;lambdas2[p]=0;lambdas3[p]=0;lambdas4[p]=0;

if(trytypes[p]==1)	//lasso-sparse - lambda is set to match the lassosum paper
{lambdas[p]=trylams[p]*pow(data_length,-.5);}

if(trytypes[p]==2)	//lasso-shrink - exp(betaj^2) = 2/lam^2
{lambdas[p]=pow(2,.5);}

if(trytypes[p]==3)	//ridge - exp(betaj^2) = lam
{lambdas[p]=1.0;}

if(trytypes[p]==4)	//bolt - exp(betaj^2) = plam + p2lam2, with p2lam2/plam = f2/(1-f2) - can have p=1 and p2=0
{
if(tryps[p]==1)
{lambdas[p]=1;lambdas2[p]=0;}
else
{
lambdas[p]=(1-tryf2s[p])/tryps[p];
lambdas2[p]=tryf2s[p]/tryp2s[p];
}
}

if(trytypes[p]==5)	//bayesr-sparse - exp(betaj^2) = p2lam2 + p3lam3 + p4lam4
{
value=tryp2s[p]/100+tryp3s[p]/10+tryp4s[p];
lambdas[p]=0.0;
lambdas2[p]=0.01/value;
lambdas3[p]=0.1/value;
lambdas4[p]=1.0/value;
}

if(trytypes[p]==6)	//bayesr-shrink - exp(betaj^2) =  plam + p2lam2 + p3lam3 + p4lam4
{
value=tryps[p]/1000+tryp2s[p]/100+tryp3s[p]/10+tryp4s[p];
lambdas[p]=0.001/value;
lambdas2[p]=0.01/value;
lambdas3[p]=0.1/value;
lambdas4[p]=1.0/value;
}

if(trytypes[p]==7)	//elastic - exp(betaj^2) = 2p/lam^2 + 2p2/lam2^2 + p3lam3 (alt 2(p+p2)/lam^2 + p3lam3)
{
if(tryp3s[p]==0)	//lasso (will have p=p2=1/2)
{lambdas[p]=1.0;lambdas2[p]=1.0;}
if(tryp3s[p]==1)	//ridge
{lambdas3[p]=1.0;}
if(tryp3s[p]>0||tryp3s[p]<1)	//elastic
{
lambdas[p]=pow(4*tryps[p]/(1-tryf2s[p]),.5);
lambdas2[p]=pow(4*tryp2s[p]/(1-tryf2s[p]),.5);
lambdas3[p]=tryf2s[p]/tryp3s[p];
}
}

if(trytypes[p]==8)	//ldpred - exp(betaj^2) = plam (special case of bolt)
{
lambdas[p]=pow(tryps[p],-1);
lambdas2[p]=0;
}
}	//end of p loop

///////////////////////////

