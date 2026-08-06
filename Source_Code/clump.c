/*
Copyright 2026 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//use conditional analysis to get Lbit (must be at least one) and load effs3

///////////////////////////

//get significance threshold
value5=pow(normal_inv(cutoff/2),2);

//use indexer to store top hits
indexer=malloc(sizeof(int)*50);

for(count=0;count<50;count++)
{
//work out which predictor to pick next
if(count==0)    //find best predictor based on raw significance (will always use)
{
best=0;
value=chis[bitstart];
value2=rhos[bitstart];
for(j=1;j<bitlength;j++)
{
if(chis[bitstart+j]>value){best=j;value=chis[bitstart+j];value2=rhos[bitstart+j];}
}
}
else    //find best predictor based on conditional significance (will only use if significant)
{
Z=malloc(sizeof(double)*bitlength);
Z2=malloc(sizeof(double)*bitlength);
ZTZ=malloc(sizeof(double)*count*count);
ZTZ2=malloc(sizeof(double)*count);
ZTdata=malloc(sizeof(double)*count*(bitlength+1));

//load up correlations for picked predictors - could reverse shrink
for(j2=0;j2<count;j2++)
{
for(j=0;j<count;j++)
{
if(indexer[j2]>=indexer[j]){ZTZ[j+j2*count]=cors[(size_t)indexer[j]*bitlength+indexer[j2]];}
else{ZTZ[j+j2*count]=cors[(size_t)indexer[j2]*bitlength+indexer[j]];}
}
}

//load up cross-correlations - again could reverse shrink
for(j2=0;j2<bitlength;j2++)
{
for(j=0;j<count;j++)
{
if(j2>=indexer[j]){ZTdata[j+j2*count]=cors[(size_t)indexer[j]*bitlength+j2];}
else{ZTdata[j+j2*count]=cors[(size_t)j2*bitlength+indexer[j]];}
}
}

//final column of ZTdata contains rhos for included predictors
for(j=0;j<count;j++){ZTdata[j+bitlength*count]=rhos[bitstart+indexer[j]];}

//get inv(ZTZ) ZTdata
eigen_invert(ZTZ, count, ZTZ2, bitlength+1, ZTdata, 1);

//residual variance is 1 - rhos_picked ZTdata_bitlength
value3=1.0;
for(j=0;j<count;j++){value3-=rhos[bitstart+indexer[j]]*ZTdata[j+bitlength*count];}

for(j2=0;j2<bitlength;j2++)
{
//top is rhoj2 - rho_picked ZTdata_j2
value=rhos[bitstart+j2];
for(j=0;j<count;j++){value-=rhos[bitstart+indexer[j]]*ZTdata[j+j2*count];}

//bottom is 1 - cor_picked ZTdata_j2 - once again could reverse shrink
value2=1.0;
for(j=0;j<count;j++)
{
if(j2>=indexer[j]){value2-=cors[(size_t)indexer[j]*bitlength+j2]*ZTdata[j+j2*count];}
else{value2-=cors[(size_t)j2*bitlength+indexer[j]]*ZTdata[j+j2*count];}
}

//store test statistic and effect size
if(value2>0.1&&value3>0.1){Z[j2]=nss[bitstart+j2]*pow(value,2)/value2/value3/(1-pow(value,2)/value2/value3);Z2[j2]=value/value2;}
else{Z[j2]=0;Z2[j2]=0;}
}

//check test statistics zero for already included predictors (should be because value2=0)
for(j=0;j<count;j++)
{
if(Z[indexer[j]]>0){printf("not work bit %d count %d pred %d; please tell Doug\n", bit+1, count+1, bitstart+indexer[j]);exit(1);}
}

//get most significant predictor and corresponding effect size
best=0;
value=Z[0];
value2=Z2[0];
for(j=1;j<bitlength;j++)
{
if(Z[j]>value){best=j;value=Z[j];value2=Z2[j];}
}

free(Z);free(Z2);free(ZTZ);free(ZTZ2);free(ZTdata);

//will stop if not significant
if(value<value5){break;}
}

//so the most significant predictor is best, which has test statistic value and effect size value2
fprintf(output3,"%s %d %.0f %.4e %.4e\n", preds[bitstart+best], chr[bitstart+best], bp[bitstart+best], erfc(pow(value,.5)*M_SQRT1_2), erfc(pow(chis[bitstart+best],.5)*M_SQRT1_2));

indexer[count]=best;
}

if(count==50){printf("Error, Window %d contains at least 50 significant predictors (P<%.4e); if you see this message, please can you tell Doug\n\n", bit+1, cutoff);} 

free(indexer);

///////////////////////////

