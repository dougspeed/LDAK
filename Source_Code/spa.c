/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Functions for performing SPA

///////////////////////////

int spa_logistic_one(double score, double *data, int ns, double *probs, double *weights, double *stats, double scal2)
{
int i, count, flag;
double value, value2, value3, diff, relax;

double kk, kkb, k0, k1, k1b, k2, ww, vv, ss, ssold;


//K0(t) is sum ( log (probs exp(tX) + (1-probs)) ) - t XT probs
//K1(t) is sum ( X probs / (probs + (1-probs) exp(-tX)) ) - XT probs
//K2(t) is sum ( X^2 probs (1-probs) exp(-tX) / (probs + (1-probs) exp(-tX))^2 )

//useful to compute sum (data probs)
value3=0;for(i=0;i<ns;i++){value3+=data[i]*probs[i];}

//flag will indicate success: 1=good, 2=approx (timeout), -1=fail

//starting knot is zero
kk=0;

//work out starting cgfs (k0 and k1 are trivial)
k0=0;k1=0;k2=0;
for(i=0;i<ns;i++){k2+=pow(data[i],2)*weights[i];}

//starting stat is non-spa stat
ssold=score/k2;

count=0;
while(1)
{
if(count==10)
{
printf("Warning, SPA did not converge within %d iterations\n\n", count);
if(flag==1){flag=2;}
break;
}

//work out full move
diff=(score-k1)/k2;

flag=-1;
relax=1;
while(relax>0.001)
{
//get proposed move and corresponding k1
kkb=kk+relax*diff;
k1b=-value3;
for(i=0;i<ns;i++)
{k1b+=data[i]*probs[i]/(probs[i]+(1-probs[i])*exp(-data[i]*kkb));}

if(fabs(k1b-score)<fabs(k1-score))	//good move - update knot and cgfs and break
{
kk=kkb;
k1=k1b;
k0=-kk*value3;k2=0;
for(i=0;i<ns;i++)
{
value=exp(-data[i]*kk);
value2=probs[i]+(1-probs[i])*value;
k0+=log(value2/value);
k2+=pow(data[i],2)*weights[i]*value*pow(value2,-2);
flag=1;
}
break;
}
else	//bad move - will try smaller move
{relax*=0.5;}
}

if(flag==-1)	//failed to move
{break;}

if(kk*score-k0>0&&k2>0)	//can compute the spa statistic and maybe break
{
ww=pow(2*(kk*score-k0),.5);
if(kk<0){ww=-ww;}
vv=kk*pow(k2,.5);
ss=(ww+pow(ww,-1)*log(vv/ww));

if(fabs(ss-ssold)<0.01){break;}

ssold=ss;
}
else	//failed this time - hopefully will work later
{flag=-1;}

count++;
}

if(flag==1||flag==2)	//load up stats (remembering to scale)
{
stats[2]=ss*scal2;
stats[1]=stats[0]/stats[2];
stats[3]=erfc(fabs(stats[2])*M_SQRT1_2);
}

return(flag);
}

////////

int spa_logistic_two(double score2, double XTCX, double *data, int ns, double *probs, double *weights, double *stats, double scal2)
{
int i, count, flag;
double score, value, value2, value3;

double kk, k0, k1, k2, ww, vv, ss, ss2, ssold, ssave, pva;


//K0(t) is sum ( log (probs exp(tX) + (1-probs)) )
//K1 is sum ( X probs / (probs + (1-probs) exp(-tX)) )
//K2 is sum ( X^2 probs (1-probs) exp(-tX) / (probs + (1-probs) exp(-tX))^2 )

//useful to compute sum (data probs)
value3=0;for(i=0;i<ns;i++){value3+=data[i]*probs[i];}

//flag will indicate success: 1=good, -1=fail (ignore timeout here)

//first solve for positive score
score=fabs(score2);

//starting knot is based on taylors expantion of K1 around zero: K1(t) = sum (X probs) + t sum (X^2 probs (1-probs)) = score2
kk=score/XTCX;

count=0;
while(1)
{
if(count==10){printf("Warning, SPA did not converge within %d iterations\n\n", count);break;}

k0=-kk*value3;k1=-value3;k2=0;
for(i=0;i<ns;i++)
{
value=exp(-data[i]*kk);
value2=probs[i]+(1-probs[i])*value;
k0+=log(value2/value);
k1+=data[i]*probs[i]/value2;
k2+=pow(data[i],2)*weights[i]*value*pow(value2,-2);
}

if(kk*score-k0>0&&k2>0)	//can compute the spa statistic
{
ww=pow(2*(kk*score-k0),.5);
if(kk<0){ww=-ww;}
vv=kk*pow(k2,.5);
ss=(ww+pow(ww,-1)*log(vv/ww));

if(count>0&&flag==1)	//the previous iteration worked
{
if(fabs(ss-ssold)<0.01){break;}
}
ssold=ss;
flag=1;
}
else	//spa has failed this iteration - hopefully it works later
{
flag=-1;
}

//update kk
kk-=(k1-score)/k2;

count++;
}

if(flag==1)	//now solve for the negative score
{
score=-fabs(score2);

//starting knot is based on taylors expantion of K1 around zero: K1(t) = sum (X probs) + t sum (X^2 probs (1-probs)) = score2
kk=score/XTCX;

count=0;
while(1)
{
if(count==10){printf("Warning, SPA did not converge within %d iterations\n\n", count);break;}

k0=-kk*value3;k1=-value3;k2=0;
for(i=0;i<ns;i++)
{
value=exp(-data[i]*kk);
value2=probs[i]+(1-probs[i])*value;
k0+=log(value2/value);
k1+=data[i]*probs[i]/value2;
k2+=pow(data[i],2)*weights[i]*value*pow(value2,-2);
}

if(kk*score-k0>0&&k2>0)	//can compute the spa statistic
{
ww=pow(2*(kk*score-k0),.5);
if(kk<0){ww=-ww;}
vv=kk*pow(k2,.5);
ss2=(ww+pow(ww,-1)*log(vv/ww));

if(count>0&&flag==1)	//the previous iteration worked
{
if(fabs(ss2-ssold)<0.01){break;}
}
ssold=ss2;
flag=1;
}
else	//will use the non-spa statistic
{
flag=-1;
break;
}

//update kk
kk-=(k1-score)/k2;

count++;
}
}

if(flag==1)	//get the average and load up stats (remembering to account for sign and scale)
{
pva=.5*(erfc(fabs(ss)*M_SQRT1_2)+erfc(fabs(ss2)*M_SQRT1_2));
ssave=normal_inv(1-.5*pva);

//check ssave between ss and ss2
if(fabs(ss)>fabs(ss2))	//ss biggest
{
if(fabs(ssave)<fabs(ss2)||fabs(ssave)>fabs(ss)){ssave=fabs(ss2);
printf("warn 1 score %f flag %d positive %f negative %f average %f orig %f pva %f\n", score2, flag, ss, ss2, ssave, stats[2]/scal2, pva);
}
}
else	//ss2 biggest
{
if(fabs(ssave)<fabs(ss)||fabs(ssave)>fabs(ss2)){ssave=fabs(ss);
printf("warn 2 score %f flag %d positive %f negative %f average %f orig %f pva %f\n", score2, flag, ss, ss2, ssave, stats[2]/scal2, pva);
}
}

if(score2>0){stats[2]=ssave*scal2;}
else{stats[2]=-ssave*scal2;}
stats[1]=stats[0]/stats[2];
stats[3]=erfc(fabs(stats[2])*M_SQRT1_2);
}

return(flag);
} 

///////////////////////////

void exact_cumulants(double *data, double centre, double *Yadj, int ns, int num_knots, double *knots)
{
int i, j, k;

int num_bins, *map;
double min, max, scale, sum, sum2, sum3, mean, mean2, mean3, value, value2;

double *bins, *K0, *K1, *K2, **CGF0, **CGF1, **CGF2;


//get max
max=fabs(data[0]);
for(i=1;i<ns;i++)
{
if(fabs(data[i])>max){max=fabs(data[i]);}
}
scale=2.0/max;

num_bins=4;

bins=malloc(sizeof(double)*4);
bins[0]=-centre*scale;
bins[1]=(1-centre)*scale;
bins[2]=(2-centre)*scale;
bins[3]=0;

map=malloc(sizeof(int)*ns);

for(i=0;i<ns;i++)
{
map[i]=-1;
if(data[i]==-centre){map[i]=0;}
if(data[i]==1-centre){map[i]=1;}
if(data[i]==2-centre){map[i]=2;}
if(data[i]==0){map[i]=3;}
if(map[i]==-1){printf("Error %d value %f centres %f\n", i, data[i], centre);exit(1);}
}

CGF0=malloc(sizeof(double*)*num_knots);
CGF1=malloc(sizeof(double*)*num_knots);
CGF2=malloc(sizeof(double*)*num_knots);

for(j=0;j<num_knots;j++)
{
CGF0[j]=malloc(sizeof(double)*num_bins);
CGF1[j]=malloc(sizeof(double)*num_bins);
CGF2[j]=malloc(sizeof(double)*num_bins);
}

//get min and max phenotype
min=Yadj[0];
max=Yadj[0];
for(i=1;i<ns;i++)
{
if(Yadj[i]>max){max=Yadj[i];}
if(Yadj[i]<min){min=Yadj[i];}
}

//compute cumulants for knots x bins
#pragma omp parallel for private(j, k, i, value, value2, sum, sum2, sum3, mean, mean2, mean3) schedule(static)
for(j=0;j<num_knots;j++)
{
for(k=0;k<num_bins;k++)
{
//what is max value of bins * knots * Yadj
if(bins[k]*knots[j]>0){value=max*bins[k]*knots[j];}
else{value=min*bins[k]*knots[j];}
sum=0;
sum2=0;
sum3=0;
for(i=0;i<ns;i++)
{
value2=exp(bins[k]*knots[j]*Yadj[i]-value);
sum+=value2;
sum2+=Yadj[i]*value2;
sum3+=pow(Yadj[i],2)*value2;
}
mean=sum/ns;
mean2=sum2/ns;
mean3=sum3/ns;

CGF0[j][k]=log(mean)+value;
CGF1[j][k]=mean2/mean;
CGF2[j][k]=mean3/mean-pow(mean2/mean,2);

if(CGF0[j][k]!=CGF0[j][k]||CGF1[j][k]!=CGF1[j][k]||CGF2[j][k]!=CGF2[j][k])
{printf("Error for knot %d %f, bin %f means %f %f %f\n", j+1, knots[j], bins[k], mean, mean2, mean3);exit(1);}
}
}

K0=malloc(sizeof(double)*num_knots);
K1=malloc(sizeof(double)*num_knots);
K2=malloc(sizeof(double)*num_knots);

for(j=0;j<num_knots;j++)
{
K0[j]=0;
K1[j]=0;
K2[j]=0;
for(i=0;i<ns;i++)
{
k=map[i];
K0[j]+=CGF0[j][k];
K1[j]+=bins[k]*CGF1[j][k];
K2[j]+=pow(bins[k],2)*CGF2[j][k];
}
}

FILE *test=fopen("true.txt","w");
for(j=0;j<num_knots;j++)
{fprintf(test, "%f %f %f %f\n", knots[j], K0[j], K1[j], K2[j]);}
fclose(test);

for(j=0;j<num_knots;j++){free(CGF0[j]);free(CGF1[j]);free(CGF2[j]);}
free(bins);free(map);free(CGF0);free(CGF1);free(CGF2);free(K0);free(K1);free(K2);

}

///////////////////////////

int spa_test(double score2, double *data, int ns, int num_knots, double *knots, int num_bins, double *bins, double **CGF0, double **CGF1, double **CGF2, double **CGF3, double *stats, double scal2)
//uses taylors expansion for bins, interpolation for knots
{
int i, j, k, mark;
double value, value2, value3, max, diff, scale;

double score, *acs, *acs2, *acs3, *dcs, *dcs2, *dcs3, *ecs, *ecs2, *fcs;

int fleft, fright, fmiddle;
double kk, k0, k2, *k1s, ww, vv;


acs=malloc(sizeof(double)*num_bins);
acs2=malloc(sizeof(double)*num_bins);
acs3=malloc(sizeof(double)*num_bins);
dcs=malloc(sizeof(double)*num_bins);
dcs2=malloc(sizeof(double)*num_bins);
dcs3=malloc(sizeof(double)*num_bins);
ecs=malloc(sizeof(double)*num_bins);
ecs2=malloc(sizeof(double)*num_bins);
fcs=malloc(sizeof(double)*num_bins);
k1s=malloc(sizeof(double)*num_knots);

//find max (absolute) genotype - will then scale all genotypes by 2/max
max=fabs(data[0]);
for(i=1;i<ns;i++)
{
if(fabs(data[i])>max){max=fabs(data[i]);}
}

//scale score by 2/max
score=score2*2/max;

//for each individual, find mark, such that bins[mark] is closest value to 2/max*data[i]
value=2.0+2*pow(num_bins-1,-1);
value2=0.25*(num_bins-1);
for(k=0;k<num_bins;k++)
{acs[k]=0;acs2[k]=0;acs3[k]=0;dcs[k]=0;dcs2[k]=0;dcs3[k]=0;ecs[k]=0;ecs2[k]=0;fcs[k]=0;}

for(i=0;i<ns;i++)
{
value3=2*data[i]/max;
mark=(int)(value2*(value3+value));
if(mark<0){mark=0;}
if(mark>num_bins-1){mark=num_bins-1;}

diff=value3-bins[mark];
acs[mark]++;
acs2[mark]+=value3;
acs3[mark]+=pow(value3,2);
dcs[mark]+=diff;
dcs2[mark]+=value3*diff;
dcs3[mark]+=pow(value3,2)*diff;
ecs[mark]+=pow(diff,2);
ecs2[mark]+=value3*pow(diff,2);
fcs[mark]+=pow(diff,3);
}

//for(k=0;k<num_bins;k++){printf("bin %f has %f %f %f and %f %f %f\n", bins[k], acs[k], acs2[k],acs3[k], dcs[k], dcs2[k], dcs3[k]);}

//compute first cumulants for all knots (bit redundant for final knot)
for(j=0;j<num_knots;j++)
{
k1s[j]=0;for(k=0;k<num_bins;k++){k1s[j]+=acs2[k]*CGF1[j][k]+knots[j]*dcs2[k]*CGF2[j][k]+pow(knots[j],2)/2*ecs2[k]*CGF3[j][k];}

//k0=0;for(k=0;k<num_bins;k++){k0+=acs[k]*CGF0[j][k]+knots[j]*dcs[k]*CGF1[j][k]+pow(knots[j],2)/2*ecs[k]*CGF2[j][k]+pow(knots[j],3)/6*fcs[k]*CGF3[j][k];}
//k2=0;for(k=0;k<num_bins;k++){k2+=acs3[k]*CGF2[j][k]+knots[j]*dcs3[k]*CGF3[j][k];}
//printf("Count %d knot %f k1 %f k0 %f k2 %f score %f\n", j, knots[j], k0, k1s[j], k2, score);
}

//find the points where first cumulant stops being monotonic
fleft=(num_knots-2)/2;
while(fleft>0)
{
if(k1s[fleft-1]>k1s[fleft]){break;}
fleft--;
}
fright=(num_knots-1)/2;
while(fright<num_knots-2)
{
if(k1s[fright+1]<k1s[fright]){break;}
fright++;
}
//if(fleft>0||fright<num_bins-2){printf("Warning, first cumulant is only monotonic from %d to %d\n", fleft+1, fright+1);}

if(score>=k1s[fleft]&&score<=k1s[fright])	//score within range - find bounding knots using divide and conquer
{
while(fright-fleft>1)
{
fmiddle=(fleft+fright)/2;
if(score<=k1s[fmiddle])	//root must be between fleft and fmiddle
{fright=fmiddle;}
else	//root must be between fmiddle and fright
{fleft=fmiddle;}
}
}
else	//score outside range - will either extend range or use just one knot
{
printf("Warning, score is %f but extreme values of (monotonic) k1 are %f and %f (knots %f %f)\n", score, k1s[fleft], k1s[fright], knots[fleft], knots[fright]);

if(fleft==0&&score/k1s[0]>1.01&&k1s[0]/k1s[1]>1.01)	//appears it will help to increase left range
{
free(acs);free(acs2);free(acs3);free(dcs);free(dcs2);free(dcs3);free(ecs);free(ecs2);free(fcs);free(k1s);
return(-2);
}
else	//has failed (might later change so use first knot)
{
stats[2]=0;
stats[1]=-9999;
stats[3]=1;
free(acs);free(acs2);free(acs3);free(dcs);free(dcs2);free(dcs3);free(ecs);free(ecs2);free(fcs);free(k1s);
return(-1);
//fright=fleft;
}

if(fright==num_bins-2&&score/k1s[num_knots-2]>1.01&&k1s[num_knots-3]/k1s[num_knots-2]>1.01)	//appears it will help to increase right range
{
free(acs);free(acs2);free(acs3);free(dcs);free(dcs2);free(dcs3);free(ecs);free(ecs2);free(fcs);free(k1s);
return(-2);
}
else	//has failed (might later change so use final knot)
{
stats[2]=0;
stats[1]=-9999;
stats[3]=1;
free(acs);free(acs2);free(acs3);free(dcs);free(dcs2);free(dcs3);free(ecs);free(ecs2);free(fcs);free(k1s);
return(-1);
//fleft=fright;
}
}

//printf("knots are %f and %f, score is %f k1s are %f %f\n", knots[fleft], knots[fright], score, k1s[fleft], k1s[fright]);

if(knots[fleft]<=0&&knots[fright]>=0)	//knots either side of zero, so compute non-spa statistic
{
k2=0;for(k=0;k<num_bins;k++){k2+=acs3[k]*CGF2[num_knots-1][k];}

if(k2>0)
{
stats[2]=score*pow(k2,-.5)*scal2;
stats[1]=stats[0]/stats[2];
stats[3]=erfc(fabs(stats[2])*M_SQRT1_2);
free(acs);free(acs2);free(acs3);free(dcs);free(dcs2);free(dcs3);free(ecs);free(ecs2);free(fcs);free(k1s);
return(2);
}
else	//has failed, and nothing we can do
{
stats[2]=0;
stats[1]=-9999;
stats[3]=1;
free(acs);free(acs2);free(acs3);free(dcs);free(dcs2);free(dcs3);free(ecs);free(ecs2);free(fcs);free(k1s);
return(-1);
}
}
else	//use spa
{
//scale indicates fraction between left and right knots (0 => left knot, 1 => right knot)
if(fleft==fright){scale=0;}
else{scale=(k1s[fleft]-score)/(k1s[fleft]-k1s[fright]);}

//printf("k1s are %f %f %f scale %f\n", k1s[fleft], score, k1s[fright], scale);

//compute location of final knot
kk=(1-scale)*knots[fleft]+scale*knots[fright];

//compute k0 and k2 for final knot
value=0;for(k=0;k<num_bins;k++){value+=acs[k]*CGF0[fleft][k]+knots[fleft]*dcs[k]*CGF1[fleft][k]+pow(knots[fleft],2)/2*ecs[k]*CGF2[fleft][k]+pow(knots[fleft],3)/6*fcs[k]*CGF3[fleft][k];}
value2=0;for(k=0;k<num_bins;k++){value2+=acs[k]*CGF0[fright][k]+knots[fright]*dcs[k]*CGF1[fright][k]+pow(knots[fright],2)/2*ecs[k]*CGF2[fright][k]+pow(knots[fright],3)/6*fcs[k]*CGF3[fright][k];}
k0=(1-scale)*value+scale*value2;

value=0;for(k=0;k<num_bins;k++){value+=acs3[k]*CGF2[fleft][k]+knots[fleft]*dcs3[k]*CGF3[fleft][k];}
value2=0;for(k=0;k<num_bins;k++){value2+=acs3[k]*CGF2[fright][k]+knots[fright]*dcs3[k]*CGF3[fright][k];}
k2=(1-scale)*value+scale*value2;

//printf("values %f %f and %f from %f %f w %f v %f\n", kk*score, k0, k2, value, value2, pow(2*(kk*score-k0),.5), kk*pow(k2,.5));

if(kk*score-k0>0&&k2>0)	//can compute the spa statistic
{
ww=pow(2*(kk*score-k0),.5);
if(kk<0){ww=-ww;}
vv=kk*pow(k2,.5);

stats[2]=(ww+pow(ww,-1)*log(vv/ww))*scal2;
stats[1]=stats[0]/stats[2];
stats[3]=erfc(fabs(stats[2])*M_SQRT1_2);

free(acs);free(acs2);free(acs3);free(dcs);free(dcs2);free(dcs3);free(ecs);free(ecs2);free(fcs);free(k1s);
return(1);
}
else	//have failed, so try to compute non-spa statistic
{
k2=0;for(k=0;k<num_bins;k++){k2+=acs3[k]*CGF2[num_knots-1][k];}

if(k2>0)
{
stats[2]=score*pow(k2,-.5)*scal2;
stats[1]=stats[0]/stats[2];
stats[3]=erfc(fabs(stats[2])*M_SQRT1_2);
free(acs);free(acs2);free(acs3);free(dcs);free(dcs2);free(dcs3);free(ecs);free(ecs2);free(fcs);free(k1s);
return(2);
}
else	//has failed, and nothing we can do
{
stats[2]=0;
stats[1]=-9999;
stats[3]=1;
free(acs);free(acs2);free(acs3);free(dcs);free(dcs2);free(dcs3);free(ecs);free(ecs2);free(fcs);free(k1s);
return(-1);
}
}
}	//end of using spa

printf("Error, should not be here\n\n");exit(1);
}

////////

int spa_test2(double score2, double *data, int ns, int num_knots, double *knots, int num_bins, double *bins, double **CGF0, double **CGF1, double **CGF2, double **CGF3, double *stats, int type)
//type 0 - will scale genotypes and score, type 1 - will not
{
int i, k, mark;
double value, value2, value3, max, scale;

double score, *acs, *acs2, *acs3;

int fleft, fright, fmiddle;
double k1left, k1right, k1middle, kk, k0, k2, ww, vv;


if(type==0)	//find max (absolute) genotype - will then scale all genotypes by 2/max
{
max=fabs(data[0]);
for(i=1;i<ns;i++)
{
if(fabs(data[i])>max){max=fabs(data[i]);}
}
}
else{max=2;}

//scale score by 2/max
score=score2*2/max;

//compute bin counts
acs=malloc(sizeof(double)*num_bins);
acs2=malloc(sizeof(double)*num_bins);
acs3=malloc(sizeof(double)*num_bins);
for(k=0;k<num_bins;k++){acs[k]=0;acs2[k]=0;acs3[k]=0;}

//for each individual, find mark, such that 2/max*data[i] is between bins[mark] and bins[mark+1]
//then find scale, which indicates fraction between left and right bin (0 => left bin, 1 => right bin)
value=.5*(num_bins-1);
value2=.5*(num_bins-1)/max;

for(i=0;i<ns;i++)
{
mark=(int)(value+value2*data[i]);
if(mark<0){mark=0;}
if(mark>=num_bins-1){mark=num_bins-2;}
scale=value2*data[i]-.5*value*bins[mark];

value3=2*data[i]/max;
acs[mark]+=1-scale;
acs[mark+1]+=scale;
acs2[mark]+=(1-scale)*value3;
acs2[mark+1]+=scale*value3;
acs3[mark]+=(1-scale)*pow(value3,2);
acs3[mark+1]+=scale*pow(value3,2);
}

//for(k=0;k<num_bins;k++){printf("bin %f has %f %f %f\n", bins[k], acs[k], acs2[k],acs3[k]);}

//first compute k1 for left and right knots, to confirm these bound score (ignore last knot, which is zero)
fleft=0;
k1left=0;for(k=0;k<num_bins;k++){k1left+=acs2[k]*CGF1[fleft][k];}
fright=num_knots-2;
k1right=0;for(k=0;k<num_bins;k++){k1right+=acs2[k]*CGF1[fright][k];}

//if(score<k1left||score>k1right){printf("Warning, score is %f but extreme values of k1 are %f and %f knots %d %d are %f %f\n", score, k1left, k1right, fleft, fright, knots[fleft], knots[fright]);}

if(score<k1left)	//see if it will help to increase left range
{
fmiddle=1;
k1middle=0;for(k=0;k<num_bins;k++){k1middle+=acs2[k]*CGF1[fmiddle][k];}

if(k1left-score>0.01&&k1left/k1middle<0.99)	//appears it will help to increase range
{
stats[1]=-1;
free(acs);free(acs2);free(acs3);
return(0);
}
else	//will compute using left knot
{fright=fleft;k1right=k1left;}
}

if(score>k1right)	//see if it will help to increase right range
{
fmiddle=num_knots-3;
k1middle=0;for(k=0;k<num_bins;k++){k1middle+=acs2[k]*CGF1[fmiddle][k];}

if(score-k1right>0.01&&k1middle/k1right<0.99)	//appears it will help to increase range
{
stats[1]=-1;
free(acs);free(acs2);free(acs3);
return(0);
}
else	//will compute using last knot
{fleft=fright;k1left=k1right;}
}

while(fright-fleft>1)	//find closest knots using divide and conquer
{
fmiddle=(fleft+fright)/2;
k1middle=0;for(k=0;k<num_bins;k++){k1middle+=acs2[k]*CGF1[fmiddle][k];}

if(score<=k1middle)	//root must be between fleft and fmiddle
{fright=fmiddle;k1right=k1middle;}
else	//root must be between fmiddle and fright
{fleft=fmiddle;k1left=k1middle;}
}

//printf("knots are %f and %f, score is %f k1s are %f %f\n", knots[fleft], knots[fright], score, k1left, k1right);

if(knots[fleft]<=0&&knots[fright]>=0)	//knots either side of zero - use standard approximation
{
k2=0;for(k=0;k<num_bins;k++){k2+=acs3[k]*CGF2[num_knots-1][k];}

if(k2>0)
{
stats[2]=score*pow(k2,-.5);
stats[1]=stats[0]/stats[2];
stats[3]=erfc(fabs(stats[2])*M_SQRT1_2);
}
else	//has failed, and nothing we can do
{
stats[2]=0;
stats[1]=-9999;
stats[3]=1;
}
}
else	//use spa
{
//scale indicates fraction between left and right knots (0 => left knot, 1 => right knot)
scale=(k1left-score)/(k1left-k1right);

//printf("k1s are %f %f %f scale %f\n", k1left, score, k1right, scale);

//compute location of final knot
kk=(1-scale)*knots[fleft]+scale*knots[fright];

//compute k0 and k2 for final knot
value=0;for(k=0;k<num_bins;k++){value+=acs[k]*CGF0[fleft][k];}
value2=0;for(k=0;k<num_bins;k++){value2+=acs[k]*CGF0[fright][k];}
k0=(1-scale)*value+scale*value2;

value=0;for(k=0;k<num_bins;k++){value+=acs3[k]*CGF2[fleft][k];}
value2=0;for(k=0;k<num_bins;k++){value2+=acs3[k]*CGF2[fright][k];}
k2=(1-scale)*value+scale*value2;

//printf(" %f %f and %f from %f %f w %f v %f\n", kk*score, k0, k2, value, value2, pow(2*(kk*score-k0),.5), kk*pow(k2,.5));

if(kk*score-k0>0&&k2>0)	//can compute the square root
{
ww=pow(2*(kk*score-k0),.5);
if(kk<0){ww=-ww;}
vv=kk*pow(k2,.5);
stats[2]=ww+pow(ww,-1)*log(vv/ww);
stats[1]=stats[0]/stats[2];
stats[3]=erfc(fabs(stats[2])*M_SQRT1_2);
}
else	//have failed, and nothing we can do
{
stats[2]=0;
stats[1]=-9999;
stats[3]=1;
}
}	//end of using spa

free(acs);free(acs2);free(acs3);

return(0);
}

///////////////////////////

void empirical_cumulants(double *Yadj, int ns, int num_knots, double *knots, int num_bins, double *bins, double **CGF0, double **CGF1, double **CGF2, double **CGF3, double spamax)
{
int i, j, k;
double min, max, sum, sum2, sum3, sum4, mean, mean2, mean3, mean4, value, value2;


//get min and max phenotype
min=Yadj[0];
max=Yadj[0];
for(i=1;i<ns;i++)
{
if(Yadj[i]>max){max=Yadj[i];}
if(Yadj[i]<min){min=Yadj[i];}
}

//set basic knots - making sure final one is zero
value=pow(num_knots,-1);
//for(j=0;j<num_knots-1;j++){knots[j]=normal_inv(value*(j+1));}
for(j=0;j<num_knots-1;j++){knots[j]=tan(M_PI*(value*(j+1)-.5));}
knots[num_knots-1]=0;

//ensure size of largest knots is spamax
value=-spamax/knots[0];
for(j=0;j<num_knots;j++){knots[j]*=value;}

//compute cumulants for knots x bins
#pragma omp parallel for private(j, k, i, value, value2, sum, sum2, sum3, sum4, mean, mean2, mean3, mean4) schedule(static)
for(j=0;j<num_knots;j++)
{
for(k=0;k<num_bins;k++)
{
//what is max value of bins * knots * Yadj
if(bins[k]*knots[j]>0){value=max*bins[k]*knots[j];}
else{value=min*bins[k]*knots[j];}
sum=0;
sum2=0;
sum3=0;
sum4=0;
for(i=0;i<ns;i++)
{
value2=exp(bins[k]*knots[j]*Yadj[i]-value);
sum+=value2;
sum2+=Yadj[i]*value2;
sum3+=pow(Yadj[i],2)*value2;
sum4+=pow(Yadj[i],3)*value2;
}
mean=sum/ns;
mean2=sum2/ns;
mean3=sum3/ns;
mean4=sum4/ns;

CGF0[j][k]=log(mean)+value;
CGF1[j][k]=mean2/mean;
CGF2[j][k]=mean3/mean-pow(mean2/mean,2);
CGF3[j][k]=mean4/mean-3*mean2*mean3*pow(mean,-2)+2*pow(mean2/mean,3);

if(CGF0[j][k]!=CGF0[j][k]||CGF1[j][k]!=CGF1[j][k]||CGF2[j][k]!=CGF2[j][k])
{printf("Error for knot %d %f, bin %f means %f %f %f\n", j+1, knots[j], bins[k], mean, mean2, mean3);exit(1);}
}
}
}

////////

/*
void logistic_cumulants(double *nullprobs, int ns, int num_knots, double *knots, int num_bins, double *bins, double **CGF0, double **CGF1, double **CGF2, double **CGF3, double spamax)
{
int i, j, k;
double min, max, sum, sum2, sum3, sum4, mean, mean2, mean3, mean4, value, value2;


//set basic knots - making sure final one is zero
value=pow(num_knots,-1);
//for(j=0;j<num_knots-1;j++){knots[j]=normal_inv(value*(j+1));}
for(j=0;j<num_knots-1;j++){knots[j]=tan(M_PI*(value*(j+1)-.5));}
knots[num_knots-1]=0;

//ensure size of largest knots is spamax
value=-spamax/knots[0];
for(j=0;j<num_knots;j++){knots[j]*=value;}

//compute cumulants for knots x bins
#pragma omp parallel for private(j, k, i, value, value2, sum, sum2, sum3, sum4, mean, mean2, mean3, mean4) schedule(static)
for(j=0;j<num_knots;j++)
{
for(k=0;k<num_bins;k++)
{




//what is max value of bins * knots * Yadj
if(bins[k]*knots[j]>0){value=max*bins[k]*knots[j];}
else{value=min*bins[k]*knots[j];}
sum=0;
sum2=0;
sum3=0;
sum4=0;
for(i=0;i<ns;i++)
{
value2=exp(bins[k]*knots[j]*Yadj[i]-value);
sum+=value2;
sum2+=Yadj[i]*value2;
sum3+=pow(Yadj[i],2)*value2;
sum4+=pow(Yadj[i],3)*value2;
}
mean=sum/ns;
mean2=sum2/ns;
mean3=sum3/ns;
mean4=sum4/ns;

CGF0[j][k]=log(mean)+value;
CGF1[j][k]=mean2/mean;
CGF2[j][k]=mean3/mean-pow(mean2/mean,2);
CGF3[j][k]=mean4/mean-3*mean2*mean3*pow(mean,-2)+2*pow(mean2/mean,3);

if(CGF0[j][k]!=CGF0[j][k]||CGF1[j][k]!=CGF1[j][k]||CGF2[j][k]!=CGF2[j][k])
{printf("Error for knot %d %f, bin %f means %f %f %f\n", j+1, knots[j], bins[k], mean, mean2, mean3);exit(1);}
}
}
}
*/

///////////////////////////

