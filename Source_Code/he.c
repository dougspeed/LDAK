/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Haseman Elston and PCGC regression - will always have at least one subset

///////////////////////////

void he_reg(int ns, int num_subs, int num_covars, int num_envs, int num_tops, int num_kins, int num_regs, int **subindex, double *Y, double *Z, float **Mkins_single, double *kintraces, double *kinsums, double *X, int Xtotal, int *Xstarts, int *Xends, double *Xsums, double prev, int np, int discenv, char *oversfile, double **ssums, int num_blocks, int memsave, char **kinstems, char **ids3, char *outfile, int type, double trun, double missingvalue)
//type=0 - he, type=1 - pcgc
{
size_t scounta, scountb, stotal, smark;
int i, i2, j, k, k2, k3, p, p2, r, s, count, kcount, start, end, flag;
double sum, sum2, sumsq, mean, mean2, var, value, factor;
float *datatemp;
char **wantids;

int num_fixed, total, *indexer, *indexer2, *subs, *block;
double *thetas, *thetasds, *thetapvas, *Yadj, covher, topher;
double *hers, *hersds, *shares, *sharesds, *cohers, *cohers2, *jacks, stats[4];
double likenull, like, lrtstat, lrtpva;

double Ysumsq, *Dtemp, *DTD, *DTS, *DTDs, *DTSs, *DTD2, **MDTD, **MDTS;
double Ysumsqa, Ysumsqb, *DTDa, *DTDb, *DTSa, *DTSb, **MDTDa, **MDTDb, **MDTSa, **MDTSb;

int *Sindexer, *Snums;
double *Sscales, *SJ, *Shers, *Shersds;

char filename[500], filename2[500], filename3[500], filename4[500], filename5[500], filename6[500];
FILE **Minputs, *output, *output2, *output3, *output4, *output5, *output6;


//set num_fixed, total and stotal
num_fixed=num_covars+num_envs+num_tops;
total=num_kins+num_regs;
stotal=(size_t)ns*(size_t)(ns-1)/2;

//allocate variables

indexer=malloc(sizeof(int)*ns);
indexer2=malloc(sizeof(int)*ns);

thetas=malloc(sizeof(double)*num_fixed);
thetasds=malloc(sizeof(double)*num_fixed);
thetapvas=malloc(sizeof(double)*num_fixed);
Yadj=malloc(sizeof(double)*ns);

hers=malloc(sizeof(double)*(total+1));
hersds=malloc(sizeof(double)*(total+1));

if(total>0)	//using kinships
{
shares=malloc(sizeof(double)*total);
sharesds=malloc(sizeof(double)*total);
cohers=malloc(sizeof(double)*total*total);
cohers2=malloc(sizeof(double)*total*total);

Dtemp=malloc(sizeof(double)*ns*total);
DTD=malloc(sizeof(double)*total*total);
DTS=malloc(sizeof(double)*total);
DTDs=malloc(sizeof(double)*total*total);
DTSs=malloc(sizeof(double)*total);
DTD2=malloc(sizeof(double)*total);
MDTD=malloc(sizeof(double*)*num_blocks);
MDTS=malloc(sizeof(double*)*num_blocks);
for(p=0;p<num_blocks;p++){MDTD[p]=malloc(sizeof(double)*total*total);MDTS[p]=malloc(sizeof(double)*total);}

DTDa=malloc(sizeof(double)*total*total);
DTDb=malloc(sizeof(double)*total*total);
DTSa=malloc(sizeof(double)*total);
DTSb=malloc(sizeof(double)*total);
MDTDa=malloc(sizeof(double*)*num_blocks);
MDTDb=malloc(sizeof(double*)*num_blocks);
MDTSa=malloc(sizeof(double*)*num_blocks);
MDTSb=malloc(sizeof(double*)*num_blocks);
for(p=0;p<num_blocks;p++)
{
MDTDa[p]=malloc(sizeof(double)*total*total);
MDTDb[p]=malloc(sizeof(double)*total*total);
MDTSa[p]=malloc(sizeof(double)*total);
MDTSb[p]=malloc(sizeof(double)*total);
}
}

subs=malloc(sizeof(int)*ns);
block=malloc(sizeof(int)*ns);

jacks=malloc(sizeof(double)*(total+1+total)*num_blocks);

////////

if(num_kins>0&&memsave==1)	//set indexers based on first id file (if nk>1, have checked all id files match) - also allocate datatemp
{
sprintf(filename, "%s.grm.id", kinstems[0]);
kcount=countrows(filename);
wantids=malloc(sizeof(char*)*kcount);
read_ids(filename, NULL, NULL, wantids, kcount, NULL, 0, 0);
count=find_strings(wantids, kcount, ids3, ns, indexer, indexer2, NULL, NULL, NULL, NULL, 3);
if(count!=ns){printf("Doug Error %d %d\n", count, ns);exit(1);}
for(i=0;i<kcount;i++){free(wantids[i]);}free(wantids);

datatemp=malloc(sizeof(float)*kcount);
}
else	//order of kinships matches that of phenotypes and regions
{
for(i=0;i<ns;i++){indexer[i]=i;indexer2[i]=i;}
}

//work out which subset each individual is in (each will be in exactly one)
for(s=0;s<num_subs;s++)
{
for(i=0;i<subindex[s][0];i++){subs[subindex[s][1+i]]=s;}
}

//divide individuals into blocks
for(p=0;p<num_blocks;p++)
{
start=(double)p/num_blocks*ns;
end=(double)(p+1)/num_blocks*ns;
for(i=start;i<end;i++){block[i]=p;}
}
permute_int(block,ns);

//solve covariates (get covher and topher) and fill Yadj
if(type==0)	//linear model - Yadj contains standardized residuals
{reg_covar_lin(Y, Z, ns, num_covars+num_envs, num_tops, thetas, thetasds, thetapvas, Yadj, 1, &covher, &topher);}
else	//logistic model - Yadj contains pcgc residuals
{reg_covar_log(Y, Z, ns, num_covars+num_envs, num_tops, NULL, thetas, thetasds, thetapvas, Yadj, 0, &covher, &topher, prev, 0.001, 100);}

if(num_covars>1){printf("Proportion of variance explained by the %d covariates: %.4f\n", num_covars, covher);}
if(num_tops==1){printf("Proportion of variance explained by the top predictor: %.4f\n", topher);}
if(num_tops>1){printf("Proportion of variance explained by the %d top predictor: %.4f\n", num_tops, topher);}
if(num_fixed>1){printf("\n");}

if(num_regs>0)	//regress covariates out of regions and update Xsums
{
reg_covar_matrix(X, Z, ns, Xtotal, num_fixed);

for(r=0;r<num_regs;r++)
{
sumsq=0;
for(j=Xstarts[r];j<Xends[r];j++)
{
for(i=0;i<ns;i++){sumsq+=pow(X[i+j*ns],2);}
}
Xsums[r]=sumsq/ns;
}
}

if(num_kins>0&&memsave==1)	//open and check kinships
{
Minputs=malloc(sizeof(FILE *)*num_kins);
for(k=0;k<num_kins;k++)
{
sprintf(filename, "%s.grm.bin", kinstems[k]);
if((Minputs[k]=fopen(filename,"rb"))==NULL)
{printf("Error opening %s\n\n", filename);exit(1);}
fseeko(Minputs[k], 0, SEEK_END);
if(ftello(Minputs[k])!=(off_t)sizeof(float)*kcount*(kcount+1)/2)
{printf("Error reading %s; should have size %jd not %jd\n\n", filename, (off_t)sizeof(float)*kcount*(kcount+1)/2, ftello(Minputs[k]));exit(1);}
}	//end of k loop
}

////////

//ready to start - S = products, D = [K1 K2 etc]

if(total==0)	//null model - just get likelihoods
{
hers[0]=0;hersds[0]=0;
Ysumsq=0;
for(i=1;i<ns;i++)
{
for(i2=0;i2<i;i2++){Ysumsq+=pow(Yadj[i]*Yadj[i2],2);}
}
likenull=-.5*stotal*(1+log(2*M_PI*Ysumsq/stotal));
like=likenull;lrtstat=0;lrtpva=1;
}
else
{
sprintf(filename,"%s.progress", outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fclose(output);

//set values to zero
Ysumsqa=0;
Ysumsqb=0;
for(k=0;k<total;k++)
{
for(k2=0;k2<total;k2++){DTDa[k+k2*total]=0;DTDb[k+k2*total]=0;}
DTSa[k]=0;DTSb[k]=0;
}

for(p=0;p<num_blocks;p++)
{
for(k=0;k<total;k++)
{
for(k2=0;k2<total;k2++){MDTDa[p][k+k2*total]=0;MDTDb[p][k+k2*total]=0;}
MDTSa[p][k]=0;MDTSb[p][k]=0;
}}

//loop through, filling a sample's worth of contributions to DTSa, DTSb, DTDa and DTDb at a time
scounta=0;scountb=0;
for(i=1;i<ns;i++)
{
//load up first i-1 rows of Dtemp with scaled kinships
for(k=0;k<num_kins;k++)
{
if(memsave==1)	//fill up datatemp - kinships for i start at i*(i+1)/2
{
smark=sizeof(float)*indexer[i]*(indexer[i]+1)/2;
fseeko(Minputs[k], smark, SEEK_SET);
fread(datatemp, sizeof(float), indexer[i], Minputs[k]);
for(i2=0;i2<i;i2++){Dtemp[i2+k*ns]=datatemp[indexer[i2]]/kintraces[k];}
}
else	//kinships already stored
{
for(i2=0;i2<i;i2++){Dtemp[i2+k*ns]=Mkins_single[k][(size_t)indexer2[i]*ns+indexer2[i2]]/kintraces[k];}
}
}

for(r=0;r<num_regs;r++)	//will compute kinships manually (probably faster to create a temporary matrix and vector, then use dgemv)
{
for(i2=0;i2<i;i2++)
{
value=0;for(j=Xstarts[r];j<Xends[r];j++){value+=X[indexer2[i]+j*ns]*X[indexer2[i2]+j*ns];}
Dtemp[i2+(num_kins+r)*ns]=value/Xsums[r];
}
}

//now loop through pairs
for(i2=0;i2<i;i2++)
{
if((scounta+scountb)%100000000==0)
{
printf("Processing Sample Pair %jd of %jd\n", (scounta+scountb)+1, stotal);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Considering Sample Pair %jd of %jd\n", (scounta+scountb)+1, stotal);
fclose(output);
}

//will add to blocks of both individuals (if the two blocks the same, add to just one)
p=block[indexer2[i]];
p2=block[indexer2[i2]];

if(trun==-9999)	//using all values
{
if(subs[indexer2[i]]==subs[indexer2[i2]])	//in same subset - update Ysumsqa, DTDa, DTSa, MDTDa, MDTSa
{
Ysumsqa+=pow(Yadj[indexer2[i]]*Yadj[indexer2[i2]],2);
for(k=0;k<total;k++)
{
for(k2=0;k2<total;k2++)
{
DTDa[k+k2*total]+=Dtemp[i2+k*ns]*Dtemp[i2+k2*ns];
MDTDa[p][k+k2*total]+=Dtemp[i2+k*ns]*Dtemp[i2+k2*ns];
if(p2!=p){MDTDa[p2][k+k2*total]+=Dtemp[i2+k*ns]*Dtemp[i2+k2*ns];}
}
DTSa[k]+=Dtemp[i2+k*ns]*Yadj[indexer2[i]]*Yadj[indexer2[i2]];
MDTSa[p][k]+=Dtemp[i2+k*ns]*Yadj[indexer2[i]]*Yadj[indexer2[i2]];
if(p2!=p){MDTSa[p2][k]+=Dtemp[i2+k*ns]*Yadj[indexer2[i]]*Yadj[indexer2[i2]];}
}
scounta++;
}
else	//in different subsets
{
Ysumsqb+=pow(Yadj[indexer2[i]]*Yadj[indexer2[i2]],2);
for(k=0;k<total;k++)
{
for(k2=0;k2<total;k2++)
{
DTDb[k+k2*total]+=Dtemp[i2+k*ns]*Dtemp[i2+k2*ns];
MDTDb[p][k+k2*total]+=Dtemp[i2+k*ns]*Dtemp[i2+k2*ns];
if(p2!=p){MDTDb[p2][k+k2*total]+=Dtemp[i2+k*ns]*Dtemp[i2+k2*ns];}
}
DTSb[k]+=Dtemp[i2+k*ns]*Yadj[indexer2[i]]*Yadj[indexer2[i2]];
MDTSb[p][k]+=Dtemp[i2+k*ns]*Yadj[indexer2[i]]*Yadj[indexer2[i2]];
if(p2!=p){MDTSb[p2][k]+=Dtemp[i2+k*ns]*Yadj[indexer2[i]]*Yadj[indexer2[i2]];}
}
scountb++;
}
}
else	//will be truncating (based on the first kinship)
{
if(Dtemp[i2]<=trun)
{
if(subs[indexer2[i]]==subs[indexer2[i2]])	//in same subset - update Ysumsqa, DTDa, DTSa, MDTDa, MDTSa
{
Ysumsqa+=pow(Yadj[indexer2[i]]*Yadj[indexer2[i2]],2);

for(k=0;k<total;k++)
{
for(k2=0;k2<total;k2++)
{
DTDa[k+k2*total]+=Dtemp[i2+k*ns]*Dtemp[i2+k2*ns];
MDTDa[p][k+k2*total]+=Dtemp[i2+k*ns]*Dtemp[i2+k2*ns];
if(p2!=p){MDTDa[p2][k+k2*total]+=Dtemp[i2+k*ns]*Dtemp[i2+k2*ns];}
}
DTSa[k]+=Dtemp[i2+k*ns]*Yadj[indexer2[i]]*Yadj[indexer2[i2]];
MDTSa[p][k]+=Dtemp[i2+k*ns]*Yadj[indexer2[i]]*Yadj[indexer2[i2]];
if(p2!=p){MDTSa[p2][k]+=Dtemp[i2+k*ns]*Yadj[indexer2[i]]*Yadj[indexer2[i2]];}
}
scounta++;
}
else	//in different subsets
{
Ysumsqb+=pow(Yadj[indexer2[i]]*Yadj[indexer2[i2]],2);
for(k=0;k<total;k++)
{
for(k2=0;k2<total;k2++)
{
DTDb[k+k2*total]+=Dtemp[i2+k*ns]*Dtemp[i2+k2*ns];
MDTDb[p][k+k2*total]+=Dtemp[i2+k*ns]*Dtemp[i2+k2*ns];
if(p2!=p){MDTDb[p2][k+k2*total]+=Dtemp[i2+k*ns]*Dtemp[i2+k2*ns];}
}
DTSb[k]+=Dtemp[i2+k*ns]*Yadj[indexer2[i]]*Yadj[indexer2[i2]];
MDTSb[p][k]+=Dtemp[i2+k*ns]*Yadj[indexer2[i]]*Yadj[indexer2[i2]];
if(p2!=p){MDTSb[p2][k]+=Dtemp[i2+k*ns]*Yadj[indexer2[i]]*Yadj[indexer2[i2]];}
}
scountb++;
}
}
}
}}	//end i2 and i loop
printf("\n");

//now add as and bs together
Ysumsq=Ysumsqa+Ysumsqb;

for(k=0;k<total;k++)
{
for(k2=0;k2<total;k2++)
{
DTD[k+k2*total]=DTDa[k+k2*total]+DTDb[k+k2*total];
for(p=0;p<num_blocks;p++){MDTD[p][k+k2*total]=MDTDa[p][k+k2*total]+MDTDb[p][k+k2*total];}
}
DTS[k]=DTSa[k]+DTSb[k];
for(p=0;p<num_blocks;p++){MDTS[p][k]=MDTSa[p][k]+MDTSb[p][k];}
}

if(num_kins>0&&memsave==1)	//can shut kinships
{
for(k=0;k<num_kins;k++){fclose(Minputs[k]);}free(Minputs);
}

////////

//save DTD and DTS
for(k=0;k<total;k++)
{
for(k2=0;k2<total;k2++){DTDs[k+k2*total]=DTD[k+k2*total];}
DTSs[k]=DTS[k];
}

//solve for all
printf("Analysing all %jd pairs of samples\n", scounta+scountb);
(void)eigen_invert(DTD, total, DTD2, 1, DTS, 1);
sum=0;for(k=0;k<total;k++){sum+=DTS[k];}
for(k=0;k<total;k++){hers[k]=DTS[k];}
hers[total]=sum;
for(k=0;k<total;k++){shares[k]=DTS[k]/sum;}

//get likelihoods and lrts
likenull=-.5*(scounta+scountb)*(1+log(2*M_PI*Ysumsq/(scounta+scountb)));
for(k=0;k<total;k++){Ysumsq-=DTSs[k]*DTS[k];}
like=-.5*(scounta+scountb)*(1+log(2*M_PI*Ysumsq/(scounta+scountb)));
lrtstat=2*(like-likenull);
lrtpva=erfc(pow(lrtstat,.5)*M_SQRT1_2);

for(p=0;p<num_blocks;p++)
{
if(p%5000==0){printf("Performing Jackknife %d out of %d\n", p+1, num_blocks);}
for(k=0;k<total;k++)
{
for(k2=0;k2<total;k2++){DTD[k+k2*total]=DTDs[k+k2*total]-MDTD[p][k+k2*total];}
DTS[k]=DTSs[k]-MDTS[p][k];
}

(void)eigen_invert(DTD, total, DTD2, 1, DTS, 1);
sum=0;for(k=0;k<total;k++){sum+=DTS[k];}
for(k=0;k<total;k++){jacks[k+p*(total+1+total)]=DTS[k];}
jacks[total+p*(total+1+total)]=sum;
for(k=0;k<total;k++){jacks[total+1+k+p*(total+1+total)]=DTS[k]/sum;}
}	//end of p loop

for(k=0;k<total+1+total;k++)	//get sds
{
sum=0;sumsq=0;
for(p=0;p<num_blocks;p++){sum+=jacks[k+p*(total+1+total)];sumsq+=pow(jacks[k+p*(total+1+total)],2);}
mean=sum/num_blocks;
var=(num_blocks-1)*(sumsq/num_blocks-pow(mean,2));
if(k<total+1){hersds[k]=pow(var,.5);}
else{sharesds[k-total-1]=pow(var,.5);}
}

//get covariances for hers
for(k=0;k<total;k++)
{
for(k2=0;k2<total;k2++)
{
sum=0;sum2=0;sumsq=0;
for(p=0;p<num_blocks;p++)
{
sum+=jacks[k+p*(total+1+total)];
sum2+=jacks[k2+p*(total+1+total)];
sumsq+=jacks[k+p*(total+1+total)]*jacks[k2+p*(total+1+total)];
}
mean=sum/num_blocks;mean2=sum2/num_blocks;
var=(num_blocks-1)*(sumsq/num_blocks-mean*mean2);
cohers[k+k2*total]=var;
}}

//and covariances for shares
for(k=0;k<total;k++)
{
for(k2=0;k2<total;k2++)
{
sum=0;sum2=0;sumsq=0;
for(p=0;p<num_blocks;p++)
{
sum+=jacks[total+1+k+p*(total+1+total)];
sum2+=jacks[total+1+k2+p*(total+1+total)];
sumsq+=jacks[total+1+k+p*(total+1+total)]*jacks[total+1+k2+p*(total+1+total)];
}
mean=sum/num_blocks;mean2=sum2/num_blocks;
var=(num_blocks-1)*(sumsq/num_blocks-mean*mean2);
cohers2[k+k2*total]=var;
}}

printf("\nTotal heritability %.4f (SE %.4f)\n\n", hers[total], hersds[total]);
}	//end of total>0

////////

//adjust for tops
for(k=0;k<total+1;k++){hers[k]*=(1-topher);hersds[k]*=(1-topher);}

//save stuff

flag=0;
for(k=0;k<num_kins;k++){flag+=(kinsums[k]==-9999);}
if(flag==0)
{
sum=0;
for(k=0;k<num_kins;k++){sum+=kinsums[k];}
for(r=0;r<num_regs;r++){sum+=Xsums[r];}
}

if(type==0){sprintf(filename2,"%s.he", outfile);}
else{sprintf(filename2,"%s.pcgc", outfile);}
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2, "Num_Kinships %d\nNum_Regions %d\nNum_Top_Predictors %d\nNum_Covariates %d\nNum_Environments %d\n", num_kins, num_regs, num_tops, num_covars, num_envs);
fprintf(output2, "Coeffsfile %s.coeff\nCovar_Heritability %.4f\n", outfile, covher);
fprintf(output2, "Total_Samples %d\nWith_Phenotypes %d\n", ns, np);
fprintf(output2, "Null_Likelihood %.6f\nAlt_Likelihood %.6f\n", likenull, like);
if(total==1){fprintf(output2, "LRT_Stat %.4f\nLRT_P %.4e\n", lrtstat, lrtpva);}
else{fprintf(output2, "LRT_Stat %.4f\nLRT_P NA\n", lrtstat);}

fprintf(output2, "Component Heritability SE Size Mega_Intensity SE\n");
if(flag==0)	//might have null model
{
for(k=0;k<num_kins;k++){fprintf(output2, "Her_K%d %.6f %.6f %.2f %.6f %.6f\n", k+1, hers[k], hersds[k], kinsums[k], hers[k]/kinsums[k]*1000000, hersds[k]/kinsums[k]*1000000);}
for(r=0;r<num_regs;r++){fprintf(output2, "Her_R%d %.6f %.6f %.2f %.6f %.6f\n", r+1, hers[num_kins+r], hersds[num_kins+r], Xsums[r], hers[num_kins+r]/Xsums[r]*1000000, hersds[num_kins+r]/Xsums[r]*1000000);}
fprintf(output2, "Her_Top %.6f NA NA NA NA\n", topher);
if(total==0){fprintf(output2, "Her_All %.6f NA NA NA NA\n", topher);}
else{fprintf(output2, "Her_All %.6f %.6f %.2f %.6f %.6f\n", hers[total]+topher, hersds[total], sum, hers[total]/sum*1000000, hersds[total]/sum*1000000);}
}
else	//must be non-null
{
for(k=0;k<num_kins;k++){fprintf(output2, "Her_K%d %.6f %.6f NA NA NA\n", k+1, hers[k], hersds[k]);}
for(r=0;r<num_regs;r++){fprintf(output2, "Her_R%d %.6f %.6f NA NA NA\n", r+1, hers[num_kins+r], hersds[num_kins+r]);}
fprintf(output2, "Her_Top %.6f NA NA NA NA\n", topher);
fprintf(output2, "Her_All %.6f %.6f NA NA NA\n", hers[total]+topher, hersds[total]);
}
fclose(output2);

sprintf(filename3,"%s.coeff", outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
if(type==0){fprintf(output3, "Component Effect SE P\n");}
else{fprintf(output3, "Component Log_Odds SE P\n");}
fprintf(output3, "Intercept %.4e %.4e %.4e\n", thetas[0], thetasds[0], thetapvas[0]);
for(j=1;j<num_covars;j++){fprintf(output3, "Covariate_%d %.4e %.4e %.4e\n",j, thetas[j], thetasds[j], thetapvas[j]);}
for(j=0;j<num_envs;j++){fprintf(output3, "Enviromental_%d %.4e %.4e %.4e\n",j, thetas[num_covars+j], thetasds[num_covars+j], thetapvas[num_covars+j]);}
fclose(output3);

sprintf(filename4,"%s.share", outfile);
if((output4=fopen(filename4,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename4);exit(1);}
fprintf(output4, "Component Share SE Expected Enrichment SE\n");
if(flag==0)
{
for(k=0;k<num_kins;k++){fprintf(output4, "Share_K%d %.6f %.6f %.6f %.6f %.6f\n", k+1, shares[k], sharesds[k], kinsums[k]/sum, shares[k]/kinsums[k]*sum, sharesds[k]/kinsums[k]*sum);}
for(r=0;r<num_regs;r++){fprintf(output4, "Share_R%d %.6f %.6f %.6f %.6f %.6f\n", r+1, shares[num_kins+r], sharesds[num_kins+r], Xsums[r]/sum, shares[num_kins+r]/Xsums[r]*sum, sharesds[num_kins+r]/Xsums[r]*sum);}
}
else
{
for(k=0;k<num_kins;k++){fprintf(output4, "Share_K%d %.6f %.6f NA NA NA\n", k+1, shares[k], sharesds[k]);}
for(r=0;r<num_regs;r++){fprintf(output4, "Share_R%d %.6f %.6f NA NA NA\n", r+1, shares[num_kins+r], sharesds[num_kins+r]);}
}
fclose(output4);

if(type==0&&prev!=-9999)	//binary he
{
factor=get_factor(Y, ns, prev, -9999, outfile);

sprintf(filename5,"%s.he.liab", outfile);
if((output5=fopen(filename5,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename5);exit(1);}
fprintf(output5, "Num_Kinships %d\nNum_Regions %d\nNum_Top_Predictors %d\nNum_Covariates %d\nNum_Environments %d\n", num_kins, num_regs, num_tops, num_covars, num_envs);
fprintf(output5, "Coeffsfile %s.coeff.liab\nCovar_Heritability %.4f\n", outfile, covher*factor);
fprintf(output5, "Total_Samples %d\nWith_Phenotypes %d\n", ns, np);
fprintf(output5, "Null_Likelihood %.6f\nAlt_Likelihood %.6f\n", likenull, like);
if(total==1){fprintf(output5, "LRT_Stat %.4f\nLRT_P %.4e\n", lrtstat, lrtpva);}
else{fprintf(output5, "LRT_Stat %.4f\nLRT_P NA\n", lrtstat);}

fprintf(output5, "Component Heritability SE Size Mega_Intensity SE\n");
if(flag==0)	//might have null model
{
for(k=0;k<num_kins;k++){fprintf(output5, "Her_K%d %.6f %.6f %.2f %.6f %.6f\n", k+1, hers[k]*factor, hersds[k]*factor, kinsums[k], hers[k]/kinsums[k]*1000000*factor, hersds[k]/kinsums[k]*1000000*factor);}
for(r=0;r<num_regs;r++){fprintf(output5, "Her_R%d %.6f %.6f %.2f %.6f %.6f\n", r+1, hers[num_kins+r]*factor, hersds[num_kins+r]*factor, Xsums[r], hers[num_kins+r]/Xsums[r]*1000000*factor, hersds[num_kins+r]/Xsums[r]*1000000*factor);}
fprintf(output5, "Her_Top %.6f NA NA NA NA\n", topher*factor);
if(total==0){fprintf(output5, "Her_All %.6f NA NA NA NA\n", topher);}
else{fprintf(output5, "Her_All %.6f %.6f %.2f %.6f %.6f\n", hers[total]*factor+topher*factor, hersds[total]*factor, sum, hers[total]/sum*1000000*factor, hersds[total]/sum*1000000*factor);}
}
else	//must be non-null
{
for(k=0;k<num_kins;k++){fprintf(output5, "Her_K%d %.6f %.6f NA NA NA\n", k+1, hers[k]*factor, hersds[k]*factor);}
for(r=0;r<num_regs;r++){fprintf(output5, "Her_R%d %.6f %.6f NA NA NA\n", r+1, hers[num_kins+r]*factor, hersds[num_kins+r]*factor);}
fprintf(output5, "Her_Top %.6f NA NA NA NA\n", topher*factor);
fprintf(output5, "Her_All %.6f %.6f NA NA NA\n", hers[total]*factor+topher*factor, hersds[total]*factor);
}
fclose(output5);
}

if(type==1)	//save marginals for pcgc
{
sprintf(filename5,"%s.pcgc.marginal", outfile);
if((output5=fopen(filename5,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename5);exit(1);}
fprintf(output5, "Num_Kinships %d\nNum_Regions %d\nNum_Top_Predictors %d\nNum_Covariates %d\nNum_Environments %d\n", num_kins, num_regs, num_tops, num_covars, num_envs);
fprintf(output5, "Coeffsfile %s.coeff\nCovar_Heritability %.4f\n", outfile, covher);
fprintf(output5, "Total_Samples %d\nWith_Phenotypes %d\n", ns, np);
fprintf(output5, "Null_Likelihood %.6f\nAlt_Likelihood %.6f\n", likenull, like);
if(total==1){fprintf(output5, "LRT_Stat %.4f\nLRT_P %.4e\n", lrtstat, lrtpva);}
else{fprintf(output5, "LRT_Stat %.4f\nLRT_P NA\n", lrtstat);}

fprintf(output5, "Component Heritability SE Size Mega_Intensity SE\n");
if(flag==0)	//might have null model
{
for(k=0;k<num_kins;k++){fprintf(output5, "Her_K%d %.6f %.6f %.2f %.6f %.6f\n", k+1, hers[k]*(1-covher), hersds[k]*(1-covher), kinsums[k], hers[k]/kinsums[k]*1000000*(1-covher), hersds[k]/kinsums[k]*1000000*(1-covher));}
for(r=0;r<num_regs;r++){fprintf(output5, "Her_R%d %.6f %.6f %.2f %.6f %.6f\n", r+1, hers[num_kins+r]*(1-covher), hersds[num_kins+r]*(1-covher), Xsums[r], hers[num_kins+r]/Xsums[r]*1000000*(1-covher), hersds[num_kins+r]/Xsums[r]*1000000*(1-covher));}
fprintf(output5, "Her_Top %.6f NA NA NA NA\n", topher*(1-covher));
if(total==0){fprintf(output5, "Her_All %.6f NA NA NA NA\n", topher*(1-covher));}
else{fprintf(output5, "Her_All %.6f %.6f %.2f %.6f %.6f\n", hers[total]*(1-covher)+topher*(1-covher), hersds[total]*(1-covher), sum, hers[total]/sum*1000000*(1-covher), hersds[total]/sum*1000000*(1-covher));}
}
else	//must be non-null
{
for(k=0;k<num_kins;k++){fprintf(output5, "Her_K%d %.6f %.6f NA NA NA\n", k+1, hers[k]*(1-covher), hersds[k]*(1-covher));}
for(r=0;r<num_regs;r++){fprintf(output5, "Her_R%d %.6f %.6f NA NA NA\n", r+1, hers[num_kins+r]*(1-covher), hersds[num_kins+r]*(1-covher));}
fprintf(output5, "Her_Top %.6f NA NA NA NA\n", topher*(1-covher));
fprintf(output5, "Her_All %.6f %.6f NA NA NA\n", hers[total]*(1-covher)+topher*(1-covher), hersds[total]*(1-covher));
}
fclose(output5);
}

sprintf(filename6,"%s.cross", outfile);
if((output6=fopen(filename6,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename6);exit(1);}
for(k=0;k<num_kins;k++){fprintf(output6, "Her_K%d\t", k+1);}
for(r=0;r<num_regs;r++){fprintf(output6, "Her_R%d\t", r+1);}
fprintf(output6, "\n");
for(k=0;k<total;k++)
{
for(k2=0;k2<total;k2++)
{fprintf(output6, "%.6f\t", cohers[k+k2*total]);}
fprintf(output6, "\n");
}
fclose(output6);

////////

if(discenv==1)	//get shares for 1+num_envs+num_regs genetic kinships
{
printf("Computing heritabilities for %d subgroups\n", num_envs);

Sindexer=malloc(sizeof(int)*ns);
Snums=malloc(sizeof(int)*(1+num_envs));
Sscales=malloc(sizeof(double)*total);
SJ=malloc(sizeof(double)*total);
Shers=malloc(sizeof(double)*(1+num_envs));
Shersds=malloc(sizeof(double)*(1+num_envs));
if(memsave==1){wantids=malloc(sizeof(char*)*ns);}

//first do for all individuals - just need to add up genetic heritabilities
for(k=0;k<1+num_envs;k++){SJ[k]=1;}
for(k=1+num_envs;k<num_kins;k++){SJ[k]=0;}
for(r=0;r<num_regs;r++){SJ[num_kins+r]=1;}

sum=0;for(k=0;k<total;k++){sum+=hers[k]*SJ[k];}
value=0;for(k=0;k<total;k++){for(k2=0;k2<total;k2++){value+=cohers[k+k2*total]*SJ[k]*SJ[k2];}}

Snums[num_envs]=ns;
Shers[num_envs]=sum;
if(value>=0){Shersds[num_envs]=pow(value,.5);}
else{Shersds[num_envs]=-9999;}

//now for each subgroup
for(j=0;j<num_envs;j++)
{
//work out individuals in subgroup
count=0;
for(i=0;i<ns;i++)
{
if(Z[i+(num_covars+num_tops+j)*ns]==1){Sindexer[count]=i;count++;}
}

if(memsave==1)	//load these into wantids
{
for(i=0;i<count;i++){wantids[i]=ids3[Sindexer[i]];}
}

//get traces across these individuals
for(k=0;k<num_kins;k++)
{
if(memsave==0)
{
sum=0;for(i=0;i<count;i++){sum+=Mkins_single[k][(size_t)Sindexer[i]*ns+Sindexer[i]];}
Sscales[k]=sum/count;
}
else{Sscales[k]=read_kin_trace(kinstems[k], count, wantids, NULL, 1, 0, 1);}
}
for(r=0;r<num_regs;r++)
{
sum=0;
for(i=0;i<count;i++)
{
for(k=Xstarts[r];k<Xends[r];k++){sum+=pow(X[Sindexer[i]+k*ns],2);}
}
Sscales[num_kins+r]=sum/count;
}

//need to add up genetic heritabilities, scaled by subset trace over full trace
for(k=0;k<1+num_envs;k++){SJ[k]=Sscales[k]/kintraces[k];}
for(k=1+num_envs;k<num_kins;k++){SJ[k]=0;}
for(r=0;r<num_regs;r++){SJ[num_kins+r]=Sscales[num_kins+r]/Xsums[r];}

sum=0;for(k=0;k<total;k++){sum+=hers[k]*SJ[k];}
value=0;for(k=0;k<total;k++){for(k2=0;k2<total;k2++){value+=cohers[k+k2*total]*SJ[k]*SJ[k2];}}

Snums[j]=count;
Shers[j]=sum;
if(value>=0){Shersds[j]=pow(value,.5);}
else{Shersds[j]=-9999;}
}	//end of j loop

sprintf(filename,"%s.subgroups", outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fprintf(output, "Component Num_Samples Heritability SE\n");
for(j=0;j<num_envs;j++){fprintf(output, "Her_Sub%d %d %.6f %.6f\n", j+1, Snums[j], Shers[j], Shersds[j]);}
fprintf(output, "Her_All %d %.6f %.6f\n", Snums[num_envs], Shers[num_envs], Shersds[num_envs]);
fclose(output);

printf("Estimates saved in %s\n\n", filename);
free(Sindexer);free(Snums);free(Shers);free(Shersds);free(SJ);
if(memsave==1){free(wantids);}
}	//end of discenv=1

if(strcmp(oversfile,"blank")!=0)	//compute and save cats and enrichments - cant have regions
{
sprintf(filename,"%s.cats", outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}

fprintf(output, "Component Heritability SE\n");
for(k=0;k<num_kins;k++)
{
sum=0;sum2=0;
for(k2=0;k2<num_kins;k2++)
{
sum+=hers[k2]*ssums[k][k2];
for(k3=0;k3<num_kins;k3++){sum2+=cohers[k2+k3*total]*ssums[k][k2]*ssums[k][k3];}
}
if(sum2>0){value=pow(sum2,.5);}
else{value=-9999;}
fprintf(output, "Cat_K%d %.6f %.6f\n", k+1, sum, value);
}
fclose(output);

sprintf(filename,"%s.enrich", outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}

fprintf(output, "Component Share SE Expected Enrichment SE\n");
for(k=0;k<num_kins;k++)
{
sum=0;sum2=0;
for(k2=0;k2<num_kins;k2++)
{
sum+=shares[k2]*ssums[k][k2];
for(k3=0;k3<num_kins;k3++){sum2+=cohers2[k2+k3*total]*ssums[k][k2]*ssums[k][k3];}
}
if(sum2>0){value=pow(sum2,.5);}
else{value=-9999;}
fprintf(output, "Enrich_K%d %.6f %.6f %.6f %.6f %.6f\n", k+1, sum, value, ssums[k][num_kins+1], sum/ssums[k][num_kins+1], value/ssums[k][num_kins+1]);
}
fclose(output);
}

////////

if(num_subs>1)	//now solve for within and across
{
//first within subsets - save DTDa and DTSa, solve and get likelihood and LRTs
for(k=0;k<total;k++)
{
for(k2=0;k2<total;k2++){DTDs[k+k2*total]=DTDa[k+k2*total];}
DTSs[k]=DTSa[k];
}

printf("Analysing %jd pairs of samples in the same subsets\n", scounta);
(void)eigen_invert(DTDa, total, DTD2, 1, DTSa, 1);
sum=0;for(k=0;k<total;k++){sum+=DTSa[k];}
for(k=0;k<total;k++){hers[k]=DTSa[k];}
hers[total]=sum;
for(k=0;k<total;k++){shares[k]=DTSa[k]/sum;}

likenull=-.5*scounta*(1+log(2*M_PI*Ysumsqa/scounta));
for(k=0;k<total;k++){Ysumsq-=DTSs[k]*DTSa[k];}
like=-.5*scounta*(1+log(2*M_PI*Ysumsqa/scounta));
lrtstat=2*(like-likenull);
lrtpva=erfc(pow(lrtstat,.5)*M_SQRT1_2);

for(p=0;p<num_blocks;p++)
{
if(p%5000==0){printf("Performing within-subset Jackknife %d out of %d\n", p+1, num_blocks);}
for(k=0;k<total;k++)
{
for(k2=0;k2<total;k2++){DTDa[k+k2*total]=DTDs[k+k2*total]-MDTDa[p][k+k2*total];}
DTSa[k]=DTSs[k]-MDTSa[p][k];
}

(void)eigen_invert(DTDa, total, DTD2, 1, DTSa, 1);
sum=0;for(k=0;k<total;k++){sum+=DTSa[k];}
for(k=0;k<total;k++){jacks[k+p*(total+1+total)]=DTSa[k];}
jacks[total+p*(total+1+total)]=sum;
for(k=0;k<total;k++){jacks[total+1+k+p*(total+1+total)]=DTSa[k]/sum;}
}	//end of p loop

for(k=0;k<total+1+total;k++)	//get sds
{
sum=0;sumsq=0;
for(p=0;p<num_blocks;p++){sum+=jacks[k+p*(total+1+total)];sumsq+=pow(jacks[k+p*(total+1+total)],2);}
mean=sum/num_blocks;
var=(num_blocks-1)*(sumsq/num_blocks-pow(mean,2));
if(k<total+1){hersds[k]=pow(var,.5);}
else{sharesds[k-total-1]=pow(var,.5);}
}

//get covariances
for(k=0;k<total;k++)
{
for(k2=0;k2<total;k2++)
{
sum=0;sum2=0;sumsq=0;
for(p=0;p<num_blocks;p++)
{
sum+=jacks[k+p*(total+1+total)];
sum2+=jacks[k2+p*(total+1+total)];
sumsq+=jacks[k+p*(total+1+total)]*jacks[k2+p*(total+1+total)];
}
mean=sum/num_blocks;mean2=sum2/num_blocks;
var=(num_blocks-1)*(sumsq/num_blocks-mean*mean2);
cohers[k+k2*total]=var;
}}

printf("\nTotal heritability %.4f (SE %.4f)\n\n", hers[total], hersds[total]);
stats[0]=hers[total];stats[1]=hersds[total];

//adjust for tops then save (already have flag, but need to recompute sum)

for(k=0;k<total+1;k++){hers[k]*=(1-topher);hersds[k]*=(1-topher);}
if(flag==0)
{
sum=0;for(k=0;k<num_kins;k++){sum+=kinsums[k];}for(r=0;r<num_regs;r++){sum+=Xsums[r];}
}

if(type==0){sprintf(filename2,"%s.he.within", outfile);}
else{sprintf(filename2,"%s.pcgc.within", outfile);}
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2, "Num_Kinships %d\nNum_Regions %d\nNum_Top_Predictors %d\nNum_Covariates %d\nNum_Environments %d\n", num_kins, num_regs, num_tops, num_covars, num_envs);
fprintf(output2, "Coeffsfile %s.coeff\nCovar_Heritability %.4f\n", outfile, covher);
fprintf(output2, "Total_Samples %d\nWith_Phenotypes %d\n", ns, np);
fprintf(output2, "Null_Likelihood %.6f\nAlt_Likelihood %.6f\n", likenull, like);
if(total==1){fprintf(output2, "LRT_Stat %.4f\nLRT_P %.4e\n", lrtstat, lrtpva);}
else{fprintf(output2, "LRT_Stat %.4f\nLRT_P NA\n", lrtstat);}

fprintf(output2, "Component Heritability SE Size Mega_Intensity SE\n");
if(flag==0)	//might have null model (although would be a very weird analysis)
{
for(k=0;k<num_kins;k++){fprintf(output2, "Her_K%d %.6f %.6f %.2f %.6f %.6f\n", k+1, hers[k], hersds[k], kinsums[k], hers[k]/kinsums[k]*1000000, hersds[k]/kinsums[k]*1000000);}
for(r=0;r<num_regs;r++){fprintf(output2, "Her_R%d %.6f %.6f %.2f %.6f %.6f\n", r+1, hers[num_kins+r], hersds[num_kins+r], Xsums[r], hers[num_kins+r]/Xsums[r]*1000000, hersds[num_kins+r]/Xsums[r]*1000000);}
fprintf(output2, "Her_Top %.6f NA NA NA NA\n", topher);
if(total==0){fprintf(output2, "Her_All %.6f NA NA NA NA\n", topher);}
else{fprintf(output2, "Her_All %.6f %.6f %.2f %.6f %.6f\n", hers[total]+topher, hersds[total], sum, hers[total]/sum*1000000, hersds[total]/sum*1000000);}
}
else	//must be non-null
{
for(k=0;k<num_kins;k++){fprintf(output2, "Her_K%d %.6f %.6f NA NA NA\n", k+1, hers[k], hersds[k]);}
for(r=0;r<num_regs;r++){fprintf(output2, "Her_R%d %.6f %.6f NA NA NA\n", r+1, hers[num_kins+r], hersds[num_kins+r]);}
fprintf(output2, "Her_Top %.6f NA NA NA NA\n", topher);
fprintf(output2, "Her_All %.6f %.6f NA NA NA\n", hers[total]+topher, hersds[total]);
}
fclose(output2);

sprintf(filename4,"%s.share.within", outfile);
if((output4=fopen(filename4,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename4);exit(1);}
fprintf(output4, "Component Share SE Expected Enrichment SE\n");
if(flag==0)
{
for(k=0;k<num_kins;k++){fprintf(output4, "Share_K%d %.6f %.6f %.6f %.6f %.6f\n", k+1, shares[k], sharesds[k], kinsums[k]/sum, shares[k]/kinsums[k]*sum, sharesds[k]/kinsums[k]*sum);}
for(r=0;r<num_regs;r++){fprintf(output4, "Share_R%d %.6f %.6f %.6f %.6f %.6f\n", r+1, shares[num_kins+r], sharesds[num_kins+r], Xsums[r]/sum, shares[num_kins+r]/Xsums[r]*sum, sharesds[num_kins+r]/Xsums[r]*sum);}
}
else
{
for(k=0;k<num_kins;k++){fprintf(output4, "Share_K%d %.6f %.6f NA NA NA\n", k+1, shares[k], sharesds[k]);}
for(r=0;r<num_regs;r++){fprintf(output4, "Share_R%d %.6f %.6f NA NA NA\n", r+1, shares[num_kins+r], sharesds[num_kins+r]);}
}
fclose(output4);

if(type==0&&prev!=-9999)	//binary he - already have factor
{
sprintf(filename5,"%s.he.within.liab", outfile);
if((output5=fopen(filename5,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename5);exit(1);}
fprintf(output5, "Num_Kinships %d\nNum_Regions %d\nNum_Top_Predictors %d\nNum_Covariates %d\nNum_Environments %d\n", num_kins, num_regs, num_tops, num_covars, num_envs);
fprintf(output5, "Coeffsfile %s.coeff.liab\nCovar_Heritability %.4f\n", outfile, covher*factor);
fprintf(output5, "Total_Samples %d\nWith_Phenotypes %d\n", ns, np);
fprintf(output5, "Null_Likelihood %.6f\nAlt_Likelihood %.6f\n", likenull, like);
if(total==1){fprintf(output5, "LRT_Stat %.4f\nLRT_P %.4e\n", lrtstat, lrtpva);}
else{fprintf(output5, "LRT_Stat %.4f\nLRT_P NA\n", lrtstat);}

fprintf(output5, "Component Heritability SE Size Mega_Intensity SE\n");
if(flag==0)
{
for(k=0;k<num_kins;k++){fprintf(output5, "Her_K%d %.6f %.6f %.2f %.6f %.6f\n", k+1, hers[k]*factor, hersds[k]*factor, kinsums[k], hers[k]/kinsums[k]*1000000*factor, hersds[k]/kinsums[k]*1000000*factor);}
for(r=0;r<num_regs;r++){fprintf(output5, "Her_R%d %.6f %.6f %.2f %.6f %.6f\n", r+1, hers[num_kins+r]*factor, hersds[num_kins+r]*factor, Xsums[r], hers[num_kins+r]/Xsums[r]*1000000*factor, hersds[num_kins+r]/Xsums[r]*1000000*factor);}
fprintf(output5, "Her_Top %.6f NA NA NA NA\n", topher*factor);
if(total==0){fprintf(output5, "Her_All %.6f NA NA NA NA\n", topher);}
else{fprintf(output5, "Her_All %.6f %.6f %.2f %.6f %.6f\n", hers[total]*factor+topher*factor, hersds[total]*factor, sum, hers[total]/sum*1000000*factor, hersds[total]/sum*1000000*factor);}
}
else
{
for(k=0;k<num_kins;k++){fprintf(output5, "Her_K%d %.6f %.6f NA NA NA\n", k+1, hers[k]*factor, hersds[k]*factor);}
for(r=0;r<num_regs;r++){fprintf(output5, "Her_R%d %.6f %.6f NA NA NA\n", r+1, hers[num_kins+r]*factor, hersds[num_kins+r]*factor);}
fprintf(output5, "Her_Top %.6f NA NA NA NA\n", topher*factor);
fprintf(output5, "Her_All %.6f %.6f NA NA NA\n", hers[total]*factor+topher*factor, hersds[total]*factor);
}
fclose(output5);
}

if(type==1)
{
sprintf(filename5,"%s.pcgc.within.marginal", outfile);
if((output5=fopen(filename5,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename5);exit(1);}
fprintf(output5, "Num_Kinships %d\nNum_Regions %d\nNum_Top_Predictors %d\nNum_Covariates %d\nNum_Environments %d\n", num_kins, num_regs, num_tops, num_covars, num_envs);
fprintf(output5, "Coeffsfile %s.coeff\nCovar_Heritability %.4f\n", outfile, covher);
fprintf(output5, "Total_Samples %d\nWith_Phenotypes %d\n", ns, np);
fprintf(output5, "Null_Likelihood %.6f\nAlt_Likelihood %.6f\n", likenull, like);
if(total==1){fprintf(output5, "LRT_Stat %.4f\nLRT_P %.4e\n", lrtstat, lrtpva);}
else{fprintf(output5, "LRT_Stat %.4f\nLRT_P NA\n", lrtstat);}

fprintf(output5, "Component Heritability SE Size Mega_Intensity SE\n");
if(flag==0)	//might have null model
{
for(k=0;k<num_kins;k++){fprintf(output5, "Her_K%d %.6f %.6f %.2f %.6f %.6f\n", k+1, hers[k]*(1-covher), hersds[k]*(1-covher), kinsums[k], hers[k]*(1-covher)/kinsums[k]*1000000, hersds[k]*(1-covher)/kinsums[k]*1000000);}
for(r=0;r<num_regs;r++){fprintf(output5, "Her_R%d %.6f %.6f %.2f %.6f %.6f\n", r+1, hers[num_kins+r]*(1-covher), hersds[num_kins+r]*(1-covher), Xsums[r], hers[num_kins+r]*(1-covher)/Xsums[r]*1000000, hersds[num_kins+r]*(1-covher)/Xsums[r]*1000000);}
fprintf(output5, "Her_Top %.6f NA NA NA NA\n", topher*(1-covher));
if(total==0){fprintf(output5, "Her_All %.6f NA NA NA NA\n", topher*(1-covher));}
else{fprintf(output5, "Her_All %.6f %.6f %.2f %.6f %.6f\n", hers[total]*(1-covher)+topher*(1-covher), hersds[total]*(1-covher), sum, hers[total]*(1-covher)/sum*1000000, hersds[total]*(1-covher)/sum*1000000);}
}
else	//must be non-null
{
for(k=0;k<num_kins;k++){fprintf(output5, "Her_K%d %.6f %.6f NA NA NA\n", k+1, hers[k]*(1-covher), hersds[k]*(1-covher));}
for(r=0;r<num_regs;r++){fprintf(output5, "Her_R%d %.6f %.6f NA NA NA\n", r+1, hers[num_kins+r]*(1-covher), hersds[num_kins+r]*(1-covher));}
fprintf(output5, "Her_Top %.6f NA NA NA NA\n", topher*(1-covher));
fprintf(output5, "Her_All %.6f %.6f NA NA NA\n", hers[total]*(1-covher)+topher*(1-covher), hersds[total]*(1-covher));
}
fclose(output5);
}

sprintf(filename6,"%s.cross.within", outfile);
if((output6=fopen(filename6,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename6);exit(1);}
for(k=0;k<num_kins;k++){fprintf(output6, "Her_K%d\t", k+1);}
for(r=0;r<num_regs;r++){fprintf(output6, "Her_R%d\t", r+1);}
fprintf(output6, "\n");
for(k=0;k<total;k++)
{
for(k2=0;k2<total;k2++)
{fprintf(output6, "%.6f\t", cohers[k+k2*total]);}
fprintf(output6, "\n");
}
fclose(output6);

////////

//now across subsets - save DTDb and DTSb
for(k=0;k<total;k++)
{
for(k2=0;k2<total;k2++){DTDs[k+k2*total]=DTDb[k+k2*total];}
DTSs[k]=DTSb[k];
}

printf("Analysing %jd pairs of samples in different subsets\n", scountb);
(void)eigen_invert(DTDb, total, DTD2, 1, DTSb, 1);
sum=0;for(k=0;k<total;k++){sum+=DTSb[k];}
for(k=0;k<total;k++){hers[k]=DTSb[k];}
hers[total]=sum;
for(k=0;k<total;k++){shares[k]=DTSb[k]/sum;}

likenull=-.5*scountb*(1+log(2*M_PI*Ysumsqb/scountb));
for(k=0;k<total;k++){Ysumsq-=DTSs[k]*DTSb[k];}
like=-.5*scountb*(1+log(2*M_PI*Ysumsqb/scountb));
lrtstat=2*(like-likenull);
lrtpva=erfc(pow(lrtstat,.5)*M_SQRT1_2);

for(p=0;p<num_blocks;p++)
{
if(p%5000==0){printf("Performing across-subset Jackknife %d out of %d\n", p+1, num_blocks);}
for(k=0;k<total;k++)
{
for(k2=0;k2<total;k2++){DTDb[k+k2*total]=DTDs[k+k2*total]-MDTDb[p][k+k2*total];}
DTSb[k]=DTSs[k]-MDTSb[p][k];
}

(void)eigen_invert(DTDb, total, DTD2, 1, DTSb, 1);
sum=0;for(k=0;k<total;k++){sum+=DTSb[k];}
for(k=0;k<total;k++){jacks[k+p*(total+1+total)]=DTSb[k];}
jacks[total+p*(total+1+total)]=sum;
for(k=0;k<total;k++){jacks[total+1+k+p*(total+1+total)]=DTSb[k]/sum;}
}	//end of p loop

for(k=0;k<total+1+total;k++)	//get sds
{
sum=0;sumsq=0;
for(p=0;p<num_blocks;p++){sum+=jacks[k+p*(total+1+total)];sumsq+=pow(jacks[k+p*(total+1+total)],2);}
mean=sum/num_blocks;
var=(num_blocks-1)*(sumsq/num_blocks-pow(mean,2));
if(k<total+1){hersds[k]=pow(var,.5);}
else{sharesds[k-total-1]=pow(var,.5);}
}

//get covariances
for(k=0;k<total;k++)
{
for(k2=0;k2<total;k2++)
{
sum=0;sum2=0;sumsq=0;
for(p=0;p<num_blocks;p++)
{
sum+=jacks[k+p*(total+1+total)];
sum2+=jacks[k2+p*(total+1+total)];
sumsq+=jacks[k+p*(total+1+total)]*jacks[k2+p*(total+1+total)];
}
mean=sum/num_blocks;mean2=sum2/num_blocks;
var=(num_blocks-1)*(sumsq/num_blocks-mean*mean2);
cohers[k+k2*total]=var;
}}

printf("\nTotal heritability %.4f (SE %.4f)\n\n", hers[total], hersds[total]);
stats[2]=hers[total];stats[3]=hersds[total];

//adjust for tops then save (already have flag, but need to recompute sum)

for(k=0;k<total+1;k++){hers[k]*=(1-topher);hersds[k]*=(1-topher);}
if(flag==0)
{
sum=0;for(k=0;k<num_kins;k++){sum+=kinsums[k];}for(r=0;r<num_regs;r++){sum+=Xsums[r];}
}

if(type==0){sprintf(filename2,"%s.he.across", outfile);}
else{sprintf(filename2,"%s.pcgc.across", outfile);}
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2, "Num_Kinships %d\nNum_Regions %d\nNum_Top_Predictors %d\nNum_Covariates %d\nNum_Environments %d\n", num_kins, num_regs, num_tops, num_covars, num_envs);
fprintf(output2, "Coeffsfile %s.coeff\nCovar_Heritability %.4f\n", outfile, covher);
fprintf(output2, "Total_Samples %d\nWith_Phenotypes %d\n", ns, np);
fprintf(output2, "Null_Likelihood %.6f\nAlt_Likelihood %.6f\n", likenull, like);
if(total==1){fprintf(output2, "LRT_Stat %.4f\nLRT_P %.4e\n", lrtstat, lrtpva);}
else{fprintf(output2, "LRT_Stat %.4f\nLRT_P NA\n", lrtstat);}

fprintf(output2, "Component Heritability SE Size Mega_Intensity SE\n");
if(flag==0)	//might have null model (although would be a very weird analysis)
{
for(k=0;k<num_kins;k++){fprintf(output2, "Her_K%d %.6f %.6f %.2f %.6f %.6f\n", k+1, hers[k], hersds[k], kinsums[k], hers[k]/kinsums[k]*1000000, hersds[k]/kinsums[k]*1000000);}
for(r=0;r<num_regs;r++){fprintf(output2, "Her_R%d %.6f %.6f %.2f %.6f %.6f\n", r+1, hers[num_kins+r], hersds[num_kins+r], Xsums[r], hers[num_kins+r]/Xsums[r]*1000000, hersds[num_kins+r]/Xsums[r]*1000000);}
fprintf(output2, "Her_Top %.6f NA NA NA NA\n", topher);
if(total==0){fprintf(output2, "Her_All %.6f NA NA NA NA\n", topher);}
else{fprintf(output2, "Her_All %.6f %.6f %.2f %.6f %.6f\n", hers[total]+topher, hersds[total], sum, hers[total]/sum*1000000, hersds[total]/sum*1000000);}
}
else	//must be non-null
{
for(k=0;k<num_kins;k++){fprintf(output2, "Her_K%d %.6f %.6f NA NA NA\n", k+1, hers[k], hersds[k]);}
for(r=0;r<num_regs;r++){fprintf(output2, "Her_R%d %.6f %.6f NA NA NA\n", r+1, hers[num_kins+r], hersds[num_kins+r]);}
fprintf(output2, "Her_Top %.6f NA NA NA NA\n", topher);
fprintf(output2, "Her_All %.6f %.6f NA NA NA\n", hers[total]+topher, hersds[total]);
}
fclose(output2);

sprintf(filename4,"%s.share.across", outfile);
if((output4=fopen(filename4,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename4);exit(1);}
fprintf(output4, "Component Share SE Expected Enrichment SE\n");
if(flag==0)
{
for(k=0;k<num_kins;k++){fprintf(output4, "Share_K%d %.6f %.6f %.6f %.6f %.6f\n", k+1, shares[k], sharesds[k], kinsums[k]/sum, shares[k]/kinsums[k]*sum, sharesds[k]/kinsums[k]*sum);}
for(r=0;r<num_regs;r++){fprintf(output4, "Share_R%d %.6f %.6f %.6f %.6f %.6f\n", r+1, shares[num_kins+r], sharesds[num_kins+r], Xsums[r]/sum, shares[num_kins+r]/Xsums[r]*sum, sharesds[num_kins+r]/Xsums[r]*sum);}
}
else
{
for(k=0;k<num_kins;k++){fprintf(output4, "Share_K%d %.6f %.6f NA NA NA\n", k+1, shares[k], sharesds[k]);}
for(r=0;r<num_regs;r++){fprintf(output4, "Share_R%d %.6f %.6f NA NA NA\n", r+1, shares[num_kins+r], sharesds[num_kins+r]);}
}
fclose(output4);

if(type==0&&prev!=-9999)	//binary he - already have factor
{
sprintf(filename5,"%s.he.across.liab", outfile);
if((output5=fopen(filename5,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename5);exit(1);}
fprintf(output5, "Num_Kinships %d\nNum_Regions %d\nNum_Top_Predictors %d\nNum_Covariates %d\nNum_Environments %d\n", num_kins, num_regs, num_tops, num_covars, num_envs);
fprintf(output5, "Coeffsfile %s.coeff.liab\nCovar_Heritability %.4f\n", outfile, covher*factor);
fprintf(output5, "Total_Samples %d\nWith_Phenotypes %d\n", ns, np);
fprintf(output5, "Null_Likelihood %.6f\nAlt_Likelihood %.6f\n", likenull, like);
if(total==1){fprintf(output5, "LRT_Stat %.4f\nLRT_P %.4e\n", lrtstat, lrtpva);}
else{fprintf(output5, "LRT_Stat %.4f\nLRT_P NA\n", lrtstat);}

fprintf(output5, "Component Heritability SE Size Mega_Intensity SE\n");
if(flag==0)
{
for(k=0;k<num_kins;k++){fprintf(output5, "Her_K%d %.6f %.6f %.2f %.6f %.6f\n", k+1, hers[k]*factor, hersds[k]*factor, kinsums[k], hers[k]/kinsums[k]*1000000*factor, hersds[k]/kinsums[k]*1000000*factor);}
for(r=0;r<num_regs;r++){fprintf(output5, "Her_R%d %.6f %.6f %.2f %.6f %.6f\n", r+1, hers[num_kins+r]*factor, hersds[num_kins+r]*factor, Xsums[r], hers[num_kins+r]/Xsums[r]*1000000*factor, hersds[num_kins+r]/Xsums[r]*1000000*factor);}
fprintf(output5, "Her_Top %.6f NA NA NA NA\n", topher*factor);
if(total==0){fprintf(output5, "Her_All %.6f NA NA NA NA\n", topher);}
else{fprintf(output5, "Her_All %.6f %.6f %.2f %.6f %.6f\n", hers[total]*factor+topher*factor, hersds[total]*factor, sum, hers[total]/sum*1000000*factor, hersds[total]/sum*1000000*factor);}
}
else
{
for(k=0;k<num_kins;k++){fprintf(output5, "Her_K%d %.6f %.6f NA NA NA\n", k+1, hers[k]*factor, hersds[k]*factor);}
for(r=0;r<num_regs;r++){fprintf(output5, "Her_R%d %.6f %.6f NA NA NA\n", r+1, hers[num_kins+r]*factor, hersds[num_kins+r]*factor);}
fprintf(output5, "Her_Top %.6f NA NA NA NA\n", topher*factor);
fprintf(output5, "Her_All %.6f %.6f NA NA NA\n", hers[total]*factor+topher*factor, hersds[total]*factor);
}
fclose(output5);
}

if(type==1)
{
sprintf(filename5,"%s.pcgc.across.marginal", outfile);
if((output5=fopen(filename5,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename5);exit(1);}
fprintf(output5, "Num_Kinships %d\nNum_Regions %d\nNum_Top_Predictors %d\nNum_Covariates %d\nNum_Environments %d\n", num_kins, num_regs, num_tops, num_covars, num_envs);
fprintf(output5, "Coeffsfile %s.coeff\nCovar_Heritability %.4f\n", outfile, covher);
fprintf(output5, "Total_Samples %d\nWith_Phenotypes %d\n", ns, np);
fprintf(output5, "Null_Likelihood %.6f\nAlt_Likelihood %.6f\n", likenull, like);
if(total==1){fprintf(output5, "LRT_Stat %.4f\nLRT_P %.4e\n", lrtstat, lrtpva);}
else{fprintf(output5, "LRT_Stat %.4f\nLRT_P NA\n", lrtstat);}

fprintf(output5, "Component Heritability SE Size Mega_Intensity SE\n");
if(flag==0)	//might have null model
{
for(k=0;k<num_kins;k++){fprintf(output5, "Her_K%d %.6f %.6f %.2f %.6f %.6f\n", k+1, hers[k]*(1-covher), hersds[k]*(1-covher), kinsums[k], hers[k]*(1-covher)/kinsums[k]*1000000, hersds[k]*(1-covher)/kinsums[k]*1000000);}
for(r=0;r<num_regs;r++){fprintf(output5, "Her_R%d %.6f %.6f %.2f %.6f %.6f\n", r+1, hers[num_kins+r]*(1-covher), hersds[num_kins+r]*(1-covher), Xsums[r], hers[num_kins+r]*(1-covher)/Xsums[r]*1000000, hersds[num_kins+r]*(1-covher)/Xsums[r]*1000000);}
fprintf(output5, "Her_Top %.6f NA NA NA NA\n", topher*(1-covher));
if(total==0){fprintf(output5, "Her_All %.6f NA NA NA NA\n", topher*(1-covher));}
else{fprintf(output5, "Her_All %.6f %.6f %.2f %.6f %.6f\n", hers[total]*(1-covher)+topher*(1-covher), hersds[total]*(1-covher), sum, hers[total]*(1-covher)/sum*1000000, hersds[total]*(1-covher)/sum*1000000);}
}
else	//must be non-null
{
for(k=0;k<num_kins;k++){fprintf(output5, "Her_K%d %.6f %.6f NA NA NA\n", k+1, hers[k]*(1-covher), hersds[k]*(1-covher));}
for(r=0;r<num_regs;r++){fprintf(output5, "Her_R%d %.6f %.6f NA NA NA\n", r+1, hers[num_kins+r]*(1-covher), hersds[num_kins+r]*(1-covher));}
fprintf(output5, "Her_Top %.6f NA NA NA NA\n", topher*(1-covher));
fprintf(output5, "Her_All %.6f %.6f NA NA NA\n", hers[total]*(1-covher)+topher*(1-covher), hersds[total]*(1-covher));
}
fclose(output5);
}

sprintf(filename6,"%s.cross.across", outfile);
if((output6=fopen(filename6,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename6);exit(1);}
for(k=0;k<num_kins;k++){fprintf(output6, "Her_K%d\t", k+1);}
for(r=0;r<num_regs;r++){fprintf(output6, "Her_R%d\t", r+1);}
fprintf(output6, "\n");
for(k=0;k<total;k++)
{
for(k2=0;k2<total;k2++)
{fprintf(output6, "%.6f\t", cohers[k+k2*total]);}
fprintf(output6, "\n");
}
fclose(output6);

////////

//test statistic appears to be (t1-t2)^2/(v1+v2) where tk is total genetic and vk its variance
lrtstat=pow(stats[0]-stats[2],2)/(pow(stats[1],2)+pow(stats[3],2));
lrtpva=erfc(pow(lrtstat,.5)*M_SQRT1_2);
printf("Test statistic: %.4f, Pvalue: %.4e\n\n", lrtstat, lrtpva);

if(type==0){sprintf(filename,"%s.he.compare", outfile);}
else{sprintf(filename,"%s.pcgc.compare", outfile);}
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fprintf(output, "Her_Within %.6f\nSE_Within %.6f\n", stats[0], stats[1]);
fprintf(output, "Her_Across %.6f\nSE_Across %.6f\n", stats[2], stats[3]);
fprintf(output, "LRT_Stat %.6f\nLRT_P %.4e\n", lrtstat, lrtpva);
fclose(output);
}	//end of num_subs>1

////////

if(type==0)
{
printf("Main results saved in %s.he", outfile);
if(num_subs>1){printf(", with subset versions saved in %s.he.within and %s.he.across", outfile, outfile);}
}
else
{
printf("Main results saved in %s.pcgc", outfile);
if(num_subs>1){printf(", with subset versions saved in %s.pcgc.within and %s.pcgc.across", outfile, outfile);}
}
printf("\n\n");

////////

free(indexer);free(indexer2);
free(thetas);free(thetasds);free(thetapvas);free(Yadj);
free(hers);free(hersds);
if(total>0)
{
free(shares);free(sharesds);free(cohers);free(cohers2);
free(Dtemp);free(DTD);free(DTS);free(DTDs);free(DTSs);free(DTD2);
for(p=0;p<num_blocks;p++){free(MDTD[p]);free(MDTS[p]);}
free(MDTD);free(MDTS);
free(DTDa);free(DTSa);free(DTDb);free(DTSb);
for(p=0;p<num_blocks;p++){free(MDTDa[p]);free(MDTSa[p]);free(MDTDb[p]);free(MDTSb[p]);}
free(MDTDa);free(MDTSa);free(MDTDb);free(MDTSb);
}
free(subs);free(block);free(jacks);
if(num_kins>0&&memsave==1){free(datatemp);}

}	//end of he_reg

///////////////////////////

