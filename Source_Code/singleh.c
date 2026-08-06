/*
Copyright 2026 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Single-predictor DRM regression - might permute (no top-preds, enviro, spa or sample weights) - must have dtype=1
//drm=1 - autosomes, drm=2 - chrx
//first compute revised phenotypes (absolute difference to genotype-specific median), then regress on genotypes
//for mmaa, maybe better to multiply residual variance by YinvVY, but this is non-trivial to compute
//if permuting, we permute predictors when reading (whereas normally we permute phenotype and covariates)

///////////////////////////

if(dougvar==1||dougvar==2){printf("Zero Her\n\n");}
if(dougvar==1||dougvar==3){printf("Raw phenotypes\n\n");}
if(dougvar==9){printf("Just her\n\n");}
if(dougvar2!=-9999){printf("Force her to %f\n", dougvar2);}

//passqc indicates whether doing any qc
passqc=(minmaf!=-9999||maxmaf!=-9999||minvar!=-9999||minobs!=-9999||mininfo!=-9999);

//set number of pedigree heritabilities
num_hers1=20;

//set parpreds - indicates which predictors are considered PAR
parpreds=malloc(sizeof(int)*data_length);
if(drm==1)  //all predictors are PRS
{
for(j=0;j<data_length;j++){parpreds[j]=1;}
}
else    //divide based on xpar
{
if(xpar==0) //no par
{
for(j=0;j<data_length;j++){parpreds[j]=0;}
}
if(xpar==1) //all par
{
for(j=0;j<data_length;j++){parpreds[j]=1;}
}
if(xpar==19) //pars contains snps outside 2699520 and 154931044
{
for(j=0;j<data_length;j++){parpreds[j]=(bp[j]<2699520||bp[j]>154931044);}
}
if(xpar==38) //pars contains snps outside 2781479 and 155701383
{
for(j=0;j<data_length;j++){parpreds[j]=(bp[j]<2781479||bp[j]>155701383);}
}

count=0;for(j=0;j<data_length;j++){count+=parpreds[j];}
printf("%d of the %d SNPs are determined to be within the PAR", count, data_length);
if(xpar==19){printf(" (defined as SNPs with basepairs below 2,699,520 or above 154,931,044)");}
if(xpar==38){printf(" (defined as SNPs with basepairs below 2,781,479 or above 155,701,383)");}
printf("\n\n");
}

if(drm==2&&xpar!=1)  //read sexes and get numbers of males and females
{
sexes=malloc(sizeof(int)*num_samples_use);
read_sexfile(sexfile, sexes, num_samples_use, ids1, ids2, ids3);
num_males=0;for(i=0;i<num_samples_use;i++){num_males+=(sexes[i]==1);}
num_females=0;for(i=0;i<num_samples_use;i++){num_females+=(sexes[i]==2);}

printf("There are %d males and %d females", num_males, num_females);
if(num_males<3||num_females<3){printf("note that the sex-specific analyses require at least three samples)");}
printf("\n\n");
}

//allocate variables

data_warn2(bitsize,2*num_samples_use);
data=malloc(sizeof(double)*num_samples_use*bitsize);
data2=malloc(sizeof(double)*num_samples_use*bitsize);

if(strcmp(relfile,"blank")==0){anal_warn2(bitsize,num_samples_use);}
else{anal_warn2(bitsize,5*num_samples_use);}

order=malloc(sizeof(int)*num_samples_use);
order2=malloc(sizeof(int)*num_samples_use);
Y=malloc(sizeof(double)*num_samples_use);
Z=malloc(sizeof(double)*num_samples_use*num_fixed);

keepsamps2=malloc(sizeof(int)*num_samples_use);
allids1=malloc(sizeof(char*)*num_samples_use);
allids2=malloc(sizeof(char*)*num_samples_use);
allids3=malloc(sizeof(char*)*num_samples_use);
retain=malloc(sizeof(int)*num_fixed);

thetas=malloc(sizeof(double)*num_fixed);
thetasds=malloc(sizeof(double)*num_fixed);
thetapvas=malloc(sizeof(double)*num_fixed);

Yadj=malloc(sizeof(double)*num_samples_use*bitsize);
Yadj2=malloc(sizeof(double)*num_samples_use);
Yadj3=malloc(sizeof(double)*num_samples_use);

stats=malloc(sizeof(double)*14*bitsize);

YTdata=malloc(sizeof(double)*2);
XTCX=malloc(sizeof(double)*4);
XTCX2=malloc(sizeof(double)*4);

nums=malloc(sizeof(int*)*bitsize);
for(j=0;j<bitsize;j++){nums[j]=malloc(sizeof(int)*3);}

svars=malloc(sizeof(double*)*bitsize);
for(j=0;j<bitsize;j++){svars[j]=malloc(sizeof(double)*3);}

Mscales=malloc(sizeof(double)*bitsize);

if(strcmp(relfile,"blank")!=0)
{
maxpairs=countrows(relfile);
total=num_hers1;
if(nmcmc>num_hers1){total=nmcmc;}

value=(double)(maxpairs+num_samples_use)/1024/1024/1024*24;
value2=(double)(nmcmc+2*total)/1024*num_samples_use/1024*8/1024;
if(value+value2){printf("Warning, to estimate heritability requires %.1f Gb; sorry, this can not be reduced\n\n", value+value2);}

firsts=malloc(sizeof(int)*(maxpairs+num_samples_use));
seconds=malloc(sizeof(int)*(maxpairs+num_samples_use));
relats=malloc(sizeof(double)*(maxpairs+num_samples_use));

tryhers=malloc(sizeof(double)*num_hers1);
toprats=malloc(sizeof(double)*num_hers1);
botrats=malloc(sizeof(double)*num_hers1);
polates=malloc(sizeof(double)*num_hers1);

gaussian=malloc(sizeof(double)*num_samples_use*nmcmc);
cX=malloc(sizeof(double)*num_samples_use*total);
cR=malloc(sizeof(double)*num_samples_use*total);
}

//deal with progress file
sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}

////////

//will perform up to three analyses (for all individual, males then females)
for(dflag=0;dflag<3;dflag++)
{
//flag indicates whether performing this analysis
if(drm==1)  //must have dflag=0 - always performing
{flag=1;}
else    //depends on number of parpreds, males, females and maybe phenotype
{
count=0;for(j=0;j<data_length;j++){count+=parpreds[j];}

if(dflag==0)    //perform if there are some parpreds
{
flag=(count>0);
if(flag==1)
{
if(drm==2){printf("-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --\n\n");}
printf("This analysis tests the %d PAR predictors using all samples\n\n", count);
fclose(output);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"This analysis tests the %d PAR predictors using all samples\n", count);
}
}

if(dflag==1)    //perform if there are some non-parpreds and males, and reduced phenotype non-trivial
{
flag=(count<data_length&&num_males>=3);

if(flag==1) //get variance of phenotype across males
{
sum=0;sumsq=0;indcount=0;
for(i=0;i<num_samples_use;i++)
{
if(sexes[i]==1){sum+=resp[i];sumsq+=pow(resp[i],2);indcount++;}
}
mean=sum/indcount;
var=sumsq/indcount-pow(mean,2);
if(var==0){printf("Warning, the phenotype is trivial when restricted to males\n\n");flag=0;}
}

if(flag==1)
{
printf("-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --\n\n");
printf("This analysis tests the %d non-PAR predictors using %d males\n\n", data_length-count, num_males);
fclose(output);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"This analysis tests the %d non-PAR predictors using %d males\n", data_length-count, num_males);
}
}

if(dflag==2)    //perform if there are some non-parpreds and females, and reduced phenotype non-trivial
{
flag=(count<data_length&&num_females>=3);

if(flag==1) //get variance of phenotype across females
{
sum=0;sumsq=0;indcount=0;
for(i=0;i<num_samples_use;i++)
{
if(sexes[i]==2){sum+=resp[i];sumsq+=pow(resp[i],2);indcount++;}
}
mean=sum/indcount;
var=sumsq/indcount-pow(mean,2);
if(var==0){printf("Warning, the phenotype is trivial when restricted to females\n\n");flag=0;}
}

if(flag==1)
{
printf("-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --\n\n");
printf("This analysis tests the %d non-PAR predictors using %d females\n\n", data_length-count, num_females);
fclose(output);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"This analysis tests the %d non-PAR predictors using %d females\n", data_length-count, num_females);
}
}
}

if(flag==1) //perform the analysis
{
//work out which individuals are being used
if(dflag==0)    //using all individuals
{
num_samples_use2=num_samples_use;
for(i=0;i<num_samples_use;i++)
{
copy_string(allids1,i,ids1[i]);
copy_string(allids2,i,ids2[i]);
copy_string(allids3,i,ids3[i]);
order[i]=i;
}
}
if(dflag==1)    //using only males
{
num_samples_use2=0;
for(i=0;i<num_samples_use;i++)
{
if(sexes[i]==1)
{
copy_string(allids1,num_samples_use2,ids1[i]);
copy_string(allids2,num_samples_use2,ids2[i]);
copy_string(allids3,num_samples_use2,ids3[i]);
order[num_samples_use2]=i;
num_samples_use2++;
}
}
}
if(dflag==2)    //using only females
{
num_samples_use2=0;
for(i=0;i<num_samples_use;i++)
{
if(sexes[i]==2)
{
copy_string(allids1,num_samples_use2,ids1[i]);
copy_string(allids2,num_samples_use2,ids2[i]);
copy_string(allids3,num_samples_use2,ids3[i]);
order[num_samples_use2]=i;
num_samples_use2++;
}
}
}

//get keepsamps2 (then maybe permute)
for(i=0;i<num_samples_use2;i++){keepsamps2[i]=keepsamps[order[i]];}
if(permute==1){permute_int(keepsamps2,num_samples_use2);}

//fill Y
for(i=0;i<num_samples_use2;i++){Y[i]=resp[order[i]];}

//load up first column of Z
for(i=0;i<num_samples_use2;i++){Z[i]=covar[order[i]];}
retain[0]=0;
num_fixed2=1;

//for subsequent columns, will only use if non-trivial
for(j=1;j<num_fixed;j++)
{
sum=0;sumsq=0;
for(i=0;i<num_samples_use2;i++){sum+=covar[order[i]+j*num_samples_use];sumsq+=pow(covar[order[i]+j*num_samples_use],2);}
mean=sum/num_samples_use2;
var=sumsq/num_samples_use2-pow(mean,2);
if(var>0)   //non-trivial
{
for(i=0;i<num_samples_use2;i++){Z[i+num_fixed2*num_samples_use2]=covar[order[i]+j*num_samples_use];}
retain[j]=num_fixed2;
num_fixed2++;
}
else
{
if(dflag==0){printf("Error BE33, please tell Doug (%d of %d)\n\n", j, num_fixed);exit(1);}
if(dflag==1){printf("Warning, Covariate %d is trivial when restricted to males, so has been excluded from the analysis\n", j);}
else{printf("Warning, Covariate %d is trivial when restricted to females, so has been excluded from the analysis\n", j);}
retain[j]=-1;
}
}
if(num_fixed2<num_fixed){printf("\n");}

//get phenotype order
dptrs=malloc(sizeof(struct sorting_double)*num_samples_use2);
for(i=0;i<num_samples_use2;i++){dptrs[i].value=Y[i];dptrs[i].index=i;}
qsort(dptrs, num_samples_use2, sizeof(struct sorting_double), compare_sorting_double);
for(i=0;i<num_samples_use2;i++){order2[i]=dptrs[i].index;}
free(dptrs);

//get median phenotype
value=.5*Y[order2[(num_samples_use2-1)/2]]+.5*Y[order2[num_samples_use2/2]];

//get absolute differences
for(i=0;i<num_samples_use2;i++){Yadj2[i]=fabs(Y[i]-value);}

if(dougvar==1||dougvar==3)  //just here for testing (reverts to original phenotypes)
{
for(i=0;i<num_samples_use2;i++){Yadj2[i]=Y[i];}
}

//solve null model, get thetas, adjust response and standardize
reg_covar_lin(Yadj2, Z, num_samples_use2, num_fixed2, 0, thetas, thetasds, thetapvas, Yadj3, 1, NULL, NULL);

//save coefficients
if(dflag==0){sprintf(filename2,"%s.coeff.all", outfile);}
if(dflag==1){sprintf(filename2,"%s.coeff.males", outfile);}
if(dflag==2){sprintf(filename2,"%s.coeff.females", outfile);}
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2, "Component\tEffect\tSE\tP\n");
fprintf(output2, "Intercept\t%.4e\t%.4e\t%.4e\n", thetas[0], thetasds[0], thetapvas[0]);
for(j=1;j<num_covars;j++)
{
if(retain[j]==-1){fprintf(output2, "Covariate_%d\tNA\tNA\tNA\n",j);}
else{fprintf(output2, "Covariate_%d\t%.4e\t%.4e\t%.4e\n",j, thetas[retain[j]], thetasds[retain[j]], thetapvas[retain[j]]);}
}
fclose(output2);

if(strcmp(relfile,"blank")!=0) //read related pairs and estimate null heritability
{
printf("Reading details for %d related pairs from %s\n", maxpairs, relfile);
num_rels=read_relfile(relfile, maxpairs, firsts, seconds, relats, NULL, num_samples_use2, allids3, 1);

//fill gaussian
for(g=0;g<nmcmc;g++)
{
for(i=0;i<num_samples_use2;i++){gaussian[(size_t)g*num_samples_use2+i]=rnorm_safe();}
}

//set heritability grid
value=pow(.8/0.001,1.0/(num_hers1-1));
for(k=0;k<num_hers1;k++){tryhers[k]=0.001*pow(value,k);}

//estimate denominator for each heritability in turn
for(k=0;k<num_hers1;k++)
{
//put noise into cR, and set cX to zero
for(g=0;g<nmcmc;g++)
{
for(i=0;i<num_samples_use2;i++){cX[(size_t)g*num_samples_use2+i]=0;}
for(i=0;i<num_samples_use2;i++){cR[(size_t)g*num_samples_use2+i]=gaussian[(size_t)g*num_samples_use2+i];}
}

//compute denominator
sparse_cgd(num_samples_use2, nmcmc, num_rels, firsts, seconds, relats, cX, cR, -9999, NULL, NULL, 0, NULL, 3, 0.000001, tryhers[k], NULL, nmcmc, gaussian, -9999, NULL, NULL, NULL, NULL, NULL, -9999, NULL, botrats+k, NULL, NULL);
}

//put num_hers1 copies of phenotype into cR, and set cX to zero
total=num_hers1;
for(k=0;k<num_hers1;k++)
{
for(i=0;i<num_samples_use2;i++){cX[(size_t)k*num_samples_use2+i]=0;}
for(i=0;i<num_samples_use2;i++){cR[(size_t)k*num_samples_use2+i]=Yadj3[i];}
}

//compute numerators
sparse_cgd(num_samples_use2, total, num_rels, firsts, seconds, relats, cX, cR, -9999, NULL, Yadj3, 0, NULL, 4, 0.000001, -9999, NULL, -9999, NULL, -9999, NULL, NULL, NULL, NULL, NULL, num_hers1, tryhers, toprats, NULL, NULL);

//get ratios, then interpolate
for(k=0;k<num_hers1;k++){polates[k]=toprats[k]/botrats[k];}
her=inter_hers(num_hers1, polates, tryhers, 1);
printf("The estimated heritability for the NULL model is %.4f\n", her);
if(her==0){printf("Warning, this is very low, so has been increased to 0.01\n");her=0.01;}
if(her>0.8){printf("Warning, this is very high, so has been reduced to 0.8\n");her=0.8;}
printf("\n");

if(dougvar==9){exit(1);}
if(dougvar==1||dougvar==2){her=0;}
if(dougvar2!=-9999){her=dougvar2;}

 //prepare for cholesky
cholmod_start (&chol_c);
chol_T=cholmod_allocate_triplet(num_samples_use2,num_samples_use2,num_rels,-1,CHOLMOD_REAL,&chol_c);
if(chol_T==NULL){printf("Failed to allocate triplet matrix with %d individuals and %d related pairs; please tell Doug\n\n", num_samples_use2, num_rels);exit(1);}

chol_B=cholmod_allocate_dense(num_samples_use2,bitsize,num_samples_use2,CHOLMOD_REAL, &chol_c);
chol_B2=cholmod_allocate_dense(num_samples_use2,bitsize,num_samples_use2,CHOLMOD_REAL, &chol_c);
chol_rhs=(double *) chol_B->x;
chol_rhs2=(double *) chol_B2->x;

//first load variance matrix into chol_T - remember we need K her + I (1-her)
chol_Ti = (int *) chol_T->i;
chol_Tj = (int *) chol_T->j;
chol_Tx = (double *) chol_T->x;

for(i=0;i<num_rels;i++)
{
if(firsts[i]==seconds[i]){value=relats[i]*her+1-her;}
else{value=relats[i]*her;}

if(firsts[i]>=seconds[i]){chol_Ti[i]=firsts[i];chol_Tj[i]=seconds[i];chol_Tx[i]=value;}
else{chol_Ti[i]=seconds[i];chol_Tj[i]=firsts[i];chol_Tx[i]=value;}
}
chol_T->nnz=num_rels;

//now transfer chol_T to chol_A (can then free chol_T)
chol_A=cholmod_triplet_to_sparse(chol_T, num_rels, &chol_c);
cholmod_free_triplet(&chol_T, &chol_c);

//prepare for decomposition
chol_L=cholmod_analyze(chol_A, &chol_c);

predicted_nnz=0;for(i=0;i<num_samples_use2;i++){predicted_nnz+=((int *) chol_L->ColCount)[i];}
value=(double)predicted_nnz*12/1024/1024/1024;
if(value>.1){printf("Warning, the Cholesky decomposition is predicted to use %.2f Gb (note that this is approximate, so the actual memory may be a few times higher)\n\n", value);}

//get cholesky (can then free chol_A)
info=cholmod_factorize(chol_A,chol_L,&chol_c);
if(info!=1){printf("Cholesky failed\n\n");exit(1);}
cholmod_free_sparse(&chol_A, &chol_c);
}

////////

//deal with on-the-fly file
if(dflag==0){sprintf(filename2,"%s.DRM.all", outfile);}
if(dflag==1){sprintf(filename2,"%s.DRM.males", outfile);}
if(dflag==2){sprintf(filename2,"%s.DRM.females", outfile);}
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2, "Chromosome\tSNP\tBasepair\tA1\tA2\tEffect_LIN\tSE_LIN\tZ_LIN\tP_LIN\t");
fprintf(output2, "Effect_ADD\tSE_ADD\tZ_ADD\tP_ADD\tEffect_DOM\tSE_DOM\tZ_DOM\tP_DOM\tCHISQ_GEN\tP_GEN\t");
fprintf(output2, "A1_Freq\tn\tn_A1A1\tn_A1A2\tn_A2A2\tn_NA\tVar_A1A1\tVar_A1A2\tVar_A2A2\n");

//ready for bit loop
count2=0;
bittotal=(data_length-1)/bitsize+1;
for(bit=0;bit<bittotal;bit++)
{
bitstart=bit*bitsize;
bitend=(bit+1)*bitsize;
if(bitend>data_length){bitend=data_length;}
bitlength=bitend-bitstart;

if(bit%10==0||exact==1)
{
printf("Performing DRM regression for Chunk %d of %d\n", bit+1, bittotal);
fclose(output);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Performing DRM regression for Chunk %d of %d\n", bit+1, bittotal);

fclose(output2);
if((output2=fopen(filename2,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename2);exit(1);}
}

//work out how many predictors in this bit
if(dflag==0)    //testing par predictors
{
count=0;for(j=bitstart;j<bitend;j++){count+=(parpreds[j]==1);}
}
else    //testing nonpar predictors
{
count=0;for(j=bitstart;j<bitend;j++){count+=(parpreds[j]==0);}
}

if(count>0) //testing this bit
{
//read data and compute statistics (but do not centre or set missing to mean) - must have dtype=1
current=read_data_fly(datafile, 1, data, NULL, num_samples_use2, keepsamps2, bitstart, bitend, keeppreds_use, NULL, -9999, num_samples, num_preds, -9999, -9999, -9999, NULL, missingvalue, -9999, -9999, -9999, maxthreads);
stand_data(data, centres+bitstart, mults+bitstart, sqdevs+bitstart, rates+bitstart, infos+bitstart, num_samples_use2, bitlength, missingvalue, 0, 0, -9999, NULL, 0);

//exclude predictors not being tested
for(j=bitstart;j<bitend;j++)
{
if(dflag==0&&parpreds[j]==0){mults[j]=-9999;}
if(dflag!=0&&parpreds[j]==1){mults[j]=-9999;}
}

//get genotype counts and variances
#pragma omp parallel for private(j, i, sum, sum2, sum3, sumsq, sumsq2, sumsq3) schedule(static)
for(j=0;j<bitlength;j++)
{
nums[j][0]=0;nums[j][1]=0;nums[j][2]=0;
sum=0;sum2=0;sum3=0;
sumsq=0;sumsq2=0;sumsq3=0;
for(i=0;i<num_samples_use2;i++)
{
if(data[(size_t)j*num_samples_use2+i]==0){nums[j][0]++;sum+=Y[i];sumsq+=pow(Y[i],2);}
if(data[(size_t)j*num_samples_use2+i]==1){nums[j][1]++;sum2+=Y[i];sumsq2+=pow(Y[i],2);}
if(data[(size_t)j*num_samples_use2+i]==2){nums[j][2]++;sum3+=Y[i];sumsq3+=pow(Y[i],2);}
}

if(nums[j][0]>0){svars[j][0]=sumsq/nums[j][0]-pow(sum/nums[j][0],2);}
else{svars[j][0]=0;}
if(nums[j][1]>0){svars[j][1]=sumsq2/nums[j][1]-pow(sum2/nums[j][1],2);}
else{svars[j][1]=0;}
if(nums[j][2]>0){svars[j][2]=sumsq3/nums[j][2]-pow(sum3/nums[j][2],2);}
else{svars[j][2]=0;}
}

if(dflag==1)    //replace heterogeneous genotypes with missing (must update counts, vars, centres and maybe mults)
{
for(j=0;j<bitlength;j++){count2+=(mults[bitstart+j]!=-9999&&nums[j][1]>0);}

#pragma omp parallel for private(j, i) schedule(dynamic)
for(j=0;j<bitlength;j++)
{
if(nums[j][1]>0)
{
for(i=0;i<num_samples_use2;i++)
{
if(data[(size_t)j*num_samples_use2+i]==1){data[(size_t)j*num_samples_use2+i]=missingvalue;}
}
nums[j][1]=0;
svars[j][1]=0;
centres[bitstart+j]=2.0*nums[j][2]/(nums[j][0]+nums[j][2]);
if(nums[j][0]==0||nums[j][2]==0){mults[bitstart+j]=-9999;}
}
}
}

if(passqc==1)	//perform additional qc
{
for(j=0;j<bitlength;j++)
{
if(mults[bitstart+j]!=-9999)
{
maf=centres[bitstart+j]/2+(centres[bitstart+j]>1)*(1-centres[bitstart+j]);
if(minmaf!=-9999&&maf<minmaf){mults[bitstart+j]=-9999;}
if(maxmaf!=-9999&&maf>maxmaf){mults[bitstart+j]=-9999;}
if(minvar!=-9999&&sqdevs[bitstart+j]<minvar){mults[bitstart+j]=-9999;}
if(minobs!=-9999&&rates[bitstart+j]<minobs){mults[bitstart+j]=-9999;}
if(mininfo!=-9999&&infos[bitstart+j]<mininfo){mults[bitstart+j]=-9999;}
}
}
}

//for each SNP, compute new phenotypes
#pragma omp parallel for private(j) schedule(static)
for(j=0;j<bitlength;j++)
{
if(mults[bitstart+j]!=-9999){abs_median(data+(size_t)j*num_samples_use2, Yadj+(size_t)j*num_samples_use2, Y, num_samples_use2, order2, missingvalue, nums[j][0], nums[j][1], nums[j][2]);}
else
{
for(i=0;i<num_samples_use2;i++){Yadj[(size_t)j*num_samples_use2+i]=0;}
}
}

if(dougvar==1||dougvar==3)  //just here for testing (reverts to original phenotypes)
{
for(j=0;j<bitlength;j++)
{
for(i=0;i<num_samples_use2;i++){Yadj[(size_t)j*num_samples_use2+i]=Y[i];}
}
}

/*
if(bit==0)  //save new phens (just for testing)
{
if(dflag==0){sprintf(filename3,"%s.check.all", outfile);}
if(dflag==1){sprintf(filename3,"%s.check.males", outfile);}
if(dflag==2){sprintf(filename3,"%s.check.females", outfile);}
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3,"FID\tIID\tRaw_Phen");
for(j=0;j<num_fixed2;j++){fprintf(output3, "\tCov_%d", j+1);}
for(j=0;j<10;j++)
{
if(j<bitlength){fprintf(output3,"\tSNP_%d\tAdjPhen_%d", j+1, j+1);}
}
fprintf(output3,"\n");

for(i=0;i<num_samples_use2;i++)
{
fprintf(output3,"%s\t%s\t%.6f", allids1[i], allids2[i], Y[i]);
for(j=0;j<num_fixed2;j++){fprintf(output3, "\t%f", Z[i+j*num_samples_use2]);}
for(j=0;j<10;j++)
{
if(j<bitlength){fprintf(output3,"\t%.6f\t%.6f", data[(size_t)j*num_samples_use2+i], Yadj[i+j*num_samples_use2]);}
}
fprintf(output3,"\n");
}
fclose(output3);
}
*/

//regress covariates out of the new phenotypes
reg_covar_matrix(Yadj, Z, num_samples_use2, bitlength, num_fixed2);

//standardize phenotypes, saving scaling
stand_matrix_nomiss_scale(Yadj, num_samples_use2, num_samples_use2, bitlength, Mscales);

//deal with missing data and fill data2
#pragma omp parallel for private(j,i,mean,mean2) schedule(static)
for(j=0;j<bitlength;j++)
{
if(mults[bitstart+j]!=-9999)
{
mean=(double)(nums[j][1]+2*nums[j][2])/(nums[j][0]+nums[j][1]+nums[j][2]);
mean2=(double)(nums[j][1])/(nums[j][0]+nums[j][1]+nums[j][2]);
for(i=0;i<num_samples_use2;i++)
{
if(data[(size_t)j*num_samples_use2+i]!=missingvalue){data2[(size_t)j*num_samples_use2+i]=(data[(size_t)j*num_samples_use2+i]==1);}
else{data[(size_t)j*num_samples_use2+i]=mean; data2[(size_t)j*num_samples_use2+i]=mean2;}
}
}
else
{
for(i=0;i<num_samples_use2;i++){data[(size_t)j*num_samples_use2+i]=0;data2[(size_t)j*num_samples_use2+i]=0;;}
}
}

//regress covariates out of data and data2
reg_covar_matrix(data, Z, num_samples_use2, bitlength, num_fixed2);
reg_covar_matrix(data2, Z, num_samples_use2, bitlength, num_fixed2);

if(strcmp(relfile,"blank")!=0)  //compute invV X and invV X2 (it is OK that data is zero for excluded predictors)
{
for(j=0;j<bitlength;j++)
{
for(i=0;i<num_samples_use2;i++){chol_rhs[(size_t)j*num_samples_use2+i]=data[(size_t)j*num_samples_use2+i];}
}
chol_B->ncol=bitlength;
chol_S=cholmod_solve(CHOLMOD_A,chol_L,chol_B, &chol_c);
chol_sol=(double *) chol_S->x;

for(j=0;j<bitlength;j++)
{
for(i=0;i<num_samples_use2;i++){chol_rhs2[(size_t)j*num_samples_use2+i]=data2[(size_t)j*num_samples_use2+i];}
}
chol_B2->ncol=bitlength;
chol_S2=cholmod_solve(CHOLMOD_A,chol_L,chol_B2, &chol_c);
chol_sol2=(double *) chol_S2->x;
}

//ready to compute test statistics
for(j=0;j<bitlength;j++)
{
mark=14*j;

if(mults[bitstart+j]!=-9999)
{
//1dof model
if(strcmp(relfile,"blank")==0)  //classical
{
//compute YTdata and XTCX
value3=ddot_(&num_samples_use2, data+(size_t)j*num_samples_use2, &one, Yadj+(size_t)j*num_samples_use2, &one);
value4=ddot_(&num_samples_use2, data+(size_t)j*num_samples_use2, &one, data+(size_t)j*num_samples_use2, &one);

//compute estimated effect YTCX/XTCX and its variance YCT(YC-XCvalue)/XTCX(n-nf-1)
value=value3/value4;
value2=(num_samples_use2/value4-pow(value,2))/(num_samples_use2-num_fixed2-1);
}
else    //mmaa
{
//compute YTdata and XTCX
value3=ddot_(&num_samples_use2, chol_sol+(size_t)j*num_samples_use2, &one, Yadj+(size_t)j*num_samples_use2, &one);
value4=ddot_(&num_samples_use2, chol_sol+(size_t)j*num_samples_use2, &one, data+(size_t)j*num_samples_use2, &one);

//compute estimated effect YTCX/XTCX and its variance 1/XTCX
value=value3/value4;
value2=pow(value4,-1);
}

//store results (remember to scale)
stats[0+mark]=value/Mscales[j];
stats[1+mark]=pow(value2,.5)/Mscales[j];
stats[2+mark]=stats[0+mark]/stats[1+mark];
stats[3+mark]=erfc(fabs(stats[2+mark])*M_SQRT1_2);

if(nums[j][0]>0&&nums[j][1]>0&&nums[j][2]>0) //can do 2dof test (note that we can use value3 and value4 from the 1dof model)
{
if(strcmp(relfile,"blank")==0)  //classical
{
//fill YTdata
YTdata[0]=value3;
YTdata[1]=ddot_(&num_samples_use2, data2+(size_t)j*num_samples_use2, &one, Yadj+(size_t)j*num_samples_use2, &one);

//fill XTCX
//XTCX[0]=ddot_(&num_samples_use2, data+(size_t)j*num_samples_use2, &one, data+(size_t)j*num_samples_use2, &one);
XTCX[0]=value4;
XTCX[1]=ddot_(&num_samples_use2, data+(size_t)j*num_samples_use2, &one, data2+(size_t)j*num_samples_use2, &one);
XTCX[3]=ddot_(&num_samples_use2, data2+(size_t)j*num_samples_use2, &one, data2+(size_t)j*num_samples_use2, &one);
}
else    //mmaa
{
//fill YTdata
YTdata[0]=value3;
YTdata[1]=ddot_(&num_samples_use2, chol_sol2+(size_t)j*num_samples_use2, &one, Yadj+(size_t)j*num_samples_use2, &one);

//fill XTCX
XTCX[0]=value4;
XTCX[1]=ddot_(&num_samples_use2, chol_sol+(size_t)j*num_samples_use2, &one, data2+(size_t)j*num_samples_use2, &one);
XTCX[3]=ddot_(&num_samples_use2, chol_sol2+(size_t)j*num_samples_use2, &one, data2+(size_t)j*num_samples_use2, &one);
}

//get determinant
value=XTCX[0]*XTCX[3]-pow(XTCX[1],2);
if(value==0)
{
printf("Error, inversion has failed for Predictor %s, please tell Doug\nGenotype counts %d %d %d (missing %d), XTCX is %f %f %f, YTdata %f %f\n\n", preds[bitstart+j], nums[j][0], nums[j][1], nums[j][2], num_samples_use2-nums[j][0]-nums[j][1]-nums[j][2], XTCX[0], XTCX[1], XTCX[3] ,YTdata[0], YTdata[1]);exit(1);
}

//get inverse
XTCX2[0]=XTCX[3]/value;
XTCX2[1]=-XTCX[1]/value;
XTCX2[3]=XTCX[0]/value;

//effect size estimates are inv(XTX) XTY
mean=XTCX2[0]*YTdata[0]+XTCX2[1]*YTdata[1];
mean2=XTCX2[1]*YTdata[0]+XTCX2[3]*YTdata[1];

if(strcmp(relfile,"blank")==0)
{
//estimate residual variance
value=(num_samples_use2-mean*YTdata[0]-mean2*YTdata[1])/(num_samples_use2-num_fixed2-2);

//effect size variances are inv(XTX) resvar
var=XTCX2[0]*value;
var2=XTCX2[3]*value;

//chisq(2) test statistic is t(means) XTY / residual variance
value4=(mean*YTdata[0]+mean2*YTdata[1])/value;
}
else
{
//effect size variances are inv(XTX)
var=XTCX2[0];
var2=XTCX2[3];

//chisq(2) test statistic is t(means) XTY
value4=(mean*YTdata[0]+mean2*YTdata[1]);
}

//convert test statistic to p-value, then to chisq(1) stat
value=pochisq(value4, 2);


if(value>0)	//have a valid p-value, so convert to chisq(1) statistic
{value2=pow(normal_inv(value/2),2);}
else	//not valid, so use approximation from wiki
//https://en.wikipedia.org/wiki/Noncentral_chi-squared_distribution#Cumulative_distribution_function
{
value3=pow(value4/2,1.0/3)-1+pow(4.5*2,-1);
value2=pow(value3,2)*4.5*2;
}

//store results
stats[4+mark]=mean/Mscales[j];
stats[5+mark]=pow(var,.5)/Mscales[j];
stats[6+mark]=stats[4+mark]/stats[5+mark];
stats[7+mark]=erfc(fabs(stats[6+mark])*M_SQRT1_2);
stats[8+mark]=mean2/Mscales[j];
stats[9+mark]=pow(var2,.5)/Mscales[j];
stats[10+mark]=stats[8+mark]/stats[9+mark];
stats[11+mark]=erfc(fabs(stats[10+mark])*M_SQRT1_2);
stats[12+mark]=value2;
stats[13+mark]=erfc(pow(value2,.5)*M_SQRT1_2);
}   //end of 2dof model
}
}   //end of j loop

//save results
for(j=0;j<bitlength;j++)
{
mark=14*j;

if(mults[bitstart+j]!=-9999)
{
fprintf(output2, "%d\t%s\t%.0f\t%s\t%s\t", chr[bitstart+j], preds[bitstart+j], bp[bitstart+j], along1[bitstart+j], along2[bitstart+j]);
fprintf(output2, "%.4e\t%.4e\t%.4f\t%.4e\t", stats[0+mark], stats[1+mark], stats[2+mark], stats[3+mark]);
if(nums[j][0]>0&&nums[j][1]>0&&nums[j][2]>0)
{
fprintf(output2, "%.4e\t%.4e\t%.4f\t%.4e\t", stats[4+mark], stats[5+mark], stats[6+mark], stats[7+mark]);
fprintf(output2, "%.4e\t%.4e\t%.4f\t%.4e\t%.4f\t%.4e\t", stats[8+mark], stats[9+mark], stats[10+mark], stats[11+mark], stats[12+mark], stats[13+mark]);
}
else{fprintf(output2, "%.4e\t%.4e\t%.4f\t%.4e\tNA\tNA\tNA\tNA\t%.4f\t%.4e\t", stats[0+mark], stats[1+mark], stats[2+mark], stats[3+mark], pow( stats[2+mark],2), stats[3+mark]);}

//add on frequency, counts and variances
fprintf(output2, "%.4f\t%d\t%d\t%d\t%d\t%d\t", centres[bitstart+j]/2, nums[j][2]+nums[j][1]+nums[j][0], nums[j][2], nums[j][1], nums[j][0], num_samples_use2-nums[j][2]-nums[j][1]-nums[j][0]);
fprintf(output2, "%.4e\t%.4e\t%.4e\n", svars[j][2], svars[j][1], svars[j][0]);
}
else	//did not test, but will include in assoc if not doing qc and was a potential test
{
if(passqc==0)
{
if((dflag==0&&parpreds[j]==1)||(dflag!=0&&parpreds[j]==0))
{
fprintf(output2, "%d\t%s\t%.0f\t%s\t%s\t", chr[bitstart+j], preds[bitstart+j], bp[bitstart+j], along1[bitstart+j], along2[bitstart+j]);
fprintf(output2, "NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t");
fprintf(output2, "%.4f\t%d\t%d\t%d\t%d\t%d\tNA\tNA\tNA\n", centres[bitstart+j]/2, nums[j][2]+nums[j][1]+nums[j][0], nums[j][2], nums[j][1], nums[j][0], num_samples_use2-nums[j][2]-nums[j][1]-nums[j][0]);
}
}
}
}	//end of j loop

if(strcmp(relfile,"blank")!=0) //must free chol_S and chol_S2
{
cholmod_free_dense(&chol_S, &chol_c);
cholmod_free_dense(&chol_S2, &chol_c);
}
}   //end of count>0
}	//end of bit loop
printf("\n");
if(count2>0){printf("Warning, %d predictors had heterogeneous genotypes for one or more samples (the heterogeneous genotypes were set to missing)\n\n", count2);}

fclose(output2);

if(dflag==0)    //tested par predictors
{
count=0;for(j=0;j<data_length;j++){count+=(parpreds[j]==1&&mults[j]==-9999);}
count2=0;for(j=0;j<data_length;j++){count2+=(parpreds[j]==1);}
}
else    //tested nonpar predictors
{
count=0;for(j=0;j<data_length;j++){count+=(parpreds[j]==0&&mults[j]==-9999);}
count2=0;for(j=0;j<data_length;j++){count2+=(parpreds[j]==0);}
}

if(count==count2)
{
if(passqc==0)
{printf("Warning, all %d predictors were excluded because they were trivial\n\n", count);}
else
{printf("Warning, all %d predictors were excluded because they were trivial or failed quality control\n\n", count);}
}
else
{
if(count>0)
{
if(passqc==0)
{printf("Warning, %d predictors were excluded because they were trivial\n\n", count);}
else
{printf("Warning, %d predictors were excluded because they were trivial or failed quality control\n\n", count);}
}
}

printf("Main results saved in %s\n\n", filename2);

for(i=0;i<num_samples_use2;i++){free(allids1[i]);free(allids2[i]);free(allids3[i]);}
}   //end of flag=1

if(drm==1){break;}
}   //end of dflag loop

fclose(output);

////////

free(parpreds);
if(drm==2&&xpar!=1){free(sexes);}

free(data);free(data2);
free(order);free(order2);free(Y);free(Z);
free(keepsamps2);free(allids1);free(allids2);free(allids3);free(retain);
free(thetas);free(thetasds);free(thetapvas);
free(Yadj);free(Yadj2);free(Yadj3);
free(stats);
free(YTdata);free(XTCX);free(XTCX2);
for(j=0;j<bitsize;j++){free(nums[j]);}free(nums);
for(j=0;j<bitsize;j++){free(svars[j]);}free(svars);
free(Mscales);
if(strcmp(relfile,"blank")!=0)
{
free(firsts);free(seconds);free(relats);
free(tryhers);free(toprats);free(botrats);free(polates);
free(gaussian);free(cX);free(cR);

chol_B->ncol=bitsize;
chol_B2->ncol=bitsize;
cholmod_free_dense(&chol_B, &chol_c);
cholmod_free_dense(&chol_B2, &chol_c);
cholmod_free_factor(&chol_L, &chol_c);
cholmod_finish (&chol_c);
}

///////////////////////////

