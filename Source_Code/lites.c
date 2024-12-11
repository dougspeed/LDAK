/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Construct prediction models from summary statistics - assumes phenotype has variance one (also standardize predictors)
//Have already filtered out predictors with mults=-9999

///////////////////////////

//neff is the average sample size for the main summaries
sum=0;for(j=0;j<data_length;j++){sum+=nss[j];}
neff=sum/data_length;

//work out which predictors remain within each block
usedpreds=malloc(sizeof(int)*num_preds);
for(j=0;j<num_preds;j++){usedpreds[j]=0;}
for(j=0;j<data_length;j++){usedpreds[keeppreds_use[j]]=1;}

blockuse=malloc(sizeof(int *)*bittotal);
count=0;
for(bit=0;bit<bittotal;bit++)
{
blockuse[bit]=malloc(sizeof(int)*(2+blocklengths[bit]));
blockuse[bit][0]=count;
count2=0;
for(j=0;j<blocklengths[bit];j++)
{
if(usedpreds[blockstarts[bit]+j]==1){blockuse[bit][2+count2]=j;count2++;}
}
blockuse[bit][1]=count2;
count+=count2;
}

free(usedpreds);

//can exclude blocks that contain no predictors
count=0;
for(bit=0;bit<bittotal;bit++)
{
if(blockuse[bit][1]>0)
{
if(count!=bit)
{
blockstarts[count]=blockstarts[bit];
blockends[count]=blockends[bit];
blocklengths[count]=blocklengths[bit];
blockindexes[count]=blockindexes[bit];
free(blockuse[count]);
blockuse[count]=malloc(sizeof(int)*(2+blockuse[bit][1]));
memcpy(blockuse[count],blockuse[bit],sizeof(int)*(2+blockuse[bit][1]));
}
count++;
}
}

if(count<bittotal){printf("Warning, due to predictor filtering, the number of windows has been reduced from %d to %d\n\n", bittotal, count);}
for(bit=count;bit<bittotal;bit++){free(blockuse[bit]);}
bittotal=count;

//get bitmax
bitmax=blockends[0]-blockstarts[0];
for(bit=1;bit<bittotal;bit++)
{
if(blockends[bit]-blockstarts[bit]>bitmax){bitmax=blockends[bit]-blockstarts[bit];}
}

//read centres, means, variances and rjksums from bin file (have already allocated the first three)
rjksums=malloc(sizeof(double)*data_length);

sprintf(filename2,"%s.cors.bin", corname);
if((input2=fopen(filename2,"rb"))==NULL)
{printf("Error opening %s\n\n",filename2);exit(1);}

//read means, allowing for filtering
fseeko(input2, 0, SEEK_SET);
current=0;
for(j=0;j<data_length;j++)
{
if(keeppreds_use[j]!=current)
{fseeko(input2, (off_t)sizeof(double)*keeppreds_use[j], SEEK_SET);}
if(fread(centres+j, sizeof(double), 1, input2)!=1)
{printf("Error reading mean for Predictor %d from %s\n\n", j+1, filename2);exit(1);}
current=keeppreds_use[j]+1;
}

//read scalings, allowing for filtering
fseeko(input2, (off_t)sizeof(double)*num_preds, SEEK_SET);
current=0;
for(j=0;j<data_length;j++)
{
if(keeppreds_use[j]!=current)
{fseeko(input2, (off_t)sizeof(double)*num_preds+sizeof(double)*keeppreds_use[j], SEEK_SET);}
if(fread(mults+j, sizeof(double), 1, input2)!=1)
{printf("Error reading scaling for Predictor %d from %s\n\n", j+1, filename2);exit(1);}
current=keeppreds_use[j]+1;
}

//read variances, allowing for filtering
fseeko(input2, (off_t)sizeof(double)*num_preds*2, SEEK_SET);
current=0;
for(j=0;j<data_length;j++)
{
if(keeppreds_use[j]!=current)
{fseeko(input2, (off_t)sizeof(double)*num_preds*2+sizeof(double)*keeppreds_use[j], SEEK_SET);}
if(fread(sqdevs+j, sizeof(double), 1, input2)!=1)
{printf("Error reading variance for Predictor %d from %s\n\n", j+1, filename2);exit(1);}
current=keeppreds_use[j]+1;
}

//read rjksums, allowing for filtering
fseeko(input2, (off_t)sizeof(double)*num_preds*3, SEEK_SET);
current=0;
for(j=0;j<data_length;j++)
{
if(keeppreds_use[j]!=current)
{fseeko(input2, (off_t)sizeof(double)*num_preds*3+sizeof(double)*keeppreds_use[j], SEEK_SET);}
if(fread(rjksums+j, sizeof(double), 1, input2)!=1)
{printf("Error reading tagging for Predictor %d from %s\n\n", j+1, filename2);exit(1);}
current=keeppreds_use[j]+1;
}

fclose(input2);

////////

if(skipcv==0)	//using pseudo cross-validation - get pseudo rhos and neff2, average size for training summaries
{
//only need rhos2 and rhos3 (not nss2 / nss3 or chis2 / chis3)
rhos2=malloc(sizeof(double)*data_length);
rhos3=malloc(sizeof(double)*data_length);

if(cvprop!=-9999)	//make them from cors.noise
{
//allocate and read datarands
datarands=malloc(sizeof(double)*data_length);
sprintf(filename,"%s.cors.noise", corname);
read_values(filename,datarands,data_length,keeppreds_use,1,0,0);

//rhos2 is rhos plus datarands*root(cvprop/(1-cvprop)/nss)
value=pow(cvprop/(1-cvprop),.5);
for(j=0;j<data_length;j++){rhos2[j]=rhos[j]+datarands[j]*value*pow(nss[j],-.5);}

//rhos3 is complement
for(j=0;j<data_length;j++){rhos3[j]=(rhos[j]-(1-cvprop)*rhos2[j])/cvprop;}

//set neff2
neff2=(1-cvprop)*neff;

free(datarands);
}
else	//read from pseudo files
{
//will also read nss2, nss3, chis2 and chis3 (then delete)
nss2=malloc(sizeof(double)*data_length);
nss3=malloc(sizeof(double)*data_length);
chis2=malloc(sizeof(double)*data_length);
chis3=malloc(sizeof(double)*data_length);

sprintf(filename2,"%s.train.summaries", pseudostem);
printf("Reading training summary statistics from %s\n", filename2);
read_sumsfile(filename2, nss2, chis2, rhos2, NULL, data_length, preds, al1, al2, bimfile, amb, -9999, 1.0, 0, -9999, 0);
printf("First few stats and ns are: %s %.3f %.1f", preds[0], chis2[0], nss2[0]);
for(j=1;j<data_length;j++){if(j<3){printf(" | %s %.3f %.1f", preds[j], chis2[j], nss2[j]);}}
printf("\n\n");

sprintf(filename3,"%s.test.summaries", pseudostem);
printf("Reading test summary statistics from %s\n", filename3);
read_sumsfile(filename3, nss3, chis3, rhos3, NULL, data_length, preds, al1, al2, bimfile, amb, -9999, 1.0, 0, -9999, 0);
printf("First few stats and ns are: %s %.3f %.1f", preds[0], chis3[0], nss3[0]);
for(j=1;j<data_length;j++){if(j<3){printf(" | %s %.3f %.1f", preds[j], chis3[j], nss3[j]);}}
printf("\n\n");

//set neff2
sum=0;for(j=0;j<data_length;j++){sum+=nss2[j];}
neff2=sum/data_length;

free(nss2);free(nss3);free(chis2);free(chis3);
}
}

if(skipcv==0)	//set highlds, which indicates which predictors to exclude when training models
{
highlds=malloc(sizeof(int)*data_length);
for(j=0;j<data_length;j++){highlds[j]=0;}

if(checkld==1)	//read predictors to exclude
{
if(strcmp(ldfile,"blank")==0)	//use file generated by calc-cors
{
sprintf(filename,"%s.cors.highld", corname);
count=countrows(filename);

head=check_head(filename,"None","None",0);

if(head==1){printf("There are no high-LD predictors (%s is empty)\n\n", filename);}
else
{
printf("Reading list of %d high-LD predictors from %s\n", count, filename);
wantpreds=malloc(sizeof(char*)*count);
indexer=malloc(sizeof(int)*data_length);
read_strings(filename, wantpreds, count, NULL, 1, 0);

count2=find_strings(preds, data_length, wantpreds, count, indexer, NULL, NULL, NULL, NULL, NULL, 3);
if(count2==0){printf("Warning, none of the predictors are in the data\n\n");}
if(count2>0&&count2<count){printf("Warning, only %d of these are in the data\n\n", count2);}
if(count2==count){printf("All of these are in the data\n\n");}

for(j=0;j<count2;j++){highlds[indexer[j]]=1;}

for(j=0;j<count;j++){free(wantpreds[j]);}free(wantpreds);free(indexer);
}
}
else	//exclude predictors in ldfile
{
count=countrows(ldfile);
printf("Reading list of %d high-LD predictors from %s\n", count, ldfile);
wantpreds=malloc(sizeof(char*)*count);
indexer=malloc(sizeof(int)*data_length);
read_strings(ldfile, wantpreds, count, NULL, 1, 0);

count2=find_strings(preds, data_length, wantpreds, count, indexer, NULL, NULL, NULL, NULL, NULL, 3);
if(count2==0){printf("Warning, none of the predictors are in the data\n\n");}
if(count2>0&&count2<count){printf("Warning, only %d of these are in the data\n\n", count2);}
if(count2==count){printf("All of these are in the data\n\n");}

for(j=0;j<count2;j++){highlds[indexer[j]]=1;}

for(j=0;j<count;j++){free(wantpreds[j]);}free(wantpreds);free(indexer);
}
}
}

if(her==-9999&&strcmp(fracfile,"blank")==0)	//obtain a sensible estimate of heritability using ldscores
{
sum=0;sum2=0;
for(j=0;j<data_length;j++){sum+=nss[j]/neff*(chis[j]-1);sum2+=pow(nss[j]/neff,2)*rjksums[j];}
her=sum/sum2*data_length/neff;

printf("Estimated heritability is %.4f\n", her);
if(her<0.01){printf("Warning, this is very low, so has been increased to 0.01\n");her=0.01;}
if(her>maxher){printf("Warning, this is very high, so has been reduced to %.4f\n", maxher);her=maxher;}
printf("\n");
}

//set and save parameters (will allocate trytypes, tryhers, trylams, tryscales, tryps, tryf2s)
#include "setparamsb.c"

if(shrink<1)	//scale tryscales
{
for(p=0;p<num_try;p++){tryscales[p]*=shrink;}
}

////////

//allocate variables

value=(double)bitmax/1024*bitmax/1024/1024*12;
if(value>1){printf("Warning, to store the correlations requires %.1f Gb\n\n", value);}

total=num_try;
if(prsvar==1&&num_blocks+1>total){total=num_blocks+1;}

if(prsvar==0){anal_warn(data_length, 3*num_try);}
else{anal_warn(data_length, 3*total+num_blocks);}

value=(double)bitmax/1024*bitmax/1024*12/1024;
if(value>1){printf("Warning, to store the correlations requires %.1f Gb\n\n", value);}

YTdata=malloc(sizeof(double)*data_length);
if(prsvar==1){YTdata2=malloc(sizeof(double)*data_length*num_blocks);}

exps=malloc(sizeof(double)*data_length);

lambdas=malloc(sizeof(double)*num_try);
lambdas2=malloc(sizeof(double)*num_try);
lambdas3=malloc(sizeof(double)*num_try);
lambdas4=malloc(sizeof(double)*num_try);

cors_single=malloc(sizeof(float)*bitmax*bitmax);
cors=malloc(sizeof(double)*bitmax*bitmax);
cors2=malloc(sizeof(double)*bitmax*total);

effs=malloc(sizeof(double)*data_length*total);
effs2=malloc(sizeof(double)*bitmax*total);
probs=malloc(sizeof(double)*data_length*num_try);
residuals=malloc(sizeof(double)*total);
ess=malloc(sizeof(double)*total);
ess2=malloc(sizeof(double)*total);

if(skipcv==0)
{
predvars=malloc(sizeof(double)*num_try);
predcors=malloc(sizeof(double)*num_try);
}

////////

if(prsvar==1)	//get jackknife statistics
{
neff3=(1-jackprop)*neff;

//XTY is neff3 x (rhos plus datarands2*root(jackprop/(1-jackprop)/nss)
value=pow(jackprop/(1-jackprop),.5);

readfloats=malloc(sizeof(float)*num_preds);

sprintf(filename,"%s.cors.jackknifes", corname);
if((input=fopen(filename,"rb"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}

for(p=0;p<num_blocks;p++)
{
if(fread(readfloats, sizeof(float), num_preds, input)!=num_preds)
{printf("Error reading jackknifes for Block %d from %s\n\n", p+1, filename);exit(1);}

for(j=0;j<data_length;j++)
{YTdata2[(size_t)p*data_length+j]=neff3*(rhos[j]+readfloats[keeppreds_use[j]]*value*pow(nss[j],-.5));}
}
fclose(input);

free(readfloats);
}

//fill exps, which contains expected per-predicted heritabilies (sum to one)

if(strcmp(indhers,"blank")!=0)	//already have
{
for(j=0;j<data_length;j++){exps[j]=weights[j];}
}
else
{
for(j=0;j<data_length;j++)
{
if(hwestand==1){exps[j]=weights[j]*pow(centres[j]*(1-centres[j]/2),1+power);}
else{exps[j]=weights[j]*pow(sqdevs[j],1+power);}
}

//make sure exps sum to one
sum=0;for(j=0;j<data_length;j++){sum+=exps[j];}
for(j=0;j<data_length;j++){exps[j]=exps[j]/sum;}
}

if(checkfreq==1)	//filter predictors based on frequency
{
printf("Calculating the differences in allele frequencies between the summary statistics and the correlations\n");
//set chisq threshold corresponding to P=1e-6
value=23.92813;

//get num_preds_use from rows 5 of root file
sprintf(filename, "%s.cors.root", corname);
read_integers(filename, &readint, 1, NULL, 2, 3, 0);

count=0;
for(j=0;j<data_length;j++)
{
//get the correlation frequency
value2=centres[j]/2;

//and the (weighted) mean frequency
value3=(nss[j]*a1freq[j]+readint*value2)/(nss[j]+readint);

//get the lrt stat
value4=4*(nss[j]*a1freq[j]*log(a1freq[j]/value3)+readint*value2*log(value2/value3)+nss[j]*(1-a1freq[j])*log((1-a1freq[j])/(1-value3))+readint*(1-value2)*log((1-value2)/(1-value3)));

if(value4>value){exps[j]=0;count++;}
}

if(count==0){printf("There are no significant differences, so no predictors are excluded\n\n");}
else{printf("Warning, the difference was significant for %d of the %d predictors, so these will be excluded\n\n", count, data_length);}
}

////////

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

if(trytypes[p]==4)	//bolt - exp(betaj^2) = plam + p2lam2, with p2lam2/plam = f2/(1-f2) - p and p2 within (0,1)
{
lambdas[p]=(1-tryf2s[p])/tryps[p];
lambdas2[p]=tryf2s[p]/tryp2s[p];
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

if(trytypes[p]==7)	//elastic - exp(betaj^2) = p/lam^2 + p2/lam2^2 + p3lam3
{
if(tryp3s[p]==0)	//lasso
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

if(trytypes[p]==8)	//ldpred (normal) - exp(betaj^2) = plam
{
lambdas[p]=pow(tryps[p],-1);
lambdas2[p]=0;
}

//if(trytypes[p]==9)	//ldpred (probs) - do not use lambdas
}	//end of p loop

//blank progress file
sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fclose(output);

///////////////////////

if(skipcv==0)	//solve using training summaries - var(Y)=1, so YTY=neff2 - then test using test summaries
{
if(num_try==1)	//trivial case - just need to set best to zero
{
printf("Can skip training phase because there is only set of parameters\n\n");

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Can skip training phase because there is only set of parameters\n");
fclose(output);

best=0;
}
else
{
//screen and file print
printf("Estimating effect sizes for %d models using pseudo training summary statistics\n", num_try);

if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Constructing %d models using pseudo training summary statistics\n", num_try);
fclose(output);

//get XTY for each predictor
for(j=0;j<data_length;j++){YTdata[j]=rhos2[j]*neff2;}

//set starting effects based on type of model
#pragma omp parallel for private(p,j,value,value2) schedule(dynamic, 1)
for(p=0;p<num_try;p++)
{
for(j=0;j<data_length;j++)
{
if(exps[j]>0&&highlds[j]==0)
{
//value and value2 are how much to scale gaussian and laplace parameters
value=exps[j]*tryhers[p];
value2=pow(exps[j]*tryhers[p],-.5);

if(trytypes[p]==1)	//lasso-sparse - set to zero
{effs[(size_t)p*data_length+j]=0;}
if(trytypes[p]==2)	//lasso - posterior mean, weighted by tagging
{effs[(size_t)p*data_length+j]=get_postmean(YTdata[j], lambdas[p]*value2, -9999, -9999, -9999, neff2, 1.0, -9999, -9999, -9999, -9999, NULL, 2, NULL)/rjksums[j];}
if(trytypes[p]==3)	//ridge - posterior mean, weighted by tagging
{effs[(size_t)p*data_length+j]=get_postmean(YTdata[j], lambdas[p]*value, -9999, -9999, -9999, neff2, 1.0, -9999, -9999, -9999, -9999, NULL, 3, NULL)/rjksums[j];}
if(trytypes[p]==4)	//bolt - posterior mean, weighted by tagging
{effs[(size_t)p*data_length+j]=get_postmean(YTdata[j], lambdas[p]*value, lambdas2[p]*value, -9999, -9999, neff2, 1.0, tryps[p], tryp2s[p], -9999, -9999, NULL, 4, NULL)/rjksums[j];}
if(trytypes[p]==5||trytypes[p]==6)	//bayesr or bayesr-shrink - posterior mean, weighted by tagging
{effs[(size_t)p*data_length+j]=get_postmean(YTdata[j], lambdas[p]*value, lambdas2[p]*value, lambdas3[p]*value, lambdas4[p]*value, neff2, 1.0, tryps[p], tryp2s[p], tryp3s[p], tryp4s[p], NULL, trytypes[p], NULL)/rjksums[j];}
if(trytypes[p]==7)	//elastic - posterior mean, weighted by tagging
{effs[(size_t)p*data_length+j]=get_postmean(YTdata[j], lambdas[p]*value2, lambdas2[p]*value2, lambdas3[p]*value, -9999, neff2, 1.0, tryps[p], tryp2s[p], tryp3s[p], -9999, NULL, 7, NULL)/rjksums[j];}

if(trytypes[p]==8)	//ldpred (normal) - posterior mean, weighted by tagging (use bolt function)
{effs[(size_t)p*data_length+j]=get_postmean(YTdata[j], lambdas[p]*value, lambdas2[p]*value, -9999, -9999, neff2, 1.0, tryps[p], tryp2s[p], -9999, -9999, NULL, 4, NULL)/rjksums[j];}
if(trytypes[p]==9)	//ldpred (probs) - bolt posterior mean, weighted by tagging
{
if(value<tryvars[p]){effs[(size_t)p*data_length+j]=get_postmean(YTdata[j], tryvars[p], 0.0, -9999, -9999, neff2, 1.0, value/tryvars[p], 1-value/tryvars[p], -9999, -9999, NULL, 4, NULL)/rjksums[j];}
else{effs[(size_t)p*data_length+j]=get_postmean(YTdata[j], tryvars[p], 0, -9999, -9999, neff2, 1.0, 1.0, 0.0, -9999, -9999, NULL, 4, NULL)/rjksums[j];}
}
}
//else	//set to zero
{effs[(size_t)p*data_length+j]=0;}

if(dougvar==1){effs[(size_t)p*data_length+j]=0;}
}
}

//set predvars to zero
for(p=0;p<num_try;p++){predvars[p]=0;}

//re-open cors
sprintf(filename2,"%s.cors.bin", corname);
if((input2=fopen(filename2,"rb"))==NULL)
{printf("Error opening %s\n\n",filename2);exit(1);}

ecount=0;
wcount=0;
for(bit=0;bit<bittotal;bit++)
{
bitstart=blockuse[bit][0];
bitend=bitstart+blockuse[bit][1];
bitlength=blockuse[bit][1];
bitlength2=blockends[bit]-blockstarts[bit];

if(bit%50==0)
{
printf("Estimating training effect sizes for Window %d of %d\n", bit+1, bittotal);

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Estimating training effect sizes for Window %d of %d\n", bit+1, bittotal);
fclose(output);
}

//read lower triangle correlations
fseeko(input2, (off_t)sizeof(double)*num_preds*4+sizeof(float)*blockindexes[bit], SEEK_SET);
for(k=0;k<bitlength2-1;k++)
{
count2=fread(cors_single+(size_t)k*bitlength2+k+1, sizeof(float), bitlength2-k-1, input2);
if(count2!=bitlength2-k-1)
{printf("Error reading correlations for Window %d from %s\n\n", bit+1, filename2);exit(1);}
}

//now extract correlations for predictors we are using
#pragma omp parallel for private(k,j) schedule(static)
for(k=0;k<bitlength;k++)
{
for(j=0;j<k;j++){cors[(size_t)k*bitlength+j]=cors_single[(size_t)blockuse[bit][2+j]*bitlength2+blockuse[bit][2+k]];}
cors[(size_t)k*bitlength+k]=1.0;
for(j=k+1;j<bitlength;j++){cors[(size_t)k*bitlength+j]=cors_single[(size_t)blockuse[bit][2+k]*bitlength2+blockuse[bit][2+j]];}
}

//save effect sizes for bit, in case fail to converge
#pragma omp parallel for private(p,j) schedule(static)
for(p=0;p<num_try;p++)
{
for(j=bitstart;j<bitend;j++){effs2[j-bitstart+p*bitmax]=effs[(size_t)p*data_length+j];}
}

//calculate ess for all models - equal to 2YTXbeta - t(beta) XTXbeta
alpha=1.0;beta=0.0;
dgemm_("N", "N", &bitlength, &num_try, &bitlength, &alpha, cors, &bitlength, effs+bitstart, &data_length, &beta, cors2, &bitlength);

#pragma omp parallel for private(p,j) schedule(static)
for(p=0;p<num_try;p++)
{
ess[p]=0;
for(j=bitstart;j<bitend;j++){ess[p]+=effs[(size_t)p*data_length+j]*(2*YTdata[j]/neff2-cors2[(j-bitstart)+p*bitlength]);}
}

for(count=0;count<maxiter;count++)
{
for(j=bitstart;j<bitend;j++)
{
if(exps[j]>0&&highlds[j]==0)
{
//get XjTXbeta for all models
alpha=1.0;beta=0.0;
dgemv_("T", &bitlength, &num_try, &alpha, effs+bitstart, &data_length, cors+(j-bitstart)*bitlength, &one, &beta, residuals, &one);

#pragma omp parallel for private(p,sum,value,value2,postmean) schedule(dynamic, 1)
for(p=0;p<num_try;p++)
{
//get XjT residuals
sum=YTdata[j]-tryscales[p]*neff2*(residuals[p]-effs[(size_t)p*data_length+j]);

//value and value2 are how much to scale gaussian and laplace parameters
value=exps[j]*tryhers[p];
value2=pow(exps[j]*tryhers[p],-.5);

//update effect size
if(trytypes[p]==1)	//lasso-sparse
{postmean=get_postmean(sum, lambdas[p]*value2*neff2, -9999, -9999, -9999, neff2, 1.0, -9999, -9999, -9999, -9999, NULL, 1, NULL);}
if(trytypes[p]==2)	//lasso
{postmean=get_postmean(sum, lambdas[p]*value2, -9999, -9999, -9999, neff2, 1.0, -9999, -9999, -9999, -9999, NULL, 2, probs+(size_t)p*data_length+j);}
if(trytypes[p]==3)	//ridge
{postmean=get_postmean(sum, lambdas[p]*value, -9999, -9999, -9999, neff2, 1.0, -9999, -9999, -9999, -9999, NULL, 3, probs+(size_t)p*data_length+j);}
if(trytypes[p]==4)	//bolt
{postmean=get_postmean(sum, lambdas[p]*value, lambdas2[p]*value, -9999, -9999, neff2, 1.0, tryps[p], tryp2s[p], -9999, -9999, NULL, 4, probs+(size_t)p*data_length+j);}
if(trytypes[p]==5||trytypes[p]==6)	//bayesr or bayesr-shrink
{postmean=get_postmean(sum, lambdas[p]*value, lambdas2[p]*value, lambdas3[p]*value, lambdas4[p]*value, neff2, 1.0, tryps[p], tryp2s[p], tryp3s[p], tryp4s[p], NULL, trytypes[p], probs+(size_t)p*data_length+j);}
if(trytypes[p]==7)	//elastic
{postmean=get_postmean(sum, lambdas[p]*value2, lambdas2[p]*value2, lambdas3[p]*value, -9999, neff2, 1.0, tryps[p], tryp2s[p], tryp3s[p], -9999, NULL, 7, probs+(size_t)p*data_length+j);}

if(trytypes[p]==8)	//ldpred (normal)
{postmean=get_postmean(sum, lambdas[p]*value, lambdas2[p]*value, -9999, -9999, neff2, 1.0, tryps[p], tryp2s[p], -9999, -9999, NULL, 4, NULL);}
if(trytypes[p]==9)	//ldpred (probs)
{
if(value<tryvars[p]){postmean=get_postmean(sum, tryvars[p], 0.0, -9999, -9999, neff2, 1.0, value/tryvars[p], 1-value/tryvars[p], -9999, -9999, NULL, 4, NULL);}
else{postmean=get_postmean(sum, tryvars[p], 0.0, -9999, -9999, neff2, 1.0, 1.0, 0.0, -9999, -9999, NULL, 4, NULL);}
}

effs[(size_t)p*data_length+j]=postmean;
}	//end of p loop
}}	//end of j loop

//save old ess
for(p=0;p<num_try;p++){ess2[p]=ess[p];}

//get new ess
alpha=1.0;beta=0.0;
dgemm_("N", "N", &bitlength, &num_try, &bitlength, &alpha, cors, &bitlength, effs+bitstart, &data_length, &beta, cors2, &bitlength);

#pragma omp parallel for private(p,j) schedule(static)
for(p=0;p<num_try;p++)
{
ess[p]=0;
for(j=bitstart;j<bitend;j++){ess[p]+=effs[(size_t)p*data_length+j]*(2*YTdata[j]/neff2-cors2[(j-bitstart)+p*bitlength]);}
}

//see whether converged
cflag=0;for(p=0;p<num_try;p++){cflag+=(fabs(ess[p]-ess2[p])<tol);}
if(cflag==num_try){break;}
}	//end of count loop

if(count==maxiter){printf("Warning, Window %d did not converge within %d iterations\n", bit+1, maxiter);}

cflag=0;for(p=0;p<num_try;p++){cflag+=(ess[p]>1);}
if(cflag>0)	//some estimates were suspect - revert all models to saved estimates
{
#pragma omp parallel for private(p,j) schedule(static)
for(p=start;p<end;p++)
{
for(j=bitstart;j<bitend;j++){effs[(size_t)p*data_length+j]=effs2[j-bitstart+p*bitmax];}
}

printf("Warning, Window %d went wild, its contribution will be excluded\n", bit+1);
ecount++;
wcount+=bitlength;
}

//recompute cors2, first setting small correlations to zero
for(k=0;k<bitlength;k++)
{
for(j=0;j<bitlength;j++)
{
if(fabs(cors[(size_t)k*bitlength+j])<0.01){cors[(size_t)k*bitlength+j]=0;}
}
}
alpha=1.0;beta=0.0;
dgemm_("N", "N", &bitlength, &num_try, &bitlength, &alpha, cors, &bitlength, effs+bitstart, &data_length, &beta, cors2, &bitlength);

//add on contribution to predvars
#pragma omp parallel for private(p,j) schedule(static)
for(p=0;p<num_try;p++)
{
for(j=bitstart;j<bitend;j++){predvars[p]+=effs[(size_t)p*data_length+j]*cors2[j-bitstart+p*bitlength];}
}
}	//end of bit loop
if(wcount==0){printf("\nCompleted: all windows converged\n\n");}
else{printf("\nCompleted:, %d windows failed to converge (in total, %d predictors)\n\n", wcount, ecount);}

fclose(input2);

////////

//test the models, and record the best (note that if exps[j]=0 or highld[j]=1, the predictor will have effect zero)
printf("Testing the models using pseudo test summary statistics\n");

//predcors = rhos3 beta / predvars
alpha=1.0;beta=0.0;
dgemv_("T", &data_length, &num_try, &alpha, effs, &data_length, rhos3, &one, &beta, predcors, &one);
for(p=0;p<num_try;p++)
{
if(predvars[p]>0){predcors[p]=predcors[p]*pow(predvars[p],-.5);}
else{predcors[p]=-9999;}
}

//find best
best=-1;
for(p=0;p<num_try;p++)
{
if(predcors[p]!=-9999)
{
if(best==-1){best=p;value=predcors[p];}
if(predcors[p]>value){best=p;value=predcors[p];}
}
}

if(best==-1)	//not possible to compute a correlation for any models
{printf("Error, it was not possible to compute a correlation for any of the models (suggesting they all have zero effect sizes)\n\n");exit(1);}

//save correlations
sprintf(filename2,"%s.cors",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2,"Model\tCorrelation\n");

for(p=0;p<num_try;p++)
{
if(predcors[p]!=-9999){fprintf(output2,"%d\t%.6f\n", p+1, predcors[p]);}
else{fprintf(output2,"%d\tNA\n", p+1);}
}
fclose(output2);

//save the best model
sprintf(filename3,"%s.best",outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3, "Model\tType\tHeritability\tlambda\tscale\tp\tf2\tp1\tp2\tp3\tp4\n");

if(trytypes[best]==1){fprintf(output3, "%d\tlasso-sparse\tNA\t%.4f\t%.4f\tNA\tNA\tNA\tNA\tNA\tNA\n", best+1, trylams[best], tryscales[best]);}
if(trytypes[best]==2){fprintf(output3, "%d\tlasso\t%.4f\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n", best+1, tryhers[best]);}
if(trytypes[best]==3){fprintf(output3, "%d\tridge\t%.4f\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n", best+1, tryhers[best]);}
if(trytypes[best]==4){fprintf(output3, "%d\tbolt\t%.4f\tNA\tNA\t%.4f\t%.4f\tNA\tNA\tNA\tNA\n", best+1, tryhers[best], tryps[best], tryf2s[best]);}
if(trytypes[best]==5){fprintf(output3, "%d\tbayesr\t%.4f\tNA\tNA\tNA\tNA\t%.4f\t%.4f\t%.4f\t%.4f\n", best+1, tryhers[best], tryps[best], tryp2s[best], tryp3s[best], tryp4s[best]);}
if(trytypes[best]==6){fprintf(output3, "%d\tbayesr-shrink\t%.4f\tNA\tNA\tNA\tNA\t%.4f\t%.4f\t%.4f\t%.4f\n", best+1, tryhers[best], tryps[best], tryp2s[best], tryp3s[best], tryp4s[best]);}
if(trytypes[best]==7){fprintf(output3, "%d\telastic\t%.4f\tNA\tNA\t%.4f\t%.4f\tNA\tNA\tNA\tNA\n", best+1, tryhers[best], 1-tryp3s[best], tryf2s[best]);}
if(trytypes[best]==8){fprintf(output3, "%d\tLDpred1\t%.4f\tNA\tNA\t%.4f\tNA\tNA\tNA\tNA\tNA\n", best+1, tryhers[best], tryps[best]);}
if(trytypes[best]==9){fprintf(output3, "%d\tLDpred2\t%.4f\tNA\tNA\t%.4f\tNA\tNA\tNA\tNA\tNA\n", best+1, tryhers[best], tryvars[best]);}
fclose(output3);

printf("The estimated accuracies are saved in %s, while the parameters corresponding to the best model (which has correlation %.2f) are saved in %s\n\n", filename2, predcors[best], filename3);
}	//end of not trivial
}	//end of using training and test summary statistics

////////

//solve for main statistics - var(Y)=1, so YTY=neff

//set start and end
if(skipcv==0)	//test only best model
{start=best;end=best+1;}
else	//test all models
{start=0;end=num_try;}

total2=end-start;

//screen and file print
if(skipcv==0){printf("Estimating effect sizes using the best parameters and summary statistics in %s\n", sumsfile);}
else{printf("Estimating effect sizes for %d models using summary statistics in %s\n", total2, sumsfile);}

if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Constructing %d models using summary statistics in %s\n", total2, sumsfile);
fclose(output);

//get XTY for each predictor
for(j=0;j<data_length;j++){YTdata[j]=rhos[j]*neff;}

//set starting effect based on type of model, and set probablities to zero
#pragma omp parallel for private(p,j,value,value2) schedule(dynamic, 1)
for(p=start;p<end;p++)
{
for(j=0;j<data_length;j++)
{
if(exps[j]>0)
{
//value and value2 are how much to scale gaussian and laplace parameters
value=exps[j]*tryhers[p];
value2=pow(exps[j]*tryhers[p],-.5);

if(trytypes[p]==1)	//lasso-sparse - set to zero
{effs[(size_t)p*data_length+j]=0;}
if(trytypes[p]==2)	//lasso - posterior mean, weighted by tagging
{effs[(size_t)p*data_length+j]=get_postmean(YTdata[j], lambdas[p]*value2, -9999, -9999, -9999, neff, 1.0, -9999, -9999, -9999, -9999, NULL, 2, NULL)/rjksums[j];}
if(trytypes[p]==3)	//ridge - posterior mean, weighted by tagging
{effs[(size_t)p*data_length+j]=get_postmean(YTdata[j], lambdas[p]*value, -9999, -9999, -9999, neff, 1.0, -9999, -9999, -9999, -9999, NULL, 3, NULL)/rjksums[j];}
if(trytypes[p]==4)	//bolt - posterior mean, weighted by tagging
{effs[(size_t)p*data_length+j]=get_postmean(YTdata[j], lambdas[p]*value, lambdas2[p]*value, -9999, -9999, neff, 1.0, tryps[p], tryp2s[p], -9999, -9999, NULL, 4, NULL)/rjksums[j];}
if(trytypes[p]==5||trytypes[p]==6)	//bayesr or bayesr-shrink - posterior mean, weighted by tagging
{effs[(size_t)p*data_length+j]=get_postmean(YTdata[j], lambdas[p]*value, lambdas2[p]*value, lambdas3[p]*value, lambdas4[p]*value, neff, 1.0, tryps[p], tryp2s[p], tryp3s[p], tryp4s[p], NULL, trytypes[p], NULL)/rjksums[j];}
if(trytypes[p]==7)	//elastic - posterior mean, weighted by tagging
{effs[(size_t)p*data_length+j]=get_postmean(YTdata[j], lambdas[p]*value2, lambdas2[p]*value2, lambdas3[p]*value, -9999, neff, 1.0, tryps[p], tryp2s[p], tryp3s[p], -9999, NULL, 7, NULL)/rjksums[j];}

if(trytypes[p]==8)	//ldpred (normal) - posterior mean, weighted by tagging (use bolt function)
{effs[(size_t)p*data_length+j]=get_postmean(YTdata[j], lambdas[p]*value, lambdas2[p]*value, -9999, -9999, neff, 1.0, tryps[p], tryp2s[p], -9999, -9999, NULL, 4, NULL)/rjksums[j];}
if(trytypes[p]==9)	//ldpred (probs) - bolt posterior mean, weighted by tagging
{
if(value<tryvars[p]){effs[(size_t)p*data_length+j]=get_postmean(YTdata[j], tryvars[p], 0.0, -9999, -9999, neff, 1.0, value/tryvars[p], 1-value/tryvars[p], -9999, -9999, NULL, 4, NULL)/rjksums[j];}
else{effs[(size_t)p*data_length+j]=get_postmean(YTdata[j], tryvars[p], 0, -9999, -9999, neff, 1.0, 1.0, 0.0, -9999, -9999, NULL, 4, NULL)/rjksums[j];}
}
}
else	//set to zero
{effs[(size_t)p*data_length+j]=0;}

if(dougvar==1){effs[(size_t)p*data_length+j]=0;}

probs[(size_t)p*data_length+j]=0;
}
}

//re-open cors
sprintf(filename2,"%s.cors.bin", corname);
if((input2=fopen(filename2,"rb"))==NULL)
{printf("Error opening %s\n\n",filename2);exit(1);}

ecount=0;
wcount=0;
for(bit=0;bit<bittotal;bit++)
{
bitstart=blockuse[bit][0];
bitend=bitstart+blockuse[bit][1];
bitlength=blockuse[bit][1];
bitlength2=blockends[bit]-blockstarts[bit];

if(bit%100==0)
{
printf("Estimating final effect sizes for Window %d of %d\n", bit+1, bittotal);

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Estimating final effect sizes for Window %d of %d\n", bit+1, bittotal);
fclose(output);
}

//read lower triangle correlations
fseeko(input2, (off_t)sizeof(double)*num_preds*4+sizeof(float)*blockindexes[bit], SEEK_SET);
for(k=0;k<bitlength2-1;k++)
{
if(fread(cors_single+(size_t)k*bitlength2+k+1, sizeof(float), bitlength2-k-1, input2)!=bitlength2-k-1)
{printf("Error reading correlations for Window %d from %s\n\n", bit+1, filename2);exit(1);}
}

//now extract correlations for predictors we are using
#pragma omp parallel for private(k,j) schedule(static)
for(k=0;k<bitlength;k++)
{
for(j=0;j<k;j++){cors[(size_t)k*bitlength+j]=cors_single[(size_t)blockuse[bit][2+j]*bitlength2+blockuse[bit][2+k]];}
cors[(size_t)k*bitlength+k]=1.0;
for(j=k+1;j<bitlength;j++){cors[(size_t)k*bitlength+j]=cors_single[(size_t)blockuse[bit][2+k]*bitlength2+blockuse[bit][2+j]];}
}

//save effect sizes for bit, in case fail to converge
#pragma omp parallel for private(p,j) schedule(static)
for(p=start;p<end;p++)
{
for(j=bitstart;j<bitend;j++){effs2[j-bitstart+p*bitmax]=effs[(size_t)p*data_length+j];}
}

//calculate ess for all models - equal to 2YTXbeta - t(beta) XTXbeta
alpha=1.0;beta=0.0;
dgemm_("N", "N", &bitlength, &total2, &bitlength, &alpha, cors, &bitlength, effs+(size_t)start*data_length+bitstart, &data_length, &beta, cors2+start*bitlength, &bitlength);

#pragma omp parallel for private(p,j) schedule(static)
for(p=start;p<end;p++)
{
ess[p]=0;
for(j=bitstart;j<bitend;j++){ess[p]+=effs[(size_t)p*data_length+j]*(2*YTdata[j]/neff-cors2[j-bitstart+p*bitlength]);}
}

for(count=0;count<maxiter;count++)
{
for(j=bitstart;j<bitend;j++)
{
if(exps[j]>0)
{
//get XjTXbeta for all models
alpha=1.0;beta=0.0;
dgemv_("T", &bitlength, &total2, &alpha, effs+(size_t)start*data_length+bitstart, &data_length, cors+(j-bitstart)*bitlength, &one, &beta, residuals+start, &one);

#pragma omp parallel for private(p,sum,value,value2,postmean) schedule(dynamic, 1)
for(p=start;p<end;p++)
{
//get XjT residuals
sum=YTdata[j]-tryscales[p]*neff*(residuals[p]-effs[(size_t)p*data_length+j]);

//value and value2 are how much to scale gaussian and laplace parameters
value=exps[j]*tryhers[p];
value2=pow(exps[j]*tryhers[p],-.5);

//update effect size
if(trytypes[p]==1)	//lasso-sparse
{postmean=get_postmean(sum, lambdas[p]*value2*neff, -9999, -9999, -9999, neff, 1.0, -9999, -9999, -9999, -9999, NULL, 1, NULL);}
if(trytypes[p]==2)	//lasso
{postmean=get_postmean(sum, lambdas[p]*value2, -9999, -9999, -9999, neff, 1.0, -9999, -9999, -9999, -9999, NULL, 2, probs+(size_t)p*data_length+j);}
if(trytypes[p]==3)	//ridge
{postmean=get_postmean(sum, lambdas[p]*value, -9999, -9999, -9999, neff, 1.0, -9999, -9999, -9999, -9999, NULL, 3, probs+(size_t)p*data_length+j);}
if(trytypes[p]==4)	//bolt
{postmean=get_postmean(sum, lambdas[p]*value, lambdas2[p]*value, -9999, -9999, neff, 1.0, tryps[p], tryp2s[p], -9999, -9999, NULL, 4, probs+(size_t)p*data_length+j);}
if(trytypes[p]==5||trytypes[p]==6)	//bayesr or bayesr-shrink
{postmean=get_postmean(sum, lambdas[p]*value, lambdas2[p]*value, lambdas3[p]*value, lambdas4[p]*value, neff, 1.0, tryps[p], tryp2s[p], tryp3s[p], tryp4s[p], NULL, trytypes[p], probs+(size_t)p*data_length+j);}
if(trytypes[p]==7)	//elastic
{postmean=get_postmean(sum, lambdas[p]*value2, lambdas2[p]*value2, lambdas3[p]*value, -9999, neff, 1.0, tryps[p], tryp2s[p], tryp3s[p], -9999, NULL, 7, probs+(size_t)p*data_length+j);}

if(trytypes[p]==8)	//ldpred (normal)
{postmean=get_postmean(sum, lambdas[p]*value, lambdas2[p]*value, -9999, -9999, neff, 1.0, tryps[p], tryp2s[p], -9999, -9999, NULL, 4, NULL);}
if(trytypes[p]==9)	//ldpred (probs)
{
if(value<tryvars[p]){postmean=get_postmean(sum, tryvars[p], 0.0, -9999, -9999, neff, 1.0, value/tryvars[p], 1-value/tryvars[p], -9999, -9999, NULL, 4, NULL);}
else{postmean=get_postmean(sum, tryvars[p], 0.0, -9999, -9999, neff, 1.0, 1.0, 0.0, -9999, -9999, NULL, 4, NULL);}
}

effs[(size_t)p*data_length+j]=postmean;
}	//end of p loop
}}	//end of j loop

//save old ess
for(p=start;p<end;p++){ess2[p]=ess[p];}

//get new ess
alpha=1.0;beta=0.0;
dgemm_("N", "N", &bitlength, &total2, &bitlength, &alpha, cors, &bitlength, effs+(size_t)start*data_length+bitstart, &data_length, &beta, cors2+start*bitlength, &bitlength);

#pragma omp parallel for private(p,j) schedule(static)
for(p=start;p<end;p++)
{
ess[p]=0;
for(j=bitstart;j<bitend;j++){ess[p]+=effs[(size_t)p*data_length+j]*(2*YTdata[j]/neff-cors2[j-bitstart+p*bitlength]);}
}

//see whether converged
cflag=0;for(p=start;p<end;p++){cflag+=(fabs(ess[p]-ess2[p])<tol);}
if(cflag==total2){break;}
}	//end of count loop

if(count==maxiter){printf("Warning, Window %d did not converge within %d iterations\n", bit+1, maxiter);}

cflag=0;for(p=start;p<end;p++){cflag+=(ess[p]>1);}
if(cflag>0)	//some estimates were suspect - revert all models to saved estimates
{
#pragma omp parallel for private(p,j) schedule(static)
for(p=start;p<end;p++)
{
for(j=bitstart;j<bitend;j++){effs[(size_t)p*data_length+j]=effs2[j-bitstart+p*bitmax];}
}

printf("Warning, Window %d went wild, its contribution will be excluded\n", bit+1);
ecount++;
wcount+=bitlength;
}
}	//end of bit loop
if(wcount==0){printf("Completed, all windows converged\n\n");}
else{printf("Warning, %d windows failed to converge (in total, %d predictors)\n\n", wcount, ecount);}

fclose(input2);

////////

if(prsvar==0)	//save all models
{
sprintf(filename2,"%s.effects",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}

fprintf(output2,"Predictor A1 A2 Centre");
for(p=start;p<end;p++){fprintf(output2," Model%d", p+1);}
fprintf(output2,"\n");
for(j=0;j<data_length;j++)
{
if(exps[j]>0)
{
fprintf(output2, "%s %c %c %.6f", preds[j], al1[j], al2[j], centres[j]);
for(p=start;p<end;p++){fprintf(output2," %.4e", effs[(size_t)p*data_length+j]*mults[j]);}
fprintf(output2,"\n");
}
}
fclose(output2);

sprintf(filename3,"%s.probs",outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3,"Predictor");
for(p=start;p<end;p++){fprintf(output3," Model%d", p+1);}
fprintf(output3,"\n");

for(j=0;j<data_length;j++)
{
fprintf(output3, "%s ", preds[j]);
for(p=start;p<end;p++){fprintf(output3," %.4e", probs[(size_t)p*data_length+j]);}
fprintf(output3,"\n");
}
fclose(output3);

if(total2==1){printf("Model saved in %s, with posterior probabilities in %s\n\n", filename2, filename3);}
else{printf("Models saved in %s, with posterior probabilities in %s\n\n", filename2, filename3);}
}

////////

if(prsvar==1)	//solve for jackknife statistics - var(Y)=1, so YTY=neff3 (must have skipcv=0 so will have found best)
{
//move best effects into end model
for(j=0;j<data_length;j++){effs[(size_t)num_blocks*data_length+j]=effs[(size_t)best*data_length+j];}

//screen and file print
printf("Estimating effect sizes using the best parameters and %d sets of jackknife summary statistics\n", num_blocks);

if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Constructing %d models using jackknife summary statistics\n",num_blocks);
fprintf(output, "Model\tType\tHeritability\tlambda\tscale\tp\tf2\tp1\tp2\tp3\tp4\tNum_Predictors_Failed\n");
fclose(output);

//set starting effect based on type of model
#pragma omp parallel for private(p,j,value,value2) schedule(dynamic, 1)
for(p=0;p<num_blocks;p++)
{
for(j=0;j<data_length;j++)
{
if(exps[j]>0)
{
//value and value2 are how much to scale gaussian and laplace parameters
value=exps[j]*tryhers[best];
value2=pow(exps[j]*tryhers[best],-.5);

if(trytypes[best]==1)	//lasso-sparse - set to zero
{effs[(size_t)p*data_length+j]=0;}
if(trytypes[best]==2)	//lasso - posterior mean, weighted by tagging
{effs[(size_t)p*data_length+j]=get_postmean(YTdata2[(size_t)p*data_length+j], lambdas[best]*value2, -9999, -9999, -9999, neff3, 1.0, -9999, -9999, -9999, -9999, NULL, 2, NULL)/rjksums[j];}
if(trytypes[best]==3)	//ridge - posterior mean, weighted by tagging
{effs[(size_t)p*data_length+j]=get_postmean(YTdata2[(size_t)p*data_length+j], lambdas[best]*value, -9999, -9999, -9999, neff3, 1.0, -9999, -9999, -9999, -9999, NULL, 3, NULL)/rjksums[j];}
if(trytypes[best]==4)	//bolt - posterior mean, weighted by tagging
{effs[(size_t)p*data_length+j]=get_postmean(YTdata2[(size_t)p*data_length+j], lambdas[best]*value, lambdas2[best]*value, -9999, -9999, neff3, 1.0, tryps[best], tryp2s[best], -9999, -9999, NULL, 4, NULL)/rjksums[j];}
if(trytypes[best]==5||trytypes[best]==6)	//bayesr or bayesr-shrink - posterior mean, weighted by tagging
{effs[(size_t)p*data_length+j]=get_postmean(YTdata2[(size_t)p*data_length+j], lambdas[best]*value, lambdas2[best]*value, lambdas3[best]*value, lambdas4[best]*value, neff3, 1.0, tryps[best], tryp2s[best], tryp3s[best], tryp4s[best], NULL, trytypes[best], NULL)/rjksums[j];}
if(trytypes[best]==7)	//elastic - posterior mean, weighted by tagging
{effs[(size_t)p*data_length+j]=get_postmean(YTdata2[(size_t)p*data_length+j], lambdas[best]*value2, lambdas2[best]*value2, lambdas3[best]*value, -9999, neff3, 1.0, tryps[best], tryp2s[best], tryp3s[best], -9999, NULL, 7, NULL)/rjksums[j];}

if(trytypes[best]==8)	//ldpred (normal) - posterior mean, weighted by tagging (use bolt function)
{effs[(size_t)p*data_length+j]=get_postmean(YTdata2[(size_t)p*data_length+j], lambdas[best]*value, lambdas2[best]*value, -9999, -9999, neff3, 1.0, tryps[best], tryp2s[best], -9999, -9999, NULL, 4, NULL)/rjksums[j];}
if(trytypes[best]==9)	//ldpred (probs) - bolt posterior mean, weighted by tagging
{
if(value<tryvars[best]){effs[(size_t)p*data_length+j]=get_postmean(YTdata2[(size_t)p*data_length+j], tryvars[best], 0.0, -9999, -9999, neff3, 1.0, value/tryvars[best], 1-value/tryvars[best], -9999, -9999, NULL, 4, NULL)/rjksums[j];}
else{effs[(size_t)p*data_length+j]=get_postmean(YTdata2[(size_t)p*data_length+j], tryvars[best], 0, -9999, -9999, neff3, 1.0, 1.0, 0.0, -9999, -9999, NULL, 4, NULL)/rjksums[j];}
}
}
else	//set to zero
{effs[(size_t)p*data_length+j]=0;}

if(dougvar==1){effs[(size_t)p*data_length+j]=0;}
}
}

//re-open cors
sprintf(filename2,"%s.cors.bin", corname);
if((input2=fopen(filename2,"rb"))==NULL)
{printf("Error opening %s\n\n",filename2);exit(1);}

ecount=0;
wcount=0;
for(bit=0;bit<bittotal;bit++)
{
bitstart=blockuse[bit][0];
bitend=bitstart+blockuse[bit][1];
bitlength=blockuse[bit][1];
bitlength2=blockends[bit]-blockstarts[bit];

if(bit%50==0)
{
printf("Estimating effect sizes for Window %d of %d\n", bit+1, bittotal);

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Estimating effect sizes for Window %d of %d\n", bit+1, bittotal);
fclose(output);
}

//read lower triangle correlations
fseeko(input2, (off_t)sizeof(double)*num_preds*4+sizeof(float)*blockindexes[bit], SEEK_SET);
for(k=0;k<bitlength2-1;k++)
{
if(fread(cors_single+(size_t)k*bitlength2+k+1, sizeof(float), bitlength2-k-1, input2)!=bitlength2-k-1)
{printf("Error reading correlations for Window %d from %s\n\n", bit+1, filename2);exit(1);}
}

//now extract correlations for predictors we are using
#pragma omp parallel for private(k,j) schedule(static)
for(k=0;k<bitlength;k++)
{
for(j=0;j<k;j++){cors[(size_t)k*bitlength+j]=cors_single[(size_t)blockuse[bit][2+j]*bitlength2+blockuse[bit][2+k]];}
cors[(size_t)k*bitlength+k]=1.0;
for(j=k+1;j<bitlength;j++){cors[(size_t)k*bitlength+j]=cors_single[(size_t)blockuse[bit][2+k]*bitlength2+blockuse[bit][2+j]];}
}

//save effect sizes for bit, in case fail to converge
#pragma omp parallel for private(p,j) schedule(static)
for(p=0;p<num_blocks;p++)
{
for(j=bitstart;j<bitend;j++){effs2[j-bitstart+p*bitmax]=effs[(size_t)p*data_length+j];}
}

//calculate ess for all models - equal to 2YTXbeta - t(beta) XTXbeta
alpha=1.0;beta=0.0;
dgemm_("N", "N", &bitlength, &num_blocks, &bitlength, &alpha, cors, &bitlength, effs+bitstart, &data_length, &beta, cors2, &bitlength);

#pragma omp parallel for private(p,j) schedule(static)
for(p=0;p<num_blocks;p++)
{
ess[p]=0;
for(j=bitstart;j<bitend;j++){ess[p]+=effs[(size_t)p*data_length+j]*(2*YTdata2[(size_t)p*data_length+j]/neff3-cors2[(j-bitstart)+p*bitlength]);}
}

for(count=0;count<maxiter;count++)
{
for(j=bitstart;j<bitend;j++)
{
if(exps[j]>0)
{
//get XjTXbeta for all models
alpha=1.0;beta=0.0;
dgemv_("T", &bitlength, &num_blocks, &alpha, effs+bitstart, &data_length, cors+(j-bitstart)*bitlength, &one, &beta, residuals, &one);

#pragma omp parallel for private(p,sum,value,value2,postmean) schedule(dynamic, 1)
for(p=0;p<num_blocks;p++)
{
//get XjT residuals
sum=YTdata2[(size_t)p*data_length+j]-tryscales[best]*neff3*(residuals[p]-effs[(size_t)p*data_length+j]);

//value and value2 are how much to scale gaussian and laplace parameters
value=exps[j]*tryhers[best];
value2=pow(exps[j]*tryhers[best],-.5);

//update effect size
if(trytypes[best]==1)	//lasso-sparse
{postmean=get_postmean(sum, lambdas[best]*value2*neff3, -9999, -9999, -9999, neff3, 1.0, -9999, -9999, -9999, -9999, NULL, 1, NULL);}
if(trytypes[best]==2)	//lasso
{postmean=get_postmean(sum, lambdas[best]*value2, -9999, -9999, -9999, neff3, 1.0, -9999, -9999, -9999, -9999, NULL, 2, NULL);}
if(trytypes[best]==3)	//ridge
{postmean=get_postmean(sum, lambdas[best]*value, -9999, -9999, -9999, neff3, 1.0, -9999, -9999, -9999, -9999, NULL, 3, NULL);}
if(trytypes[best]==4)	//bolt
{postmean=get_postmean(sum, lambdas[best]*value, lambdas2[best]*value, -9999, -9999, neff3, 1.0, tryps[best], tryp2s[best], -9999, -9999, NULL, 4, NULL);}
if(trytypes[best]==5||trytypes[best]==6)	//bayesr or bayesr-shrink
{postmean=get_postmean(sum, lambdas[best]*value, lambdas2[best]*value, lambdas3[best]*value, lambdas4[best]*value, neff3, 1.0, tryps[best], tryp2s[best], tryp3s[best], tryp4s[best], NULL, trytypes[best], NULL);}
if(trytypes[best]==7)	//elastic
{postmean=get_postmean(sum, lambdas[best]*value2, lambdas2[best]*value2, lambdas3[best]*value, -9999, neff3, 1.0, tryps[best], tryp2s[best], tryp3s[best], -9999, NULL, 7, NULL);}

if(trytypes[best]==8)	//ldpred (normal)
{postmean=get_postmean(sum, lambdas[best]*value, lambdas2[best]*value, -9999, -9999, neff3, 1.0, tryps[best], tryp2s[best], -9999, -9999, NULL, 4, NULL);}
if(trytypes[best]==9)	//ldpred (probs)
{
if(value<tryvars[best]){postmean=get_postmean(sum, tryvars[best], 0.0, -9999, -9999, neff3, 1.0, value/tryvars[best], 1-value/tryvars[best], -9999, -9999, NULL, 4, NULL);}
else{postmean=get_postmean(sum, tryvars[best], 0.0, -9999, -9999, neff3, 1.0, 1.0, 0.0, -9999, -9999, NULL, 4, NULL);}
}

effs[(size_t)p*data_length+j]=postmean;
}	//end of p loop
}}	//end of j loop

//save old ess
for(p=0;p<num_blocks;p++){ess2[p]=ess[p];}

//get new ess
alpha=1.0;beta=0.0;
dgemm_("N", "N", &bitlength, &num_blocks, &bitlength, &alpha, cors, &bitlength, effs+bitstart, &data_length, &beta, cors2, &bitlength);

#pragma omp parallel for private(p,j) schedule(static)
for(p=0;p<num_blocks;p++)
{
ess[p]=0;
for(j=bitstart;j<bitend;j++){ess[p]+=effs[(size_t)p*data_length+j]*(2*YTdata2[(size_t)p*data_length+j]/neff3-cors2[(j-bitstart)+p*bitlength]);}
}

//see whether converged
cflag=0;for(p=0;p<num_blocks;p++){cflag+=(fabs(ess[p]-ess2[p])<tol);}
if(cflag==num_blocks){break;}
}	//end of count loop

if(count==maxiter){printf("Warning, Window %d did not converge within %d iterations\n", bit+1, maxiter);}

cflag=0;for(p=0;p<num_blocks;p++){cflag+=(ess[p]>1);}
if(cflag>0)	//some estimates were suspect - revert all models to saved estimates
{
#pragma omp parallel for private(p,j) schedule(static)
for(p=0;p<num_blocks;p++)
{
for(j=bitstart;j<bitend;j++){effs[(size_t)p*data_length+j]=effs2[j-bitstart+p*bitmax];}
}

printf("Warning, Window %d went wild, its contribution will be excluded\n", bit+1);
ecount++;
wcount+=bitlength;
}
}	//end of bit loop
if(wcount==0){printf("Completed, all windows converged\n\n");}
else{printf("Warning, %d windows failed to converge (in total, %d predictors)\n\n", wcount, ecount);}

fclose(input2);

//save models
sprintf(filename2,"%s.effects",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}

fprintf(output2,"Predictor A1 A2 Centre Real");
for(p=0;p<num_blocks;p++){fprintf(output2," Jackknife_%.2f", jackprop);}
fprintf(output2,"\n");
for(j=0;j<data_length;j++)
{
if(exps[j]>0)
{
fprintf(output2, "%s %c %c %.6f %.4e", preds[j], al1[j], al2[j], centres[j], effs[(size_t)num_blocks*data_length+j]*mults[j]);
for(p=0;p<num_blocks;p++){fprintf(output2," %.4e", effs[(size_t)p*data_length+j]*mults[j]);}
fprintf(output2,"\n");
}
}
fclose(output2);

sprintf(filename3,"%s.probs",outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3,"Predictor Model%d\n", best+1);
fprintf(output3,"\n");

for(j=0;j<data_length;j++){fprintf(output3, "%s %.4e\n", preds[j], probs[(size_t)best*data_length+j]);}
fclose(output3);

printf("Models saved in %s, with posterior probabilities in %s\n\n", filename2, filename3);
}	//end of prsvar=1

////////

for(bit=0;bit<bittotal;bit++){free(blockuse[bit]);}free(blockuse);
free(rjksums);
if(skipcv==0){free(rhos2);free(rhos3);}
if(skipcv==0){free(highlds);}
free(trytypes);free(tryhers);free(trylams);free(tryscales);free(tryps);free(tryp2s);free(tryp3s);free(tryp4s);free(tryf2s);free(tryvars);
free(YTdata);
if(prsvar==1){free(YTdata2);}
free(exps);
free(lambdas);free(lambdas2);free(lambdas3);free(lambdas4);
free(cors_single);free(cors);free(cors2);
free(effs);free(effs2);free(probs);free(residuals);free(ess);free(ess2);
if(skipcv==0){free(predvars);free(predcors);}

///////////////////////////

