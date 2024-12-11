/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//LOCO regression - might permute, have top-preds and spa (no enviro or sample weights)
//performing linear regression (mode==131), quasi-logistic regression (mode==132&fastgwa==0) or logistic regression (mode==131&fastgwa!=0)
//for quasi-logistic spa, found best to regress (Y-mu) on X (instead of (Y-mu)/W on WX or alternatives)

///////////////////////////

//read loco.details for all phenotypes - want sample sizes, scaling and root of scaling - and check prs files

Pincs=malloc(sizeof(double)*num_resps_use);
Plambdas=malloc(sizeof(double)*num_resps_use);
Proots=malloc(sizeof(double)*num_resps_use);

for(m=0;m<num_resps_use;m++)
{
if(kvikparity==0){sprintf(filename,"%s.loco.details", prsfile);}
else{sprintf(filename,"%s.pheno%d.loco.details", prsfile, keepresps[m]+1);}
if(just_check(filename)!=0)
{
if(kvikstep==-9999&&gctastep==-9999&&faststep==-9999){printf("Error reading %s; this file would have been created using \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n", filename);}
else{printf("Error reading %s; this file would have been created using \"--kvik-step1\", \"--GCTA-LOCO-step1\" or \"--fastGWA-step1\"\n\n", filename);}
exit(1);
}

count=countrows(filename);
if(count!=6)
{printf("Error, %s should have six rows (not %d), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n", filename, count);exit(1);}

//open and get sample size
if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}
if(fscanf(input, "%s %d ", readstring, &readint)!=2)
{printf("Error reading Row 1 of %s\n\n", filename);exit(1);}
if(strcmp(readstring,"Actual_Sample_Size")!=0)
{printf("Error reading %s; Row 1 should begin \"Actual_Sample_Size\" (not %s), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n", filename, readstring);exit(1);}

//read effective sample size, and divide by actual sample size
if(fscanf(input, "%s %lf ", readstring, &readdouble)!=2)
{printf("Error reading Row 2 of %s\n\n", filename);exit(1);}
if(strcmp(readstring,"Approx_Effective_Sample_Size")!=0)
{printf("Error reading %s; Row 2 should begin \"Approx_Effective_Sample_Size\" (not %s), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n", filename, readstring);exit(1);}
Pincs[m]=readdouble/readint;

//read scaling, and also get square root
if(fscanf(input, "%s %lf ", readstring, Plambdas+m)!=2)
{printf("Error reading Row 3 of %s\n\n", filename);exit(1);}
if(strcmp(readstring,"Scaling_Estimate")!=0)
{printf("Error reading %s; Row 3 should begin \"Scaling_Estimate\" (not %s), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n", filename, readstring);exit(1);}
Proots[m]=pow(Plambdas[m],.5);

fclose(input);

if(kvikparity==0){sprintf(filename,"%s.loco.prs", prsfile);}
else{sprintf(filename,"%s.pheno%d.loco.prs", prsfile, keepresps[m]+1);}
if(just_check(filename)!=0)
{
if(kvikstep==-9999&&gctastep==-9999&&faststep==-9999){printf("Error reading %s; this file would have been created using \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n", filename);}
else{printf("Error reading %s; this file would have been created using \"--kvik-step1\", \"--GCTA-LOCO-step1\" or \"--fastGWA-step1\"\n\n", filename);}
exit(1);
}

count=countcols(filename);
if(fastgwa==0&&count!=2+num_chr2){printf("Error, %s should have %d columns (not %d), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n", filename, 2+num_chr2, count);exit(1);}
if(fastgwa!=0&&count!=3){printf("Error, %s should have three columns (not %d), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n", filename, count);exit(1);}

count2=countrows(filename);
if(count2!=readint+1)
{printf("Error, %s should have %d rows (not %d), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n", filename, readint+1, count2);exit(1);}
}	//end of m loop

////////

value=(double)num_samples_use/1024/1024/1024*num_chr2*num_resps_use;
if(value>1){printf("Warning, to store the PRS requires %.1f Gb; sorry, this can not be reduced\n\n", value);}

//deal with prs files

if(fastgwa==0)
{
chrindex=malloc(sizeof(int)*num_chr2);
prs=malloc(sizeof(double)*num_samples_use*num_chr2*num_resps_use);
}
else
{
pedprs=malloc(sizeof(double)*num_samples_use*num_resps_use);
}

if(fastgwa==0)	//open first prs file and check chromosomes
{
if(kvikparity==0){sprintf(filename,"%s.loco.prs", prsfile);}
else{sprintf(filename,"%s.pheno1.loco.prs", prsfile);}
if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}
if(fscanf(input, "%s %s ", readstring, readstring2)!=2)
{printf("Error reading first two elements of %s\n\n", filename);exit(1);}

//read chr numbers
for(j=0;j<num_chr2;j++)
{
if(fscanf(input, "Chr%d ", chrindex+j)!=1)
{printf("Error reading Element %d of %s\n\n", j+3, filename);exit(1);}
}
fclose(input);

//check all chr covered
usedpreds=malloc(sizeof(int)*(chrindex[num_chr2-1]+1));
for(j=0;j<chrindex[num_chr2-1]+1;j++){usedpreds[j]=0;}
for(j=0;j<num_chr2;j++){usedpreds[chrindex[j]]=1;}

for(j=0;j<data_length;j++)
{
if(chr[j]>chrindex[num_chr2-1])
{printf("Error, will be analyzing predictors on Chromosome %d, but this chromosome was not included when making the PRS\n\n", chr[j]);exit(1);}
if(usedpreds[chr[j]]==0)
{printf("Error, will be analyzing predictors on Chromosome %d, but this chromosome was not included when making the PRS\n\n", chr[j]);exit(1);}
}
free(usedpreds);
}

for(m=0;m<num_resps_use;m++)	//set prs values to zero
{
if(fastgwa==0)	//have values for each chromosome
{
for(j=0;j<num_chr2;j++)
{
for(i=0;i<num_samples_use;i++){prs[i+(j+m*num_chr2)*num_samples_use]=0;}
}
}
else	//have just one set of values
{
for(i=0;i<num_samples_use;i++){pedprs[i+m*num_samples_use]=0;}
}
}

for(m=0;m<num_resps_use;m++)	//check then read each PRS in turn
{
if(kvikparity==0){sprintf(filename,"%s.loco.prs", prsfile);}
else{sprintf(filename,"%s.pheno%d.loco.prs", prsfile, keepresps[m]+1);}

count=countrows(filename)-1;
wantids=malloc(sizeof(char*)*count);
indexer=malloc(sizeof(int)*num_samples_use);
indexer2=malloc(sizeof(int)*num_samples_use);
usedids=malloc(sizeof(int)*num_samples_use);

printf("Reading LOCO PRS for %d samples from %s\n", count, filename);

read_ids(filename, NULL, NULL, wantids, count, NULL, 1, 0);
count2=find_strings(wantids, count, ids3, num_samples_use, indexer, indexer2, NULL, NULL, NULL, NULL, 3);
if(count2==0){printf("Error, none of the %d samples in %s are in %s\n\n", count, filename, idsfile);exit(1);}

//check we have prs for all samples with non-missing phenotypes
for(i=0;i<num_samples_use;i++){usedids[i]=(resp[i+m*num_samples_use]!=missingvalue);}
for(i=0;i<count2;i++){usedids[indexer2[i]]+=2;}

//check overlap equals respcounts[m]
count3=0;for(i=0;i<num_samples_use;i++){count3+=(usedids[i]==3);}
if(count3<respcounts[m])
{
if(respcounts[m]==num_samples_use){printf("Error, %s contains PRS for only %d of the %d samples\n\n", filename, count3, num_samples_use);exit(1);}
else{printf("Error, %s contains PRS for only %d of the %d samples with non-missing phenotypes\n\n", filename, count3, respcounts[m]);exit(1);}
}

//check whether extra samples
count3=0;for(i=0;i<num_samples_use;i++){count3+=(usedids[i]==2);}
if(count3>0){printf("Warning, %s contains PRS for %d extra samples (ideally, the individuals used now should be the same as those used with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n", filename, count3);}

//open PRS and skip header row
if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}

found=0;
for(i=0;i<count;i++)
{
if(indexer[found]==i)	//using this row
{
if(fastgwa==0)	//have values for each chromosome
{
if(fscanf(input, "%s %s ", readstring, readstring2)!=2)
{printf("Error reading first two elements of Row %d of %s\n\n", i+2, filename);exit(1);}
for(j=0;j<num_chr2;j++)
{
if(fscanf(input, "%lf ", prs+indexer2[found]+(j+m*num_chr2)*num_samples_use)!=1)
{printf("Error reading PRS %d from Row %d of %s\n\n", j+1, i+1, filename);exit(1);}
}
}
else	//have just one value
{
if(fscanf(input, "%s %s %lf ", readstring, readstring2, pedprs+indexer2[found]+m*num_samples_use)!=3)
{printf("Error reading Row %d of %s\n\n", i+2, filename);exit(1);}
}
}
else	//skip this row
{
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}
}

found++;
if(found==count2){break;}
}

fclose(input);

for(i=0;i<count;i++){free(wantids[i]);}free(wantids);free(indexer);free(indexer2);free(usedids);
}
printf("\n");

////////

if(dougvar2!=-9999)
{
for(m=0;m<num_resps_use;m++)
{
printf("%d Original scaling is %f divide by %f\n", m+1, Plambdas[m], 1.0/Plambdas[m]);
Plambdas[m]=dougvar2;
Proots[m]=pow(Plambdas[m],.5);
printf("%d New scaling is %f divide by %f\n", m+1, Plambdas[m], 1.0/Plambdas[m]);
}
}

if(dougvar2==99)	//blank prs and set scalings to one
{
printf("blanking prs\n");
for(m=0;m<num_resps_use;m++)
{
if(fastgwa==0)	//have values for each chromosome
{
for(j=0;j<num_chr2;j++)
{
for(i=0;i<num_samples_use;i++){prs[i+(j+m*num_chr2)*num_samples_use]=0;}
}
}
else	//have just one set of values
{
for(i=0;i<num_samples_use;i++){pedprs[i+m*num_samples_use]=0;}
}

Plambdas[m]=1;Proots[m]=1;
}

printf("Revised scaling is one - this does not work for binary\n");	
}

////////

//passqc indicates whether doing any qc
passqc=(minmaf!=-9999||maxmaf!=-9999||minvar!=-9999||minobs!=-9999||mininfo!=-9999);

//threshold is the p-value threshold used with adjpreds
threshold=0.05;
if(spatest==1&&spathresh<threshold){threshold=spathresh;}

//allocate variables

data_warn2(bitsize,num_samples_use);
data=malloc(sizeof(double)*num_samples_use*bitsize);
datasqs=malloc(sizeof(double)*bitsize);

order=malloc(sizeof(int)*num_samples_use);
Y=malloc(sizeof(double)*num_samples_use*num_resps_use);
Z=malloc(sizeof(double)*num_samples_use*num_fixed);

Pscales=malloc(sizeof(double)*num_resps_use);
Pvarphens=malloc(sizeof(double)*num_resps_use);

thetas=malloc(sizeof(double)*num_fixed*num_resps_use);
thetasds=malloc(sizeof(double)*num_fixed*num_resps_use);
thetapvas=malloc(sizeof(double)*num_fixed*num_resps_use);
Yadj=malloc(sizeof(double)*num_samples_use*num_resps_use);
Yadj2=malloc(sizeof(double)*num_samples_use*num_resps_use);

if(mode==132)
{
nullprobs=malloc(sizeof(double)*num_samples_use*num_resps_use);
nullweights=malloc(sizeof(double)*num_samples_use*num_resps_use);
}

tindex=malloc(sizeof(int)*data_length);
stats=malloc(sizeof(double)*4*bitsize);
spastatus=malloc(sizeof(int)*bitsize);

YTdata=malloc(sizeof(double)*bitsize*num_resps_use);
XTCX=malloc(sizeof(double)*bitsize);

if(adjpreds==1)
{
Z2=malloc(sizeof(double)*num_samples_use*num_fixed*num_resps_use);
Z3=malloc(sizeof(double)*num_samples_use*num_fixed*num_resps_use);
ZTdata=malloc(sizeof(double)*num_fixed);
}

if(spatest==1)
{
usedpreds=malloc(sizeof(int)*bitsize);

if(mode==131||fastgwa==0)
{
Pspamax=malloc(sizeof(double)*num_resps_use);
knots=malloc(sizeof(double)*num_knots*num_resps_use);
bins=malloc(sizeof(double)*num_bins);
CGF0=malloc(sizeof(double*)*num_knots*num_resps_use);
CGF1=malloc(sizeof(double*)*num_knots*num_resps_use);
CGF2=malloc(sizeof(double*)*num_knots*num_resps_use);
CGF3=malloc(sizeof(double*)*num_knots*num_resps_use);

for(j=0;j<num_knots*num_resps_use;j++)
{
CGF0[j]=malloc(sizeof(double)*num_bins);
CGF1[j]=malloc(sizeof(double)*num_bins);
CGF2[j]=malloc(sizeof(double)*num_bins);
CGF3[j]=malloc(sizeof(double)*num_bins);
}
}
else
{
nullprobs2=malloc(sizeof(double)*num_samples_use*num_resps_use);
nullweights2=malloc(sizeof(double)*num_samples_use*num_resps_use);
}
}

//prepare for reading data
if(binary==0){open_datagz(&datainputgz, datafile, num_samples, genskip, genheaders, genprobs);}
current=0;

//will divide thetas and SDs by prop non-missing (this is on top of scaling due to PRS)
for(m=0;m<num_resps_use;m++)
{Pscales[m]=pow((double)respcounts[m]/num_samples_use,-1);}

//deal with order
for(i=0;i<num_samples_use;i++){order[i]=i;}
if(permute==1){permute_int(order,num_samples_use);}

//fill Y (maybe permuted)
for(m=0;m<num_resps_use;m++)
{
for(i=0;i<num_samples_use;i++){Y[i+m*num_samples_use]=resp[order[i]+m*num_samples_use];}
}

if(pad==1)	//pad missing
{impute_matrix_missing(Y, num_samples_use, num_samples_use, num_resps_use, missingvalue);}

//fill Z (maybe permuted)
for(j=0;j<num_fixed;j++)
{
for(i=0;i<num_samples_use;i++){Z[i+j*num_samples_use]=covar[order[i]+j*num_samples_use];}
}

for(m=0;m<num_resps_use;m++)	//solve null models
{
if(mode==131)	//get thetas, adjust response and get variance
{
reg_covar_lin(Y+m*num_samples_use, Z, num_samples_use, num_covars, num_tops, thetas+m*num_fixed, thetasds+m*num_fixed, thetapvas+m*num_fixed, Yadj+m*num_samples_use, 0, NULL, NULL);
}
else	//solve null model, get null weights, and compute W x adjusted response
{
reg_covar_log(Y+m*num_samples_use, Z, num_samples_use, num_covars, num_tops, offsets, thetas+m*num_fixed, thetasds+m*num_fixed, thetapvas+m*num_fixed, nullprobs+m*num_samples_use, 2, NULL, NULL, -9999, 0.001, 100);
for(i=0;i<num_samples_use;i++){nullweights[i+m*num_samples_use]=nullprobs[i+m*num_samples_use]*(1-nullprobs[i+m*num_samples_use]);}
for(i=0;i<num_samples_use;i++){Yadj[i+m*num_samples_use]=(Y[i+m*num_samples_use]-nullprobs[i+m*num_samples_use]);}
}

if(respcounts[m]<num_samples_use)	//adjust estimates for padding and scaling
{
for(j=0;j<num_fixed;j++)
{
thetas[j+m*num_fixed]=thetas[j+m*num_fixed]*Pscales[m]*Plambdas[m];
thetasds[j+m*num_fixed]=thetasds[j+m*num_fixed]*Pscales[m]*Proots[m];
thetapvas[j+m*num_fixed]=erfc(fabs(thetas[j+m*num_fixed]/thetasds[j+m*num_fixed])*M_SQRT1_2);
}
}

//save coefficients
if(kvikparity==0){sprintf(filename,"%s.coeff", outfile);}
else{sprintf(filename,"%s.pheno%d.coeff", outfile, keepresps[m]+1);}
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to\t%s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
if(mode==131){fprintf(output, "Component\tEffect\tSE\tP\n");}
else{fprintf(output, "Component\tLogOR\tSE\tP\n");}
fprintf(output, "Intercept\t%.4e\t%.4e\t%.4e\n", thetas[0+m*num_fixed], thetasds[0+m*num_fixed], thetapvas[0+m*num_fixed]);
for(j=1;j<num_covars;j++){fprintf(output, "Covariate_%d\t%.4e\t%.4e\t%.4e\n",j, thetas[j+m*num_fixed], thetasds[j+m*num_fixed], thetapvas[j+m*num_fixed]);}
fclose(output);
}	//end of m loop

//store original Yadj (saves needing to recompute for each chromosome, and used for spa initialization)
for(m=0;m<num_resps_use;m++)
{
for(i=0;i<num_samples_use;i++){Yadj2[i+m*num_samples_use]=Yadj[i+m*num_samples_use];}
}

//set tindex - indicates whether each predictor is a top preds (and if so, which)
for(j=0;j<data_length;j++){tindex[j]=-9999;}
if(num_tops>0)	//find the tops
{
mark=0;
for(j=0;j<num_tops;j++)
{
while(tkeeppreds[j]>keeppreds_use[mark]){mark++;}
tindex[mark]=j;
}
}

if(adjpreds==1)	//fill Z3 with either Z or WZ, then get Z2 = Z inv(ZTZ3) - per-phenotype values unnecessary for linear regression
{
for(m=0;m<num_resps_use;m++)
{
if(mode==131){copy_matrix(num_samples_use, num_fixed, Z, Z3+m*num_samples_use*num_fixed, 0, NULL);}
else{copy_matrix(num_samples_use, num_fixed, Z, Z3+m*num_samples_use*num_fixed, 1, nullweights+m*num_samples_use);}

get_Z2(Z2+m*num_samples_use*num_fixed, Z, Z3+m*num_samples_use*num_fixed, num_samples_use, num_fixed);
}
}

if(fastgwa==0)	//adjust responses for first chromosome PRS and get variances
{
mark2=0;while(chr[0]!=chrindex[mark2]){mark2++;}

for(m=0;m<num_resps_use;m++)
{
if(mode==131)	//adjusted phenotype is Yadj2 - PRS, get variance of Yadj2
{
for(i=0;i<num_samples_use;i++){Yadj[i+m*num_samples_use]=Yadj2[i+m*num_samples_use]-prs[i+(mark2+m*num_chr2)*num_samples_use];}
sum=0;sumsq=0;for(i=0;i<num_samples_use;i++){sum+=Yadj[i+m*num_samples_use];sumsq+=pow(Yadj[i+m*num_samples_use],2);}
mean=sum/num_samples_use;Pvarphens[m]=sumsq/num_samples_use-pow(mean,2);
}
else	//adjusted phenotype is Yadj2 - W PRS, get variance of 1/root(W) Yadj2
{
for(i=0;i<num_samples_use;i++){Yadj[i+m*num_samples_use]=Yadj2[i+m*num_samples_use]-prs[i+(mark2+m*num_chr2)*num_samples_use]*nullweights[i+m*num_samples_use];}
sum=0;sumsq=0;for(i=0;i<num_samples_use;i++){sum+=Yadj[i+m*num_samples_use]*pow(nullweights[i+m*num_samples_use],-.5);sumsq+=pow(Yadj[i+m*num_samples_use],2)/nullweights[i+m*num_samples_use];}
mean=sum/num_samples_use;Pvarphens[m]=sumsq/num_samples_use-pow(mean,2);
}
}
}
else	//similar but with pedigree PRS (and do not need variance for logistic regression)
{
for(m=0;m<num_resps_use;m++)
{
if(mode==131)
{
for(i=0;i<num_samples_use;i++){Yadj[i+m*num_samples_use]=Yadj2[i+m*num_samples_use]-pedprs[i+m*num_samples_use];}
sum=0;sumsq=0;for(i=0;i<num_samples_use;i++){sum+=Yadj[i+m*num_samples_use];sumsq+=pow(Yadj[i+m*num_samples_use],2);}
mean=sum/num_samples_use;Pvarphens[m]=sumsq/num_samples_use-pow(mean,2);
}
else
{
for(i=0;i<num_samples_use;i++){Yadj[i+m*num_samples_use]=Yadj2[i+m*num_samples_use]-pedprs[i+m*num_samples_use]*nullweights[i+m*num_samples_use];}
}
}
}

if(spatest==0)	//set spastatus to zero
{
for(j=0;j<bitsize;j++){spastatus[j]=0;}
}
else
{
if(mode==131||fastgwa==0)	//using empirical SPA
{
//set bins - evenly spaced from -2 to 2
for(k=0;k<num_bins;k++){bins[k]=-2+(double)k/(num_bins-1)*4;}

for(m=0;m<num_resps_use;m++)
{
if(spamax==-9999)	//set spamax based on variance of Yadj
{
sum=0;sumsq=0;for(i=0;i<num_samples_use;i++){sum+=Yadj2[i+m*num_samples_use];sumsq+=pow(Yadj2[i+m*num_samples_use],2);}
mean=sum/num_samples_use;var=sumsq/num_samples_use-pow(mean,2);
Pspamax[m]=2000*pow(num_samples_use*var,-.5);
}
else{Pspamax[m]=spamax;}

if(mpheno!=-1){printf("Computing empirical CGF (%d knots between %f and %f, and %d bins between -2 and 2)\n", num_knots, -Pspamax[m], Pspamax[m], num_bins);}
else{printf("Phenotype %d: computing empirical CGF (%d knots between %f and %f, and %d bins between -2 and 2)\n", m+1, num_knots, -Pspamax[m], Pspamax[m], num_bins);}

empirical_cumulants(Yadj+m*num_samples_use, num_samples_use, num_knots, knots+m*num_knots, num_bins, bins, CGF0+m*num_knots, CGF1+m*num_knots, CGF2+m*num_knots, CGF3+m*num_knots, Pspamax[m]);
}
printf("\n");
}
else	//using logistic spa - convenient to get some extra probs and weights
{
for(m=0;m<num_resps_use;m++)
{
for(i=0;i<num_samples_use;i++)
{
nullprobs2[i+m*num_samples_use]=nullprobs[i+m*num_samples_use]+pedprs[i+m*num_samples_use]*nullweights[i+m*num_samples_use];
if(nullprobs2[i+m*num_samples_use]<=0){nullprobs2[i+m*num_samples_use]=1e-10;}
if(nullprobs2[i+m*num_samples_use]>=1){nullprobs2[i+m*num_samples_use]=1-1e-10;}
nullweights2[i+m*num_samples_use]=nullprobs2[i+m*num_samples_use]*(1-nullprobs2[i+m*num_samples_use]);
}
}
}
}

////////

//deal with progress and on-the-fly files

if(kvikparity==0){printf("Main results saved in %s.assoc, with a summary version in %s.summaries and p-values in %s.pvalues\n\n", outfile, outfile, outfile);}
else
{
if(mpheno!=-1){printf("Main results saved in %s.pheno%d.assoc, with a summary version in %s.pheno%d.summaries and p-values in %s.pheno%d.pvalues\n\n", outfile, keepresps[0]+1, outfile, keepresps[0]+1, outfile, keepresps[0]+1);}
else{printf("Main results saved in %s.phenoX.assoc, with a summary version in %s.phenoX.summaries and p-values in %s.phenoX.pvalues, where X is the phenotype number\n\n", outfile, outfile, outfile);}
}

if(kvikparity==0||mpheno==-1){sprintf(filename,"%s.progress", outfile);}
else{sprintf(filename,"%s.pheno%d.progress", outfile, keepresps[0]+1);}
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}

for(m=0;m<num_resps_use;m++)
{
if(kvikparity==0){sprintf(filename2,"%s.assoc", outfile);}
else{sprintf(filename2,"%s.pheno%d.assoc", outfile, keepresps[m]+1);}
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2, "Chromosome\tPredictor\tBasepair\tA1\tA2\tWald_Stat\tWald_P\tEffect\tSE\tMAF\tCallRate\tMachR2\tSPA_Status\n");
fclose(output2);

if(kvikparity==0){sprintf(filename3,"%s.summaries", outfile);}
else{sprintf(filename3,"%s.pheno%d.summaries", outfile, keepresps[m]+1);}
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3, "Predictor\tA1\tA2\tZ\tn\tA1Freq\n");
fclose(output3);

if(kvikparity==0){sprintf(filename4,"%s.pvalues", outfile);}
else{sprintf(filename4,"%s.pheno%d.pvalues", outfile, keepresps[m]+1);}
if((output4=fopen(filename4,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename4);exit(1);}
fprintf(output4, "Predictor\tP\n");
fclose(output4);
}

//work out how many bits
bittotal=0;
bitend=0;
while(bitend<data_length)
{
bitstart=bitend;
bitend=bitstart+bitsize;
if(bitend>data_length){bitend=data_length;}
while(chr[bitend-1]>chr[bitstart]){bitend--;}

bittotal++;
}

//ready for bit loop
bit=0;
bitend=0;
while(bitend<data_length)
{
bitstart=bitend;
bitend=bitstart+bitsize;
if(bitend>data_length){bitend=data_length;}
while(chr[bitend-1]>chr[bitstart]){bitend--;}
bitlength=bitend-bitstart;

if(bit%10==0)
{
if(mode==131)
{
printf("Performing linear regression for Chunk %d of %d (Chromosome %d)\n", bit+1, bittotal, chr[bitstart]);
fclose(output);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Performing linear regression for Chunk %d of %d (Chromosome %d)\n", bit+1, bittotal, chr[bitstart]);
}
else
{
if(fastgwa==0)
{
printf("Performing quasi-logistic regression for Chunk %d of %d (Chromosome %d)\n", bit+1, bittotal, chr[bitstart]);
fclose(output);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Performing quasi-logistic regression for Chunk %d of %d (Chromosome %d)\n", bit+1, bittotal, chr[bitstart]);
}
else
{
printf("Performing logistic regression for Chunk %d of %d (Chromosome %d)\n", bit+1, bittotal, chr[bitstart]);
fclose(output);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Performing logistic regression for Chunk %d of %d (Chromosome %d)\n", bit+1, bittotal, chr[bitstart]);
}
}
}

if(fastgwa==0)	//see if changing prs (not necessary when using pedigree PRS)
{
if(chr[bitstart]!=chrindex[mark2])	//must update phenotype and variance, and maybe SPA
{
while(chr[bitstart]!=chrindex[mark2]){mark2++;}

for(m=0;m<num_resps_use;m++)
{
if(mode==131)
{
for(i=0;i<num_samples_use;i++){Yadj[i+m*num_samples_use]=Yadj2[i+m*num_samples_use]-prs[i+(mark2+m*num_chr2)*num_samples_use];}
sum=0;sumsq=0;for(i=0;i<num_samples_use;i++){sum+=Yadj[i+m*num_samples_use];sumsq+=pow(Yadj[i+m*num_samples_use],2);}
mean=sum/num_samples_use;Pvarphens[m]=sumsq/num_samples_use-pow(mean,2);
}
else
{
for(i=0;i<num_samples_use;i++){Yadj[i+m*num_samples_use]=Yadj2[i+m*num_samples_use]-prs[i+(mark2+m*num_chr2)*num_samples_use]*nullweights[i+m*num_samples_use];}
sum=0;sumsq=0;for(i=0;i<num_samples_use;i++){sum+=Yadj[i+m*num_samples_use]*pow(nullweights[i+m*num_samples_use],-.5);sumsq+=pow(Yadj[i+m*num_samples_use],2)/nullweights[i+m*num_samples_use];}
mean=sum/num_samples_use;Pvarphens[m]=sumsq/num_samples_use-pow(mean,2);
}

if(spatest==1)	//must be using canoical spa
{
if(mpheno!=-1){printf("Recomputing empirical CGF (%d knots between %f and %f, and %d bins between -2 and 2)\n", num_knots, -Pspamax[m], Pspamax[m], num_bins);}
else{printf("Phenotype %d: recomputing empirical CGF (%d knots between %f and %f, and %d bins between -2 and 2)\n", m+1, num_knots, -Pspamax[m], Pspamax[m], num_bins);}

empirical_cumulants(Yadj+m*num_samples_use, num_samples_use, num_knots, knots+m*num_knots, num_bins, bins, CGF0+m*num_knots, CGF1+m*num_knots, CGF2+m*num_knots, CGF3+m*num_knots, Pspamax[m]);
}
}
}
}

//read data, compute statistics, centre and set missing to zero
if(dtype==1)	//fast way
{
read_bed_wrapper(datafile, data, centres+bitstart, mults+bitstart, sqdevs+bitstart, rates+bitstart, infos+bitstart, num_samples_use, keepsamps, bitlength, keeppreds_use+bitstart, num_samples, num_preds, missingvalue, bedzeros, bedones, bedtwos, 1, maxthreads);
}
else	//slow way
{
current=read_data_fly(datafile, dtype, data, NULL, num_samples_use, keepsamps, bitstart, bitend, keeppreds_use, datainputgz, current, num_samples, num_preds, genskip, genheaders, genprobs, bgen_indexes, missingvalue, -9999, -9999, nonsnp, maxthreads);
stand_data(data, centres+bitstart, mults+bitstart, sqdevs+bitstart, rates+bitstart, infos+bitstart, num_samples_use, bitlength, missingvalue, 0, 0, -9999, NULL, 1);
}

if(passqc==1)	//perform qc (will already have mults=-9999 for trivial predictors)
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

//get XTY
alpha=1.0;beta=0.0;
dgemm_("T", "N", &bitlength, &num_resps_use, &num_samples_use, &alpha, data, &num_samples_use, Yadj, &num_samples_use, &beta, YTdata, &bitlength);

//ready to test
for(m=0;m<num_resps_use;m++)
{
if(m==0||mode==132)	//for linear regression, only need to get XTCX for the first phenotype
{
if(adjpreds==2)	//regress out covariates (its ok if have already done this for previous phenotypes)
{
if(mode==131)
{reg_covar_matrix(data, Z, num_samples_use, bitlength, num_fixed);}
else
{reg_covar_weighted(data, Z, num_samples_use, bitlength, num_fixed, nullweights+m*num_samples_use);}
}

//get XTCX
if(mode==131)
{
if(adjpreds==2)	//must compute
{
#pragma omp parallel for private(j, i) schedule(static)
for(j=0;j<bitlength;j++)
{
XTCX[j]=0;for(i=0;i<num_samples_use;i++){XTCX[j]+=pow(data[(size_t)j*num_samples_use+i],2);}
}
}
else	//already have
{
for(j=0;j<bitlength;j++){XTCX[j]=num_samples_use*sqdevs[bitstart+j];}
}
}
else	//must compute
{
#pragma omp parallel for private(j, i) schedule(static)
for(j=0;j<bitlength;j++)
{
XTCX[j]=0;for(i=0;i<num_samples_use;i++){XTCX[j]+=pow(data[(size_t)j*num_samples_use+i],2)*nullweights[i+m*num_samples_use];}
}
}
}

for(j=0;j<bitlength;j++)
{
mark=4*j;

if(mults[bitstart+j]!=-9999)	//will be testing - remember to scale theta and SD by Pscales and Plambdas / Proots
{
if(tindex[bitstart+j]!=-9999)	//a top predictor, so have already tested (and scaled)
{
stats[0+mark]=thetas[num_covars+tindex[bitstart+j]+m*num_fixed];
stats[1+mark]=thetasds[num_covars+tindex[bitstart+j]+m*num_fixed];
stats[2+mark]=stats[0+mark]/stats[1+mark];
stats[3+mark]=thetapvas[num_covars+tindex[bitstart+j]+m*num_fixed];
}
else	//not a top, so must test
{
if(mode==131||fastgwa==0)
{
value=YTdata[j+m*bitlength]/XTCX[j];
value2=(Pvarphens[m]*num_samples_use/XTCX[j]-pow(value,2))/(num_samples_use-num_fixed-1);
}
else
{
value=YTdata[j+m*bitlength]/XTCX[j];
value2=pow(XTCX[j],-1);
}

stats[0+mark]=value*Pscales[m]*Plambdas[m];
stats[1+mark]=pow(value2,.5)*Pscales[m]*Proots[m];
stats[2+mark]=stats[0+mark]/stats[1+mark];
stats[3+mark]=erfc(fabs(stats[2+mark])*M_SQRT1_2);
}
}	//end of testing
}	//end of j loop

if(adjpreds==1)	//revisit the most significant predictors
{
for(j=0;j<bitlength;j++)
{
mark=4*j;

if(mults[bitstart+j]!=-9999)	//might consider
{
if(tindex[bitstart+j]==-9999&&stats[3+mark]<threshold)	//recalculate XTCX, then test statistics
{
//save original value
value4=XTCX[j];

//need to subtract Z2 x Z3TX from data (its ok if have already done this for previous phenotypes)
alpha=1.0;beta=0.0;
dgemv_("T", &num_samples_use, &num_fixed, &alpha, Z3+m*num_samples_use*num_fixed, &num_samples_use, data+(size_t)j*num_samples_use, &one, &beta, ZTdata, &one);

alpha=-1.0;beta=1.0;
dgemv_("N", &num_samples_use, &num_fixed, &alpha, Z2+m*num_samples_use*num_fixed, &num_samples_use, ZTdata, &one, &beta, data+(size_t)j*num_samples_use, &one);

if(mode==131)
{
XTCX[j]=0;for(i=0;i<num_samples_use;i++){XTCX[j]+=pow(data[(size_t)j*num_samples_use+i],2);}
}
else	//must compute
{
XTCX[j]=0;for(i=0;i<num_samples_use;i++){XTCX[j]+=pow(data[(size_t)j*num_samples_use+i],2)*nullweights[i+m*num_samples_use];}
}

if(mode==131||fastgwa==0)
{
value=YTdata[j+m*bitlength]/XTCX[j];
value2=(Pvarphens[m]*num_samples_use/XTCX[j]-pow(value,2))/(num_samples_use-num_fixed-1);
}
else
{
value=YTdata[j+m*bitlength]/XTCX[j];
value2=pow(XTCX[j],-1);
}

stats[0+mark]=value*Pscales[m]*Plambdas[m];
stats[1+mark]=pow(value2,.5)*Pscales[m]*Proots[m];
stats[2+mark]=stats[0+mark]/stats[1+mark];
stats[3+mark]=erfc(fabs(stats[2+mark])*M_SQRT1_2);

//restore original value
XTCX[j]=value4;
}}	//end of revisiting
}	//end of j loop
}	//end of adjpreds=1

if(spatest==1)	//compute spa test statistic
{
//work out which predictors to test
count=0;
for(j=0;j<bitlength;j++)
{
mark=4*j;

if(mults[bitstart+j]!=-9999)	//might test
{
if(tindex[bitstart+j]==-9999&&stats[3+mark]<spathresh)	//will compute SPA
{usedpreds[count]=j;count++;}
else	//will report non-SPA results
{spastatus[j]=0;}
}
else	//trivial, so will not test, but helps to set spastatus to missing
{spastatus[j]=-9999;}
}

while(1)
{
#pragma omp parallel for private(j2, j, mark) schedule(static)
for(j2=0;j2<count;j2++)
{
j=usedpreds[j2];
mark=4*j;

if(mode==131||fastgwa==0){spastatus[j]=spa_test(YTdata[j+m*bitlength], data+(size_t)j*num_samples_use, num_samples_use, num_knots, knots+m*num_knots, num_bins, bins, CGF0+m*num_knots, CGF1+m*num_knots, CGF2+m*num_knots, CGF3+m*num_knots, stats+mark, Proots[m]);}
else
{
if(spaside==1){spastatus[j]=spa_logistic_one(YTdata[j+m*bitlength], data+(size_t)j*num_samples_use, num_samples_use, nullprobs2+m*num_samples_use, nullweights2+m*num_samples_use, stats+mark, Proots[m]);}
else{spastatus[j]=spa_logistic_two(YTdata[j+m*bitlength], XTCX[j], data+(size_t)j*num_samples_use, num_samples_use, nullprobs2+m*num_samples_use, nullweights2+m*num_samples_use, stats+mark, Proots[m]);}
}
}

if(mode==132&&fastgwa!=0)	//finished
{break;}

//must have empirical spa solver - see if necessary to increase range for any predictors (and recompute SPA)
count=0;
for(j=0;j<bitlength;j++)
{
if(spastatus[j]==-2){usedpreds[count]=j;count++;}
}
if(count==0){break;}

Pspamax[m]*=5;
if(mpheno!=-1){printf("Increasing SPA range to %f and recomputing empirical CGF\n", Pspamax[m]);}
else{printf("Phenotype %d: increasing SPA range to %f and recomputing empirical CGF\n", m+1, Pspamax[m]);}

empirical_cumulants(Yadj+m*num_samples_use, num_samples_use, num_knots, knots+m*num_knots, num_bins, bins, CGF0+m*num_knots, CGF1+m*num_knots, CGF2+m*num_knots, CGF3+m*num_knots, Pspamax[m]);
}
}

//reopen output files
if(kvikparity==0){sprintf(filename2,"%s.assoc", outfile);}
else{sprintf(filename2,"%s.pheno%d.assoc", outfile, keepresps[m]+1);}
if((output2=fopen(filename2,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename2);exit(1);}
if(kvikparity==0){sprintf(filename3,"%s.summaries", outfile);}
else{sprintf(filename3,"%s.pheno%d.summaries", outfile, keepresps[m]+1);}
if((output3=fopen(filename3,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename3);exit(1);}
if(kvikparity==0){sprintf(filename4,"%s.pvalues", outfile);}
else{sprintf(filename4,"%s.pheno%d.pvalues", outfile, keepresps[m]+1);}
if((output4=fopen(filename4,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename4);exit(1);}

//save results
for(j=0;j<bitlength;j++)
{
mark=4*j;

if(mults[bitstart+j]!=-9999)	//tested predictor
{
//print assoc
fprintf(output2, "%d\t%s\t%.0f\t%s\t%s\t", chr[bitstart+j], preds[bitstart+j], bp[bitstart+j], along1[bitstart+j], along2[bitstart+j]);
fprintf(output2, "%.4f\t%.4e\t%.4e\t%.4e\t", stats[2+mark], stats[3+mark], stats[0+mark], stats[1+mark]);
if(nonsnp==0){fprintf(output2, "%.6f\t", centres[bitstart+j]/2+(centres[bitstart+j]>1)*(1-centres[bitstart+j]));}
else{fprintf(output2, "NA\t");}
if(genprobs<2){fprintf(output2, "%.4f\tNA\t", rates[bitstart+j]);}
else{fprintf(output2, "%.4f\t%.2f\t", rates[bitstart+j], infos[bitstart+j]);}
if(spastatus[j]==0){fprintf(output2, "NOT_USED\n");}
if(spastatus[j]==1){fprintf(output2, "SUCCESS\n");}
if(spastatus[j]==2){fprintf(output2, "APPROX\n");}
if(spastatus[j]==-1){fprintf(output2, "FAILED\n");}

//print summaries and pvalues - effective sample size allows for Pincs
fprintf(output3, "%s\t%s\t%s\t%.4f\t%.1f\t%.4f\n", preds[bitstart+j], along1[bitstart+j], along2[bitstart+j], stats[2+mark], rates[bitstart+j]*respcounts[m]*Pincs[m], centres[bitstart+j]/2);
fprintf(output4, "%s\t%.4e\n", preds[bitstart+j], stats[3+mark]);
}
else	//did not test, but will include in assoc if not doing qc
{
if(passqc==0)
{
fprintf(output2, "%d\t%s\t%.0f\t%s\t%s\t", chr[bitstart+j], preds[bitstart+j], bp[bitstart+j], along1[bitstart+j], along2[bitstart+j]);
fprintf(output2, "NA\tNA\tNA\tNA\t");
if(nonsnp==0){fprintf(output2, "%.6f\t", centres[bitstart+j]/2+(centres[bitstart+j]>1)*(1-centres[bitstart+j]));}
else{fprintf(output2, "NA\t");}
if(genprobs<2){fprintf(output2, "%.4f\tNA\t", rates[bitstart+j]);}
else{fprintf(output2, "%.4f\t%.2f\t", rates[bitstart+j], infos[bitstart+j]);}
if(spastatus[j]==0){fprintf(output2, "NOT_USED\n");}
}
}
}	//end of j loop

fclose(output2);
fclose(output3);
fclose(output4);
}	//end of m loop

bit++;
}	//end of while loop
printf("\n");

fclose(output);

count=0;for(j=0;j<num_preds_use;j++){count+=(mults[j]==-9999);}
if(count==num_preds_use)
{
if(passqc==0)
{printf("Warning, all %d predictors were excluded because they were trivial (showed no variation)\n\n", count);}
else
{printf("Warning, all %d predictors were excluded because they were trivial or failed quality control\n\n", count);}
}
else
{
if(count>0)
{
if(passqc==0)
{printf("Warning, %d predictors were excluded because they were trivial (showed no variation)\n\n", count);}
else
{printf("Warning, %d predictors were excluded because they were trivial or failed quality control\n\n", count);}
}
}

if(kvikparity==0){printf("Main results saved in %s.assoc, with a summary version in %s.summaries and p-values in %s.pvalues\n\n", outfile, outfile, outfile);}
else
{
if(mpheno!=-1){printf("Main results saved in %s.pheno%d.assoc, with a summary version in %s.pheno%d.summaries and p-values in %s.pheno%d.pvalues\n\n", outfile, keepresps[0]+1, outfile, keepresps[0]+1, outfile, keepresps[0]+1);}
else{printf("Main results saved in %s.phenoX.assoc, with a summary version in %s.phenoX.summaries and p-values in %s.phenoX.pvalues, where X is the phenotype number\n\n", outfile, outfile, outfile);}
}

free(Pincs);free(Plambdas);free(Proots);
if(fastgwa==0){free(chrindex);free(prs);}
else{free(pedprs);}
free(data);free(datasqs);
free(order);free(Y);free(Z);
free(Pscales);free(Pvarphens);
free(thetas);free(thetasds);;free(thetapvas);free(Yadj);free(Yadj2);
if(mode==132){free(nullprobs);free(nullweights);}
free(tindex);free(stats);free(spastatus);
free(YTdata);free(XTCX);
if(adjpreds==1){free(Z2);free(Z3);free(ZTdata);}
if(spatest==1)
{
free(usedpreds);
if(mode==131||fastgwa==0)
{
for(j=0;j<num_knots*num_resps_use;j++){free(CGF0[j]);free(CGF1[j]);free(CGF2[j]);free(CGF3[j]);}
free(Pspamax);free(knots);free(bins);free(CGF0);free(CGF1);free(CGF2);free(CGF3);
}
else{free(nullprobs2);free(nullweights2);}
}
if(binary==0){gzclose(datainputgz);}

///////////////////////////

