/*
Copyright 2026 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//LOCO regression - might permute, have spa and sample weights (no top-preds or enviro)
//performing linear regression (mode==131), quasi-logistic regression (mode==132&fastgwa==0) or logistic regression (mode==131&fastgwa!=0)
//for quasi-logistic spa, found best to regress (Y-mu) on X (instead of (Y-mu)/W on WX or alternatives)
//if using sample weights, will have adjpreds=2; if multi=1, must have fastgwa=0
//if mode==132, will have set sandwich=0

///////////////////////////

//check and read loco.details for all phenotypes - want increase in mse, scaling and root of scaling, and effect calibration - and check prs files

Plambdas=malloc(sizeof(double)*num_resps_use);
Proots=malloc(sizeof(double)*num_resps_use);
Pcals=malloc(sizeof(double)*num_resps_use);
Pcals2=malloc(sizeof(double)*num_resps_use);

for(m=0;m<num_resps_use;m++)
{
if(mpheno!=-1){sprintf(filename,"%s.loco.details", prsfile);}
else{sprintf(filename,"%s.pheno%d.loco.details", prsfile, keepresps[m]+1);}
if(just_check(filename)!=0)
{
if(kvikstep==-9999&&gctastep==-9999&&faststep==-9999){printf("Error reading %s; this file should have been created using \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n", filename);}
else{printf("Error reading %s; this file should have been created using \"--kvik-step1\", \"--GCTA-LOCO-step1\" or \"--fastGWA-step1\"\n\n", filename);}
exit(1);
}

count=countrows(filename);
if(count!=7)
{printf("Error, %s should have seven rows (not %d), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n", filename, count);exit(1);}

//open check first row
if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}
if(fscanf(input, "%s %d ", readstring, &readint)!=2)
{printf("Error reading Row 1 of %s\n\n", filename);exit(1);}
if(strcmp(readstring,"Actual_Sample_Size")!=0)
{printf("Error reading %s; Row 1 should begin \"Actual_Sample_Size\" (not %s), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n", filename, readstring);exit(1);}

//read effective sample size, and divide by actual sample size - used to assign this to Pincs, then scale sample sizes by this
if(fscanf(input, "%s %lf ", readstring, &readdouble)!=2)
{printf("Error reading Row 2 of %s\n\n", filename);exit(1);}
if(strcmp(readstring,"Approx_Effective_Sample_Size")!=0)
{printf("Error reading %s; Row 2 should begin \"Approx_Effective_Sample_Size\" (not %s), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n", filename, readstring);exit(1);}
//Pincs[m]=readdouble/num_samples_use;

//read scaling, and also get square root (but set to one if using sandwich)
if(fscanf(input, "%s %lf ", readstring, Plambdas+m)!=2)
{printf("Error reading Row 3 of %s\n\n", filename);exit(1);}
if(strcmp(readstring,"Scaling_Estimate")!=0)
{printf("Error reading %s; Row 3 should begin \"Scaling_Estimate\" (not %s), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n", filename, readstring);exit(1);}
Proots[m]=pow(Plambdas[m],.5);

if(sandwich==1){Plambdas[m]=1;Proots[m]=1;}

//skip one row, and get effect calibration, then set SD calibration to effect calibration / square root
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}
if(fscanf(input, "%s %lf ", readstring, Pcals+m)!=2)
{printf("Error reading Row 5 of %s\n\n", filename);exit(1);}
if(strcmp(readstring,"Effect_Calibration")!=0)
{printf("Error reading %s; Row 5 should begin \"Effect_Calibration\" (not %s), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n", filename, readstring);exit(1);}
Pcals2[m]=Pcals[m]/Proots[m];

fclose(input);
}	//end of m loop

////////

//prepare for reading prs files

value=(double)num_samples_use/1024/1024/1024*num_chr2*num_resps_use;
if(value>1){printf("Warning, to store the PRS requires %.1f Gb; sorry, this can not be reduced\n\n", value);}

chrindex2=malloc(sizeof(int)*num_chr2);
prs=malloc(sizeof(double)*num_samples_use*num_chr2*num_resps_use);

if(fastgwa==0)	//open first prs file and check chromosomes
{
if(mpheno!=-1){sprintf(filename,"%s.loco.prs", prsfile);}
else{sprintf(filename,"%s.pheno1.loco.prs", prsfile);}
if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}
if(fscanf(input, "%s %s ", readstring, readstring2)!=2)
{printf("Error reading first two elements of %s\n\n", filename);exit(1);}

//read chr numbers
for(j=0;j<num_chr2;j++)
{
if(fscanf(input, "Chr%d ", chrindex2+j)!=1)
{printf("Error reading Element %d of %s\n\n", j+3, filename);exit(1);}
}
fclose(input);

//check all chr covered
usedpreds=malloc(sizeof(int)*(chrindex2[num_chr2-1]+1));
for(j=0;j<chrindex2[num_chr2-1]+1;j++){usedpreds[j]=0;}
for(j=0;j<num_chr2;j++){usedpreds[chrindex2[j]]=1;}

for(j=0;j<data_length;j++)
{
if(chr[j]>chrindex2[num_chr2-1])
{printf("Error, will be analyzing predictors on Chromosome %d, but this chromosome was not included when making the PRS\n\n", chr[j]);exit(1);}
if(usedpreds[chr[j]]==0)
{printf("Error, will be analyzing predictors on Chromosome %d, but this chromosome was not included when making the PRS\n\n", chr[j]);exit(1);}
}
free(usedpreds);
}
else	//set chrindex2 artificially (know that num_chr<=num_chr2)
{
chrindex2[0]=chr[0];
count=1;
for(j=1;j<data_length;j++)
{
if(chr[j]>chrindex2[count-1]){chrindex2[count]=chr[j];count++;}
}
while(count<num_chr2){chrindex2[count]=chrindex2[count-1];count++;}
}

//can now read (will have values for all samples)

for(m=0;m<num_resps_use;m++)	//check then read each PRS file in turn
{
if(mpheno!=-1){sprintf(filename,"%s.loco.prs", prsfile);}
else{sprintf(filename,"%s.pheno%d.loco.prs", prsfile, keepresps[m]+1);}

wantids=malloc(sizeof(char*)*num_samples_use2);
indexer=malloc(sizeof(int)*num_samples_use2);
indexer2=malloc(sizeof(int)*num_samples_use2);

printf("Reading LOCO PRS for %d samples from %s\n", num_samples_use2, filename);

read_ids(filename, NULL, NULL, wantids, num_samples_use2, NULL, 1, 0);
count=find_strings(wantids, num_samples_use2, ids3, num_samples_use, indexer, indexer2, NULL, NULL, NULL, NULL, 3);
if(count!=num_samples_use){printf("Error 6XX; please tell Doug %d %d\n\n", count, num_samples_use);exit(1);}

//open PRS and skip header row
if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}

found=0;
for(i=0;i<num_samples_use2;i++)
{
if(i==indexer[found])	//will be using
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
else	//have just one value - so copy into remaining chromosomes
{
if(fscanf(input, "%s %s %lf ", readstring, readstring2, prs+indexer2[found]+m*num_chr2*num_samples_use)!=3)
{printf("Error reading Row %d of %s\n\n", i+2, filename);exit(1);}
for(j=1;j<num_chr2;j++){prs[indexer2[found]+(j+m*num_chr2)*num_samples_use]=prs[indexer2[found]+m*num_chr2*num_samples_use];}
}

found++;
if(found==num_samples_use){break;}
}
else	//skip row
{
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}
}
}

fclose(input);

for(i=0;i<num_samples_use2;i++){free(wantids[i]);}free(wantids);free(indexer);free(indexer2);
}

////////

if(dougvar2==99)	//blank prs and set scalings to one
{
printf("blanking prs\n");
for(m=0;m<num_resps_use;m++)
{
for(j=0;j<num_chr2;j++)
{
for(i=0;i<num_samples_use;i++){prs[i+(j+m*num_chr2)*num_samples_use]=0;}
}
Plambdas[m]=1;Proots[m]=1;Pcals[m]=1;Pcals2[m]=1;
}

printf("Revised scaling is one - this does not work for binary\n\n");	
}

if(dougvar2==98)	//set scalings to one
{
printf("changing calibrations\n");
for(m=0;m<num_resps_use;m++){Plambdas[m]=1;Proots[m]=1;Pcals[m]=1;Pcals2[m]=1;}
}

////////

if(multi==1)	//read transformation and inverse
{
Umat=malloc(sizeof(double)*num_resps_use*num_resps_use);
Umat2=malloc(sizeof(double)*num_resps_use*num_resps_use);

sprintf(filename,"%s.multivariate.transformation", prsfile);
for(m=0;m<num_resps_use;m++){read_values(filename, Umat+m*num_resps_use, num_resps_use, NULL, m+1, 0, 0);}

sprintf(filename,"%s.multivariate.inverse", prsfile);
for(m=0;m<num_resps_use;m++){read_values(filename, Umat2+m*num_resps_use, num_resps_use, NULL, m+1, 0, 0);}
}

//passqc indicates whether doing any qc
passqc=(minmaf!=-9999||maxmaf!=-9999||minvar!=-9999||minobs!=-9999||mininfo!=-9999);

//threshold is the p-value threshold used with adjpreds
threshold=0.05;
if(spatest==1&&spathresh>threshold){threshold=spathresh;}

//allocate variables

data_warn2(bitsize,num_samples_use);
data=malloc(sizeof(double)*num_samples_use*bitsize);
datasqs=malloc(sizeof(double)*bitsize);

order=malloc(sizeof(int)*num_samples_use);
sweights=malloc(sizeof(double)*num_samples_use);
Y=malloc(sizeof(double)*num_samples_use*num_resps_use);
Z=malloc(sizeof(double)*num_samples_use*num_fixed);

Pscales=malloc(sizeof(double)*num_resps_use);
Pscales2=malloc(sizeof(double)*num_resps_use);
Pscales3=malloc(sizeof(double)*num_resps_use);
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

stats=malloc(sizeof(double)*4*bitsize*num_resps_use);
spastatus=malloc(sizeof(int)*bitsize*num_resps_use);

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

if(multi==1){stats2=malloc(sizeof(double)*4*bitsize*num_resps_use);}

//prepare for reading data
if(binary==0){open_datagz(&datainputgz, datafile, num_samples, genskip, genheaders, genprobs);}
current=0;

//deal with order
for(i=0;i<num_samples_use;i++){order[i]=i;}
if(permute==1){permute_int(order,num_samples_use);}

if(strcmp(sampwfile,"blank")!=0)	//get square roots of sample weights (W^1/2)
{
readdoubles=malloc(sizeof(double)*num_samples_use);
read_sampwfile(sampwfile, readdoubles, num_samples_use, ids3, ids1, ids2);

//scale roots so the mean inverse is one (might help stability)
sum=0;for(i=0;i<num_samples_use;i++){sum+=pow(readdoubles[i],-1);}
mean=sum/num_samples_use;
for(i=0;i<num_samples_use;i++){readdoubles[i]*=mean;}

for(i=0;i<num_samples_use;i++){sweights[i]=pow(readdoubles[order[i]],0.5);}
free(readdoubles);
}
else
{
for(i=0;i<num_samples_use;i++){sweights[i]=1;}
}

//fill Y (maybe permuted and multiplying by W^1/2)
for(m=0;m<num_resps_use;m++)
{
for(i=0;i<num_samples_use;i++){Y[i+m*num_samples_use]=resp[order[i]+m*num_samples_use]*sweights[i];}
}

//fill Z (maybe permuted and multiplying by W^1/2)
for(j=0;j<num_fixed;j++)
{
for(i=0;i<num_samples_use;i++){Z[i+j*num_samples_use]=covar[order[i]+j*num_samples_use]*sweights[i];}
}

//decide how much to scale effect sizes, sds and spa test statistics (separate to lambdas and roots)
for(m=0;m<num_resps_use;m++)
{
if(mode==131||fastgwa==0)	//using linear (or quasi linear) model
{
//will scale by 1/p, 1/p and 1, where p is prop non-missing
Pscales[m]=pow((double)respcounts[m]/num_samples_use,-1);
Pscales2[m]=pow((double)respcounts[m]/num_samples_use,-1);
Pscales3[m]=1;
}
else	//using logistic model
{
//will scale by 1/p, root(1/p) and ratio, where p is prop non-missing
Pscales[m]=pow((double)respcounts[m]/num_samples_use,-1);
Pscales2[m]=pow((double)respcounts[m]/num_samples_use,-.5);
Pscales3[m]=Pscales[m]/Pscales2[m];
}
}

for(m=0;m<num_resps_use;m++)	//solve null models (and pad any missing values)
{
if(mode==131)	//get thetas, adjusted response and pad missing
{
reg_covar_lin_missing(Y+m*num_samples_use, Z, num_samples_use, num_fixed, thetas+m*num_fixed, thetasds+m*num_fixed, thetapvas+m*num_fixed, Yadj+m*num_samples_use, missingvalue);
}
else		//get thetas, nullprobs, nullweights and adjusted response (Y-mu), padding missing
{
reg_covar_log_missing(Y+m*num_samples_use, Z, num_samples_use, num_fixed, offsets, thetas+m*num_fixed, thetasds+m*num_fixed, thetapvas+m*num_fixed, nullprobs+m*num_samples_use, nullweights+m*num_samples_use, Yadj+m*num_samples_use, 0.001, 100, missingvalue, 0);
}

//save coefficients
if(mpheno!=-1){sprintf(filename,"%s.coeff", outfile);}
else{sprintf(filename,"%s.pheno%d.coeff", outfile, keepresps[m]+1);}
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to\t%s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
if(mode==131){fprintf(output, "Component\tEffect\tSE\tP\n");}
else{fprintf(output, "Component\tLogOR\tSE\tP\n");}
fprintf(output, "Intercept\t%.4e\t%.4e\t%.4e\n", thetas[0+m*num_fixed], thetasds[0+m*num_fixed], thetapvas[0+m*num_fixed]);
for(j=1;j<num_covars;j++){fprintf(output, "Covariate_%d\t%.4e\t%.4e\t%.4e\n",j, thetas[j+m*num_fixed], thetasds[j+m*num_fixed], thetapvas[j+m*num_fixed]);}
fclose(output);
}	//end of m loop

if(multi==1)	//replace Yadj with transformed phenotypes (and regress out covariates to be safe)
{
for(m=0;m<num_resps_use;m++)
{
for(i=0;i<num_samples_use;i++){Yadj2[i+m*num_samples_use]=Yadj[i+m*num_samples_use];}
}
alpha=1.0;beta=0.0;
dgemm_("N", "T", &num_samples_use, &num_resps_use, &num_resps_use, &alpha, Yadj2, &num_samples_use, Umat, &num_resps_use, &beta, Yadj, &num_samples_use);

reg_covar_matrix(Yadj, Z, num_samples_use, num_resps_use, num_fixed);

//compute Pvarphens (will be updated), and get number of non trivial phenotypes
total=0;
for(m=0;m<num_resps_use;m++)
{
sum=0;sumsq=0;for(i=0;i<num_samples_use;i++){sum+=Yadj[i+m*num_samples_use];sumsq+=pow(Yadj[i+m*num_samples_use],2);}
mean=sum/num_samples_use;Pvarphens[m]=sumsq/num_samples_use-pow(mean,2);
total+=(Pvarphens[m]>0);
}
}

if(spatest==0)	//convenient to set spastatus to zero
{
for(m=0;m<num_resps_use;m++)
{
for(j=0;j<bitsize;j++){spastatus[j+m*bitsize]=0;}
}
}
else
{
if(mode==131||fastgwa==0)	//using empirical spa - can set bins and spamax
{
//set bins - evenly spaced from -2 to 2
for(k=0;k<num_bins;k++){bins[k]=-2+(double)k/(num_bins-1)*4;}

for(m=0;m<num_resps_use;m++)
{
if(spamax==-9999)	//set spamax based on variance of Yadj
{
sum=0;sumsq=0;for(i=0;i<num_samples_use;i++){sum+=Yadj[i+m*num_samples_use];sumsq+=pow(Yadj[i+m*num_samples_use],2);}
mean=sum/num_samples_use;var=sumsq/num_samples_use-pow(mean,2);
if(var>0){Pspamax[m]=2000*pow(num_samples_use*var,-.5);}
else{Pspamax[m]=0;}
}
else{Pspamax[m]=spamax;}
}
}
}

if(adjpreds==1)	//fill Z3 with either copies of Z or WZ for each phenotype, then get Z2 = Z inv(ZTZ3) - per-phenotype values unnecessary for linear regression
{
for(m=0;m<num_resps_use;m++)
{
if(mode==131){copy_matrix(num_samples_use, num_fixed, Z, Z3+m*num_samples_use*num_fixed, 0, NULL);}
else{copy_matrix(num_samples_use, num_fixed, Z, Z3+m*num_samples_use*num_fixed, 1, nullweights+m*num_samples_use);}

get_Z2(Z2+m*num_samples_use*num_fixed, Z, Z3+m*num_samples_use*num_fixed, num_samples_use, num_fixed);
}
}

//adjust responses for first chromosome prs and get variances
mark2=0;while(chr[0]!=chrindex2[mark2]){mark2++;}

for(m=0;m<num_resps_use;m++)	//note that for logistic regression, sweights must be one
{
if(mode==131)	//get adjusted phenotype Yadj2 = Yadj - PRS W^1/2, and its variance
{
for(i=0;i<num_samples_use;i++){Yadj2[i+m*num_samples_use]=Yadj[i+m*num_samples_use]-prs[order[i]+(mark2+m*num_chr2)*num_samples_use]*sweights[i];}
sum=0;sumsq=0;for(i=0;i<num_samples_use;i++){sum+=Yadj2[i+m*num_samples_use];sumsq+=pow(Yadj2[i+m*num_samples_use],2);}
mean=sum/num_samples_use;Pvarphens[m]=sumsq/num_samples_use-pow(mean,2);
}
else	//get adjusted phenotype Yadj2 = Yadj - W PRS, and the variance of 1/root(W) Yadj2 (variance redundant if fastgwa!=0)
{
for(i=0;i<num_samples_use;i++){Yadj2[i+m*num_samples_use]=Yadj[i+m*num_samples_use]-prs[order[i]+(mark2+m*num_chr2)*num_samples_use]*nullweights[i+m*num_samples_use];}
sum=0;sumsq=0;for(i=0;i<num_samples_use;i++){sum+=Yadj2[i+m*num_samples_use]*pow(nullweights[i+m*num_samples_use],-.5);sumsq+=pow(Yadj2[i+m*num_samples_use],2)/nullweights[i+m*num_samples_use];}
mean=sum/num_samples_use;Pvarphens[m]=sumsq/num_samples_use-pow(mean,2);
}
}

if(spatest==1)
{
for(m=0;m<num_resps_use;m++)
{
if(mode==131||fastgwa==0)	//using empirical SPA
{
if(mpheno!=-1){printf("Computing empirical CGF (%d knots between %f and %f, and %d bins between -2 and 2)\n", num_knots, -Pspamax[m], Pspamax[m], num_bins);}
else{printf("Phenotype %d: computing empirical CGF (%d knots between %f and %f, and %d bins between -2 and 2)\n", m+1, num_knots, -Pspamax[m], Pspamax[m], num_bins);}

empirical_cumulants(Yadj2+m*num_samples_use, num_samples_use, num_knots, knots+m*num_knots, num_bins, bins, CGF0+m*num_knots, CGF1+m*num_knots, CGF2+m*num_knots, CGF3+m*num_knots, Pspamax[m]);
}
else	//using logistic spa - convenient to get some extra probs and weights
{
for(i=0;i<num_samples_use;i++)
{
nullprobs2[i+m*num_samples_use]=nullprobs[i+m*num_samples_use]+prs[order[i]+(mark2+m*num_chr2)*num_samples_use]*nullweights[i+m*num_samples_use];
if(nullprobs2[i+m*num_samples_use]<=0){nullprobs2[i+m*num_samples_use]=1e-10;}
if(nullprobs2[i+m*num_samples_use]>=1){nullprobs2[i+m*num_samples_use]=1-1e-10;}
nullweights2[i+m*num_samples_use]=nullprobs2[i+m*num_samples_use]*(1-nullprobs2[i+m*num_samples_use]);
}
}
}
if(mode==131||fastgwa==0){printf("\n");}
}

////////

//deal with progress and on-the-fly files

sprintf(filename,"%s.progress", outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}

for(m=0;m<num_resps_use;m++)
{
if(mpheno!=-1){sprintf(filename2,"%s.assoc", outfile);}
else{sprintf(filename2,"%s.pheno%d.assoc", outfile, keepresps[m]+1);}
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
if(mode==131){fprintf(output2, "Chromosome\tPredictor\tBasepair\tA1\tA2\tWald_Stat\tWald_P\tEffect\tSE\tMAF\tCallRate\tMachR2\tSPA_Status\n");}
else{fprintf(output2, "Chromosome\tPredictor\tBasepair\tA1\tA2\tWald_Stat\tWald_P\tApprox_Log_OR\tApprox_SE\tMAF\tCallRate\tMachR2\tSPA_Status\n");}
fclose(output2);

if(mpheno!=-1){sprintf(filename3,"%s.summaries", outfile);}
else{sprintf(filename3,"%s.pheno%d.summaries", outfile, keepresps[m]+1);}
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3, "Predictor\tA1\tA2\tZ\tn\tA1Freq\n");
fclose(output3);

if(mpheno!=-1){sprintf(filename4,"%s.pvalues", outfile);}
else{sprintf(filename4,"%s.pheno%d.pvalues", outfile, keepresps[m]+1);}
if((output4=fopen(filename4,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename4);exit(1);}
fprintf(output4, "Predictor\tP\n");
fclose(output4);
}

if(multi==1)
{
sprintf(filename5,"%s.factors", outfile);
if((output5=fopen(filename5,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename5);exit(1);}
fprintf(output5, "Predictor");
for(m=0;m<num_resps_use;m++){fprintf(output5,"\tStats_%d", m+1);}
fprintf(output5,"\n");
fclose(output5);

sprintf(filename6,"%s.combined", outfile);
if((output6=fopen(filename6,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename6);exit(1);}
fprintf(output6, "Predictor\tChiSq_%d_Stat\tP\tChiSq_1_Stat\n", total);
fclose(output6);
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

if(chr[bitstart]!=chrindex2[mark2])	//new chromosome
{
while(chr[bitstart]!=chrindex2[mark2]){mark2++;}

for(m=0;m<num_resps_use;m++)
{
if(mode==131)	//get adjusted phenotype Yadj2 = Yadj - PRS W^1/2, and its variance
{
for(i=0;i<num_samples_use;i++){Yadj2[i+m*num_samples_use]=Yadj[i+m*num_samples_use]-prs[order[i]+(mark2+m*num_chr2)*num_samples_use]*sweights[i];}
sum=0;sumsq=0;for(i=0;i<num_samples_use;i++){sum+=Yadj2[i+m*num_samples_use];sumsq+=pow(Yadj2[i+m*num_samples_use],2);}
mean=sum/num_samples_use;Pvarphens[m]=sumsq/num_samples_use-pow(mean,2);
}
else	//get adjusted phenotype Yadj2 = Yadj - W PRS, and the variance of 1/root(W) Yadj2 (variance redundant if fastgwa!=0)
{
for(i=0;i<num_samples_use;i++){Yadj2[i+m*num_samples_use]=Yadj[i+m*num_samples_use]-prs[order[i]+(mark2+m*num_chr2)*num_samples_use]*nullweights[i+m*num_samples_use];}
sum=0;sumsq=0;for(i=0;i<num_samples_use;i++){sum+=Yadj2[i+m*num_samples_use]*pow(nullweights[i+m*num_samples_use],-.5);sumsq+=pow(Yadj2[i+m*num_samples_use],2)/nullweights[i+m*num_samples_use];}
mean=sum/num_samples_use;Pvarphens[m]=sumsq/num_samples_use-pow(mean,2);
}
}

if(spatest==1)
{
if(mode==131||fastgwa==0)	//using empirical SPA
{
if(mpheno!=-1){printf("Recomputing empirical CGF\n");}
else{printf("Recomputing empirical CGF for each phenotype\n");}
for(m=0;m<num_resps_use;m++)
{
empirical_cumulants(Yadj2+m*num_samples_use, num_samples_use, num_knots, knots+m*num_knots, num_bins, bins, CGF0+m*num_knots, CGF1+m*num_knots, CGF2+m*num_knots, CGF3+m*num_knots, Pspamax[m]);
}
}
else	//using logistic spa - convenient to get some extra probs and weights
{
for(m=0;m<num_resps_use;m++)
{
for(i=0;i<num_samples_use;i++)
{
nullprobs2[i+m*num_samples_use]=nullprobs[i+m*num_samples_use]+prs[order[i]+(mark2+m*num_chr2)*num_samples_use]*nullweights[i+m*num_samples_use];
if(nullprobs2[i+m*num_samples_use]<=0){nullprobs2[i+m*num_samples_use]=1e-10;}
if(nullprobs2[i+m*num_samples_use]>=1){nullprobs2[i+m*num_samples_use]=1-1e-10;}
nullweights2[i+m*num_samples_use]=nullprobs2[i+m*num_samples_use]*(1-nullprobs2[i+m*num_samples_use]);
}
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

if(strcmp(sampwfile,"blank")!=0)	//multiply by W^1/2
{copy_matrix(num_samples_use, bitlength, data, data, 1, sweights);}

//get XTY
alpha=1.0;beta=0.0;
dgemm_("T", "N", &bitlength, &num_resps_use, &num_samples_use, &alpha, data, &num_samples_use, Yadj2, &num_samples_use, &beta, YTdata, &bitlength);

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
if(adjpreds!=2)	//already have
{
for(j=0;j<bitlength;j++){XTCX[j]=num_samples_use*sqdevs[bitstart+j];}
}
else	//must compute
{
#pragma omp parallel for private(j, i) schedule(static)
for(j=0;j<bitlength;j++)
{
XTCX[j]=0;for(i=0;i<num_samples_use;i++){XTCX[j]+=pow(data[(size_t)j*num_samples_use+i],2);}
}
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

if(Pvarphens[m]>0)
{
#pragma omp parallel for private(j, mark, value, value2, value3, i) schedule(static)
for(j=0;j<bitlength;j++)
{
mark=4*j+m*4*bitlength;

if(mults[bitstart+j]!=-9999)	//will be testing - remember Pscales and Pscales2 and Pcals and Pcals2
{
if(sandwich==0)	//use standard estimates (can be linear or logistic)
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
}
else	//use sandwich estimates (must be linear)
{
//effect size remains YTCX/XTCX
value=YTdata[j+m*bitlength]/XTCX[j];

//variance is sum(Xi^2 ri^2)/XTCX^2, where ri are residuals
value3=0;
for(i=0;i<num_samples_use;i++)
{value3+=pow(data[(size_t)j*num_samples_use+i]*(Yadj2[i+m*num_samples_use]-value*data[(size_t)j*num_samples_use+i]),2);}
value2=value3*pow(XTCX[j],-2);
}

if(multi==0)	//will correct for padding now
{
stats[0+mark]=value*Pscales[m]*Pcals[m];
stats[1+mark]=pow(value2,.5)*Pscales2[m]*Pcals2[m];
}
else	//will correct for padding after transforming
{
stats[0+mark]=value*Pcals[m];
stats[1+mark]=pow(value2,.5)*Pcals2[m];
}
stats[2+mark]=stats[0+mark]/stats[1+mark];
stats[3+mark]=erfc(fabs(stats[2+mark])*M_SQRT1_2);
}	//end of testing
}	//end of j loop

if(adjpreds==1)	//revisit the most significant predictors
{
for(j=0;j<bitlength;j++)
{
mark=4*j+m*4*bitlength;

if(mults[bitstart+j]!=-9999)	//might consider
{
if(stats[3+mark]<threshold)	//recalculate XTCX, then test statistics
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

if(sandwich==0)	//use standard estimates (can be linear or logistic)
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
}
else	//use sandwich estimates (must be linear)
{
//effect size remains YTCX/XTCX
value=YTdata[j+m*bitlength]/XTCX[j];

//variance is sum(Xi^2 ri^2)/XTCX^2, where ri are residuals
value3=0;
for(i=0;i<num_samples_use;i++)
{value3+=pow(data[(size_t)j*num_samples_use+i]*(Yadj2[i+m*num_samples_use]-value*data[(size_t)j*num_samples_use+i]),2);}
value2=value3*pow(XTCX[j],-2);
}

if(multi==0)	//will correct for padding now
{
stats[0+mark]=value*Pscales[m]*Pcals[m];
stats[1+mark]=pow(value2,.5)*Pscales2[m]*Pcals2[m];
}
else	//will correct for padding after transforming
{
stats[0+mark]=value*Pcals[m];
stats[1+mark]=pow(value2,.5)*Pcals2[m];
}
stats[2+mark]=stats[0+mark]/stats[1+mark];
stats[3+mark]=erfc(fabs(stats[2+mark])*M_SQRT1_2);

//restore original value
XTCX[j]=value4;
}}	//end of revisiting
}	//end of j loop
}	//end of adjpreds=1

if(spatest==1)	//compute spa test statistic - remember Pscales3 and Proots (Pscales3=1 for multi, so no need to exclude)
{
//work out which predictors to test
count=0;
for(j=0;j<bitlength;j++)
{
mark=4*j+m*4*bitlength;

if(mults[bitstart+j]!=-9999)	//might test
{
if(stats[3+mark]<spathresh)	//will compute SPA
{usedpreds[count]=j;count++;}
else	//will report non-SPA results
{spastatus[j+m*bitlength]=0;}
}
else	//trivial, so will not test, but helps to set spastatus to missing
{spastatus[j+m*bitlength]=-9999;}
}

while(1)
{
#pragma omp parallel for private(j2, j, mark) schedule(static)
for(j2=0;j2<count;j2++)
{
j=usedpreds[j2];
mark=4*j+m*4*bitlength;

if(mode==131||fastgwa==0){spastatus[j+m*bitlength]=spa_test(YTdata[j+m*bitlength], data+(size_t)j*num_samples_use, num_samples_use, num_knots, knots+m*num_knots, num_bins, bins, CGF0+m*num_knots, CGF1+m*num_knots, CGF2+m*num_knots, CGF3+m*num_knots, stats+mark, Pscales3[m]*Proots[m]);}
else
{
if(spaside==1){spastatus[j+m*bitlength]=spa_logistic_one(YTdata[j+m*bitlength], data+(size_t)j*num_samples_use, num_samples_use, nullprobs2+m*num_samples_use, nullweights2+m*num_samples_use, stats+mark, Pscales3[m]*Proots[m]);}
else{spastatus[j+m*bitlength]=spa_logistic_two(YTdata[j+m*bitlength], XTCX[j], data+(size_t)j*num_samples_use, num_samples_use, nullprobs2+m*num_samples_use, nullweights2+m*num_samples_use, stats+mark, Pscales3[m]*Proots[m]);}
}
}

if(mode==132&&fastgwa!=0)	//finished
{break;}

//must have empirical spa solver - see if necessary to increase range for any predictors (and recompute SPA)
count=0;
for(j=0;j<bitlength;j++)
{
if(spastatus[j+m*bitlength]==-2){usedpreds[count]=j;count++;}
}
if(count==0){break;}

Pspamax[m]*=5;
if(mpheno!=-1){printf("Increasing SPA range to %f and recomputing empirical CGF\n", Pspamax[m]);}
else{printf("Phenotype %d: increasing SPA range to %f and recomputing empirical CGF\n", m+1, Pspamax[m]);}

empirical_cumulants(Yadj2+m*num_samples_use, num_samples_use, num_knots, knots+m*num_knots, num_bins, bins, CGF0+m*num_knots, CGF1+m*num_knots, CGF2+m*num_knots, CGF3+m*num_knots, Pspamax[m]);
}
}
}
else	//trivial (transformed) phenotype
{
for(j=0;j<bitlength;j++)
{
mark=4*j+m*4*bitlength;
stats[0+mark]=0;stats[1+mark]=0;stats[2+mark]=0;stats[3+mark]=1;
}
}
}	//end of m loop

if(multi==1)	//get combined pvalues then transform statistics
{
//save current values
sprintf(filename5,"%s.factors", outfile);
if((output5=fopen(filename5,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename5);exit(1);}

for(j=0;j<bitlength;j++)
{
if(mults[bitstart+j]!=-9999)	//tested predictor
{
fprintf(output5,"%s", preds[bitstart+j]);
for(m=0;m<num_resps_use;m++)
{
mark=4*j+m*4*bitlength;
fprintf(output5, "\t%.6f", stats[2+mark]);
}
fprintf(output5,"\n");
}
else
{
fprintf(output5,"%s", preds[bitstart+j]);
for(m=0;m<num_resps_use;m++){fprintf(output5, "\tNA");}
fprintf(output5,"\n");
}
}
fclose(output5);

sprintf(filename6,"%s.combined", outfile);
if((output6=fopen(filename6,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename6);exit(1);}

for(j=0;j<bitlength;j++)
{
if(mults[bitstart+j]!=-9999)	//tested predictor
{
sum=0;
for(m=0;m<num_resps_use;m++)
{
mark=4*j+m*4*bitlength;
sum+=pow(stats[2+mark],2);
}
value=pochisq(sum, total);
if(value>0)	//have a valid p-value, so convert to chisq(1) statistic
{value2=pow(normal_inv(value/2),2);}
else	//not valid, so use approximation from wiki
//https://en.wikipedia.org/wiki/Noncentral_chi-squared_distribution#Cumulative_distribution_function
{
value3=pow(sum/total,1.0/3)-1-pow(4.5*total,-1);
value2=pow(value3,2)*4.5*total;
}
fprintf(output6, "%s\t%.4f\t%.4e\t%.4f\n", preds[bitstart+j], sum, value, value2);
}
else{fprintf(output6, "%s\tNA\tNA\tNA\n", preds[bitstart+j]);}
}
fclose(output6);

#pragma omp parallel for private(j, m, m2, mark, mark3) schedule(static)
for(j=0;j<bitlength;j++)
{
if(mults[bitstart+j]!=-9999)	//tested predictor
{
//first compute estimate and variance (save in stats2)
for(m=0;m<num_resps_use;m++)
{
mark=4*j+m*4*bitlength;
stats2[0+mark]=0;
stats2[1+mark]=0;
for(m2=0;m2<num_resps_use;m2++)
{
mark3=4*j+m2*4*bitlength;
stats2[0+mark]+=stats[0+mark3]*Umat2[m+m2*num_resps_use];
stats2[1+mark]+=pow(stats[1+mark3]*Umat2[m+m2*num_resps_use],2);
}
}

//now copy back into stats, remembering to scale
for(m=0;m<num_resps_use;m++)
{
mark=4*j+m*4*bitlength;
stats[0+mark]=stats2[0+mark]*Pscales[m];
stats[1+mark]=pow(stats2[1+mark],.5)*Pscales2[m];
stats[2+mark]=stats[0+mark]/stats[1+mark];
stats[3+mark]=erfc(fabs(stats[2+mark])*M_SQRT1_2);
}
}}
}

for(m=0;m<num_resps_use;m++)	//save
{
//reopen output files
if(mpheno!=-1){sprintf(filename2,"%s.assoc", outfile);}
else{sprintf(filename2,"%s.pheno%d.assoc", outfile, keepresps[m]+1);}
if((output2=fopen(filename2,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename2);exit(1);}
if(mpheno!=-1){sprintf(filename3,"%s.summaries", outfile);}
else{sprintf(filename3,"%s.pheno%d.summaries", outfile, keepresps[m]+1);}
if((output3=fopen(filename3,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename3);exit(1);}
if(mpheno!=-1){sprintf(filename4,"%s.pvalues", outfile);}
else{sprintf(filename4,"%s.pheno%d.pvalues", outfile, keepresps[m]+1);}
if((output4=fopen(filename4,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename4);exit(1);}

//save results
for(j=0;j<bitlength;j++)
{
mark=4*j+m*4*bitlength;

if(mults[bitstart+j]!=-9999)	//tested predictor
{
//print assoc
fprintf(output2, "%d\t%s\t%.0f\t%s\t%s\t", chr[bitstart+j], preds[bitstart+j], bp[bitstart+j], along1[bitstart+j], along2[bitstart+j]);
fprintf(output2, "%.4f\t%.4e\t%.4e\t%.4e\t", stats[2+mark], stats[3+mark], stats[0+mark], stats[1+mark]);
if(nonsnp==0){fprintf(output2, "%.6f\t", centres[bitstart+j]/2+(centres[bitstart+j]>1)*(1-centres[bitstart+j]));}
else{fprintf(output2, "NA\t");}
if(genprobs<2){fprintf(output2, "%.4f\tNA\t", rates[bitstart+j]);}
else{fprintf(output2, "%.4f\t%.2f\t", rates[bitstart+j], infos[bitstart+j]);}
if(spastatus[j+m*bitlength]==0){fprintf(output2, "NOT_USED\n");}
if(spastatus[j+m*bitlength]==1){fprintf(output2, "SUCCESS\n");}
if(spastatus[j+m*bitlength]==2){fprintf(output2, "APPROX\n");}
if(spastatus[j+m*bitlength]==-1){fprintf(output2, "FAILED\n");}

//print summaries and pvalues
fprintf(output3, "%s\t%s\t%s\t%.4f\t%.0f\t%.6f\n", preds[bitstart+j], along1[bitstart+j], along2[bitstart+j], stats[2+mark], rates[bitstart+j]*respcounts[m], centres[bitstart+j]/2);
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
fprintf(output2, "NOT_USED\n");
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

count=0;for(j=0;j<data_length;j++){count+=(mults[j]==-9999);}
if(count==data_length)
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

if(mpheno!=-1){printf("Main results saved in %s.assoc, with a summary version in %s.summaries and p-values in %s.pvalues\n\n", outfile, outfile, outfile);}
else
{printf("Main results saved in %s.phenoX.assoc, with a summary version in %s.phenoX.summaries and p-values in %s.phenoX.pvalues, where X is the phenotype number\n\n", outfile, outfile, outfile);}

free(Plambdas);free(Proots);free(Pcals);free(Pcals2);
free(chrindex2);free(prs);
if(multi==1){free(Umat);free(Umat2);}
free(data);free(datasqs);
free(order);free(sweights);free(Y);free(Z);
free(Pscales);free(Pscales2);free(Pscales3);free(Pvarphens);
free(thetas);free(thetasds);free(thetapvas);free(Yadj);free(Yadj2);
if(mode==132){free(nullprobs);free(nullweights);}
free(stats);free(spastatus);
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
if(multi==1){free(stats2);}
if(binary==0){gzclose(datainputgz);}

///////////////////////////

