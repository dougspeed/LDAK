/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Make phenotypes

///////////////////////////

//allocate variables
data_warn2(bitsize,num_samples_use);
data=malloc(sizeof(double)*num_samples_use*bitsize);

unders=malloc(sizeof(double*)*num_phenos);
phens=malloc(sizeof(double*)*num_phenos);
for(m=0;m<num_phenos;m++)
{
unders[m]=malloc(sizeof(double)*num_samples_use);
phens[m]=malloc(sizeof(double)*num_samples_use);
}

//set unders to zero
for(m=0;m<num_phenos;m++)
{
for(i=0;i<num_samples_use;i++){unders[m][i]=0;}
}

if(strcmp(covarfile,"blank")!=0)	//using covariates
{
vsds=malloc(sizeof(double)*num_covars);
veffs=malloc(sizeof(double)*num_covars);
vadds=malloc(sizeof(double)*num_samples_use);
}

if(prev!=-9999)	//will use liabilities
{
liabs=malloc(sizeof(double*)*num_phenos);
for(m=0;m<num_phenos;m++){liabs[m]=malloc(sizeof(double)*num_samples_use);}
}

//prepare for reading data
if(binary==0){open_datagz(&datainputgz, datafile, num_samples, genskip, genheaders, genprobs);}
current=0;

if(her>0)
{
//deal with progress file
sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fclose(output);

//ready for bit loop
bittotal=(data_length-1)/bitsize+1;
ecount=0;
wcount=0;
for(bit=0;bit<bittotal;bit++)
{
bitstart=bit*bitsize;
bitend=(bit+1)*bitsize;
if(bitend>data_length){bitend=data_length;}
bitlength=bitend-bitstart;

printf("Making phenotypes for Chunk %d of %d\n", bit+1, bittotal);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Making phenotypes for Chunk %d of %d\n", bit+1, bittotal);
fclose(output);

current=read_data_fly(datafile, dtype, data, NULL, num_samples_use, keepsamps, bitstart, bitend, keeppreds_use, datainputgz, current, num_samples, num_preds, genskip, genheaders, genprobs, bgen_indexes, missingvalue, -9999, -9999, nonsnp, maxthreads);
stand_data(data, centres+bitstart, mults+bitstart, sqdevs+bitstart, rates+bitstart, infos+bitstart, num_samples_use, bitlength, missingvalue, power, 0, hwestand, weights+bitstart, 1);

for(j=0;j<bitlength;j++)
{
for(m=0;m<num_phenos;m++)
{
if(effects[m][bitstart+j]!=0)
{
if(mults[j]==-9999||weights[j]==0)
{
if(weights[j]==0)
{
if(ecount<5){printf("Warning, Predictor %s was picked to be causal but has weight zero\n", preds[j]);}
ecount++;
}
else
{
if(wcount<5){printf("Warning, Predictor %s was picked to be causal but is trivial\n", preds[j]);}
wcount++;
}
}
else
{
for(i=0;i<num_samples_use;i++){unders[m][i]+=data[(size_t)j*num_samples_use+i]*effects[m][bitstart+j];}
}
}}	//end of effect!=0 and m loop
}	//end of j loop
}	//end of bit loop
if(ecount>=5){printf("In total, %d weight-zero predictors were picked to be causal\n", ecount);}
if(wcount>=5){printf("In total, %d trivial predictors were picked to be causal\n", wcount);}
printf("\n");
}	//end of her>0

//have breeding values, now add noise

for(m=0;m<num_phenos;m++)
{
//work out how much noise to add
sum=0;sumsq=0;
for(i=0;i<num_samples_use;i++){sum+=unders[m][i];sumsq+=pow(unders[m][i],2);}
mean=sum/num_samples_use;
var=sumsq/num_samples_use-pow(mean,2);
if(var==0&&her!=0)
{
printf("Error, phenotype %d has variance zero, indicating that all its causal predictors are trivial", m+1);
if(strcmp(weightsfile,"blank")!=0){printf(" or have weight zero");}
printf("\n\n");exit(1);
}

if(her==0){value=1;}
else{value=pow(var*(1-her)/her,.5);}

//add on noise
for(i=0;i<num_samples_use;i++){phens[m][i]=unders[m][i]+rnorm_safe()*value;}

//will scale everything so variance is one
sum=0;sumsq=0;
for(i=0;i<num_samples_use;i++){sum+=phens[m][i];sumsq+=pow(phens[m][i],2);}
mean=sum/num_samples_use;
var=sumsq/num_samples_use-pow(mean,2);

value=pow(var,-.5);
for(i=0;i<num_samples_use;i++){unders[m][i]*=value;phens[m][i]*=value;}
for(j=0;j<data_length;j++){effects[m][j]*=value;}
}

if(strcmp(covarfile,"blank")!=0)	//add covariate contribution - remember first column of covars is intercept
{
//get sds of covariates
vsds[0]=0;
for(j=1;j<num_covars;j++)
{
sum=0;sumsq=0;
for(i=0;i<num_samples_use;i++){sum+=covar[i+j*num_samples_use];sumsq+=pow(covar[i+j*num_samples_use],2);}
mean=sum/num_samples_use;
var=sumsq/num_samples_use-pow(mean,2);
vsds[j]=pow(var,.5);
}

for(m=0;m<num_phenos;m++)
{
//sample effects
veffs[0]=0;
for(j=1;j<num_covars;j++){veffs[j]=rnorm_safe()*pow(vsds[j],-1);}

//get contributions
alpha=1.0;beta=0.0;
dgemv_("N", &num_samples_use, &num_covars, &alpha, covar, &num_samples_use, veffs, &one, &beta, vadds, &one);

//scale contributions to have mean zero and variance cher/(1-cher)
sum=0;sumsq=0;
for(i=0;i<num_samples_use;i++){sum+=vadds[i];sumsq+=pow(vadds[i],2);}
mean=sum/num_samples_use;
var=sumsq/num_samples_use-pow(mean,2);

value=pow(cher/(1-cher)/var,.5);
for(i=0;i<num_samples_use;i++){vadds[i]=(vadds[i]-mean)*value;}

//add contributions to phens
for(i=0;i<num_samples_use;i++){phens[m][i]+=vadds[i];}

//scale everything so variance is one
sum=0;sumsq=0;
for(i=0;i<num_samples_use;i++){sum+=phens[m][i];sumsq+=pow(phens[m][i],2);}
mean=sum/num_samples_use;
var=sumsq/num_samples_use-pow(mean,2);

value=pow(var,-.5);
for(i=0;i<num_samples_use;i++){unders[m][i]*=value;phens[m][i]*=value;}
for(j=0;j<data_length;j++){effects[m][j]*=value;}
}
}

if(prev!=-9999)	//convert to binary
{
value=normal_inv(1-prev);
for(m=0;m<num_phenos;m++)
{
for(i=0;i<num_samples_use;i++){liabs[m][i]=phens[m][i];phens[m][i]=(liabs[m][i]>value);}
}
}

//now save
sprintf(filename2, "%s.pheno", outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n", filename2);exit(1);}
fprintf(output2, "FID IID");
for(m=0;m<num_phenos;m++){fprintf(output2, " Phen%d", m+1);}
fprintf(output2,"\n");

for(i=0;i<num_samples_use;i++)
{
fprintf(output2, "%s %s", ids1[i], ids2[i]);
for(m=0;m<num_phenos;m++)
{
if(prev!=-9999){fprintf(output2, " %d", (int)phens[m][i]);}
else{fprintf(output2, " %.6f", phens[m][i]);}
}
fprintf(output2, "\n");
}
fclose(output2);

sprintf(filename3, "%s.breed", outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n", filename3);exit(1);}
fprintf(output3, "FID IID");
for(m=0;m<num_phenos;m++){fprintf(output3, " Phen%d", m+1);}
fprintf(output3,"\n");

for(i=0;i<num_samples_use;i++)
{
fprintf(output3, "%s %s", ids1[i], ids2[i]);
for(m=0;m<num_phenos;m++){fprintf(output3, " %.6f", unders[m][i]);}
fprintf(output3, "\n");
}
fclose(output3);

sprintf(filename4, "%s.effects", outfile);
if((output4=fopen(filename4,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n", filename4);exit(1);}
fprintf(output4,"Phenotype Predictor A1 A2 Centre Effect Var_Explained\n");

if(her>0)
{
for(m=0;m<num_phenos;m++)
{
for(j=0;j<data_length;j++)
{
if(effects[m][j]!=0&&mults[j]!=-9999)
fprintf(output4, "%d %s %c %c %.6f %.6f %.6f\n", m+1, preds[j], al1[j], al2[j], centres[j], mults[j]*effects[m][j]*pow(weights[j],.5), pow(mults[j]*effects[m][j],2)*weights[j]*sqdevs[j]);
}
}
}
fclose(output4);

if(prev!=-9999)
{
sprintf(filename5, "%s.liab", outfile);
if((output5=fopen(filename5,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n", filename5);exit(1);}
fprintf(output5, "FID IID");
for(m=0;m<num_phenos;m++){fprintf(output5, " Phen%d", m+1);}
fprintf(output5,"\n");

for(i=0;i<num_samples_use;i++)
{
fprintf(output5, "%s %s", ids1[i], ids2[i]);
for(m=0;m<num_phenos;m++){fprintf(output5, " %.6f", liabs[m][i]);}
fprintf(output5, "\n");
}
fclose(output5);
}

if(prev!=-9999){printf("Phenotypes saved in %s, with liabilities in %s, breeding values in %s and effect sizes in %s\n\n", filename2, filename5, filename3, filename4);}
else{printf("Phenotypes saved in %s, with breeding values in %s and effect sizes in %s\n\n", filename2, filename3, filename4);}

//free allocation from setdl.c
for(m=0;m<num_phenos;m++){free(effects[m]);}free(effects);

//frees from above
free(data);
for(m=0;m<num_phenos;m++){free(unders[m]);free(phens[m]);}free(unders);free(phens);
if(strcmp(covarfile,"blank")!=0){free(vsds);free(veffs);free(vadds);}
if(prev!=-9999){for(m=0;m<num_phenos;m++){free(liabs[m]);}free(liabs);}
if(binary==0){gzclose(datainputgz);}

///////////////////////////

