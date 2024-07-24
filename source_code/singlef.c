/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Single-predictor linear regression with trios (no top-preds, enviro, spa, sample weights, bootstrapping)

///////////////////////////

//allocate variables

data_warn2(bitsize,4*num_samples_use);

data=malloc(sizeof(double)*3*num_samples_use*bitsize);
data2=malloc(sizeof(double)*num_samples_use*bitsize);

order=malloc(sizeof(int)*num_samples_use);
Y=malloc(sizeof(double)*num_samples_use);
Z=malloc(sizeof(double)*num_samples_use*num_fixed);

thetas=malloc(sizeof(double)*num_fixed);
thetasds=malloc(sizeof(double)*num_fixed);
thetapvas=malloc(sizeof(double)*num_fixed);
Yadj=malloc(sizeof(double)*num_samples_use);

stats=malloc(sizeof(double)*4*bitsize);
stats2=malloc(sizeof(double)*12*bitsize);
stats3=malloc(sizeof(double)*8*bitsize);
stats4=malloc(sizeof(double)*8*bitsize);

YTdata=malloc(sizeof(double)*3*bitsize);
XTCX=malloc(sizeof(double)*9*bitsize);
XTCX2=malloc(sizeof(double)*9*bitsize);
XTCX3=malloc(sizeof(double)*4*bitsize);
XTCX4=malloc(sizeof(double)*4*bitsize);

//prepare for reading data
if(binary==0){open_datagz(&datainputgz, datafile, num_samples, genskip, genheaders, genprobs);}
current=0;

//deal with order (not bootstrapping)
for(i=0;i<num_samples_use;i++){order[i]=i;}
if(permute==1){permute_int(order,num_samples_use);}

//fill Y and Z
for(i=0;i<num_samples_use;i++)
{
Y[i]=resp[order[i]];
for(j=0;j<num_fixed;j++){Z[i+j*num_samples_use]=covar[order[i]+j*num_samples_use];}
}

//solve null model, get thetas and adjust response
reg_covar_lin(Y, Z, num_samples_use, num_covars, num_tops, thetas, thetasds, thetapvas, Yadj, 0, NULL, NULL);

//save
sprintf(filename,"%s.coeff", outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fprintf(output, "Component Effect SD P\n");
fprintf(output, "Intercept %.6f %.6f %.4e\n", thetas[0], thetasds[0], thetapvas[0]);
for(j=1;j<num_covars;j++){fprintf(output, "Covariate_%d %.6f %.6f %.4e\n",j, thetas[j], thetasds[j], thetapvas[j]);}
fclose(output);

////////

//deal with progress and on-the-fly files

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}

sprintf(filename2,"%s.basic",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2, "Chromosome Predictor Basepair A1 A2 Wald_Stat Wald_P Effect SD A1_Mean MAF\n");

sprintf(filename3,"%s.trios",outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3, "Chromosome Predictor Basepair A1 A2 Offspring_Effect Offspring_SD Offspring_Z Offspring_P Father_Effect Father_SD Father_Z Father_P Mother_Effect Mother_SD Mother_Z Mother_P\n");

sprintf(filename4,"%s.offspring.fathers",outfile);
if((output4=fopen(filename4,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename4);exit(1);}
fprintf(output4, "Chromosome Predictor Basepair A1 A2 Offspring_Effect Offspring_SD Offspring_Z Offspring_P Father_Effect Father_SD Father_Z Father_P\n");

sprintf(filename5,"%s.offspring.mothers",outfile);
if((output5=fopen(filename5,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename5);exit(1);}
fprintf(output5, "Chromosome Predictor Basepair A1 A2 Offspring_Effect Offspring_SD Offspring_Z Offspring_P Mother_Effect Mother_SD Mother_Z Mother_P\n");

//ready for bit loop
bittotal=(data_length-1)/bitsize+1;
for(bit=0;bit<bittotal;bit++)
{
bitstart=bit*bitsize;
bitend=(bit+1)*bitsize;
if(bitend>data_length){bitend=data_length;}
bitlength=bitend-bitstart;

if(bit%10==0)
{
printf("Performing linear regression for Chunk %d of %d\n", bit+1, bittotal);
fclose(output);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Performing linear regression for Chunk %d of %d\n", bit+1, bittotal);
}

fclose(output2);
if((output2=fopen(filename2,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename2);exit(1);}
fclose(output3);
if((output3=fopen(filename3,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename3);exit(1);}
fclose(output4);
if((output4=fopen(filename4,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename4);exit(1);}
fclose(output5);
if((output5=fopen(filename5,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename5);exit(1);}

//read three sets of data for chunk
current=read_data_fly(datafile, dtype, data, NULL, 3*num_samples_use, keepsamps, bitstart, bitend, keeppreds_use, datainputgz, current, num_samples, num_preds, genskip, genheaders, genprobs, bedbytes, missingvalue, -9999, -9999, nonsnp, maxthreads);

//set missing to mean for three sets of data individually (for first set, record means)
impute_matrix_missing(data, num_samples_use, 3*num_samples_use, bitlength, missingvalue, centres+bitstart);
impute_matrix_missing(data+num_samples_use, num_samples_use, 3*num_samples_use, bitlength, missingvalue, NULL);
impute_matrix_missing(data+2*num_samples_use, num_samples_use, 3*num_samples_use, bitlength, missingvalue, NULL);

//regress covariates out of the three sets of data
reg_covar_matrix(data, Z, num_samples_use, 3*bitlength, num_fixed);

//get YTdata for all three sets of data
token=3*bitlength;
alpha=1.0;beta=0.0;
dgemv_("T", &num_samples_use, &token, &alpha, data, &num_samples_use, Yadj, &one, &beta, YTdata, &one);

//get XTXC and inverse for all predictors (set centres to -9999 if cholesky fails)
wcount=0;
for(j=0;j<bitlength;j++)
{
//get XTCX
alpha=1.0;beta=0.0;
dgemm_("T", "N", &three, &three, &num_samples_use, &alpha, data+(size_t)3*j*num_samples_use, &num_samples_use, data+(size_t)3*j*num_samples_use, &num_samples_use, &beta, XTCX+9*j, &three);

//copy then invert
for(k=0;k<9;k++){XTCX2[k+9*j]=XTCX[k+9*j];}
(void)cholesky_invert(XTCX2+9*j, 3, -1, NULL, &info, 0);

if(info!=0)
{
if(wcount<5){printf("Warning, Predictor %s will be ignored (due to lack of variation and/or colinearity with parental predictors)\n", preds[bitstart+j]);}
centres[bitstart+j]=-9999;
wcount++;
}

//fill XTCX3 (offspring-fathers) and invert
XTCX3[0+4*j]=XTCX[0+9*j];
XTCX3[1+4*j]=XTCX[1+9*j];
XTCX3[2+4*j]=XTCX[3+9*j];
XTCX3[3+4*j]=XTCX[4+9*j];
(void)cholesky_invert(XTCX3+4*j, 2, -1, NULL, &info, 0);

//fill XTCX4 (offspring-fathers) and invert
XTCX4[0+4*j]=XTCX[0+9*j];
XTCX4[1+4*j]=XTCX[2+9*j];
XTCX4[2+4*j]=XTCX[6+9*j];
XTCX4[3+4*j]=XTCX[8+9*j];
(void)cholesky_invert(XTCX4+4*j, 2, -1, NULL, &info, 0);
}
if(wcount>5){printf("In total, %d predictors will be ignored\n", wcount);}
if(wcount>0){printf("\n");}

////////

//ready to test
#pragma omp parallel for private(j, i, mark, mean, mean2, mean3, var, var2, var3) schedule(static)
for(j=0;j<bitlength;j++)
{
if(centres[bitstart+j]!=-9999)	//cholesky worked
{
//test basic model (just offspring)
mark=4*j;

//effect size estimate is inv(XTX) XTY
mean=YTdata[3*j]/XTCX[9*j];

//fill residuals
for(i=0;i<num_samples_use;i++){data2[(size_t)j*num_samples_use+i]=Yadj[i]-mean*data[(size_t)3*j*num_samples_use+i];}

//compute variance
var=club_sandwich_one(XTCX[9*j], data+(size_t)3*j*num_samples_use, data2+(size_t)j*num_samples_use, num_fams, famindex, num_samples_use, num_fixed);

stats[0+mark]=mean;
stats[1+mark]=pow(var,.5);
stats[2+mark]=stats[0+mark]/stats[1+mark];
stats[3+mark]=erfc(fabs(stats[2+mark])*M_SQRT1_2);

//test full model (offspring, father and mother)
mark=12*j;

//effect size estimates are inv(XTX) XTY
mean=XTCX2[0+9*j]*YTdata[3*j]+XTCX2[1+9*j]*YTdata[3*j+1]+XTCX2[2+9*j]*YTdata[3*j+2];
mean2=XTCX2[3+9*j]*YTdata[3*j]+XTCX2[4+9*j]*YTdata[3*j+1]+XTCX2[5+9*j]*YTdata[3*j+2];
mean3=XTCX2[6+9*j]*YTdata[3*j]+XTCX2[7+9*j]*YTdata[3*j+1]+XTCX2[8+9*j]*YTdata[3*j+2];

//fill residuals
for(i=0;i<num_samples_use;i++)
{data2[(size_t)j*num_samples_use+i]=Yadj[i]-mean*data[(size_t)3*j*num_samples_use+i]-mean2*data[(size_t)(3*j+1)*num_samples_use+i]-mean3*data[(size_t)(3*j+2)*num_samples_use+i];}

//compute variances
club_sandwich_three(XTCX2+9*j, data+(size_t)3*j*num_samples_use, data2+(size_t)j*num_samples_use, &var, &var2, &var3, num_fams, famindex, num_samples_use, num_fixed);

stats2[0+mark]=mean;
stats2[1+mark]=pow(var,.5);
stats2[2+mark]=stats2[0+mark]/stats2[1+mark];
stats2[3+mark]=erfc(fabs(stats2[2+mark])*M_SQRT1_2);

stats2[4+mark]=mean2;
stats2[5+mark]=pow(var2,.5);
stats2[6+mark]=stats2[4+mark]/stats2[5+mark];
stats2[7+mark]=erfc(fabs(stats2[6+mark])*M_SQRT1_2);

stats2[8+mark]=mean3;
stats2[9+mark]=pow(var3,.5);
stats2[10+mark]=stats2[8+mark]/stats2[9+mark];
stats2[11+mark]=erfc(fabs(stats2[10+mark])*M_SQRT1_2);

//test offspring-father model
mark=8*j;

//effect size estimates are inv(XTX) XTY
mean=XTCX3[0+4*j]*YTdata[3*j]+XTCX3[1+4*j]*YTdata[3*j+1];
mean2=XTCX3[2+4*j]*YTdata[3*j]+XTCX3[3+4*j]*YTdata[3*j+1];

//fill residuals
for(i=0;i<num_samples_use;i++)
{data2[(size_t)j*num_samples_use+i]=Yadj[i]-mean*data[(size_t)3*j*num_samples_use+i]-mean2*data[(size_t)(3*j+1)*num_samples_use+i];}

//compute variances
club_sandwich_two(XTCX3+4*j, data+(size_t)3*j*num_samples_use, data2+(size_t)j*num_samples_use, &var, &var2, num_fams, famindex, num_samples_use, num_fixed, 0);

stats3[0+mark]=mean;
stats3[1+mark]=pow(var,.5);
stats3[2+mark]=stats3[0+mark]/stats3[1+mark];
stats3[3+mark]=erfc(fabs(stats3[2+mark])*M_SQRT1_2);

stats3[4+mark]=mean2;
stats3[5+mark]=pow(var2,.5);
stats3[6+mark]=stats3[4+mark]/stats3[5+mark];
stats3[7+mark]=erfc(fabs(stats3[6+mark])*M_SQRT1_2);

//test offspring-mother model
mark=8*j;

//effect size estimates are inv(XTX) XTY
mean=XTCX4[0+4*j]*YTdata[3*j]+XTCX4[1+4*j]*YTdata[3*j+2];
mean2=XTCX4[2+4*j]*YTdata[3*j]+XTCX4[3+4*j]*YTdata[3*j+2];

//fill residuals
for(i=0;i<num_samples_use;i++)
{data2[(size_t)j*num_samples_use+i]=Yadj[i]-mean*data[(size_t)3*j*num_samples_use+i]-mean2*data[(size_t)(3*j+2)*num_samples_use+i];}

//compute variances
club_sandwich_two(XTCX4+4*j, data+(size_t)3*j*num_samples_use, data2+(size_t)j*num_samples_use, &var, &var2, num_fams, famindex, num_samples_use, num_fixed, 1);

stats4[0+mark]=mean;
stats4[1+mark]=pow(var,.5);
stats4[2+mark]=stats4[0+mark]/stats4[1+mark];
stats4[3+mark]=erfc(fabs(stats4[2+mark])*M_SQRT1_2);

stats4[4+mark]=mean2;
stats4[5+mark]=pow(var2,.5);
stats4[6+mark]=stats4[4+mark]/stats4[5+mark];
stats4[7+mark]=erfc(fabs(stats4[6+mark])*M_SQRT1_2);
}
}	//end of j loop

//save results
for(j=0;j<bitlength;j++)
{
if(centres[bitstart+j]!=-9999)	//have results
{
//basic model
mark=4*j;
fprintf(output2, "%d %s %.0f %c %c ", chr[bitstart+j], preds[bitstart+j], cmbp[bitstart+j], al1[bitstart+j], al2[bitstart+j]);
fprintf(output2, "%.4f %.4e %.4e %.4e ", stats[2+mark], stats[3+mark], stats[0+mark], stats[1+mark]);
fprintf(output2, "%.6f ", centres[bitstart+j]);
if(nonsnp==0){fprintf(output2, "%.6f\n", centres[bitstart+j]/2+(centres[bitstart+j]>1)*(1-centres[bitstart+j]));}
else{fprintf(output2, "NA\n");}

//full model
mark=12*j;
fprintf(output3, "%d %s %.0f %c %c ", chr[bitstart+j], preds[bitstart+j], cmbp[bitstart+j], al1[bitstart+j], al2[bitstart+j]);
fprintf(output3, "%.4f %.4e %.4e %.4e ", stats2[0+mark], stats2[1+mark], stats2[2+mark], stats2[3+mark]);
fprintf(output3, "%.4f %.4e %.4e %.4e ", stats2[4+mark], stats2[5+mark], stats2[6+mark], stats2[7+mark]);
fprintf(output3, "%.4f %.4e %.4e %.4e\n", stats2[8+mark], stats2[9+mark], stats2[10+mark], stats2[11+mark]);

//first offspring model
mark=8*j;
fprintf(output4, "%d %s %.0f %c %c ", chr[bitstart+j], preds[bitstart+j], cmbp[bitstart+j], al1[bitstart+j], al2[bitstart+j]);
fprintf(output4, "%.4f %.4e %.4e %.4e ", stats3[0+mark], stats3[1+mark], stats3[2+mark], stats3[3+mark]);
fprintf(output4, "%.4f %.4e %.4e %.4e\n", stats3[4+mark], stats3[5+mark], stats3[6+mark], stats3[7+mark]);

//second offspring model
mark=8*j;
fprintf(output5, "%d %s %.0f %c %c ", chr[bitstart+j], preds[bitstart+j], cmbp[bitstart+j], al1[bitstart+j], al2[bitstart+j]);
fprintf(output5, "%.4f %.4e %.4e %.4e ", stats4[0+mark], stats4[1+mark], stats4[2+mark], stats4[3+mark]);
fprintf(output5, "%.4f %.4e %.4e %.4e\n", stats4[4+mark], stats4[5+mark], stats4[6+mark], stats4[7+mark]);
}
else	//trivial
{
fprintf(output2, "%d %s %.0f %c %c ", chr[bitstart+j], preds[bitstart+j], cmbp[bitstart+j], al1[bitstart+j], al2[bitstart+j]);
fprintf(output2, "NA NA NA NA NA NA\n");

fprintf(output3, "%d %s %.0f %c %c ", chr[bitstart+j], preds[bitstart+j], cmbp[bitstart+j], al1[bitstart+j], al2[bitstart+j]);
fprintf(output3, "NA NA NA NA NA NA NA NA NA NA NA NA NA NA\n");

fprintf(output4, "%d %s %.0f %c %c ", chr[bitstart+j], preds[bitstart+j], cmbp[bitstart+j], al1[bitstart+j], al2[bitstart+j]);
fprintf(output4, "NA NA NA NA NA NA NA NA NA NA\n");

fprintf(output5, "%d %s %.0f %c %c ", chr[bitstart+j], preds[bitstart+j], cmbp[bitstart+j], al1[bitstart+j], al2[bitstart+j]);
fprintf(output5, "NA NA NA NA NA NA NA NA NA NA\n");
}
}	//end of j loop
}	//end of bit loop
if(wcount==0){printf("\n");}

fclose(output);
fclose(output2);
fclose(output3);

printf("Results saved in %s, %s, %s and %s\n\n", filename2, filename3, filename4, filename5);

//write details
count=0;for(j=0;j<data_length;j++){count+=(centres[j]!=-9999);}

sprintf(filename,"%s.details", outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fprintf(output,"Analysis Trios\nNum_Samples %d\n", num_samples_use);
fprintf(output,"Num_Predictors %d\nNum_Tested %d\n", data_length, count);
fclose(output);

////////

free(data);free(data2);
free(order);free(Y);free(Z);
free(thetas);free(thetasds);free(thetapvas);free(Yadj);
free(stats);free(stats2);
free(YTdata);free(XTCX);free(XTCX2);
if(binary==0){gzclose(datainputgz);}

///////////////////////////

