/*
Copyright 2026 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Join results from running sub-step

///////////////////////////

//allocations

Umat2=malloc(sizeof(double)*num_resps_use*num_resps_use);
effs=malloc(sizeof(double)*data_length*num_resps_use);
effs2=malloc(sizeof(double)*data_length*num_resps_use);

Y=malloc(sizeof(double)*num_samples_use*num_resps_use);
Yadj=malloc(sizeof(double)*num_samples_use*num_resps_use);
Yadj2=malloc(sizeof(double)*num_samples_use*num_resps_use);
residuals=malloc(sizeof(double)*num_samples_use*num_resps_use);
residuals2=malloc(sizeof(double)*num_samples_use*num_resps_use);

Mscales=malloc(sizeof(double)*num_resps_use);
Mneffs=malloc(sizeof(double)*num_resps_use);

readdoubles=malloc(sizeof(double)*7);

//read inverse transformation

sprintf(filename,"%s.multivariate.inverse", outfile);
if(just_check(filename)!=0)
{printf("Error reading %s; this file should have been created using \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\" with \"--sub-step CUT\"\n\n", filename);exit(1);}

count=countrows(filename);
count2=countcols(filename);
if(count!=num_resps_use||count2!=num_resps_use)
{printf("Error, %s should have %d rows and %d columns (not %d and %d), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\" with \"--sub-step CUT\"\n\n", filename, num_resps_use, num_resps_use, count, count2);exit(1);}

for(m=0;m<num_resps_use;m++){read_values(filename, Umat2+m*num_resps_use, num_resps_use, NULL, m+1, 0, 0);}

//check all files complete and correct size (and get num_test and Mscales)

for(m=0;m<num_resps_use;m++)
{
sprintf(filename,"%s.pheno%d.effects.temp", outfile, m+1);
if(just_check(filename)!=0)
{printf("Error reading %s; this file should have been created using \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\" with \"--sub-step %d\"\n\n", filename, m+1);exit(1);}

count=countcols(filename);
if(count!=6)
{printf("Error, %s should have six columns (not %d), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\" with \"--sub-step %d\"\n\n", filename, count, m+1);exit(1);}

count2=countrows(filename)-1;
if(count2!=data_length){printf("Error, the number of predictors in %s (%d) does not match the number of predictors being analyzed (%d)\n\n", filename, count2, data_length);exit(1);}

sprintf(filename,"%s.pheno%d.predictions.temp", outfile, m+1);
if(just_check(filename)!=0)
{printf("Error reading %s; this file should have been created using \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\" with \"--sub-step %d\"\n\n", filename, m+1);exit(1);}

count=countcols(filename);
if(count!=3)
{printf("Error, %s should have three columns (not %d), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\" with \"--sub-step %d\"\n\n", filename, count, m+1);exit(1);}

count2=countrows(filename);
if(m==0){num_test=count2;}
else
{
if(count2!=num_test){printf("Error, the number of samples in %s (%d) does not match the number in %s.pheno1.predictions.temp (%d)\n\n", filename, count2, outfile, num_test);exit(1);}
}

if(loco==1)
{
sprintf(filename,"%s.pheno%d.loco.details.temp", outfile, m+1);
if(just_check(filename)!=0)
{printf("Error reading %s; this file should have been created using \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\" with \"--sub-step %d\"\n\n", filename, m+1);exit(1);}

count=countcols(filename);
if(count!=2)
{printf("Error, %s should have two columns (not %d), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\" with \"--sub-step %d\"\n\n", filename, count, m+1);exit(1);}

count=countrows(filename);
if(count!=8)
{printf("Error, %s should have eight rows (not %d), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\" with \"--sub-step %d\"\n\n", filename, count, m+1);exit(1);}

read_values(filename, Mscales+m, 1, NULL, 2, 7, 0);
}
}

////////

//untransform genome-wide effects, and correct for padded values
printf("Transforming effect sizes for %d phenotypes\n\n", num_resps_use);

for(m=0;m<num_resps_use;m++)
{
sprintf(filename,"%s.pheno%d.effects.temp", outfile, m+1);
read_values(filename, effs2+m*data_length, data_length, NULL, 5, 1, 0);
}

//read in centres
sprintf(filename,"%s.pheno1.effects.temp", outfile);
read_values(filename, centres, data_length, NULL, 4, 1, 0);

//read in mults (well, qc)
sprintf(filename,"%s.pheno1.effects.temp", outfile);
read_values(filename, mults, data_length, NULL, 6, 1, 0);

//transform effects
alpha=1.0;beta=0.0;
dgemm_("N", "T", &data_length, &num_resps_use, &num_resps_use, &alpha, effs2, &data_length, Umat2, &num_resps_use, &beta, effs, &data_length);

//correct for padded values
for(m=0;m<num_resps_use;m++)
{
if(respcounts[m]<num_samples_use)
{
//divide effects by p
value=(double)num_samples_use/respcounts[m];
for(j=0;j<data_length;j++){effs[m*data_length+j]*=value;}
}
}

//save effects, remembering to scale
for(m=0;m<num_resps_use;m++)
{
sprintf(filename2,"%s.pheno%d.effects", outfile, m+1);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output2,"Predictor A1 A2 Centre Effect\n");

for(j=0;j<data_length;j++)
{
if(mults[j]==1){fprintf(output2, "%s %c %c %.6f %.4e\n", preds[j], al1[j], al2[j], centres[j], effs[m*data_length+j]/Mscales[m]);}
}
fclose(output2);
}

//get mse of untransformed models, and then update effective sample size 
printf("Measuring the accuracy of the best models for the original phenotypes\n");

//read in phenotypes and predictions
for(m=0;m<num_resps_use;m++)
{
sprintf(filename,"%s.pheno%d.predictions.temp", outfile, m+1);
read_values(filename, Y+m*num_test, num_test, NULL, 1, 0, 1);
read_values(filename, Yadj+m*num_test, num_test, NULL, 2, 0, 0);
read_values(filename, residuals+m*num_test, num_test, NULL, 3, 0, 0);
}

//untransform adjusted phenotypes and predictions
alpha=1.0;beta=0.0;
dgemm_("N", "T", &num_test, &num_resps_use, &num_resps_use, &alpha, Yadj, &num_test, Umat2, &num_resps_use, &beta, Yadj2, &num_test);
dgemm_("N", "T", &num_test, &num_resps_use, &num_resps_use, &alpha, residuals, &num_test, Umat2, &num_resps_use, &beta, residuals2, &num_test);

for(m=0;m<num_resps_use;m++)
{
if(respcounts[m]<num_samples_use)	//correct for padded values - get Y - Xbeta/p = Y - (Y-R)/p = (1-1/p) Y + R/p
{
value=1-(double)num_samples_use/respcounts[m];
value2=(double)num_samples_use/respcounts[m];
for(i=0;i<num_test;i++){residuals2[(size_t)m*num_test+i]=value*Yadj2[i+m*num_test]+value2*residuals2[(size_t)m*num_test+i];}
}

//compute relative mse (across non-missing individuals)
sumsq=0;sumsq2=0;
for(i=0;i<num_test;i++)
{
if(Y[i+m*num_test]!=missingvalue)
{
sumsq+=pow(residuals2[(size_t)m*num_test+i],2);
sumsq2+=pow(Yadj2[i+m*num_test],2);
}
}

value=sumsq/sumsq2;
printf("Phenotype %d: mean squared error %.4f\n", m+1, value);
Mneffs[m]=respcounts[m]/value;
}
printf("\n");

if(loco==1)	//save loco.details, updating effective size
{
for(m=0;m<num_resps_use;m++)
{
sprintf(filename,"%s.pheno%d.loco.details.temp", outfile, m+1);
read_values(filename, readdoubles, 7, NULL, 2, 0, 1);

sprintf(filename2,"%s.pheno%d.loco.details", outfile, m+1);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2,"Actual_Sample_Size %d\n", (int)readdoubles[0]);
fprintf(output2,"Approx_Effective_Sample_Size %.1f\n", Mneffs[m]);
fprintf(output2,"Scaling_Estimate %.4f\n", readdoubles[2]);
fprintf(output2,"Scaling_SE %.4f\n", readdoubles[3]);
fprintf(output2,"Effect_Calibration %.4f\n", readdoubles[4]);
fprintf(output2,"Power %.4f\n", readdoubles[5]);
fprintf(output2,"Heritability %.4f\n", readdoubles[6]);
fclose(output2);
}
}

printf("Best-fitting models saved in %s.phenoX.effects, where X is the phenotype number\n\n", outfile);

if(loco==1)
{printf("LOCO results saved in %s.phenoX.loco.prs and %s.phenoX.loco.prs, where X is the phenotype number\n\n", outfile, outfile);}

printf("Note that you can now delete the files with suffix \".temp\"\n\n");

free(Umat2);free(effs);free(effs2);
free(Y);free(Yadj);free(Yadj2);free(residuals);free(residuals2);
free(Mscales);free(Mneffs);
free(readdoubles);

///////////////////////////

