/*
Copyright 2026 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Jackknife to get variances

///////////////////////////

//cflag indicates whether using cv
cflag=(cvprop!=-9999||strcmp(bvsfile,"blank")!=0);

//get number of datapoints and set num_blocks (have already checked jackfile has at least two columns)
if(strcmp(jackfile,"blank")!=0)
{count=countcols(jackfile);count2=countrows(jackfile);}
else
{count=(countcols(proffile)-4)/2;count2=countrows(proffile)-1;}

if(count2<3){printf("Error, there are only %d datapoints, it is not possible to continue\n\n", count2);exit(1);}
if(count2<200){printf("Warning, there are usually thousands of datapoints (not %d)\n\n", count2);}

if(strcmp(jackfile,"blank")!=0){num_covars=count-1;}
else{num_covars=1;}

if(num_blocks==-1){num_blocks=count2;}

//allocate variables
order=malloc(sizeof(int)*count2);
sX=malloc(sizeof(double)*count2*3);
sZ=malloc(sizeof(double)*count2*(num_covars+1));
sT=malloc(sizeof(double)*count2);
sXTX=malloc(sizeof(double)*9);
sXTX2=malloc(sizeof(double)*9);
stats=malloc(sizeof(double)*num_blocks*4);
likes=malloc(sizeof(double)*num_blocks);

if(strcmp(proffile,"blank")!=0)	//using a profile file, so can be quite big
{
guesses=malloc(sizeof(double*)*(count+cflag));
for(k=0;k<count+cflag;k++){guesses[k]=malloc(sizeof(double)*count2);}
}

if(dichot==1)
{
dptrs=malloc(sizeof(struct sorting_double)*count2);
sR1=malloc(sizeof(int)*(count2+1));
sR2=malloc(sizeof(int)*(count2+1));
}

//open results file
sprintf(filename,"%s.jack",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s - check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
if(strcmp(jackfile,"blank")!=0){fprintf(output, "Measure Estimate SE\n");}
else{fprintf(output, "Profile Measure Estimate SE\n");}

////////

if(strcmp(jackfile,"blank")!=0)	//using a jack file
{
//load up predicted and observed values
if(count==2){printf("Reading predicted and observed values for %d samples from %s\n\n", count2, jackfile);}
else{printf("Reading predicted and observed values for %d samples, as well as %d covariates, from %s\n\n", count2, num_covars-1, jackfile);}
read_values(jackfile,sX,count2,NULL,1,0,0);
read_values(jackfile,sX+count2,count2,NULL,2,0,0);

//read covariates
for(j=0;j<count2;j++){sZ[j]=1;}
if(count>2)
{
for(k=1;k<num_covars;k++){read_values(jackfile,sZ+k*count2,count2,NULL,2+k,0,0);}
}

//shuffle the values
for(j=0;j<count2;j++){order[j]=j;}
permute_int(order,count2);

for(k=0;k<2;k++)
{
for(j=0;j<count2;j++){sT[j]=sX[j+k*count2];}
for(j=0;j<count2;j++){sX[j+k*count2]=sT[order[j]];}
}
for(k=1;k<num_covars;k++)
{
for(j=0;j<count2;j++){sT[j]=sZ[j+k*count2];}
for(j=0;j<count2;j++){sZ[j+k*count2]=sT[order[j]];}
}

//set third column of sX to one and final column of sZ to predictions
for(j=0;j<count2;j++){sX[j+2*count2]=1;sZ[j+num_covars*count2]=sX[j];}

if(dichot==1||prev!=-9999)	//check observed values are zero or one (but not all same) and sort
{
total=0;
for(j=0;j<count2;j++)
{
if(sX[j+count2]!=0&&sX[j+count2]!=1)
{printf("Error reading %s; values in Column 2 should be either one or zero (corresponding to cases and controls)\n\n", jackfile);exit(1);}
total+=sX[j+count2];
}
if(total==0){printf("Error reading %s; all values in Column 2 are zero\n\n", jackfile);exit(1);}
if(total==count2){printf("Error reading %s; all values in Column 2 are one\n\n", jackfile);exit(1);}

printf("The phenotype is binary, with %d cases and %d controls\n\n", total, count2-total);

if(prev!=-9999){factor=get_factor(sX+count2, count2, prev, -9999, outfile);}
}

//ready to analyse

if(num_blocks>count2){num_blocks=count2;}

//get correlation across all datapoints
alpha=1.0;beta=0.0;
dgemm_("T", "N", &three, &three, &count2, &alpha, sX, &count2, sX, &count2, &beta, sXTX, &three);
sum=sXTX[2];sum2=sXTX[5];
sumsq=sXTX[0];sumsq2=sXTX[4];sumsq3=sXTX[1];
value=(count2*sumsq3-sum*sum2)*pow(count2*sumsq-sum*sum,-.5)*pow(count2*sumsq2-sum2*sum2,-.5);

//store var of observed values in value5 - will use this to scale mse and mae
value5=sumsq2/count2-pow(sum2/count2,2);

//now mse and mae
sum=0;for(j=0;j<count2;j++){sum+=pow(sX[j]-sX[j+count2],2);}
value2=sum/count2;
sum=0;for(j=0;j<count2;j++){sum+=fabs(sX[j]-sX[j+count2]);}
value3=sum/count2;

if(dichot==1)	//now auc and likelihood (must first sort)
{
for(j=0;j<count2;j++){dptrs[j].value=sX[j];dptrs[j].index=j;}
qsort(dptrs, count2, sizeof(struct sorting_double), compare_sorting_double);

sR1[0]=0;sR2[0]=0;found=1;
last=dptrs[0].value-1;  //ensures first value visited always treated as new
for(j=0;j<count2;j++)
{
j2=dptrs[j].index;
if(dptrs[j].value==last)	//same as previous value
{sR1[found-1]+=sX[j2+count2];sR2[found-1]+=1-sX[j2+count2];}
else	//new
{sR1[found]=sR1[found-1]+sX[j2+count2];sR2[found]=sR2[found-1]+1-sX[j2+count2];found++;}
last=dptrs[j].value;
}

sum=0;for(j=1;j<found;j++){sum+=.5*(sR1[j]-sR1[j-1])*(sR2[j]+sR2[j-1]);}
value4=sum/sR1[found-1]/sR2[found-1];

value6=like_diff_log(sX+count2, sZ, count2, num_covars+1, 0.001, 100, 0, 0);
}

//now jackknife
for(p=0;p<num_blocks;p++)
{
start=(double)(p)*count2/num_blocks;
end=(double)(p+1)*count2/num_blocks;
if(p%500000==0){printf("Performing Jackknife %d out of %d\n", p+1, num_blocks);}
count3=end-start;

//first correlation
alpha=1.0;beta=0.0;
dgemm_("T", "N", &three, &three, &count3, &alpha, sX+start, &count2, sX+start, &count2, &beta, sXTX2, &three);
sum=sXTX[2]-sXTX2[2];sum2=sXTX[5]-sXTX2[5];
sumsq=sXTX[0]-sXTX2[0];sumsq2=sXTX[4]-sXTX2[4];sumsq3=sXTX[1]-sXTX2[1];
stats[p]=((count2-count3)*sumsq3-sum*sum2)*pow((count2-count3)*sumsq-sum*sum,-.5)*pow((count2-count3)*sumsq2-sum2*sum2,-.5);

//now mse and mae
sum=0;for(j=start;j<end;j++){sum+=pow(sX[j]-sX[j+count2],2);}
stats[p+num_blocks]=(value2*count2-sum)/(count2-count3);
sum=0;for(j=start;j<end;j++){sum+=fabs(sX[j]-sX[j+count2]);}
stats[p+2*num_blocks]=(value3*count2-sum)/(count2-count3);

if(dichot==1)	//now auc and likelihood
{
//for auc, no need to resort
sR1[0]=0;sR2[0]=0;found=1;
last=dptrs[0].value-1;
for(j=0;j<count2;j++)
{
j2=dptrs[j].index;
if(j2<start||j2>=end)
{
if(dptrs[j].value==last)	//same as previous value
{sR1[found-1]+=sX[j2+count2];sR2[found-1]+=1-sX[j2+count2];}
else	//new
{sR1[found]=sR1[found-1]+sX[j2+count2];sR2[found]=sR2[found-1]+1-sX[j2+count2];found++;}
last=dptrs[j].value;
}
}

sum=0;for(j=1;j<found;j++){sum+=.5*(sR1[j]-sR1[j-1])*(sR2[j]+sR2[j-1]);}
stats[p+3*num_blocks]=sum/sR1[found-1]/sR2[found-1];

likes[p]=like_diff_log(sX+count2, sZ, count2, num_covars+1, 0.001, 100, start, end);
}
}	//end of p loop
printf("\n");

sum=0;sumsq=0;
for(p=0;p<num_blocks;p++){sum+=stats[p];sumsq+=pow(stats[p],2);}
mean=sum/num_blocks;
var=(num_blocks-1)*(sumsq/num_blocks-pow(mean,2));
printf("Correlation is %.4f (SE %.4f)\n", value, pow(var,.5));
fprintf(output, "Correlation %.6f %.6f\n", value, pow(var,.5));

sum=0;sumsq=0;
for(p=0;p<num_blocks;p++){sum+=pow(stats[p],2);sumsq+=pow(stats[p],4);}
mean=sum/num_blocks;
var=(num_blocks-1)*(sumsq/num_blocks-pow(mean,2));
printf("Correlation squared is %.4f (SE %.4f)\n", pow(value,2), pow(var,.5));
fprintf(output, "Squared_correlation %.6f %.6f\n", pow(value,2), pow(var,.5));

sum=0;sumsq=0;
for(p=0;p<num_blocks;p++){sum+=stats[p+num_blocks];sumsq+=pow(stats[p+num_blocks],2);}
mean=sum/num_blocks;
var=(num_blocks-1)*(sumsq/num_blocks-pow(mean,2));
printf("Mean squared error is %.4f (SE %.4f)\n", value2/value5, pow(var,.5)/value5);
fprintf(output, "Mean_squared_error %.6f %.6f\n", value2/value5, pow(var,.5)/value5);

sum=0;sumsq=0;
for(p=0;p<num_blocks;p++){sum+=stats[p+2*num_blocks];sumsq+=pow(stats[p+2*num_blocks],2);}
mean=sum/num_blocks;
var=(num_blocks-1)*(sumsq/num_blocks-pow(mean,2));
printf("Mean absolute error is %.4f (SE %.4f)\n", value3*pow(value5,-.5), pow(var,.5)*pow(value5,-.5));
fprintf(output, "Mean_absolute_error %.6f %.6f\n", value3*pow(value5,-.5), pow(var,.5)*pow(value5,-.5));

if(dichot==1)	//now auc
{
sum=0;sumsq=0;
for(p=0;p<num_blocks;p++){sum+=stats[p+3*num_blocks];sumsq+=pow(stats[p+3*num_blocks],2);}
mean=sum/num_blocks;
var=(num_blocks-1)*(sumsq/num_blocks-pow(mean,2));
printf("Area under curve is %.4f (SE %.4f)\n", value4, pow(var,.5));
fprintf(output, "Area_under_curve %.6f %.6f\n", value4, pow(var,.5));

sum=0;sumsq=0;
for(p=0;p<num_blocks;p++){sum+=likes[p];sumsq+=pow(likes[p],2);}
mean=sum/num_blocks;
var=(num_blocks-1)*(sumsq/num_blocks-pow(mean,2));
printf("Nagelkerke R2 is %.4f (SE %.4f)\n", value6, pow(var,.5));
fprintf(output, "Nagelkerke_R2 %.6f %.6f\n", value6, pow(var,.5));
}

if(prev!=-9999)	//now liability correlation
{
sum=0;sumsq=0;
for(p=0;p<num_blocks;p++){sum+=pow(stats[p],2);sumsq+=pow(stats[p],4);}
mean=sum/num_blocks;
var=(num_blocks-1)*(sumsq/num_blocks-pow(mean,2));
printf("Correlation squared on the liability scale is %.4f (SE %.4f)\n", pow(value,2)*factor, pow(var,.5)*factor);
fprintf(output, "Liability_squared_correlation %.6f %.6f\n", pow(value,2)*factor, pow(var,.5)*factor);
}
printf("\n");
}
else	//using profile file - might have missing phenotypes
{
//read observed values (minus covariates) into second column of sX, and predictions into first columns of guesses
printf("Reading observed and predicted phenotypes for %d samples from %s\n", count2, proffile);
if((input=fopen(proffile,"r"))==NULL)
{printf("Error opening %s\n\n",proffile);exit(1);}
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}

count4=0;
wcount=0;
for(j=0;j<count2;j++)
{
if(fscanf(input, "%s %s %s %s", readstring, readstring2, readstring3, readstring4)!=4)
{printf("Error reading first four values of Row %d of %s\n\n", j+1, proffile);exit(1);}
if(strcmp(readstring3,"NA")!=0)	//have phenotype - will use sample
{
sX[count4+count2]=atof(readstring3)-atof(readstring4); 
for(k=0;k<count;k++)
{
if(fscanf(input, "%lf %s ", guesses[k]+count4, readstring)!=2)
{printf("Error reading Profile %d from Row %d of %s\n\n", k+1, j+1, proffile);exit(1);}
}
count4++;
}
else	//missing phenotype - will skip sample
{
if(wcount<5){printf("Warning, phenotype is missing for Sample %s %s\n", readstring, readstring2);}
wcount++;
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}
}
}
fclose(input);

if(wcount>5){printf("In total, phenotypes are missing for %d samples\n", wcount);}
printf("\n");

if(count4==0){printf("Error, no samples have phenotypes\n\n");exit(1);}
if(count4<3){printf("Error, only %d samples have phenotypes\n\n", count4);exit(1);}

if(dichot==1||prev!=-9999)	//check observed values are zero or one (but not all same) and sort
{
total=0;
for(j=0;j<count4;j++)
{
if(sX[j+count2]!=0&&sX[j+count2]!=1)
{printf("Error reading %s; phenotypes should be either one or zero (corresponding to cases and controls)\n\n", proffile);exit(1);}
total+=sX[j+count2];
}
if(total==0){printf("Error reading %s; all phenotypes are zero\n\n", jackfile);exit(1);}
if(total==count4){printf("Error reading %s; all phenotypes are one\n\n", jackfile);exit(1);}

printf("The phenotype is binary, with %d cases and %d controls\n\n", total, count4-total);

if(prev!=-9999){factor=get_factor(sX+count2, count4, prev, -9999, outfile);}
}

if(cflag==1)  //regress phenotypes on profiles (plus intercept), then save prediction in final column of guesses
{
sP=malloc(sizeof(double)*count4*(count+1));
sP2=malloc(sizeof(double)*(count+1));

if(cvprop!=-9999)   //simply use first cvprop individuals to estimate weights
{
count3=(int)(cvprop*count4);
if(count3==0)   //a trick to allow doug to use all individuals for fitting and estimating
{count3=count4;}
}
else    //work out how many cv samples, then put these first
{
//read cvsamples
count5=countrows(bvsfile);
printf("Reading list of %d cross-validation samples from %s\n", count5, bvsfile);
wantids=malloc(sizeof(char*)*count5);
read_ids(bvsfile, NULL, NULL, wantids, count5, NULL, 0, 0);

//read profile samples and phenotypes
wantids2=malloc(sizeof(char*)*count2);
readdoubles=malloc(sizeof(double)*count2);
read_ids(proffile, NULL, NULL, wantids2, count2, NULL, 1, 0);
read_values(proffile, readdoubles, count2, NULL, 3, 1, 1);

//remove samples missing phenotypes
found=0;
for(j=0;j<count2;j++)
{
if(readdoubles[j]!=-9999)
{
if(j!=found){free(wantids2[found]);copy_string(wantids2,found,wantids2[j]);}
found++;
}
}
if(found!=count4){printf("Error AZZ7, please tell Doug\n\n");exit(1);}

//find cvsamples
indexer=malloc(sizeof(int)*count4);
indexer2=malloc(sizeof(int)*count4);
count3=find_strings(wantids2, count4, wantids, count5, indexer, NULL, proffile, bvsfile, NULL, NULL, 3);
if(count3==0){printf("Error, no overlap between the %d samples in %s and the %d samples with phenotypes in %s\n\n", count5, bvsfile, count4, proffile);exit(1);}
if(count3==count4){printf("Error, %s contains all %d samples with phenotypes in %s\n\n", bvsfile, count4, proffile);exit(1);}

//now rearrange so that cvsamples are first
for(j=0;j<count4;j++){indexer2[j]=1;}
for(j=0;j<count3;j++){indexer2[indexer[j]]=0;}

found=count3;
for(j=0;j<count4;j++)
{
if(indexer2[j]==1){indexer[found]=j;found++;}
}
if(found!=count4){printf("Error BZZ7, please tell Doug\n\n");exit(1);}

for(j=0;j<count4;j++){sT[j]=sX[j+count2];}
for(j=0;j<count4;j++){sX[j+count2]=sT[indexer[j]];}
for(k=0;k<count;k++)
{
for(j=0;j<count4;j++){sT[j]=guesses[k][j];}
for(j=0;j<count4;j++){guesses[k][j]=sT[indexer[j]];}
}

for(j=0;j<count5;j++){free(wantids[j]);}free(wantids);
for(j=0;j<count4;j++){free(wantids2[j]);}free(wantids2);
free(indexer);free(indexer2);
}

//now do the regression
for(k=0;k<count;k++)
{
for(j=0;j<count3;j++){sP[j+k*count3]=guesses[k][j];}
}
for(j=0;j<count3;j++){sP[j+count*count3]=1;}

if(dichot==0)   //estimate coefficients
{reg_covar_lin(sX+count2, sP, count3, count+1, 0, sP2, NULL, NULL, NULL, -9999, NULL, NULL);}
else{reg_covar_log(sX+count2, sP, count3, count+1, 0, NULL, sP2, NULL, NULL, NULL, -9999, NULL, NULL, -9999, 0.001, 100);}

//combine across all datapoints (had considered converting to probabilities for binary traits)
for(j=0;j<count4;j++)
{
guesses[count][j]=sP2[count];
for(k=0;k<count;k++){guesses[count][j]+=sP2[k]*guesses[k][j];}
//if(dichot==1){guesses[count][j]=pow(1+exp(-guesses[count][j]),-1);}
}

//save scaled weights
sum=0;for(k=0;k<count;k++){sum+=sP2[k];}
sprintf(filename2,"%s.weights",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s - check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2, "Profile Weight\n");
for(k=0;k<count;k++){fprintf(output2, "%d %.4f\n", k+1, sP2[k]/sum);}
fclose(output2);

if(count3<count4)   //reduce datapoints and count4
{
printf("Reducing the number of datapoints from %d to %d\n\n", count4, count4-count3);
count4=count4-count3;
for(j=0;j<count4;j++)
{
sX[count2+j]=sX[count2+count3+j];
for(k=0;k<count+1;k++){guesses[k][j]=guesses[k][count3+j];}
}
}

free(sP);free(sP2);
}

//shuffle the values
for(j=0;j<count4;j++){order[j]=j;}
permute_int(order,count4);

for(j=0;j<count4;j++){sT[j]=sX[j+count2];}
for(j=0;j<count4;j++){sX[j+count2]=sT[order[j]];}
for(k=0;k<count+cflag;k++)
{
for(j=0;j<count4;j++){sT[j]=guesses[k][j];}
for(j=0;j<count4;j++){guesses[k][j]=sT[order[j]];}
}

//set third column of sX to one
for(j=0;j<count4;j++){sX[j+2*count2]=1;}

//ready to analyse - remember, we only have count4 datapoints (but sX has size count2)

if(num_blocks>count4){num_blocks=count4;}

for(k=0;k<count+cflag;k++)
{
if(k<count){printf("Analysing Profile %d\n", k+1);}
else{printf("Analysing combined profile\n");}

//load predicted values for Profile k into first column of sX
for(j=0;j<count4;j++){sX[j]=guesses[k][j];}

//get correlation across all datapoints
alpha=1.0;beta=0.0;
dgemm_("T", "N", &three, &three, &count4, &alpha, sX, &count2, sX, &count2, &beta, sXTX, &three);
sum=sXTX[2];sum2=sXTX[5];
sumsq=sXTX[0];sumsq2=sXTX[4];sumsq3=sXTX[1];
value=(count4*sumsq3-sum*sum2)*pow(count4*sumsq-sum*sum,-.5)*pow(count4*sumsq2-sum2*sum2,-.5);

//store variance of observed values in value5 - will use this to scale mse and mae
value5=sumsq2/count4-pow(sum2/count4,2);

//now mse and mae
sum=0;for(j=0;j<count4;j++){sum+=pow(sX[j]-sX[j+count2],2);}
value2=sum/count4;
sum=0;for(j=0;j<count4;j++){sum+=fabs(sX[j]-sX[j+count2]);}
value3=sum/count4;

if(dichot==1)	//now auc (must first sort)
{
for(j=0;j<count4;j++){dptrs[j].value=sX[j];dptrs[j].index=j;}
qsort(dptrs, count4, sizeof(struct sorting_double), compare_sorting_double);

sR1[0]=0;sR2[0]=0;found=1;
last=dptrs[0].value-1;
for(j=0;j<count4;j++)
{
j2=dptrs[j].index;
if(dptrs[j].value==last)	//same as previous value
{sR1[found-1]+=sX[j2+count2];sR2[found-1]+=1-sX[j2+count2];}
else	//new
{sR1[found]=sR1[found-1]+sX[j2+count2];sR2[found]=sR2[found-1]+1-sX[j2+count2];found++;}
last=dptrs[j].value;
}

sum=0;for(j=1;j<found;j++){sum+=.5*(sR1[j]-sR1[j-1])*(sR2[j]+sR2[j-1]);}
value4=sum/sR1[found-1]/sR2[found-1];
}	//end of dichot

//now jackknife

for(p=0;p<num_blocks;p++)
{
start=(double)(p)*count4/num_blocks;
end=(double)(p+1)*count4/num_blocks;
if(p%500000==0){printf("Performing Jackknife %d out of %d\n", p+1, num_blocks);}
count3=end-start;

//first correlation
alpha=1.0;beta=0.0;
dgemm_("T", "N", &three, &three, &count3, &alpha, sX+start, &count2, sX+start, &count2, &beta, sXTX2, &three);
sum=sXTX[2]-sXTX2[2];sum2=sXTX[5]-sXTX2[5];
sumsq=sXTX[0]-sXTX2[0];sumsq2=sXTX[4]-sXTX2[4];sumsq3=sXTX[1]-sXTX2[1];
stats[p]=((count4-count3)*sumsq3-sum*sum2)*pow((count4-count3)*sumsq-sum*sum,-.5)*pow((count4-count3)*sumsq2-sum2*sum2,-.5);

//now mse and mae
sum=0;for(j=start;j<end;j++){sum+=pow(sX[j]-sX[j+count2],2);}
stats[p+num_blocks]=(value2*count4-sum)/(count4-count3);
sum=0;for(j=start;j<end;j++){sum+=fabs(sX[j]-sX[j+count2]);}
stats[p+2*num_blocks]=(value3*count4-sum)/(count4-count3);

if(dichot==1)	//now auc
{
sR1[0]=0;sR2[0]=0;found=1;
last=dptrs[0].value-1;
for(j=0;j<count4;j++)
{
j2=dptrs[j].index;
if(j2<start||j2>=end)
{
if(dptrs[j].value==last)	//same as previous value
{sR1[found-1]+=sX[j2+count2];sR2[found-1]+=1-sX[j2+count2];}
else	//new
{sR1[found]=sR1[found-1]+sX[j2+count2];sR2[found]=sR2[found-1]+1-sX[j2+count2];found++;}
last=dptrs[j].value;
}
}

sum=0;for(j=1;j<found;j++){sum+=.5*(sR1[j]-sR1[j-1])*(sR2[j]+sR2[j-1]);}
stats[p+3*num_blocks]=sum/sR1[found-1]/sR2[found-1];
}
}	//end of p loop
printf("\n");

if(k<count){sprintf(readstring, "%d", k+1);}
else{strcpy(readstring,"Combined");}

sum=0;sumsq=0;
for(p=0;p<num_blocks;p++){sum+=stats[p];sumsq+=pow(stats[p],2);}
mean=sum/num_blocks;
var=(num_blocks-1)*(sumsq/num_blocks-pow(mean,2));
printf("Correlation is %.4f (SE %.4f)\n", value, pow(var,.5));
fprintf(output, "%s Correlation %.6f %.6f\n", readstring, value, pow(var,.5));

sum=0;sumsq=0;
for(p=0;p<num_blocks;p++){sum+=pow(stats[p],2);sumsq+=pow(stats[p],4);}
mean=sum/num_blocks;
var=(num_blocks-1)*(sumsq/num_blocks-pow(mean,2));
printf("Correlation squared is %.4f (SE %.4f)\n", pow(value,2), pow(var,.5));
fprintf(output, "%s Squared_correlation %.6f %.6f\n", readstring, pow(value,2), pow(var,.5));

sum=0;sumsq=0;
for(p=0;p<num_blocks;p++){sum+=stats[p+num_blocks];sumsq+=pow(stats[p+num_blocks],2);}
mean=sum/num_blocks;
var=(num_blocks-1)*(sumsq/num_blocks-pow(mean,2));
printf("Mean squared error is %.4f (SE %.4f)\n", value2/value5, pow(var,.5)/value5);
fprintf(output, "%s Mean_squared_error %.6f %.6f\n", readstring, value2/value5, pow(var,.5)/value5);

sum=0;sumsq=0;
for(p=0;p<num_blocks;p++){sum+=stats[p+2*num_blocks];sumsq+=pow(stats[p+2*num_blocks],2);}
mean=sum/num_blocks;
var=(num_blocks-1)*(sumsq/num_blocks-pow(mean,2));
printf("Mean absolute error is %.4f (SE %.4f)\n", value3*pow(value5,-.5), pow(var,.5)*pow(value5,-.5));
fprintf(output, "%s Mean_absolute_error %.6f %.6f\n", readstring, value3*pow(value5,-.5), pow(var,.5)*pow(value5,-.5));

if(dichot==1)	//now auc
{
sum=0;sumsq=0;
for(p=0;p<num_blocks;p++){sum+=stats[p+3*num_blocks];sumsq+=pow(stats[p+3*num_blocks],2);}
mean=sum/num_blocks;
var=(num_blocks-1)*(sumsq/num_blocks-pow(mean,2));
printf("Area under curve is %.4f (SE %.4f)\n", value4, pow(var,.5));
fprintf(output, "%s Area_under_curve %.6f %.6f\n", readstring, value4, pow(var,.5));
}

if(prev!=-9999)	//now liability correlation
{
sum=0;sumsq=0;
for(p=0;p<num_blocks;p++){sum+=pow(stats[p],2);sumsq+=pow(stats[p],4);}
mean=sum/num_blocks;
var=(num_blocks-1)*(sumsq/num_blocks-pow(mean,2));
printf("Correlation squared on the liability scale is %.4f (SE %.4f)\n", pow(value,2)*factor, pow(var,.5)*factor);
fprintf(output, "%s Liability_squared_correlation %.6f %.6f\n", readstring, pow(value,2)*factor, pow(var,.5)*factor);
}
printf("\n");
}
}

fclose(output);
printf("Results saved in %s\n\n", filename);

free(order);free(sX);free(sZ);free(sT);free(sXTX);free(sXTX2);free(stats);free(likes);
if(strcmp(proffile,"blank")!=0){for(k=0;k<count+cflag;k++){free(guesses[k]);}free(guesses);}
if(dichot==1){free(dptrs);free(sR1);free(sR2);}

///////////////////////////

