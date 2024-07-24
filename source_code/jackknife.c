/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
f
    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Jackknife to get variances

///////////////////////////

//get number of datapoints and set num_blocks
if(strcmp(jackfile,"blank")!=0)
{count=countcols(jackfile);count2=countrows(jackfile);}
else
{count=(countcols(proffile)-4)/2;count2=countrows(proffile)-1;}

if(count2<3){printf("Error, there are only %d datapoints, it is not possible to continue\n\n", count2);exit(1);}
if(count2<200){printf("Warning, there are usually thousands of datapoints (not %d)\n\n", count2);}

if(num_blocks==-1||num_blocks>count2){num_blocks=count2;}

//allocate variables
order=malloc(sizeof(int)*count2);
sX=malloc(sizeof(double)*count2*3);
sXTX=malloc(sizeof(double)*9);
sXTX2=malloc(sizeof(double)*9);
sW=malloc(sizeof(double)*count2);
sT=malloc(sizeof(double)*count2);
stats=malloc(sizeof(double)*num_blocks*4);

if(strcmp(proffile,"blank")!=0)	//using a profile file, so can be quite big
{
guesses=malloc(sizeof(double*)*count);
for(k=0;k<count;k++){guesses[k]=malloc(sizeof(double)*count2);}
}

if(auc==1)
{
dptrs=malloc(sizeof(struct sorting_double)*count2);
sZ1=malloc(sizeof(int)*(count2+1));
sZ2=malloc(sizeof(int)*(count2+1));
}

//open results file
sprintf(filename,"%s.jack",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s - check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
if(strcmp(jackfile,"blank")!=0){fprintf(output, "Measure Estimate SD\n");}
else{fprintf(output, "Profile Measure Estimate SD\n");}

////////

if(strcmp(jackfile,"blank")!=0)	//using a jack file
{
//load up predicted and observed values, and if present, weights
if(count==2){printf("Reading predicted and observed values for %d samples from %s\n", count2, jackfile);}
else{printf("Reading predicted and observed values for %d samples, as well as weightings, from %s\n", count2, jackfile);}
read_values(jackfile,sX,count2,NULL,1,0,0);
read_values(jackfile,sX+count2,count2,NULL,2,0,0);
for(j=0;j<count2;j++){sX[j+2*count2]=1;}

if(count==3)
{
read_values(jackfile,sW,count2,NULL,3,0,0);
for(j=0;j<count2;j++)
{
if(sW[j]<=0){printf("Error reading %s; weight %d is %f (all weights should be positive)\n\n", jackfile, j+1, sW[j]);exit(1);} 
}
}
else
{
for(j=0;j<count2;j++){sW[j]=1;}
}

if(permute==1)	//permute first two columns of sX and then sW (latter might be redundant)
{
for(j=0;j<count2;j++){order[j]=j;}
permute_int(order,count2);

for(k=0;k<2;k++)
{
for(j=0;j<count2;j++){sT[j]=sX[j+k*count2];}
for(j=0;j<count2;j++){sX[j+k*count2]=sT[order[j]];}
}
for(j=0;j<count2;j++){sT[j]=sW[j];}
for(j=0;j<count2;j++){sW[j]=sT[order[j]];}
}

if(auc==1||prev!=-9999)	//check observed values are zero or one (but not all same) and sort
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

//sort based on first column
for(j=0;j<count2;j++){dptrs[j].value=sX[j];dptrs[j].index=j;}
qsort(dptrs, count2, sizeof(struct sorting_double), compare_sorting_double);

if(prev!=-9999){factor=get_factor(sX+count2, count2, prev, -9999, outfile);}
}

printf("\n");

//first jackknife without weights
if(count==3){printf("First analyzing without weights\n");}

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

if(auc==1)	//now auc
{
sZ1[0]=0;sZ2[0]=0;found=1;
last=dptrs[0].value-1;
for(j=0;j<count2;j++)
{
j2=dptrs[j].index;
if(dptrs[j].value==last)	//same as previous value
{sZ1[found-1]+=sX[j2+count2];sZ2[found-1]+=1-sX[j2+count2];}
else	//new
{sZ1[found]=sZ1[found-1]+sX[j2+count2];sZ2[found]=sZ2[found-1]+1-sX[j2+count2];found++;}
last=dptrs[j].value;
}

sum=0;for(j=1;j<found;j++){sum+=.5*(sZ1[j]-sZ1[j-1])*(sZ2[j]+sZ2[j-1]);}
value4=sum/sZ1[found-1]/sZ2[found-1];
}	//end of auc

//now jackknife
for(p=0;p<num_blocks;p++)
{
start=(double)p/num_blocks*count2;
end=(double)(p+1)/num_blocks*count2;
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

if(auc==1)	//now auc
{
sZ1[0]=0;sZ2[0]=0;found=1;
last=dptrs[0].value-1;
for(j=0;j<count2;j++)
{
j2=dptrs[j].index;
if(j2<start||j2>=end)
{
if(dptrs[j].value==last)	//same as previous value
{sZ1[found-1]+=sX[j2+count2];sZ2[found-1]+=1-sX[j2+count2];}
else	//new
{sZ1[found]=sZ1[found-1]+sX[j2+count2];sZ2[found]=sZ2[found-1]+1-sX[j2+count2];found++;}
last=dptrs[j].value;
}
}

sum=0;for(j=1;j<found;j++){sum+=.5*(sZ1[j]-sZ1[j-1])*(sZ2[j]+sZ2[j-1]);}
stats[p+3*num_blocks]=sum/sZ1[found-1]/sZ2[found-1];
}
}	//end of p loop
printf("\n");

sum=0;sumsq=0;
for(p=0;p<num_blocks;p++){sum+=stats[p];sumsq+=pow(stats[p],2);}
mean=sum/num_blocks;
var=(num_blocks-1)*(sumsq/num_blocks-pow(mean,2));
printf("Correlation is %.4f (SD %.4f)\n", value, pow(var,.5));
fprintf(output, "Correlation %.6f %.6f\n", value, pow(var,.5));

sum=0;sumsq=0;
for(p=0;p<num_blocks;p++){sum+=pow(stats[p],2);sumsq+=pow(stats[p],4);}
mean=sum/num_blocks;
var=(num_blocks-1)*(sumsq/num_blocks-pow(mean,2));
printf("Correlation squared is %.4f (SD %.4f)\n", pow(value,2), pow(var,.5));
fprintf(output, "Squared_correlation %.6f %.6f\n", pow(value,2), pow(var,.5));

sum=0;sumsq=0;
for(p=0;p<num_blocks;p++){sum+=stats[p+num_blocks];sumsq+=pow(stats[p+num_blocks],2);}
mean=sum/num_blocks;
var=(num_blocks-1)*(sumsq/num_blocks-pow(mean,2));
printf("Mean squared error is %.4f (SD %.4f)\n", value2/value5, pow(var,.5)/value5);
fprintf(output, "Mean_squared_error %.6f %.6f\n", value2/value5, pow(var,.5)/value5);

sum=0;sumsq=0;
for(p=0;p<num_blocks;p++){sum+=stats[p+2*num_blocks];sumsq+=pow(stats[p+2*num_blocks],2);}
mean=sum/num_blocks;
var=(num_blocks-1)*(sumsq/num_blocks-pow(mean,2));
printf("Mean absolute error is %.4f (SD %.4f)\n", value3*pow(value5,-.5), pow(var,.5)*pow(value5,-.5));
fprintf(output, "Mean_absolute_error %.6f %.6f\n", value3*pow(value5,-.5), pow(var,.5)*pow(value5,-.5));

if(auc==1)	//now auc
{
sum=0;sumsq=0;
for(p=0;p<num_blocks;p++){sum+=stats[p+3*num_blocks];sumsq+=pow(stats[p+3*num_blocks],2);}
mean=sum/num_blocks;
var=(num_blocks-1)*(sumsq/num_blocks-pow(mean,2));
printf("Area under curve is %.4f (SD %.4f)\n", value4, pow(var,.5));
fprintf(output, "Area_under_curve %.6f %.6f\n", value4, pow(var,.5));
}

if(prev!=-9999)	//now liability correlation
{
sum=0;sumsq=0;
for(p=0;p<num_blocks;p++){sum+=pow(stats[p],2);sumsq+=pow(stats[p],4);}
mean=sum/num_blocks;
var=(num_blocks-1)*(sumsq/num_blocks-pow(mean,2));
printf("Correlation squared on the liability scale is %.4f (SD %.4f)\n", pow(value,2)*factor, pow(var,.5)*factor);
fprintf(output, "Liability_squared_correlation %.6f %.6f\n", pow(value,2)*factor, pow(var,.5)*factor);
}
printf("\n");

if(count==3)	//jackknife with weights
{
printf("Now analyzing with weights\n");

//get weighted correlation across all datapoints
for(j=0;j<count2;j++){sX[j]*=pow(sW[j],.5);sX[j+count2]*=pow(sW[j],.5);sX[j+2*count2]=pow(sW[j],.5);}
alpha=1.0;beta=0.0;
dgemm_("T", "N", &three, &three, &count2, &alpha, sX, &count2, sX, &count2, &beta, sXTX, &three);
sum=sXTX[2];sum2=sXTX[5];
sumsq=sXTX[0];sumsq2=sXTX[4];sumsq3=sXTX[1];
value=(sXTX[8]*sumsq3-sum*sum2)*pow(sXTX[8]*sumsq-sum*sum,-.5)*pow(sXTX[8]*sumsq2-sum2*sum2,-.5);

//store weighted variance of observed values in value5 - will use this to scale mse and mae
value5=sumsq2/sXTX[8]-pow(sum2/sXTX[8],2);

//now weighted mse and mae
sum=0;for(j=0;j<count2;j++){sum+=pow(sX[j]-sX[j+count2],2);}
value2=sum/sXTX[8];
sum=0;for(j=0;j<count2;j++){sum+=pow(sW[j],.5)*fabs(sX[j]-sX[j+count2]);}
value3=sum/sXTX[8];

for(p=0;p<num_blocks;p++)
{
start=(double)p/num_blocks*count2;
end=(double)(p+1)/num_blocks*count2;
if(p%500000==0){printf("Performing Jackknife %d out of %d\n", p+1, num_blocks);}
count3=end-start;

//first correlation
alpha=1.0;beta=0.0;
dgemm_("T", "N", &three, &three, &count3, &alpha, sX+start, &count2, sX+start, &count2, &beta, sXTX2, &three);
sum=sXTX[2]-sXTX2[2];sum2=sXTX[5]-sXTX2[5];value4=sXTX[8]-sXTX2[8];
sumsq=sXTX[0]-sXTX2[0];sumsq2=sXTX[4]-sXTX2[4];sumsq3=sXTX[1]-sXTX2[1];
stats[p]=(value4*sumsq3-sum*sum2)*pow(value4*sumsq-sum*sum,-.5)*pow(value4*sumsq2-sum2*sum2,-.5);

//now mse and mae
sum=0;for(j=start;j<end;j++){sum+=pow(sX[j]-sX[j+count2],2);}
stats[p+num_blocks]=(value2*sXTX[8]-sum)/value4;
sum=0;for(j=start;j<end;j++){sum+=pow(sW[j],.5)*fabs(sX[j]-sX[j+count2]);}
stats[p+2*num_blocks]=(value3*sXTX[8]-sum)/value4;
}

sum=0;sumsq=0;
for(p=0;p<num_blocks;p++){sum+=stats[p];sumsq+=pow(stats[p],2);}
mean=sum/num_blocks;
var=(num_blocks-1)*(sumsq/num_blocks-pow(mean,2));
printf("\nWeighted correlation is %.4f (SD %.4f)\n", value, pow(var,.5));
fprintf(output, "Weighted_correlation %.6f %.6f\n", value, pow(var,.5));

sum=0;sumsq=0;
for(p=0;p<num_blocks;p++){sum+=pow(stats[p],2);sumsq+=pow(stats[p],4);}
mean=sum/num_blocks;
var=(num_blocks-1)*(sumsq/num_blocks-pow(mean,2));
printf("Weighted correlation squared is %.4f (SD %.4f)\n", pow(value,2), pow(var,.5));
fprintf(output, "Weighted_squared_correlation %.6f %.6f\n", pow(value,2), pow(var,.5));

sum=0;sumsq=0;
for(p=0;p<num_blocks;p++){sum+=stats[p+num_blocks];sumsq+=pow(stats[p+num_blocks],2);}
mean=sum/num_blocks;
var=(num_blocks-1)*(sumsq/num_blocks-pow(mean,2));
printf("Weighted mean squared error is %.4f (SD %.4f)\n", value2/value5, pow(var,.5)/value5);
fprintf(output, "Weighted_mean_squared_error %.6f %.6f\n", value2/value5, pow(var,.5)/value5);

sum=0;sumsq=0;
for(p=0;p<num_blocks;p++){sum+=stats[p+2*num_blocks];sumsq+=pow(stats[p+2*num_blocks],2);}
mean=sum/num_blocks;
var=(num_blocks-1)*(sumsq/num_blocks-pow(mean,2));
printf("Weighted mean absolute error is %.4f (SD %.4f)\n", value3*pow(value5,-.5), pow(var,.5)*pow(value5,-.5));
fprintf(output, "Weighted_mean_absolute_error %.6f %.6f\n", value3*pow(value5,-.5), pow(var,.5)*pow(value5,-.5));
printf("\n");
}	//end of count==3
}
else	//using profile file - might have missing phenotypes
{
//read observed values (minus covariates) into second column of sX, and predictions into guesses
guesses=malloc(sizeof(double*)*count);
for(k=0;k<count;k++){guesses[k]=malloc(sizeof(double)*count2);}

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

//ready to analyse - remember, we only have count4 datapoints (but sX has size count2)

//set third column of sX to one
for(j=0;j<count4;j++){sX[j+2*count2]=1;}

if(permute==1)	//permute second columns of sX and then all columns of guesses
{
for(j=0;j<count4;j++){order[j]=j;}
permute_int(order,count4);

for(j=0;j<count4;j++){sT[j]=sX[j+count2];}
for(j=0;j<count4;j++){sX[j+count2]=sT[order[j]];}

for(k=0;k<count;k++)
{
for(j=0;j<count4;j++){sT[j]=guesses[k][j];}
for(j=0;j<count4;j++){guesses[k][j]=sT[order[j]];}
}
}

if(auc==1||prev!=-9999)	//check observed values are zero or one (but not all same)
{
total=0;
for(j=0;j<count4;j++)
{
if(sX[j+count2]!=0&&sX[j+count2]!=1)
{printf("Error reading %s; phenotypes should be either one or zero (corresponding to cases and controls)\n\n", proffile);exit(1);}
total+=sX[j+count2];
}
if(total==0){printf("Error reading %s; all phenotypes are zero\n\n", proffile);exit(1);}
if(total==count4){printf("Error reading %s; all phenotypes are one\n\n", proffile);exit(1);}

if(prev!=-9999){factor=get_factor(sX+count2, count2, prev, -9999, outfile);}
}

for(k=0;k<count;k++)
{
printf("Analysing Phenotype %d\n", k+1);
//load predicted values for Profile k into first column of sX
for(j=0;j<count4;j++){sX[j]=guesses[k][j];}

//get correlation across all datapoints
alpha=1.0;beta=0.0;
dgemm_("T", "N", &three, &three, &count4, &alpha, sX, &count2, sX, &count2, &beta, sXTX, &three);
sum=sXTX[2];sum2=sXTX[5];
sumsq=sXTX[0];sumsq2=sXTX[4];sumsq3=sXTX[1];
value=(count4*sumsq3-sum*sum2)*pow(count4*sumsq-sum*sum,-.5)*pow(count4*sumsq2-sum2*sum2,-.5);

//store variance of observed values in value5 - will use this to scale mse and mae
value5=sumsq2/count2-pow(sum2/count2,2);

//now mse and mae
sum=0;for(j=0;j<count4;j++){sum+=pow(sX[j]-sX[j+count2],2);}
value2=sum/count4;
sum=0;for(j=0;j<count4;j++){sum+=fabs(sX[j]-sX[j+count2]);}
value3=sum/count4;

if(auc==1)	//now auc (must first sort)
{
for(j=0;j<count4;j++){dptrs[j].value=sX[j];dptrs[j].index=j;}
qsort(dptrs, count4, sizeof(struct sorting_double), compare_sorting_double);

sZ1[0]=0;sZ2[0]=0;found=1;
last=dptrs[0].value-1;
for(j=0;j<count4;j++)
{
j2=dptrs[j].index;
if(dptrs[j].value==last)	//same as previous value
{sZ1[found-1]+=sX[j2+count2];sZ2[found-1]+=1-sX[j2+count2];}
else	//new
{sZ1[found]=sZ1[found-1]+sX[j2+count2];sZ2[found]=sZ2[found-1]+1-sX[j2+count2];found++;}
last=dptrs[j].value;
}

sum=0;for(j=1;j<found;j++){sum+=.5*(sZ1[j]-sZ1[j-1])*(sZ2[j]+sZ2[j-1]);}
value4=sum/sZ1[found-1]/sZ2[found-1];
}	//end of auc

//now jackknife

for(p=0;p<num_blocks;p++)
{
start=(double)p/num_blocks*count4;
end=(double)(p+1)/num_blocks*count4;
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

if(auc==1)	//now auc
{
sZ1[0]=0;sZ2[0]=0;found=1;
last=dptrs[0].value-1;
for(j=0;j<count4;j++)
{
j2=dptrs[j].index;
if(j2<start||j2>=end)
{
if(dptrs[j].value==last)	//same as previous value
{sZ1[found-1]+=sX[j2+count2];sZ2[found-1]+=1-sX[j2+count2];}
else	//new
{sZ1[found]=sZ1[found-1]+sX[j2+count2];sZ2[found]=sZ2[found-1]+1-sX[j2+count2];found++;}
last=dptrs[j].value;
}
}

sum=0;for(j=1;j<found;j++){sum+=.5*(sZ1[j]-sZ1[j-1])*(sZ2[j]+sZ2[j-1]);}
stats[p+3*num_blocks]=sum/sZ1[found-1]/sZ2[found-1];
}
}	//end of p loop
printf("\n");

sum=0;sumsq=0;
for(p=0;p<num_blocks;p++){sum+=stats[p];sumsq+=pow(stats[p],2);}
mean=sum/num_blocks;
var=(num_blocks-1)*(sumsq/num_blocks-pow(mean,2));
printf("Correlation is %.4f (SD %.4f)\n", value, pow(var,.5));
fprintf(output, "%d Correlation %.6f %.6f\n", k+1, value, pow(var,.5));

sum=0;sumsq=0;
for(p=0;p<num_blocks;p++){sum+=pow(stats[p],2);sumsq+=pow(stats[p],4);}
mean=sum/num_blocks;
var=(num_blocks-1)*(sumsq/num_blocks-pow(mean,2));
printf("Correlation squared is %.4f (SD %.4f)\n", pow(value,2), pow(var,.5));
fprintf(output, "%d Squared_correlation %.6f %.6f\n", k+1, pow(value,2), pow(var,.5));

sum=0;sumsq=0;
for(p=0;p<num_blocks;p++){sum+=stats[p+num_blocks];sumsq+=pow(stats[p+num_blocks],2);}
mean=sum/num_blocks;
var=(num_blocks-1)*(sumsq/num_blocks-pow(mean,2));
printf("Mean squared error is %.4f (SD %.4f)\n", value2/value5, pow(var,.5)/value5);
fprintf(output, "%d Mean_squared_error %.6f %.6f\n", k+1, value2/value5, pow(var,.5)/value5);

sum=0;sumsq=0;
for(p=0;p<num_blocks;p++){sum+=stats[p+2*num_blocks];sumsq+=pow(stats[p+2*num_blocks],2);}
mean=sum/num_blocks;
var=(num_blocks-1)*(sumsq/num_blocks-pow(mean,2));
printf("Mean absolute error is %.4f (SD %.4f)\n", value3*pow(value5,-.5), pow(var,.5)*pow(value5,-.5));
fprintf(output, "%d Mean_absolute_error %.6f %.6f\n", k+1, value3*pow(value5,-.5), pow(var,.5)*pow(value5,-.5));

if(auc==1)	//now auc
{
sum=0;sumsq=0;
for(p=0;p<num_blocks;p++){sum+=stats[p+3*num_blocks];sumsq+=pow(stats[p+3*num_blocks],2);}
mean=sum/num_blocks;
var=(num_blocks-1)*(sumsq/num_blocks-pow(mean,2));
printf("Area under curve is %.4f (SD %.4f)\n", value4, pow(var,.5));
fprintf(output, "%d Area_under_curve %.6f %.6f\n", k+1, value4, pow(var,.5));
}

if(prev!=-9999)	//now liability correlation
{
sum=0;sumsq=0;
for(p=0;p<num_blocks;p++){sum+=pow(stats[p],2);sumsq+=pow(stats[p],4);}
mean=sum/num_blocks;
var=(num_blocks-1)*(sumsq/num_blocks-pow(mean,2));
printf("Correlation squared on the liability scale is %.4f (SD %.4f)\n", pow(value,2)*factor, pow(var,.5)*factor);
fprintf(output, "%d Liability_squared_correlation %.6f %.6f\n", k+1, pow(value,2)*factor, pow(var,.5)*factor);
}
printf("\n");
}
}

fclose(output);
printf("Results saved in %s\n\n", filename);

free(order);free(sX);free(sXTX);free(sXTX2);free(sW);free(sT);free(stats);
if(strcmp(proffile,"blank")!=0){for(k=0;k<count;k++){free(guesses[k]);}free(guesses);}
if(auc==1){free(dptrs);free(sZ1);free(sZ2);}

///////////////////////////

