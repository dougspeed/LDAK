/*
Copyright 2026 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//deal with extra traits - non-target traits can have missing summary statistics

///////////////////////////

sigma=malloc(sizeof(double)*num_sums3*num_sums3);
omega=malloc(sizeof(double)*num_sums3*num_sums3);

rjksums2=malloc(sizeof(double)*data_length*num_cors);
rjksums3=malloc(sizeof(double)*data_length*num_cors);
rjksums4=malloc(sizeof(double)*data_length*num_cors);
comat=malloc(sizeof(double)*num_sums3*num_sums3);
comat2=malloc(sizeof(double)*num_sums3);
comat3=malloc(sizeof(double)*num_sums3);
comat4=malloc(sizeof(double)*num_sums3);

//rjksums2 contains power -1 taggings for all correlations (can be missing values)
//rjksums3 contains best power taggings for all correlations (can be missing values)
for(s=0;s<num_cors;s++)
{
sprintf(filename2,"%s.cors.bin", corstems[s]);
if((input2=fopen(filename2,"rb"))==NULL)
{printf("Error re-opening %s\n\n",filename2);exit(1);}

for(j=0;j<data_length;j++){rjksums2[j+s*data_length]=0;}
for(j=0;j<Dnp2[s];j++)
{
fseeko(input2, sizeof(double)*Dnp[s]*3+sizeof(double)*Dkp[s][j], SEEK_SET);
if(fread(rjksums2+Dkp2[s][j]+s*data_length, sizeof(double), 1, input2)!=1)
{printf("Error reading Tagging %d for Predictor %d from %s\n\n", k+1, j+1, filename2);exit(1);}
}

for(j=0;j<data_length;j++){rjksums3[j+s*data_length]=0;}
for(j=0;j<Dnp2[s];j++)
{
fseeko(input2, sizeof(double)*Dnp[s]*(3+best)+sizeof(double)*Dkp[s][j], SEEK_SET);
if(fread(rjksums3+Dkp2[s][j]+s*data_length, sizeof(double), 1, input2)!=1)
{printf("Error reading Tagging %d for Predictor %d from %s\n\n", k+1, j+1, filename2);exit(1);}
}

fclose(input2); 
}

//rjksums4 contains cross-ancestry taggings (can ingore target populations and those not being used)
indexer=malloc(sizeof(int)*num_cors);

//set all values to one
for(s=0;s<num_cors;s++)
{
for(j=0;j<data_length;j++){rjksums4[j+s*data_length]=1;}
}

//are we using any non-target correlations
for(s=0;s<num_cors;s++){indexer[s]=0;}
for(q=1;q<num_sums3;q++)
{
if(sumpops[q]!=sumpops[0]){indexer[sumpops[q]]=1;}
}
count=0;for(s=0;s<num_cors;s++){count+=indexer[s];}

if(count>0)
{
printf("Calculating cross-ancestry correlations for %d ancestries\n", count);

for(bit=0;bit<bittotal;bit++)
{
bitstart=blockstarts[bit];
bitend=blockends[bit];
bitlength=blockends[bit]-blockstarts[bit];

if(bit%50==0)
{
printf("Calculating cross-ancestry taggings for Window %d of %d\n", bit+1, bittotal);

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Calculating cross-ancestry taggings for Window %d of %d\n", bit+1, bittotal);
fclose(output);
}

//re-open target correlations
s=sumpops[0];
sprintf(filename2,"%s.cors.bin", corstems[s]);
if((input2=fopen(filename2,"rb"))==NULL)
{printf("Error re-opening %s\n\n",filename2);exit(1);}

//read lower triangle of target correlations
fseeko(input2, Dindexes[s][bit], SEEK_SET);
for(k=0;k<Dsizes[s][bit][0]-1;k++)
{
count2=fread(cors_short+(size_t)k*Dsizes[s][bit][0]+k+1, sizeof(unsigned short), Dsizes[s][bit][0]-k-1, input2);
if(count2!=Dsizes[s][bit][0]-k-1)
{printf("Error reading correlations for Window %d from %s\n\n", bit+1, filename2);exit(1);}
}
fclose(input2);

//extract correlations (no need for signs
#pragma omp parallel for private(k,j) schedule(static)
for(k=0;k<bitlength;k++)
{
for(j=k+1;j<bitlength;j++)
{cors[(size_t)k*bitlength+j]=0.00005*cors_short[(size_t)Duse[s][bit][k]*Dsizes[s][bit][0]+Duse[s][bit][j]]-1;}
}

for(s=0;s<num_cors;s++)
{
if(indexer[s]==1)   //using this population
{
sprintf(filename2,"%s.cors.bin", corstems[s]);
if((input2=fopen(filename2,"rb"))==NULL)
{printf("Error re-opening %s\n\n",filename2);exit(1);}

//read lower triangle of non-target correlations
fseeko(input2, Dindexes[s][bit], SEEK_SET);
for(k=0;k<Dsizes[s][bit][0]-1;k++)
{
count2=fread(cors_short+(size_t)k*Dsizes[s][bit][0]+k+1, sizeof(unsigned short), Dsizes[s][bit][0]-k-1, input2);
if(count2!=Dsizes[s][bit][0]-k-1)
{printf("Error reading correlations for Window %d from %s\n\n", bit+1, filename2);exit(1);}
}
fclose(input2);

//add on sums (remember signs)
#pragma omp parallel for private(k,j,value) schedule(static)
for(k=0;k<Dsizes[s][bit][1];k++)
{
for(j=0;j<k;j++)
{
value=(0.00005*cors_short[(size_t)Duse[s][bit][j]*Dsizes[s][bit][0]+Duse[s][bit][k]]-1)*Dsigns[s][bit][j]*Dsigns[s][bit][k];
rjksums4[bitstart+Duse2[s][bit][k]+s*data_length]+=cors[(size_t)Duse2[s][bit][j]*bitlength+Duse2[s][bit][k]]*value;
}
for(j=k+1;j<Dsizes[s][bit][1];j++)
{
value=(0.00005*cors_short[(size_t)Duse[s][bit][k]*Dsizes[s][bit][0]+Duse[s][bit][j]]-1)*Dsigns[s][bit][j]*Dsigns[s][bit][k];
rjksums4[bitstart+Duse2[s][bit][k]+s*data_length]+=cors[(size_t)Duse2[s][bit][k]*bitlength+Duse2[s][bit][j]]*value;
}
}

}}   //end of using this population and s loop
}   //end of bit loop

printf("\n");
}   //end of count>0

free(indexer);

////////

printf("Computing omega and sigma, the MTAG covariance matrices\n");

//start with diagonals
#pragma omp parallel for private(q,j,sum,sumsq) schedule(static)
for(q=0;q<num_sums3;q++)
{
//get same-trait intercept
if(mtagcopy==0) //set to one
{sigma[q+q*num_sums3]=1;}
else    //estimate its value
{
sigma[q+q*num_sums3]=solve_intercept_wrapper1(rjksums2+sumpops[q]*data_length, Mnss[q], Mchis[q], data_length);
if(sigma[q+q*num_sums3]<0.5)
{
printf("Warning, trait %d has a very low estimated intercept (%.4f); has been increased to 0.5", q+1, sigma[q+q*num_sums3]);
sigma[q+q*num_sums3]=0.5;
}
}

//get same-trait correlation, ensuring values not too small
sum=0;sumsq=0;
for(j=0;j<data_length;j++)
{
if(Mnss[q][j]>0)	//have summaries
{
if(mtagcopy==0){sumsq+=(pow(Mrhos[q][j],2)-sigma[q+q*num_sums3]/Mnss[q][j])/rjksums[j];sum+=pow(rjksums[j],-1);}
else{sumsq+=pow(Mrhos[q][j],2)-sigma[q+q*num_sums3]/Mnss[q][j];sum++;}
}
}
if(sumsq<0.01){sumsq=0.01;}
omega[q+q*num_sums3]=sumsq/sum;
}

//now target vs other traits
#pragma omp parallel for private(q,j,sum,sumsq) schedule(static)
for(q=1;q<num_sums3;q++)
{
if(sumpops[0]!=sumpops[q]||pritraits[q]==1) //set intercept to zero
{
sigma[q]=0;
sigma[q*num_sums3]=0;
}
else    //estimate its value
{
if(mtagcopy==0){sigma[q]=solve_intercept_wrapper2(rjksums2+sumpops[0]*data_length, rjksums3+sumpops[0]*data_length, Mnss[0], Mnss[q], Mchis[0], Mchis[q], Mrhos[0], Mrhos[q], data_length, q);}
else{sigma[q]=solve_intercept_wrapper2(rjksums2+sumpops[0]*data_length, rjksums2+sumpops[0]*data_length, Mnss[0], Mnss[q], Mchis[0], Mchis[q], Mrhos[0], Mrhos[q], data_length, q);}
sigma[q*num_sums3]=sigma[q];
}

sumsq=0;sum=0;
for(j=0;j<data_length;j++)
{
if(Mnss[0][j]>0&&Mnss[q][j]>0)	//have both sets of summaries
{
if(mtagcopy==0){sumsq+=(Mrhos[0][j]*Mrhos[q][j]-sigma[q]*pow(Mnss[0][j],-.5)*pow(Mnss[q][j],-.5))/rjksums[j];sum+=pow(rjksums[j],-1);}
else{sumsq+=Mrhos[0][j]*Mrhos[q][j]-sigma[q]*pow(Mnss[0][j],-.5)*pow(Mnss[q][j],-.5);sum++;}
}
}
omega[q]=sumsq/sum;
omega[q*num_sums3]=sumsq/sum;
}

if(mtagforce==0)    //blank secondary traits with low correlation
{
count=0;for(q=0;q<num_sums3;q++){count+=(pritraits[q]==0);}

count2=0;
for(q=1;q<num_sums3;q++)
{
value=omega[q]*pow(omega[0]*omega[q+q*num_sums3],-.5);
if(pritraits[q]==0&&fabs(value)<mtagcor)
{
for(j=0;j<data_length;j++){Mnss[q][j]=0;Mchis[q][j]=0;Mrhos[q][j]=0;}
pritraits[q]=2;
count2++;
}
}

if(count2>0)
{
if(count2<count){printf("Only %d of the %d secondary summary statistics have correlation above %.4f (the remainder will be ignored)\n\n", count-count2, count, mtagcor);}
else{printf("None of the %d secondary summary statistics have correlation above %.4f (so all will be ignored)\n\n", count, mtagcor);}
}
}

//get remaining off-diagonal values (if required)
for(q=1;q<num_sums3;q++)
{
//set values to zero, then update for non-excluded traits
for(q2=q+1;q2<num_sums3;q2++){sigma[q+q2*num_sums3]=0;sigma[q2+q*num_sums3]=0;omega[q+q2*num_sums3]=0;omega[q2+q*num_sums3]=0;}

if(pritraits[q]!=2)
{
#pragma omp parallel for private(q2,j,sum,sumsq) schedule(dynamic)
for(q2=q+1;q2<num_sums3;q2++)
{
if(pritraits[q2]!=2)
{
//get cross-trait intercept
if(sumpops[q]!=sumpops[q2]||(pritraits[q]==1&&pritraits[q2]==1)) //set to zero
{
sigma[q+q2*num_sums3]=0;
sigma[q2+q*num_sums3]=0;
}
else    //estimate its value
{
if(mtagcopy==0){sigma[q+q2*num_sums3]=solve_intercept_wrapper2(rjksums2+sumpops[q]*data_length, rjksums3+sumpops[q]*data_length, Mnss[q], Mnss[q2], Mchis[q], Mchis[q2], Mrhos[q], Mrhos[q2], data_length, q);}
else{sigma[q+q2*num_sums3]=solve_intercept_wrapper2(rjksums2+sumpops[q]*data_length, rjksums2+sumpops[q]*data_length, Mnss[q], Mnss[q2], Mchis[q], Mchis[q2], Mrhos[q], Mrhos[q2], data_length, q);}
sigma[q2+q*num_sums3]=sigma[q+q2*num_sums3];
}

//get cross-trait correlations
sum=0;sumsq=0;
for(j=0;j<data_length;j++)
{
if(Mnss[q][j]>0&&Mnss[q2][j]>0)	//have both sets of summaries
{
if(mtagcopy==0){sumsq+=(Mrhos[q][j]*Mrhos[q2][j]-sigma[q+q2*num_sums3]*pow(Mnss[q][j],-.5)*pow(Mnss[q2][j],-.5))/rjksums[j];sum+=pow(rjksums[j],-1);}
else{sumsq+=Mrhos[q][j]*Mrhos[q2][j]-sigma[q+q2*num_sums3]*pow(Mnss[q][j],-.5)*pow(Mnss[q2][j],-.5);sum++;}
}
}
omega[q+q2*num_sums3]=sumsq/sum;
omega[q2+q*num_sums3]=sumsq/sum;
}}
}
}

////////

//comat contains average covariance matrix for all (non-trivial) traits, comat2 contains omega_1 / omega_11 (or zero if trait trivial)
for(q=0;q<num_sums3;q++)
{
sum=0;for(j=0;j<data_length;j++){sum+=Mnss[q][j];}
mean=sum/data_length;

if(mean>0)
{
comat[q+q*num_sums3]=omega[q+q*num_sums3]-omega[q]*omega[q]/omega[0]+sigma[q+q*num_sums3]/mean;
comat2[q]=omega[q]/omega[0];
}
else{comat[q+q*num_sums3]=1;comat2[q]=0;}

for(q2=q+1;q2<num_sums3;q2++)
{
sum2=0;for(j=0;j<data_length;j++){sum2+=Mnss[q2][j];}
mean2=sum2/data_length;
if(mean>0&&mean2>0)
{
comat[q+q2*num_sums3]=omega[q+q2*num_sums3]-omega[q]*omega[q2]/omega[0]+sigma[q+q2*num_sums3]*pow(mean,-.5)*pow(mean2,-.5);
comat[q2+q*num_sums3]=comat[q+q2*num_sums3];
}
else{comat[q+q2*num_sums3]=0;comat[q2+q*num_sums3]=0;}
}
}

if(mtagforce==1)  //will use all (non-trivial) traits
{mtag_select(comat, comat2, comat3, comat4, num_sums3, pritraits, -9999);}
else    //use target trait, primary traits and selected secondary traits
{mtag_select(comat, comat2, comat3, comat4, num_sums3, pritraits, mtagpercent/100);}

printf("Here are the estimated genetic correlations, intercepts, importances and weightings\n");
for(q=0;q<num_sums3;q++)
{printf("Trait %d (%s): correlation %.4f, intercept %.4f, importance %.4f and weighting %.4f\n", q+1, sumstems[q], omega[q]*pow(omega[0]*omega[q+q*num_sums3],-.5), sigma[q], comat4[q]/comat4[0], comat3[q]);}
printf("\n");




sprintf(filename3,"%s.mtag.weights",outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3,"Summary_Statistics File Correlation Effective_Sample_Size Importance Weighting\n");
for(q=0;q<num_sums3;q++){fprintf(output3,"%d %s %.4f %.2f %.4f %.4f\n", q+1, sumstems[q], omega[q]*pow(omega[0]*omega[q+q*num_sums3],-.5), comat4[q], comat4[q]/comat4[0], comat3[q]);}
fclose(output3);

//how many primary and secondary traits are being used?
count=0;for(q=1;q<num_sums3;q++){count+=(pritraits[q]==1);}
count2=0;for(q=1;q<num_sums3;q++){count2+=(pritraits[q]==0&&comat3[q]!=0);}

sprintf(filename3,"%s.mtag.summary",outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3,"Using the target summary statistics\nall %d primary summary statistics\nand %d secondary summary statistics\n", count, count2);
fclose(output3);

//now use pritrait to indicate which traits being used
for(q=0;q<num_sums3;q++){pritraits[q]=(comat3[q]!=0);}

////////

//set summary statistics to zero for extra traits with weight zero
for(q=1;q<num_sums3;q++)
{
if(pritraits[q]==0)
{
for(j=0;j<data_length;j++){Mnss[q][j]=0;Mchis[q][j]=0;Mrhos[q][j]=0;}
}
}

if(impsums==1)  //impute summary statistic for extra trait with non-zero weight (might not have correlations for all predictors)
{
time(&starttime);

indexer=malloc(sizeof(int)*bitmax);
indexer2=malloc(sizeof(int)*bitmax);
cors3=malloc(sizeof(double)*bitmax);

for(q=1;q<num_sums3;q++)
{
if(pritraits[q]!=0)
{
count3=0;for(j=0;j<data_length;j++){count3+=(Mnss[q][j]==0&&rjksums2[j+sumpops[q]*data_length]>0);}
if(count3>0)
{
printf("Imputing summary statistics for %d predictors for Trait %d (%s)\n", count3, q+1, sumstems[q]);

//re-open corresponding correlations
s=sumpops[q];

sprintf(filename2,"%s.cors.bin", corstems[s]);
if((input2=fopen(filename2,"rb"))==NULL)
{printf("Error re-opening %s\n\n",filename2);exit(1);}

for(bit=0;bit<bittotal;bit++)
{
bitstart=blockstarts[bit];
bitend=blockends[bit];
bitlength=blockends[bit]-blockstarts[bit];

if(bit%50==0)
{
printf("Imputing for Window %d of %d\n",bit+1, bittotal);

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Trait %d: imputing summary statistics for Window %d of %d\n", q+1, bit+1, bittotal);
fclose(output);
}

//get numbers and indexes of predictors with and without summary statistics (indexes are relative to predictors with correlations)
count=0;
for(j=0;j<Dsizes[s][bit][1];j++)
{
if(Mnss[q][bitstart+Duse2[s][bit][j]]>0){indexer[count]=j;count++;}
else{indexer2[j-count]=j;}
}

if(count<Dsizes[s][bit][1]&&count>Dsizes[s][bit][1]/2) //some missing in this window (but not too many)
{
//read lower triangle of target correlations
fseeko(input2, Dindexes[s][bit], SEEK_SET);
for(k=0;k<Dsizes[s][bit][0]-1;k++)
{
count2=fread(cors_short+(size_t)k*Dsizes[s][bit][0]+k+1, sizeof(unsigned short), Dsizes[s][bit][0]-k-1, input2);
if(count2!=Dsizes[s][bit][0]-k-1)
{printf("Error reading correlations for Window %d from %s\n\n", bit+1, filename2);exit(1);}
}

//extract correlations for predictors with summary statistics (remember signs and shrink)
#pragma omp parallel for private(k,k2,j,j2) schedule(static)
for(k=0;k<count;k++)
{
k2=indexer[k];
cors[(size_t)k*count+k]=1.0;
for(j=k+1;j<count;j++)
{
j2=indexer[j];
cors[(size_t)k*count+j]=(0.00005*cors_short[(size_t)Duse[s][bit][k2]*Dsizes[s][bit][0]+Duse[s][bit][j2]]-1)*Dsigns[s][bit][j2]*Dsigns[s][bit][k2]*shrink;
}
}

//put Z scores for present predictors in cors3, and get mean sample size
sum=0;
for(j=0;j<count;j++)
{
j2=indexer[j];
if(Mrhos[q][bitstart+Duse2[s][bit][j2]]>0){cors3[j]=pow(Mchis[q][bitstart+Duse2[s][bit][j2]],.5);}
else{cors3[j]=-pow(Mchis[q][bitstart+Duse2[s][bit][j2]],.5);}
sum+=Mnss[q][bitstart+Duse2[s][bit][j2]];
}
mean=sum/count;

//get cholesky
dpotrf_("L", &count, cors, &count, &info);
if(info!=0){printf("Error, unable to decompose correlations for %d predictors in Window %d\n\n", count, bit+1);exit(1);}

//compute inv(cors) cors3
dpotrs_("L", &count, &one, cors, &count, cors3, &count, &info);
if(info!=0){printf("Error, unable to invert correlations for %d predictors in Window %d\n\n", count, bit+1);exit(1);}

//imputed Z scores are C_j cors3 (remember signs and shrink)
#pragma omp parallel for private(j,j2,k,k2,sum) schedule(static)
for(j=0;j<Dsizes[s][bit][1]-count;j++)
{
j2=indexer2[j];
sum=0;
for(k=0;k<count;k++)
{
k2=indexer[k];
if(j2>k2){sum+=(0.00005*cors_short[(size_t)Duse[s][bit][k2]*Dsizes[s][bit][0]+Duse[s][bit][j2]]-1)*Dsigns[s][bit][j2]*Dsigns[s][bit][k2]*shrink*cors3[k];}
else{sum+=(0.00005*cors_short[(size_t)Duse[s][bit][j2]*Dsizes[s][bit][0]+Duse[s][bit][k2]]-1)*Dsigns[s][bit][j2]*Dsigns[s][bit][k2]*shrink*cors3[k];}
}

Mnss[q][bitstart+Duse2[s][bit][j2]]=mean;
Mchis[q][bitstart+Duse2[s][bit][j2]]=pow(sum,2);
if(sum>0){Mrhos[q][bitstart+Duse2[s][bit][j2]]=pow(Mchis[q][bitstart+Duse2[s][bit][j2]]/(Mchis[q][bitstart+Duse2[s][bit][j2]]+mean),.5);}
else{Mrhos[q][bitstart+Duse2[s][bit][j2]]=-pow(Mchis[q][bitstart+Duse2[s][bit][j2]]/(Mchis[q][bitstart+Duse2[s][bit][j2]]+mean),.5);}
}
}
}   //end of bit loop
printf("\n");

fclose(input2);
}}}   //end of q loop and count3>0 and pritraits!=0

time(&endtime);
printf("Time check: have so far spent %.2f hours\n\n", (double)(endtime-starttime)/60/60);

free(indexer);free(indexer2);
}

////////

if(mtagreg==1)
{
sprintf(filename3,"%s.mtag.regions",outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
for(q=0;q<num_sums3;q++)
{
if(pritraits[q]!=0){fprintf(output3,"%.4f ", omega[q]*pow(omega[0]*omega[q+q*num_sums3],-.5));}
}
fprintf(output3,"%d\n", data_length);
}

//ready to combine summary statistics for target trait
printf("Computing the MTAG summary statistics\n\n");

for(bit=0;bit<bittotal;bit++)
{
bitstart=blockstarts[bit];
bitend=blockends[bit];
bitlength=blockends[bit]-blockstarts[bit];

if(mtagreg==1)  //recompute omega (for contributing traits) - must have mtagcopy=0
{
#pragma omp parallel for private(q,j,sum,sumsq) schedule(dynamic)
for(q=0;q<num_sums3;q++)
{
if(pritraits[q]==1)
{
sum=0;sumsq=0;
for(j=bitstart;j<bitend;j++)
{
if(Mnss[q][j]>0){sumsq+=(pow(Mrhos[q][j],2)-sigma[q+q*num_sums3]/Mnss[q][j])/rjksums[j];sum+=pow(rjksums[j],-1);}
}
if(sumsq<0.000001){sumsq=0.000001;}
omega[q+q*num_sums3]=sumsq/sum;
}}

#pragma omp parallel for private(q,q2,j,sum,sumsq,value) schedule(dynamic)
for(q=0;q<num_sums3;q++)
{
if(pritraits[q]==1)
{
for(q2=q+1;q2<num_sums3;q2++)
{
if(pritraits[q2]==1)
{
sumsq=0;sum=0;
for(j=bitstart;j<bitend;j++)
{
if(Mnss[q][j]>0&&Mnss[q2][j]>0){sumsq+=(Mrhos[q][j]*Mrhos[q2][j]-sigma[q+q2*num_sums3]*pow(Mnss[q][j],-.5)*pow(Mnss[q2][j],-.5))/rjksums[j];sum+=pow(rjksums[j],-1);}
}
value=pow(omega[q+q*num_sums3]*omega[q2+q2*num_sums3],.5);
if(sumsq/sum>value){sumsq=sum*value;}
if(sumsq/sum<-value){sumsq=-sum*value;}
omega[q+q2*num_sums3]=sumsq/sum;
omega[q2+q*num_sums3]=sumsq/sum;
}}
}}

for(q=0;q<num_sums3;q++)
{
if(pritraits[q]==1){fprintf(output3,"%.4f ", omega[q]*pow(omega[0]*omega[q+q*num_sums3],-.5));}
}
fprintf(output3,"%d\n",bitlength);
}

for(s=0;s<num_cors;s++){cohers[s+bit*num_cors]=0;}
for(j=bitstart;j<bitend;j++) 
{
//get covariance matrix and omega_1 / omega_11 (blank rows and columns corresponding to missing traits)
for(q=0;q<num_sums3;q++)
{
if(Mnss[q][j]>0)
{
comat[q+q*num_sums3]=omega[q+q*num_sums3]-omega[q]*omega[q]/omega[0]+sigma[q+q*num_sums3]/Mnss[q][j];
comat2[q]=omega[q]/omega[0];
}
else{comat[q+q*num_sums3]=1;comat2[q]=0;}

for(q2=q+1;q2<num_sums3;q2++)
{
if(Mnss[q][j]>0&&Mnss[q2][j]>0)
{
comat[q+q2*num_sums3]=omega[q+q2*num_sums3]-omega[q]*omega[q2]/omega[0]+sigma[q+q2*num_sums3]*pow(Mnss[q][j],-.5)*pow(Mnss[q2][j],-.5);
comat[q2+q*num_sums3]=comat[q+q2*num_sums3];
}
else{comat[q+q2*num_sums3]=0;comat[q2+q*num_sums3]=0;}
}
}

//compute inverse comat times comat2
for(q=0;q<num_sums3;q++){comat3[q]=comat2[q];}
eigen_invert(comat, num_sums3, comat4, 1, comat3, 1);

//get denominator, which is sum(comat2 * comat3)
sum=0;for(q=0;q<num_sums3;q++){sum+=comat2[q]*comat3[q];}

if(sum<=0)
{
printf("Error, a predictor has non-positive variance, please can you send Doug the screen output\n\n");
printf("Predictor %d %s\n", j+1, preds[j]);
for(q=0;q<num_sums3;q++){printf("%d size %f weight %f | ", q+1, Mnss[q][j], comat3[q]);}
printf(" --- %f\n", sum);
exit(1);
}

//compute meta rho, set effective sample size to sum (inverse of variance), and get chi
rhos4[j]=0;
for(q=0;q<num_sums3;q++)
{
value=comat3[q]/sum;
rhos4[j]+=value*Mrhos[q][j];
cohers[sumpops[q]+bit*num_cors]+=pow(value,2);
}
nss4[j]=sum;
}

//scale cohers so it sums to one
sum=0;for(s=0;s<num_cors;s++){sum+=cohers[s+bit*num_cors];}
for(s=0;s<num_cors;s++){cohers[s+bit*num_cors]=cohers[s+bit*num_cors]/sum;}

if(dougvar==3)  //reset 
{
for(bit=0;bit<bittotal;bit++)
{
for(s=0;s<num_cors;s++){cohers[s+bit*num_cors]=0;}
cohers[sumpops[0]+bit*num_cors]=1;
}
}
}   //end of bit loop

if(mtagreg==1){fclose(output3);}

if(num_cors>1)
{
printf("When estimating the final effect sizes, will use the following (squared) ancestry proportions\n");
for(s=0;s<num_cors;s++){printf("Correlations with stem %s have weight %.4f\n", corstems[s], cohers[s]);}
printf("\n");
}

sprintf(filename3,"%s.mtag",outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3,"Predictor A1 A2 Z n\n");
for(j=0;j<data_length;j++){fprintf(output3,"%s %c %c %.4f %.2f\n", preds[j], al1[j], al2[j], pow(nss4[j],.5)*rhos4[j], nss4[j]);}
fclose(output3);

free(sigma);free(omega);
free(rjksums2);free(rjksums3);free(rjksums4);free(comat);free(comat2);free(comat3);free(comat4);

///////////////////////////

