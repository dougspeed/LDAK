/*
Copyright 2026 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//deal with extra traits - non-focal traits can have missing summary statistics

///////////////////////////

indexer=malloc(sizeof(int)*num_cors);
rjksums2=malloc(sizeof(double)*data_length*num_cors);
rjksums3=malloc(sizeof(double)*data_length*num_cors*num_cors);
sigma=malloc(sizeof(double)*num_sums3*num_sums3);
omega=malloc(sizeof(double)*num_sums3*num_sums3);
comat=malloc(sizeof(double)*num_sums3*num_sums3);
comat2=malloc(sizeof(double)*num_sums3);
comat3=malloc(sizeof(double)*num_sums3);
comat4=malloc(sizeof(double)*num_sums3);

//which correlations are we using
for(s=0;s<num_cors;s++){indexer[s]=0;}
for(q=0;q<num_sums3;q++){indexer[sumpops[q]]=1;}
count=0;
for(s=0;s<num_cors;s++)
{
if(indexer[s]==1){indexer[count]=s;count++;}
}

if(count==1){printf("Recomputing taggings using estimated per-predictor heritabilities\n");}
else{printf("Calculating cross-ancestry taggings for %d ancestries\n", count);}

value=(double)bitmax/1024*bitmax/1024/1024*2*count;
if(value>1){printf("Warning, to store the correlations requires %.1f Gb\n\n", value);}
printf("\n");

Dcors_short=malloc(sizeof(unsigned short*)*count);
Mfilename=malloc(sizeof(char*)*count);
Minput=malloc(sizeof(FILE *)*count);
for(s=0;s<count;s++){Dcors_short[s]=malloc(sizeof(unsigned short)*bitmax*bitmax);Mfilename[s]=malloc(sizeof(char)*500);}

//rjksums2 and rjksums3 start at zero, then update when necessary
for(s=0;s<num_cors;s++)
{
for(j=0;j<data_length;j++){rjksums2[j+s*data_length]=0;}
for(s2=0;s2<num_cors;s2++)
{
for(j=0;j<data_length;j++){rjksums3[j+(s2+s*num_cors)*data_length]=0;}
}
}

for(bit=0;bit<bittotal;bit++)
{
bitstart=blockstarts[bit];
bitend=blockends[bit];
bitlength=blockends[bit]-blockstarts[bit];

if(bit%50==0)
{
printf("Calculating taggings for Window %d of %d\n", bit+1, bittotal);

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Calculating taggings for Window %d of %d\n", bit+1, bittotal);
fclose(output);
}

//read in all correlations and compute self values
#pragma omp parallel for private(s,s2,j,j2,k,k2,count2,value,value2,value3) schedule(static)
for(s2=0;s2<count;s2++)
{
s=indexer[s2];

sprintf(Mfilename[s2],"%s.cors.bin", corstems[s]);
if((Minput[s2]=fopen(Mfilename[s2],"rb"))==NULL)
{printf("Error re-opening %s\n\n",Mfilename[s2]);exit(1);}

//read lower triangle of focal correlations
fseeko(Minput[s2], Dindexes[s][bit], SEEK_SET);
for(k=0;k<Dsizes[s][bit][0]-1;k++)
{
count2=fread(Dcors_short[s2]+(size_t)k*Dsizes[s][bit][0]+k+1, sizeof(unsigned short), Dsizes[s][bit][0]-k-1, Minput[s2]);
if(count2!=Dsizes[s][bit][0]-k-1)
{printf("Error reading correlations for Window %d from %s\n\n", bit+1, Mfilename[s2]);exit(1);}
}
fclose(Minput[s2]);

//update rjksums2 and rjksums3 (correcting for sample size)
value=(double)(Dns[s]-1)/(Dns[s]-2);
value2=pow(Dns[s]-2,-1);
for(k=0;k<Dsizes[s][bit][1];k++)
{
k2=Duse2[s][bit][k];
rjksums2[bitstart+k2+s*data_length]++;
rjksums3[bitstart+k2+(s+s*num_cors)*data_length]+=exps[bitstart+k2];
for(j=k+1;j<Dsizes[s][bit][1];j++)
{
j2=Duse2[s][bit][j];
value3=pow(0.00005*Dcors_short[s2][(size_t)Duse[s][bit][k]*Dsizes[s][bit][0]+Duse[s][bit][j]]-1,2)*value-value2;
rjksums2[bitstart+k2+s*data_length]+=value3;
rjksums2[bitstart+j2+s*data_length]+=value3;
rjksums3[bitstart+k2+(s+s*num_cors)*data_length]+=value3*exps[bitstart+j2];
rjksums3[bitstart+j2+(s+s*num_cors)*data_length]+=value3*exps[bitstart+k2];
}
}
}

//get cross correlations
for(s2=0;s2<count-1;s2++)
{
s=indexer[s2];

//set cors to zero
#pragma omp parallel for private(k,j) schedule(static)
for(k=0;k<bitlength;k++)
{
for(j=k;j<bitlength;j++){cors[(size_t)k*bitlength+j]=0;}
}

//extract correlations for ancestry s (remember signs)
#pragma omp parallel for private(k,j) schedule(static)
for(k=0;k<Dsizes[s][bit][1];k++)
{
cors[(size_t)Duse2[s][bit][k]*bitlength+Duse2[s][bit][k]]=1;
for(j=k+1;j<Dsizes[s][bit][1];j++)
{cors[(size_t)Duse2[s][bit][k]*bitlength+Duse2[s][bit][j]]=0.00005*Dcors_short[s2][(size_t)Duse[s][bit][k]*Dsizes[s][bit][0]+Duse[s][bit][j]]-1;}
}

//get cross-correlations with later ancestries
#pragma omp parallel for private(s,s3,j,j2,k,k2,value,value2,value3,value4) schedule(static)
for(s3=s2+1;s3<count;s3++)
{
s=indexer[s3];

value3=pow((Dns[indexer[s2]]-2)*(Dns[s]-2),.5);
value=(value3+1)/value3;
value2=pow(value3,-1);
for(k=0;k<Dsizes[s][bit][1];k++)
{
k2=Duse2[s][bit][k];
if(cors[(size_t)k2*bitlength+k2]==1){rjksums3[bitstart+k2+(indexer[s2]+s*num_cors)*data_length]+=exps[bitstart+k2];}
for(j=k+1;j<Dsizes[s][bit][1];j++)
{
j2=Duse2[s][bit][j];
value3=(0.00005*Dcors_short[s3][(size_t)Duse[s][bit][k]*Dsizes[s][bit][0]+Duse[s][bit][j]]-1)*Dsigns[s][bit][j]*Dsigns[s][bit][k];
value4=cors[(size_t)k2*bitlength+j2]*value3*value;
if(value4>0){value4-=value2;}
else{value4+=value2;}
rjksums3[bitstart+k2+(indexer[s2]+s*num_cors)*data_length]+=value4*exps[bitstart+j2];
rjksums3[bitstart+j2+(indexer[s2]+s*num_cors)*data_length]+=value4*exps[bitstart+k2];
}
}
}
}
}   //end of bit loop
printf("\n");

//ensure rjksums2 are not below 1
#pragma omp parallel for private(s,s2,j,count2,count3) schedule(static)
for(s2=0;s2<count;s2++)
{
s=indexer[s2];
for(j=0;j<Dnp2[s];j++)
{
if(rjksums2[Dkp2[s][j]+s*data_length]<1){rjksums2[Dkp2[s][j]+s*data_length]=1;}
}
}

//copy off-diagonals of rjksums3
for(s2=0;s2<count-1;s2++)
{
for(s3=s2+1;s3<count;s3++)
{copy_matrix(1, data_length, rjksums3+(indexer[s2]+indexer[s3]*num_cors)*data_length, rjksums3+(indexer[s3]+indexer[s2]*num_cors)*data_length, 0, NULL);}
}

for(s=0;s<count;s++){free(Dcors_short[s]);free(Mfilename[s]);}free(Dcors_short);free(Mfilename);free(Minput);

////////

printf("Computing omega and sigma, the MTAG covariance matrices\n");

//start with diagonals, checking values not too wild
#pragma omp parallel for private(q) schedule(dynamic)
for(q=0;q<num_sums3;q++)
{
if(mtagcopy==0) //fix intercept at one, estimate slope
{solve_sums_lite(sigma+q+q*num_sums3, omega+q+q*num_sums3, rjksums2+sumpops[q]*data_length, rjksums3+(sumpops[q]+sumpops[q]*num_cors)*data_length, Mnss[q], Mchis[q], data_length, 0, q);}
else    //estimate intercept and slope
{solve_sums_lite(sigma+q+q*num_sums3, omega+q+q*num_sums3, rjksums2+sumpops[q]*data_length, rjksums3+(sumpops[q]+sumpops[q]*num_cors)*data_length, Mnss[q], Mchis[q], data_length, 1, q);}

if(sigma[q+q*num_sums3]<0.5)
{
printf("Warning, trait %d has a very low estimated intercept (%.4f); has been increased to 0.5\n", q+1, sigma[q+q*num_sums3]);
sigma[q+q*num_sums3]=0.5;
}
if(omega[q+q*num_sums3]<0.01)
{
printf("Warning, trait %d has a very low estimated heritability (%.4f); has been increased to 0.01\n", q+1, omega[q+q*num_sums3]);
omega[q+q*num_sums3]=0.01;
}
if(omega[q+q*num_sums3]>0.8)
{printf("Warning, trait %d has a very high estimated heritability (%.4f)\n", q+1, omega[q+q*num_sums3]);}
}

//now focal vs other traits, checking values not too wild
#pragma omp parallel for private(q,value) schedule(dynamic)
for(q=1;q<num_sums3;q++)
{
if(sumpops[0]==sumpops[q])  //estimate intercept and slope
{solve_cors_lite(sigma+q, omega+q, rjksums2+sumpops[0]*data_length, rjksums2+sumpops[q]*data_length, rjksums3+(sumpops[0]+sumpops[q]*num_cors)*data_length, Mnss[0], Mnss[q], Mchis[0], Mchis[q], Mrhos[0], Mrhos[q], data_length, 1, 0, q);}
else //fix intercept at zero, estimate slope
{solve_cors_lite(sigma+q, omega+q, rjksums2+sumpops[0]*data_length, rjksums2+sumpops[q]*data_length, rjksums3+(sumpops[0]+sumpops[q]*num_cors)*data_length, Mnss[0], Mnss[q], Mchis[0], Mchis[q], Mrhos[0], Mrhos[q], data_length, 0, 0, q);}

if(sigma[q]>0.99)
{
printf("Warning, trait pair (1-%d) has estimated intercept %.4f; has been reduced to 0.99\n", q+1, sigma[q]);
sigma[q]=0.99;
}
if(sigma[q]<-0.99)
{
printf("Warning, trait pair (1-%d) has estimated intercept %.4f; has been increased to -0.99\n", q+1, sigma[q]);
sigma[q]=-0.99;
}

value=omega[q]*pow(omega[0]*omega[q+q*num_sums3],-.5);
if(value>0.99)
{
printf("Warning, trait pair (1-%d) has estimated correlation %.4f; has been reduced to 0.99\n", q+1, value);
omega[q]=0.99*pow(omega[0]*omega[q+q*num_sums3],.5);
}
if(value<-0.99)
{
printf("Warning, trait pair (1-%d) has estimated correlation %.4f; has been increased to -0.99\n", q+1, value);
omega[q]=-0.99*pow(omega[0]*omega[q+q*num_sums3],.5);
}

sigma[q*num_sums3]=sigma[q];
omega[q*num_sums3]=omega[q];
}

if(mtagforce==0)    //blank secondary traits with low correlation (record this by setting pritraits=2)
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

//get remaining off-diagonal values (if required), checking values not too wild
for(q=1;q<num_sums3;q++)
{
//set values to zero, then update for non-excluded traits
for(q2=q+1;q2<num_sums3;q2++){sigma[q+q2*num_sums3]=0;sigma[q2+q*num_sums3]=0;omega[q+q2*num_sums3]=0;omega[q2+q*num_sums3]=0;}

if(pritraits[q]!=2)
{
#pragma omp parallel for private(q2,value) schedule(dynamic)
for(q2=q+1;q2<num_sums3;q2++)
{
if(pritraits[q2]!=2)
{
if(sumpops[q]==sumpops[q2])    //estimate intercept and slope
{solve_cors_lite(sigma+q+q2*num_sums3, omega+q+q2*num_sums3, rjksums2+sumpops[q]*data_length, rjksums2+sumpops[q2]*data_length, rjksums3+(sumpops[q]+sumpops[q2]*num_cors)*data_length, Mnss[q], Mnss[q2], Mchis[q], Mchis[q2], Mrhos[q], Mrhos[q2], data_length, 1, q, q2);}
else //fix intercept at zero, estimate slope
{solve_cors_lite(sigma+q+q2*num_sums3, omega+q+q2*num_sums3, rjksums2+sumpops[q]*data_length, rjksums2+sumpops[q2]*data_length, rjksums3+(sumpops[q]+sumpops[q2]*num_cors)*data_length, Mnss[q], Mnss[q2], Mchis[q], Mchis[q2], Mrhos[q], Mrhos[q2], data_length, 0, q, q2);}

if(sigma[q+q2*num_sums3]>0.99)
{
printf("Warning, trait pair (%d-%d) has estimated intercept %.4f; has been reduced to 0.99\n", q+1, q+2, sigma[q+q2*num_sums3]);
sigma[q+q2*num_sums3]=0.99;
}
if(sigma[q+q2*num_sums3]<-0.99)
{
printf("Warning, trait pair (%d-%d) has estimated intercept %.4f; has been increased to -0.99", q+1, q+2, sigma[q+q2*num_sums3]);
sigma[q+q2*num_sums3]=-0.99;
}

value=omega[q+q2*num_sums3]*pow(omega[q+q*num_sums3]*omega[q2+q2*num_sums3],-.5);
if(value>0.99)
{
printf("Warning, trait pair (%d-%d) has estimated correlation %.4f; has been reduced to 0.99\n", q+1, q2+1, value);
omega[q+q2*num_sums3]=0.99*pow(omega[q+q*num_sums3]*omega[q2+q2*num_sums3],.5);
}
if(value<-0.99)
{
printf("Warning, trait pair (%d-%d) has estimated correlation %.4f; has been increased to -0.99", q+1, q2+1, value);
omega[q+q2*num_sums3]=-0.99*pow(omega[q+q*num_sums3]*omega[q2+q2*num_sums3],.5);
}

sigma[q2+q*num_sums3]=sigma[q+q2*num_sums3];
omega[q2+q*num_sums3]=omega[q+q2*num_sums3];
}
}
}
}

if(mtagmeta==1) //all primary traits should have same heritability as focal trait, and have no sample overlap
{
for(q=1;q<num_sums3;q++)
{
if(pritraits[q]==1)
{
omega[q+q*num_sums3]=omega[0];
sigma[q]=0;omega[q]=omega[0];
sigma[q*num_sums3]=0;omega[q*num_sums3]=omega[0];

for(q2=q+1;q2<num_sums3;q2++)      
{
if(pritraits[q2]==1)
{
sigma[q+q2*num_sums3]=0;omega[q+q2*num_sums3]=omega[0];
sigma[q2+q*num_sums3]=0;omega[q2+q*num_sums3]=omega[0];
}}
}}
}

/*
printf("\nhere is omega\n");
print_top(omega,num_sums3,num_sums3,num_sums3);
printf("here is sigma\n");
print_top(sigma,num_sums3,num_sums3,num_sums3);
printf("here is correlation matrix\n");
for(q=0;q<num_sums3;q++)
{
for(q2=0;q2<num_sums3;q2++)
{
printf("%.4f ", omega[q+q2*num_sums3]*pow(omega[q+q*num_sums3]*omega[q2+q2*num_sums3],-.5));
}
printf("\n");
}
printf("\n");
*/

////////

//note that for individual snps, will use omega*exps[j] or omega*exps[j]*rjksums3[j], and will use sigma/nss[j]
//comat contains (approx) average covariance matrix for all (non-trivial) traits, comat2 contains omega_1 / omega_11 (or zero if trait trivial)
for(q=0;q<num_sums3;q++)
{
sum=0;for(j=0;j<data_length;j++){sum+=Mnss[q][j];}
mean=sum/data_length;

if(mean>0)
{
comat[q+q*num_sums3]=omega[q+q*num_sums3]/data_length-omega[q]*omega[q]/omega[0]/data_length+sigma[q+q*num_sums3]/mean;
comat2[q]=omega[q]/omega[0];
}
else{comat[q+q*num_sums3]=1;comat2[q]=0;}

for(q2=q+1;q2<num_sums3;q2++)
{
sum2=0;for(j=0;j<data_length;j++){sum2+=Mnss[q2][j];}
mean2=sum2/data_length;
if(mean>0&&mean2>0)
{
comat[q+q2*num_sums3]=omega[q+q2*num_sums3]/data_length-omega[q]*omega[q2]/omega[0]/data_length+sigma[q+q2*num_sums3]*pow(mean,-.5)*pow(mean2,-.5);
comat[q2+q*num_sums3]=comat[q+q2*num_sums3];
}
else{comat[q+q2*num_sums3]=0;comat[q2+q*num_sums3]=0;}
}
}

if(mtagforce==1)  //will use all (non-trivial) traits
{mtag_select(comat, comat2, comat3, comat4, num_sums3, pritraits, -9999);}
else    //use focal trait, primary traits and selected secondary traits
{mtag_select(comat, comat2, comat3, comat4, num_sums3, pritraits, mtagpercent/100);}

printf("Here are the estimated genetic correlations, intercepts, importances and weightings\n");
for(q=0;q<num_sums3;q++)
{printf("Trait %d (%s): correlation %.4f, intercept %.4f, importance %.4f and weighting %.4f\n", q+1, sumstems[q], omega[q]*pow(omega[0]*omega[q+q*num_sums3],-.5), sigma[q], comat4[q]/comat4[0], comat3[q]);}
printf("\n");

if(num_focals==1){sprintf(filename3,"%s.mtag.weights",outfile);}
else{sprintf(filename3,"%s.focal%d.mtag.weights", outfile, cur_focal+1);}
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3,"Summary_Statistics File Correlation Effective_Sample_Size Importance Weighting\n");
for(q=0;q<num_sums3;q++){fprintf(output3,"%d %s %.4f %.2f %.4f %.4f\n", q+1, sumstems[q], omega[q]*pow(omega[0]*omega[q+q*num_sums3],-.5), comat4[q], comat4[q]/comat4[0], comat3[q]);}
fclose(output3);

//how many primary and secondary traits are being used?
count=0;for(q=1;q<num_sums3;q++){count+=(pritraits[q]==1);}
count2=0;for(q=1;q<num_sums3;q++){count2+=(pritraits[q]==0&&comat3[q]!=0);}

if(num_focals==1){sprintf(filename3,"%s.mtag.summary",outfile);}
else{sprintf(filename3,"%s.focal%d.mtag.summary", outfile, cur_focal+1);}
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3,"Using the focal summary statistics\nall %d primary summary statistics\nand %d secondary summary statistics\n", count, count2);
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

if(impsums==1)  //impute missing summary statistic for extra trait with non-zero weight (if necessary)
{
for(q=1;q<num_sums3;q++)
{
if(pritraits[q]==1)
{
count=0;for(j=0;j<data_length;j++){count+=(Mnss[q][j]==0&&rjksums2[j+sumpops[q]*data_length]>0);}
if(count>0)
{
printf("Imputing summary statistics for (up to) %d predictors for Trait %d (%s)\n", count, q+1, sumstems[q]);
impute_sums(q, Mnss, Mrhos, Mchis, data_length, bittotal, blockstarts, blockends, sumpops[q], NULL, corstems, -9999, Dindexes, Dsizes, Duse, Duse2, Dsigns, shrink, outfile, cors_short, cors, preds);
}
}
}
}

////////

if(mtagreg==1)
{
if(num_focals==1){sprintf(filename3,"%s.mtag.regions",outfile);}
else{sprintf(filename3,"%s.focal%d.mtag.regions", outfile, cur_focal+1);}
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
for(q=0;q<num_sums3;q++)
{
if(pritraits[q]==1){fprintf(output3,"%.4f ", omega[q]*pow(omega[0]*omega[q+q*num_sums3],-.5));}
}
fprintf(output3,"%d\n", data_length);
}

//ready to combine summary statistics for focal trait
printf("Computing the MTAG summary statistics\n\n");

wcount=0;
for(bit=0;bit<bittotal;bit++)
{
bitstart=blockstarts[bit];
bitend=blockends[bit];
bitlength=blockends[bit]-blockstarts[bit];

if(mtagreg==1)  //recompute omega (for contributing traits) - must have mtagcopy=0
{
#pragma omp parallel for private(q) schedule(dynamic)
for(q=0;q<num_sums3;q++)
{
if(pritraits[q]==1)
{
solve_sums_lite(sigma+q+q*num_sums3, omega+q+q*num_sums3, rjksums2+sumpops[q]*data_length+bitstart, rjksums3+(sumpops[q]+sumpops[q]*num_cors)*data_length+bitstart, Mnss[q]+bitstart, Mchis[q]+bitstart, bitlength, 2, q);
if(omega[q+q*num_sums3]<=0){omega[q+q*num_sums3]=1e-10;}
}}

#pragma omp parallel for private(q,q2,value) schedule(dynamic)
for(q=0;q<num_sums3;q++)
{
if(pritraits[q]==1)
{
for(q2=q+1;q2<num_sums3;q2++)
{
if(pritraits[q2]==1)
{
solve_cors_lite(sigma+q+q2*num_sums3, omega+q+q2*num_sums3, rjksums2+sumpops[q]*data_length+bitstart, rjksums2+sumpops[q2]*data_length+bitstart, rjksums3+(sumpops[q]+sumpops[q2]*num_cors)*data_length+bitstart, Mnss[q]+bitstart, Mnss[q2]+bitstart, Mchis[q]+bitstart, Mchis[q2]+bitstart, Mrhos[q]+bitstart, Mrhos[q2]+bitstart, bitlength, 2, q, q2);

value=omega[q+q2*num_sums3]*pow(omega[q+q*num_sums3]*omega[q2+q2*num_sums3],-.5);
if(value>0.99){omega[q+q2*num_sums3]=0.99*pow(omega[q+q*num_sums3]*omega[q2+q2*num_sums3],.5);}
if(value<-0.99){omega[q+q2*num_sums3]=-0.99*pow(omega[q+q*num_sums3]*omega[q2+q2*num_sums3],.5);}

omega[q2+q*num_sums3]=omega[q+q2*num_sums3];
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
if(dougvar!=1)
{
comat[q+q*num_sums3]=omega[q+q*num_sums3]*exps[j]-omega[q]*omega[q]/omega[0]*exps[j]+sigma[q+q*num_sums3]/Mnss[q][j];
comat2[q]=omega[q]/omega[0];
}
else
{
value=rjksums3[j+(sumpops[q]+sumpops[q]*num_cors)*data_length];
value2=rjksums3[j+(sumpops[0]+sumpops[q]*num_cors)*data_length];
value3=rjksums3[j+(sumpops[0]+sumpops[0]*num_cors)*data_length];
comat[q+q*num_sums3]=omega[q+q*num_sums3]*value-omega[q]*omega[q]/omega[0]*value2*value2/value3+sigma[q+q*num_sums3]/Mnss[q][j];
comat2[q]=omega[q]/omega[0]*value2/value3;
}
}
else{comat[q+q*num_sums3]=1;comat2[q]=0;}

for(q2=q+1;q2<num_sums3;q2++)
{
if(Mnss[q][j]>0&&Mnss[q2][j]>0)
{
if(dougvar!=1)
{
comat[q+q2*num_sums3]=omega[q+q2*num_sums3]*exps[j]-omega[q]*omega[q2]/omega[0]*exps[j]+sigma[q+q2*num_sums3]*pow(Mnss[q][j],-.5)*pow(Mnss[q2][j],-.5);
}
else
{
value=rjksums3[j+(sumpops[q]+sumpops[q2]*num_cors)*data_length];
value2=rjksums3[j+(sumpops[0]+sumpops[q]*num_cors)*data_length];
value3=rjksums3[j+(sumpops[0]+sumpops[0]*num_cors)*data_length];
value4=rjksums3[j+(sumpops[0]+sumpops[q2]*num_cors)*data_length];
comat[q+q2*num_sums3]=omega[q+q2*num_sums3]*value-omega[q]*omega[q2]/omega[0]*value2*value4/value3+sigma[q+q2*num_sums3]*pow(Mnss[q][j],-.5)*pow(Mnss[q2][j],-.5);
}
comat[q2+q*num_sums3]=comat[q+q2*num_sums3];
}
else{comat[q+q2*num_sums3]=0;comat[q2+q*num_sums3]=0;}
}
}

//compute inverse comat times comat2
for(q=0;q<num_sums3;q++){comat3[q]=comat2[q];}
info=eigen_solve_check(comat, num_sums3, comat4, comat3);
if(info>0)
{
if(wcount<5){printf("Warning, the matrix for predictor %d (%s) has %d non-positive eigenvectors\n", j+1, preds[j], info);}
wcount++;
}

//get denominator, which is sum(comat2 * comat3)
sum=0;for(q=0;q<num_sums3;q++){sum+=comat2[q]*comat3[q];}

if(sum<=0)
{
printf("Error, a predictor has non-positive variance, please can you send Doug the screen output\n\n");
printf("Predictor %d %s\n", j+1, preds[j]);
for(q=0;q<num_sums3;q++){printf("%d size %f orig %f weight %f | ", q+1, Mnss[q][j], comat2[q], comat3[q]);}
printf(" --- %f\n", sum);

printf("exps %e\n", exps[j]);
print_top(rjksums3+j,1,num_cors*num_cors,data_length);
exit(1);
}

//compute meta rho and set effective sample size to sum (inverse of variance)
rhos4[j]=0;
for(q=0;q<num_sums3;q++){rhos4[j]+=comat3[q]/sum*Mrhos[q][j];}
nss4[j]=sum;
chis4[j]=nss4[j]*pow(rhos4[j],2)/(1-pow(rhos4[j],2));

//add on contributions to cohers - be careful for focal trait when using metasum
if(metasum==0){cohers[sumpops[0]+bit*num_cors]+=comat3[0]/sum;}
else
{
for(s=0;s<num_cors;s++){cohers[s+bit*num_cors]+=cohers2[s]*comat3[0]/sum;}
}

for(q=1;q<num_sums3;q++){cohers[sumpops[q]+bit*num_cors]+=comat3[q]/sum;}
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
if(wcount>5){printf("In total, %d matrices had non-positive eigenvalues\n", count2);}
if(wcount>0){printf("\n");}

if(mtagreg==1){fclose(output3);}

if(num_cors>1)
{
printf("Here is the ancestry composition of the combined summary statistics\n");
for(s=0;s<num_cors;s++){printf("Correlations with stem %s have average contribution %.4f\n", corstems[s], cohers[s]);}
printf("\n");
}

if(num_focals==1){sprintf(filename3,"%s.mtag",outfile);}
else{sprintf(filename3,"%s.focal%d.mtag", outfile, cur_focal+1);}
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3,"Predictor\tA1\tA2\tZ\tn\n");
for(j=0;j<data_length;j++){fprintf(output3,"%s\t%c\t%c\t%.4f\t%.2f\n", preds[j], al1[j], al2[j], pow(nss4[j]/(1-pow(rhos4[j],2)),.5)*rhos4[j], nss4[j]);}
fclose(output3);

free(indexer);free(rjksums2);free(rjksums3);free(sigma);free(omega);free(comat);free(comat2);free(comat3);free(comat4);

///////////////////////////

