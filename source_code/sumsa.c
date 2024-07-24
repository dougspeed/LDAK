/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Summary heritability

///////////////////////////

if(strcmp(powfile,"blank")!=0)	//deal with powers (will do first, as will find out faster if invalid)
{
printf("Reading %d powers from %s\n\n",num_pows,powfile);
powers=malloc(sizeof(double)*num_pows);
read_values(powfile,powers,num_pows,NULL,1,0,0);
}

//read in tag details
count=countrows(tagfile)-2-num_parts;
model_warn(count, num_parts);

catlabels=malloc(sizeof(char *)*num_parts);
spreds=malloc(sizeof(char*)*count);
sal1=malloc(sizeof(char)*count);
sal2=malloc(sizeof(char)*count);
stags=malloc(sizeof(double)*count);
svars=malloc(sizeof(double*)*num_parts);
ssums=malloc(sizeof(double*)*num_parts);
for(q=0;q<num_parts;q++)
{
catlabels[q]=malloc(sizeof(char)*500);
svars[q]=malloc(sizeof(double)*count);
ssums[q]=malloc(sizeof(double)*(num_parts+2));
}

//read tagfile - if extracting, count will be reset
count=read_tagfile(tagfile, catlabels, spreds, sal1, sal2, stags, svars, ssums, num_parts, count, parttype, altfile, bpredfile, cpredfile, tagone);

//get summary stats
snss=malloc(sizeof(double)*count);
schis=malloc(sizeof(double)*count);
srhos=malloc(sizeof(double)*count);

printf("Reading summary statistics from %s\n", sumsfile);
read_sumsfile(sumsfile, snss, schis, srhos, count, spreds, sal1, sal2, tagfile, amb, fixn, scaling, 1, checksums);
printf("First few stats and ns are: %s %.3f %.1f", spreds[0], schis[0], snss[0]);
for(j=1;j<count;j++){if(j<3){printf(" | %s %.3f %.1f", spreds[j], schis[j], snss[j]);}}
printf("\n\n");

//print out overlap
count2=0;for(j=0;j<count;j++){count2+=(snss[j]!=-9999);}
sprintf(filename,"%s.overlap",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fprintf(output,"Tagging_File_Predictors %d\nSummary_Statistic_Predictors %d\nOverlap_Proportion %.4f\n", count, count2, (double)count2/count);
fclose(output);

if(cutoff!=-9999)	//look for big predictors
{
count2=0;
for(j=0;j<count;j++)
{
if(snss[j]!=-9999&&pow(srhos[j],2)>=cutoff){snss[j]=-9999;count2++;}
}
if(count2==0){printf("There are no predictors that explain at least %.6f of phenotypic variance\n\n", cutoff);}
if(count2==1){printf("There is 1 predictor that explains at least %.6f of phenotypic variance; this will be excluded\n\n", cutoff);}
if(count2>1){printf("There are %d predictors that explain at least %.6f of phenotypic variance; these will be excluded\n\n", count2, cutoff);}
}

//count number of predictors we will use
count2=0;for(j=0;j<count;j++){count2+=(snss[j]!=-9999);}
if(count2==0){printf("Error, there are no predictors\n\n");exit(1);}
//if(count2<3){printf("Error, there are only %d predictors, can not continue\n\n", count2);exit(1);}
if(count2<20000){printf("Warning, there are usually tens of thousands of predictors (not %d)\n\n", count2);}

if(count2<count)	//squeeze down and reset count
{
count2=0;
for(j=0;j<count;j++)
{
if(snss[j]!=-9999)
{
if(count2!=j)
{
free(spreds[count2]);copy_string(spreds,count2,spreds[j]);
sal1[count2]=sal1[j];
sal2[count2]=sal2[j];
stags[count2]=stags[j];
for(q=0;q<num_parts;q++){svars[q][count2]=svars[q][j];}
snss[count2]=snss[j];
schis[count2]=schis[j];
srhos[count2]=srhos[j];
}
count2++;
}}
for(j=count2;j<count;j++){free(spreds[j]);}
count=count2;
}

count2=0;count3=0;
for(j=0;j<count;j++)
{
if(stags[j]<1){count2++;count3+=(stags[j]<0);stags[j]=1;}
}
if(count2>0){printf("Warning, %d SNPs have tagging less than 1 (%d less than 0); these have been replaced with 1\n\n", count2, count3);}

count2=0;
for(j=0;j<count;j++)
{
if(schis[j]==0){count2++;schis[j]=1e-6;}
}
if(count2>0){printf("Warning, %d SNPs have test statistic 0; these have been replaced with 1e-6\n\n", count2);}

gif=get_sum_stats(stags, snss, schis, srhos, count, 0, cutoff);

if(strcmp(cvsfile,"blank")!=0)	//have cv predictors
{
count2=countrows(cvsfile);
printf("Reading list of %d cross-validation predictors from %s\n", count2, cvsfile);
wantpreds=malloc(sizeof(char*)*count2);
read_strings(cvsfile, wantpreds, count2, NULL, 1, 0);
cvindex=malloc(sizeof(int)*count2);
ncv=find_strings(spreds, count, wantpreds, count2, cvindex, NULL, NULL, NULL, NULL, NULL, 3);
if(ncv==0){printf("Error, none of these are in the tagging file and have summary statistics\n\n");exit(1);}
if(ncv<count2){printf("%d of these are in the tagging file and have summary statistics\n\n",ncv);}
else{printf("All of these are in the tagging file and have summary statistics\n\n");}
for(j=0;j<count2;j++){free(wantpreds[j]);}free(wantpreds);

cvexps=malloc(sizeof(double)*ncv);
}
else{ncv=0;}

////////

if(prev!=-9999)	//get factor - must have ascer
{factor=get_factor(NULL, -9999, prev, ascer, outfile);}

//set num_parts_use and num_anals
if(divide==-9999){num_parts_use=num_parts;num_anals=1;}
else{num_parts_use=divide;num_anals=num_parts/divide;}

//set total and allocate
total=num_parts_use+gcon+cept;
stats=malloc(sizeof(double)*(total+1+num_parts_use)*3);
likes=malloc(sizeof(double)*11);
cohers=malloc(sizeof(double)*num_parts_use*num_parts_use);
cohers2=malloc(sizeof(double)*num_parts_use);
influs=malloc(sizeof(double)*num_parts_use);

//following two options are mutually exclusive
if(strcmp(catfile,"blank")!=0)	//first analysis will contain num_reds categories
{
svars2=malloc(sizeof(double*)*num_reds);ssums2=malloc(sizeof(double*)*num_reds);
for(q=0;q<num_reds;q++){svars2[q]=malloc(sizeof(double)*count);ssums2[q]=malloc(sizeof(double)*(num_reds+2));}
stats2=malloc(sizeof(double)*(num_reds+gcon+cept+1+num_reds));
}
if(divide!=-9999)	//each analysis will contain divide categories
{
svars2=malloc(sizeof(double*)*divide);ssums2=malloc(sizeof(double*)*divide);
for(q=0;q<divide;q++){svars2[q]=malloc(sizeof(double)*count);ssums2[q]=malloc(sizeof(double)*(divide+2));}
if(strcmp(powfile,"blank")!=0){powlikes=malloc(sizeof(double)*num_pows*6);}
}

//ready to perform analyses

for(m=0;m<num_anals;m++)
{
//deal with progress file
if(divide==-9999){sprintf(filename,"%s.progress",outfile);}
else{sprintf(filename,"%s.%d.progress",outfile,m+1);}
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fclose(output);

if(divide==-9999)	//might be using catfile, taufile or ldsc (but not more than one)
{
if(strcmp(catfile,"blank")==0&&strcmp(taufile,"blank")==0&&ldsc==0)	//the normal solver
{solve_sums(stats, likes, cohers, cohers2, influs, num_parts, gcon, cept, num_blocks, count, ncv, cvindex, cvexps, stags, svars, ssums, snss, schis, parttype, tol, maxiter, chisol, 0, filename);}

if(strcmp(catfile,"blank")!=0)
{
//solve for num_reds categories (indexed in keepparts2) - will not be jackknifing
for(q=0;q<num_reds;q++)
{
for(j=0;j<count;j++){svars2[q][j]=svars[keepparts2[q]][j];}
for(q2=0;q2<num_reds;q2++){ssums2[q][q2]=ssums[keepparts2[q]][keepparts2[q2]];}
ssums2[q][num_reds]=ssums[keepparts2[q]][num_parts];
ssums2[q][num_reds+1]=ssums[keepparts2[q]][num_parts+1];
}

if(num_reds==1){printf("Solving for only one category\n");}
else{printf("Solving for only %d categories\n", num_reds);}
solve_sums(stats2, NULL, NULL, NULL, NULL, num_reds, gcon, cept, -9999, count, ncv, cvindex, NULL, stags, svars2, ssums2, snss, schis, parttype, tol, maxiter, chisol, 1, filename);

//load coefficients into stats and solve for full model
for(q=0;q<num_parts;q++){stats[q]=0;}
for(q=0;q<num_reds;q++){stats[keepparts2[q]]=stats2[q];}
if(gcon==1){stats[num_parts]=stats2[num_reds];}
if(cept==1){stats[num_parts+gcon]=stats2[num_reds+gcon];}

if(num_parts==2){printf("Solving for both categories\n");}
else{printf("Solving for all %d categories\n", num_parts);}
solve_sums(stats, likes, cohers, cohers2, influs, num_parts, gcon, cept, num_blocks, count, ncv, cvindex, cvexps, stags, svars, ssums, snss, schis, parttype, tol, maxiter, chisol, 2, filename);
}

if(strcmp(taufile,"blank")!=0)	//have already checked taufile is correct length
{
//stats should contain heritabilities (tau x sums), gc and intercept-1
read_values(taufile,stats,num_parts+gcon+cept,NULL,1,0,0);
for(q=0;q<num_parts;q++){stats[q]*=ssums[q][q];}
if(cept==1){stats[num_parts+gcon]-=1;}

solve_sums(stats, likes, NULL, NULL, influs, num_parts, gcon, cept, -9999, count, ncv, cvindex, cvexps, stags, svars, ssums, snss, schis, parttype, tol, maxiter, chisol, 3, filename);

//convenient to set sds and variances to zero
for(q=0;q<total+1+num_parts;q++){stats[q+total+1+num_parts]=0;stats[q+2*(total+1+num_parts)]=0;}
for(q=0;q<num_parts*num_parts;q++){cohers[q]=0;}
for(q=0;q<num_parts;q++){cohers2[q]=0;}
}

if(ldsc==1)	//try to recreate ldsc
{solve_sums(stats, likes, cohers, cohers2, influs, num_parts, gcon, cept, num_blocks, count, ncv, cvindex, cvexps, stags, svars, ssums, snss, schis, parttype, tol, maxiter, chisol, 4, filename);}
}
else	//using divide - can not be using catfile, taufile or ldsc (but can be using cvpredictors)
{
for(q=0;q<divide;q++)
{
for(j=0;j<count;j++){svars2[q][j]=svars[q+m*divide][j];}
for(q2=0;q2<divide;q2++){ssums2[q][q2]=ssums[q+m*divide][q2+m*divide];}
ssums2[q][divide]=ssums[q+m*divide][num_parts];
}
sum=0;for(q=0;q<divide;q++){sum+=ssums2[q][divide];}
for(q=0;q<divide;q++){ssums2[q][divide+1]=ssums2[q][divide]/sum;}

printf("Analysing Category %d of %d\n", m+1, num_anals);
if(m==0||uptaus==0)	//solve first category normally
{solve_sums(stats, likes, cohers, cohers2, influs, divide, gcon, cept, num_blocks, count, ncv, cvindex, cvexps, stags, svars2, ssums2, snss, schis, 1, tol, maxiter, chisol, 0, filename);}
else	//start values based on coefficients from previous run
{
for(q=0;q<divide;q++){stats[q]*=ssums[q+m*divide][q+m*divide]/ssums[q+(m-1)*divide][q+(m-1)*divide];}
solve_sums(stats, likes, cohers, cohers2, influs, divide, gcon, cept, num_blocks, count, ncv, cvindex, cvexps, stags, svars2, ssums2, snss, schis, 1, tol, maxiter, chisol, 5, filename);}

if(strcmp(powfile,"blank")!=0)
{
powlikes[m+2*num_anals]=stats[total];
powlikes[m+3*num_anals]=stats[total+total+1+divide];
powlikes[m+4*num_anals]=likes[0];
powlikes[m+5*num_anals]=likes[1];
}
}

////////

if(strcmp(cvsfile,"blank")!=0)	//will only save cv predictions and correlation
{
if(divide==-9999){sprintf(filename,"%s.cv.exps",outfile);}
else{sprintf(filename,"%s.%d.cv.exps",outfile,m+1);}
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fprintf(output,"Predictor Stat Expectation Weighting\n");
for(j=0;j<ncv;j++){fprintf(output,"%s %.4f %.4f %.4f\n", spreds[cvindex[j]], schis[cvindex[j]], cvexps[j], 1.0/stags[cvindex[j]]);}
fclose(output);

sum=0;sum2=0;sumsq=0;sumsq2=0;sumsq3=0;
for(j=0;j<ncv;j++)
{
sum+=schis[cvindex[j]];sum2+=cvexps[j];
sumsq+=pow(schis[cvindex[j]],2);sumsq2+=pow(cvexps[j],2);sumsq3+=schis[cvindex[j]]*cvexps[j];
}
mean=sum/ncv;mean2=sum2/ncv;
value=(sumsq3-ncv*mean*mean2)/pow(sumsq-ncv*mean*mean,.5)/pow(sumsq2-ncv*mean2*mean2,.5);

if(divide==-9999){sprintf(filename2,"%s.cv.cor",outfile);}
else{sprintf(filename2,"%s.%d.cv.cor",outfile,m+1);}
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2,"Num_Predictors %d\nCorrelation %.4f\n", ncv, value);
fclose(output2);

printf("Across the %d predictors, the correlation between the true and predicted test statistics is %.4f (predicted test are saved in %s)\n\n", ncv, value, filename);
}
else	//saving most things
{
//set ssums_use
if(divide==-9999){ssums_use=ssums;}
else{ssums_use=ssums2;}

//get sum of influences, and its sd
sum=0;sum2=0;
for(q=0;q<num_parts_use;q++)
{
sum+=influs[q]*stats[q];
for(q2=0;q2<num_parts_use;q2++){sum2+=influs[q]*influs[q2]*cohers[q+q2*num_parts_use];}
}

if(divide==-9999){sprintf(filename,"%s.hers",outfile);}
else{sprintf(filename,"%s.%d.hers",outfile,m+1);}
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fprintf(output,"Component Heritability SD Influence SD\n");
for(q=0;q<num_parts_use;q++)
{
if(parttype==0&&q<num_parts_use-1){fprintf(output,"Her_A%d ",q+1);}
if(parttype==0&&q==num_parts_use-1){fprintf(output,"Her_Base ");}
if(parttype==1){fprintf(output,"Her_P%d ",q+1);}
fprintf(output,"%.6f %.6f %.6f %.6f\n", stats[q], stats[q+total+1+num_parts_use], influs[q]*stats[q], fabs(influs[q])*stats[q+total+1+num_parts_use]);
}
if(sum2>0){fprintf(output,"Her_All %.6f %.6f %.6f %.6f\n", stats[total], stats[total+total+1+num_parts_use], sum, pow(sum2,.5));}
else{fprintf(output,"Her_All %.6f %.6f %.6f NA\n", stats[total], stats[total+total+1+num_parts_use], sum);}
fclose(output);

if(divide==-9999){sprintf(filename2,"%s.cats",outfile);}
else{sprintf(filename2,"%s.%d.cats",outfile,m+1);}
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2,"Component Heritability SD\n");
for(q=0;q<num_parts_use;q++)
{
if(parttype==0&&q<num_parts_use-1){fprintf(output2,"Cat_A%d ",q+1);}
if(parttype==0&&q==num_parts_use-1){fprintf(output2,"Cat_Base ");}
if(parttype==1){fprintf(output2,"Cat_P%d ",q+1);}
fprintf(output2, "%.6f %.6f\n", stats[total+1+q], stats[total+1+q+total+1+num_parts_use]);
}
fclose(output2);

//shares (at bottom, code for trying to work out influences)
if(divide==-9999){sprintf(filename3,"%s.share",outfile);}
else{sprintf(filename3,"%s.%d.share",outfile,m+1);}
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3,"Component Share SD\n");
for(q=0;q<num_parts_use;q++)
{
if(parttype==0&&q<num_parts_use-1){fprintf(output3,"Share_A%d ",q+1);}
if(parttype==0&&q==num_parts_use-1){fprintf(output3,"Share_Base ");}
if(parttype==1){fprintf(output3,"Share_P%d ",q+1);}
fprintf(output3, "%.6f %.6f\n", stats[q]/stats[total], stats[q+2*(total+1+num_parts_use)]);
}
fclose(output3);

if(divide==-9999){sprintf(filename4,"%s.enrich",outfile);}
else{sprintf(filename4,"%s.%d.enrich",outfile,m+1);}
if((output4=fopen(filename4,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename4);exit(1);}
fprintf(output4,"Component Share SD Expected Enrichment SD Z-Stat1 Z-Stat2\n");
for(q=0;q<num_parts_use;q++)
{
//compute share and get sd
value=stats[total+1+q]/stats[total];
value2=stats[total+1+q+2*(total+1+num_parts_use)];
//print category names
if(parttype==0&&q<num_parts_use-1){fprintf(output4,"Enrich_A%d ",q+1);}
if(parttype==0&&q==num_parts_use-1){fprintf(output4,"Enrich_Base ");}
if(parttype==1){fprintf(output4,"Enrich_P%d ",q+1);}
//print share, sd, expected, enrichment and sd
fprintf(output4,"%.6f %.6f %.6f %.6f %.6f ", value, value2, ssums_use[q][num_parts_use+1], value/ssums_use[q][num_parts_use+1], value2/ssums_use[q][num_parts_use+1]);
//print test statistics
if(value2>0){fprintf(output4, "%6f ", (value-ssums_use[q][num_parts_use+1])/value2);}
else{fprintf(output4, "NA ");}
if(cohers2[q]>0){fprintf(output4, "%6f\n", stats[total]*(value-ssums_use[q][num_parts_use+1])/cohers2[q]);}
else{fprintf(output4, "NA\n");}
}
fclose(output4);

if(prev!=-9999)
{
if(divide==-9999){sprintf(filename5,"%s.hers.liab",outfile);}
else{sprintf(filename5,"%s.%d.hers.liab",outfile,m+1);}
if((output5=fopen(filename5,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename5);exit(1);}
fprintf(output5,"Component Heritability SD Influence SD\n");
for(q=0;q<num_parts_use;q++)
{
if(parttype==0&&q<num_parts_use-1){fprintf(output5,"Her_A%d ",q+1);}
if(parttype==0&&q==num_parts_use-1){fprintf(output5,"Her_Base ");}
if(parttype==1){fprintf(output5,"Her_P%d ",q+1);}
fprintf(output5,"%.6f %.6f %.6f %.6f\n", stats[q]*factor, stats[q+total+1+num_parts_use]*factor, influs[q]*stats[q], fabs(influs[q])*stats[q+total+1+num_parts_use]);
}
if(sum2>0){fprintf(output5,"Her_All %.6f %.6f %.6f %.6f\n", stats[total]*factor, stats[total+total+1+num_parts_use]*factor, sum, pow(sum2,.5));}
else{fprintf(output5,"Her_All %.6f %.6f %.6f NA\n", stats[total]*factor, stats[total+total+1+num_parts_use]*factor, sum);}
fclose(output5);

if(divide==-9999){sprintf(filename6,"%s.cats.liab",outfile);}
else{sprintf(filename6,"%s.%d.cats.liab",outfile,m+1);}
if((output6=fopen(filename6,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename6);exit(1);}
fprintf(output6,"Component Heritability SD\n");
for(q=0;q<num_parts_use;q++)
{
if(parttype==0&&q<num_parts_use-1){fprintf(output6,"Cat_A%d ",q+1);}
if(parttype==0&&q==num_parts_use-1){fprintf(output6,"Cat_Base ");}
if(parttype==1){fprintf(output6,"Cat_P%d ",q+1);}
fprintf(output6, "%.6f %.6f\n", stats[total+1+q]*factor, stats[total+1+q+total+1+num_parts_use]*factor);
}
fclose(output6);
}

///////

//for extra, want effective size and mean (already have weighted gif)
sum=0;for(j=0;j<count;j++){sum+=pow(stags[j],-1);}
sum2=0;for(j=0;j<count;j++){sum2+=snss[j]/stags[j];}
value=sum2/sum;
sum2=0;for(j=0;j<count;j++){sum2+=schis[j]/stags[j];}
mean=sum2/sum;

if(divide==-9999){sprintf(filename7,"%s.extra",outfile);}
else{sprintf(filename7,"%s.%d.extra",outfile,m+1);}
if((output7=fopen(filename7,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename7);exit(1);}
fprintf(output7,"Weighted_Sample_Size %.2f\n", value);
fprintf(output7,"Num_Datapoints %d\nEffective_Size %.2f\n", count, sum);
fprintf(output7,"Weighted_Mean_Chisq %.4f\nWeighted_GIF %.4f\n", mean, gif);
fprintf(output7,"Null_logl %.6f\nAlt_logl %.6f\n", likes[0], likes[1]);

if(gcon==1){fprintf(output7,"Scaling_Estimate %.6f\nScaling_SD %.6f\n", stats[num_parts_use], stats[num_parts_use+total+1+num_parts_use]);}
if(cept==1){fprintf(output7,"Intercept_Estimate %.6f\nIntercept_SD %.6f\n", 1+stats[num_parts_use+gcon], stats[num_parts_use+gcon+total+1+num_parts_use]);}
if(likes[6]==1){fprintf(output7,"Converged YES\n");}
else{fprintf(output7,"Converged NO\n");}
fprintf(output7,"Num_Neg_Exp_Her %d\nEffective_Neg_Exp_Her %.2f\n", (int)likes[7], likes[8]);
fprintf(output7,"Num_Non_Pos_Stat %d\nEffective_Non_Pos_Stat %.2f\n", (int)likes[9], likes[10]);
fclose(output7);

if(divide==-9999){sprintf(filename8,"%s.cross",outfile);}
else{sprintf(filename8,"%s.%d.cross",outfile,m+1);}
if((output8=fopen(filename8,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename8);exit(1);}
fprintf(output8, "Component Estimate");
for(q=0;q<num_parts_use;q++){fprintf(output8, " Covar_%d", q+1);}
fprintf(output8, "\n");
for(q=0;q<num_parts_use;q++)
{
if(parttype==0&&q<num_parts_use-1){fprintf(output8,"Her_A%d ",q+1);}
if(parttype==0&&q==num_parts_use-1){fprintf(output8,"Her_Base ");}
if(parttype==1){fprintf(output8,"Her_P%d ",q+1);}
fprintf(output8, "%.6f", stats[q]);
for(q2=0;q2<num_parts_use;q2++)
{fprintf(output8, " %.6f", cohers[q+q2*num_parts_use]);}
fprintf(output8, "\n");
}
fclose(output8);

if(divide==-9999){sprintf(filename9,"%s.taus",outfile);}
else{sprintf(filename9,"%s.%d.taus",outfile,m+1);}
if((output9=fopen(filename9,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename9);exit(1);}
fprintf(output9, "Component Estimate");
for(q=0;q<num_parts_use;q++){fprintf(output9, " Covar_%d", q+1);}
fprintf(output9, "\n");

for(q=0;q<num_parts_use;q++)
{
if(parttype==0&&q<num_parts_use-1){fprintf(output9,"Tau_A%d ",q+1);}
if(parttype==0&&q==num_parts_use-1){fprintf(output9,"Tau_Base ");}
if(parttype==1){fprintf(output9,"Tau_P%d ",q+1);}
fprintf(output9, "%.6e", stats[q]/ssums_use[q][q]);
for(q2=0;q2<num_parts_use;q2++)
{fprintf(output9, " %.6e", cohers[q+q2*num_parts_use]/ssums_use[q][q]/ssums_use[q2][q2]);}
fprintf(output9, "\n");
}
fclose(output9);

if(divide==-9999){sprintf(filename10,"%s.labels",outfile);}
else{sprintf(filename10,"%s.%d.labels",outfile,m+1);}
if((output10=fopen(filename10,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename10);exit(1);}
fprintf(output10,"Category Label\n");
if(divide==-9999)
{
for(q=0;q<num_parts_use;q++)
{
if(parttype==0&&q<num_parts_use-1){fprintf(output10,"A%d %s\n",q+1, catlabels[q]);}
if(parttype==0&&q==num_parts_use-1){fprintf(output10,"Base %s", catlabels[q]);}
if(parttype==1){fprintf(output10,"P%d %s\n",q+1, catlabels[q]);}
}
}
else	//must have partitions (currently not possible to have labels)
{
for(q=0;q<num_parts_use;q++){fprintf(output10,"P%d %s\n",q+1, catlabels[q+m*divide]);}
}
fclose(output10);

printf("Main results saved in %s, %s, %s, %s and %s", filename, filename2, filename3, filename4, filename7);
if(prev!=-9999){printf(", with liability versions saved in %s and %s", filename5, filename6);}
printf("\n\n");
}	//end of saving most things
}	//end of m loop

if(strcmp(powfile,"blank")!=0)	//compute parameters
{
printf("Using a grid search to estimate mean and variance of power parameter\n");

//will try m from -1.25 to .75 (201) and s from .01 to .5 (50)
num_means=201;
mvals=malloc(sizeof(double)*num_means);
for(j=0;j<num_means;j++){mvals[j]=-1.25+.01*j;}
num_sds=50;
svals=malloc(sizeof(double)*num_sds);
for(k=0;k<num_sds;k++){svals[k]=.01+.01*k;}

perfs=malloc(sizeof(double*)*num_means);
for(j=0;j<num_means;j++){perfs[j]=malloc(sizeof(double)*num_sds);}

//find max powlikes, then subtract max and get exponent
value=powlikes[5*num_anals];
for(m=1;m<num_anals;m++)
{
if(powlikes[m+5*num_anals]>value){value=powlikes[m+5*num_anals];}
}
for(m=0;m<num_anals;m++){powlikes[m]=exp(powlikes[m+5*num_anals]-value);}

//loop through all pairs of means and sds
for(j=0;j<num_means;j++)
{
//need to subtract something to avoid massive values
value=fabs(powers[0]-mvals[j]);
for(m=1;m<num_anals;m++)
{
if(fabs(powers[m]-mvals[j])<value){value=fabs(powers[m]-mvals[j]);}
}
for(k=0;k<num_sds;k++)
{
//get sums and sumsqs for all data points
sum=0;sum2=0;sumsq=0;sumsq2=0;sumsq3=0;
for(m=0;m<num_anals;m++)
{
value2=exp(-.5*(pow(powers[m]-mvals[j],2)-pow(value,2))*pow(svals[k],-2));
sum+=powlikes[m];sum2+=value2;sumsq+=pow(powlikes[m],2);sumsq2+=pow(value2,2);sumsq3+=powlikes[m]*value2;
}

//now get max correlation with one value removed
perfs[j][k]=-2;
for(m=0;m<num_anals;m++)
{
value2=exp(-.5*(pow(powers[m]-mvals[j],2)-pow(value,2))*pow(svals[k],-2));
sum-=powlikes[m];sum2-=value2;sumsq-=pow(powlikes[m],2);sumsq2-=pow(value2,2);sumsq3-=powlikes[m]*value2;
value3=((num_anals-1)*sumsq3-sum*sum2)*pow((num_anals-1)*sumsq-sum*sum,-.5)*pow((num_anals-1)*sumsq2-sum2*sum2,-.5);
if(isinf(value3)!=1&&value3>perfs[j][k]){perfs[j][k]=value3;}
sum+=powlikes[m];sum2+=value2;sumsq+=pow(powlikes[m],2);sumsq2+=pow(value2,2);sumsq3+=powlikes[m]*value2;
}
}}	//end of j and k loops

//get max
topm=0;tops=0;value=perfs[0][0];
for(j=0;j<num_means;j++)
{
for(k=0;k<num_sds;k++)
{
if(perfs[j][k]>value){topm=j;tops=k;value=perfs[j][k];}
}}
printf("Best power is %.4f (SD %.4f)\n\n", mvals[topm], svals[tops]);

//column two of powlikes is expected value (ensure same sum as column one)
value=fabs(powers[0]-mvals[topm]);
for(m=1;m<num_anals;m++)
{
if(fabs(powers[m]-mvals[topm])<value){value=fabs(powers[m]-mvals[topm]);}
}
for(m=0;m<num_anals;m++)
{powlikes[num_anals+m]=value2=exp(-.5*(pow(powers[m]-mvals[topm],2)-pow(value,2))*pow(svals[tops],-2));}
sum=0;sum2=0;
for(m=0;m<num_anals;m++){sum+=powlikes[m];sum2+=powlikes[num_anals+m];}
for(m=0;m<num_anals;m++){powlikes[num_anals+m]*=sum/sum2;}

//also get mode
found=0;
for(m=1;m<num_anals;m++)
{
if(powlikes[m]>powlikes[found]){found=m;}
}

//save
sprintf(filename,"%s.power",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fprintf(output,"Component Mode Mean SD Correlation\nPower %.4f %.4f %.4f %.4f\n", powers[found], mvals[topm], svals[tops], perfs[topm][tops]);
fclose(output);

sprintf(filename2,"%s.power.cors",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2,"SD_Mean");
for(j=0;j<num_means;j++){fprintf(output2," %.3f", mvals[j]);}
fprintf(output2,"\n");
for(k=0;k<num_sds;k++)
{
fprintf(output2,"%.3f", svals[k]);
for(j=0;j<num_means;j++){fprintf(output2," %.3f", perfs[j][k]);}
fprintf(output2,"\n");
}
fclose(output2);

sprintf(filename3,"%s.power.plot",outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3,"Power Her_All SD Null_Chisq_Likelihood Alt_Chisq_Likelihood Observed Expected\n");
for(m=0;m<num_anals;m++){fprintf(output3,"%.4f %.6f %.6f %.6f %.6f %.6f %.6f\n", powers[m], powlikes[m+2*num_anals], powlikes[m+3*num_anals], powlikes[m+4*num_anals], powlikes[m+5*num_anals], powlikes[m], powlikes[m+num_anals]);}
fclose(output3);

free(mvals);free(svals);
for(j=0;j<num_means;j++){free(perfs[j]);}free(perfs);
}	//end of have powfile

for(q=0;q<num_parts;q++){free(catlabels[q]);}free(catlabels);
for(j=0;j<count;j++){free(spreds[j]);}free(spreds);free(sal1);free(sal2);free(stags);
for(q=0;q<num_parts;q++){free(svars[q]);free(ssums[q]);}free(svars);free(ssums);
free(snss);free(schis);free(srhos);
if(strcmp(cvsfile,"blank")!=0){free(cvindex);free(cvexps);}
free(stats);free(likes);free(cohers);free(cohers2);free(influs);

if(strcmp(catfile,"blank")!=0)
{
for(q=0;q<num_reds;q++){free(svars2[q]);free(ssums2[q]);}free(svars2);free(ssums2);free(stats2);
}
if(divide!=-9999)
{
for(q=0;q<divide;q++){free(svars2[q]);free(ssums2[q]);}free(svars2);free(ssums2);
}
if(strcmp(powfile,"blank")!=0){free(powers);free(powlikes);}

///////////////////////////

//code for multiplicative solver 
/*
else	//multiplicative solver - must have parttype=1, can not be using catfile or taufile
{
//solve for first num_bases categories + will not be jackknifing
for(q=0;q<num_bases;q++)
{
for(j=0;j<count;j++){svars2[q][j]=svars[q][j];}
for(q2=0;q2<num_bases;q2++){ssums2[q][q2]=ssums[q][q2];}
ssums2[q][num_bases]=ssums[q][num_parts];
ssums2[q][num_bases+1]=ssums[q][num_parts+1];
}

printf("Solving for first %d categories\n", num_bases);
solve_sums(stats2, NULL, NULL, num_bases, gcon, cept, -9999, count, ncv, cvindex, NULL, stags, svars2, ssums2, snss, schis, 1, tol, sumiter, chisol, 4, ldsc, filename);

//get weightings (equal to taus)
for(k=0;k<num_bases;k++){sweights[k]=stats2[k]/ssums2[k][num_bases];}
sum=0;for(k=0;k<num_bases;k++){sum+=sweights[k];}
if(sum<0){printf("Warning, the required multiplier is negative (%f), which may cause problems; please tell Doug\n\n", sum);}

//now collapse svars into svars3 and ssums into ssums3
for(q=0;q<num_parts_use;q++)
{
for(j=0;j<count;j++)
{
svars3[q][j]=0;
for(k=0;k<num_bases;k++){svars3[q][j]+=svars[k+q*num_bases][j]*sweights[k];}
}

for(q2=0;q2<num_parts_use;q2++)
{
ssums3[q][q2]=0;
for(k=0;k<num_bases;k++){ssums3[q][q2]+=ssums[k+q*num_bases][k+q2*num_bases]*sweights[k];}
}
ssums3[q][num_parts_use]=0;
for(k=0;k<num_bases;k++){ssums3[q][num_parts_use]+=ssums[k+q*num_bases][num_parts]*sweights[k];}

//using partitions, so last value of ssums3[q] is sum of first num_parts_use values
ssums3[q][num_parts_use+1]=0;
for(q2=0;q2<num_parts_use;q2++){ssums3[q][num_parts_use+1]+=ssums3[q][q2];}
value=0;
for(k=0;k<num_bases;k++){value+=ssums[k+q*num_bases][num_parts+1]*sweights[k];}
}

//load coefficients into stats and solve for multiplicative model
for(q=0;q<num_parts_use;q++){stats[q]=0;}
for(k=0;k<num_bases;k++){stats[0]+=stats2[k];}
if(gcon==1){stats[num_parts_use]=stats2[num_bases];}
if(cept==1){stats[num_parts_use+gcon]=stats2[num_bases+gcon];}

printf("Now including remaining %d categories multiplicatively\n", num_parts-num_bases);
solve_sums(stats, likes, cohers, num_parts_use, gcon, cept, num_blocks, count, ncv, cvindex, cvexps, stags, svars3, ssums3, snss, schis, 1, tol, sumiter, chisol, 5, ldsc, filename);

//now recover taus (and covariances) for each original category
for(q=0;q<num_parts_use;q++)
{
for(k=0;k<num_bases;k++)
{
stats3[k+q*num_bases]=stats[q]/ssums3[q][num_parts_use]*sweights[k];
for(q2=0;q2<num_parts_use;q2++)
{
for(k2=0;k2<num_bases;k2++){stats3[num_parts+k+q*num_bases+(k2+q2*num_bases)*num_parts]=cohers[q+q2*num_parts_use]/ssums3[q][num_parts_use]/ssums3[q2][num_parts_use]*sweights[k]*sweights[k2];}
}
}}
}	//end of num_bases>1

//////////////////////////

//code for old influences - equal to intensity/av_intensity minus its expected value
if(divide==-9999){sprintf(filename3,"%s.share",outfile);}
else{sprintf(filename3,"%s.%d.share",outfile,m+1);}
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3,"Component Share SD\n");
for(q=0;q<num_parts_use;q++)
{
//value is divisor (her_base or her_all) and value 2 is expected share (so expected influence is I(value2>0)) 
//if(parttype==0&&q<num_parts_use-1){fprintf(output3,"Share_A%d ",q+1);value=stats[num_parts_use-1];value2=0;}
//if(parttype==0&&q==num_parts_use-1){fprintf(output3,"Share_Base ");value=stats[num_parts_use-1];value2=1;}
//if(parttype==1){fprintf(output3,"Share_P%d ",q+1);value=stats[total];value2=ssums_use[q][num_parts_use]/sum;}
//fprintf(output3, "%.6f %.6f %.6f %.6f %.6f\n", stats[q]/value, stats[q+2*(total+1+num_parts_use)], value2, stats[q]/value/ssums_use[q][num_parts_use]*sum-(value2>0), stats[q+2*(total+1+num_parts_use)]/ssums_use[q][num_parts_use]*sum);

//value is expected share, value2 is E(intensity/av_intensity)
if(parttype==0&&q<num_parts_use-1){fprintf(output3,"Share_A%d ",q+1);value=0;value2=0;}
if(parttype==0&&q==num_parts_use-1){fprintf(output3,"Share_Base ");value=1;value2=1;}
if(parttype==1){fprintf(output3,"Share_P%d ",q+1);value=ssums_use[q][num_parts_use]/sum;value2=1;}
fprintf(output3, "%.6f %.6f %.6f %.6f %.6f\n", stats[q]/stats[total], stats[q+2*(total+1+num_parts_use)], value, stats[q]/stats[total]/ssums_use[q][num_parts_use]*sum-value2, stats[q+2*(total+1+num_parts_use)]/ssums_use[q][num_parts_use]*sum);
}
fclose(output3);

*/

