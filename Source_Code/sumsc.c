/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Pathways heritability

///////////////////////////

//read in base details

sprintf(filename,"%s.base.tagging", pathfile);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--calc-taggings\"\n\n", filename);exit(1);}
count4=countrows(filename)-1;

spreds=malloc(sizeof(char*)*count4);
sal1=malloc(sizeof(char)*count4);
sal2=malloc(sizeof(char)*count4);
stags=malloc(sizeof(double)*count4);
svars=malloc(sizeof(double*)*3);
ssums=malloc(sizeof(double*)*3);
for(q=0;q<3;q++)
{
svars[q]=malloc(sizeof(double)*count4);
ssums[q]=malloc(sizeof(double)*5);
}

//read tagfile - if extracting, count will be reset
count=read_tagfile_base(filename, spreds, sal1, sal2, stags, svars, count4, bpredfile, cpredfile);

//get summary stats
snss=malloc(sizeof(double)*count);
schis=malloc(sizeof(double)*count);
srhos=malloc(sizeof(double)*count);

printf("Reading summary statistics from %s\n", sumsfile);
read_sumsfile(sumsfile, snss, schis, srhos, NULL, count, spreds, sal1, sal2, tagfile, amb, fixn, scaling, 1, checksums, sformat);
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

//get indexes of those we will use, and squeeze down
indexer=malloc(sizeof(int)*count);
count2=0;
for(j=0;j<count;j++)
{
if(snss[j]!=-9999)
{
indexer[count2]=j;
if(count2!=j)
{
free(spreds[count2]);copy_string(spreds,count2,spreds[j]);
sal1[count2]=sal1[j];
sal2[count2]=sal2[j];
stags[count2]=stags[j];
svars[0][count2]=svars[0][j];
svars[1][count2]=svars[1][j];
snss[count2]=snss[j];
schis[count2]=schis[j];
srhos[count2]=srhos[j];
}
count2++;
}}
for(j=count2;j<count;j++){free(spreds[j]);}
count=count2;

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

if(dougvar==1)	//for testing, set genic taggings to zero
{
printf("Setting genic values to zero (for testing only)\n\n");
for(j=0;j<count;j++){svars[1][j]=0;}
}

////////

//allocate

stats=malloc(sizeof(double)*(3+gcon+cept+1+3)*3);
stats2=malloc(sizeof(double)*4);
likes=malloc(sizeof(double)*11);
cohers=malloc(sizeof(double)*2*2);
cohers2=malloc(sizeof(double)*3*3);
influs=malloc(sizeof(double)*2);
readfloats=malloc(sizeof(float)*count4);

//read ssums for base and genic - after pathway name, need elements 6, 7, 8, 10, 11, 12 of pathway.sums
sprintf(filename2,"%s.pathway.sums", pathfile);
if((input2=fopen(filename2,"r"))==NULL)
{printf("Error opening %s\n\n",filename2);exit(1);}
if(fscanf(input2, "%s %lf %lf %s %lf %lf %lf %s %lf ", readstring, ssums[0]+0, ssums[0]+1, readstring, ssums[0]+2, ssums[1]+0, ssums[1]+1, readstring, ssums[1]+2)!=9)
{printf("Error reading first row of %s\n\n", filename2);exit(1);}
fclose(input2);

//put proportions as last elements (probably not used)
ssums[0][3]=1;
ssums[1][3]=ssums[1][2]/ssums[0][2];

//test base model
num_parts_use=2;
total=2+gcon+cept;

//screen and file print
printf("Solving for base category and genic predictors\n\n");

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fprintf(output,"Solving for base category and genic predictors\n");
fclose(output);

solve_sums(stats, likes, cohers, NULL, influs, num_parts_use, gcon, cept, num_blocks, count, 0, NULL, NULL, stags, svars, ssums, snss, schis, tol, maxiter, chisol, 0, filename);

//get sum of influences, and its sd
sum=0;sum2=0;
for(q=0;q<num_parts_use;q++)
{
sum+=influs[q]*stats[q];
for(q2=0;q2<num_parts_use;q2++){sum2+=influs[q]*influs[q2]*cohers[q+q2*num_parts_use];}
}

sprintf(filename2,"%s.hers",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2,"Component Heritability SE Influence SE\n");
for(q=0;q<2;q++)
{
if(q==0){fprintf(output2,"Her_Base ");}
else{fprintf(output2,"Her_Genic ");}
fprintf(output2,"%.6f %.6f %.6f %.6f\n", stats[q], stats[q+total+1+num_parts_use], influs[q]*stats[q], fabs(influs[q])*stats[q+total+1+num_parts_use]);
}
if(sum2>0){fprintf(output2,"Her_All %.6f %.6f %.6f %.6f\n", stats[total], stats[total+total+1+num_parts_use], sum, pow(sum2,.5));}
else{fprintf(output2,"Her_All %.6f %.6f %.6f NA\n", stats[total], stats[total+total+1+num_parts_use], sum);}
fclose(output2);

sprintf(filename3,"%s.cats",outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3,"Component Heritability SE\n");
for(q=0;q<num_parts_use;q++)
{
if(q==0){fprintf(output3,"Cat_Base ");}
else{fprintf(output3,"Cat_Genic ");}
fprintf(output3, "%.6f %.6f\n", stats[total+1+q], stats[total+1+q+total+1+num_parts_use]);
}
fclose(output3);

//for extra, want effective size and mean (already have weighted gif)
sum=0;for(j=0;j<count;j++){sum+=pow(stags[j],-1);}
sum2=0;for(j=0;j<count;j++){sum2+=snss[j]/stags[j];}
value=sum2/sum;
sum2=0;for(j=0;j<count;j++){sum2+=schis[j]/stags[j];}
mean=sum2/sum;

sprintf(filename4,"%s.extra",outfile);
if((output4=fopen(filename4,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename4);exit(1);}
fprintf(output4,"Weighted_Sample_Size %.2f\n", value);
fprintf(output4,"Num_Datapoints %d\nEffective_Size %.2f\n", count, sum);
fprintf(output4,"Weighted_Mean_Chisq %.4f\nWeighted_GIF %.4f\n", mean, gif);
fprintf(output4,"Null_logl %.6f\nAlt_logl %.6f\n", likes[0], likes[1]);

if(gcon==1){fprintf(output4,"Scaling_Estimate %.6f\nScaling_SE %.6f\n", stats[num_parts_use], stats[num_parts_use+total+1+num_parts_use]);}
if(cept==1){fprintf(output4,"Intercept_Estimate %.6f\nIntercept_SE %.6f\n", 1+stats[num_parts_use+gcon], stats[num_parts_use+gcon+total+1+num_parts_use]);}
if(likes[6]==1){fprintf(output4,"Converged YES\n");}
else{fprintf(output4,"Converged NO\n");}
fprintf(output4,"Num_Neg_Exp_Her %d\nEffective_Neg_Exp_Her %.2f\n", (int)likes[7], likes[8]);
fprintf(output4,"Num_Non_Pos_Stat %d\nEffective_Non_Pos_Stat %.2f\n", (int)likes[9], likes[10]);
fclose(output4);

////////

//save coefficients to use as starting values
stats2[0]=stats[0];
stats2[1]=stats[1];
if(gcon==1){stats2[2]=stats[2];}
if(cept==1){stats2[2+gcon]=stats[2+gcon];}

//re-open progress
sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n",filename);exit(1);}

//re-open pathway.sums
sprintf(filename2,"%s.pathway.sums", pathfile);
if((input2=fopen(filename2,"r"))==NULL)
{printf("Error opening %s\n\n",filename2);exit(1);}

//open pathway.tagging
sprintf(filename3,"%s.pathway.tagging", pathfile);
if((input3=fopen(filename3,"r"))==NULL)
{printf("Error opening %s\n\n",filename3);exit(1);}

//open results file
sprintf(filename4,"%s.enrich",outfile);
if((output4=fopen(filename4,"w"))==NULL)
{printf("Error opening %s\n\n",filename4);exit(1);}
fprintf(output4,"Component Share SE Expected Enrichment SE Z-Stat1 Z-Stat2\n");

for(q=0;q<num_parts;q++)	//test annotation q
{
if(q%10==0)
{
printf("Testing pathway %d of %d\n", q+1, num_parts);

fclose(output);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n",filename);exit(1);}
fprintf(output, "Testing pathway %d of %d\n", q+1, num_parts);

fclose(output4);
sprintf(filename4,"%s.enrich",outfile);
if((output4=fopen(filename4,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename4);exit(1);}
}

//read ssums for base, genic and annotation q
if(fscanf(input2, "%s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf ", readstring, ssums[0]+0, ssums[0]+1, ssums[0]+2, ssums[0]+3, ssums[1]+0, ssums[1]+1, ssums[1]+2, ssums[1]+3, ssums[2]+0, ssums[2]+1, ssums[2]+2, ssums[2]+3)!=13)
{printf("Error reading Row %d of %s\n\n", q+1, filename2);exit(1);}

//put proportions as last elements
ssums[0][4]=1;
ssums[1][4]=ssums[1][3]/ssums[0][3];
ssums[2][4]=ssums[2][3]/ssums[0][3];

//read label
if(fscanf(input3, "%s ", readstring)!=1)
{printf("Error reading Row %d of %s\n\n", q+1, filename3);exit(1);}

//read taggings
for(j=0;j<count4;j++)
{
if(fscanf(input3, "%f ", readfloats+j)!=1)
{printf("Error reading Row %d of %s\n\n", q+1, filename3);exit(1);}
}
for(j=0;j<count;j++){svars[2][j]=readfloats[indexer[j]];}

//load coefficients into stats
stats[0]=stats2[0];
stats[1]=stats2[1];
stats[2]=0;
if(gcon==1){stats[3]=stats2[2];}
if(cept==1){stats[3+gcon]=stats2[2+gcon];}

//test model
num_parts_use=3;
total=3+gcon+cept;

solve_sums(stats, NULL, NULL, cohers2, NULL, num_parts_use, gcon, cept, num_blocks, count, 0, NULL, NULL, stags, svars, ssums, snss, schis, tol, maxiter, chisol, 6, NULL);

//compute share for annotation and get sd
q2=2;
value=stats[total+1+q2]/stats[total];
value2=stats[total+1+q2+2*(total+1+num_parts_use)];

fprintf(output4,"%s ",readstring);
if(value2!=-9999){fprintf(output4,"%.6f %.6f %.6f %.6f %.6f ", value, value2, ssums[q2][num_parts_use+1], value/ssums[q2][num_parts_use+1], value2/ssums[q2][num_parts_use+1]);}
else{fprintf(output4,"%.6f NA %.6f %.6f NA ", value, ssums[q2][num_parts_use+1], value/ssums[q2][num_parts_use+1]);}

//print test statistics
if(value2>0){fprintf(output4, "%6f ", (value-ssums[q2][num_parts_use+1])/value2);}
else{fprintf(output4, "NA ");}
if(cohers2[q2]>0){fprintf(output4, "%6f\n", stats[total]*(value-ssums[q2][num_parts_use+1])/cohers2[q2]);}
else{fprintf(output4, "NA\n");}
}	//end of q loop
printf("\n");
fclose(output);
fclose(input2);
fclose(input3);
fclose(output4);

for(j=0;j<count;j++){free(spreds[j]);}free(spreds);free(sal1);free(sal2);free(stags);
for(q=0;q<3;q++){free(svars[q]);free(ssums[q]);}free(svars);free(ssums);
free(snss);free(schis);free(srhos);
free(indexer);
free(stats);free(stats2);free(likes);free(cohers);free(cohers2);free(influs);free(readfloats);

///////////////////////////

