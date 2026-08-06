/*
Copyright 2026 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Summary correlations - if using xagfile, have set num_parts=3

///////////////////////////

//read in tag details
if(strcmp(tagfile,"blank")!=0)
{
count=countrows(tagfile)-2-num_parts;
model_warn(count, num_parts);
}
else
{
count=countrows(xagfile)-5;
model_warn(count, 1);
}

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
if(strcmp(tagfile,"blank")!=0){count=read_tagfile(tagfile, catlabels, spreds, sal1, sal2, stags, svars, ssums, num_parts, count, parttype, altfile, bpredfile, cpredfile, tagone);}
else{count=read_tagfile(xagfile, catlabels, spreds, sal1, sal2, stags, svars, ssums, num_parts, count, parttype, altfile, bpredfile, cpredfile, tagone);}

//get summary stats
snss=malloc(sizeof(double)*count);
snss2=malloc(sizeof(double)*count);
schis=malloc(sizeof(double)*count);
schis2=malloc(sizeof(double)*count);
srhos=malloc(sizeof(double)*count);
srhos2=malloc(sizeof(double)*count);

printf("Reading summary statistics for Trait 1 from %s\n", sumsfile);
if(strcmp(tagfile,"blank")!=0){read_sumsfile(sumsfile, snss, schis, srhos, NULL, count, spreds, sal1, sal2, tagfile, amb, scaling, checksums);}
else{read_sumsfile(sumsfile, snss, schis, srhos, NULL, count, spreds, sal1, sal2, xagfile, amb, scaling, checksums);}
printf("First few stats and ns are: %s %.3f %.1f", spreds[0], schis[0], snss[0]);
for(j=1;j<count;j++){if(j<3){printf(" | %s %.3f %.1f", spreds[j], schis[j], snss[j]);}}
printf("\n\n");

printf("Reading summary statistics for Trait 2 from %s\n", sums2file);
if(strcmp(tagfile,"blank")!=0){read_sumsfile(sums2file, snss2, schis2, srhos2, NULL, count, spreds, sal1, sal2, tagfile, amb, scaling2, checksums);}
else{read_sumsfile(sums2file, snss2, schis2, srhos2, NULL, count, spreds, sal1, sal2, xagfile, amb, scaling2, checksums);}
printf("First few stats and ns are: %s %.3f %.1f", spreds[0], schis2[0], snss2[0]);
for(j=1;j<count;j++){if(j<3){printf(" | %s %.3f %.1f", spreds[j], schis2[j], snss2[j]);}}
printf("\n\n");

//print out overlap
count2=0;for(j=0;j<count;j++){count2+=(snss[j]>0&&snss2[j]>0);}
sprintf(filename,"%s.overlap",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fprintf(output,"Tagging_File_Predictors %d\nSummary_Statistic_Predictors %d\nOverlap_Proportion %.4f\n", count, count2, (double)count2/count);
fclose(output);

if(count2<count){printf("In total, %d of the %d predictors have summary statistics for both traits\n\n", count2, count);} 

if(cutoff!=-9999)	//look for big predictors
{
count2=0;
for(j=0;j<count;j++)
{
if(snss[j]>0&&pow(srhos[j],2)>=cutoff){snss[j]=0;count2++;}
}
if(count2==0){printf("There are no predictors that explain at least %.6f of phenotypic variance for Trait 1\n\n", cutoff);}
if(count2==1){printf("There is 1 predictor that explains at least %.6f of phenotypic variance for Trait 1; this will be excluded\n\n", cutoff);}
if(count2>1){printf("There are %d predictors that explain at least %.6f of phenotypic variance for Trait 1; these will be excluded\n\n", count2, cutoff);}

count2=0;
for(j=0;j<count;j++)
{
if(snss2[j]>0&&pow(srhos2[j],2)>=cutoff){snss2[j]=0;count2++;}
}
if(count2==0){printf("There are no predictors that explain at least %.6f of phenotypic variance for Trait 2\n\n", cutoff);}
if(count2==1){printf("There is 1 predictor that explains at least %.6f of phenotypic variance for Trait 2; this will be excluded\n\n", cutoff);}
if(count2>1){printf("There are %d predictors that explain at least %.6f of phenotypic variance for Trait 2; these will be excluded\n\n", count2, cutoff);}
}

if(cutoff2!=-9999)	//truncate big predictors
{
count2=0;
for(j=0;j<count;j++)
{
if(snss[j]>0&&schis[j]>=cutoff2)
{
schis[j]=cutoff2;
if(srhos[j]>0){srhos[j]=pow(cutoff2/(cutoff2+snss[j]),.5);}
else{srhos[j]=-pow(cutoff2/(cutoff2+snss[j]),.5);}
count2++;
}
}
if(count2==0){printf("There are no predictors with chi-squared statistic over %.2f for Trait 1\n\n", cutoff2);}
if(count2==1){printf("There is 1 predictor with chi-squared statistic over %.2f for Trait 1; this will be truncated\n\n", cutoff2);}
if(count2>1){printf("There are %d predictors with chi-squared statistic over %.2f for Trait 1; these will be truncated\n\n", count2, cutoff2);}

count2=0;
for(j=0;j<count;j++)
{
if(snss2[j]>0&&schis2[j]>=cutoff2)
{
schis2[j]=cutoff2;
if(srhos2[j]>0){srhos2[j]=pow(cutoff2/(cutoff2+snss2[j]),.5);}
else{srhos2[j]=-pow(cutoff2/(cutoff2+snss2[j]),.5);}
count2++;
}
}
if(count2==0){printf("There are no predictors with chi-squared statistic over %.2f for Trait 2\n\n", cutoff2);}
if(count2==1){printf("There is 1 predictor with chi-squared statistic over %.2f for Trait 2; this will be truncated\n\n", cutoff2);}
if(count2>1){printf("There are %d predictors with chi-squared statistic over %.2f for Trait 2; these will be truncated\n\n", count2, cutoff2);}
}

//count number of predictors we will use
count2=0;for(j=0;j<count;j++){count2+=(snss[j]>0&&snss2[j]>0);}
if(count2==0){printf("Error, there are no predictors\n\n");exit(1);}
//if(count2<3){printf("Error, there are only %d predictors, can not continue\n\n", count2);exit(1);}
if(count2<20000){printf("Warning, there are usually tens of thousands of predictors (not %d)\n\n", count2);}

if(count2<count)	//squeeze down and reset count
{
count2=0;
for(j=0;j<count;j++)
{
if(snss[j]>0&&snss2[j]>0)
{
if(count2!=j)
{
free(spreds[count2]);copy_string(spreds,count2,spreds[j]);
sal1[count2]=sal1[j];
sal2[count2]=sal2[j];
stags[count2]=stags[j];
for(q=0;q<num_parts;q++){svars[q][count2]=svars[q][j];}
snss[count2]=snss[j];snss2[count2]=snss2[j];
schis[count2]=schis[j];schis2[count2]=schis2[j];
srhos[count2]=srhos[j];srhos2[count2]=srhos2[j];
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
if(schis2[j]==0){count2++;schis2[j]=1e-6;}
}
if(count2>0){printf("Warning, %d SNPs have test statistic 0; these have been replaced with 1e-6\n\n", count2);}

(void)get_sum_stats(stags, snss, schis, srhos, count, 1, cutoff);
(void)get_sum_stats(stags, snss2, schis2, srhos2, count, 2, cutoff);

if(strcmp(xagfile,"blank")!=0)	//copy svars into svars2, and get reduced version of ssums
{
svars2=malloc(sizeof(double*)*num_parts);
svars2[0]=malloc(sizeof(double)*count);
svars2[1]=malloc(sizeof(double)*count);
svars2[2]=malloc(sizeof(double)*count);
ssums2=malloc(sizeof(double*));
ssums2[0]=malloc(sizeof(double)*3);

for(j=0;j<count;j++){svars2[0][j]=svars[0][j];}
for(j=0;j<count;j++){svars2[1][j]=svars[1][j];}
for(j=0;j<count;j++){svars2[2][j]=svars[2][j];}
ssums2[0][0]=ssums[0][0];
ssums2[0][1]=ssums[0][0];
ssums2[0][2]=1;
}

////////

//set totals and allocate
if(strcmp(tagfile,"blank")!=0)
{
total=2*(num_parts+gcon+cept)+num_parts+oversamp;
total2=num_parts+gcon+cept;
total3=total+4+num_parts;
stats=malloc(sizeof(double)*total3*2);
stats2=malloc(sizeof(double)*4);
}
else    //have cross-tagging - will only have one category
{
total=2*(1+gcon+cept)+1+oversamp;
total2=1+gcon+cept;
total3=total+5;
stats=malloc(sizeof(double)*total3*2);
stats2=malloc(sizeof(double)*4);
}

if(strcmp(tagfile,"blank")!=0)
{
sprintf(filename,"%s.progress",outfile);
stats2[0]=prev;stats2[1]=prev2;stats2[2]=ascer;stats2[3]=ascer2;
solve_cors(stats, num_parts, gcon, cept, oversamp, num_blocks, count, stags, svars, ssums, snss, schis, srhos, snss2, schis2, srhos2, tol, maxiter, filename, 0, stats2);

//save
sprintf(filename,"%s.cors",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fprintf(output,"Component Value SE\n");
fprintf(output,"Her1_All %.6f %.6f\n", stats[total], stats[total+total3]);
fprintf(output,"Her2_All %.6f %.6f\n", stats[total+1], stats[total+1+total3]);
fprintf(output,"Coher_All %.6f %.6f\n", stats[total+2], stats[total+2+total3]);
fprintf(output,"Cor_All %.6f %.6f\n", stats[total+3], stats[total+3+total3]);
if(gcon==1)
{
fprintf(output,"Scaling1 %.6f %.6f\n", stats[num_parts], stats[num_parts+total3]);
fprintf(output,"Scaling2 %.6f %.6f\n", stats[total2+num_parts], stats[total2+num_parts+total3]);
}
if(cept==1)
{
fprintf(output,"Intercept1 %.6f %.6f\n", stats[num_parts+gcon], stats[num_parts+gcon+total3]);
fprintf(output,"Intercept2 %.6f %.6f\n", stats[total2+num_parts+gcon], stats[total2+num_parts+gcon+total3]);
}
if(oversamp==1){fprintf(output,"Overlap %.6f %.6f\n", stats[2*total2+num_parts], stats[2*total2+num_parts+total3]);}
else{fprintf(output,"Overlap 0 0\n");}
fclose(output);

sprintf(filename2,"%s.cors.full",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2,"Category Trait1_Her SE Trait2_Her SE Both_Coher SE Correlation SE\n");
for(q=0;q<num_parts;q++)
{
if(parttype==0&&q<num_parts-1){fprintf(output2,"A%d ",q+1);}
if(parttype==0&&q==num_parts-1){fprintf(output2,"Base ");}
if(parttype==1){fprintf(output2,"P%d ",q+1);}
fprintf(output2,"%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n", stats[q], stats[q+total3], stats[total2+q], stats[total2+q+total3], stats[2*total2+q], stats[2*total2+q+total3], stats[total+4+q], stats[total+4+q+total3]);
}
fprintf(output2,"ALL %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n", stats[total], stats[total+total3], stats[total+1], stats[total+1+total3], stats[total+2], stats[total+2+total3], stats[total+3], stats[total+3+total3]);
fclose(output2);

sprintf(filename3,"%s.labels",outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3,"Category Label\n");
for(q=0;q<num_parts;q++)
{
if(parttype==0&&q<num_parts-1){fprintf(output3,"A%d %s\n",q+1, catlabels[q]);}
if(parttype==0&&q==num_parts-1){fprintf(output3,"Base %s\n", catlabels[q]);}
if(parttype==1){fprintf(output3,"P%d %s\n",q+1, catlabels[q]);}
}
fclose(output3);

if(prev!=-9999)
{
factor=get_factor(NULL, -9999, prev, ascer, NULL);
factor2=get_factor(NULL, -9999, prev2, ascer2, NULL);

sprintf(filename4,"%s.cors.liab",outfile);
if((output4=fopen(filename4,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename4);exit(1);}
fprintf(output4,"Component Value SE\n");
fprintf(output4,"Her1_All %.6f %.6f\n", stats[total]*factor, stats[total+total3]*factor);
fprintf(output4,"Her2_All %.6f %.6f\n", stats[total+1]*factor2, stats[total+1+total3]*factor2);
fprintf(output4,"Coher_All %.6f %.6f\n", stats[total+2]*pow(factor*factor2,.5), stats[total+2+total3]*pow(factor*factor2,.5));
fprintf(output4,"Cor_All %.6f %.6f\n", stats[total+3], stats[total+3+total3]);

fprintf(output4,"Genetic_Separation %.6f %.6f\n", stats2[0], stats2[1]);
fprintf(output4,"Max_AUC %.6f %.6f\n", stats2[2], stats2[3]);
fclose(output4);
}

if(prev==-9999){printf("Results saved in %s and %s\n\n", filename, filename2);}
else{printf("Results saved in %s, %s and %s\n\n", filename, filename2, filename4);}
}
else	//cross tagging
{
sprintf(filename,"%s.progress",outfile);
stats2[0]=prev;stats2[1]=prev2;stats2[3]=ascer;stats2[4]=ascer2;
solve_cors(stats, 1, gcon, cept, oversamp, num_blocks, count, stags, svars2, ssums2, snss, schis, srhos, snss2, schis2, srhos2, tol, maxiter, filename, 1, NULL);

//save
sprintf(filename,"%s.cors",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fprintf(output,"Component Value SE\n");
fprintf(output,"Her1_All %.6f %.6f\n", stats[total], stats[total+total3]);
fprintf(output,"Her2_All %.6f %.6f\n", stats[total+1], stats[total+1+total3]);
fprintf(output,"Coher_All %.6f %.6f\n", stats[total+2], stats[total+2+total3]);
fprintf(output,"Cor_All %.6f %.6f\n", stats[total+3], stats[total+3+total3]);
if(gcon==1)
{
fprintf(output,"Scaling1 %.6f %.6f\n", stats[1], stats[1+total3]);
fprintf(output,"Scaling2 %.6f %.6f\n", stats[total2+1], stats[total2+1+total3]);
}
if(cept==1)
{
fprintf(output,"Intercept1 %.6f %.6f\n", stats[1+gcon], stats[1+gcon+total3]);
fprintf(output,"Intercept2 %.6f %.6f\n", stats[total2+1+gcon], stats[total2+1+gcon+total3]);
}
if(oversamp==1){fprintf(output,"Overlap %.6f %.6f\n", stats[2*total2+1], stats[2*total2+1+total3]);}
else{fprintf(output,"Overlap 0 0\n");}
fclose(output);

printf("Results saved in %s\n\n", filename);
}

////////

for(q=0;q<num_parts;q++){free(catlabels[q]);}free(catlabels);
for(j=0;j<count;j++){free(spreds[j]);}free(spreds);free(sal1);free(sal2);free(stags);
for(q=0;q<num_parts;q++){free(svars[q]);free(ssums[q]);}free(svars);free(ssums);
free(snss);free(snss2);free(schis);free(schis2);free(srhos);free(srhos2);
if(strcmp(xagfile,"blank")!=0){free(svars2[0]);free(svars2[1]);free(svars2[2]);free(svars2);free(ssums2[0]);free(ssums2);}
free(stats);free(stats2);

///////////////////////////

