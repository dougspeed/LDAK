/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Summary correlations

///////////////////////////

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
snss2=malloc(sizeof(double)*count);
schis=malloc(sizeof(double)*count);
schis2=malloc(sizeof(double)*count);
srhos=malloc(sizeof(double)*count);
srhos2=malloc(sizeof(double)*count);

printf("Reading summary statistics for Trait 1 from %s\n", sumsfile);
read_sumsfile(sumsfile, snss, schis, srhos, count, spreds, sal1, sal2, tagfile, amb, fixn, scaling, 2, checksums);
printf("First few stats and ns are: %s %.3f %.1f", spreds[0], schis[0], snss[0]);
for(j=1;j<count;j++){if(j<3){printf(" | %s %.3f %.1f", spreds[j], schis[j], snss[j]);}}
printf("\n\n");

printf("Reading summary statistics for Trait 2 from %s\n", sums2file);
read_sumsfile(sums2file, snss2, schis2, srhos2, count, spreds, sal1, sal2, tagfile, amb, fixn2, scaling2, 2, checksums);
printf("First few stats and ns are: %s %.3f %.1f", spreds[0], schis2[0], snss2[0]);
for(j=1;j<count;j++){if(j<3){printf(" | %s %.3f %.1f", spreds[j], schis2[j], snss2[j]);}}
printf("\n\n");

//print out overlap
count2=0;for(j=0;j<count;j++){count2+=(snss[j]!=-9999&&snss2[j]!=-9999);}
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
if(count2==0){printf("There are no predictors that explain at least %.6f of phenotypic variance for Trait 1\n\n", cutoff);}
if(count2==1){printf("There is 1 predictor that explains at least %.6f of phenotypic variance for Trait 1; this will be excluded\n\n", cutoff);}
if(count2>1){printf("There are %d predictors that explain at least %.6f of phenotypic variance for Trait 1; these will be excluded\n\n", count2, cutoff);}

count2=0;
for(j=0;j<count;j++)
{
if(snss2[j]!=-9999&&pow(srhos2[j],2)>=cutoff){snss2[j]=-9999;count2++;}
}
if(count2==0){printf("There are no predictors that explain at least %.6f of phenotypic variance for Trait 2\n\n", cutoff);}
if(count2==1){printf("There is 1 predictor that explains at least %.6f of phenotypic variance for Trait 2; this will be excluded\n\n", cutoff);}
if(count2>1){printf("There are %d predictors that explain at least %.6f of phenotypic variance for Trait 2; these will be excluded\n\n", count2, cutoff);}
}

//count number of predictors we will use
count2=0;for(j=0;j<count;j++){count2+=(snss[j]!=-9999&&snss2[j]!=-9999);}
if(count2==0){printf("Error, there are no predictors\n\n");exit(1);}
//if(count2<3){printf("Error, there are only %d predictors, can not continue\n\n", count2);exit(1);}
if(count2<20000){printf("Warning, there are usually tens of thousands of predictors (not %d)\n\n", count2);}

if(count2<count)	//squeeze down and reset count
{
count2=0;
for(j=0;j<count;j++)
{
if(snss[j]!=-9999&&snss2[j]!=-9999)
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

////////

//set totals and allocate
total=2*(num_parts+gcon+cept)+num_parts+1;
total2=num_parts+gcon+cept;
stats=malloc(sizeof(double)*(total+4+num_parts)*2);

sprintf(filename,"%s.progress",outfile);
solve_cors(stats, num_parts, gcon, cept, num_blocks, count, stags, svars, ssums, snss, schis, srhos, snss2, schis2, srhos2, parttype, tol, maxiter, filename);

//save
sprintf(filename,"%s.cors",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fprintf(output,"Component Value SD\n");
fprintf(output,"Her1_All %.6f %.6f\n", stats[total], stats[total+total+4+num_parts]);
fprintf(output,"Her2_All %.6f %.6f\n", stats[total+1], stats[total+1+total+4+num_parts]);
fprintf(output,"Coher_All %.6f %.6f\n", stats[total+2], stats[total+2+total+4+num_parts]);
fprintf(output,"Cor_All %.6f %.6f\n", stats[total+3], stats[total+3+total+4+num_parts]);
if(gcon==1)
{
fprintf(output,"Scaling1 %.6f %.6f\n", stats[num_parts], stats[num_parts+total+4+num_parts]);
fprintf(output,"Scaling2 %.6f %.6f\n", stats[total2+num_parts], stats[total2+num_parts+total+4+num_parts]);
}
if(cept==1)
{
fprintf(output,"Intercept1 %.6f %.6f\n", stats[num_parts+gcon], stats[num_parts+gcon+total+4+num_parts]);
fprintf(output,"Intercept2 %.6f %.6f\n", stats[total2+num_parts+gcon], stats[total2+num_parts+gcon+total+4+num_parts]);
}
fprintf(output,"Overlap %.6f %.6f\n", stats[2*total2+num_parts], stats[2*total2+num_parts+total+4+num_parts]);
fclose(output);

sprintf(filename2,"%s.cors.full",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2,"Category Trait1_Her SD Trait2_Her SD Both_Coher SD Correlation SD\n");
for(q=0;q<num_parts;q++)
{
if(parttype==0&&q<num_parts-1){fprintf(output2,"A%d ",q+1);}
if(parttype==0&&q==num_parts-1){fprintf(output2,"Base ");}
if(parttype==1){fprintf(output2,"P%d ",q+1);}
fprintf(output2,"%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n", stats[q], stats[q+total+4+num_parts], stats[total2+q], stats[total2+q+total+4+num_parts], stats[2*total2+q], stats[2*total2+q+total+4+num_parts], stats[total+4+q], stats[total+4+q+total+4+num_parts]);
}
fprintf(output2,"ALL %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n", stats[total], stats[total+total+4+num_parts], stats[total+1], stats[total+1+total+4+num_parts], stats[total+2], stats[total+2+total+4+num_parts], stats[total+3], stats[total+3+total+4+num_parts]);
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

printf("Results saved in %s and %s\n\n", filename, filename2);

////////

for(q=0;q<num_parts;q++){free(catlabels[q]);}free(catlabels);
for(j=0;j<count;j++){free(spreds[j]);}free(spreds);free(sal1);free(sal2);free(stags);
for(q=0;q<num_parts;q++){free(svars[q]);free(ssums[q]);}free(svars);free(ssums);
free(snss);free(snss2);free(schis);free(schis2);free(srhos);free(srhos2);
free(stats);

///////////////////////////

