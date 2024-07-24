/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Adjust summary statistics for winners curse

///////////////////////////

//read predictor details
count=countrows(sumsfile)-1;
printf("Reading details for %d predictors from %s\n\n", count, sumsfile);

preds=malloc(sizeof(char*)*count);
al1=malloc(sizeof(char)*count);
al2=malloc(sizeof(char)*count);

count2=countcols(sumsfile);
if(find_head("Predictor", sumsfile, count2)!=-1){cols[0]=find_head("Predictor", sumsfile, count2);}
else{cols[0]=find_head("SNP", sumsfile, count2);}
cols[1]=find_head("A1", sumsfile, count2);
cols[2]=find_head("A2", sumsfile, count2);

//get predictor names
read_strings(sumsfile, preds, count, NULL, cols[0]+1, 1);

//get alleles
fullal1=malloc(sizeof(char*)*count);
fullal2=malloc(sizeof(char*)*count);
read_strings(sumsfile, fullal1, count, NULL, cols[1]+1, 1);
read_strings(sumsfile, fullal2, count, NULL, cols[2]+1, 1);

for(j=0;j<count;j++)
{
if(strlen(fullal1[j])+strlen(fullal2[j])>2)
{printf("Error reading %s; Predictor %s has multi-character alleles (%s and %s)\n\n", sumsfile, preds[j], fullal1[j], fullal2[j]);exit(1);}
al1[j]=fullal1[j][0];
al2[j]=fullal2[j][0];
}
for(j=0;j<count;j++){free(fullal1[j]);free(fullal2[j]);}free(fullal1);free(fullal2);

if(extract==1)	//extract
{
predorder=malloc(sizeof(int)*count);
check_dups(preds,count,sumsfile,predorder,1);

usedpreds=malloc(sizeof(int)*count);
for(j=0;j<count;j++){usedpreds[j]=1;}

count2=extraction(usedpreds, count, preds, chr, predorder, bpredfile, cpredfile, onechr, onesnp, sumsfile);
if(count2<count){printf("Will use %d of the predictors in %s\n\n", count2, sumsfile);}
else{printf("Will use all %d predictors in %s\n\n", count, sumsfile);}

if(count2<count)	//squeeze down
{
count2=0;
for(j=0;j<count;j++)
{
if(usedpreds[j]==1)
{
if(count2!=j)
{
free(preds[count2]);copy_string(preds,count2,preds[j]);
al1[count2]=al1[j];
al2[count2]=al2[j];
}
count2++;
}}
for(j=count2;j<count;j++){free(preds[j]);}
count=count2;
}

free(predorder);free(usedpreds);
}

//read in summaries - must have values for all predictors
nss=malloc(sizeof(double)*count);
chis=malloc(sizeof(double)*count);
rhos=malloc(sizeof(double)*count);

printf("Reading summary statistics from %s\n", sumsfile);
read_sumsfile(sumsfile, nss, chis, rhos, count, preds, al1, al2, NULL, amb, fixn, scaling, 2, 0);
printf("First few stats and ns are: %s %.3f %.1f", preds[0], chis[0], nss[0]);
for(j=1;j<count;j++){if(j<3){printf(" | %s %.3f %.1f", preds[j], chis[j], nss[j]);}}
printf("\n\n");

//get normal threshold and chisq thresholds
value=-normal_inv(cutoff/2);
value2=pow(value,2);

//count how many exceed chisq threshold
count2=0;for(j=0;j<count;j++){count2+=(chis[j]>=value2);}
if(count2==0){printf("Error, no predictors have p-value below %.2e\n\n", cutoff);exit(1);}
if(count2==1){printf("There is one predictor with p-value below %.2e\n\n", cutoff);}
else{printf("There are %d predictors with p-value below %.2e\n\n", count2, cutoff);}

//save details of significant and non-significant predictors
sprintf(filename,"%s.sig",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
sprintf(filename2,"%s.nonsig",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}

for(j=0;j<count;j++)
{
if(chis[j]>=value2){fprintf(output,"%s\n", preds[j]);}
else{fprintf(output2,"%s\n", preds[j]);}
}
fclose(output);
fclose(output2);

printf("Lists of significant and non-significant predictors saved in %s and %s\n\n", filename, filename2);

////////

//loop through predictors - if significant, update test statistic allowing for winners curse
for(j=0;j<count;j++)
{
if(chis[j]>=value2)	//update
{
//get sd (note that var = (1-rho^2)/n = 1/(chi+n)
value3=pow(chis[j]+nss[j],-.5);

//find max likelihood value
value4=secant_trun(fabs(rhos[j]), value3, value*value3, tol, maxiter);

//update effect and stat
if(rhos[j]>0){rhos[j]=value4;}
else{rhos[j]=-value4;}
chis[j]=1+nss[j]*pow(value4,2)/(1-pow(value4,2));
}
}	//end of j loop

//save
sprintf(filename,"%s.summaries.adjusted",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fprintf(output, "Predictor A1 A2 Direction Stat n\n");
for(j=0;j<count;j++)
{fprintf(output, "%s %c %c %d %.4f %d\n", preds[j], al1[j], al2[j], (rhos[j]>=0)-(rhos[j]<0), chis[j], (int)nss[j]);}
fclose(output);

printf("Adjusted summary statistics saved in %s\n\n", filename);

for(j=0;j<count;j++){free(preds[j]);}free(preds);
free(al1);free(al2);
free(nss);free(chis);free(rhos);

///////////////////////////

