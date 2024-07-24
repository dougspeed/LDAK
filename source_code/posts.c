/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Calculate posterior effects - first with per-predictor variance, then pre-predictor prob 

///////////////////////////

if(countcols(expfile)!=4)
{printf("Error, %s should have 4 columns (not %d), suggesting the file has been changed since creation with \"--calc-exps\"\n\n", expfile, countcols(expfile));exit(1);}

//read in expectations
count=countrows(expfile)-1;
preds=malloc(sizeof(char*)*count);
al1=malloc(sizeof(char)*count);
al2=malloc(sizeof(char)*count);
exps=malloc(sizeof(double)*count);

rs=malloc(sizeof(char)*10000000);

//open and skip first line
printf("Reading expectations for %d predictors from %s\n", count, expfile);
if((input=fopen(expfile,"r"))==NULL)
{printf("Error opening %s\n\n",expfile);exit(1);}
if(fscanf(input, "%s %s %s %s ", readstring, readstring2, readstring2, readstring2)!=4)
{printf("Error reading Row 1 of %s\n\n", expfile);exit(1);}

if(strcmp(readstring,"Predictor")!=0)
{printf("Error reading %s; first element should be \"Predictor\" (not %s), suggesting the file has been changed since creation with \"--calc-exps\"\n\n", expfile, readstring);exit(1);}

for(j=0;j<count;j++)
{
if(fscanf(input, "%s %c %c %lf ", rs, al1+j, al2+j, exps+j)!=4)
{printf("Error reading Row %d of %s\n\n", j+2, expfile);exit(1);}
copy_string(preds,j,rs);
}
fclose(input);

if(extract==1)	//extract
{
predorder=malloc(sizeof(int)*count);
check_dups(preds,count,expfile,predorder,1);

usedpreds=malloc(sizeof(int)*count);
for(j=0;j<count;j++){usedpreds[j]=1;}

count2=extraction(usedpreds, count, preds, chr, predorder, bpredfile, cpredfile, onechr, onesnp, expfile);
if(count2<count){printf("Will use %d of the predictors in %s\n\n", count2, expfile);}
else{printf("Will use all %d predictors in %s\n\n", count, expfile);}

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
exps[count2]=exps[j];
}
count2++;
}}
for(j=count2;j<count;j++){free(preds[j]);}
count=count2;
}

free(predorder);free(usedpreds);
}

//read in summaries (will require all summaries present)
nss=malloc(sizeof(double)*count);
chis=malloc(sizeof(double)*count);
rhos=malloc(sizeof(double)*count);

printf("Reading summary statistics from %s\n", sumsfile);
read_sumsfile(sumsfile, nss, chis, rhos, count, preds, al1, al2, expfile, amb, fixn, scaling, 0, -9999);
printf("First few stats and ns are: %s %.3f %.1f", preds[0], chis[0], nss[0]);
for(j=1;j<count;j++){if(j<3){printf(" | %s %.3f %.1f", preds[j], chis[j], nss[j]);}}
printf("\n\n");

//get average sample size
sum=0;for(j=0;j<count;j++){sum+=nss[j];}
neff=sum/count;
printf("Average sample size is %.1f\n\n", neff);

if(cvar==-9999)	//set prior variance based on average sample size
{cvar=29.72/neff;}
printf("The standardized effect sizes of causal predictors are assumed to come from a Gaussian distribution with mean zero and variance %.4e; change the latter using \"--causal-variance\"\n\n", cvar);

//for each SNP compute posterior mean, variance and odds (assume Var(X)=Var(Y)=1)
pmeans=malloc(sizeof(double)*count);
pvars=malloc(sizeof(double)*count);
bfs=malloc(sizeof(double)*count);
podds=malloc(sizeof(double)*count);

wcount=0;
for(j=0;j<count;j++)
{
//first compute (log) bf assuming prior variance is exps[j]
//mean is value rho, var is value V, bf is value rho^2/2V + log(1-value)/2
//where value=C/(C+V), C=exps[j], V=(1-rho^2)/n = 1/(chi+n) (and rho^2/V=chi)
value=exps[j]/(exps[j]+pow(chis[j]+nss[j],-1));
pmeans[j]=value*rhos[j];
pvars[j]=value/(chis[j]+nss[j]);
bfs[j]=.5*value*chis[j]+.5*log(1-value);

//now compute post prob assuming prior variance is cvar and prior prob is exps[j]/cvar
value3=exps[j]/cvar;
if(value3>.99)
{
if(wcount<5){printf("Warning, %s has prior probability of association %.1f, will be set to 0.99\n\n", preds[j], value3);}
value3=.99;
wcount++;
}

//bf same as above, except now C=cvar
value=cvar/(cvar+pow(chis[j]+nss[j],-1));
value2=.5*value*chis[j]+.5*log(1-value);

//post odds are bf x prior odds
podds[j]=exp(value2+log(value3/(1-value3)));
if(value2+log(value3/(1-value3))>300*log(10)){podds[j]=1e300;}
if(value2+log(value3/(1-value3))<-300*log(10)){podds[j]=1e-300;}
}

if(wcount>5){printf("In total, %d predictors have prior probability greater than 0.99\n", wcount);}
if(wcount>0){printf("\n");}

//save
sprintf(filename,"%s.posts",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fprintf(output,"Predictor A1 A2 Posterior_Mean Posterior_Variance Log10_Bayes_Factor Test_Stat P Posterior_Odds\n");
for(j=0;j<count;j++)
{
fprintf(output, "%s %c %c ", preds[j], al1[j], al2[j]);
fprintf(output, "%.4e %.4e %.4f %.4f %.4e %.4e\n", pmeans[j], pvars[j], bfs[j]/log(10), pow(pmeans[j],2)/pvars[j], erfc(fabs(pmeans[j]*pow(pvars[j],-.5))*M_SQRT1_2), podds[j]);
}
fclose(output);

sprintf(filename2,"%s.score",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2,"Predictor A1 A2 Centre Effect\n");
for(j=0;j<count;j++){fprintf(output2, "%s %c %c NA %.6f\n", preds[j], al1[j], al2[j], pmeans[j]);}
fclose(output2);

sprintf(filename3,"%s.pvalues",outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3,"Predictor P\n");
for(j=0;j<count;j++){fprintf(output3, "%s %.4e\n", preds[j], erfc(fabs(pmeans[j]*pow(pvars[j],-.5))*M_SQRT1_2));}
fclose(output3);

printf("Posterior estimates saved in %s, with a score file in %s (note that these are standardized scores, so if running \"--calc-scores\" you should use \"--power -1)\n\n", filename, filename2);

for(j=0;j<count;j++){free(preds[j]);}free(preds);
free(al1);free(al2);free(exps);free(rs);
free(nss);free(chis);free(rhos);
free(pmeans);free(pvars);free(bfs);free(podds);

///////////////////////////

