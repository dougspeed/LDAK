/*
Copyright 2026 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.f

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//compare frequencies between summary statistics and correlations, maybe infer ancestry for bulks, and blank statistics

///////////////////////////

//read centres for all correlations (allowing for filtering, missing values and flips)
centres2=malloc(sizeof(double)*data_length*num_cors);

for(s=0;s<num_cors;s++)
{
sprintf(filename2,"%s.cors.bin", corstems[s]);
if((input2=fopen(filename2,"rb"))==NULL)
{printf("Error re-opening %s\n\n",filename2);exit(1);}

for(j=0;j<data_length;j++){centres2[j+s*data_length]=0;}
for(j=0;j<Dnp2[s];j++)
{
fseeko(input2, sizeof(double)*Dkp[s][j], SEEK_SET);
if(fread(centres2+Dkp2[s][j]+s*data_length, sizeof(double), 1, input2)!=1)
{printf("Error reading mean for Predictor %d from %s\n\n", j+1, filename2);exit(1);}

if(Dal1[s][Dkp[s][j]]!=al1[Dkp2[s][j]]){centres2[Dkp2[s][j]+s*data_length]=2-centres2[Dkp2[s][j]+s*data_length];}
}

fclose(input2);
}

if(metasum==1)  //infer cohers2 for focal summaries (must be mode=159 and will have frequencies)
{
printf("Inferring the ancestry of the summary statistics in %s\n", sumstems[0]);

STS=malloc(sizeof(double)*num_cors*num_cors);
STS2=malloc(sizeof(double)*num_cors);
STM=malloc(sizeof(double)*num_cors);
STM2=malloc(sizeof(double)*num_cors);

for(s=0;s<num_cors;s++)
{
for(s2=0;s2<num_cors;s2++){STS[s+s2*num_cors]=0;}
STM[s]=0;
}

count=0;
for(j=0;j<data_length;j++)
{
flag=1;
for(s=0;s<num_cors;s++)
{
if(centres2[j+s*data_length]==0){flag=0;break;}
}
if(flag==1)
{
alpha=1.0;beta=1.0;
dgemm_("T", "N", &num_cors, &num_cors, &one, &alpha, centres2+j, &data_length, centres2+j, &data_length, &beta, STS, &num_cors);

for(s=0;s<num_cors;s++){STM[s]+=centres2[j+s*data_length]*2*a1freq[j];}
count++;
}
}
if(count==0){printf("Error, there are no predictors common to all correlations\n\n");exit(1);}

//compute proportions
eigen_invert(STS, num_cors, STS2, 1, STM, 1);

//save original values, then truncate values outside of (0,1) and scale so sum to one, save results in cohers2
for(s=0;s<num_cors;s++){STM2[s]=STM[s];}
sum=0;
for(s=0;s<num_cors;s++)
{
if(STM[s]<0){STM[s]=0;}
if(STM[s]>1){STM[s]=1;}
sum+=STM[s];
}
for(s=0;s<num_cors;s++){cohers2[s]=STM[s]/sum;}

printf("The estimated ancestry proportions are\n");
for(s=0;s<num_cors;s++){printf("Correlations with stem %s: %.4f (estimate before truncation / scaling %.4f)\n", corstems[s], cohers2[s], STM2[s]);}
printf("\n");

if(cohers2[sumpops[0]]==0){printf("Error, the assigned correlations (those with stem %s) have estimated proportion zero; please reassign\n\n", corstems[sumpops[0]]);exit(1);}

free(STS);free(STS2);free(STM);free(STM2);
}

////////

if(gotfreq==1)
{
if(num_cors==1){printf("Measuring the similarity between the summary statistics and the predictor-predictor correlations\n");}

flag=0;
for(q=0;q<num_sums1;q++)	//find mean difference in frequency across correlations (ignoring missing predictors and those with trivial frequency)
{
if(num_cors>1)
{
if(num_sums1==1){printf("Measuring the similarity between the summary statistics and the %d sets of predictor-predictor correlations\n", num_cors);}
else{printf("Measuring the similarity between the summary statistics in %s and the %d sets of predictor-predictor correlations\n", sumstems[q], num_cors);}
}

best=-1;
for(s=0;s<num_cors;s++)
{
sum=0;count=0;
for(j=0;j<data_length;j++)
{
if(Mnss[q][j]>0&&centres2[j+s*data_length]!=0)	//present in summary statistics and correlations
{sum+=fabs(Ma1freq[q][j]-centres2[j+s*data_length]/2);count++;}
}

if(count>0)
{
mean=sum/count;
if(num_cors==1)
{
if(num_sums1==1){printf("The mean difference in frequencies is %.4f\n", mean);}
else{printf("For the summary statistics in %s, the mean difference in frequencies is %.4f\n", sumstems[q], mean);}
if(mean>0.01){printf("Warning, this is quite high (ideally the difference is below 0.01)\n");}
}
else{printf("When compared with Set %d, the mean difference in frequencies is %.4f\n", s+1, mean);}

if(best==-1){best=s;value=mean;}
if(mean<value){best=s;value=mean;}
}
else
{
if(sumpops[q]==s){printf("Error, there is no overlap between the predictors in %s and the correlations with stem %s\n\n", sumstems[q], corstems[s]);exit(1);}
printf("Warning, there is no overlap between the predictors in %s and the correlations with stem %s\n\n", sumstems[q], corstems[s]);
}
}

if(num_cors>1)
{
printf("The best-fitting correlations are Set %d (those with stem %s)\n\n", best+1, corstems[best]);
if(best!=sumpops[q])
{
if(checkpops==1){printf("Error, the summary statistics have instead been assigned to Set %d (those with stem %s)\n\n", sumpops[q]+1, corstems[sumpops[q]]);flag++;}
else{printf("Warning, the summary statistics have instead been assigned to Set %d (those with stem %s)\n\n", sumpops[q]+1, corstems[sumpops[q]]);}
}
}

if(q==0){value2=value;}
}
if(num_cors==1){printf("\n");}

if(flag>0)
{
if(flag<num_sums1){printf("In total, %d of the summary statistics appear to have been incorrectly assigned\n", flag);}
printf("If you are confident that the summary statistics are assigned to the correct predictor-predictor correlations, you can use \"--check-assignments NO\" to override this error\n\n");
exit(1);
}

flag2=0;
for(q=0;q<num_sums1;q++)	//compute per-predictor differences (maybe filtering)
{
if(q==0)    //will save and measure overall similarity
{
if(num_focals==1){sprintf(filename2,"%s.frequency.check",outfile);}
else{sprintf(filename2,"%s.focal%d.frequency.check", outfile, cur_focal+1);}
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2, "Predictor A1 A2 A1Freq_Cor A1Freq_SS Difference Kept\n");

sum2=0;
value2=0;
count2=0;
}

count=0;
for(j=0;j<data_length;j++)
{
if(Mnss[q][j]>0&&centres2[j+sumpops[q]*data_length]!=0)	//present in summary statistics and correlations
{
if(metasum==0||q>0) //easy case - one ancestry
{value=centres2[j+sumpops[q]*data_length]/2;}
else    //hard case - need to get average frequency across contributing ancestries
{
value3=0;
sum=0;  //fact that cohers2[sumpops[0]]>0 ensures sum>0 
for(s=0;s<num_cors;s++)
{
if(centres2[j+s*data_length]!=0){value3+=cohers2[s]*centres2[j+s*data_length]/2;sum+=cohers2[s];}
}
value=value3/sum;
}

flag=-1;
value3=Ma1freq[q][j];
if(checkfreq==1)
{
if(fabs(value-Ma1freq[q][j])>maxfreq)   //exclude predictor
{Mnss[q][j]=0;Mchis[q][j]=0;Mrhos[q][j]=0;Ma1freq[q][j]=-9999;count++;flag=0;}
else    //keep predictor
{flag=1;}
}

if(q==0)    //save and measure overall similarity
{
if(flag==-1){fprintf(output2, "%s %c %c %.6f %.6f %.6f NA\n", preds[j], al1[j], al2[j], value, value3, fabs(value-value3));}
if(flag==0){fprintf(output2, "%s %c %c %.6f %.6f %.6f 0\n", preds[j], al1[j], al2[j], value, value3, fabs(value-value3));}
if(flag==1){fprintf(output2, "%s %c %c %.6f %.6f %.6f 1\n", preds[j], al1[j], al2[j], value, value3, fabs(value-value3));}

sum2+=fabs(value-value3);
if(fabs(value-value3)>value2){value2=fabs(value-value3);}
count2++;
}
}}

if(q==0){fclose(output2);}

if(count>0)
{
if(num_sums3==1){printf("Warning, %d predictors are excluded due to frequency differences\n", count);}
else{printf("Warning, %d of the predictors in %s are excluded due to frequency differences\n", count, sumstems[q]);}
flag2=1;
}
}
if(flag2==1){printf("\n");}

if(mcmcfinal==0&&finemap==0)    //look at focal similarity
{
mean2=sum2/count2;
if(mean2<0.005&&value2<0.1)
{
if(checkmcmc==1){printf("Error, the predictor-predictor correlations with stem %s are a very good match to the summary statistics in %s (the mean difference in frequencies is %.4f, while the maximum difference is %.4f), so we recommend using \"--MCMC-solve YES\" to switch from variational Bayes to MCMC (you can override this error using \"--check-MCMC NO\")\n\n", corstems[sumpops[0]], sumstems[0], mean2, value2);exit(1);}
else{printf("Warning, the predictor-predictor correlations with stem %s are a very good match to the summary statistics in %s (the mean difference in frequencies is %.4f, while the maximum difference is %.4f), so it is normally better to switch from variational Bayes to MCMC\n\n", corstems[sumpops[0]], sumstems[0], mean2, value2);}
}
}
}
else    //do not have frequencies
{
if(num_focals==1){sprintf(filename2,"%s.frequency.check",outfile);}
else{sprintf(filename2,"%s.focal%d.frequency.check", outfile, cur_focal+1);}
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2,"The summary statistic file %s does not contain a column called A1Freq\n", sumstems[0]);
fclose(output2);
}

//blank summary statistics if not present in the corresponding correlations (probably no need to do for q=0)
flag=0;
for(q=0;q<num_sums1;q++)	
{
count=0;
for(j=0;j<data_length;j++)
{
if(Mnss[q][j]>0&&centres2[j+sumpops[q]*data_length]==0){Mnss[q][j]=0;Mchis[q][j]=0;Mrhos[q][j]=0;count++;}
}
if(count>0){printf("Warning, %d of the predictors in %s are not in the corresponding predictor-predictor correlations (so are excluded)\n", count, sumstems[q]);flag=1;}
}
if(flag==1){printf("\n");}

////////

if(num_sums2>0) //must be mode=159 - maybe decide bulk summaries correlations, filter based on frequencies and blank if not present in corresponding correlations
{
if(num_cors==1){printf("Measuring the similarity between the bulk summary statistics and the predictor-predictor correlations\n");}
else{printf("Measuring the similarity between the bulk summary statistics and the %d sets of predictor-predictor correlations\n", num_cors);}

best=-1;
for(s=0;s<num_cors;s++)
{
sum=0;count=0;
for(j=0;j<data_length;j++)
{
if(Ma1freq[num_sums1][j]>0&&centres2[j+s*data_length]!=0)	//present in summary statistics and correlations
{sum+=fabs(Ma1freq[num_sums1][j]-centres2[j+s*data_length]/2);count++;}
}

if(count>0)
{
mean=sum/count;

if(num_cors==1)
{
printf("The mean difference in frequencies is %.4f\n", mean);
if(mean>0.01)
{
if(checkpops==1){printf("Error, this is quite high (ideally the difference is below 0.01), indicating you should include predictor-predictor correlations from an additional ancestry (note that the bulk summary statistics provided at www.dougspeed.com are from UK Biobank GWAS); if you are confident the correlations are correct, you can override this error using \"--check-assignments NO\")\n\n");exit(1);}
printf("Warning, this is quite high (ideally the difference is below 0.01)\n");
}
}
else{printf("When compared with Set %d, the mean difference in frequencies is %.4f\n", s+1, mean);}

if(best==-1){best=s;value=mean;}
if(mean<value){best=s;value=mean;}
}
else{printf("Warning, there is no overlap between the bulk summary statistics and the correlations with stem %s\n", corstems[s]);}
}
if(best==-1){printf("Error, there is no overlap with any correlations\n\n");exit(1);}

if(num_cors>1)
{
printf("The best-fitting correlations are Set %d (those with stem %s)\n", best+1, corstems[best]);
if(value>0.01)
{
if(checkpops==1){printf("Error, the corresponding mean is quite high (ideally the difference is below 0.01), indicating you should include predictor-predictor correlations from an additional ancestry (note that the bulk summary statistics provided at www.dougspeed.com are from UK Biobank GWAS); if you are confident the correlations are correct, you can override this error using \"--check-assignments NO\")\n\n");exit(1);}
printf("Warning, the corresponding mean is quite high (ideally the difference is below 0.01)\n");
}
}

for(q=num_sums1;q<num_sums3;q++){sumpops[q]=best;}

if(checkfreq==1)
{
count=0;
for(j=0;j<data_length;j++)
{
if(Ma1freq[num_sums1][j]>0&&centres2[j+sumpops[num_sums1]*data_length]!=0)	//present in summary statistics and correlations
{
value=centres2[j+sumpops[num_sums1]*data_length]/2;
if(fabs(value-Ma1freq[num_sums1][j])>maxfreq)   //exclude predictor
{
for(q=num_sums1;q<num_sums3;q++){Mnss[q][j]=0;Mchis[q][j]=0;Mrhos[q][j]=0;}
count++;
}
}
}
if(count>0){printf("Warning, %d of the predictors in %s are excluded due to frequency differences\n", count, bulkfile);}
}

count2=0;
for(j=0;j<data_length;j++)
{
if(Ma1freq[num_sums1][j]>0&&centres2[j+sumpops[num_sums1]*data_length]==0)
{
for(q=num_sums1;q<num_sums3;q++){Mnss[q][j]=0;Mchis[q][j]=0;Mrhos[q][j]=0;count2++;}
}
}
if(count2>0){printf("Warning, %d of the predictors in %s are not in the corresponding predictor-predictor correlations (so are excluded)\n", count2, bulkfile);}

printf("\n");
}

if(mode==180&&gotfreq==1)   //copy in frequencies for predictors we might impute
{
q=0;
for(j=0;j<data_length;j++)
{
if(Mnss[q][j]==0&&centres2[j+sumpops[q]*data_length]!=0)	//might impute
{Ma1freq[q][j]=centres2[j+sumpops[q]*data_length]/2;}
}
}

free(centres2);

///////////////////////////

