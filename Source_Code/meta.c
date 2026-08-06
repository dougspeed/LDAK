/*
Copyright 2026 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.
*/

///////////////////////////

//Construct prediction models from summary statistics - must allow for case that target summary statistics missing for some predictors (happens if unable to impute a window) 
//if num_cors>1, must have power=-9999, or -1, -.75, -.5, -.25 or 0

///////////////////////////

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fclose(output);

//allocate and load multitrait summary statistics
value=(double)data_length/1024*num_sums3*4/1024/1024*8;
if(value>1){printf("Warning, to store the summary statistics requires %.1f Gb\n\n", value);}

Mnss=malloc(sizeof(double*)*num_sums3);
Mchis=malloc(sizeof(double*)*num_sums3);
Mrhos=malloc(sizeof(double*)*num_sums3);
Ma1freq=malloc(sizeof(double*)*num_sums3);

for(q=1;q<num_sums3;q++)
{
Mnss[q]=malloc(sizeof(double)*data_length);
Mchis[q]=malloc(sizeof(double)*data_length);
Mrhos[q]=malloc(sizeof(double)*data_length);
Ma1freq[q]=malloc(sizeof(double)*data_length);
}

//can just reassign for target trait
Mnss[0]=nss;Mchis[0]=chis;Mrhos[0]=rhos;Ma1freq[0]=a1freq;

if(num_sums3>1)  //read in other summary statistics (missing will get ns, rhos, chis and freq zero)
{
if(num_sums3==2){printf("Reading summary statistics from %s\n", sumstems[1]);}
else{printf("Reading %d sets of summary statistics (if using multiple threads, these will be read in a random order)\n", num_sums3-1);}

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
if(num_sums3==2){fprintf(output, "Reading summary statistics from %s\n", sumstems[1]);}
else{fprintf(output, "Reading %d sets of summary statistics (if using multiple threads, these will be read in a random order)\n", num_sums3-1);}
fclose(output);

#pragma omp parallel for private(q) schedule(static)
for(q=1;q<num_sums3;q++)
{read_sumsfile(sumstems[q], Mnss[q], Mchis[q], Mrhos[q], Ma1freq[q], data_length, preds, al1, al2, bimfile, amb, 1.0, checksums);}
printf("\n");

for(q=1;q<num_sums3;q++)
{
max=0;for(j=0;j<data_length;j++){max+=(Mnss[q][j]>max)*(Mnss[q][j]-max);}
count=0;for(j=0;j<data_length;j++){count+=(Mnss[q][j]>0&&Mnss[q][j]<0.5*max);}
sum=0;for(j=0;j<data_length;j++){sum+=Mnss[q][j];}
mean=sum/data_length;
printf("For the summary statistics in %s, the maximum (average) sample size is %.0f (%.1f)\n", sumstems[q], max, mean);
if(count>0){printf("Warning, %d predictors have sample size less than half the maximum (%.1f)\n", count, 0.5*max);}
}
}

////////

printf("Meta-analyzing summary statistics\n\n");

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Meta-analyzing summary statistics\n");
fclose(output);

//do we have freq for all summaries?
flag=1;
for(q=0;q<num_sums3;q++)
{
if(find_head("A1Freq", sumstems[q], -1)==-1){flag=0;break;}
}

sprintf(filename2,"%s.meta",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}

if(flag==1){fprintf(output2,"Predictor\tA1\tA2\tZ\tn\tA1Freq\tP\n");}
else{fprintf(output2,"Predictor\tA1\tA2\tZ\tn\tP\n");}

count=0;
for(j=0;j<data_length;j++)
{
value=0;
value3=0;
sum=0;
for(q=0;q<num_sums3;q++)
{
value+=Mrhos[q][j]*Mnss[q][j];
value3+=Ma1freq[q][j]*Mnss[q][j];
sum+=Mnss[q][j];
}
if(sum>0)
{
value2=pow(sum/(1-pow(value/sum,2)),.5)*value/sum;
if(flag==1){fprintf(output2,"%s\t%c\t%c\t%.4f\t%0.f\t%.4f\t%.4e\n", preds[j], al1[j], al2[j], value2, sum, value3/sum, erfc(fabs(value2)*M_SQRT1_2));}
else{fprintf(output2,"%s\t%c\t%c\t%.4f\t%0.f\t%.4e\n", preds[j], al1[j], al2[j], value2, sum, erfc(fabs(value2)*M_SQRT1_2));}
count++;
}
}
fclose(output2);

if(count==0){printf("Error, none of the %d predictors have valid summary statistics\n\n", data_length);exit(1);}
if(count==data_length){printf("All %d predictors have valid summary statistics\n\n", data_length);}
else{printf("Warning, only %d of the %d predictors have valid summary statistics\n\n", count, data_length);}

printf("The meta-analyzed summary statistics are saved in %s\n\n", filename2);

////////

for(q=1;q<num_sums3;q++){free(Mnss[q]);free(Mchis[q]);free(Mrhos[q]);free(Ma1freq[q]);}
free(Mnss);free(Mchis);free(Mrhos);free(Ma1freq);

///////////////////////////

