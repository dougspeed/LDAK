/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Calculate expected heritability contributed or tagged by predictors

///////////////////////////

if(mode==146)
{
//read in taus from file
stats=malloc(sizeof(double)*num_parts);
sprintf(filename,"%s.taus",outfile);
read_values(filename,stats,num_parts,NULL,2,1,0);

count=countrows(matfile)-1;
printf("Calculating per-predictor heritabilities for %d predictors\n", count);

if((input=fopen(matfile,"r"))==NULL)
{printf("Error opening %s\n\n",matfile);exit(1);}
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}

sprintf(filename2,"%s.ind.hers",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2, "Predictor Heritability\n");

sprintf(filename3,"%s.ind.hers.positive",outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output3, "Predictor Heritability\n");

for(j=0;j<count;j++)
{
if(fscanf(input, "%s ", readstring)!=1)
{printf("Error reading first element of Row %d of %s\n\n", j+2, matfile);exit(1);}
value=0;
for(q=0;q<num_parts;q++)
{
if(fscanf(input, "%lf ", &value2)!=1)
{printf("Error reading Element %d of Row %d of %s, suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", q+2, j+2, tagfile);exit(1);}
value+=stats[q]*value2;
}

fprintf(output2, "%s %.6e\n", readstring, value);
if(value>0){fprintf(output3, "%s %.6e\n", readstring, value);}
}

fclose(input);
fclose(output2);
fclose(output3);

printf("\nEstimates of per-predictor heritabilities saved in %s and %s\n\n", filename2, filename3);

free(stats);
}

////////

if(mode==149)	//calc-exps
{
//read in tag details - can not be extracting
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

count=read_tagfile(tagfile, catlabels, spreds, sal1, sal2, stags, svars, ssums, num_parts, count, parttype, "blank", "blank", "blank", 0);

//allocations
stats=malloc(sizeof(double)*num_parts);

//stats should contain the tau coefficients
read_values(taufile,stats,num_parts,NULL,2,1,0);

printf("Computing estimates for %d predictors\n", count);

sprintf(filename,"%s.exps",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fprintf(output,"Predictor A1 A2 Heritability\n");

for(j=0;j<count;j++)
{
value=0;
for(q=0;q<num_parts;q++){value+=stats[q]*svars[q][j];}
fprintf(output, "%s %c %c %.4e\n", spreds[j], sal1[j], sal2[j], value);
}

fclose(output);

printf("\nEstimates of heritabilities tagged by predictors saved in %s\n\n", filename);

for(q=0;q<num_parts;q++){free(catlabels[q]);}free(catlabels);
for(j=0;j<count;j++){free(spreds[j]);}free(spreds);free(sal1);free(sal2);free(stags);
for(q=0;q<num_parts;q++){free(svars[q]);free(ssums[q]);}free(svars);free(ssums);
free(stats);
}

///////////////////////////

