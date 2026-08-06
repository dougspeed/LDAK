/*
Copyright 2026 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Join cross-tagging files

///////////////////////////

count=countcols(tagstems[0]);
if(count!=12){printf("Error, %s should have at 12 columns (not %d), suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", tagstems[0], count);exit(1);}

//read top of first file to get flag (0=MAF, 1=Variance)
if((input=fopen(tagstems[0],"r"))==NULL)
{printf("Error opening %s\n\n",tagstems[0]);exit(1);}

//get flag based on element seven
if(fscanf(input, "%s %s %s %s %s %s %s %s %s ", readstring, readstring2, readstring2, readstring2, readstring2, readstring2, readstring3, readstring2, readstring2)!=9)
{printf("Error reading Row 1 of %s\n\n", tagstems[0]);exit(1);}
if(strcmp(readstring3,"MAF")!=0&&strcmp(readstring3,"Variance")!=0)
{printf("Error, the seventh element of %s should be \"MAF\" or \"Variance\" (not %s), suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", tagstems[0], readstring3);exit(1);}
flag=(strcmp(readstring3,"Variance")==0);

fclose(input);

//now do some checks for all files
for(k=0;k<num_tags;k++)
{
count=countcols(tagstems[k]);
if(count!=12){printf("Error, %s should have at 12 columns (not %d), suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", tagstems[k], count);exit(1);}

if((input=fopen(tagstems[k],"r"))==NULL)
{printf("Error opening %s\n\n",tagstems[k]);exit(1);}
if(fscanf(input, "%s %s %s %s %s %s %s %s %s ", readstring, readstring2, readstring2, readstring2, readstring2, readstring2, readstring3, readstring2, readstring2)!=9)
{printf("Error reading Row 1 of %s\n\n", tagstems[k]);exit(1);}

if(strcmp(readstring,"Predictor")!=0)
{printf("Error reading %s; first element should be \"Predictor\" (not %s), suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", tagstems[k], readstring);exit(1);}
if(strcmp(readstring3,"MAF")!=0&&strcmp(readstring3,"Variance")!=0)
{printf("Error, the seventh element of %s should be \"MAF\" or \"Variance\" (not %s), suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", tagstems[k], readstring3);exit(1);}

if(strcmp(readstring2,"Variance")==0&&flag==0)
{printf("Error, %s was created using \"--hwe-stand NO\", but %s was not\n\n", tagstems[k], tagstems[0]);exit(1);}
if(strcmp(readstring2,"MAF")==0&&flag==1)
{printf("Error, %s was created using \"--hwe-stand NO\", but %s was not\n\n", tagstems[0], tagstems[k]);exit(1);}

if(fscanf(input, "%s %s %s ", readstring, readstring2, readstring3)!=3)
{printf("Error reading Elements 10, 11 and 12 of %s, suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", tagstems[k]);exit(1);}
if(strcmp(readstring,"Subset_1")!=0||strcmp(readstring2,"Subset_2")!=0||strcmp(readstring3,"Cross_Term")!=0)
{printf("Error reading %s; the last three elements of Row 1 should be \"Subset_1\", \"Subset_2\" and \"Cross_Term\" (not %s, %s and %s), suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", tagstems[k], readstring, readstring2, readstring3);exit(1);}

fclose(input);
}

////////

//allocate and set ssums to 0
ssums=malloc(sizeof(double*)*3);
ssums[0]=malloc(sizeof(double)*4);
ssums[1]=malloc(sizeof(double)*4);
ssums[2]=malloc(sizeof(double)*4);
ssums[0][0]=0;ssums[0][1]=0;ssums[0][2]=0;ssums[0][3]=0;
ssums[1][0]=0;ssums[1][1]=0;ssums[1][2]=0;ssums[1][3]=0;
ssums[2][0]=0;ssums[2][1]=0;ssums[2][2]=0;ssums[2][3]=0;

//get lengths of each file
pcounts=malloc(sizeof(int)*num_tags);
printf("Counting the number of rows in each tagging file\n\n");
for(k=0;k<num_tags;k++)
{
pcounts[k]=countrows(tagstems[k])-5;
if(pcounts[k]<1)
{printf("Error, %s should have at least six rows (not %d), suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", tagstems[k], pcounts[k]+5);exit(1);}
}

//totals will record numbers of reference, regression and heritability predictors
total=0;
total2=0;
total3=0;

if(checkdups==1)	//check no overlap in predictors
{
count=0;for(k=0;k<num_tags;k++){count+=pcounts[k];}
printf("Checking none of the %d predictors are duplicates; if you are sure this is the case, you can skip this using \"--check-dups NO\"\n", count);

wantpreds=malloc(sizeof(char*)*count);
count=0;
for(k=0;k<num_tags;k++)
{
read_strings(tagstems[k], wantpreds+count, pcounts[k], NULL, 1, 1);
count+=pcounts[k];
}
check_dups(wantpreds, count, "the tagging files", NULL, 1);
printf("\n");

for(j=0;j<count;j++){free(wantpreds[j]);}free(wantpreds);
}

//sort out top line of tagging file
sprintf(filename2,"%s.cross.tagging",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}

fprintf(output2,"Predictor A1 A2 Neighbours Tagging ");
if(flag==0){fprintf(output2,"Weight MAF Categories Exp_Heritability");}
else{fprintf(output2,"Weight Variance Categories Exp_Heritability");}
fprintf(output2," Subset_1 Subset_2 Cross_Term\n");

//now ready to read and print

for(k=0;k<num_tags;k++)
{
//open and skip first line
if((input=fopen(tagstems[k],"r"))==NULL)
{printf("Error opening %s\n\n",tagstems[k]);exit(1);}
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}
printf("Reading tagging details for %d predictors from %s\n", pcounts[k], tagstems[k]);

//now read and print next pcounts[k] lines
for(j=0;j<pcounts[k];j++)
{
if(fscanf(input, "%s ", readstring)!=1)
{printf("Error reading first element of Row %d of %s\n\n", j+2, tagstems[k]);exit(1);}
fprintf(output2, "%s", readstring);
for(q=1;q<12;q++)
{
if(fscanf(input, "%s ", readstring)!=1)
{printf("Error reading Element %d of Row %d of %s\n\n", q+1, j+2, tagstems[k]);exit(1);}
fprintf(output2, " %s", readstring);
}
fprintf(output2, "\n");
}	//end of j loop

for(q2=0;q2<3;q2++)	//get the sums from the partition lines
{
if(fscanf(input, "The %s %s %s %s %s %s %s %s ", readstring, readstring, readstring, readstring, readstring, readstring, readstring, readstring)!=8)
{printf("Error reading first nine elements of Row %d of %s, suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", pcounts[k]+q2+2, tagstems[k]);exit(1);}
if(fscanf(input, "%lf %lf %lf ", &value, &value2, &value3)!=3)
{printf("Error reading Elements 10, 11 and 12 of Row %d of %s, suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", pcounts[k]+q2+2, tagfile);exit(1);}
ssums[0][q2]+=value;
ssums[1][q2]+=value2;
ssums[2][q2]+=value3;
}

//and read the final line
if(fscanf(input, "%s ", readstring)!=1)
{printf("Error reading the first element of Row %d of %s, suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", pcounts[k]+5, tagstems[k]);exit(1);}
if(strcmp(readstring,"The")==0)
{printf("Error, the tagging file %s was made using an older version of LDAK; please remake using this version\n\n", tagfile);exit(1);}
if(strcmp(readstring,"There")!=0)
{printf("Error reading the first element of Row %d of %s, suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", pcounts[k]+5, tagstems[k]);exit(1);}

if(fscanf(input, "%s %d %s %d %s %d %s %s ", readstring, &readint, readstring, &readint2, readstring, &readint3, readstring, readstring)!=8)
{printf("Error reading first 9 elements of Row %d of %s, suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", pcounts[k]+5, tagstems[k]);exit(1);}
total+=readint;
total2+=readint2;
total3+=readint3;

if(fscanf(input, "%lf %lf %lf ", &value, &value2, &value3)!=3)
{printf("Error reading Element 10, 11 and 12 of Row %d of %s, suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", pcounts[k]+5, tagfile);exit(1);}
ssums[0][3]+=value;
ssums[1][3]+=value2;
ssums[2][3]+=value3;

fclose(input);
}	//end of k loop
printf("\n");

//print the sums
for(q2=0;q2<3;q2++)
{
if(q2==2){fprintf(output2, "The relative contribution of cross term to each category");}
else{fprintf(output2, "The relative contribution of Subset %d to each category", q2+1);}
for(q=0;q<3;q++){fprintf(output2, " %.4f", ssums[q][q2]);}
fprintf(output2,"\n");
}

fprintf(output2, "There are %d reference %d regression %d heritability predictors", total, total2, total3);
for(q=0;q<3;q++){fprintf(output2, " %d", (int)ssums[q][3]);}
fprintf(output2, "\n");

fclose(output2);

printf("New cross-ancestry tagging file saved in %s\n\n", filename2);

for(q=0;q<3;q++){free(ssums[q]);}free(ssums);
free(pcounts);

///////////////////////////

