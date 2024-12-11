/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Join tagging files (and maybe matrices)

///////////////////////////

count=countcols(tagstems[0])-9;
if(count<1){printf("Error, %s should have at least 10 columns (not %d), suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", tagstems[0], 9+count);exit(1);}

//read top of first file to get flag (0=MAF, 1=Variance), parttype (0=annotations, 1=partitions) and num_parts
if((input=fopen(tagstems[0],"r"))==NULL)
{printf("Error opening %s\n\n",tagstems[0]);exit(1);}

//get flag based on element seven
if(fscanf(input, "%s %s %s %s %s %s %s %s %s ", readstring, readstring2, readstring2, readstring2, readstring2, readstring2, readstring3, readstring2, readstring2)!=9)
{printf("Error reading Row 1 of %s\n\n", tagstems[0]);exit(1);}
if(strcmp(readstring3,"MAF")!=0&&strcmp(readstring3,"Variance")!=0)
{printf("Error, the seventh element of %s should be \"MAF\" or \"Variance\" (not %s), suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", tagstems[0], readstring3);exit(1);}
flag=(strcmp(readstring3,"Variance")==0);

//get parttype based on final element
for(q=0;q<count;q++)
{
if(fscanf(input, "%s ", readstring)!=1)
{printf("Error reading element %d of %s, suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", 10+q, tagstems[0]);exit(1);}
}
parttype=(strcmp(readstring,"Base")!=0);

//num_parts is count or count-1, depending on whether have a background
num_parts=count-(strcmp(readstring,"Background")==0);

fclose(input);

//get labels, allowing for possibility of background
catlabels=malloc(sizeof(char *)*(num_parts+1));
for(q=0;q<num_parts+1;q++){catlabels[q]=malloc(sizeof(char)*500);}
strcpy(catlabels[num_parts],"Background");

if((input=fopen(tagstems[0],"r"))==NULL)
{printf("Error opening %s\n\n",tagstems[0]);exit(1);}
if(fscanf(input, "%s %s %s %s %s %s %s %s %s ", readstring, readstring2, readstring2, readstring2, readstring2, readstring2, readstring3, readstring2, readstring2)!=9)
{printf("Error reading Row 1 of %s\n\n", tagstems[0]);exit(1);}

for(q=0;q<count;q++)
{
if(fscanf(input, "%s ", catlabels[q])!=1)
{printf("Error reading element %d of %s, suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", 10+q, tagstems[0]);exit(1);}
}
fclose(input);

//now read tops of all files to get addpart, plus many checks
addpart=0;
for(k=0;k<num_tags;k++)
{
count=countcols(tagstems[k])-9;
if(count<1){printf("Error, %s should have at least 10 columns (not %d), suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", tagstems[k], 9+count);exit(1);}

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

for(q=0;q<count;q++)	//now category names
{
if(fscanf(input, "%s ", readstring)!=1)
{printf("Error reading element %d of %s, suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", 10+q, tagstems[k]);exit(1);}
if(strcmp(readstring,catlabels[q])!=0)
{printf("Error, the label of Category %d in %s (%s) does not match that in %s (%s)\n\n", q+1, tagstems[k], readstring, tagstems[0], catlabels[q]);exit(1);}
}

if(strcmp(readstring,"Base")==0&&parttype==1)
{printf("Error, %s was created using partitions, but %s was not\n\n", tagstems[0], tagstems[k]);exit(1);}
if(strcmp(readstring,"Base")!=0&&parttype==0)
{printf("Error, %s was created using partitions, but %s was not\n\n", tagstems[k], tagstems[0]);exit(1);}

if(strcmp(readstring,"Background")==0){addpart2=1;addpart=1;}
else{addpart2=0;}

if(count-addpart2!=num_parts)
{
if(parttype==0){printf("Error, the number of annotations used to construct %s (%d) does not match the number used for %s (%d)\n\n", tagstems[k], count-1, tagstems[0], num_parts-1);exit(1);}
else{printf("Error, the number of partitions used to construct %s (%d) does not match the number used for %s (%d)\n\n", tagstems[k], count-addpart2, tagstems[0], num_parts);exit(1);}
}

if(parttype==0)
{
if(count==1){printf("%s provides taggings for the base\n",tagstems[k]);}
if(count==2){printf("%s provides taggings for 1 annotation plus the base\n",tagstems[k]);}
if(count>2){printf("%s provides taggings for %d annotation plus the base\n",tagstems[k], count-1);}
}
else
{
if(count==1){printf("%s provides taggings for 1 partition\n", tagstems[k]);}
else
{
printf("%s provides taggings for %d partitions", tagstems[k], count);
if(addpart2==1){printf(" (including a background partition)");}
printf("\n");
}
}

fclose(input);
}
printf("\n");

if(strcmp(matlist,"blank")!=0)	//check sizes and labels of matrices
{
for(k=0;k<num_tags;k++)
{
if(countcols(matstems[k])!=countcols(tagstems[k])-8)
{printf("Error, %s should have %d columns (not %d), indicating that it does not provide the heritability matrix corresponding to %s\n\n", matstems[k], countcols(tagstems[k])-8, countcols(matstems[k]), tagstems[k]);exit(1);}

if((input=fopen(matstems[k],"r"))==NULL)
{printf("Error opening %s\n\n",matstems[k]);exit(1);}
if(fscanf(input, "%s ", readstring)!=1)
{printf("Error reading Row 1 of %s\n\n", matstems[k]);exit(1);}

for(q=0;q<count;q++)
{
if(fscanf(input, "%s ", readstring)!=1)
{printf("Error reading Row 1 of %s\n\n", matstems[k]);exit(1);}
if(strcmp(readstring,catlabels[q])!=0)
{printf("Error, the label of Category %d in %s (%s) does not match that in %s (%s)\n\n", q+1, matstems[k], readstring, tagstems[k], catlabels[q]);exit(1);}
}
fclose(input);
}
}

////////

//allocate and set ssums to 0
ssums=malloc(sizeof(double*)*(num_parts+addpart));
for(q=0;q<num_parts+addpart;q++)
{
ssums[q]=malloc(sizeof(double)*(num_parts+addpart+1));
for(q2=0;q2<num_parts+addpart+1;q2++){ssums[q][q2]=0;}
}

//get lengths of each file
pcounts=malloc(sizeof(int)*num_tags);
printf("Counting the number of rows in each tagging file\n\n");
for(k=0;k<num_tags;k++)
{
count=countcols(tagstems[k])-9;
pcounts[k]=countrows(tagstems[k])-count-2;
if(pcounts[k]<1)
{printf("Error, %s should have at least %d rows (not %d), suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", tagstems[k], count+3, pcounts[k]+count+2);exit(1);}
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
sprintf(filename2,"%s.tagging",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}

fprintf(output2,"Predictor A1 A2 Neighbours Tagging ");
if(flag==0){fprintf(output2,"Weight MAF Categories Exp_Heritability");}
else{fprintf(output2,"Weight Variance Categories Exp_Heritability");}
for(q=0;q<num_parts+addpart;q++){fprintf(output2," %s", catlabels[q]);}
fprintf(output2,"\n");

if(strcmp(matlist,"blank")!=0)	//and also of heritability matrix
{
sprintf(filename4,"%s.matrix",outfile);
if((output4=fopen(filename4,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename4);exit(1);}

fprintf(output4,"Predictor");
for(q=0;q<num_parts+addpart;q++){fprintf(output4," %s", catlabels[q]);}
fprintf(output4,"\n");
}

//now ready to read and print

for(k=0;k<num_tags;k++)
{
//start with tagging file
count=countcols(tagstems[k])-9;
if(count==num_parts+addpart){addpart2=0;}
else{addpart2=1;}

//open and skip first line
if((input=fopen(tagstems[k],"r"))==NULL)
{printf("Error opening %s\n\n",tagstems[k]);exit(1);}
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}
printf("Reading tagging details for %d predictors from %s\n", count2, tagstems[k]);

//now read and print next pcounts[k] lines, adding on zero if addpart2=1
for(j=0;j<pcounts[k];j++)
{
if(fscanf(input, "%s ", readstring)!=1)
{printf("Error reading first element of Row %d of %s\n\n", j+2, tagstems[k]);exit(1);}
fprintf(output2, "%s", readstring);
for(q=1;q<9+count;q++)
{
if(fscanf(input, "%s ", readstring)!=1)
{printf("Error reading Element %d of Row %d of %s\n\n", q+1, j+2, tagstems[k]);exit(1);}
fprintf(output2, " %s", readstring);
}
if(addpart2==1){fprintf(output2, "0.0000");}
fprintf(output2, "\n");
}	//end of j loop

for(q2=0;q2<count;q2++)	//get the sums from the partition lines
{
if(fscanf(input, "The %s %s %s %s %s %s %s %s ", readstring, readstring, readstring, readstring, readstring, readstring, readstring, readstring)!=8)
{printf("Error reading first 9 elements of Row %d of %s, suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", pcounts[k]+q2+2, tagstems[k]);exit(1);}
for(q=0;q<count;q++)
{
if(fscanf(input, "%lf ", &value)!=1)
{printf("Error reading Element %d of Row %d of %s, suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", q+10, pcounts[k]+q2+2, tagfile);exit(1);}
ssums[q][q2]+=value;
}
}

//and read the final line (do separately, as must go into num_parts+addpart, and record num preds)
if(fscanf(input, "%s ", readstring)!=1)
{printf("Error reading the first element of Row %d of %s, suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", pcounts[k]+count+2, tagstems[k]);exit(1);}
if(strcmp(readstring,"The")==0)
{printf("Error, the tagging file %s was made using an older version of LDAK; please remake using this version\n\n", tagfile);exit(1);}
if(strcmp(readstring,"There")!=0)
{printf("Error reading the first element of Row %d of %s, suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", pcounts[k]+count+2, tagstems[k]);exit(1);}

if(fscanf(input, "%s %d %s %d %s %d %s %s ", readstring, &readint, readstring, &readint2, readstring, &readint3, readstring, readstring)!=8)
{printf("Error reading first 9 elements of Row %d of %s, suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", pcounts[k]+count+2, tagstems[k]);exit(1);}
total+=readint;
total2+=readint2;
total3+=readint3;

for(q=0;q<count;q++)
{
if(fscanf(input, "%lf ", &value)!=1)
{printf("Error reading Element %d of Row %d of %s, suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", q+10, pcounts[k]+count+2, tagfile);exit(1);}
ssums[q][num_parts+addpart]+=value;
}

fclose(input);

if(strcmp(matlist,"blank")!=0)	//now read and print heritability matrix
{
count=countrows(matstems[k])-1;
if(count!=pcounts[k])
{printf("Error, %s should have %d rows (not %d), suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", matstems[k], pcounts[k]+1, count+1);exit(1);}

printf("Reading heritability details for %d predictors from %s\n", pcounts[k] , matstems[k]);

//open and skip first line
if((input3=fopen(matstems[k],"r"))==NULL)
{printf("Error opening %s\n\n",matstems[k]);exit(1);}
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input3, "%c", &readchar);}

//now read and print next pcounts[k] lines
for(j=0;j<pcounts[k];j++)
{
for(q=0;q<1+count;q++)
{
if(fscanf(input3, "%s ", readstring)!=1)
{printf("Error reading Element %d of Row %d of %s\n\n", q+1, j+2, matstems[k]);exit(1);}
fprintf(output4, "%s ", readstring);
}
fprintf(output4, "\n");
}

fclose(input3);
}
}	//end of k loop
printf("\n");

//print the sums
for(q2=0;q2<num_parts;q2++)
{
if(parttype==0)
{
if(q2==num_parts-1){fprintf(output2, "The relative contribution of the Base to each category");}
else{fprintf(output2, "The relative contribution of Annotation %d to each category", q2+1);}
}
else{fprintf(output2, "The relative contribution of Partition %d to each category", q2+1);}
for(q=0;q<num_parts+addpart;q++){fprintf(output2, " %.4f", ssums[q][q2]);}
fprintf(output2, "\n");
}

if(addpart==1)
{
fprintf(output2, "The relative contribution of the Background to each category");
for(q=0;q<num_parts+addpart;q++){fprintf(output2, " %.4f", ssums[q][num_parts]);}
fprintf(output2, "\n");
}

fprintf(output2, "There are %d reference %d regression %d heritability predictors", total, total2, total3);
for(q=0;q<num_parts+addpart;q++){fprintf(output2, " %d", (int)ssums[q][num_parts+addpart]);}
fprintf(output2, "\n");

fclose(output2);
if(strcmp(matlist,"blank")!=0){fclose(output4);}

printf("New tagging file saved in %s", filename2);
if(strcmp(matlist,"blank")!=0){printf(", with new heritability matrix saved in %s", filename4);}
printf("\n\n");

for(q=0;q<num_parts+1;q++){free(catlabels[q]);}free(catlabels);
for(q=0;q<num_parts;q++){free(ssums[q]);}free(ssums);
free(pcounts);

///////////////////////////

