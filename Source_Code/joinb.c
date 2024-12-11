/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Join pathway tagging files

///////////////////////////

//check sums files exist and get num_parts
sprintf(filename,"%s.pathway.sums",tagstems[0]);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--calc-taggings\"\n\n", filename);exit(1);}
num_parts=countrows(filename);

for(k=1;k<num_tags;k++)
{
sprintf(filename,"%s.pathway.sums",tagstems[k]);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--calc-taggings\"\n\n", filename);exit(1);}
count=countrows(filename);
if(count!=num_parts){printf("Error, %s was constructed using %d annotations while %s was constructed using %d annotations\n\n", tagstems[k], num_parts, tagstems[0], count);exit(1);}
}

//allocations
pcounts=malloc(sizeof(int)*num_tags);
Minput=malloc(sizeof(FILE *)*num_tags);
readfloats=malloc(sizeof(float)*12);
writefloats=malloc(sizeof(float)*12);

////////

//start with base.tagging files

sprintf(filename,"%s.base.tagging",tagstems[0]);
if(countcols(filename)!=11){printf("Error, %s should have 11 columns (not %d), suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", filename, countcols(filename));exit(1);}

//read top of first file to get flag (0=MAF, 1=Variance)
if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}

//get flag based on element seven
if(fscanf(input, "%s %s %s %s %s %s %s %s %s ", readstring, readstring2, readstring2, readstring2, readstring2, readstring2, readstring3, readstring2, readstring2)!=9)
{printf("Error reading Row 1 of %s\n\n", filename);exit(1);}
if(strcmp(readstring3,"MAF")!=0&&strcmp(readstring3,"Variance")!=0)
{printf("Error, the seventh element of %s should be \"MAF\" or \"Variance\" (not %s), suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", filename, readstring3);exit(1);}
flag=(strcmp(readstring3,"Variance")==0);

fclose(input);

//check tops of remaining files
for(k=1;k<num_tags;k++)
{
sprintf(filename,"%s.base.tagging",tagstems[k]);
if(countcols(filename)!=11){printf("Error, %s should have 11 columns (not %d), suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", filename, countcols(filename));exit(1);}

if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}
if(fscanf(input, "%s %s %s %s %s %s %s %s %s ", readstring, readstring2, readstring2, readstring2, readstring2, readstring2, readstring3, readstring2, readstring2)!=9)
{printf("Error reading Row 1 of %s\n\n", filename);exit(1);}

if(strcmp(readstring,"Predictor")!=0)
{printf("Error reading %s; first element should be \"Predictor\" (not %s), suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", filename, readstring);exit(1);}
if(strcmp(readstring3,"MAF")!=0&&strcmp(readstring3,"Variance")!=0)
{printf("Error, the seventh element of %s should be \"MAF\" or \"Variance\" (not %s), suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", filename, readstring3);exit(1);}

if(strcmp(readstring2,"Variance")==0&&flag==0)
{printf("Error, %s was created using \"--hwe-stand NO\", but %s was not\n\n", tagstems[k], tagstems[0]);exit(1);}
if(strcmp(readstring2,"MAF")==0&&flag==1)
{printf("Error, %s was created using \"--hwe-stand NO\", but %s was not\n\n", tagstems[0], tagstems[k]);exit(1);}

fclose(input);
}

//get lengths of each file
printf("Counting the number of rows in each base tagging file\n\n");
for(k=0;k<num_tags;k++)
{
sprintf(filename,"%s.base.tagging",tagstems[k]);
pcounts[k]=countrows(filename)-1;
if(pcounts[k]<1)
{printf("Error, %s should have at least 2 rows (not %d), suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", filename, pcounts[k]+1);exit(1);}
}

if(checkdups==1)	//check no overlap in predictors
{
count=0;for(k=0;k<num_tags;k++){count+=pcounts[k];}
printf("Checking none of the %d predictors are duplicates; if you are sure this is the case, you can skip this using \"--check-dups NO\"\n", count);

wantpreds=malloc(sizeof(char*)*count);
count=0;
for(k=0;k<num_tags;k++)
{
sprintf(filename,"%s.base.tagging",tagstems[k]);
read_strings(filename, wantpreds+count, pcounts[k], NULL, 1, 1);
count+=pcounts[k];
}
check_dups(wantpreds, count, "the tagging files", NULL, 1);
printf("\n");

for(j=0;j<count;j++){free(wantpreds[j]);}free(wantpreds);
}

//sort out top line of base.tagging file
sprintf(filename2,"%s.base.tagging",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}

fprintf(output2,"Predictor A1 A2 Neighbours Tagging ");
if(flag==0){fprintf(output2,"Weight MAF Categories Exp_Heritability Base Genic\n");}
else{fprintf(output2,"Weight Variance Categories Exp_Heritability Base Genic\n");}

//now ready to read and print

for(k=0;k<num_tags;k++)
{
//open and skip first line
sprintf(filename,"%s.base.tagging",tagstems[k]);
if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}
printf("Reading tagging details for %d predictors from %s\n", pcounts[k], filename);

//now read and print next pcounts[k] lines
for(j=0;j<pcounts[k];j++)
{
if(fscanf(input, "%s ", readstring)!=1)
{printf("Error reading Element %d of Row %d of %s\n\n", q+1, j+2, filename);exit(1);}
fprintf(output2, "%s", readstring);
for(q=1;q<11;q++)
{
if(fscanf(input, "%s ", readstring)!=1)
{printf("Error reading Element %d of Row %d of %s\n\n", q+1, j+2, filename);exit(1);}
fprintf(output2, " %s", readstring);
}
fprintf(output2, "\n");
}	//end of j loop

fclose(input);
}	//end of k loop
printf("\n");

////////

//now pathway tagging and sums files

//check pathway tagging files exist and right sizes, then open
for(k=0;k<num_tags;k++)
{
sprintf(filename,"%s.pathway.tagging",tagstems[k]);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--calc-taggings\"\n\n", filename);exit(1);}

count=countcols(filename)-1;
if(count!=pcounts[k]){printf("Error, number of predictors in %s (%d) does not match number in %s.base.tagging (%d)\n\n", filename, count, tagstems[k], pcounts[k]);exit(1);}

if((Minput[k]=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}
}

//open pathway.tagging file
sprintf(filename2,"%s.pathway.tagging",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}

//now merge
for(q=0;q<num_parts;q++)
{
for(k=0;k<num_tags;k++)
{
if(k==0)	//read and print first element
{
if(fscanf(Minput[0], "%s ", readstring)!=1)
{printf("Error reading first element of Row %d of %s.pathway.tagging\n\n", q+1, tagstems[0]);}
fprintf(output2,"%s", readstring);
}
else	//read and check first element
{
if(fscanf(Minput[k], "%s ", readstring2)!=1)
{printf("Error reading first element of Row %d of %s.pathway.tagging\n\n", q+1, tagstems[k]);}
if(strcmp(readstring,readstring2)!=0)
{printf("Error, the first elements of Row %d of %s.pathway.tagging and %s.pathway.tagging (%s and %s) do not match\n\n", q+1, tagstems[0], tagstems[k], readstring, readstring2);exit(1);}
}

//read and print remaining elements
for(j=0;j<pcounts[k];j++)
{
if(fscanf(Minput[k], "%s ", readstring2)!=1)
{printf("Error reading Element %d of Row %d of %s.pathway.tagging\n\n", j+2, q+1, tagstems[k]);}
fprintf(output2," %s", readstring2);
}
}
fprintf(output2,"\n");
}

for(k=0;k<num_tags;k++){fclose(Minput[k]);}

////////

//deal with sums files

//open sums
for(k=0;k<num_tags;k++)
{
sprintf(filename,"%s.pathway.sums",tagstems[k]);
if((Minput[k]=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}
}

//open pathway.sums file
sprintf(filename2,"%s.pathway.sums",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}

//now merge
for(q=0;q<num_parts;q++)
{
if(fscanf(Minput[0], "%s %f %f %f %f %f %f %f %f %f %f %f %f ", readstring, writefloats, writefloats+1, writefloats+2, writefloats+3, writefloats+4, writefloats+5, writefloats+6, writefloats+7, writefloats+8, writefloats+9, writefloats+10, writefloats+11)!=13)
{printf("Error reading Row %d of %s.pathway.sums\n\n", q+1, tagstems[0]);}

for(k=1;k<num_tags;k++)
{
if(fscanf(Minput[k], "%s %f %f %f %f %f %f %f %f %f %f %f %f ", readstring2, readfloats, readfloats+1, readfloats+2, readfloats+3, readfloats+4, readfloats+5, readfloats+6, readfloats+7, readfloats+8, readfloats+9, readfloats+10, readfloats+11)!=13)
{printf("Error reading Row %d of %s.pathway.sums\n\n", q+1, tagstems[k]);}

if(strcmp(readstring,readstring2)!=0)
{printf("Error, the first elements of Row %d of %s.pathway.sums and %s.pathway.sums (%s and %s) do not match\n\n", q+1, tagstems[0], tagstems[k], readstring, readstring2);exit(1);}

for(j=0;j<12;j++){writefloats[j]+=readfloats[j];}
}

fprintf(output2,"%s", readstring);
fprintf(output2," %.4f %.4f %.4f %d", writefloats[0], writefloats[1], writefloats[2], (int)writefloats[3]);
fprintf(output2," %.4f %.4f %.4f %d", writefloats[4], writefloats[5], writefloats[6], (int)writefloats[7]);
fprintf(output2," %.4f %.4f %.4f %d\n", writefloats[8], writefloats[9], writefloats[10], (int)writefloats[11]);
}
fclose(output2);

for(k=0;k<num_tags;k++){fclose(Minput[k]);}

free(pcounts);
free(Minput);free(readfloats);free(writefloats);

///////////////////////////

