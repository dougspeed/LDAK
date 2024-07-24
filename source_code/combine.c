/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Merge tagging files - this will delete category labels, but that is not a problem

///////////////////////////

//check each file and get number of partitions
pcounts=malloc(sizeof(int)*num_tags);
for(k=0;k<num_tags;k++)
{
pcounts[k]=countcols(tagstems[k])-9;	
if(pcounts[k]<1)
{printf("Error, %s should have at least 10 columns (not %d), suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", tagstems[k], 9+pcounts[k]);exit(1);}

if((input=fopen(tagstems[k],"r"))==NULL)
{printf("Error opening %s\n\n",tagstems[k]);exit(1);}
if(fscanf(input, "%s %s %s %s %s %s %s %s %s ", readstring, readstring2, readstring2, readstring2, readstring2, readstring2, readstring3, readstring2, readstring2)!=9)
{printf("Error reading Row 1 of %s\n\n", tagstems[k]);exit(1);}

if(strcmp(readstring,"Predictor")!=0)
{printf("Error reading %s; first element should be \"Predictor\" (not %s), suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", tagstems[k], readstring);exit(1);}
//if(strcmp(readstring3,"MAF")!=0&&strcmp(readstring3,"Variance")!=0)
//{printf("Error, the seventh element of %s should be \"MAF\" or \"Variance\" (not %s), suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", tagstems[k], readstring3);exit(1);}

fclose(input);
}

//and lengths should match
count=countrows(tagstems[0])-2-pcounts[0];
for(k=1;k<num_tags;k++)
{
count2=countrows(tagstems[k])-2-pcounts[k];	
if(count2!=count)
{printf("Error, %s contains %d predictors, but %s has %d (all tagging files should contain the same predictors)\n\n", tagstems[k], count2, tagstems[0], count);exit(1);}
}

//get final number of partitions
num_parts=0;for(k=0;k<num_tags;k++){num_parts+=pcounts[k];}

////////

//read details, checking they match those in first
preds=malloc(sizeof(char*)*count);
al1=malloc(sizeof(char)*count);
al2=malloc(sizeof(char)*count);
tally1=malloc(sizeof(double)*count);
tally2=malloc(sizeof(double)*count);
tally3=malloc(sizeof(double)*count*num_parts);
usedpreds=malloc(sizeof(int)*count);
ssums=malloc(sizeof(double*)*num_parts);
for(q=0;q<num_parts;q++)
{
ssums[q]=malloc(sizeof(double)*(num_parts+1));
for(q2=0;q2<num_parts;q2++){ssums[q][q2]=0;}
}

rs=malloc(sizeof(char)*10000000);

mark=0;
for(k=0;k<num_tags;k++)
{
//open tagging file and skip to end of first row
printf("Reading tagging details from %s\n", tagstems[k]);
if((input=fopen(tagstems[k],"r"))==NULL)
{printf("Error opening %s\n\n",tagstems[k]);exit(1);}
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}

for(j=0;j<count;j++)
{
if(fscanf(input, "%s %c %c %lf %lf %s %s %d %s ", rs, &readchar, &readchar2, &readdouble, &readdouble2, readstring, readstring, &readint, readstring)!=9)
{printf("Error reading Row %d of %s\n\n", j+2, tagstems[k]);exit(1);}

for(q=0;q<pcounts[k];q++)
{
if(fscanf(input, "%lf ", tally3+(size_t)(mark+q)*count+j)!=1)
{printf("Error reading Element %d of Row %d of %s\n\n", 10+q, j+2, tagstems[k]);exit(1);}
}

if(k==0)	//load up preds, al1, al2, neighbours (tally1), tagging (tally2) and numcats
{
copy_string(preds,j,rs);
al1[j]=readchar;
al2[j]=readchar2;
tally1[j]=readdouble;
tally2[j]=readdouble2;
usedpreds[j]=readint;
}
else	//check them against loaded values (for numcats add on)
{
if(strcmp(rs,preds[j])!=0)
{printf("Error, Row %d of %s corresponds to %s, while Row %d of %s corresponds to %s (tagging files must be created from the same predictors)\n\n", j+2, tagstems[k], rs, j+2, tagstems[0], preds[j]);exit(1);}

if((readchar!=al1[j]&&readchar!=al2[j])||(readchar2!=al1[j]&&readchar2!=al2[j]))
{printf("Error, %s has alleles %c and %c in %s, but %c and %c in %s\n\n", preds[j], readchar, readchar2, tagstems[k], al1[j], al2[j], tagstems[0]);exit(1);}

if(readdouble!=tally1[j])
{printf("Error, %s has %.6f neighbours in %s, but %.6f in %s (make sure you use the same value for \"--window-kb\" or \"--window-cm\" when creating tagging files)\n\n", preds[j], readdouble, tagstems[k], tally1[j], tagstems[0]);exit(1);}

if(readdouble2!=tally2[j])
{printf("Error, %s has tagging %.6f in %s, but %.6f in %s (make sure you use the same value for \"--window-kb\" or \"--window-cm\" when creating tagging files)\n\n", preds[j], readdouble2, tagstems[k], tally2[j], tagstems[0]);exit(1);}

usedpreds[j]+=readint;
}
}	//end of j loop

//now penultimate rows
for(q2=0;q2<pcounts[k];q2++)
{
if(fscanf(input, "The %s %s %s %s %s %s %s %s ", readstring, readstring, readstring, readstring, readstring, readstring, readstring, readstring)!=8)
{printf("Error reading Row %d of %s, suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", count+2+q2, tagstems[k]);exit(1);}
for(q=0;q<pcounts[k];q++)
{
if(fscanf(input, "%lf ", ssums[mark+q]+mark+q2)!=1)
{printf("Error reading Row %d of %s\n\n", count+2+q2, tagstems[k]);exit(1);}
}
}

//and final row
if(fscanf(input, "%s ", readstring)!=1)
{printf("Error reading the first element of Row %d of %s, suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", count+2+pcounts[k], tagstems[k]);exit(1);}
if(strcmp(readstring,"The")==0)
{printf("Error, the tagging file %s was made using an older version of LDAK; please remake using this version\n\n", tagstems[k]);exit(1);}
if(strcmp(readstring,"There")!=0)
{printf("Error reading the first element of Row %d of %s, suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", count+2+pcounts[k], tagstems[k]);exit(1);}

if(fscanf(input, "%s %d %s %d %s %d %s %s ", readstring, &readint, readstring, &readint2, readstring, &readint3, readstring, readstring)!=8)
{printf("Error reading first 9 elements of Row %d of %s, suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", count+2+pcounts[k], tagstems[k]);exit(1);}

if(k==0)
{
total=readint;
total2=readint2;
total3=readint3;
}
else
{
if(readint!=total){printf("Error, the number of reference predictors in %s (%d) does not match the number in %s (%d)\n\n", tagstems[k], readint, tagstems[0], total);exit(1);}
if(readint2!=total2){printf("Error, the number of regression predictors in %s (%d) does not match the number in %s (%d)\n\n", tagstems[k], readint2, tagstems[0], total2);exit(1);}
if(readint3!=total3){printf("Error, the number of heritability predictors in %s (%d) does not match the number in %s (%d)\n\n", tagstems[k], readint3, tagstems[0], total3);exit(1);}
}

for(q=0;q<pcounts[k];q++)
{
if(fscanf(input, "%lf ", ssums[mark+q]+num_parts)!=1)
{printf("Error reading Row %d of %s\n\n", count+2+pcounts[k], tagstems[k]);exit(1);}
}

fclose(input);
mark+=pcounts[k];
}	//end of k loop
printf("\n");

//ready to print
sprintf(filename2,"%s.tagging",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}

//top line
fprintf(output2, "Predictor A1 A2 Neighbours Tagging BLANK BLANK Categories BLANK");
for(q=0;q<num_parts;q++){fprintf(output2," Model_%d", q+1);}
fprintf(output2,"\n");

//now rest
for(j=0;j<count;j++)
{
fprintf(output2, "%s %c %c %d %.3f -1 -1 %d -1", preds[j], al1[j], al2[j], (int)tally1[j], tally2[j], usedpreds[j]);
for(q=0;q<num_parts;q++){fprintf(output2, " %.4f", tally3[(size_t)q*count+j]);}
fprintf(output2,"\n");
}

for(q2=0;q2<num_parts;q2++)
{
fprintf(output2, "The relative contribution of Category %d to each category", q2+1);
for(q=0;q<num_parts;q++){fprintf(output2, " %.4f", ssums[q][q2]);}
fprintf(output2, "\n");
}

fprintf(output2, "There are %d reference %d regression %d heritability predictors", total, total2, total3);
for(q=0;q<num_parts;q++){fprintf(output2, " %d", (int)ssums[q][num_parts]);}
fprintf(output2, "\n");

fclose(output2);

printf("New tagging file saved in %s\n\n", filename2);

free(pcounts);
for(j=0;j<count;j++){free(preds[j]);}free(preds);
free(al1);free(al2);free(tally1);free(tally2);free(tally3);free(usedpreds);
for(k=0;k<num_parts;k++){free(ssums[k]);}free(ssums);
free(rs);

///////////////////////////

