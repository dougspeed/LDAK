/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Read in details for sections, partitions, genes or windows

///////////////////////////

if(mode==102||mode==103||mode==104)	//read sections
{
sprintf(filename,"%ssection.details",folder);
if(just_check(filename)!=0)
{printf("Error reading %s; this should have been created using \"--cut-weights\"\n\n", filename);exit(1);}

num_sections=countrows(filename)-4;
if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}

//read top line
if(fscanf(input, "%s %s %s ", readstring, readstring2, readstring3)!=3)
{printf("Error reading Row 1 of %s, suggesting the file has been changed since creation with \"--cut-weights\"\n\n", filename);exit(1);}
if(strcmp(readstring,"Datafiles")!=0)
{printf("Error, %s should begin \"Datafiles\" (not %s), suggesting the file has been changed since creation with \"--cut-weights\"\n\n", filename, readstring);exit(1);}

//read second line (Using Thinned Predictors / Filtered and Thinned / Using All Predictors / Using Filtered Predictors)
if(fscanf(input, "%s %s %s ", readstring, readstring2, readstring3)!=3)
{printf("Error reading Row 2 of %s, suggesting the file has been changed since creation with \"--cut-weights\"\n\n", filename);exit(1);}
if(strcmp(readstring2,"Thinned")!=0&&strcmp(readstring2,"and")!=0&&strcmp(readstring2,"All")!=0&&strcmp(readstring2,"Filtered")!=0)
{printf("Error reading Row 2 of %s, suggesting the file has been changed since creation with \"--cut-weights\"\n\n", filename);exit(1);}

//set nothin
if(strcmp(readstring2,"Thinned")==0||strcmp(readstring2,"and")==0){nothin=0;}
if(strcmp(readstring2,"All")==0||strcmp(readstring2,"Filtered")==0){nothin=1;}

if(nothin==1&&spread==1){printf("Error, it is not possible to use \"--spread YES\", because you used \"--no-thin YES\" with \"--cut-weights\"\n\n");exit(1);}

//check extract
if((strcmp(readstring2,"and")==0||strcmp(readstring2,"Filtered")==0)&&extract==0)
{printf("Error, the predictor filterings (\"--extract\", \"--exclude\" or \"--chr\") used with \"--cut-weights\" should also be provided now\n\n");exit(1);}
if((strcmp(readstring2,"Thinned")==0||strcmp(readstring2,"All")==0)&&extract==1)
{printf("Error, the predictor filterings (\"--extract\", \"--exclude\" or \"--chr\") provided now were not used with \"--cut-weights\"\n\n");exit(1);}

//read window_kb, window_length or window_cm
if(fscanf(input, "%s %s %s ", readstring, readstring2, readstring3)!=3)
{printf("Error reading Row 3 of %s\n\n", filename);exit(1);}
if(strcmp(readstring,"Window")!=0||(strcmp(readstring2,"kb")!=0&&strcmp(readstring2,"length")!=0&&strcmp(readstring2,"cM")!=0))
{printf("Error reading Row 3 of %s, suggesting the file has been changed since creation with \"--cut-weights\"\n\n", filename);exit(1);}

if(window_kb==-9999&&window_length==-9999&&window_cm==-9999)	//use these values
{
if(strcmp(readstring2,"kb")==0){window_kb=atof(readstring3);}
if(strcmp(readstring2,"length")==0){window_length=atoi(readstring3);}
if(strcmp(readstring2,"cM")==0){window_cm=atof(readstring3);}
}

//skip header row
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}

sstarts=malloc(sizeof(int)*num_sections);
sends=malloc(sizeof(int)*num_sections);
sstarts2=malloc(sizeof(int)*num_sections);
sends2=malloc(sizeof(int)*num_sections);
for(j=0;j<num_sections;j++)
{
if(fscanf(input,"%d %d-%d %d-%d ", &readint, sstarts+j, sends+j, sstarts2+j, sends2+j)!=5)
{printf("Error reading Row %d of %s\n\n", 5+j, filename);exit(1);}
sstarts[j]--;sstarts2[j]--;
}	//end of j loop

if(section>num_sections)
{printf("Error, the section provided (%d) is higher than the total number (%d)\n\n", section, num_sections);exit(1);}

if(section_start>num_sections)
{printf("Error, the start-section provided (%d) is higher than the total number (%d)\n\n", section_start, num_sections);exit(1);}

fclose(input);
}	//end of reading sections

///////////////////////////

if(mode==112||mode==113)	//read partitions (for mode=113, only need num_parts)
{
sprintf(filename,"%spartition.details",folder);
if(just_check(filename)!=0)
{printf("Error reading %s; this should have been created using \"--cut-kins\"\n\n", filename);exit(1);}

num_parts=countrows(filename)-3;
if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}

if(mode==112)
{
//read top line
if(fscanf(input, "%s %s %s ", readstring, readstring2, readstring3)!=3)
{printf("Error reading Row 1 of %s, suggesting the file has been changed since creation with \"--cut-kins\"\n\n", filename);exit(1);}
if(strcmp(readstring,"Datafiles")!=0)
{printf("Error, %s should begin \"Datafiles\" (not %s), suggesting the file has been changed since creation with \"--cut-kins\"\n\n", filename, readstring);exit(1);}

//read second line (Using All Predictors / Using Filtered Predictors)
if(fscanf(input, "%s %s %s ", readstring, readstring2, readstring3)!=3)
{printf("Error reading Row 2 of %s, suggesting the file has been changed since creation with \"--cut-kins\"\n\n", filename);exit(1);}
if(strcmp(readstring2,"All")!=0&&strcmp(readstring2,"Filtered")!=0)
{printf("Error reading Row 2 of %s, suggesting the file has been changed since creation with \"--cut-kins\"\n\n", filename);exit(1);}

//check extract
if(strcmp(readstring2,"Filtered")==0&&extract==0)
{printf("Error, the predictor filterings (\"--extract\", \"--exclude\" or \"--chr\") used with \"--cut-kins\" should also be provided now\n\n");exit(1);}
if(strcmp(readstring2,"All")==0&&extract==1)
{printf("Error, the predictor filterings (\"--extract\", \"--exclude\" or \"--chr\") provided now were not used with \"--cut-kins\"\n\n");exit(1);}

//skip header row
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}

pstarts=malloc(sizeof(int)*num_parts);
pends=malloc(sizeof(int)*num_parts);
for(q=0;q<num_parts;q++)
{
if(fscanf(input,"%d %s %s ", &readint, readstring, readstring2)!=3)
{printf("Error reading Row %d of %s\n\n", 4+q, filename);exit(1);}
if(strcmp(readstring,"-1")==0)	//using partitions
{
if(q==0)	//set partpref then check
{
count=strlen(readstring2)-1;
if(readstring2[count]!='1')
{printf("Error reading filename on Row %d of %s\n\n", 4, filename);exit(1);}
strcpy(partpref2,readstring2);partpref2[count]='\0';
(void)append_check(partpref,partpref2,workdir);
}

sprintf(readstring3,"%s%d",partpref2,q+1);
if(strcmp(readstring2,readstring3)!=0)
{printf("Error reading filename on Row %d of %s\n\n", 4+q, filename);exit(1);}
pstarts[q]=-9999;pends[q]=-9999;
}
else	//so using predictor numbers
{pstarts[q]=atoi(readstring)-1;pends[q]=atoi(readstring2);}
}	//end of q loop
}

if(partition>num_parts)
{printf("Error, the partition provided (%d) is higher than the total number (%d)\n\n", partition, num_parts);exit(1);}

fclose(input);
}	//end of reading partitions

///////////////////////////

if(mode==137||mode==138||mode==139||mode==140)	//read genes
{
if(mode!=140){sprintf(filename,"%sgenes.details",folder);}
else{sprintf(filename,"%s.genes.details",outfile);}
if(just_check(filename)!=0)
{printf("Error reading %s; this should have been created using \"--cut-genes\"\n\n", filename);exit(1);}

num_genes=countrows(filename)-3;
if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}

//read top line
if(fscanf(input, "%s %s %s ", readstring, readstring2, readstring3)!=3)
{printf("Error reading Row 1 of %s, suggesting the file has been changed since creation with \"--cut-kins\"\n\n", filename);exit(1);}
if(strcmp(readstring,"Datafiles")!=0)
{printf("Error, %s should begin \"Datafiles\" (not %s), suggesting the file has been changed since creation with \"--cut-genes\"\n\n", filename, readstring);exit(1);}

//read second line (Using All Predictors / Using Filtered Predictors)
if(fscanf(input, "%s %s %s ", readstring, readstring2, readstring3)!=3)
{printf("Error reading Row 2 of %s, suggesting the file has been changed since creation with \"--cut-genes\"\n\n", filename);exit(1);}
if(strcmp(readstring2,"All")!=0&&strcmp(readstring2,"Filtered")!=0)
{printf("Error reading Row 2 of %s, suggesting the file has been changed since creation with \"--cut-genes\"\n\n", filename);exit(1);}

if(mode!=139||cut1!=-9999)	//check extract
{
if(strcmp(readstring2,"Filtered")==0&&extract==0)
{printf("Error, the predictor filterings (\"--extract\", \"--exclude\" or \"--chr\") used with \"--cut-genes\" should also be provided now\n\n");exit(1);}
if(strcmp(readstring2,"All")==0&&extract==1)
{printf("Error, the predictor filterings (\"--extract\", \"--exclude\" or \"--chr\") provided now were not used with \"--cut-genes\"\n\n");exit(1);}
}

//skip header row
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}

gnames=malloc(sizeof(char*)*num_genes);
gstarts=malloc(sizeof(int)*num_genes);
gends=malloc(sizeof(int)*num_genes);
gparts=malloc(sizeof(int)*num_genes);
gchr=malloc(sizeof(int)*num_genes);
gbp1=malloc(sizeof(double)*num_genes);
gbp2=malloc(sizeof(double)*num_genes);
gpvas=malloc(sizeof(double)*num_genes);

for(j=0;j<num_genes;j++)
{
if(fscanf(input,"%s %d %d %d %s %d %lf %lf ", readstring, gstarts+j, gends+j, gparts+j, readstring2, gchr+j, gbp1+j, gbp2+j)!=8)
{printf("Error reading Row %d of %s\n\n", 4+j, filename);exit(1);}
copy_string(gnames,j,readstring);
gstarts[j]--;
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}
}	//end of j loop

num_parts=gparts[num_genes-1];

if(partition>num_parts)
{printf("Error, the partition provided (%d) is higher than the total number (%d)\n\n", partition, num_parts);exit(1);}

fclose(input);
}	//end of reading genes

///////////////////////////

if(mode==192||mode==193||mode==194)	//read partitions
{
sprintf(filename,"%spartition.details",folder);
if(just_check(filename)!=0)
{printf("Error reading %s; this should have been created using \"--cut-gre\"\n\n", filename);exit(1);}

num_parts=countrows(filename)-4;
if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}

//read top line
if(fscanf(input, "%s %s %s ", readstring, readstring2, readstring3)!=3)
{printf("Error reading Row 1 of %s, suggesting the file has been changed since creation with \"--cut-gre\"\n\n", filename);exit(1);}
if(strcmp(readstring,"Datafiles")!=0)
{printf("Error, %s should begin \"Datafiles\" (not %s), suggesting the file has been changed since creation with \"--cut-gre\"\n\n", filename, readstring);exit(1);}

//read second line (Using All Predictors / Using Filtered Predictors)
if(fscanf(input, "%s %s %s ", readstring, readstring2, readstring3)!=3)
{printf("Error reading Row 2 of %s, suggesting the file has been changed since creation with \"--cut-gre\"\n\n", filename);exit(1);}
if(strcmp(readstring2,"All")!=0&&strcmp(readstring2,"Filtered")!=0)
{printf("Error reading Row 2 of %s, suggesting the file has been changed since creation with \"--cut-gre\"\n\n", filename);exit(1);}

//check extract
if(strcmp(readstring2,"Filtered")==0&&extract==0)
{printf("Error, the predictor filterings (\"--extract\", \"--exclude\" or \"--chr\") used with \"--cut-gre\" should also be provided now\n\n");exit(1);}
if(strcmp(readstring2,"All")==0&&extract==1)
{printf("Error, the predictor filterings (\"--extract\", \"--exclude\" or \"--chr\") provided now were not used with \"--cut-gre\"\n\n");exit(1);}

//read third line (Using All Samples / Using Filtered Samples)
if(fscanf(input, "%s %s %s ", readstring, readstring2, readstring3)!=3)
{printf("Error reading Row 3 of %s, suggesting the file has been changed since creation with \"--cut-gre\"\n\n", filename);exit(1);}
if(strcmp(readstring2,"All")!=0&&strcmp(readstring2,"Filtered")!=0)
{printf("Error reading Row 3 of %s, suggesting the file has been changed since creation with \"--cut-gre\"\n\n", filename);exit(1);}

//check filtering
if(strcmp(readstring2,"Filtered")==0&&strcmp(bsampfile,"blank")==0&&strcmp(csampfile,"blank")==0)
{printf("Error, the sample filterings (\"--keep\", \"--remove\" or \"--pheno\") used with \"--cut-gre\" should also be provided now\n\n");exit(1);}
if(strcmp(readstring2,"All")==0&&(strcmp(bsampfile,"blank")!=0||strcmp(csampfile,"blank")!=0))
{printf("Error, the predictor filterings (\"--keep\", \"--remove\" or \"--pheno\") provided now were not used with \"--cut-gre\"\n\n");exit(1);}

//skip header row
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}

pstarts=malloc(sizeof(int)*num_parts);
pends=malloc(sizeof(int)*num_parts);
for(q=0;q<num_parts;q++)
{
if(fscanf(input,"%d %d %d ", &readint, pstarts+q, pends+q)!=3)
{printf("Error reading Row %d of %s\n\n", 5+q, filename);exit(1);}
pstarts[q]--;
}	//end of q loop

if(partition>num_parts)
{printf("Error, the partition provided (%d) is higher than the total number (%d)\n\n", partition, num_parts);exit(1);}

fclose(input);
}	//end of reading partitions

///////////////////////////

