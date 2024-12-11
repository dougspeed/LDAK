/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Check and read details from the prs root file - here for modes 131, 132, 134 and 140 (134 becomes either 131 or 132)

///////////////////////////

//check and open root file

sprintf(filename,"%s.root", prsfile);
if(just_check(filename)!=0)
{
if(mode==131||mode==132){printf("Error reading %s; this file would have been created using \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\" with \"--LOCO YES\"\n\n", filename);}
else{printf("Error reading %s; this file would have been created using \"--kvik-step1\", \"--GCTA-LOCO-step1\" or \"--fastGWA-step1\"\n\n", filename);}
exit(1);
}

count=countrows(filename);
if(count!=11)
{printf("Error, %s should have 11 rows (not %d), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n", filename, count);exit(1);}

if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}

//read analysis type

if(fscanf(input, "%s %s ", readstring, readstring2)!=2)
{printf("Error reading Row 1 of %s\n\n", filename);exit(1);}
if(strcmp(readstring,"Analysis")!=0)
{printf("Error reading %s; Row 1 should begin \"Analysis\" (not %s), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n",filename, readstring);exit(1);}

fastgwa=-1;
if(strcmp(readstring2,"KVIK")==0){fastgwa=0;}
if(strcmp(readstring2,"fastGWA")==0){fastgwa=1;}
if(strcmp(readstring2,"KVIK_Pedigree")==0){fastgwa=2;}

if(fastgwa==-1){printf("Error reading %s; Element 2 of Row 1 should be \"KVIK\", \"KVIK_Pedigree\" or \"fastGWA\" (not %s), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n", filename, readstring2);exit(1);}

//read regression model

if(fscanf(input, "%s %s ", readstring, readstring2)!=2)
{printf("Error reading Row 2 of %s\n\n", filename);exit(1);}
if(strcmp(readstring,"Regression")!=0)
{printf("Error reading %s; Row 2 should begin \"Regression\" (not %s), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n",filename, readstring);exit(1);}

if(mode==131)	//check consistent
{
if(strcmp(readstring2,"Logistic")==0){printf("Error, \"--binary YES\" was used when making the PRS, so you should now use \"--logistic\" (instead of \"--linear\")\n\n");exit(1);}
}
if(mode==132)	//check consistent
{
if(strcmp(readstring2,"Linear")==0){printf("Error, \"--binary YES\" was not used when making the PRS, so you should now use \"--linear\" (instead of \"--logistic\")\n\n");exit(1);}
}
if(mode==134)	//switch to 131 or 132
{
if(strcmp(readstring2,"Linear")==0){mode=131;}
else{mode=132;}
}

//read data, phenotypes, covariates, tops and factors

if(fscanf(input, "%s %s ", readstring, readstring2)!=2)
{printf("Error reading Row 3 of %s\n\n", filename);exit(1);}
if(strcmp(readstring,"Datafile")!=0)
{printf("Error reading %s; Row 3 should begin \"Datafile\" (not %s), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n",filename, readstring);exit(1);}

if(fscanf(input, "%s %s ", readstring, readstring3)!=2)
{printf("Error reading Row 4 of %s\n\n", filename);exit(1);}
if(strcmp(readstring,"Phenotypes")!=0)
{printf("Error reading %s; Row 4 should begin \"Phenotypes\" (not %s), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n",filename, readstring);exit(1);}

if(fscanf(input, "%s %s ", readstring, readstring4)!=2)
{printf("Error reading Row 5 of %s\n\n", filename);exit(1);}
if(strcmp(readstring,"Covariates")!=0)
{printf("Error reading %s; Row 5 should begin \"Covariates\" (not %s), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n",filename, readstring);exit(1);}

if(fscanf(input, "%s %s ", readstring, readstring5)!=2)
{printf("Error reading Row 6 of %s\n\n", filename);exit(1);}
if(strcmp(readstring,"Top_Predictors")!=0)
{printf("Error reading %s; Row 6 should begin \"Top_Predictors\" (not %s), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n",filename, readstring);exit(1);}

if(fscanf(input, "%s %s ", readstring, readstring6)!=2)
{printf("Error reading Row 7 of %s\n\n", filename);exit(1);}
if(strcmp(readstring,"Factors")!=0)
{printf("Error reading %s; Row 7 should begin \"Factors\" (not %s), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n",filename, readstring);exit(1);}

if(mode==131||mode==132)	//check covariates are consistent
{
if(strcmp(covarfile,"blank")==0&&strcmp(readstring4,"none")!=0)
{printf("Error, \"--covar\" was used when making the PRS, but is not used now\n\n");exit(1);}
if(strcmp(covarfile,"blank")!=0&&strcmp(readstring4,"none")==0)
{printf("Error, \"--covar\" is used now, but was not used when making the PRS\n\n");exit(1);}

if(strcmp(topfile,"blank")==0&&strcmp(readstring5,"none")!=0)
{printf("Error, \"--top-preds\" was used when making the PRS, but is not used now\n\n");exit(1);}
if(strcmp(topfile,"blank")!=0&&strcmp(readstring5,"none")==0)
{printf("Error, \"--top-preds\" is used now, but was not used when making the PRS\n\n");exit(1);}

if(strcmp(factorfile,"blank")==0&&strcmp(readstring6,"none")!=0)
{printf("Error, \"--factors\" was used when making the PRS, but is not used now\n\n");exit(1);}
if(strcmp(factorfile,"blank")!=0&&strcmp(readstring6,"none")==0)
{printf("Error, \"--factors\" is used now, but was not used when making the PRS\n\n");exit(1);}
}

if(checkroot==1)	//check names consistent (most checks only required for linear / logistic)
{
flag=0;
if(strcmp(datafile,readstring2)!=0)
{printf("Error, the data file used when making the PRS (%s) does not match that used now (%s)\n\n", readstring2, datafile);flag=1;}

if(mode==131||mode==132)
{
if(strcmp(respfile,readstring3)!=0)
{printf("Error, the phenotypes file used when making the PRS (%s) does not match that used now (%s)\n\n", readstring3, respfile);flag=1;}
if(strcmp(covarfile,"blank")!=0&&strcmp(covarfile,readstring4)!=0)
{printf("Error, the covariates file used when making the PRS (%s) does not match that used now (%s)\n\n", readstring4, covarfile);flag=1;
if(strcmp(topfile,"blank")!=0&&strcmp(topfile,readstring5)!=0)
{printf("Error, the top predictors file used when making the PRS (%s) does not match that used now (%s)\n\n", readstring5, topfile);flag=1;}
if(strcmp(factorfile,"blank")!=0&&strcmp(factorfile,readstring6)!=0)
{printf("Error, the factors file used when making the PRS (%s) does not match that used now (%s)\n\n", readstring6, factorfile);flag=1;}
}
}
if(flag==1){printf("If you are sure you have provided the correct files, use \"--check-root NO\"\n\n");exit(1);}
}

//read phenotype number

if(fscanf(input, "%s %s ", readstring, readstring2)!=2)
{printf("Error reading Row 8 of %s\n\n", filename);exit(1);}
if(strcmp(readstring,"Phenotype_Number")!=0&&strcmp(readstring,"Num_Phenotypes")!=0)
{printf("Error reading %s; Row 8 should begin \"Phenotype_Number\" or \"Num_Phenotypes\" (not %s), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n",filename, readstring);exit(1);}

if(strcmp(readstring,"Phenotype_Number")==0)	//set kvikparity to zero and check mpheno
{
kvikparity=0;

if(mpheno!=-9999&&mpheno!=atoi(readstring2))
{printf("Error, the phenotype number provided now (%d) does not match that used with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\" (%d)\n\n", mpheno, atoi(readstring2));exit(1);}
}
else	//set kvikparity to the number of phenotypes and check mpheno
{
kvikparity=atoi(readstring2);

if(mpheno>0)	//have specified a phenotype
{
if(mpheno>atoi(readstring2))
{printf("Error, the phenotype number provided now (%d) is higher than the number of phenotypes analyzed with  \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\" (%d)\n\n", mpheno, kvikparity);exit(1);}
}
}

//skip two rows
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}

//read num_chr (save in num_chr2)

if(fscanf(input, "%s %d ", readstring, &num_chr2)!=2)
{printf("Error reading Row 11 of %s\n\n", filename);exit(1);}
if(strcmp(readstring,"Num_Chromosomes")!=0)
{printf("Error reading %s; Row 11 should begin \"Num_Chromosomes\" (not %s), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n", filename, readstring);exit(1);}

fclose(input);

////////

if(mode==140)	//get power (stored in Row 5 of details file), set pvafile and sumsfile and check files exist
{
if(kvikparity==0)	//previously analyzed just one phenotype
{
sprintf(filename,"%s.loco.details", prsfile);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--kvik-step1\", \"--GCTA-LOCO-step1\" or \"--fastGWA-step1\"\n\n", filename);exit(1);}

count=countrows(filename);
if(count!=6)
{printf("Error, %s should have six rows (not %d), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n", filename, count);exit(1);}

//open and skip four rows
if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}

//read power and close file
if(fscanf(input, "%s %lf ", readstring, &power)!=2)
{printf("Error reading Row 5 of %s\n\n", filename);exit(1);}
if(strcmp(readstring,"Power")!=0)
{printf("Error reading %s; Row 5 should begin \"Power\" (not %s), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n", filename, readstring);exit(1);}
fclose(input);

//outfile will end step3 - change this to step2
strcpy(filename,outfile);
filename[strlen(filename)-1]='2';

sprintf(pvafile,"%s.pvalues",filename);
sprintf(sumsfile,"%s.summaries",filename);
}
else	//previously analyzed multiple phenotypes
{
if(mpheno==-9999)
{
if(kvikparity==1){mpheno=1;}
else{printf("Error, you must use \"--mpheno\" to specify which of the %d phenotypes to analyze\n\n", kvikparity);exit(1);}
}

sprintf(filename,"%s.pheno%d.loco.details", prsfile, mpheno);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--kvik-step1\", \"--GCTA-LOCO-step1\" or \"--fastGWA-step1\"\n\n", filename);exit(1);}

count=countrows(filename);
if(count!=6)
{printf("Error, %s should have six rows (not %d), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n", filename, count);exit(1);}

//open and skip four rows
if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}

//read power and close file
if(fscanf(input, "%s %lf ", readstring, &power)!=2)
{printf("Error reading Row 5 of %s\n\n", filename);exit(1);}
if(strcmp(readstring,"Power")!=0)
{printf("Error reading %s; Row 5 should begin \"Power\" (not %s), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n", filename, readstring);exit(1);}
fclose(input);

//outfile will end step3 - change this to step2
strcpy(filename,outfile);
printf("filename %s\n",filename);
filename[strlen(filename)-1]='2';
printf("now %s\n", filename);

sprintf(pvafile,"%s.pheno%d.pvalues",filename, mpheno);
sprintf(sumsfile,"%s.pheno%d.summaries",filename, mpheno);
}

if(just_check(pvafile)!=0)
{printf("Error reading %s; this file should have been created using \"--kvik-step2\"\n\n", pvafile);exit(1);}

if(just_check(sumsfile)!=0)
{printf("Error reading %s; this file should have been created using \"--kvik-step2\"\n\n", sumsfile);exit(1);}
}

///////////////////////////

