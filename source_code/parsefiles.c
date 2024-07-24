/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Process and check files

///////////////////////////

//data - convenient to set use_data, which states whether using datafiles
//0 - no data, 1 - "normal", 2 - only through regions, 3 - only through top preds
//4 - have data, but only using ids, 5 - multiple datafiles, 6 - have bim (not using ids)
use_data=0;

if(mode==101||mode==102||mode==103||mode==104||mode==105||mode==106||mode==107||mode==108||mode==109||mode==111||mode==112||mode==114||(mode==117&&extract==1)||mode==122||mode==127||mode==128||mode==131||mode==132||mode==136||mode==137||mode==138||mode==140||(mode==139&&cut1!=-9999)||mode==141||mode==145||mode==151||mode==152||mode==153||mode==154||mode==156||mode==158||mode==160||mode==162||mode==171||mode==172||mode==173||mode==175||mode==186||mode==187||mode==188||mode==189||mode==191||mode==192||mode==193||mode==194)
{use_data=1;}

if(use_data==0&&num_regs>0)
{use_data=2;}

if(use_data==0&&strcmp(topfile,"blank")!=0)
{use_data=3;}

if(mode==177&&dtype!=-9999)
{use_data=4;}

if(mode==181||mode==182||mode==183||mode==184||mode==185||mode==190)
{use_data=5;}

if(mode==159)
{use_data=6;}

if((use_data==1||use_data==2||use_data==3)&&dtype==-9999)
{
if(mode==151||mode==152||mode==153||mode==154)
{printf("Error, you must provide a set of genetic data files using \"--bfile\" or \"--speed\"\n\n");exit(1);}

printf("Error, you must provide a set of genetic data files using \"--bfile\", \"--sp\", \"--sped\", \"--speed\" or \"--gen\"\n\n");exit(1);
}

if(use_data==5&&dtype==-9999)
{
if(mode==190){printf("Error, you must provide two sets of genetic data files using \"--mbfile\", \"--msp\", \"--msped\" or \"--mspeed\"\n\n");exit(1);}
printf("Error, you must provide one or more sets of genetic data files using \"--bfile\", \"--sp\", \"--sped\", \"--speed\", \"--gen\", \"--mbfile\", \"--msp\", \"--msped\" or \"--mgen\"\n\n");exit(1);
}

if(use_data==0&&dtype!=-9999)
{
if(strcmp(weightsfile,"blank")!=0)
{printf("Warning, genetic data files and weightings are provided but will not be used\n\n");}
else
{printf("Warning, genetic data files are provided but will not be used\n\n");}
dtype=-9999;strcpy(weightsfile,"blank");
}

if((use_data==3||use_data==4)&&strcmp(weightsfile,"blank")!=0)
{printf("Warning, weightings are provided but will not be used\n\n");strcpy(weightsfile,"blank");}

if(extract==1)
{
if(strcmp(topfile,"blank")!=0)
{printf("Warning, the predictor filterings will not apply to top predictors\n\n");}
if(strcmp(targetfile,"blank")!=0)
{printf("Warning, the predictor filterings will not apply to target predictors\n\n");}
}

////////

if(famhead==-9999)
{
if(dtype==5)	//then must be providing a list (else famhead would be set)
{
count=countcols(datalist);
{printf("Warning, it is assumed that Column %d of %s provides sample files; if these are instead fam files, please use \"--fam-files YES\"\n\n", count, datalist);}
famhead=2;
}
else{famhead=0;}
}

if(use_data==1||use_data==2||use_data==3||use_data==4||use_data==5)	//get stems - dtype must be 1-5 (dtypes 11-15 will have become dtype 5)
{
if(strcmp(datalist,"blank")==0)	//single dataset
{
num_files=1;
datastems=malloc(sizeof(char *));datastems[0]=malloc(sizeof(char)*500);
bimstems=malloc(sizeof(char *));bimstems[0]=malloc(sizeof(char)*500);
famstems=malloc(sizeof(char *));famstems[0]=malloc(sizeof(char)*500);

if(dtype!=5){fill_names(datastems, bimstems, famstems, 0, udatafile, dtype);}
else
{
strcpy(datastems[0],udatafile);
if(strcmp(ubimfile,"blank")!=0){strcpy(bimstems[0],ubimfile);}
else{strcpy(bimstems[0],udatafile);}
strcpy(famstems[0],ufamfile);
}
}
else	//reading from list
{
num_files=countrows(datalist);

count=countcols(datalist);
if(countels(datalist)!=num_files*count)
{printf("Error, each row of %s should contain the same number of files\n\n", datalist);exit(1);}

if(dtype==1&&count!=1&&count!=3)
{printf("Error, each row of %s must contain either one or three elements (not %d), providing either the data stem or bed, bim and fam files\n\n",datalist,count);exit(1);}
if(dtype==3&&count!=1&&count!=3)
{printf("Error, each row of %s must contain either one or three elements (not %d), providing either the data stem or sped, bim and fam files\n\n",datalist,count);exit(1);}
if(dtype==4&&count!=1&&count!=3)
{printf("Error, each row of %s must contain either one or three elements (not %d), providing either the data stem or speed, bim and fam files\n\n",datalist,count);exit(1);}

if(dtype==5)
{
if(genheaders<4&&count!=3)
{printf("Error, each row of %s must contain three elements (not %d), providing gen, bim and sample/fam files\n\n", datalist, count);exit(1);}

if(count!=2&&count!=3)
{printf("Error, each row of %s must contain two or three elements (not %d), providing gen and sample/fam files or gen, bim and sample/fam files\n\n", datalist, count);exit(1);}

if(count==2&&oxchr==-9999)
{printf("Error, either %s should include bim files, or you should use \"--oxford-single-chr\"\n\n", datalist);exit(1);}
}

datastems=malloc(sizeof(char *)*num_files);
bimstems=malloc(sizeof(char *)*num_files);
famstems=malloc(sizeof(char *)*num_files);

if((input=fopen(datalist,"r"))==NULL)
{printf("Error opening %s\n\n",datalist);exit(1);}

for(k=0;k<num_files;k++)
{
datastems[k]=malloc(sizeof(char)*500);
bimstems[k]=malloc(sizeof(char)*500);
famstems[k]=malloc(sizeof(char)*500);

if(count==1)	//cant be dtype=5
{
if(fscanf(input, "%s ", readstring2)!=1)
{printf("Error reading stem for Dataset %d from %s\n\n", k+1, datalist);exit(1);}
(void)append_check(readstring,readstring2,workdir);
fill_names(datastems, bimstems, famstems, k, readstring, dtype);
}
if(count==2)	//must be dtype=5
{
if(fscanf(input, "%s ", readstring2)!=1)
{printf("Error reading gen file for Dataset %d from %s\n\n", k+1, datalist);exit(1);}
(void)append_check(readstring,readstring2,workdir);
strcpy(datastems[k],readstring);
strcpy(bimstems[k],datastems[k]);

if(fscanf(input, "%s ", readstring2)!=1)
{printf("Error reading fam/sample file for Dataset %d from %s\n\n", k+1, datalist);exit(1);}
(void)append_check(readstring,readstring2,workdir);
strcpy(famstems[k],readstring);
}
if(count==3)	//can be any
{
if(fscanf(input, "%s ", readstring2)!=1)
{printf("Error reading predictor file for Dataset %d from %s\n\n", k+1, datalist);exit(1);}
(void)append_check(readstring,readstring2,workdir);
strcpy(datastems[k],readstring);
if(fscanf(input, "%s ", readstring2)!=1)
{printf("Error reading bim file for Dataset %d from %s\n\n", k+1, datalist);exit(1);}
(void)append_check(readstring,readstring2,workdir);
strcpy(bimstems[k],readstring);
if(fscanf(input, "%s ", readstring2)!=1)
{printf("Error reading fam/sample file for Dataset %d from %s\n\n", k+1, datalist);exit(1);}
(void)append_check(readstring,readstring2,workdir);
strcpy(famstems[k],readstring);
}

for(k2=0;k2<k;k2++)	//check predictor file names different to previous
{
if(strcmp(datastems[k],datastems[k2])==0)
{printf("Error reading %s; the file %s appears twice\n\n", datalist, datastems[k]);exit(1);}
}
}	//end of k loop

fclose(input);
}	//end of using datalist
}	//end of use_data!=0

////////

if(mode==190&&num_files!=2)
{printf("Error, when using \"--calc-sim-data\", you must provide a pair of genetic data files using \"--mbfile\", \"--msped\", \"--mspeed\" or \"--mgen\"\n\n");exit(1);}

if(mode==181||mode==182||mode==183||mode==184||mode==185)	//check names different to outfile
{
if(mode==181){sprintf(filename,"%s.bed",outfile);}
if(mode==182){sprintf(filename,"%s.sp",outfile);}
if(mode==183){sprintf(filename,"%s.sped",outfile);}
if(mode==184){sprintf(filename,"%s.speed",outfile);}
if(mode==185){sprintf(filename,"%s.gen",outfile);}
for(k=0;k<num_files;k++)
{
if(strcmp(datastems[k],filename)==0)
{
if(num_files==1)
{printf("Error, the name of the new predictor file (%s) matches that provided\n\n", filename);exit(1);}
else
{printf("Error, the name of the new predictor file (%s) matches that of Dataset %d\n\n", filename, k+1);exit(1);}
}}
}

////////

//check datafiles exist (and are not folders)

for(k=0;k<num_files;k++)
{
if(just_check(datastems[k])!=0)
{printf("Error reading %s; check file exists and is not empty\n\n", datastems[k]);exit(1);}
if(just_check(bimstems[k])!=0)
{printf("Error reading bim file %s; check file exists and is not empty\n\n", bimstems[k]);exit(1);}
if(just_check(famstems[k])!=0)
{printf("Error reading fam/sample file %s; check file exists and is not empty\n\n", famstems[k]);exit(1);}
}

//set datafile, bimfile and famfile (used when there is one set of data files)
if(use_data==0||use_data==5)
{strcpy(datafile,"blank");strcpy(bimfile,"blank");strcpy(famfile,"blank");}
if(use_data==1||use_data==2||use_data==3)
{strcpy(datafile,datastems[0]);strcpy(bimfile,bimstems[0]);strcpy(famfile,famstems[0]);}
if(use_data==4)
{strcpy(datafile,"blank");strcpy(bimfile,"blank");strcpy(famfile,famstems[0]);}
if(use_data==6)	//do not have fam file - currently must be mode 159
{
strcpy(datafile,"blank");sprintf(bimfile,"%s.cors.bim", corname);strcpy(famfile,"blank");
}

////////

if(dtype==5)	//give some advice
{
if(num_files==1){printf("Will be expecting the predictor file to have ");}
else{printf("Will be expecting each predictor file to have ");}
if(genskip==0){printf("no header row, ");}
if(genskip==1){printf("one header row, ");}
if(genskip>1){printf("%d header rows, ", genskip);}
if(genheaders==0){printf("then for each row to have no header columns ");}
if(genheaders==1){printf("then for each row to have one header column (which is ignored), ");}
if(genheaders==2){printf("then for each row to have two header columns (which are ignored), ");}
if(genheaders==3){printf("then for each row to have three header columns (Predictor, A1, A2), ");}
if(genheaders==4){printf("then for each row to have four header columns (Predictor, BP, A1, A2), ");}
if(genheaders==5){printf("then for each row to have five header columns (IGNORED, Predictor, BP, A1, A2), ");}
if(genprobs==0){printf("and contain two binary values (\"0 0\", \"1 0\", \"0 1\" or \"1 1\") per sample");}
if(genprobs==1&&nonsnp==0){printf("and contain one dosage value (expected allele count) per sample");}
if(genprobs==1&&nonsnp==1){printf("and contain one value per sample");}
if(genprobs==2){printf("and contain two probabilities (for A1A1 and A1A2) for each sample");}
if(genprobs==3){printf("and contain three probabilities (for A1A1, A1A2 and A2A2) for each sample");}
if(genprobs==4){printf("and contain four probabilities (for A1A1, A1A2, A2A2 and missing) for each sample; ");}
printf("; these can be changed using \"--gen-skip\", \"--gen-headers\" and \"--gen-probs\"");
if(num_files==1){printf(" (if a bim file is provided, predictor details will be read from this instead)\n\n");}
else{printf(" (if bim files are provided, predictor details will be read from these instead)\n\n");}
}

if((mode==151||mode==152||mode==153||mode==154)&&dtype==4)
{
if((input=fopen(datafile,"rb"))==NULL)
{printf("Error opening %s\n\n",datafile);exit(1);}
if(read_speed_size(filename)==2)
{printf("Error, you can not use long SPEED format with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\" (when you make %s, you should not add \"--speed-long YES\")\n\n", datafile);exit(1);}
fclose(input);
}

///////////////////////////

//data filtering

if(num_subs>0)	//check exist
{
for(s=0;s<num_subs;s++)
{
sprintf(filename,"%s%d",subpref, s+1);
if(just_check(filename)!=0)
{printf("Error reading %s; check file exists and is not empty\n\n", filename);exit(1);}
}
}

///////////////////////////

//data scaling (and pvalues) and coding

if(strcmp(centresfile,"blank")!=0)	//check correct size
{
if(countcols(centresfile)!=4)
{printf("Error, %s should have four columns, providing the predictor names, the A1 and A2 alleles, then the mean values for the A1 alleles (not %d)\n\n", centresfile, countcols(centresfile));exit(1);}
}
	
if(use_data==1||use_data==2||use_data==6)	//will use the weights vector
{
if(strcmp(weightsfile,"blank")==0&&ignoreweights!=1)	//weights not provided
{
if(mode==105||mode==112||mode==114||mode==127||mode==128||mode==136||mode==138||mode==140||mode==141||mode==145||mode==151||mode==152||mode==153||mode==154||mode==159||mode==173||mode==186||mode==187||mode==188||mode==189)	//weights are used - but fine to set to one
{
//printf("Warning, the predictor weightings have been set to one (equivalent to adding \"--ignore-weights YES\"); if you wish to specify different weightings, use \"--weights\"\n\n");
ignoreweights=1;
}
else	//weightings are not used, so set to one
{ignoreweights=1;}
}
}

if(strcmp(weightsfile,"blank")!=0)	//check correct size
{
count=countcols(weightsfile);
if(count!=2&&count!=6)
{printf("Error, %s should either have two or six columns (not %d), with Columns 1 and 2 providing predictor names then weightings\n\n", weightsfile, count);exit(1);}

if(count==6)
{
if((input=fopen(weightsfile,"r"))==NULL)
{printf("Error opening %s\n\n",weightsfile);exit(1);}
if(fscanf(input, "Predictor %s ", readstring)!=1)
{printf("Error reading Row 1 of %s; should begin \"Predictor\" (not %s), suggesting the file has been changed since creation with \"--join-weights\" or \"--calc-weights-all\")\"\n\n", weightsfile, readstring);exit(1);}
fclose(input);
}
}

if(strcmp(pvafile,"blank")!=0)	//check correct size
{
if(countcols(pvafile)!=2)
{printf("Error, %s should have two columns (not %d), providing Predictor names then p-values\n\n", pvafile, countcols(pvafile));exit(1);}
}

if(strcmp(impfile,"blank")!=0)	//check correct size
{
if(countcols(impfile)!=2)
{printf("Error, %s should have two columns (not %d), providing Predictor names then importances\n\n", impfile, countcols(impfile));exit(1);}
}

if(power==-9999)
{
if(mode==112||mode==114||mode==127||mode==128||mode==137||mode==138||mode==141||mode==145||mode==159||mode==173||num_regs>0)
{
if(hwestand==1){printf("Error, you must use \"--power\" to specify the assumed relationship between heritability and allele frequency; for human data we recommend -0.25\n\n");exit(1);}
else{printf("Error, you must use \"--power\" to specify the assumed relationship between heritability and predictor variance; for human data we recommend -0.25\n\n");exit(1);}
}

if(mode==172)
{printf("Error, you must use \"--power\" to specify how to standardize predictors; use 0 if your score file contains raw effect sizes (e.g., those from \"--calc-blups\", \"--linear\", \"--ridge\", \"--bolt\", \"--bayesr\", \"--elastic\", \"--mega-prs\" or \"--calc-pca-loads\"), but use -1 if it contains standardized effect sizes (e.g., those from \"--calc-posts\")\n\n");exit(1);}

if(mode==186||mode==187||mode==188||mode==189)
{printf("Error, you must use \"--power\" to specify how to standardize predictors; for no standardization (i.e., to simply sum allele counts) use \"--power 0\"\n\n");exit(1);}
}
else
{
if(mode!=112&&mode!=114&&mode!=127&&mode!=128&&mode!=137&&mode!=138&&mode!=140&&mode!=141&&mode!=145&&mode!=151&&mode!=152&&mode!=153&&mode!=154&&mode!=159&&mode!=172&&mode!=173&&mode!=186&&mode!=187&&mode!=188&&mode!=189&&num_regs==0)
{printf("Warning, \"--power\" was provided but will not be used\n\n");power=-9999;}
}

///////////////////////////

//kinships - set num_kins and kinship stems, then do some checks

if(strcmp(kinname,"blank")==0&&strcmp(kinlist,"blank")==0){num_kins=0;}
else
{
if(strcmp(kinname,"blank")!=0)	//will have already been prefixed by folder
{
num_kins=1;
kinstems=malloc(sizeof(char *));
kinstems[0]=malloc(sizeof(char)*500);
strcpy(kinstems[0], kinname);
}
else	//read from file
{
num_kins=countrows(kinlist);
kinstems=malloc(sizeof(char *)*num_kins);

if((input=fopen(kinlist,"r"))==NULL)
{printf("Error opening %s\n\n",kinlist);exit(1);}
for(k=0;k<num_kins;k++)
{
kinstems[k]=malloc(sizeof(char)*500);
if(fscanf(input, "%s ", filename)!=1)
{printf("Error reading Row %d of %s\n\n", k+1, kinlist);exit(1);}
(void)append_check(kinstems[k],filename,workdir);

for(k2=0;k2<k;k2++)	//check different to previous
{
if(strcmp(kinstems[k],kinstems[k2])==0)
{printf("Error reading %s; the stem %s appears twice\n\n", kinlist, kinstems[k]);exit(1);}
}
}
fclose(input);
}	//end of reading from file
}

////////

if((mode==115||mode==118||mode==119||mode==161||mode==162||mode==163||mode==164||mode==166||mode==167||mode==168||mode==169||mode==170||mode==201||mode==202||mode==203)&&num_kins>1)
{printf("Error, you can only provide one kinship matrix (not %d)\n\n", num_kins);exit(1);}

if(mode==117&&extract==1&&num_kins>1)
{printf("Error, when using \"--extract\" or \"--exclude\" to subtract Predictors from a kinship matrix, you can only provide one kinship matrix (not %d)\n\n", num_kins);exit(1);}

if(mode==120&&num_kins==1)
{printf("Error, to use \"--calc-sim-grms\" you must provide two or more kinship matrices\n\n");exit(1);}

if(mode==123&&(minrel!=-9999||maxrel!=-9999)&&num_kins!=1)
{printf("Error, you can only use \"--min-rel\" or \"--max-rel\" with \"--pcgc\" if using one kinship matrix (not %d)\n\n", num_kins);exit(1);}

if(mode==124&&(minrel!=-9999||maxrel!=-9999)&&num_kins!=1)
{printf("Error, you can only use \"--min-rel\" or \"--max-rel\" with \"--he\" if using one kinship matrix (not %d)\n\n", num_kins);exit(1);}

if(mode==131&&num_kins>1)
{printf("Error, with \"--linear\" you can provide at most one kinship matrix (not %d); for complex mixed-model linear regression you should instead perform the analysis in two steps: first use\"--solve-null\" to construct a combined kinship matrix, then use this matrix with \"--linear\"\n\n", num_kins);exit(1);}

if(mode==133&&num_kins==1)
{printf("Error, to use \"--solve-null\" you must provide multiple kinship matrices (otherwise, use \"--linear\")\n\n");exit(1);}

if(mode==138&&num_kins>1)
{printf("Error, with \"--calc-genes-reml\" you can provide at most one kinship matrix (not %d); for complex gene/chunk-based REML you should instead perform the analysis in two steps: first use\"--solve-null\" to construct a combined kinship matrix, then use this matrix with \"--calc-genes-reml\"\n\n", num_kins);exit(1);}

if(strcmp(eigenfile,"blank")!=0&&num_kins>1)
{printf("Error, you can only use \"--eigen\" if providing one kinship matrix (not %d)\n\n", num_kins);exit(1);}

if(mode==116||mode==117||mode==118||mode==119||mode==133||mode==164||mode==167||mode==168)	//check names different to outfile
{
for(k=0;k<num_kins;k++)
{
if(strcmp(kinstems[k],outfile)==0)
{
if(num_kins==1){printf("Error, the stem of the new kinship matrix (%s) matches that provided\n\n", outfile);exit(1);}
else{printf("Error, the stem of the new kinship matrix (%s) matches that of Kinship %d\n\n", outfile, k+1);exit(1);}
}
}
}

////////

if(num_kins>0)	//check required kinship files exists and set kinnums and kinsums
{
kinnums=malloc(sizeof(int)*num_kins);
kinsums=malloc(sizeof(double)*num_kins);

flag=0;	//records how many low-rank kinship matrices
for(k=0;k<num_kins;k++)
{
sprintf(filename,"%s.grm.id",kinstems[k]);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created when calculating kinships\n\n", filename);exit(1);}

if(mode!=118&&mode!=119)	//not converting
{
sprintf(filename,"%s.grm.bin",kinstems[k]);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created when calculating kinships\n\n", filename);exit(1);}
}

if(mode==118)	//convert raw
{
sprintf(filename,"%s.grm.gz",kinstems[k]);
if(just_check(filename)!=0)
{printf("Error reading %s; check file exists and is not empty\n\n", filename);exit(1);}
}

if(mode==119)	//convert raw
{
sprintf(filename,"%s.grm.raw",kinstems[k]);
if(just_check(filename)!=0)
{printf("Error reading %s; check file exists and is not empty\n\n", filename);exit(1);}
}

if(kindetails==1&&mode!=118&&mode!=119)	//check details and adjust exist, then read things from adjust
{
sprintf(filename,"%s.grm.details",kinstems[k]);
if(just_check(filename)!=0)
{
printf("Error reading %s; this file would have been created when calculating kinships; ", filename);
if(mode==121||mode==126)
{printf("you can ignore this warning by adding \"--kinship-details NO\" (however this file will be required to subsequently calculate BLUPs)");}
if(mode!=113&&mode!=116&&mode!=117&&mode!=121&&mode!=122&&mode!=126&&mode!=162)
{printf("you can ignore this warning by adding \"--kinship-details NO\"");}
printf("\n\n");exit(1);
}

if(countcols(filename)==7)
{
printf("Error, %s has seven columns (not nine), indicating it was created using an old version of LDAK; you should remake the corresponding kinship matrix using the latest version of LDAK", filename);
if(mode==121||mode==126)
{printf(", or ignore this warning by adding \"--kinship-details NO\" (however this file will be required to subsequently calculate BLUPs)");}
if(mode!=113&&mode!=116&&mode!=117&&mode!=121&&mode!=122&&mode!=126&&mode!=162)
{printf(", or ignore this warning by adding \"--kinship-details NO\"");}
printf("\n\n");exit(1);
}

sprintf(filename,"%s.grm.adjust",kinstems[k]);
if(just_check(filename)!=0)
{
printf("Error reading %s; this file would have been created when calculating kinships; ", filename);
if(mode==121||mode==126)
{printf("you can ignore this warning by adding \"--kinship-details NO\" (however this file will be required to subsequently calculate BLUPs)");}
if(mode!=113&&mode!=116&&mode!=117&&mode!=121&&mode!=122&&mode!=126&&mode!=162)
{printf("you can ignore this warning by adding \"--kinship-details NO\"");}
printf("\n\n");exit(1);
}

//read top line
if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}
if(fscanf(input, "%s %s ", readstring, readstring2)!=2)
{printf("Error reading Row 1 of %s\n\n", filename);exit(1);}
if(strcmp(readstring,"Datafile")!=0&&strcmp(readstring,"Multiple")!=0)
{printf("Error, %s should begin \"Datafile\" or \"Multiple\" (not %s), suggesting it was created using an old version of LDAK\n\n", filename, readstring);exit(1);}

if((mode==122||mode==162)&&checkroot==1)	//check kinship matrix corresponds to the datafile
{
if(strcmp(readstring,"Multiple")==0)
{printf("Error, the kinship matrix %s uses predictors from multiple datasets; you must ensure all the predictors are contained within the data file %s, then add \"--check-root NO\"\n\n", kinstems[k], datafile);exit(1);}
(void)append_check(readstring3,readstring2,workdir);
if(strcmp(datafile,readstring2)!=0&&strcmp(datafile,readstring3)!=0)
{printf("Error, the kinship matrix %s corresponds to the data file %s, which appears to be different to that provided now (%s); if you are sure the data file is correct, add \"--check-root NO\"\n\n", kinstems[k], readstring2, datafile);exit(1);}
}

//skip power line, read kinnums, skip weights line, read kinsums
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input,"%c", &readchar);}
if(fscanf(input, "Num_Preds %d ", kinnums+k)!=1)
{printf("Error reading Row 3 of %s\n\n", filename);exit(1);}
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input,"%c", &readchar);}
if(fscanf(input, "Denominator %lf ", kinsums+k)!=1)
{printf("Error reading Row 5 of %s\n\n", filename);exit(1);}
fclose(input);

if(kinnums[k]<num_samples_use){flag++;}
}	//end of kindetails=1
else
{kinnums[k]=-9999;kinsums[k]=-9999;}

if(diagonal==1)
{
sprintf(filename,"%s.grm.diag",kinstems[k]);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created when calculating kinships using \"--gxemm-free\"\n\n", filename);exit(1);}
}
}	//end of k loop

if((mode==121||mode==123||mode==124)&&flag>0)
{
if(num_kins==1){printf("Warning, the kinship matrix was computed from fewer predictors (%d) than the number of samples (%d); it might be more efficient to replace this with a region\n\n", kinnums[0], num_samples_use);}
else{printf("Warning, %d of the %d kinship matrices were computed from fewer predictors than the number of samples (%d); it might be more efficient to replace these with regions\n\n", flag, num_kins, num_samples_use);}
}
}	//end of num_kins>0

///////////////////////////

//regions, responses, summaries and fixed

if(num_regs>0)	//check regions exist and contain elements, and set rprune
{
for(r=0;r<num_regs;r++)
{
sprintf(filename,"%s%d",regpref, r+1);
if(just_check(filename)!=0)
{printf("Error reading %s; check file exists and is not empty\n\n", filename);exit(1);}
}

if(rprune==-9999)
{
if(strcmp(weightsfile,"blank")!=0)	//with weights, no pruning should be required
{rprune=2;}
else	//without weights, pruning depends on whether using summaries or not
{
if(strcmp(sumsfile,"blank")!=0){rprune=0.5;}
else{rprune=0.98;}
}
}
}

////////

if(strcmp(respfile,"blank")!=0)	//have phenotypes
{
count=countcols(respfile);
if(count<3)
{printf("Error, %s should have at least three columns (not %d)\n\n", respfile, count);exit(1);}
num_resps=count-2;

if(mpheno>num_resps)
{printf("Error, %s contains only %d phenotypes (so you can not use \"--mpheno %d\")\n\n", respfile, num_resps, mpheno);exit(1);}

if(mpheno2>num_resps)
{printf("Error, %s contains only %d phenotypes (so you can not use \"--mpheno2 %d\")\n\n", respfile, num_resps, mpheno2);exit(1);}

if((mode==229||mode==230)&&(mpheno==-9999||mpheno2==-9999)){printf("Error, you must specify two phenotypes using \"--mpheno\" and \"--mpheno2\"\n\n");exit(1);}

/*
if(strcmp(npheno,"blank")!=0)	//set mpheno
{
if(check_head_ids(respfile,0)==0)
{printf("Error, you can not use \"--pheno-name\" because %s does not have a header row\n\n", respfile);exit(1);}

if((input=fopen(respfile,"r"))==NULL)
{printf("Error opening %s\n\n", respfile);exit(1);}
if(fscanf(input, "%s %s ", readstring, readstring2)!=2)
{printf("Error reading first two elements of %s\n\n", respfile);exit(1);}

mpheno=-9999;
for(j=0;j<num_resps;j++)
{
if(fscanf(input,"%s ", readstring)!=1)
{printf("Error reading element %d of %s\n\n", 3+j, respfile);exit(1);}
if(strcmp(readstring,npheno)==0)
{
if(mpheno!=-9999){printf("Error, %s contains at least two phenotypes called \"%s\"\n\n", respfile, npheno);exit(1);}
mpheno=j+1;
}
}
if(mpheno==-9999){printf("Error, %s does not contain a phenotype called \"%s\"\n\n", respfile, npheno);exit(1);}
fclose(input);
}
*/

if(mpheno==-9999)
{
if(num_resps>1)
{
printf("Error, %s contains multiple phenotypes, so you must specify one using \"--mpheno\"", respfile);
if(mode==121||mode==123||mode==124||mode==126||mode==129||mode==130||mode==151||mode==152||mode==153||mode==154||mode==163)
{printf(", or use \"--mpheno ALL\" to analyse all phenotypes");}
printf("\n\n");exit(1);
}
else{mpheno=1;}
}

keepresps=malloc(sizeof(int)*num_resps);
if(mpheno==-1)
{
num_resps_use=num_resps;
for(m=0;m<num_resps;m++){keepresps[m]=m;}
}
else{num_resps_use=1;keepresps[0]=mpheno-1;}

if(mode==151||mode==152||mode==153||mode==154)	//check whether the phenotypes are binary, and that dichot is consistent
{
head=check_head_ids(respfile,0);
count=countrows(respfile)-head;
readdoubles=malloc(sizeof(double)*count);

//see how many binary
count4=0;
for(m=0;m<num_resps_use;m++)
{
read_values(respfile, readdoubles, count, NULL, 3+keepresps[m], head, 1);

//see if all values are 0/1/missing
count2=0;for(i=0;i<count;i++){count2+=(readdoubles[i]!=-9999&&readdoubles[i]!=0&&readdoubles[i]!=1&&readdoubles[i]!=missingvalue);}

//or if all values are 1/2/missing
count3=0;for(i=0;i<count;i++){count3+=(readdoubles[i]!=-9999&&readdoubles[i]!=1&&readdoubles[i]!=2&&readdoubles[i]!=missingvalue);}

count4+=(count2==0||count3==0);
}

if(num_resps_use==1)
{
if(dichot==-9999&&count4==1)
{printf("Error, the response in %s appears to be binary, so you must use either \"--binary YES\" or \"--binary NO\" to tell LDAK whether or not to use a quasi-logistic model\n\n", respfile);exit(1);}
if(dichot==1&&count4==0)
{printf("Error, the response in %s does not appear to be binary, so you can not use \"--binary YES\"\n\n", respfile);exit(1);}
}
else
{
if(dichot==-9999&&count4>0)
{printf("Error, some of the responses in %s appear to be binary, so you must use either \"--binary YES\" or \"--binary NO\" to tell LDAK whether or not to use a quasi-logistic model\n\n", respfile);exit(1);}
if(dichot==1&&count4<num_resps_use)
{printf("Error, some of the responses in %s appear not to be binary, so you can not use \"--binary YES\"\n\n", respfile);exit(1);}
}

free(readdoubles);
}
}	//end of got respfile
else	//do not have phenotype
{
if(mode==115||mode==121||mode==122||mode==125||mode==138||mode==140||mode==172)	//will fake a response
{num_resps_use=1;keepresps=malloc(sizeof(int));keepresps[0]=0;}
else
{num_resps_use=0;}
}

////////

if(strcmp(sumsfile,"blank")!=0)	//check correct columns and set amb and scaling
{
count=countcols(sumsfile);

if(find_head("Predictor", sumsfile, count)==-1&&find_head("SNP", sumsfile, count)==-1)
{printf("Error, %s should have a column named \"Predictor\" or \"SNP\"\n\n", sumsfile);exit(1);}
if(find_head("A1", sumsfile, count)==-1)
{printf("Error, %s should have a column named \"A1\"\n\n", sumsfile);exit(1);}
if(find_head("A2", sumsfile, count)==-1)
{printf("Error, %s should have a column named \"A1\"\n\n", sumsfile);exit(1);}

if(mode==146||(mode==147&&plet==1))	//only need unsigned statistics (Z, Stat or P)
{
if(find_head("Z", sumsfile, count)==-1&&find_head("Stat", sumsfile, count)==-1&&find_head("P", sumsfile, count)==-1)
{printf("Error, %s should have a column named \"Z\", \"Stat\" or \"P\"\n\n", sumsfile);exit(1);}
}
else	//will need signed statistics (either Z, Stat+Direction or P+Direction)
{
if(find_head("Z", sumsfile, count)==-1&&(find_head("Stat", sumsfile, count)==-1||find_head("Direction", sumsfile, count)==-1)&&(find_head("P", sumsfile, count)==-1||find_head("Direction", sumsfile, count)==-1))
{printf("Error, %s should have a column named \"Z\", two columns named \"Stat\" and \"Direction\", or two columns named \"P\" and \"Direction\"\n\n", sumsfile);exit(1);}
}

if(fixn==-9999)
{
if(find_head("n", sumsfile, count)==-1&&find_head("N", sumsfile, count)==-1)
{printf("Error, %s should have a column named \"n\" or \"N\" (or use \"--fixed-n\" if sample size is constant)\n\n", sumsfile);exit(1);}
}
else
{
if(find_head("n", sumsfile, count)==-1&&find_head("N", sumsfile, count)==-1)
{printf("Warning, you have used \"--fixed-n\", so the sample sizes in %s will be ignored\n\n", sumsfile);}
}

if(amb==-9999)
{
if(mode==146||(mode==147&&plet==1)||mode==150){amb=1;}
else{amb=0;}
}
if(scaling==-9999){scaling=1;}
}

////////

if(strcmp(sums2file,"blank")!=0)	//check correct columns and set scaling (amb already set)
{
count=countcols(sums2file);

if(find_head("Predictor", sums2file, count)==-1&&find_head("SNP", sums2file, count)==-1)
{printf("Error, %s should have a column named \"Predictor\" or \"SNP\"\n\n", sums2file);exit(1);}
if(find_head("A1", sums2file, count)==-1)
{printf("Error, %s should have a column named \"A1\"\n\n", sums2file);exit(1);}
if(find_head("A2", sums2file, count)==-1)
{printf("Error, %s should have a column named \"A1\"\n\n", sums2file);exit(1);}

if(mode==146||(mode==147&&plet==1))	//only need unsigned statistics (Z, Stat or P)
{
if(find_head("Z", sums2file, count)==-1&&find_head("Stat", sums2file, count)==-1&&find_head("P", sums2file, count)==-1)
{printf("Error, %s should have a column named \"Z\", \"Stat\" or \"P\"\n\n", sums2file);exit(1);}
}
else	//will need signed statistics (either Z, Stat+Direction or P+Direction)
{
if(find_head("Z", sums2file, count)==-1&&(find_head("Stat", sums2file, count)==-1||find_head("Direction", sums2file, count)==-1)&&(find_head("P", sums2file, count)==-1||find_head("Direction", sums2file, count)==-1))
{printf("Error, %s should have a column named \"Z\", columns named \"Stat\" and \"Direction\", or columns named \"P\" and \"Direction\"\n\n", sums2file);exit(1);}
}

if(fixn2==-9999)
{
if(find_head("n", sums2file, count)==-1&&find_head("N", sums2file, count)==-1)
{printf("Error, %s should have a column named \"n\" or \"N\" (or use \"--fixed-n2\" if sample size is constant)\n\n", sums2file);exit(1);}
}
else
{
if(find_head("n", sums2file, count)==-1&&find_head("N", sums2file, count)==-1)
{printf("Warning, you have used \"--fixed-n2\", so the sample sizes in %s will be ignored\n\n", sums2file);}
}

if(scaling2==-9999){scaling2=1;}
}

////////

if(strcmp(covarfile,"blank")!=0)	//have covariates (will add an intercept)
{
count=countcols(covarfile);
if(count<3){printf("Error, %s should have at least three columns (not %d)\n\n", covarfile, count);exit(1);}

num_covars=count-1;
}
else	//do not have covariates
{
if(mode==121||mode==122||mode==123||mode==124||mode==125||mode==126||mode==127||mode==128||mode==129||mode==130||mode==229||mode==230||mode==131||mode==132||mode==133||mode==138||mode==140||mode==151||mode==152||mode==153||mode==154||mode==156||mode==172||mode==173||mode==175||mode==194||mode==229||mode==230)	//will be using covar (these are all cases except adjust-grm, where covariates required)
{num_covars=1;}
else
{num_covars=0;}
}

if(strcmp(envfile,"blank")!=0)	//have environmental variables
{
count=countcols(envfile);
if(count<3){printf("Error, %s should have at least three columns (not %d)\n\n", envfile, count);exit(1);}
num_envs=count-2;

if(discenv==1&&num_envs<2)
{printf("Error, when using \"--subgroups YES\" there must be at least 2 environmental variables\n\n");exit(1);}

if(discenv==1&&(mode==121||mode==123||mode==124))	//estimating variances - check num_kins
{
if(num_kins!=2*num_envs){printf("Error, %s you must provide %d kinship matrices, created using \"--gxemm-free\" with \"--subgroups YES\" (not %d)\n\n", kinlist, 2*num_envs, num_kins);exit(1);}
}

if(mode==131&&num_envs!=1)
{printf("Error, when using \"--linear\" there must be exactly one environmental variable (not %d)\n\n", num_envs);exit(1);}
}
else
{num_envs=0;}

if(strcmp(topfile,"blank")!=0)	//have top preds
{num_tops=countrows(topfile);}
else
{num_tops=0;}

//set total number of fixed effects
num_fixed=num_covars+num_envs+num_tops;

if(strcmp(offsetfile,"blank")!=0)	//have offsets
{
count=countcols(offsetfile);
if(count!=3){printf("Error, %s should have three columns (not %d)\n\n", offsetfile, count);exit(1);}
}

///////////////////////////

//calculating weights

if(strcmp(infosfile,"blank")!=0)	//check correct size
{
if(countcols(infosfile)!=2)
{printf("Error, %s should have two columns, providing predictor names then info score (not %d)\n\n", infosfile, countcols(infosfile));exit(1);}
}

///////////////////////////

//calculating kinships

if(strcmp(partpref,"blank")!=0)	//check partitions exist and contain elements
{
for(q=0;q<num_parts;q++)
{
sprintf(filename,"%s%d",partpref, q+1);
if(just_check(filename)!=0)
{printf("Error reading %s; check file exists and is not empty\n\n", filename);exit(1);}
}
}

if(strcmp(invsfile,"blank")!=0)	//check inverse files exist and correct size
{
sprintf(filename,"%s.inverse.predictors", invsfile);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--solve-gre\"\n\n", filename);exit(1);}
count=countrows(filename);

sprintf(filename,"%s.inverse.bin", invsfile);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--solve-gre\"\n\n", filename);exit(1);}

if((input=fopen(filename,"rb"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}
fseeko(input, 0, SEEK_END);
if(ftello(input)!=(off_t)sizeof(float)*(count+1)*count/2)
{printf("Error reading %s; should have size %jd not %jd\n\n", filename, (off_t)sizeof(float)*(count+1)*count/2, ftello(input));exit(1);}
fclose(input);
}

///////////////////////////

//reml, blup and he/pcgc

if(strcmp(hersfile,"blank")!=0)	//check size and that values add to one
{
if(countcols(hersfile)!=1)
{printf("Error, %s should have one column (not %d)\n\n", hersfile, countcols(hersfile));exit(1);}

if(countrows(hersfile)!=num_kins+num_regs)
{printf("Error, %s should have %d rows (not %d)\n\n", hersfile, num_kins+num_regs, countrows(hersfile));exit(1);}

vstarts=malloc(sizeof(double)*(num_kins+num_regs));
read_values(hersfile, vstarts, num_kins+num_regs, NULL, 1, 0, 0);
sum=0;for(k=0;k<num_kins+num_regs;k++){sum+=vstarts[k];}
if(sum<0||sum>1)
{printf("Warning, the heritabilities in %s sum to %.4f (normally the sum is between zero and one)\n\n", hersfile, sum);}
free(vstarts);
}

if(strcmp(oversfile,"blank")!=0)	//check size and read proportions into ssums
{
count=countcols(oversfile);
if(count!=num_kins){printf("Error, the number of categories in %s (%d) does not match the number of kinship matrices (%d)\n\n", oversfile, count, num_kins);exit(1);}
if(countrows(oversfile)!=count+3){printf("Error, %s should have %d rows (not %d), suggesting the file has been changed since creation with \"--calc-overlaps\"\n\n", oversfile, count+3, countrows(oversfile));exit(1);}

ssums=malloc(sizeof(double*)*num_kins);
for(k=0;k<num_kins;k++)
{
ssums[k]=malloc(sizeof(double)*(num_kins+2));
read_values(oversfile, ssums[k], num_kins+2, NULL, k+1, 1, 0);
}
}

////////

if(strcmp(remlfile,"blank")!=0)	//must be blup or pred - check numbers match those provided (will read covariates in later) 
{
if((input=fopen(remlfile,"r"))==NULL)
{printf("Error opening %s\n\n",remlfile);exit(1);}

printf("Reading details from %s\n", remlfile);

if(fscanf(input, "Num_Kinships %d ", &readint)!=1)
{printf("Error reading Row 1 of %s (should begin \"Num_Kinships\"), suggesting the file has been changed since creation when performing REML\n\n",remlfile);exit(1);}
if(readint>0&&num_kins==0)
{
if(readint==1){printf("Error, you must use \"--grm\" to provide the kinship matrix used when performing REML\n\n");exit(1);}
else{printf("Error, you must use \"--mgrm\" to provide the %d kinship matrices used when performing REML\n\n", readint);exit(1);}
}
if(readint!=num_kins)
{printf("Error, the number of kinship matrices provided now (%d) is different to the number in %s (%d)\n\n", num_kins, remlfile, readint);exit(1);}

if(fscanf(input, "Num_Regions %d ", &num_scores)!=1)
{printf("Error reading Row 2 of %s (should begin \"Num_Regions\"), suggesting the file has been changed since creation when performing REML\n\n",remlfile);exit(1);}

if(num_scores>0&&mode==125){printf("Error, it is not possible to use \"--reml-pred\" if regions were used when performing REML\n\n");exit(1);}

if(fscanf(input, "Num_Top_Predictors %d ", &num_tops)!=1)
{printf("Error reading Row 3 of %s (should begin \"Num_Top_Predictors\"), suggesting the file has been changed since creation when performing REML\n\n",remlfile);exit(1);}

if(num_tops>0&&mode==125){printf("Error, it is not possible to use \"--reml-pred\" if top predictors were when performing REML\n\n");exit(1);}

if(fscanf(input, "Num_Covariates %d ", &readint)!=1)
{printf("Error reading Row 4 of %s (should begin \"Num_Covariates\"), suggesting the file has been changed since creation when performing REML\n\n",remlfile);exit(1);}
if(strcmp(covarfile,"blank")!=0)	//covariates provided
{
if(readint==1){printf("Error, the covariates provided now were not used when performing REML\n\n");exit(1);}
if(num_covars!=readint){printf("Error, the number of covariates in %s (%d) does not match the number in %s (%d)\n\n", covarfile, num_covars, remlfile, readint);exit(1);}
}
else	//not provided
{
if(readint>1){printf("Error, you must use \"--covar\" to provide the covariates used when performing REML\n\n");exit(1);}
}

if(fscanf(input, "Num_Environments %d ", &readint)!=1)
{printf("Error reading Row 5 of %s (should begin \"Num_Environments\"), suggesting the file has been changed since creation when performing REML\n\n",remlfile);exit(1);}
if(strcmp(envfile,"blank")!=0)	//environments provided
{
if(readint==0){printf("Error, the environmental variables provided now were not used when performing REML\n\n");exit(1);}
if(num_envs!=readint){printf("Error, the number of environmental variables in %s (%d) does not match the number in %s (%d)\n\n", envfile, num_envs, remlfile, readint);exit(1);}
}
else	//not provided
{
if(readint>0){printf("Error, you must use \"--enviro\" to provide the environmental variables used when performing REML\n\n");exit(1);}
}

if(fscanf(input, "Blupfile %s ", readstring)!=1)
{printf("Error reading Row 6 of %s (should begin \"Blupfile\"), suggesting the file has been changed since creation when performing REML\n\n",remlfile);exit(1);}
if(strcmp(readstring,"none")!=0)
{
strcpy(blupfile2,readstring);
if(append_check(blupfile,blupfile2,workdir)!=0)
{printf("Error reading blupfile %s; this file would have been created when performing REML\n\n", blupfile);exit(1);}
if(countcols(blupfile)!=2+2*(num_kins+num_scores))
{printf("Error, %s should have %d columns (not %d), suggesting the file has been changed since creation when performing REML\n\n ", blupfile, 2*(num_kins+num_scores), countcols(blupfile));exit(1);}
}
else{strcpy(blupfile,"blank");}

if(fscanf(input, "Regfile %s ", readstring)!=1)
{printf("Error reading Row 7 of %s (should begin \"Regfile\"), suggesting the file has been changed since creation when performing REML\n\n",remlfile);exit(1);}
if(strcmp(readstring,"none")!=0)
{
strcpy(regfile2,readstring);
if(append_check(regfile,regfile2,workdir)!=0)
{printf("Error reading regfile %s; this file would have been created when performing REML\n\n", regfile);exit(1);}
if(countcols(regfile)!=4+num_scores+(num_tops>0))
{printf("Error, %s should have %d columns (not %d), suggesting the file has been changed since creation when performing REML\n\n", regfile, 4+num_scores+(num_tops>0), countcols(regfile));exit(1);}
}

if(fscanf(input, "Coeffsfile %s ", cofile2)!=1)
{printf("Error reading Row 8 of %s (should begin \"Coeffsfile\"), suggesting the file has been changed since creation when performing REML\n\n",remlfile);exit(1);}
if(append_check(cofile,cofile2,workdir)!=0)
{printf("Error reading coeffsfile %s; this file would have been created when performing REML\n\n", cofile);exit(1);}
fclose(input);

printf("\n");
}	//end of using remlfile

////////

if(mode==123||mode==124)	//he or pcgc
{
//first set adjusted
if(adjusted==-9999)
{
if(num_fixed==1){adjusted=0;}
else{adjusted=1;}
}

if(adjusted==1&&num_kins>0)
{
//first check all kinships match
sprintf(filename,"%s.grm.root", kinstems[0]);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--adjust-grm\"; if you are not using adjusted kinships, use \"--adjusted NO\"\n\n", filename);exit(1);}
if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}
if(fscanf(input, "Kinship %s ", readstring)!=1)
{printf("Error reading Row 1 of %s (should begin \"Kinship\"), suggesting the file has been changed since creation with \"--adjust-grm\"\n\n",filename);exit(1);}
if(fscanf(input, "Covariates %s ", readstring2)!=1)
{printf("Error reading Row 2 of %s (should begin \"Covariates\"), suggesting the file has been changed since creation with \"--adjust-grm\"\n\n",filename);exit(1);}
if(fscanf(input, "Top_Predictors %s ", readstring3)!=1)
{printf("Error reading Row 3 of %s (should begin \"Top_Predictors\"), suggesting the file has been changed since creation with \"--adjust-grm\"\n\n",filename);exit(1);}
if(fscanf(input, "Environments %s ", readstring4)!=1)
{printf("Error reading Row 4 of %s (should begin \"Environments\"), suggesting the file has been changed since creation with \"--adjust-grm\"\n\n",filename);exit(1);}
fclose(input);

for(k=1;k<num_kins;k++)
{
sprintf(filename2,"%s.grm.root", kinstems[k]);
if(just_check(filename2)!=0)
{printf("Error reading %s; this file would have been created using \"--adjust-grm\"; if you are not using adjusted kinships, use \"--adjusted NO\"\n\n", filename2);exit(1);}
if((input2=fopen(filename2,"r"))==NULL)
{printf("Error opening %s\n\n",filename2);exit(1);}
if(fscanf(input2, "Kinship %s ", readstring)!=1)
{printf("Error reading Row 1 of %s (should begin \"Kinship\"), suggesting the file has been changed since creation with \"--adjust-grm\"\n\n",filename2);exit(1);}
if(fscanf(input2, "Covariates %s ", readstring)!=1)
{printf("Error reading Row 2 of %s (should begin \"Covariates\"), suggesting the file has been changed since creation with \"--adjust-grm\"\n\n",filename2);exit(1);}
if(strcmp(readstring,readstring2)!=0)
{printf("Error, the covariate file listed in %s (%s) does not match that listed in %s (%s)\n\n", filename2, readstring, filename, readstring2);exit(1);}
if(fscanf(input2, "Top_Predictors %s ", readstring)!=1)
{printf("Error reading Row 3 of %s (should begin \"Top_Predictors\"), suggesting the file has been changed since creation with \"--adjust-grm\"\n\n",filename2);exit(1);}
if(strcmp(readstring,readstring3)!=0)
{printf("Error, the top predictors file listed in %s (%s) does not match that listed in %s (%s)\n\n", filename2, readstring, filename, readstring3);exit(1);}
if(fscanf(input2, "Environments %s ", readstring)!=1)
{printf("Error reading Row 4 of %s (should begin \"Environments\"), suggesting the file has been changed since creation with \"--adjust-grm\"\n\n",filename2);exit(1);}
if(strcmp(readstring,readstring4)!=0)
{printf("Error, the environmental variables file listed in %s (%s) does not match that listed in %s (%s)\n\n", filename2, readstring, filename, readstring4);exit(1);}
fclose(input2);
}

if(strcmp(covarfile,"blank")==0&&strcmp(readstring2,"none")!=0)
{printf("Error, you must use \"--covar\" to provide the covariates used with \"--adjust-grm\"\n\n");exit(1);}
if(strcmp(covarfile,"blank")!=0&&strcmp(readstring2,"none")==0)
{printf("Error, covariates are provided now, but were not used with \"--adjust-grm\"\n\n");exit(1);}

if(strcmp(topfile,"blank")==0&&strcmp(readstring3,"none")!=0)
{printf("Error, you must use \"--top-preds\" to provide the top predictors used with \"--adjust-grm\"\n\n");exit(1);}
if(strcmp(topfile,"blank")!=0&&strcmp(readstring3,"none")==0)
{printf("Error, top predictors are provided now, but were not used with \"--adjust-grm\"\n\n");exit(1);}

if(strcmp(envfile,"blank")==0&&strcmp(readstring4,"none")!=0)
{printf("Error, you must use \"--enviro\" to provide the environmental variables used with \"--adjust-grm\"\n\n");exit(1);}
if(strcmp(envfile,"blank")!=0&&strcmp(readstring4,"none")==0)
{printf("Error, environmental variables are provided now, but were not used with \"--adjust-grm\"\n\n");exit(1);}

if(checkroot==1&&strcmp(readstring2,"none")!=0)	//check covarfiles match
{
(void)append_check(readstring,readstring2,workdir);
if(strcmp(covarfile,readstring)!=0&&strcmp(covarfile,readstring2)!=0)
{printf("Error, the kinship matrix %s corresponds to the covariate file %s, which appears to be different to that provided now (%s); if you are sure the kinship matrix is correct, use \"--check-root NO\"\n\n", kinstems[0], readstring2, covarfile);exit(1);}
}
if(checkroot==1&&strcmp(readstring3,"none")!=0)	//check top predictor files match
{
(void)append_check(readstring,readstring3,workdir);
if(strcmp(topfile,readstring)!=0&&strcmp(topfile,readstring3)!=0)
{printf("Error, the kinship matrix %s corresponds to the top predictors file %s, which appears to be different to that provided now (%s); if you are sure the kinship matrix is correct, use \"--check-root NO\"\n\n", kinstems[0], readstring3, topfile);exit(1);}
}
if(checkroot==1&&strcmp(readstring4,"none")!=0)	//check envfiles match
{
(void)append_check(readstring,readstring4,workdir);
if(strcmp(envfile,readstring)!=0&&strcmp(envfile,readstring4)!=0)
{printf("Error, the kinship matrix %s corresponds to the environmental variables file %s, which appears to be different to that provided now (%s); if you are sure the kinship matrix is correct, use \"--check-root NO\"\n\n", kinstems[0], readstring4, envfile);exit(1);}
}
}
}

////////

if(strcmp(relfile,"blank")!=0)
{
if(countcols(relfile)!=5&&countcols(relfile)!=6){printf("Error, %s should have either five or six columns (not %d), containing pairs of samples, their estimated relatedness, and perhaps also their environmental similarity\n\n", relfile, countcols(relfile));exit(1);}
}

///////////////////////////

//association analysis

if(mode==131&&trios==1)
{
count=countcols(famfile);
if(count<4){printf("Error, %s should have at least four columns (not %d), with the 3rd and 4th columns containing paternal and maternal IDs (PID and MID)\n\n", famfile, count);exit(1);}
}

if(strcmp(sampwfile,"blank")!=0)
{
count=countcols(sampwfile);
if(count!=3)
{printf("Error, %s should have three columns (not %d), containing the two sample IDs then a positive weight\n\n", sampwfile, count);exit(1);}
}

if(strcmp(locofile,"blank")!=0||strcmp(bocofile,"blank")!=0)	//do some basic checks for root, prs and details (more checks in getnums.c)
{
//first root

if(strcmp(locofile,"blank")!=0){sprintf(filename,"%s.root", locofile);}
else{sprintf(filename,"%s.root", bocofile);}
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--kvik-step1\" (or \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\" with \"--LOCO YES\")\n\n", filename);exit(1);}

count=countrows(filename);
if(count!=10)
{printf("Error, %s should have ten rows (not %d), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n", filename, count);exit(1);}

if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}

//read data, phenotypes and analysis model

if(fscanf(input, "%s %s ", readstring, readstring2)!=2)
{printf("Error reading Row 1 of %s\n\n", filename);exit(1);}
if(strcmp(readstring,"Datafile")!=0)
{printf("Error, %s should begin \"Datafile\" (not %s), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n",filename, readstring);exit(1);}

if(fscanf(input, "%s %s ", readstring, readstring3)!=2)
{printf("Error reading Row 2 of %s\n\n", filename);exit(1);}
if(strcmp(readstring,"Phenotypes")!=0)
{printf("Error, reading Row 2 of %s (should begin \"Phenotypes\"), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n",filename);exit(1);}

if(fscanf(input, "%s %s ", readstring, readstring4)!=2)
{printf("Error reading Row 3 of %s\n\n", filename);exit(1);}
if(strcmp(readstring,"Analysis")!=0)
{printf("Error, reading Row 3 of %s (should begin \"Analysis\"), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n",filename);exit(1);}

fastgwa=-1;
if(strcmp(readstring4,"KVIK")==0){fastgwa=0;}
if(strcmp(readstring4,"fastGWA")==0){fastgwa=1;}
if(strcmp(readstring4,"KVIK_Pedigree")==0){fastgwa=2;}
if(strcmp(readstring4,"KVIK_BLANK")==0){fastgwa=3;}

if(fastgwa==-1){printf("Error, Element 2 of Row 3 of %s should be \"KVIK\", \"KVIK_Pedigree\", \"KVIK_BLANK\" or \"fastGWA\" (not %s), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n", filename, readstring4);exit(1);}

//read regression model

if(fscanf(input, "%s %s ", readstring, readstring4)!=2)
{printf("Error reading Row 4 of %s\n\n", filename);exit(1);}
if(strcmp(readstring,"Regression")!=0)
{printf("Error, reading Row 4 of %s (should begin \"Regression\"), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n",filename);exit(1);}

dichot=-1;
if(strcmp(readstring4,"Linear")==0){dichot=0;}
if(strcmp(readstring4,"Logistic")==0){dichot=1;}

if(dichot==-1){printf("Error, Element 2 of Row 4 of %s should be \"Linear\" or \"Logistic\" (not %s), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n", filename, readstring4);exit(1);}

if(mode==131&&dichot==1)
{printf("Error, \"--binary YES\" was used when making the prediction models, so you should now use \"--logistic\" (instead of \"--linear\")\n\n");exit(1);}
if(mode==132&&dichot==0)
{printf("Error, \"--binary YES\" was NOT used when making the prediction models, so you should now use \"--linear\" (instead of \"--logistic\")\n\n");exit(1);}

//skip one row
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}

//read covariates and tops

if(fscanf(input, "%s %s ", readstring, readstring4)!=2)
{printf("Error reading Row 6 of %s\n\n", filename);exit(1);}
if(strcmp(readstring,"Covariates")!=0)
{printf("Error, reading Row 6 of %s (should begin \"Covariates\"), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n",filename);exit(1);}

if(fscanf(input, "%s %s ", readstring, readstring5)!=2)
{printf("Error reading Row 7 of %s\n\n", filename);exit(1);}
if(strcmp(readstring,"Top_Predictors")!=0)
{printf("Error, reading Row 7 of %s (should begin \"Top_Predictors\"), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n",filename);exit(1);}

if(mode==131||mode==132)	//check consistent
{
if(strcmp(covarfile,"blank")==0&&strcmp(readstring4,"none")!=0)
{printf("Error, \"--covar\" was used when making the prediction models, but is not used now\n\n");exit(1);}
if(strcmp(covarfile,"blank")!=0&&strcmp(readstring4,"none")==0)
{printf("Error, \"--covar\" is used now, but was not used when making the prediction models\n\n");exit(1);}

if(strcmp(topfile,"blank")==0&&strcmp(readstring5,"none")!=0)
{printf("Error, \"--top-preds\" was used when making the prediction models, but is not used now\n\n");exit(1);}
if(strcmp(topfile,"blank")!=0&&strcmp(readstring5,"none")==0)
{printf("Error, \"--top-preds\" is used now, but was not used when making the prediction models\n\n");exit(1);}
}

if(checkroot==1)	//check names consistent (most checks only required for linear / logistic)
{
flag=0;
if(strcmp(datafile,readstring2)!=0)
{printf("Error, the data file used when making the prediction models (%s) does not match that used now (%s)\n\n", readstring2, datafile);flag=1;}

if(mode==131||mode==132)
{
if(strcmp(respfile,readstring3)!=0)
{printf("Error, the phenotypes file used when making the prediction models (%s) does not match that used now (%s)\n\n", readstring3, respfile);flag=1;}
if(strcmp(covarfile,"blank")!=0&&strcmp(covarfile,readstring4)!=0)
{printf("Error, the covariates file used when making the prediction models (%s) does not match that used now (%s)\n\n", readstring4, covarfile);flag=1;if(strcmp(topfile,"blank")!=0&&strcmp(topfile,readstring5)!=0)
{printf("Error, the top predictors file used when making the prediction models (%s) does not match that used now (%s)\n\n", readstring5, topfile);flag=1;}
}
}
if(flag==1){printf("If you are sure you have provided the correct files, use \"--check-root NO\"\n\n");exit(1);}
}

//skip two rows
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}

//read num_chr (save in num_chr2)
if(fscanf(input, "%s %d ", readstring, &num_chr2)!=2)
{printf("Error reading Row 10 of %s\n\n", filename);exit(1);}
if(strcmp(readstring,"Num_Chromosomes")!=0)
{printf("Error, reading Row 10 of %s (should begin \"Num_Chromosomes\"), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n", filename);exit(1);}

fclose(input);

//now details - get sample size, effective size, scaling and maybe power

if(strcmp(locofile,"blank")!=0){sprintf(filename,"%s.loco.details", locofile);}
else{sprintf(filename,"%s.loco.details", bocofile);}
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--kvik-step1\" (or \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\" with \"--LOCO YES\")\n\n", filename);exit(1);}

count=countrows(filename);
if(count!=7)
{printf("Error, %s should have seven rows (not %d), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n", filename, count);exit(1);}

//open and skip one row
if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}

//read sample size
if(fscanf(input, "%s %d ", readstring, &readint)!=2)
{printf("Error reading Row 2 of %s\n\n", filename);exit(1);}
if(strcmp(readstring,"Actual_Sample_Size")!=0)
{printf("Error, reading Row 2 of %s (should begin \"Actual_Sample_Size\"), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n", filename);exit(1);}

//read neff
if(fscanf(input, "%s %s ", readstring, readstring2)!=2)
{printf("Error reading Row 3 of %s\n\n", filename);exit(1);}
if(strcmp(readstring,"Approx_Effective_Sample_Size")!=0)
{printf("Error, reading Row 3 of %s (should begin \"Approx_Effective_Sample_Size\"), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n", filename);exit(1);}
neff=atof(readstring2);

//read scaling (also compute scaling^.5)
if(fscanf(input, "%s %s ", readstring, readstring2)!=2)
{printf("Error reading Row 4 of %s\n\n", filename);exit(1);}
if(strcmp(readstring,"Scaling_Estimate")!=0)
{printf("Error, reading Row 4 of %s (should begin \"Scaling_Estimate\"), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n", filename);exit(1);}
scal1=atof(readstring2);
scal2=pow(scal1,.5);

if(fastgwa!=1)	//skip a row and read power
{
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}

if(fscanf(input, "%s %s ", readstring, readstring2)!=2)
{printf("Error reading Row 6 of %s\n\n", filename);exit(1);}
if(strcmp(readstring,"Power")!=0)
{printf("Error, reading Row 6 of %s (should begin \"Power\"), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n", filename);exit(1);}
power=atof(readstring2);
}
fclose(input);

if(strcmp(locofile,"blank")!=0)	//lastly prs
{
sprintf(filename,"%s.loco.prs", locofile);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--kvik-step1\" (or \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\" with \"--LOCO YES\")\n\n", filename);exit(1);}

count=countcols(filename);
if(fastgwa==0&&count!=2+num_chr2){printf("Error, %s should have %d columns (not %d), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n", filename, 2+num_chr2, count);exit(1);}
if(fastgwa!=0&&count!=3){printf("Error, %s should have three columns (not %d), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n", filename, count);exit(1);}

count=countrows(filename);
if(count!=readint+1)
{printf("Error, %s should have %d rows (not %d), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n", filename, readint+1, count);exit(1);}
}
}

////////

if(strcmp(genefile,"blank")!=0)	//check correct size
{
if(countcols(genefile)!=4&&countcols(genefile)!=5)
{printf("Error, %s should either four or five columns (not %d), providing the name, chromosome, start and end basepairs, and perhaps the orientation (\"+\" or \"-\") of each genomic region\n\n", genefile, countcols(genefile));exit(1);}

if(countcols(genefile)!=5&&(up_buffer>0||down_buffer>0))
{printf("Error, to use \"--up-buffer\" and \"--down-buffer\", %s should have five columns (not %d), providing the name, chromosome, start and end basepairs, and the orientation (\"+\" or \"-\") of each genomic region\n\n", genefile, countcols(genefile));exit(1);}
}

///////////////////////////

//sumher

if(strcmp(annpref,"blank")!=0)	//check annotations exist and contain elements
{
for(q=0;q<num_anns;q++)
{
sprintf(filename,"%s%d",annpref, q+1);
if(just_check(filename)!=0)
{printf("Error reading %s; check file exists and is not empty\n\n", filename);exit(1);}
}
}

if(strcmp(labfile,"blank")!=0&&num_parts>0)
{
if(countcols(labfile)!=1)
{printf("Error, %s should have one column (not %d)\n\n", labfile, countcols(labfile));exit(1);}

if(countrows(labfile)!=num_parts)
{printf("Error, the number of rows of %s (%d) does not match the number of partitions (%d)\n\n", labfile, countrows(labfile), num_parts);exit(1);}
}

if(strcmp(labfile,"blank")!=0&&num_anns>0)
{
if(countcols(labfile)!=1)
{printf("Error, %s should have one column (not %d)\n\n", labfile, countcols(labfile));exit(1);}

if(countrows(labfile)!=num_anns)
{printf("Error, the number of rows of %s (%d) does not match the number of annotations (%d)\n\n", labfile, countrows(labfile), num_anns);exit(1);}
}

if(strcmp(taglist,"blank")!=0)	//have taglist (might also have matlist - can not have matlist without taglist)
{
num_tags=countrows(taglist);
if(num_tags==1){printf("Error, %s must contain at least two tagging files (not one)\n\n", taglist);exit(1);}

tagstems=malloc(sizeof(char *)*num_tags);
if((input=fopen(taglist,"r"))==NULL)
{printf("Error opening %s\n\n",taglist);exit(1);}
for(k=0;k<num_tags;k++)
{
tagstems[k]=malloc(sizeof(char)*500);
if(fscanf(input, "%s ", filename)!=1)
{printf("Error reading Row %d of %s\n\n", k+1, taglist);exit(1);}
if(append_check(tagstems[k],filename,workdir)!=0)
{printf("Error reading tagging file %s; check file exists and is not empty\n\n", tagstems[k]);exit(1);}

for(k2=0;k2<k;k2++)	//check different to previous
{
if(strcmp(tagstems[k],tagstems[k2])==0)
{printf("Error reading %s; the file %s appears twice\n\n", taglist, tagstems[k]);exit(1);}
}
}
fclose(input);

if(strcmp(matlist,"blank")!=0)
{
if(countrows(matlist)!=num_tags)
{printf("Error, the number of rows of %s (%d) does not match the number of rows of %s (%d)\n\n", matlist, countrows(matlist), taglist, num_tags);exit(1);}

matstems=malloc(sizeof(char *)*num_tags);
if((input=fopen(matlist,"r"))==NULL)
{printf("Error opening %s\n\n",matlist);exit(1);}
for(k=0;k<num_tags;k++)
{
matstems[k]=malloc(sizeof(char)*500);
if(fscanf(input, "%s ", filename)!=1)
{printf("Error reading Row %d of %s\n\n", k+1, matlist);exit(1);}
if(append_check(matstems[k],filename,workdir)!=0)
{printf("Error reading heritability matrix %s; check file exists and is not empty\n\n", matstems[k]);exit(1);}

for(k2=0;k2<k;k2++)	//check different to previous
{
if(strcmp(matstems[k],matstems[k2])==0)
{printf("Error reading %s; the file %s appears twice\n\n", matlist, matstems[k]);exit(1);}
}
}
fclose(input);
}
}

////////

if(strcmp(tagfile,"blank")!=0)	//check correct size and set num_parts and parttype (0=annotations, 1=partitions)
{
count=countcols(tagfile);
if(count<10)
{printf("Error, %s should have at least 10 columns (not %d), suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", tagfile, count);exit(1);}
num_parts=count-9;

count=countrows_min(tagfile,num_parts+3);
if(count<num_parts+3)
{printf("Error, %s should have at least %d rows, suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", tagfile, num_parts+3);exit(1);}

if((input=fopen(tagfile,"r"))==NULL)
{printf("Error opening %s\n\n",tagfile);exit(1);}
if(fscanf(input, "%s %s %s %s %s %s %s %s %s ", readstring, readstring2, readstring2, readstring2, readstring2, readstring2, readstring2, readstring2, readstring2)!=9)
{printf("Error reading Row 1 of %s\n\n", tagfile);exit(1);}
if(strcmp(readstring,"Predictor")!=0)
{printf("Error reading %s; first element should be \"Predictor\" (not %s), suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", tagfile, readstring);exit(1);}

for(q=0;q<num_parts;q++)
{
if(fscanf(input, "%s ", readstring)!=1)
{printf("Error reading element %d of %s, suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", j+1, tagfile);exit(1);}
}
parttype=(strcmp(readstring,"Base")!=0);

fclose(input);
}

////////

if(strcmp(altfile,"blank")!=0)	//check correct size
{
count=countcols(altfile);
if(count!=2)
{printf("Error, %s should have two columns, providing predictor names then taggings (not %d)\n\n", altfile, count);exit(1);}
}

//dont check cvsfile

if(strcmp(catfile,"blank")!=0)	//set num_reds, keeppart and keeppart2
{
keepparts=malloc(sizeof(int)*num_parts);
for(j=0;j<num_parts;j++){keepparts[j]=0;}
if(parttype==0){keepparts[num_parts-1]=1;}

count=countrows(catfile);
if((input=fopen(catfile,"r"))==NULL)
{printf("Error opening %s\n\n",catfile);exit(1);}

for(j=0;j<count;j++)
{
if(fscanf(input, "%d ", &readint)!=1)
{printf("Error reading Element %d of %s\n\n", j+1, catfile);exit(1);}

if(readint<0||(readint==0&&(parttype==1||count>1)))
{printf("Error reading Element %d of %s; all values must be positive integers (not %d)", j+1, catfile, readint);exit(1);}

if(readint==0){break;}	//shorthand for only keep/use base

if(parttype==0&&readint>num_parts-1){printf("Error reading %s; Row %d specifies Annotation %d, which is greater than the number of annotations (%d)\n\n", catfile, j+1, readint, num_parts-1);exit(1);}
if(parttype==1&&readint>num_parts){printf("Error reading %s; Row %d specifies Partition %d, which is greater than the number of partitions (%d)\n\n", catfile, j+1, readint, num_parts);exit(1);}

if(keepparts[readint-1]==1){printf("Error, %d appears twice in %s\n\n", readint, catfile);exit(1);}
keepparts[readint-1]=1;
}
fclose(input);

keepparts2=malloc(sizeof(int)*num_parts);
num_reds=0;
for(j=0;j<num_parts;j++)
{
if(keepparts[j]==1){keepparts2[num_reds]=j;num_reds++;}
}
}

//check taufile in required.c

////////

if(strcmp(powfile,"blank")!=0)	//check contains numbers
{
count=countrows(powfile);
printf("Checking the %d values in %s\n", count, powfile);

if((input=fopen(powfile,"r"))==NULL)
{printf("Error opening %s\n\n",powfile);exit(1);}
for(k=0;k<count;k++)
{
if(fscanf(input, "%[0-9eE.+-] ", readstring)==1)
{
if(atof(readstring)<-1.25||atof(readstring)>0.25)
{printf("Warning, the power on Row %d is %.4f (usually values are between -1.25 and 0.25)\n", k+1, atof(readstring));}
}
else
{
if(fscanf(input, "%s ", readstring)!=1)
{printf("Error reading Row %d of %s\n\n", k+1, powfile);exit(1);}
printf("Error reading %s; the first value on Row %d is unrecognisable (%s)\n\n", powfile, k+1, readstring);exit(1);
}
}
fclose(input);

printf("\n");
}

////////

if(strcmp(matfile,"blank")!=0)	//check correct size (will have tagfile, and have set num_parts)
{
count=countcols(matfile);
if(count<2)
{printf("Error, %s should have at least two columns (not %d), suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", matfile, count);exit(1);}
if(count-1!=num_parts){printf("Error, the number of categories in %s (%d) does not match the number in %s (%d)\n\n", matfile, count-1, tagfile, num_parts);exit(1);}
}

///////////////////////////

//individual-level data prediction, then megaprs

if(strcmp(indhers,"blank")!=0)	//check correct size
{
if(countcols(indhers)!=2)
{printf("Error, %s should have two columns, providing predictor names then info score (not %d)\n\n", indhers, countcols(indhers));exit(1);}
}

////////

//dont check bvsfile

////////

if(strcmp(corslist,"blank")!=0)
{
num_cors=countrows(corslist);
if(mode==157&&num_cors==1){printf("Error, %s must contain at least two correlation stems (not one)\n\n", corslist);exit(1);}

corstems=malloc(sizeof(char *)*num_cors);
if((input=fopen(corslist,"r"))==NULL)
{printf("Error opening %s\n\n",corslist);exit(1);}
for(k=0;k<num_cors;k++)
{
corstems[k]=malloc(sizeof(char)*500);
if(fscanf(input, "%s ", filename)!=1)
{printf("Error reading Row %d of %s\n\n", k+1, corslist);exit(1);}
(void)append_check(corstems[k],filename,workdir);

if(mode==157)	//check different to previous
{
for(k2=0;k2<k;k2++)
{
if(strcmp(corstems[k],corstems[k2])==0)
{printf("Error reading %s; the file %s appears twice\n\n", corslist, corstems[k]);exit(1);}
}
}
}
fclose(input);
}

////////

if(strcmp(corname,"blank")!=0)	//check files
{
//first read root - only need num_preds_used and num_pairs (for checking sizes of bin and bim files)
sprintf(filename,"%s.cors.root", corname);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--calc-cors\"\n\n", filename);exit(1);}
if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}

if(fscanf(input, "%s %s ", readstring, readstring2)!=2)
{printf("Error reading Row 1 of %s\n\n", filename);exit(1);}
if(strcmp(readstring,"Datafile")!=0&&strcmp(readstring,"Multiple")!=0)
{printf("Error, %s should begin \"Datafile\" or \"Multiple\" (not %s), suggesting the file has been changed since creation with \"--calc-cors\"\n\n",filename, readstring);exit(1);}

if(fscanf(input, "Num_Samples %s ", readstring)!=1)
{printf("Error reading Row 2 of %s (should begin \"Num_Samples\"), suggesting the file has been changed since creation with \"--calc-cors\"\n\n",filename);exit(1);}
if(fscanf(input, "Num_Predictors %s ", readstring)!=1)
{printf("Error reading Row 3 of %s (should begin \"Num_Predictors\"), suggesting the file has been changed since creation with \"--calc-cors\"\n\n",filename);exit(1);}
if(fscanf(input, "Num_Samples_Used %s ", readstring)!=1)
{printf("Error reading Row 4 of %s (should begin \"Num_Samples_Used\"), suggesting the file has been changed since creation with \"--calc-cors\"\n\n",filename);exit(1);}

if(fscanf(input, "Num_Predictors_Used %d ", &readint)!=1)
{printf("Error reading Row 5 of %s (should begin \"Num_Predictors_Used\"), suggesting the file has been changed since creation with \"--calc-cors\"\n\n",filename);exit(1);}
if(fscanf(input, "Num_Pairs %jd ", &scount)!=1)
{printf("Error reading Row 6 of %s (should begin \"Num_Pairs\"), suggesting the file has been changed since creation with \"--calc-cors\"\n\n",filename);exit(1);}
fclose(input);

//now check bin file exists and correct size
sprintf(filename,"%s.cors.bin", corname);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--calc-cors\"\n\n", filename);exit(1);}
if((input=fopen(filename,"rb"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}
fseeko(input, 0, SEEK_END);
if(ftello(input)!=(off_t)(sizeof(int)+sizeof(double)*3)*readint+(sizeof(int)+sizeof(float))*scount)
{printf("Error reading %s; should have size %jd, but instead has size %jd\n\n", filename, (off_t)(sizeof(int)+sizeof(double)*3)*readint+(sizeof(int)+sizeof(float))*scount, ftello(input));exit(1);}
fclose(input);

//check bim file exists and correct size
sprintf(filename,"%s.cors.bim", corname);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--calc-cors\"\n\n", filename);exit(1);}

count=countrows(filename);
if(count!=readint)
{printf("Error reading %s; should have %d rows (not %d), suggesting the file has been changed since creation with \"--calc-cors\"\n\n",filename, readint, count);exit(1);}

//check highLD file exists
sprintf(filename,"%s.cors.highld", corname);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--calc-cors\"\n\n", filename);exit(1);}

//check noise file exists and correct size
sprintf(filename,"%s.cors.noise", corname);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--calc-cors\"\n\n", filename);exit(1);}

count=countrows(filename);
if(count!=readint)
{printf("Error reading %s; should have %d rows (not %d), suggesting the file has been changed since creation with \"--calc-cors\"\n\n",filename, readint, count);exit(1);}

}

////////

if(strcmp(pseudostem,"blank")!=0)	//check training and test summaries exist and look ok
{
sprintf(filename,"%s.train.summaries", pseudostem);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--pseudo-summaries\"\n\n", filename);exit(1);}
if(countcols(filename)!=6){printf("Error, %s should have six columns (not %d), suggesting the file has been changed since creation with \"--pseudo-summaries\"\n\n", filename, countcols(filename));exit(1);}

if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}
if(fscanf(input, "Predictor %s ", readstring)!=1)
{printf("Error reading Row 1 of %s (should begin \"Predictor\"), suggesting the file has been changed since creation with \"--pseudo-summaries\"\n\n",filename);exit(1);}
fclose(input);

sprintf(filename,"%s.test.summaries", pseudostem);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--pseudo-summaries\"\n\n", filename);exit(1);}
if(countcols(filename)!=6){printf("Error, %s should have six columns (not %d), suggesting the file has been changed since creation with \"--pseudo-summaries\"\n\n", filename, countcols(filename));exit(1);}

if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}
if(fscanf(input, "Predictor %s ", readstring)!=1)
{printf("Error reading Row 1 of %s (should begin \"Predictor\"), suggesting the file has been changed since creation with \"--pseudo-summaries\"\n\n",filename);exit(1);}
fclose(input);
}

////////

if(strcmp(bestfile,"blank")!=0)	//check has correct size
{
if(countrows(bestfile)!=2){printf("Error, %s should have two rows (not %d), suggesting the file has been changed since creation with \"--mega-prs\"\n\n", bestfile, countrows(bestfile));exit(1);}
if(countcols(bestfile)!=11){printf("Error, %s should have 11 columns (not %d), suggesting the file has been changed since creation with \"--mega-prs\"\n\n", bestfile, countcols(bestfile));exit(1);}

if((input=fopen(bestfile,"r"))==NULL)
{printf("Error opening %s\n\n",bestfile);exit(1);}
if(fscanf(input, "Model %s ", readstring)!=1)
{printf("Error reading Row 1 of %s (should begin \"Model\"), suggesting the file has been changed since creation with \"--mega-prs\"\n\n", bestfile);exit(1);}
fclose(input);
}

////////

if(strcmp(fracfile,"blank")!=0)	//either bolt, bayesr, elastic or mega with ptype 2, 3, 4, 5 or 6
{
//if(mode==151)	//ridge - not allowed fracfile
if(mode==152)	//bolt
{
if(countcols(fracfile)!=2)
{printf("Error, %s should have two columns (not %d), providing p and f2\n\n", fracfile, countcols(fracfile));exit(1);}
}
if(mode==153)	//bayesr
{
if(pointmass==1)
{
if(countcols(fracfile)!=4)
{printf("Error, %s should have four columns (not %d), providing  the fractions for the point mass, then the small, medium and large Gaussian distributions\n\n", fracfile, countcols(fracfile));exit(1);}
}
else
{
if(countcols(fracfile)!=4)
{printf("Error, %s should have four columns (not %d), providing the fractions for the tiny, small, medium and large Gaussian distributions\n\n", fracfile, countcols(fracfile));exit(1);}
}
}
if(mode==154)	//elastic
{
if(countcols(fracfile)!=2)
{printf("Error, %s should have two columns (not %d), providing p and f2\n\n", fracfile, countcols(fracfile));exit(1);}
}

if(ptype==2||ptype==3)	//lasso and ridge
{
if(countcols(fracfile)!=1)
{printf("Error, %s should have one column (not %d), providing the heritability\n\n", fracfile, countcols(fracfile));exit(1);}
}
if(ptype==4)	//bolt
{
if(countcols(fracfile)!=3)
{printf("Error, %s should have three columns (not %d), providing the heritability, p then f2\n\n", fracfile, countcols(fracfile));exit(1);}
}
if(ptype==5)	//bayes
{
if(countcols(fracfile)!=5)
{printf("Error, %s should have five columns (not %d), providing the heritability, then fractions for the point mass and the small, medium and large Gaussian distributions\n\n", fracfile, countcols(fracfile));exit(1);}
}
if(ptype==6)	//bayes-shrink
{
if(countcols(fracfile)!=5)
{printf("Error, %s should have five columns (not %d), providing the heritability, then fractions for the tiny, small, medium and large Gaussian distributions\n\n", fracfile, countcols(fracfile));exit(1);}
}
}

///////////////////////////

//pca, decompose and adjust-grm

if(strcmp(pcastem,"blank")!=0)
{
sprintf(filename,"%s.vect", pcastem);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created when using \"--pca\"\n\n", filename);exit(1);}
count=countcols(filename);
if(count<3)
{printf("Error, %s should have at least three columns (not %d), suggesting the file has been changed since creation with \"--pca\"\n\n",filename, count);exit(1);}
axes=count-2;

sprintf(filename,"%s.values", pcastem);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--pca\"\n\n", filename);exit(1);}

sprintf(filename,"%s.root", pcastem);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--pca\"\n\n", filename);exit(1);}
if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}
if(fscanf(input, "Kinship %s ", readstring)!=1)
{printf("Error reading Row 1 of %s (should begin \"Kinship\"), suggesting the file has been changed since creation with \"--pca\"\n\n",filename);exit(1);}
if(checkroot==1)	//check matches kinship
{
(void)append_check(readstring2,readstring,workdir);
if(strcmp(kinstems[0],readstring)!=0&&strcmp(kinstems[0],readstring2)!=0)
{printf("Error, the pca %s corresponds to the kinship matrix %s, which appears to be different to that provided now (%s); if you are sure the kinship matrix is correct, use \"--check-root NO\"\n\n", pcastem, readstring, kinstems[0]);exit(1);}
}
fclose(input);
}

////////

if(strcmp(eigenfile,"blank")!=0)
{
sprintf(filename,"%s.eigen.id", eigenfile);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--decompose\"\n\n", filename);exit(1);}
sprintf(filename,"%s.eigen.bin", eigenfile);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--decompose\"\n\n", filename);exit(1);}

sprintf(filename,"%s.eigen.root", eigenfile);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--decompose\"\n\n", filename);exit(1);}
if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}
if(fscanf(input, "Kinship %s ", readstring)!=1)
{printf("Error reading Row 1 of %s (should begin \"Kinship\"), suggesting the file has been changed since creation with \"--decompose\"\n\n",filename);exit(1);}
if(checkroot==1)	//check matches kinship
{
(void)append_check(readstring2,readstring,workdir);
if(strcmp(kinstems[0],readstring)!=0&&strcmp(kinstems[0],readstring2)!=0)
{printf("Error, the eigen-decomposition %s corresponds to the kinship matrix %s, which appears to be different to that provided now (%s); if you are sure the eigen-decomposition is correct, use \"--check-root NO\"\n\n", eigenfile, readstring, kinstems[0]);exit(1);}
}
fclose(input);
}

///////////////////////////

//stats, scores, making phenotypes, snps, jackknifing, folds, find gaussian

if(strcmp(scorefile,"blank")!=0)	//check correct size
{
count=countcols(scorefile);
if(count<5)
{printf("Error, %s should have least five columns: \"Predictor\", \"A1\", \"A2\" and \"Centre\" (the average allele count; set to NA if unknown), followed by a column of effect sizes for each profile (not %d)\n\n", scorefile, count);exit(1);}

if((input=fopen(scorefile,"r"))==NULL)
{printf("Error opening %s\n\n",scorefile);exit(1);}
if(fscanf(input, "%s %s %s %s", readstring, readstring2, readstring3, readstring4)!=4)
{printf("Error reading first four elements of %s\n\n", scorefile);exit(1);}

if(strcmp(readstring,"Predictor")!=0)
{printf("Error, the first column of %s should be named \"Predictor\" (not %s)", scorefile, readstring);exit(1);}
if(strcmp(readstring2,"A1")!=0)
{printf("Error, the second column of %s should be named \"A1\" (not %s)", scorefile, readstring2);exit(1);}
if(strcmp(readstring3,"A2")!=0)
{printf("Error, the third column of %s should be named \"A2\" (not %s)", scorefile, readstring3);exit(1);}
if(strcmp(readstring4,"Centre")!=0)
{printf("Error, the fourth column of %s should be named \"Centre\" (not %s)", scorefile, readstring4);exit(1);}
fclose(input);

num_scores=count-4;
}

if(strcmp(cofile,"blank")!=0)	//mode must be 122, 125 or 172
{
count=countrows(cofile);
count2=countcols(cofile);
count3=count-(count2==4);

if(mode==122||mode==125)	//blup or pred - have already checked number of fixed effects 
{
if(count!=1+num_fixed&&count2!=4)
{printf("Error, %s should have %d rows and 4 columns (not %d and %d), suggesting the file has been changed since creation with \"--reml\"\n\n", cofile, 1+num_fixed, count, count2);exit(1);}
}
else	//so calc-scores
{
if(count2!=1&&count2!=4)
{printf("Error, %s should have either a header row and four columns, with Column 2 providing the coefficients (the format produced by \"--reml\"), or one column with no header (not %d)\n\n", cofile, count2);exit(1);}

if(count3>1&&strcmp(covarfile,"blank")==0&&strcmp(envfile,"blank")==0)
{printf("Error, %s contains %d coefficients, so you must use \"--covars\" and/or \"--enviro\" to provide covariates and/or environmental variables\n\n", cofile, count3);exit(1);}
if(count3!=num_fixed)
{printf("Error, the number of coefficients in %s (%d) does not match the total number of covariates and environmental variables (%d)\n\n", cofile, count3, num_fixed);exit(1);}
}
}

if(strcmp(finalfile,"blank")!=0)	//check correct number of scores
{
count=countcols(finalfile);
if(count<5)
{printf("Error, %s should have least five columns: \"Predictor\", \"A1\", \"A2\" and \"Centre\" (the average allele count; set to NA if unknown), followed by a column of effect sizes for each profile (not %d)\n\n", finalfile, count);exit(1);}
if(count!=num_scores+4){printf("Error, %s should contain %d scores (not %d), to match the number in %s\n", finalfile, num_scores, count-4, scorefile);exit(1);} 
}

///////////////////////////

//old code for checking megaprs jackknifes
/*
if(strcmp(jackstem,"blank")!=0)	//check files and set num_blocks
{
sprintf(filename,"%s.jacks.root", jackstem);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--pseudo-jackknife\"\n\n", filename);exit(1);}
if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}
if(fscanf(input, "Summaries %s ", readstring)!=1)
{printf("Error reading Row 1 of %s (should begin \"Summaries\"), suggesting the file has been changed since creation with \"--pseudo-jackknife\"\n\n",filename);exit(1);}
if(checkroot==1)	//check matches summaries
{
(void)append_check(readstring2,readstring,workdir);
if(strcmp(sumsfile,readstring)!=0&&strcmp(sumsfile,readstring2)!=0)
{printf("Error, the jacks file %s corresponds to the summary statistics %s, which appears to be different to that provided now (%s); if you are sure the summary statistics are correct, use \"--check-root NO\"\n\n", jackstem, readstring, sumsfile);exit(1);}
}
if(fscanf(input, "Training_Proportion %s ", readstring)!=1)
{printf("Error reading Row 2 of %s (should begin \"Training_Proportion\"), suggesting the file has been changed since creation with \"--pseudo-jackknife\"\n\n",filename);exit(1);}
if(fscanf(input, "Num_Blocks %d ", &num_blocks)!=1)
{printf("Error reading Row 3 of %s (should begin \"Num_Blocks\"), suggesting the file has been changed since creation with \"--pseudo-jackknife\"\n\n",filename);exit(1);}
if(fscanf(input, "Num_Predictors %d ", &count)!=1)
{printf("Error reading Row 4 of %s (should begin \"Num_Predictors\"), suggesting the file has been changed since creation with \"--pseudo-jackknife\"\n\n",filename);exit(1);}
if(fscanf(input, "Num_Predictors_Used %d ", &count)!=1)
{printf("Error reading Row 5 of %s (should begin \"Num_Predictors_Used\"), suggesting the file has been changed since creation with \"--pseudo-jackknife\"\n\n",filename);exit(1);}

sprintf(filename,"%s.jacks.predictors", jackstem);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created when using \"--pseudo-jackknife\"\n\n", filename);exit(1);}
if(countcols(filename)!=2)
{printf("Error, %s should have two columns (not %d), suggesting the file has been changed since creation with \"--pseudo-jackknife\"\n\n",filename, countcols(filename));exit(1);exit(1);}
if(countrows(filename)!=count)
{printf("Error, %s should have %d rows (not %d), suggesting the file has been changed since creation with \"--pseudo-jackknife\"\n\n",filename, count, countrows(filename));exit(1);}
fclose(input);

sprintf(filename,"%s.jacks.bin", jackstem);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--pseudo-jackknife\"\n\n", filename);exit(1);}

if((input=fopen(filename,"rb"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}
fseeko(input, 0, SEEK_END);
if(ftello(input)!=(off_t)sizeof(double)*count*num_blocks)
{printf("Error reading %s; should have size %jd not %jd\n\n", filename, (off_t)sizeof(double)*count*num_blocks, ftello(input));exit(1);}
fclose(input);
}
*/

