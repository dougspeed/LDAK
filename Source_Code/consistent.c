/*
Copyright 2026 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

//////////////////////////

//Mainly check arguments provided are consistent (also sets a few variables)
//Note that in code, modes always in numerical order, but when printing use order from readargs

///////////////////////////

//data input - set values corresponding to dtypes (plus some checks for dtypes 2 and 5)

if(dtype==1||dtype==3||dtype==4)
{
if(genskip!=-9999||genheaders!=-9999||genprobs!=-9999)
{printf("Error, you can only use \"--gen-skip\", \"--gen-headers\" or \"--gen-probs\" with \"--gen\" or \"--mgen\"\n\n");exit(1);}
genprobs=1;
}

if(dtype==2)
{
if(genskip!=-9999||genheaders!=-9999||genprobs!=-9999)
{printf("Error, you can only use \"--gen-skip\", \"--gen-headers\" or \"--gen-probs\" with \"--gen\" or \"--mgen\"\n\n");exit(1);}

read_bgen_headers(udatafile, &bgen_preds, &bgen_samples, &bgen_comp, &bgen_layout, &bgen_ids);

if(strcmp(ufamfile,"blank")==0)
{
if(bgen_ids==0){printf("Error, %s does not contain sample IDs, so you must provide them using \"--sample\"\n\n", udatafile);exit(1);}
}
else
{
count=countrows(ufamfile)-famhead;
if(count!=bgen_samples){printf("Error, the number of rows of %s (%d) does not match the number of samples in %s (%d)\n\n", ufamfile, count, datafile, bgen_samples);exit(1);}
if(bgen_ids==1){printf("Warning, will read sample IDs from %s, and therefore ignore those provided in %s\n\n", ufamfile, udatafile);}
}

genprobs=3;
}

if(dtype==5)
{
if(genskip==-9999){genskip=0;}
if(genheaders==-9999){genheaders=5;}
if(genprobs==-9999){genprobs=3;}

if(strcmp(datalist,"blank")==0)	//some checks for --gen (do --mgen checks when reading datalist)
{
if(genheaders<4&&strcmp(ubimfile,"blank")==0)
{printf("Error, you must provide SNP details using \"--bim\"\n\n");exit(1);}

if(strcmp(ubimfile,"blank")==0&&oxchr==-9999)
{printf("Error, when using \"--gen\" you must also use \"--bim\" or \"--oxford-single-chr\"\n\n");exit(1);}

if(strcmp(ufamfile,"blank")==0)
{printf("Error, when using \"--gen\" you must also use \"--sample\"\n\n");exit(1);}
}
}

if(dtype==11)	//set values corresponding to --sp
{
if(strcmp(ubimfile,"blank")!=0)
{printf("Error, you can not use \"--bim\" with \"--sp\"\n\n");exit(1);}
if(strcmp(ufamfile,"blank")!=0)
{printf("Error, you can not use \"--fam\" or \"--sample\" with \"--sp\"\n\n");exit(1);}

sprintf(udatafile,"%s.sp",unamefile);
if(strcmp(ubimfile,"blank")==0){sprintf(ubimfile,"%s.bim",unamefile);}
if(strcmp(ufamfile,"blank")==0){sprintf(ufamfile,"%s.fam",unamefile);}
famhead=0;
if(genskip==-9999){genskip=0;}
if(genheaders==-9999){genheaders=0;}
if(genprobs==-9999){genprobs=1;}
dtype=5;
}

if(dtype==12)	//set values corresponding to --sp-gz
{
if(strcmp(ubimfile,"blank")!=0)
{printf("Error, you can not use \"--bim\" with \"--sp-gz\"\n\n");exit(1);}
if(strcmp(ufamfile,"blank")!=0)
{printf("Error, you can not use \"--fam\" or \"--sample\" with \"--sp-gz\"\n\n");exit(1);}

sprintf(udatafile,"%s.sp.gz",unamefile);
if(strcmp(ubimfile,"blank")==0){sprintf(ubimfile,"%s.bim",unamefile);}
if(strcmp(ufamfile,"blank")==0){sprintf(ufamfile,"%s.fam",unamefile);}
famhead=0;
if(genskip==-9999){genskip=0;}
if(genheaders==-9999){genheaders=0;}
if(genprobs==-9999){genprobs=1;}
dtype=5;
}

if(dtype==13)	//set values corresponding to --beagle-dose
{
if(genskip==-9999){genskip=1;}
if(genheaders==-9999){genheaders=3;}
if(genprobs==-9999){genprobs=1;}
dtype=5;
}

if(dtype==14)	//set values corresponding to --beagle-probs
{
if(genskip==-9999){genskip=1;}
if(genheaders==-9999){genheaders=3;}
if(genprobs==-9999){genprobs=3;}
dtype=5;
}

if(dtype==15)	//set values corresponding to --haps
{
if(genskip==-9999){genskip=0;}
if(genheaders==-9999){genheaders=5;}
if(genprobs==-9999){genprobs=0;}
dtype=5;
}

////////

if(strcmp(datalist,"blank")!=0&&mode!=181&&mode!=182&&mode!=183&&mode!=184&&mode!=185&&mode!=190)
{printf("Error, you can only use \"--mbfile\", \"--msp\", \"--msped\" or \"--mspeed\" when making data or with \"--calc-sim-data\"\n\n");exit(1);}

if((mode==151||mode==152||mode==153||mode==154)&&(dtype==2||dtype==3||dtype==5))
{printf("Error, you can only use \"--bfile\" or \"--speed\" with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"; you should first remake the data using either \"--make-bed\" or \"--make-speed\"\n\n");exit(1);}

if(mode==185&&genprobs<2)
{printf("Error, you can only use \"--make-gen\" when providing genotype probabilities (\"--bgen\", \"--gen\" or \"--mgen\"); consider instead using \"--make-bed\" or \"--make-speed\"\n\n");exit(1);}

////////

if(strcmp(ubimfile,"blank")!=0&&(dtype!=5||strcmp(datalist,"blank")!=0)&&mode!=135)
{printf("Error, you can only use \"--bim\" with \"--gen\" or \"--meta-analyse\"\n\n");exit(1);}

if(strcmp(ufamfile,"blank")!=0&&dtype!=2&&(dtype!=5||strcmp(datalist,"blank")!=0))
{printf("Error, you can only use \"--fam\" or \"--sample\" with \"--bgen\" or \"--gen\"\n\n");exit(1);}

if(oxchr!=-9999&&dtype!=5)
{printf("Error, you can only use \"--oxford-single-chr\" with \"--gen\" or \"--mgen\"\n\n");exit(1);}

if(oxchr!=-9999&&strcmp(ubimfile,"blank")!=0)
{printf("Error, you can not use \"--oxford-single-chr\" if using \"--bim\"\n\n");exit(1);}

if(nonsnp!=-9999&&dtype==-9999)
{printf("Error, you can only use \"--SNP-data\" when providing data\n\n");exit(1);}

if(nonsnp==1&&dtype==1)
{printf("Error, you can not use \"--SNP-data NO\" with \"--bfile\" or \"--mbfile\"; binary PLINK format handles only SNP data\n\n");exit(1);}

if(nonsnp==1&&genprobs==0)
{printf("Error, you can not use \"--SNP-data NO\" when providing haplotypes\n\n");exit(1);}

if(nonsnp==1&&genprobs>=2)
{printf("Error, you can not use \"--SNP-data NO\" when providing genotype probabilities\n\n");exit(1);}

if(nonsnp==1&&mode==156)
{printf("Error, you can not use \"--SNP-data NO\" with \"--calc-cors\"\n\n");exit(1);}

///////////////////////////

//data filtering

if((strcmp(bsampfile,"blank")!=0||strcmp(csampfile,"blank")!=0)&&mode==113)
{printf("Error, you can not use \"--keep\" or \"--remove\" with \"--join-kins\"\n\n");exit(1);}

if((strcmp(bsampfile,"blank")!=0||strcmp(csampfile,"blank")!=0)&&mode==162)
{printf("Error, you can not use \"--keep\" or \"--remove\" with \"--calc-pca-loads\" (the samples will be read from the .vect file)\n\n");exit(1);}

if(num_subs!=-9999||strcmp(subpref,"blank")!=0)
{
if(mode!=102&&mode!=104&&mode!=123&&mode!=124&&mode!=141)
{printf("Error, you can only use sample subsets with \"--calc-weights\", \"--calc-weights-all\", \"--he\" or \"--pcgc\"\n");exit(1);}

if(strcmp(subpref,"blank")==0)
{printf("Error, when using \"--subset-number\", you must also use \"--subset-prefix\"\n\n");exit(1);}

if(num_subs==-9999)
{printf("Error, when using \"--subset-prefix\", you must also use \"--subset-number\"\n\n");exit(1);}

if((mode==102||mode==104)&&num_subs==1)
{printf("Warning, using \"--subset-number 1\" and \"--subset-prefix %s\" is equivalent to (simply) using \"--keep %s1\"\n\n", subpref, subpref);}
if((mode==123||mode==124)&&num_subs==1)
{printf("Error, to estimate heritability across and within cohorts requires at least two subsets (not 1)\n\n");exit(1);}
if(mode==141&&num_subs!=2)
{printf("Error, to compute taggings across ancestries requires exactly two subsets (not %d)\n\n", num_subs);exit(1);}
}	//end of using subsets

////////

if(mode==108&&(onechr!=-9999||strcmp(onesnp,"blank")!=0))
{printf("Error, you can not use \"--chr\" or \"--snp\" with \"--find-tags\"; you should instead use \"--extract\" and/or \"--exclude\" to specify which predictors can be used as tags\n\n");exit(1);}

if(mode==117&&(onechr!=-9999||strcmp(onesnp,"blank")!=0))
{printf("Error, you can not use \"--chr\" or \"--snp\" with \"--sub-grm\"; to subtract predictors from a kinship matrix use either \"--extract\" or \"--exclude\"\n\n");exit(1);}

if(mode==117&&strcmp(bpredfile,"blank")!=0&&strcmp(cpredfile,"blank")!=0)
{printf("Error, you can not use both \"--extract\" and \"--exclude\" with \"--sub-grm\"\n\n");exit(1);}

if(mode==122&&extract==1)
{printf("Error, you can not use \"--extract\", \"--exclude\", \"--chr\" or \"--snp\" with \"--calc-blups\" (the predictors will be read from the grm.details file)\n\n");exit(1);}

if((mode==146||mode==147)&&(onechr!=-9999||strcmp(onesnp,"blank")!=0))
{printf("Error, you can not use \"--chr\" or \"--snp\" with \"--sum-hers\" or \"--sum-cors\")\n\n");exit(1);}

if(mode==149&&extract==1)
{printf("Error, you can not use \"--extract\", \"--exclude\", \"--chr\" or \"--snp\" with \"--calc-exps\" (expectations will be calculated for all predictors in the tagging file)\n\n");exit(1);}

if(mode==157&&extract==1)
{printf("Error, you can not use \"--extract\", \"--exclude\", \"--chr\" or \"--snp\" with \"--join-cors\"\n\n");exit(1);}

if(mode==162&&extract==1)
{printf("Error, you can not use \"--extract\", \"--exclude\", \"--chr\" or \"--snp\" with \"--calc-pca-loads\" (the predictors will be read from the grm.details file)\n\n");exit(1);}

if(strcmp(onesnp,"blank")!=0&&(strcmp(bpredfile,"blank")!=0||strcmp(cpredfile,"blank")!=0||onechr!=-9999))
{printf("Error, when using \"--snp\", you can not use \"--extract\", \"--exclude\" or \"--chr\"\n\n");exit(1);}

if(strcmp(onesnp,"blank")!=0&&mode!=131&&mode!=132&&mode!=134&&mode!=171&&mode!=172&&mode!=173&&mode!=181&&mode!=182&&mode!=183&&mode!=184&&mode!=185&&mode!=190)
{printf("Warning, I think it only makes sense to use \"--snp\" with \"--linear\", \"--logistic\", \"--calc-stats\", \"--calc-scores\", \"--make-phenos\", when making data or with \"--calc-sim-data\"\nBut if you can think of another use, please write in :)\n\n");}

if(mode==174&&(minmaf!=-9999||maxmaf!=-9999))
{printf("Error, when using \"--make-snps\", you should use \"--maf-low\" and \"--maf-high\" (rather than \"--min-maf\" and \"--max-maf\")\n\n");exit(1);}

if(minmaf!=-9999&&mode!=112&&mode!=114&&mode!=131&&mode!=132&&mode!=134&&mode!=151&&mode!=152&&mode!=153&&mode!=154&&mode!=156&&mode!=181&&mode!=182&&mode!=183&&mode!=184&&mode!=185&&mode!=186&&mode!=187&&mode!=188&&mode!=189)
{printf("Error, you can only use \"--min-maf\" with \"--calc-kins\", \"--calc-kins-direct\", \"--linear\", \"--logistic\", \"--ridge\", \"--bolt\", \"--bayesr\", \"--elastic\", \"--calc-cors\" or when making or condensing data; either first remake the data, or use \"--extract\" or \"--exclude\" to specify which predictors to analyse (you can use \"--calc-stats\" to calculate predictor allele frequencies)\n\n");exit(1);}

if(maxmaf!=-9999&&mode!=114&&mode!=131&&mode!=132&&mode!=134&&mode!=151&&mode!=152&&mode!=153&&mode!=154&&mode!=181&&mode!=182&&mode!=183&&mode!=184&&mode!=185&&mode!=186&&mode!=187&&mode!=188&&mode!=189)
{printf("Error, you can only use \"--max-maf\" with \"--calc-kins-direct\", \"--linear\", \"--logistic\", \"--ridge\", \"--bolt\", \"--bayesr\", \"--elastic\" or when making or condensing data; either first remake the data, or use \"--extract\" or \"--exclude\" to specify which predictors to analyse (you can use \"--calc-stats\" to calculate predictor allele frequencies)\n\n");exit(1);}

if((minmaf!=-9999||maxmaf!=-9999)&&nonsnp==1)
{printf("Error, you can not use \"--min-maf\" or \"--max-maf\" with \"--SNP-data NO\"\n\n");exit(1);}

if(minmaf!=-9999&&maxmaf!=-9999&&minmaf>=maxmaf)
{printf("Error, \"--min-maf\" (%.6f) must be lower than \"--max-maf\" (%.6f)\n\n", minmaf, maxmaf);exit(1);}

if(minvar!=-9999&&mode!=114&&mode!=131&&mode!=132&&mode!=134&&mode!=151&&mode!=152&&mode!=153&&mode!=154&&mode!=181&&mode!=182&&mode!=183&&mode!=184&&mode!=185&&mode!=186&&mode!=187&&mode!=188&&mode!=189)
{printf("Error, you can only use \"--min-var\" with \"--calc-kins-direct\", \"--linear\", \"--logistic\", \"--ridge\", \"--bolt\", \"--bayesr\", \"--elastic\" or when making or condensing data; either first remake the data, or use \"--extract\" or \"--exclude\" to specify which predictors to analyse\n\n");exit(1);}

if(minobs!=-9999&&mode!=114&&mode!=131&&mode!=132&&mode!=134&&mode!=151&&mode!=152&&mode!=153&&mode!=154&&mode!=156&&mode!=181&&mode!=182&&mode!=183&&mode!=184&&mode!=185&&mode!=186&&mode!=187&&mode!=188&&mode!=189)
{printf("Error, you can only use \"--min-obs\" with \"--calc-kins-direct\", \"--linear\", \"--logistic\", \"--ridge\", \"--bolt\", \"--bayesr\", \"--elastic\", \"--calc-cors\" or when making or condensing data; either first remake the data, or use \"--extract\" or \"--exclude\" to specify which predictors to analyse (you can use \"--calc-stats\" to calculate call-rates)\n\n");exit(1);}

if(mininfo!=-9999&&mode!=114&&mode!=131&&mode!=132&&mode!=134&&mode!=151&&mode!=152&&mode!=153&&mode!=154&&mode!=181&&mode!=182&&mode!=183&&mode!=184&&mode!=185)
{printf("Error, you can only use \"--min-info\" with \"--calc-kins-direct\", \"--linear\", \"--logistic\", \"--ridge\", \"--bolt\", \"--bayesr\", \"--elastic\" or when making data (and when providing genotype probabilities)\n\n");exit(1);}

if(mininfo!=-9999&&genprobs<2)
{printf("Error, you can only use \"--min-info\" when providing genotype probabilities (\"--bgen\", \"--gen\" or \"--mgen\"), as these are required to calculate the info score for each SNP\n\n");exit(1);}

///////////////////////////

//data scaling and coding

if(strcmp(centresfile,"blank")!=0&&mode!=112&&mode!=114)
{printf("Error, you can only use \"--predictor-means\" with \"--calc-kins\" or \"--calc-kins-direct\"\n\n");exit(1);}

if((mode==101||mode==102||mode==103||mode==104)&&strcmp(weightsfile,"blank")!=0)
{printf("Sorry, it is no longer possible to calculate weightings twice; a better solution is to first thin predictors to remove duplicates (which LDAK now does by default when cutting predictors into sections)\n\n");exit(1);}

if((mode==101||mode==102||mode==103||mode==104)&&(power!=-9999||hwestand!=-9999))
{printf("Error, you can not use \"--power\" or \"--hwe-stand\" when calculating weightings n\n");exit(1);}

if((mode==117||mode==122||mode==162)&&(strcmp(weightsfile,"blank")!=0||power!=-9999||hwestand!=-9999))
{printf("Error, you can not use \"--weights\", \"--power\" or \"--hwe-stand\" with \"--calc-blups\", \"--sub-grm\" or \"--calc-pca-loads\" (the predictor scalings will be read from the grm.details file)\n\n");exit(1);}

if(mode==160&&(strcmp(weightsfile,"blank")!=0||power!=-9999))
{printf("Error, you can not use \"--weights\" or \"--power\" with \"--validate\"\n\n");exit(1);}

if(mode==172&&strcmp(weightsfile,"blank")!=0)
{printf("Error, you can not use \"--weights\" with \"--calc-scores\" (if you wish to include weightings, you should incorporate them in the predictor effect sizes)\n\n");exit(1);}

if(strcmp(weightsfile,"blank")!=0&&ignoreweights==1){printf("Error, it is contradictory to use both \"--weights\" and \"--ignore-weights YES\"\n\n");exit(1);}

if(hwestand==1&&nonsnp==1)
{printf("Error, you can not use \"--hwe-stand YES\" with \"--SNP-data NO\"\n\n");exit(1);}

if(strcmp(pvafile,"blank")!=0&&mode!=106&&mode!=107&&mode!=136&&mode!=186&&mode!=187&&mode!=188&&mode!=189)
{printf("Error, you can only use \"--pvalues\" with \"--cut-genes\", \"--thin\", \"--thin-tops\" or when condensing data\n\n");exit(1);}

if(strcmp(impfile,"blank")!=0&&mode!=106)
{printf("Error, you can only use \"--importances\" with \"--thin\"\n\n");exit(1);}

if(strcmp(pvafile,"blank")!=0&&strcmp(impfile,"blank")!=0)
{printf("Error, you can not use both \"--pvalues\" and \"--importances\"\n\n");exit(1);}

////////

if(encoding!=-9999&&mode!=181&&mode!=182&&mode!=183&&mode!=184&&mode!=185)
{printf("Error, you can only use \"--encoding\" when making data\n\n");exit(1);}

if(encoding!=-9999&&(dtype==2||dtype==3||dtype==4||dtype==5))
{printf("Error, you can only use \"--encoding\" with \"--bfile\" or \"--mbfile\"; if your data are in a different format, you should first remake using \"--make-bed\"\n\n");exit(1);}

if(threshold!=-9999&&minprob!=-9999)
{printf("Error, you can not use both \"--threshold\" and \"--min-prob\"\n\n");exit(1);}

if((threshold!=-9999||minprob!=-9999)&&mode!=171&&mode!=181&&mode!=182&&mode!=183&&mode!=184&&mode!=185)
{printf("Error, you can only use \"--threshold\" or \"--min-prob\" with \"--calc-stats\" or when making data\n\n");exit(1);}

if(threshold!=-9999&&nonsnp==1)
{printf("Error, you can not use \"--threshold\" with \"--SNP-data NO\"\n\n");exit(1);}

if(threshold!=-9999&&dtype==1)
{printf("Error, it does not make sense to use \"--threshold\" with \"--bfile\" (because genotypes are already 0/1/2/NA)\n\n");exit(1);}

if(threshold!=-9999&&genprobs==0)
{printf("Error, it does not make sense to use \"--threshold\" when providing haplotypes (because genotypes are already 0/1/2/NA)\n\n");exit(1);}

if(minprob!=-9999&&genprobs<2)
{printf("Error, you can only use \"--min-prob\" when providing genotype probabilities (\"--bgen\", \"--gen\" or \"--mgen\")\n\n");exit(1);}

//////////////////////////

//kinships, regions, responses, summaries and fixed

if(strcmp(kinname,"blank")!=0&&strcmp(kinlist,"blank")!=0)
{printf("Error, you can not use both \"--grm\" and \"--mgrm\"\n\n");exit(1);}

if(strcmp(kinname,"blank")!=0||strcmp(kinlist,"blank")!=0)	//kinships provided (will check number correct in parsefiles.c)
{
if(mode==113)
{printf("Error, you can not use \"--grm\" or \"--mgrm\" with \"--join-kins\"; if you are looking to add kinship matrices, please use \"--add-grm\"\n\n");exit(1);}

if(mode!=115&&mode!=116&&mode!=117&&mode!=118&&mode!=119&&mode!=120&&mode!=121&&mode!=122&&mode!=123&&mode!=124&&mode!=125&&mode!=126&&mode!=131&&mode!=133&&mode!=138&&mode!=161&&mode!=162&&mode!=163&&mode!=164&&mode!=166&&mode!=167&&mode!=168&&mode!=169&&mode!=170&&mode!=177&&mode!=201&&mode!=202&&mode!=203)
{printf("Warning, kinships are provided but will not be used\n");
strcpy(kinname,"blank");strcpy(kinlist,"blank");}
}
else	//kinships not provided
{
if(mode==115||mode==118||mode==119||mode==161||mode==163||mode==164||mode==166||mode==167||mode==168||mode==169||mode==170||mode==201||mode==202||mode==203)
{printf("Error, you must provide a kinship matrix using \"--grm\"\n\n");exit(1);}
if(mode==116)
{printf("Error, you must provide one or more kinship matrices using \"--mgrm\"\n\n");exit(1);}
if(mode==117)
{printf("Error, you must provide either multiple kinship matrices using \"--mgrm\", or one kinship matrix and a list of predictors to extract or exclude using \"--grm\" with \"--extract\" or \"--exclude\"\n\n");exit(1);}
if(mode==120)
{printf("Error, you must provide two or more kinship matrices using \"--mgrm\"\n\n");exit(1);}
if(mode==122||mode==125)
{printf("Error, you must use \"--grm\" or \"--mgrm\" to provide the kinship matrix or matrices used with \"--reml\"\n\n");exit(1);}
if(mode==126)
{printf("Error, you must use \"--grm\" or \"--mgrm\" to provide one or more kinship matrices\n\n");exit(1);}
if(mode==133)
{printf("Error, you must provide two or more kinship matrices using \"--mgrm\"\n\n");exit(1);}
if(mode==162)
{printf("Error, you must use \"--grm\" to provide the kinship matrix used with \"--pca\"\n\n");exit(1);}
}

if(mode==113){sprintf(kinlist,"%spartition.list", folder);}

////////

if(kindetails==0&&(mode==113||mode==116||mode==117||mode==122||mode==162))
{printf("Error, you can not use \"--kinship-details NO\" with \"--join-kins\", \"--calc-blups\", \"--add-grm\", \"--sub-grm\" or \"--calc-pca-loads\"\n\n");exit(1);}

if(kindetails==1&&(mode==118||mode==119))
{printf("Error, you can not use \"--kinship-details YES\" with \"--convert-gz\" or \"--convert-raw\"\n\n");exit(1);}

////////

if(num_regs!=-9999||strcmp(regpref,"blank")!=0)	//regions provided
{
if(mode==117)
{printf("Sorry, you can no longer use regions with \"--sub-grm\"; instead use \"--extract\" or \"--exclude\" to specify which predictors to extract or exclude from the kinship matrix\n\n");exit(1);}

if(mode==122)
{printf("Error, you can not use \"--region-number\" or \"--region-prefix\" with \"--calc-blups\" (region details are linked to in the remlfile)\n\n");exit(1);}

if(mode==133)
{printf("Error, you can not use regions with \"--solve-null\"; you should instead convert each region into a kinship matrix\n\n");exit(1);}

if(mode!=121&&mode!=123&&mode!=124)
{printf("Error, you can only use regions with \"--reml\", \"--he\" or \"--pcgc\"\n\n");exit(1);}

if(strcmp(regpref,"blank")==0)
{printf("Error, when using \"--region-number\" you must also use \"--region-prefix\"\n\n");exit(1);}

if(num_regs==-9999)
{printf("Error, when using \"--region-prefix\" you must also use \"--region-number\"\n\n");exit(1);}

if(dtype==-9999)
{printf("Error, when using regions, you must provide a set of genetic data files using \"--bfile\", \"--bgen\", \"--sp\", \"--sped\", \"--speed\" or \"--gen\"\n\n");exit(1);}
}

if(rprune!=-9999&&strcmp(regpref,"blank")==0)
{printf("Error, you can only use \"--region-prune\" when using regions\"\n\n");exit(1);}

////////

if(strcmp(respfile,"blank")!=0&&strcmp(sumsfile,"blank")!=0)
{printf("Error, you can not use both \"--pheno\" and \"--summary\"\n\n");exit(1);}

if(strcmp(respfile,"blank")==0)	//phenotypes not provided
{
if(mode==121&&strcmp(sumsfile,"blank")==0)
{printf("Error, you must use \"--pheno\" or \"--summary\" to provide phenotypes or summary statistics\n\n");exit(1);}

if(mode==138&&strcmp(sumsfile,"blank")==0)
{printf("Error, you must use \"--pheno\" or \"--summary\" to provide phenotypes or summary statistics\n\n");exit(1);}

if((mode==172&&strcmp(finalfile,"blank")!=0)&&strcmp(sumsfile,"blank")==0)
{printf("Error, you must use \"--pheno\" or \"--summary\" to provide phenotypes or summary statistics\n\n");exit(1);}

if(mode==123||mode==124||mode==126||mode==127||mode==128||mode==129||mode==130||mode==229||mode==230||mode==131||mode==132||mode==133||mode==134||mode==151||mode==152||mode==153||mode==154||mode==160||mode==194)
{printf("Error, you must use \"--pheno\" to provide phenotypes\n\n");exit(1);}
}

if(mpheno!=-9999&&strcmp(respfile,"blank")==0&&mode!=140)
{printf("Error, you can only use \"--mpheno\" if providing phenotypes\n\n");exit(1);}

if(mpheno2!=-9999&&strcmp(respfile,"blank")==0)
{printf("Error, you can only use \"--mpheno2\" if providing phenotypes\n\n");exit(1);}

if(mpheno==-1&&mode!=121&&mode!=123&&mode!=124&&mode!=126&&mode!=129&&mode!=130&&mode!=131&&mode!=132&&mode!=134&&mode!=140&&mode!=151&&mode!=152&&mode!=153&&mode!=154&&mode!=163)
{printf("Error, you can only use \"--mpheno ALL\" with \"--reml\", \"--he\", \"--pcgc\", \"--fast-reml\", \"--quant-her\", \"--tetra-her\", \"--linear\", \"--logistic\", \"--kvik-step3\", \"--ridge\", \"--bolt\", \"--bayesr\", \"--elastic\" or \"--decompose\"\n\n");exit(1);}

if(mpheno==-1&&mode==131&&(strcmp(kinname,"blank")!=0||strcmp(kinlist,"blank")!=0))
{printf("Error, you can not use \"--mpheno ALL\" if providing kinships\n\n");exit(1);}

if(mpheno2!=-9999&&mode!=229&&mode!=230)
{printf("Error, you can only use \"--mpheno2\" with \"--quant-bivar\" or \"--tetra-bivar\"\n\n");exit(1);}

if(mpheno2!=-9999&&mpheno==-1)
{printf("Error, you can not use \"--mpheno2\" with \"--mpheno ALL\"\n\n");exit(1);}

if(mpheno!=-9999&&mpheno2!=-9999&&mpheno2==mpheno)
{printf("Error, \"--mpheno\" and \"--mpheno2\" should specify different phenotypes\n\n");exit(1);}

if(strcmp(npheno,"blank")!=0&&mpheno!=-9999)
{printf("Error, you can not use both \"--pheno-name\" and \"--mpheno\"\n\n");exit(1);}

if(pad!=-9999&&mode!=121&&mode!=123&&mode!=124&&mode!=126&&mode!=127&&mode!=128&&mode!=131&&mode!=132&&mode!=134&&mode!=138&&mode!=151&&mode!=152&&mode!=153&&mode!=154&&mode!=163&&mode!=194)
{printf("Error, you can only use \"--dentist\" with \"--reml\", \"--he\", \"--pcgc\", \"--fast-reml\", \"--fast-he\", \"--fast-pcgc\", \"--linear\", \"--logistic\", \"--calc-genes-reml\", \"--ridge\", \"--bolt\", \"--bayesr\", \"--elastic\", \"--decompose\" or \"--solve-gre\"\n\n");exit(1);}

if(pad==1&&mode==131&&(strcmp(kinname,"blank")!=0||strcmp(kinlist,"blank")!=0))
{printf("Error, you can not use \"--dentist YES\" if providing kinships\n\n");exit(1);}

////////

if(strcmp(sumsfile,"blank")!=0)	//summaries provided
{
if(mode!=121&&mode!=138&&mode!=146&&mode!=147&&mode!=150&&mode!=158&&mode!=159&&mode!=172&&mode!=179&&mode!=180)
{printf("Error, you can only use \"--summary\" with \"--reml\", \"--calc-genes-reml\", \"--sum-hers\", \"--sum-cors\", \"--mega-prs\", \"--giga-prs\", \"--pseudo-summaries\",  \"--calc-scores\" or \"--impute-summaries\"\n\n");exit(1);}

if(strcmp(kinname,"blank")!=0||strcmp(kinlist,"blank")!=0)
{printf("Error, you can not use \"--summary\" if providing kinships\n\n");exit(1);}

if(mode==121&&num_regs==-9999)
{printf("Error, you can only use \"--summary\" when using regions\n\n");exit(1);}

if(pad==1)
{printf("Error, you can not use \"--dentist YES\" with \"--summary\"\n\n");exit(1);}

if(nonsnp==1)
{printf("Error, you can not \"--summary\" with \"--SNP-data NO\"\n\n");exit(1);}
}

if(strcmp(sums2file,"blank")!=0&&mode==159)
{printf("Sorry, you can no longer use \"--summary2\" with \"--mega-prs\"; you should instead use \"--pseudos\" to provide training and test summary statistics\n\n");exit(1);}

if(strcmp(sums2file,"blank")!=0&&mode!=147)
{printf("Error, you can only use \"--summary2\" with \"--sum-cors\"\n\n");exit(1);}

if(strcmp(sumslist,"blank")!=0&&mode!=135&&mode!=159)
{printf("Error, you can only use \"--sumslist\" with \"--meta-analyse\", \"--mega-prs\" or \"--giga-prs\"\n\n");exit(1);}

if(strcmp(sumslist,"blank")!=0&&strcmp(sumsfile,"blank")!=0)
{printf("Error, you can not use both \"--summary\" and \"--sumslist\"\n\n");exit(1);}

if(strcmp(bulkfile,"blank")!=0&&mode!=159)
{printf("Error, you can only use \"--bulk-summaries\" with \"--mega-prs\" or \"--giga-prs\"\n\n");exit(1);}

////////

if(amb!=-9999&&strcmp(sumsfile,"blank")==0&&strcmp(sumslist,"blank")==0&&mode!=140)
{printf("Error, you can only use \"--allow-ambiguous\" with \"--summary\" or \"--sumslist\"\n\n");exit(1);}

if(scaling!=-9999&&strcmp(sumsfile,"blank")==0)
{printf("Error, you can only use \"--scaling\" with \"--summary\"\n\n");exit(1);}

if(scaling2!=-9999&&strcmp(sumsfile,"blank")==0)
{printf("Error, you can only use \"--scaling2\" with \"--summary2\"\n\n");exit(1);}

////////

if(prev!=-9999&&strcmp(wayfile,"blank")!=0)
{printf("Error, you can not use \"--prevalence\" with \"--pathways\"\n\n");exit(1);}

if(prev!=-9999&&mode!=121&&mode!=123&&mode!=124&&mode!=126&&mode!=127&&mode!=128&&mode!=129&&mode!=130&&mode!=229&&mode!=230&&mode!=138&&mode!=140&&mode!=146&&mode!=147&&mode!=173&&mode!=176&&mode!=194)
{printf("Error, you can only use \"--prevalence\" with \"--reml\", \"--he\", \"--pcgc\", \"--fast-reml\", \"--fast-he\", \"--fast-pcgc\", \"--quant-her\", \"--tetra-her\", \"--calc-genes-reml\", \"--sum-hers\", \"--sum-cors\", \"--make-phenos\", \"--jackknife\" or \"--solve-gre\"\n\n");exit(1);}

if(prev2!=-9999&&mode!=147&&mode!=230)
{printf("Error, you can only use \"--prevalence2\" with \"--sum-cors\"\n\n");exit(1);}

if(ascer!=-9999&&strcmp(sumsfile,"blank")==0)
{printf("Error, you can only use \"--ascertainment\" with \"--summary\"\n\n");exit(1);}

if(ascer!=-9999&&prev==-9999)
{printf("Error, you can only use \"--ascertainment\" with \"--prevalence\"\n\n");exit(1);}

if(ascer2!=-9999&&strcmp(sums2file,"blank")==0)
{printf("Error, you can only use \"--ascertainment2\" with \"--summary2\"\n\n");exit(1);}

if(ascer2!=-9999&&prev2==-9999)
{printf("Error, you can only use \"--ascertainment2\" with \"--prevalence2\"\n\n");exit(1);}

if(prev!=-9999&&mode==147&&prev2==-9999)
{printf("Error, if using \"--prevalence\", you must also use \"--prevalence2\"\n\n");exit(1);}

if(prev2!=-9999&&mode==147&&prev==-9999)
{printf("Error, if using \"--prevalence2\", you must also use \"--prevalence\"\n\n");exit(1);}

if(prev!=-9999&&strcmp(sumsfile,"blank")!=0&&ascer==-9999)
{printf("Error, you must use \"--ascertainment\" to specify the proportion of samples who were cases when generating the summary statistics in %s\n\n", sumsfile);exit(1);}

if(prev2!=-9999&&strcmp(sums2file,"blank")!=0&&ascer2==-9999)
{printf("Error, you must use \"--ascertainment2\" to specify the proportion of samples who were cases when generating the summary statistics in %s\n\n", sums2file);exit(1);}

////////

if((mode==169||mode==170)&&strcmp(covarfile,"blank")!=0)
{printf("Error, when using \"--gxemm-iid\" or \"--gxemm-free\" you must use \"--enviro\" instead of \"--covar\"\n\n");exit(1);}

if(strcmp(covarfile,"blank")!=0&&mode!=121&&mode!=122&&mode!=123&&mode!=124&&mode!=125&&mode!=126&&mode!=127&&mode!=128&&mode!=129&&mode!=130&&mode!=229&&mode!=230&&mode!=131&&mode!=132&&mode!=133&&mode!=134&&mode!=138&&mode!=151&&mode!=152&&mode!=153&&mode!=154&&mode!=156&&mode!=164&&mode!=172&&mode!=173&&mode!=175&&mode!=194)
{printf("Error, you can only use \"--covar\" with \"--reml\", \"--calc-blup\", \"--he\", \"--pcgc\", \"--reml-predict\", \"--fast-reml\", \"--fast-he\", \"--fast-pcgc\", \"--quant-her\", \"--tetra-her\", \"--linear\", \"--logistic\", \"--solve-null\", \"--calc-genes-reml\", \"--ridge\", \"--bolt\", \"--bayesr\", \"--elastic\", \"--calc-cors\", \"--adjust-grm\", \"--calc-scores\", \"--make-phenos\", \"--calc-inflation\" or \"--solve-gre\"\n\n");exit(1);}

if(strcmp(covarfile,"blank")!=0&&strcmp(sumsfile,"blank")!=0)
{printf("Error, you can not use \"--covar\" with \"--summary\"\n\n");exit(1);}

if(strcmp(covarnums,"blank")!=0&&strcmp(covarfile,"blank")==0)
{printf("Error, you can only use \"--covar-numbers\" if using \"--covar\"\n\n");exit(1);}

if(strcmp(covarnames,"blank")!=0&&strcmp(covarfile,"blank")==0)
{printf("Error, you can only use \"--covar-names\" if using \"--covar\"\n\n");exit(1);}

if(strcmp(covarnums,"blank")!=0&&strcmp(covarnames,"blank")!=0)
{printf("Error, you can not use both \"--covar-cols\" and \"--covar-names\"\n\n");exit(1);}

if(strcmp(envfile,"blank")!=0&&mode!=121&&mode!=122&&mode!=123&&mode!=124&&mode!=131&&mode!=164&&mode!=169&&mode!=170&&mode!=173)
{printf("Error, you can only use \"--enviro\" with \"--reml\", \"--calc-blup\", \"--he\", \"--pcgc\", \"--linear\", \"--adjust-grm\", \"--gxemm-iid\", \"--gxemm-free\", \"--calc-scores\" or \"--make-phenos\"\n\n");exit(1);}

if(strcmp(envfile,"blank")!=0&&strcmp(sumsfile,"blank")!=0)
{printf("Error, you can not use \"--enviro\" with \"--summary\"\n\n");exit(1);}

if(strcmp(envfile,"blank")!=0&&mode==131&&mpheno==-1)
{printf("Error, you can not use \"--enviro\" with \"--mpheno ALL\"\n\n");exit(1);}

if(strcmp(envfile,"blank")!=0&&mode==131&&pad==1)
{printf("Error, you can not use \"--enviro\" with \"--dentist YES\"\n\n");exit(1);}

if(strcmp(envfile,"blank")!=0&&mode==131&&(strcmp(kinname,"blank")!=0||strcmp(kinlist,"blank")!=0))
{printf("Error, you can not use \"--enviro\" if providing a kinship matrix\n\n");exit(1);}

if(strcmp(topfile,"blank")!=0&&mode==122)
{printf("Error, you can not use \"--top-preds\" with \"--calc-blups\" (top predictors details are linked to in the remlfile)\n\n");exit(1);}

if(strcmp(topfile,"blank")!=0&&mode==172)
{printf("Error, you can not use \"--top-preds\" with \"--calc-scores\"; the top predictors should instead be included in the score file\n\n");exit(1);}

if(strcmp(topfile,"blank")!=0&&mode!=121&&mode!=123&&mode!=124&&mode!=126&&mode!=127&&mode!=128&&mode!=131&&mode!=132&&mode!=133&&mode!=134&&mode!=138&&mode!=164)
{printf("Error, you can only use \"--top-preds\" with \"--reml\", \"--he\", \"--pcgc\", \"--fast-reml\", \"--fast-he\", \"--fast-pcgc\", \"--linear\", \"--logistic\", \"--solve-null\", \"--calc-genes-reml\" or \"--adjust-grm\"\n\n");exit(1);}

if(strcmp(topfile,"blank")!=0&&mode==131&&(strcmp(kinname,"blank")!=0||strcmp(kinlist,"blank")!=0))
{printf("Error, you can not use \"--top-preds\" if providing a kinship matrix\n\n");exit(1);}

if(strcmp(topfile,"blank")!=0&&strcmp(sumsfile,"blank")!=0)
{printf("Error, you can not use \"--top-preds\" with \"--summary\"\n\n");exit(1);}

if(strcmp(topfile,"blank")!=0&&mode==131&&strcmp(envfile,"blank")!=0)
{printf("Error, you can not use \"--top-preds\" with \"--enviro\"\n\n");exit(1);}

if(strcmp(topfile,"blank")!=0&&strcmp(prsfile,"blank")!=0)
{printf("Error, you can not use \"--top-preds\" with \"--PRS\"\n\n");exit(1);}

if(strcmp(topfile,"blank")!=0&&strcmp(transfile,"blank")!=0)
{printf("Error, you can not use \"--top-preds\" with \"--transform\"\n\n");exit(1);}

if(strcmp(topfile,"blank")!=0&&dtype==-9999)
{printf("Error, when using \"--top-preds\" you must provide a set of genetic data files using \"--bfile\", \"--bgen\", \"--sp\", \"--sped\", \"--speed\" or \"--gen\"\n\n");exit(1);}

if(strcmp(factorfile,"blank")!=0&&(mode==122||mode==125))
{printf("Error, you can not use \"--factors\" with \"--calc-blups\" or \"--reml-predict\"; you should instead use \"--covar\" to provide the combined covariate file produced when running REML\n\n");exit(1);}

if(strcmp(factorfile,"blank")!=0&&mode!=121&&mode!=123&&mode!=124&&mode!=126&&mode!=127&&mode!=128&&mode!=129&&mode!=130&&mode!=229&&mode!=230&&mode!=131&&mode!=132&&mode!=134&&mode!=138&&mode!=151&&mode!=152&&mode!=153&&mode!=154&&mode!=156&&mode!=173&&mode!=175&&mode!=194)
{printf("Error, you can only use \"--factors\" with \"--reml\", \"--he\", \"--pcgc\", \"--fast-reml\", \"--fast-he\", \"--fast-pcgc\", \"--quant-her\", \"--tetra-her\", \"--linear\", \"--logistic\", \"--calc-genes-reml\", \"--ridge\", \"--bolt\", \"--bayesr\", \"--elastic\", \"--calc-cors\", \"--make-phenos\", \"--calc-inflation\" or \"--solve-gre\"\n\n");exit(1);}

if(strcmp(factorfile,"blank")!=0&&strcmp(sumsfile,"blank")!=0)
{printf("Error, you can not use \"--factors\" with \"--summary\"\n\n");exit(1);}

////////

if(strcmp(rovarfile,"blank")!=0&&mode!=121)
{printf("Error, you can only use \"--random-covar\" with \"--reml\"\n\n");exit(1);}

if(strcmp(rovarfile,"blank")!=0&&strcmp(sumsfile,"blank")!=0)
{printf("Error, you can not use \"--random-covar\" with \"--summary\"\n\n");exit(1);}

if(strcmp(ractorfile,"blank")!=0&&mode!=121)
{printf("Error, you can only use \"--random-factors\" with \"--reml\"\n\n");exit(1);}

if(strcmp(ractorfile,"blank")!=0&&strcmp(sumsfile,"blank")!=0)
{printf("Error, you can not use \"--random-factors\" with \"--summary\"\n\n");exit(1);}

////////

if(strcmp(offsetfile,"blank")!=0&&mode!=132)
{printf("Error, you can only use \"--offset\" with \"--logistic\"\n\n");exit(1);}

if(strcmp(offsetfile,"blank")!=0&&mpheno==-1)
{printf("Error, you can not use \"--offset\" with \"--mpheno ALL\"\n\n");exit(1);}

//////////////////////////

//calculating weights, thinning and finding/removing tags

if(nothin!=-9999&&mode!=101)
{printf("Error, you can only use \"--no-thin\" with \"--cut-weights\"\n\n");exit(1);}

if(wprune!=-9999&&mode!=101&&mode!=106&&mode!=107&&mode!=110)
{printf("Error, you can only use \"--window-prune\" with \"--cut-weights\", \"--thin\", \"--thin-tops\" or \"--thin-common\"\n\n");exit(1);}

if(wprune!=-9999&&(nothin==1||nothin==2))
{printf("Error, it does not make sense to use \"--window-prune\" with \"--no-thin YES\" or \"--no-thin DONE\"\n\n");exit(1);}

if(window_kb!=-9999&&mode!=101&&mode!=102&&mode!=104&&mode!=106&&mode!=107&&mode!=108&&mode!=109&&mode!=110&&mode!=141&&mode!=152&&mode!=153&&mode!=154&&mode!=156)
{printf("Error, you can only use \"--window-kb\" with \"--cut-weights\", \"--calc-weights\", \"--calc-weights-all\", \"--calc-tagging\", \"--bolt\", \"--bayesr\", \"--elastic\", \"--calc-cors\", \"--thin\", \"--thin-tops\", \"--find-tags\", \"--remove-tags\" or \"--thin-common\"\n\n");exit(1);}

if(window_length!=-9999&&mode!=101&&mode!=102&&mode!=104&&mode!=106&&mode!=107&&mode!=110)
{printf("Error, you can only use \"--window-length\" with \"--cut-weights\", \"--calc-weights\", \"--calc-weights-all\", \"--thin\", \"--thin-tops\" or \"--thin-common\"\n\n");exit(1);}

if(window_cm!=-9999&&mode!=101&&mode!=102&&mode!=104&&mode!=106&&mode!=107&&mode!=108&&mode!=109&&mode!=110&&mode!=141&&mode!=152&&mode!=153&&mode!=154&&mode!=156)
{printf("Error, you can only use \"--window-cm\" with \"--cut-weights\", \"--calc-weights\", \"--calc-weights-all\", \"--calc-tagging\", \"--bolt\", \"--bayesr\", \"--elastic\", \"--calc-cors\", \"--thin\", \"--thin-tops\", \"--find-tags\", \"--remove-tags\" or \"--thin-common\"\n\n");exit(1);}

if(window_kb!=-9999&&window_length!=-9999)
{printf("Error, you can not use both \"--window-kb\" and \"--window-length\"\n\n");exit(1);}

if(window_kb!=-9999&&window_cm!=-9999)
{printf("Error, you can not use both \"--window-kb\" and \"--window-cm\"\n\n");exit(1);}

if(window_length!=-9999&&window_cm!=-9999)
{printf("Error, you can not use both \"--window-length\" and \"--window-cm\"\n\n");exit(1);}

if(window_cm!=-9999&&dtype==2)
{printf("Error, you can not use \"--window-cm\" with \"--bgen\" (the bgen file does not contain genetic distances)\n\n");exit(1);}

if(window_cm!=-9999&&dtype==5&&strcmp(ubimfile,"blank")==0)
{printf("Error, when using \"--window-cm\" with \"--gen\", you must also use \"--bim\" to provide genetic distances\n\n");exit(1);}

////////

if(section_kb!=-9999&&mode!=101)
{printf("Error, you can only use \"--section-kb\" with \"--cut-weights\"\n\n");exit(1);}

if(section_length!=-9999&&mode!=101)
{printf("Error, you can only use \"--section-length\" with \"--cut-weights\"\n\n");exit(1);}

if(section_cm!=-9999&&mode!=101)
{printf("Error, you can only use \"--section-cm\" with \"--cut-weights\"\n\n");exit(1);}

if(section_kb!=-9999&&section_length!=-9999)
{printf("Error, you can not use both \"--section-kb\" and \"--section-length\"\n\n");exit(1);}

if(section_kb!=-9999&&section_cm!=-9999)
{printf("Error, you can not use both \"--section-kb\" and \"--section-cm\"\n\n");exit(1);}

if(section_length!=-9999&&section_cm!=-9999)
{printf("Error, you can not use both \"--section-length\" and \"--section-cm\"\n\n");exit(1);}

if(section_cm!=-9999&&dtype==2)
{printf("Error, you can not use \"--section-cm\" with \"--bgen\" (the bgen file does not contain genetic distances)\n\n");exit(1);}

if(section_cm!=-9999&&dtype==5&&strcmp(ubimfile,"blank")==0)
{printf("Error, when using \"--section-cm\" with \"--gen\", you must also use \"--bim\" to provide genetic distances\n\n");exit(1);}

if(section_kb!=-9999&&window_cm!=-9999)
{printf("Error, you can not use both \"--window-cm\" and \"--section-kb\" (generally, you should use either \"--window-kb\" and \"--section-kb\" or \"--window-cm\" and \"--section-cm\")\n\n");exit(1);}

if(section_cm!=-9999&&window_kb!=-9999)
{printf("Error, you can not use both \"--window-kb\" and \"--section-cm\" (generally, you should use either \"--window-kb\" and \"--section-kb\" or \"--window-cm\" and \"--section-cm\")\n\n");exit(1);}

if(section_cm!=-9999&&window_cm==-9999)
{printf("Error, when using \"--section-cm\" you must also use \"--window-cm\"\n\n");exit(1);}

////////

if(buffer_kb!=-9999&&mode!=101)
{printf("Error, you can only use \"--buffer-kb\" with \"--cut-weights\"\n\n");exit(1);}

if(buffer_length!=-9999&&mode!=101)
{printf("Error, you can only use \"--buffer-length\" with \"--cut-weights\"\n\n");exit(1);}

if(buffer_cm!=-9999&&mode!=101)
{printf("Error, you can only use \"--buffer-cm\" with \"--cut-weights\"\n\n");exit(1);}

if(buffer_kb!=-9999&&buffer_length!=-9999)
{printf("Error, you can not use both \"--buffer-kb\" and \"--buffer-length\"\n\n");exit(1);}

if(buffer_kb!=-9999&&buffer_cm!=-9999)
{printf("Error, you can not use both \"--buffer-kb\" and \"--buffer-cm\"\n\n");exit(1);}

if(buffer_length!=-9999&&buffer_cm!=-9999)
{printf("Error, you can not use both \"--buffer-length\" and \"--buffer-cm\"\n\n");exit(1);}

if(buffer_cm!=-9999&&dtype==2)
{printf("Error, you can not use \"--buffer-cm\" with \"--bgen\" (the bgen file does not contain genetic distances)\n\n");exit(1);}

if(buffer_cm!=-9999&&dtype==5&&strcmp(ubimfile,"blank")==0)
{printf("Error, when using \"--buffer-cm\" with \"--gen\", you must also use \"--bim\" to provide genetic distances\n\n");exit(1);}

if(buffer_kb!=-9999&&window_cm!=-9999)
{printf("Error, you can not use \"--buffer-kb\" with \"--window-cm\"\n\n");exit(1);}

if(buffer_cm!=-9999&&window_kb!=-9999)
{printf("Error, you can not use \"--buffer-cm\" with \"--window-kb\"\n\n");exit(1);}

if(buffer_cm!=-9999&&window_cm==-9999)
{printf("Error, when using \"--buffer-cm\" you must also use \"--window-cm\"\n\n");exit(1);}

////////

if(section!=-9999&&mode==104)
{printf("Error, you can not use \"--section\" with \"--calc-weights-all\" (LDAK will loop through all sections)\n\n");exit(1);}

if(section!=-9999&&mode!=102)
{printf("Error, you can only use \"--section\" with \"--calc-weights\"\n\n");exit(1);}

if(section_start!=-9999&&mode!=104)
{printf("Error, you can only use \"--start-section\" with \"--calc-weights-all\"\n\n");exit(1);}

////////

if((lddecay!=-9999||halflife!=-9999)&&mode!=102&&mode!=104)
{printf("Error, you can only use \"--decay\" or \"--half-life\" with \"--calc-weights\" or \"--calc-weights-all\"\n\n");exit(1);}

if(halflife!=-9999&&lddecay==0)
{printf("Error, you can only use \"--half-life\" with \"--decay YES\"\n\n");exit(1);}

////////

if(fudge!=-9999&&mode!=102&&mode!=104)
{printf("Error, you can only use \"--quick-weights\" with \"--calc-weights\" or \"--calc-weights-all\"\n\n");exit(1);}

if(simplex!=-9999&&mode!=102&&mode!=104)
{printf("Error, you can only use \"--simplex\" with \"--calc-weights\" or \"--calc-weights-all\"\n\n");exit(1);}

if(maxtime!=-9999&&mode!=102&&mode!=104)
{printf("Error, you can only use \"--max-time\" with \"--calc-weights\" or \"--calc-weights-all\"\n\n");exit(1);}

if(maxtime!=-9999&&simplex!=1)
{printf("Error, you can only use \"--max-time\" with \"--simplex YES\"\n\n");exit(1);}

if(spread!=-9999&&mode!=103&&mode!=104)
{printf("Error, you can only use \"--spread\" with \"--join-weights\" or \"--calc-weights-all\"\n\n");exit(1);}

////////

if(strcmp(targetfile,"blank")!=0&&mode!=109)
{printf("Error, you can only use \"--targets\" with \"--remove-tags\"\n\n");exit(1);}

//////////////////////////

//calculating and manipulating kinships

if(part_length!=-9999&&mode!=111&&mode!=136&&mode!=191)
{printf("Error, you can only use \"--partition-length\" with \"--cut-kins\", \"--cut-genes\" or \"--cut-gre\"\n\n");exit(1);}

if(bychr!=-9999&&mode!=111&&mode!=136&&kvikstep!=2&&gctastep!=2&&faststep!=2)
{printf("Error, you can only use \"--by-chr\" with \"--cut-kins\" or \"--cut-genes\", or when using \"--kvik-step2\", \"--GCTA-LOCO-step2\" or \"--fastGWA-step2\"\n\n");exit(1);}

if(part_length!=-9999&&bychr==1)
{printf("Error, you can not use both \"--partition-length\" and \"--by-chr YES\"\n\n");exit(1);}

if(num_parts!=-9999||strcmp(partpref,"blank")!=0)	//partitions provided
{
if(gigaprs==1)
{printf("Error, you can not use partitions with \"--giga-prs\"\n\n");}

if(mode!=105&&mode!=111&&mode!=127&&mode!=128&&mode!=141&&mode!=145&&mode!=159)
{printf("Error, you can only use partitions with \"--cut-kins\", \"--fast-he\", \"--fast-pcgc\", \"--calc-tagging\", \"--calc-overlaps\" or \"--mega-prs\"\n\n");exit(1);}

if(num_subs!=-9999)
{printf("Error, you can not use partitions when using sample subsets\n\n");exit(1);}

if(part_length!=-9999||bychr==1)
{printf("Error, you can not use \"--partition-length\" or \"--by-chr YES\" when using partitions\n\n");exit(1);}

if(strcmp(partpref,"blank")==0)
{printf("Error, when using \"--partition-number\" you must also use \"--partition-prefix\"\n\n");exit(1);}

if(num_parts==-9999)
{printf("Error, when using \"--partition-prefix\" you must also use \"--partition-number\"\n\n");exit(1);}

if(num_parts==1&&mode==111)
{printf("Warning, using \"--partition-number 1\" and \"--partition-prefix %s\" is equivalent to (simply) using \"--calc-kins-direct\" with \"--extract %s1\"\n\n", partpref, partpref);}
}

if(checkpart!=-9999&&mode!=111)
{printf("Error, you can only use \"--check-partitions\" with \"--cut-kins\"\n\n");exit(1);}

if(checkpart!=-9999&&strcmp(partpref,"blank")==0)
{printf("Error, you can only use \"--check-partitions\" when using partitions\n\n");exit(1);}

////////

if(partition!=-9999&&mode!=112&&mode!=137&&mode!=138&&mode!=141&&mode!=192)
{printf("Error, you can only use \"--partition\" with \"--calc-kins\", \"--calc-genes-kins\", \"--calc-genes-reml\", \"--calc-tagging\" or \"--calc-gre\"\n\n");exit(1);}

if((kingz!=-9999||kinraw!=-9999)&&mode!=112&&mode!=113&&mode!=114&&mode!=116&&mode!=117&&mode!=118&&mode!=119&&mode!=133&&mode!=137&&mode!=164&&mode!=166&&mode!=167&&mode!=168&&mode!=169&&mode!=170)
{printf("Error, you can only use \"--kinship-gz\" or \"--kinship-raw\" when producing kinship matrices\n\n");exit(1);}

if(single!=-9999&&mode!=112&&mode!=114)
{printf("Error, you can only use \"--single\" with \"--calc-kins\" or \"--calc-kins-direct\"\n\n");exit(1);}

////////

if(dosage!=-9999&&mode!=112&&mode!=114)
{printf("Error, you can only use \"--male-dosage\" with \"--calc-kins\" or \"--calc-kins-direct\"\n\n");exit(1);}

if(strcmp(malesfile,"blank")!=0&&mode!=112&&mode!=114)
{printf("Error, you can only use \"--males\" with \"--calc-kins\" or \"--calc-kins-direct\"\n\n");exit(1);}

if(strcmp(malesfile,"blank")!=0&&dosage==-9999)
{printf("Error, you can only use \"--males\" if also using \"--male-dosage\"\n\n");exit(1);}

if(onlydets!=-9999&&mode!=112&&mode!=114)
{printf("Error, you can only use \"--only-details\" with \"--calc-kins\" or \"--calc-kins-direct\"\n\n");exit(1);}

if(strcmp(invsfile,"blank")!=0&&mode!=114)
{printf("Error, you can only use \"--inverse\" with \"--calc-kins-direct\"\n\n");exit(1);}

if(strcmp(invsfile,"blank")!=0&&dosage!=-9999)
{printf("Error, you can not use \"--inverse\" with \"--male-dosage\"\n\n");exit(1);}

if(square!=-9999&&mode!=112&&mode!=114)
{printf("Error, you can only use \"--square\" with \"--calc-kins\" or \"--calc-kins-direct\"\n\n");exit(1);}

////////

if((maxrel!=-9999||minrel!=-9999)&&mode!=115)
{printf("Error, you can only use \"--max-rel\" or \"--min-rel\" with \"--filter\"\n\n");exit(1);}

if(kinstand!=-9999&&mode!=115)
{printf("Error, you can only use \"--kin-stand\" with \"--filter\"\n\n");exit(1);}

////////

if(partial!=-9999&&mode!=118)
{printf("Error, you can only use \"--partial\" with \"--convert-gz\"\n\n");exit(1);}

///////////////////////////

//reml, blup, he and sum-hers

if(diagonal!=-9999&&mode!=121)
{printf("Error, you can only use \"--diagonal\" with \"--reml\"\n\n");exit(1);}

if(diagonal==1&&strcmp(kinname,"blank")==0&&strcmp(kinlist,"blank")==0)
{printf("Error, you can only use \"--diagonal YES\" when providing one or more kinship matrices\n\n");exit(1);}

if(diagonal==1&&(strcmp(rovarfile,"blank")!=0||strcmp(ractorfile,"blank")!=0))
{printf("Error, you can not use \"--diagonal YES\" with \"--random-covars\" or \"--random-factors\"\n\n");exit(1);}

if(diagonal==1&&num_regs!=-9999)
{printf("Error, you can not use \"--diagonal YES\" when using regions\n\n");exit(1);}

if(strcmp(hersfile,"blank")!=0&&mode!=121&&mode!=126&&mode!=133)
{printf("Error, you can only use \"--starts\" with \"--reml\", \"--fast-reml\" or \"--solve-null\"\n\n");exit(1);}

if(strcmp(hersfile,"blank")!=0&&(strcmp(rovarfile,"blank")!=0||strcmp(ractorfile,"blank")!=0))
{printf("Error, you can not use \"--starts\" with \"--random-covars\" or \"--random-factors\"\n\n");exit(1);}

if(strcmp(hersfile,"blank")!=0&&strcmp(kinname,"blank")==0&&strcmp(kinlist,"blank")==0&&num_regs==-9999)
{printf("Error, you can only use \"--starts\" when providing kinship matrices and/or regions\n\n");exit(1);}

if(hestart!=-9999&&mode!=121&&mode!=126)
{printf("Error, you can only use \"--he-starts\" with \"--reml\" or \"--fast-reml\"\n\n");exit(1);}

if(hestart==0&&strcmp(hersfile,"blank")!=0)
{printf("Error, you can not use \"--he-starts YES\" with \"--starts\"\n\n");exit(1);}

if(shortcut==0&&mode!=121)
{printf("Error, you can only use \"--shortcut NO\" with \"--reml\"\n\n");exit(1);}

if(shortcut==0&&strcmp(sumsfile,"blank")!=0)
{printf("Error, you can not use \"--shortcut NO\" with \"--summary\"\n\n");exit(1);}

if(shortcut==0&&diagonal==1)
{printf("Error, you can not use \"--shortcut NO\" with \"--diagonal YES\"\n\n");exit(1);}

if(discenv!=-9999&&mode!=121&&mode!=123&&mode!=124&&mode!=170)
{printf("Error, you can only use \"--subgroups\" with \"--reml\", \"--he\", \"--pcgc\" or \"--gxemm-free\"\n\n");exit(1);}

if(discenv==1&&(strcmp(rovarfile,"blank")!=0||strcmp(ractorfile,"blank")!=0))
{printf("Error, you can not use \"--subgroups YES\" with \"--random-covars\" or \"--random-factors\"\n\n");exit(1);}

if(discenv==1&&strcmp(envfile,"blank")==0)
{printf("Error, you can only use \"--subgroups YES\" with \"--enviro\"\n\n");exit(1);}

if(strcmp(oversfile,"blank")!=0&&mode!=121&&mode!=123&&mode!=124&&mode!=126)
{printf("Error, you can only use \"--overlaps\" with \"--reml\", \"--he\", \"--pcgc\" or \"--fast-reml\"\n\n");exit(1);}

if(strcmp(oversfile,"blank")!=0&&strcmp(kinname,"blank")==0&&strcmp(kinlist,"blank")==0)
{printf("Error, you can only use \"--overlaps\" when providing one or more kinship matrices\n\n");exit(1);}

if(strcmp(oversfile,"blank")!=0&&num_regs!=-9999)
{printf("Error, you can not use \"--overlaps\" when using regions\n\n");exit(1);}

////////

if(strcmp(remlfile,"blank")!=0&&mode!=122&&mode!=125)
{printf("Error, you can only use \"--remlfile\" with \"--calc-blups\" or \"--reml-predict\"\n\n");exit(1);}

////////

if(adjusted!=-9999&&mode!=123&&mode!=124)
{printf("Error, you can only use \"--adjusted\" with \"--he\" or \"--pcgc\"\n\n");exit(1);}

if(trun!=-9999&&mode!=123&&mode!=124)
{printf("Error, you can only use \"--truncate-fly\" with \"--he\" or \"--pcgc\"\n\n");exit(1);}

if(trun!=-9999&&strcmp(kinname,"blank")==0&&strcmp(kinlist,"blank")==0)
{printf("Error, you can only use \"--truncate-fly\" when providing a kinship matrix\n\n");exit(1);}

////////

if(num_vects!=-9999&&mode!=126&&mode!=127&&mode!=128)
{printf("Error, you can only use \"--repetitions\" with \"--fast-reml\", \"--fast-he\" and \"--fast-pcgc\"\n\n");exit(1);}

if(ldlt!=-9999&&mode!=126)
{printf("Error, you can only use \"--LDLT\" with \"--fast-reml\"\n\n");exit(1);}

////////

if(strcmp(relfile,"blank")!=0&&mode!=129&&mode!=130&&mode!=229&&mode!=230&&drm==-9999&&fastgwa!=1)
{printf("Error, you can only use \"--relatives\" with \"--quant-her\", \"--tetra-her\", \"--fastgwa YES\" or with \"--linear\" and \"--DRM\"\n\n");exit(1);}

///////////////////////////

//association analysis

if(strcmp(prsfile,"blank")!=0&&mode!=131&&mode!=132&&mode!=134&&mode!=140)
{printf("Error, you can only use \"--PRS\" with \"--linear\" or \"--logistic\"\n\n");exit(1);}

if(strcmp(prsfile,"blank")!=0&&strcmp(topfile,"blank")!=0)
{printf("Error, you can not use \"--PRS\" with \"--top-preds\"\n\n");exit(1);}

if(strcmp(prsfile,"blank")!=0&&mode==131&&(strcmp(kinname,"blank")!=0||strcmp(kinlist,"blank")!=0))
{printf("Error, you can not use \"--PRS\" if providing a kinship matrix\n\n");exit(1);}

if(strcmp(prsfile,"blank")!=0&&mode==131&&strcmp(envfile,"blank")!=0)
{printf("Error, you can not use \"--PRS\" with \"--enviro\"\n\n");exit(1);}

if(checkprs!=-9999&&strcmp(prsfile,"blank")==0)
{printf("Error, you can only use \"--check-PRS\" with \"--PRS\"\n\n");exit(1);}

////////

if(spatest!=-9999&&mode!=131&&mode!=132&&mode!=134)
{printf("Error, you can only use \"--spa-test\" with \"--linear\" or \"--logistic\"\n\n");exit(1);}

if(spatest==1&&mode==131&&(strcmp(kinname,"blank")!=0||strcmp(kinlist,"blank")!=0))
{printf("Error, you can not use \"--spa-test YES\" if providing a kinship matrix\n\n");exit(1);}

if(spatest==1&&mode==131&&strcmp(envfile,"blank")!=0)
{printf("Error, you can not use \"--spa-test YES\" with \"--enviro\"\n\n");exit(1);}

if(num_knots!=-9999&&spatest!=1)
{printf("Error, you can only use \"--num-knots\" with \"--spa-test YES\"\n\n");exit(1);}

if(num_bins!=-9999&&spatest!=1)
{printf("Error, you can only use \"--num-bins\" with \"--spa-test YES\"\n\n");exit(1);}

if(spathresh!=-9999&&spatest!=1)
{printf("Error, you can only use \"--spa-threshold\" with \"--spa-test YES\"\n\n");exit(1);}

if(spamax!=-9999&&spatest!=1)
{printf("Error, you can only use \"--spa-range\" with \"--spa-test YES\"\n\n");exit(1);}

////////

if(families!=-9999&&mode!=131)
{printf("Error, you can only use \"--families\" with \"--linear\"\n\n");exit(1);}

if(families==1&&(strcmp(kinname,"blank")!=0||strcmp(kinlist,"blank")!=0))
{printf("Error, you can not use \"--families YES\" if providing a kinship matrix\n\n");exit(1);}

if(families==1&&mpheno==-1)
{printf("Error, you can not use \"--families YES\" with \"--mpheno ALL\"\n\n");exit(1);}

if(families==1&&pad==1)
{printf("Error, you can not use \"--families YES\" with \"--dentist YES\"\n\n");exit(1);}

if(families==1&&strcmp(topfile,"blank")!=0)
{printf("Error, you can not use \"--families YES\" with \"--top-preds\"\n\n");exit(1);}

if(families==1&&strcmp(envfile,"blank")!=0)
{printf("Error, you can not use \"--families YES\" with \"--enviro\"\n\n");exit(1);}

if(families==1&&strcmp(prsfile,"blank")!=0)
{printf("Error, you can not use \"--families YES\" with \"--PRS\"\n\n");exit(1);}

if(families==1&&spatest==1)
{printf("Error, you can not use \"--families YES\" with \"--spa-test YES\"\n\n");exit(1);}

if(trios!=-9999&&mode!=131&&mode!=174)
{printf("Error, you can only use \"--trios\" with \"--linear\" or \"--make-snps\"\n\n");exit(1);}

if(trios==1&&(strcmp(kinname,"blank")!=0||strcmp(kinlist,"blank")!=0))
{printf("Error, you can not use \"--trios YES\" if providing a kinship matrix\n\n");exit(1);}

if(trios==1&&mpheno==-1)
{printf("Error, you can not use \"--trios YES\" with \"--mpheno ALL\"\n\n");exit(1);}

if(trios==1&&pad==1)
{printf("Error, you can not use \"--trios YES\" with \"--dentist YES\"\n\n");exit(1);}

if(trios==1&&strcmp(topfile,"blank")!=0)
{printf("Error, you can not use \"--trios YES\" with \"--top-preds\"\n\n");exit(1);}

if(trios==1&&strcmp(envfile,"blank")!=0)
{printf("Error, you can not use \"--trios YES\" with \"--enviro\"\n\n");exit(1);}

if(trios==1&&strcmp(prsfile,"blank")!=0)
{printf("Error, you can not use \"--trios YES\" with \"--PRS\"\n\n");exit(1);}

if(trios==1&&spatest==1)
{printf("Error, you can not use \"--trios YES\" with \"--spa-test YES\"\n\n");exit(1);}

if(duos!=-9999&&mode!=131)
{printf("Error, you can only use \"--duos\" with \"--linear\"\n\n");exit(1);}

if(duos!=-9999&&(strcmp(kinname,"blank")!=0||strcmp(kinlist,"blank")!=0))
{printf("Error, you can not use \"--duos\" if providing a kinship matrix\n\n");exit(1);}

if(duos==1&&mpheno==-1)
{printf("Error, you can not use \"--duos\" with \"--mpheno ALL\"\n\n");exit(1);}

if(duos!=-9999&&pad==1)
{printf("Error, you can not use \"--duos\" with \"--dentist YES\"\n\n");exit(1);}

if(duos!=-9999&&strcmp(topfile,"blank")!=0)
{printf("Error, you can not use \"--duos\" with \"--top-preds\"\n\n");exit(1);}

if(duos!=-9999&&strcmp(envfile,"blank")!=0)
{printf("Error, you can not use \"--duos\" with \"--enviro\"\n\n");exit(1);}

if(duos!=-9999&&strcmp(prsfile,"blank")!=0)
{printf("Error, you can not use \"--duos\" with \"--PRS\"\n\n");exit(1);}

if(duos!=-9999&&spatest==1)
{printf("Error, you can not use \"--duos\" with \"--spa-test YES\"\n\n");exit(1);}

if(families==1&&trios==1)
{printf("Error, you can not use both \"--families YES\" and \"--trios YES\"\n\n");exit(1);}

if(families==1&&duos!=-9999)
{printf("Error, you can not use both \"--families YES\" and \"--duos\"\n\n");exit(1);}

if(trios==1&&duos!=-9999)
{printf("Error, you can not use both \"--trios YES\" and \"--duos\"\n\n");exit(1);}

////////

if(drm!=-9999&&mode!=131)
{printf("Error, you can only use \"--DRM\" with \"--linear\"\n\n");exit(1);}

if(drm!=-9999&&(dtype==2||dtype==3||dtype==4||dtype==5||strcmp(datalist,"blank")!=0))
{printf("Error, you can only use \"--DRM\" with \"--bfile\"; if your data are in a different format, you should first remake using \"--make-bed\"\n\n");exit(1);}

if(drm!=-9999&&(strcmp(kinname,"blank")!=0||strcmp(kinlist,"blank")!=0))
{printf("Error, you can not use \"--DRM\" if providing a kinship matrix\n\n");exit(1);}

if(drm!=-9999&&mpheno==-1)
{printf("Error, you can not use \"--DRM\" with \"--mpheno ALL\"\n\n");exit(1);}

if(drm!=-9999&&pad==1)
{printf("Error, you can not use \"--DRM\" with \"--dentist YES\"\n\n");exit(1);}

if(drm!=-9999&&strcmp(topfile,"blank")!=0)
{printf("Error, you can not use \"--DRM\" with \"--top-preds\"\n\n");exit(1);}

if(drm!=-9999&&strcmp(envfile,"blank")!=0)
{printf("Error, you can not use \"--DRM\" with \"--enviro\"\n\n");exit(1);}

if(drm!=-9999&&strcmp(prsfile,"blank")!=0)
{printf("Error, you can not use \"--DRM\" with \"--PRS\"\n\n");exit(1);}

if(drm!=-9999&&spatest==1)
{printf("Error, you can not use \"--DRM\" with \"--spa-test YES\"\n\n");exit(1);}

if(drm!=-9999&&families==1)
{printf("Error, you can not use \"--DRM\" with \"--families YES\"\n\n");exit(1);}

if(drm!=-9999&&trios==1)
{printf("Error, you can not use \"--DRM\" with \"--trios YES\"\n\n");exit(1);}

if(drm!=-9999&&duos!=-9999)
{printf("Error, you can not use \"--DRM\" with \"--duos\"\n\n");exit(1);}

if(drm==1&&onechr==23)
{printf("Error, it is contradictory to use \"--DRM AUTOSOMES\" and \"--chr X\" (or \"--chr 23\")\n\n");exit(1);}

if(drm==2&&onechr!=-9999&&onechr!=23)
{printf("Error, you can not use \"--chr\" with \"--DRM CHRX\" (LDAK will automatically reduce to CHR X SNPs)\n\n");exit(1);}

if(xpar!=-9999&&drm!=2)
{printf("Error, you can only use \"--split-par\" with \"--DRM CHRX\"\n\n");exit(1);}

if(strcmp(sexfile,"blank")!=0&&drm!=2)
{printf("Error, you can only use \"--sexfile\" with \"--DRM CHRX\"\n\n");exit(1);}

////////

if(adjpreds!=-9999&&mode!=131&&mode!=132&&mode!=134)
{printf("Error, you can only use \"--adjust-predictors\" with \"--linear\" or \"--logistic\"\n\n");exit(1);}

if(adjpreds!=-9999&&strcmp(covarfile,"blank")==0&&strcmp(topfile,"blank")==0&&strcmp(factorfile,"blank")==0)
{printf("Error, you can only use \"--adjust-predictors\" if using \"--covar\", \"--top-preds\" and/or \"--factors\"\n\n");exit(1);}

if(adjpreds!=-9999&&mode==131&&(strcmp(kinname,"blank")!=0||strcmp(kinlist,"blank")!=0))
{printf("Error, you can not use \"--adjust-predictors\" if providing a kinship matrix\n\n");exit(1);}

if(adjpreds!=-9999&&mode==131&&strcmp(envfile,"blank")!=0)
{printf("Error, you can not use \"--adjust-predictors\" with \"--enviro\"\n\n");exit(1);}

if(adjpreds!=-9999&&families==1)
{printf("Error, you can not use \"--adjust-predictors\" with \"--families YES\"\n\n");exit(1);}

if(adjpreds!=-9999&&trios==1)
{printf("Error, you can not use \"--adjust-predictors\" with \"--trios YES\"\n\n");exit(1);}

if(adjpreds!=-9999&&duos!=-9999)
{printf("Error, you can not use \"--adjust-predictors\" with \"--duos\"\n\n");exit(1);}

if(strcmp(sampwfile,"blank")!=0&&(mode==151||mode==152||mode==153||mode==154))
{printf("Sorry, you can not use \"--sample-weights\" with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\" (I briefly added the option, but decided it was not sufficiently robust)\n\n");exit(1);}

if(strcmp(sampwfile,"blank")!=0&&mode!=114&&mode!=121&&mode!=126&&mode!=131&&mode!=134&&mode!=173)
{printf("Error, you can only use \"--sample-weights\" with \"--calc-kins-direct\", \"--reml\", \"--fast-reml\", \"--linear\" or \"--make-phenos\"\n\n");exit(1);}

if(strcmp(sampwfile,"blank")!=0&&mode==131&&(strcmp(kinname,"blank")!=0||strcmp(kinlist,"blank")!=0))
{printf("Error, you can not use \"--sample-weights\" if providing a kinship matrix\n\n");exit(1);}

if(strcmp(sampwfile,"blank")!=0&&strcmp(sumsfile,"blank")!=0)
{printf("Error, you can not use \"--sample-weights\" with \"--summary\"\n\n");exit(1);}

if(strcmp(sampwfile,"blank")!=0&&mode==131&&strcmp(envfile,"blank")!=0)
{printf("Error, you can not use \"--sample-weights\" with \"--enviro\"\n\n");exit(1);}

if(strcmp(sampwfile,"blank")!=0&&spatest==1)
{printf("Error, you can not use \"--sample-weights\" with \"--spa-test YES\"\n\n");exit(1);}

if(strcmp(sampwfile,"blank")!=0&&families==1)
{printf("Error, you can not use \"--sample-weights\" with \"--families YES\"\n\n");exit(1);}

if(strcmp(sampwfile,"blank")!=0&&trios==1)
{printf("Error, you can not use \"--sample-weights\" with \"--trios YES\"\n\n");exit(1);}

if(strcmp(sampwfile,"blank")!=0&&duos!=-9999)
{printf("Error, you can not use \"--sample-weights\" with \"--duos\"\n\n");exit(1);}

if(strcmp(sampwfile,"blank")!=0&&adjpreds!=-9999)
{printf("Error, you can not use \"--sample-weights\" with \"--adjust-predictors\"\n\n");exit(1);}

if(sandwich!=-9999&&mode!=131&&mode!=134)
{printf("Error, you can only use \"--sandwich\" with \"--linear\"\n\n");exit(1);}

if(sandwich==1&&(strcmp(kinname,"blank")!=0||strcmp(kinlist,"blank")!=0))
{printf("Error, you can not use \"--sandwich YES\" if providing a kinship matrix\n\n");exit(1);}

if(sandwich==1&&strcmp(envfile,"blank")!=0)
{printf("Error, you can not use \"--sandwich YES\" with \"--enviro\"\n\n");exit(1);}

if(sandwich==0&&families==1)
{printf("Error, you can not use \"--sandwich NO\" with \"--families YES\"\n\n");exit(1);}

if(sandwich==0&&trios==1)
{printf("Error, you can not use \"--sandwich NO\" with \"--trios YES\"\n\n");exit(1);}

if(sandwich==0&&duos!=-9999)
{printf("Error, you can not use \"--sandwich NO\" with \"--duos\"\n\n");exit(1);}

if(exact!=-9999&&mode!=131)
{printf("Error, you can only use \"--exact\" with \"--linear\"\n\n");exit(1);}

if(exact!=-9999&&strcmp(kinname,"blank")==0&&strcmp(kinlist,"blank")==0)
{printf("Error, you can only use \"--exact\" when providing a kinship matrix\n\n");exit(1);}

if(scoretest!=-9999&&mode!=132)
{printf("Error, you can only use \"--score-test\" with \"--logistic\"\n\n");exit(1);}

if(scoretest==0&&mpheno==-1)
{printf("Error, you can not use \"--score-test NO\" with \"--mpheno ALL\"\n\n");exit(1);}

if(scoretest==0&&pad==1)
{printf("Error, you can not use \"--score-test NO\" with \"--dentist YES\"\n\n");exit(1);}

if(scoretest==0&&strcmp(prsfile,"blank")!=0)
{printf("Error, you can not use \"--score-test NO\" with \"--PRS\"\n\n");exit(1);}

if(scoretest==0&&spatest==1)
{printf("Error, you can not use \"--score-test NO\" with \"--spa-test YES\"\n\n");exit(1);}

if(scoretest==0&&adjpreds!=-9999)
{printf("Error, you can not use \"--score-test NO\" with \"--adjust-predictors\"\n\n");exit(1);}

if(strcmp(transfile,"blank")!=0&&mode!=131)
{printf("Error, you can only use \"--transform\" with \"--linear\"\n\n");exit(1);}

if(strcmp(transfile,"blank")!=0&&mpheno!=-1)
{printf("Error, you can only use \"--transform\" with \"--mpheno ALL\"\n\n");exit(1);}

if(strcmp(transfile,"blank")!=0&&(strcmp(kinname,"blank")!=0||strcmp(kinlist,"blank")!=0))
{printf("Error, you can not use \"--transform\" if providing a kinship matrix\n\n");exit(1);}

if(strcmp(transfile,"blank")!=0&&strcmp(envfile,"blank")!=0)
{printf("Error, you can not use \"--transform\" with \"--enviro\"\n\n");exit(1);}

if(strcmp(transfile,"blank")!=0&&strcmp(prsfile,"blank")!=0)
{printf("Error, you can not use \"--transform\" with \"--PRS\"\n\n");exit(1);}

if(strcmp(transfile,"blank")!=0&&families==1)
{printf("Error, you can not use \"--transform\" with \"--families YES\"\n\n");exit(1);}

if(strcmp(transfile,"blank")!=0&&trios==1)
{printf("Error, you can not use \"--transform\" with \"--trios YES\"\n\n");exit(1);}

if(strcmp(transfile,"blank")!=0&&duos!=-9999)
{printf("Error, you can not use \"--transform\" with \"--duos\"\n\n");exit(1);}

////////

if((strcmp(genefile,"blank")!=0||chunks!=-9999||chunksbp!=-9999)&&mode!=136&&mode!=140&&mode!=185&&mode!=186&&mode!=187&&mode!=188&&mode!=189)
{printf("Error, you can only use \"--genefile\", \"--chunks\" or \"--chunks-bp\" with \"--cut-genes\" or when condensing data\n\n");exit(1);}

if((gene_buffer!=-9999||up_buffer!=-9999||down_buffer!=-9999)&&mode!=136&&mode!=140&&mode!=185&&mode!=186&&mode!=187&&mode!=188&&mode!=189)
{printf("Error, you can only use \"--gene-buffer\", \"--up-buffer\" or \"--down-buffer\" with \"--cut-genes\" or when condensing data\n\n");exit(1);}

if((gene_buffer!=-9999||up_buffer!=-9999||down_buffer!=-9999)&&strcmp(genefile,"blank")==0)
{printf("Error, you can only use \"--gene-buffer\", \"--up-buffer\" or \"--down-buffer\" when using \"--genefile\"\n\n");exit(1);}

if((up_buffer!=-9999||down_buffer!=-9999)&&gene_buffer!=-9999)
{printf("Error, you can not use \"--up-buffer\" or \"--down-buffer\" if using \"--gene-buffer\"\n\n");exit(1);}

if(up_buffer!=-9999&&down_buffer==-9999)
{printf("Error, when using \"--up-buffer\", you must also use \"--down-buffer\"\n\n");exit(1);}

if(down_buffer!=-9999&&up_buffer==-9999)
{printf("Error, when using \"--down-buffer\", you must also use \"--up-buffer\"\n\n");exit(1);}

if(minweight!=-9999&&mode!=136&&mode!=140&&mode!=185&&mode!=186&&mode!=187&&mode!=188&&mode!=189)
{printf("Error, you can only use \"--min-weights\" with \"--cut-genes\" or when condensing data\n\n");exit(1);}

if(minweight!=-9999&&chunks!=-9999)
{printf("Error, you can not use \"--min-weight\" when using \"--chunks\" (the genome will be divided so that all chunks have equal weight)\n\n");exit(1);}

if(overlap!=-9999&&mode!=136&&mode!=140&&mode!=185&&mode!=186&&mode!=187&&mode!=188&&mode!=189)
{printf("Error, you can only use \"--overlap\" with \"--cut-genes\" or when condensing data\n\n");exit(1);}

if(strcmp(pathfile,"blank")!=0&&mode!=136)
{printf("Error, you can only use \"--pathfile\" with \"--cut-genes\"\n\n");exit(1);}

if(strcmp(pathfile,"blank")!=0&&part_length!=-9999){printf("Error, you can not use \"--partition-length\" with \"--pathfile\"\n\n");exit(1);}

if(strcmp(pathfile,"blank")!=0&&bychr==1){printf("Error, you can not use \"--by-chr YES\" with \"--pathfile\"\n\n");exit(1);}

if(strcmp(pathfile,"blank")!=0&&strcmp(genefile,"blank")==0)
{printf("Error, you can only use \"--pathfile\" when using \"--genefile\"\n\n");exit(1);}

////////

if(gene_perms!=-9999&&mode!=138&&mode!=140)
{printf("Error, you can only use \"--gene-permutations\" with \"--calc-genes-reml\"\n\n");exit(1);}

if(gprune!=-9999&&mode!=138&&mode!=140)
{printf("Error, you can only use \"--gene-prune\" with \"--calc-genes-reml\"\n\n");exit(1);}

if(saveall!=-9999&&mode!=138&&mode!=140)
{printf("Error, you can only use \"--save-all\" with \"--calc-genes-reml\"\n\n");exit(1);}

if(limit!=-9999&&mode!=138&&mode!=140)
{printf("Error, you can only use \"--her-limit\" with \"--calc-genes-reml\"\n\n");exit(1);}

////////

if(magma!=-9999&&mode!=139&&mode!=140)
{printf("Error, you can only use \"--MAGMA\" with \"--join-genes-reml\"\n\n");exit(1);}

if(cut1!=-9999||cut2!=-9999)
{
if(mode!=139&&mode!=140)
{printf("Error, you can only use \"--sig1\" and \"--sig2\" with \"--join-genes-reml\"\n\n");exit(1);}

if(cut1==-9999||cut2==-9999)
{printf("Error, to identify significant genes/chunks, you must use both \"--sig1\" and \"--sig2 (typical values are 0.00001 and 0.01)\n\n");exit(1);}

if(cut2<cut1){printf("Warning, \"--sig2\" (%f) is lower than \"--sig1\" (%f), which can result in overlapping regions\n\n", cut2, cut2);}
}

if(gamp!=-9999||gam1!=-9999||gam2!=-9999)
{
if(mode!=139&&mode!=140)
{printf("Error, you can only use \"--gamma-fraction\", \"--gamma-alpha\" or \"--gamma-beta\" with \"--join-genes-reml\"\n\n");exit(1);}
if(gamp==-9999||gam1==-9999||gam2==-9999)
{printf("Error, to specify the NULL distribution of LRT statistics, you must provide \"--gamma-fraction\", \"--gamma-alpha\" and \"--gamma-beta\"\n\n");exit(1);}
}

///////////////////////////

//sumher

if(num_anns!=-9999||strcmp(annpref,"blank")!=0)	//annotations provided
{
if(gigaprs==1)
{printf("Error, you can not use annotations with \"--giga-prs\"\n\n");}

if(mode!=127&&mode!=128&&mode!=141&&mode!=145&&mode!=156&&mode!=159)
{printf("Error, you can only use annotations with \"--fast-he\", \"--fast-pcgc\", \"--calc-tagging\", \"--calc-overlaps\", \"--calc-cors\" or \"--mega-prs\"\n\n");exit(1);}

if(num_subs!=-9999)
{printf("Error, you can not use annotations when using sample subsets\n\n");exit(1);}

if(num_parts!=-9999||strcmp(partpref,"blank")!=0)
{printf("Error, you can not use both partitions and annotations; generally, you should use partitions when the lists of predictors are disjoint, and annotations when there is overlap\n\n");exit(1);}

if(strcmp(annpref,"blank")==0)
{printf("Error, when using \"--annotation-number\" you must also use \"--annotation-prefix\"\n\n");exit(1);}

if(num_anns==-9999)
{printf("Error, when using \"--annotation-prefix\" you must also use \"--annotation-number\"\n\n");exit(1);}
}

if(strcmp(genpreds,"blank")!=0&&mode!=141)
{printf("Error, you can only use \"--genic-predictors\" with \"--calc-taggings\"\n\n");exit(1);}

if(strcmp(genpreds,"blank")==0&&mode==141&&partition!=-9999)
{printf("Error, you can only use \"--partition\" when also using \"--genic-predictors\"\n\n");exit(1);}

if(strcmp(genpreds,"blank")!=0&&num_subs!=-9999)
{printf("Error, you can not use \"--genic-predictors\" when using sample subsets\n\n");exit(1);}

if(strcmp(genpreds,"blank")!=0&&strcmp(annpref,"blank")==0)
{printf("Error, you can only use \"--genic-predictors\" when using annotations\n\n");exit(1);}

if(strcmp(labfile,"blank")!=0&&mode!=127&&mode!=128&&mode!=141&&mode!=145&&mode!=156)
{printf("Error, you can only use \"--labels\" with \"--fast-he\", \"--fast-pcgc\", \"--calc-taggings\", \"--calc-overlaps\" or \"--calc-cors\"\n\n");exit(1);}

if(strcmp(labfile,"blank")!=0&&mode==156&&strcmp(annpref,"blank")==0)
{printf("Error, you can only use \"--labels\" when using annotations\n\n");exit(1);}

if(strcmp(labfile,"blank")!=0&&strcmp(partpref,"blank")==0&&strcmp(annpref,"blank")==0)
{printf("Error, you can only use \"--labels\" when using partitions or annotations\n\n");exit(1);}

if(mode==141&&strcmp(genpreds,"blank")!=0&&strcmp(labfile,"blank")==0)
{printf("Error, when using \"--genic-predictors\" you must also use \"--labels\" to provide the pathway names\n\n");exit(1);}

if(mode==156&&strcmp(annpref,"blank")!=0&&strcmp(labfile,"blank")==0)
{printf("Error, when using \"--annotations\" you must also use \"--labels\" to provide the annotation names\n\n");exit(1);}

if(backpart!=-9999&&strcmp(partpref,"blank")==0)
{printf("Error, you can only use \"--background\" when using partitions\n\n");exit(1);}

if(backpart!=-9999&&mode!=105&&mode!=127&&mode!=128&&mode!=141)
{printf("Error, you can only use \"--background\" with \"--fast-he\", \"--fast-pcgc\" or \"--calc-tagging\"\n\n");exit(1);}

if(allone==1&&strcmp(partpref,"blank")==0)
{printf("Error, you can only use \"--all-one YES\" when using partitions\n\n");exit(1);}

if(allone!=-9999&&mode!=105&&mode!=127&&mode!=128&&mode!=141&&mode!=159)
{printf("Error, you can only use \"--all-one YES\" with \"--fast-he\", \"--fast-pcgc\", \"--calc-tagging\" or \"--mega-prs\"\n\n");exit(1);}

if(reduce!=-9999&&mode!=105&&mode!=141)
{printf("Error, you can only use \"--reduce\" with \"--calc-tagging\"\n\n");exit(1);}

if(reduce==1&&strcmp(genpreds,"blank")!=0)
{printf("Error, you can not use \"--reduce YES\" with \"--genic-predictors\"\n\n");exit(1);}

////////

if(strcmp(printfile,"blank")!=0&&mode!=141)
{printf("Error, you can only use \"--regression-predictors\" with \"--calc-tagging\"\n\n");exit(1);}

if(strcmp(herfile,"blank")!=0&&mode!=141)
{printf("Error, you can only use \"--heritability-predictors\" with \"--calc-tagging\"\n\n");exit(1);}

if(unbias!=-9999&&mode!=141)
{printf("Error, you can only use \"--unbiased\" with \"--calc-tagging\"\n\n");exit(1);}

if(unbias==1&&num_subs!=-9999)
{printf("Error, you can not use \"--unbiased YES\" when using sample subsets\n\n");exit(1);}

if(cover!=-9999&&mode!=141)
{printf("Error, you can only use \"--coverage\" with \"--calc-tagging\"\n\n");exit(1);}

if(cover==1&&num_subs!=-9999)
{printf("Error, you can not use \"--coverage YES\" when using sample subsets\n\n");exit(1);}

if(cover==1&&strcmp(genpreds,"blank")!=0)
{printf("Error, you can not use \"--coverage YES\" with \"--genic-predictors\"\n\n");exit(1);}

if(fourdp!=-9999&&mode!=141)
{printf("Error, you can only use \"--full-accuracy\" with \"--calc-tagging\"\n\n");exit(1);}

if(cover==1&&strcmp(genpreds,"blank")==0)
{printf("Error, you can only use \"--full-accuracy YES\" with \"--genic-predictors\"\n\n");exit(1);}

////////

if(strcmp(taglist,"blank")!=0&&mode!=142&&mode!=143)
{printf("Error, you can only use \"--taglist\" with \"--join-tagging\" or \"--merge-tagging\"\n\n");exit(1);}

if(strcmp(matlist,"blank")!=0&&mode!=142)
{printf("Error, you can only use \"--matlist\" with \"--join-tagging\"\n\n");exit(1);}

if(strcmp(matlist,"blank")!=0&&strcmp(taglist,"blank")==0)
{printf("Error, you can only use \"--matlist\" with \"--taglist\"\n\n");exit(1);}

if(strcmp(xaglist,"blank")!=0&&mode!=142)
{printf("Error, you can only use \"--cross-taglist\" with \"--join-tagging\"\n\n");exit(1);}

if(strcmp(taglist,"blank")!=0&&strcmp(xaglist,"blank")!=0)
{printf("Error, you can not use both \"--taglist\" with \"--cross-taglist\"\n\n");exit(1);}

if(strcmp(pathlist,"blank")!=0&&mode!=142&&mode!=143)
{printf("Error, you can only use \"--pathlist\" with \"--join-tagging\" or \"--merge-tagging\"\n\n");exit(1);}

if(strcmp(taglist,"blank")!=0&&strcmp(pathlist,"blank")!=0)
{printf("Error, you can not use both \"--taglist\" with \"--pathlist\"\n\n");exit(1);}

if(checkdups!=-9999&&mode!=142)
{printf("Error, you can only use \"--check-dups\" with \"--join-tagging\"\n\n");exit(1);}

////////

if(strcmp(tagfile,"blank")!=0&&mode!=144&&mode!=146&&mode!=147&&mode!=149)
{printf("Error, you can only use \"--tagfile\" with \"--reduce-tagging\", \"--sum-hers\", \"--sum-cors\" or \"--calc-exps\"\n\n");exit(1);}

if(strcmp(xagfile,"blank")!=0&&mode!=147)
{printf("Error, you can only use \"--cross-tagfile\" with \"--sum-cors\"\n\n");exit(1);}

if(strcmp(wayfile,"blank")!=0&&mode!=146)
{printf("Error, you can only use \"--pathways\" with \"--sum-hers\"\n\n");exit(1);}

if(strcmp(altfile,"blank")!=0&&mode!=146&&mode!=147)
{printf("Error, you can only use \"--alternative-tags\" with \"--sum-hers\" or \"--sum-cors\"\n\n");exit(1);}

if(strcmp(cvsfile,"blank")!=0&&mode!=146)
{printf("Error, you can only use \"--cv-predictors\" with \"--sum-hers\"\n\n");exit(1);}

if(strcmp(cvsfile,"blank")!=0&&strcmp(wayfile,"blank")!=0)
{printf("Error, you can not use \"--cv-predictors\" with \"--pathways\"\n\n");exit(1);}

if(strcmp(catfile,"blank")!=0&&mode!=144&&mode!=146)
{printf("Error, you can only use \"--categories\" with \"--reduce-tagging\" or \"--sum-hers\"\n\n");exit(1);}

if(strcmp(catfile,"blank")!=0&&strcmp(wayfile,"blank")!=0)
{printf("Error, you can not use \"--categories\" with \"--pathways\"\n\n");exit(1);}

if(strcmp(taufile,"blank")!=0&&mode!=146&&mode!=149)
{printf("Error, you can only use \"--taus\" with \"--sum-hers\" or \"--calc-exps\"\n\n");exit(1);}

if(strcmp(taufile,"blank")!=0&&strcmp(wayfile,"blank")!=0)
{printf("Error, you can not use \"--taus\" with \"--pathways\"\n\n");exit(1);}

if(strcmp(catfile,"blank")!=0&&strcmp(taufile,"blank")!=0)
{printf("Error, you can not use both \"--categories\" and \"--taus\"\n\n");exit(1);}

////////

if(checksums!=-9999&&mode!=135&&mode!=138&&mode!=146&&mode!=147&&mode!=158)
{printf("Error, you can only use \"--check-sums\" with \"--meta-analyse\", \"--calc-genes-reml\", \"--sum-hers\", \"--sum-cors\" or \"--pseudo-summaries\"\n\n");exit(1);}

if(checksize!=-9999&&mode!=159)
{printf("Error, you can only use \"--check-sums\" with \"--mega-prs\" or \"--giga-prs\"\n\n");exit(1);}

if(checksize!=-9999&&strcmp(sumslist,"blank")==0)
{printf("Error, you can only use \"--check-size\" with \"--sumslist\"\n\n");exit(1);}

if(checkmcmc!=-9999&&mode!=159)
{printf("Error, you can only use \"--check-MCMC\" with \"--mega-prs\" or \"--giga-prs\"\n\n");exit(1);}

if(gcon!=-9999&&mode!=146&&mode!=147)
{printf("Error, you can only use \"--genomic-control\" with \"--sum-hers\" or \"--sum-cors\"\n\n");exit(1);}

if(cept!=-9999&&mode!=146&&mode!=147)
{printf("Error, you can only use \"--intercept\" with \"--sum-hers\" or \"--sum-cors\"\n\n");exit(1);}

if(oversamp!=-9999&&mode!=147)
{printf("Error, you can only use \"--overlapping-samples\" with \"--sum-cors\"\n\n");exit(1);}

////////

if(ldsc!=-9999&&mode!=146)
{printf("Error, you can only use \"--LDSC\" with \"--sum-hers\"\n\n");exit(1);}

if(ldsc==1&&strcmp(wayfile,"blank")!=0)
{printf("Error, you can not use \"--LDSC YES\" with \"--pathways\"\n\n");exit(1);}

if(ldsc==1&&(strcmp(catfile,"blank")!=0||strcmp(taufile,"blank")!=0))
{printf("Error, you can not use \"--LDSC YES\" with \"--categories\" or \"--taus\"\n\n");exit(1);}

if(chisol!=-9999&&mode!=146)
{printf("Error, you can only use \"--chisq-solver\" with \"--sum-hers\"\n\n");exit(1);}

if(ldsc==1&&chisol==1)
{printf("Error, you can not both \"--LDSC YES\" and \"--chisq-solver YES\"\n\n");exit(1);}

if(tagone!=-9999&&mode!=146&&mode!=147)
{printf("Error, you can only use \"--tag-one\" with \"--sum-hers\" or \"--sum-cors\"\n\n");exit(1);}

if(tagone==1&&strcmp(wayfile,"blank")!=0)
{printf("Error, you can not use \"--tag-one YES\" with \"--pathways\"\n\n");exit(1);}

////////

if(divide!=-9999&&mode!=146)
{printf("Error, you can only use \"--divisions\" with \"--sum-hers\"\n\n");exit(1);}

if(divide!=-9999&&strcmp(wayfile,"blank")!=0)
{printf("Error, you can not use \"--divisions\" with \"--pathways\"\n\n");exit(1);}

if(divide!=-9999&&(strcmp(catfile,"blank")!=0||strcmp(taufile,"blank")!=0))
{printf("Error, you can not use \"--divisions\" with \"--categories\" or \"--taus\"\n\n");exit(1);}

if(divide!=-9999&&ldsc==1)
{printf("Error, you can not use \"--divisions\" with \"--LDSC YES\"\n\n");exit(1);}

if(uptaus!=-9999&&mode!=146)
{printf("Error, you can only use \"--update-taus\" with \"--sum-hers\"\n\n");exit(1);}

if(uptaus==1&&strcmp(wayfile,"blank")!=0)
{printf("Error, you can not use \"--update-taus YES\" with \"--pathways\"\n\n");exit(1);}

if(uptaus==1&&divide==-9999)
{printf("Error, you can only use \"--update-taus YES\" with \"--divisions\"\n\n");exit(1);}

if(strcmp(powfile,"blank")!=0&&mode!=146&&mode!=151&&mode!=152&&mode!=153&&mode!=154)
{printf("Error, you can only use \"--powerfile\" with \"--sum-hers\", \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n");exit(1);}

if(strcmp(powfile,"blank")!=0&&strcmp(wayfile,"blank")!=0)
{printf("Error, you can not use \"--powerfile\" with \"--pathways\"\n\n");exit(1);}

if(strcmp(powfile,"blank")!=0&&mode==146&&divide==-9999)
{printf("Error, you can only use \"--powerfile\" with \"--divisions\"\n\n");exit(1);}

////////

if(plet!=-9999&&mode!=147)
{printf("Error, you can only use \"--absolutes\" with \"--sum-cors\"\n\n");exit(1);}

////////

if(strcmp(matfile,"blank")!=0&&mode!=144&&mode!=146)
{printf("Error, you can only use \"--matrix\" with \"--reduce-tagging\" or \"--sum-hers\"\n\n");exit(1);}

if(strcmp(matfile,"blank")!=0&&strcmp(wayfile,"blank")!=0)
{printf("Error, you can not use \"--matrix\" with \"--pathways\"\n\n");exit(1);}

if(strcmp(matfile,"blank")!=0&&divide!=-9999)
{printf("Error, you can not use \"--matrix\" with \"--divisions\"\n\n");exit(1);}

////////

if(strcmp(expfile,"blank")!=0&&mode!=150)
{printf("Error, you can only use \"--expectations\" with \"--calc-posts\"\n\n");exit(1);}

if(cvar!=-9999&&mode!=150)
{printf("Error, you can only use \"--causal-variance\" with \"--calc-posts\"\n\n");exit(1);}

//////////////////////////

//individual-level data prediction, then megaprs

if(strcmp(indhers,"blank")&&gigaprs==1)
{printf("Error, you can not use \"--ind-hers\" with \"--giga-prs\"\n\n");}

if(strcmp(indhers,"blank")!=0&&mode!=151&&mode!=152&&mode!=153&&mode!=154&&mode!=159)
{printf("Error, you can only use \"--ind-hers\" with \"--ridge\", \"--bolt\", \"--bayesr\", \"--elastic\" or \"--mega-prs\"\n\n");exit(1);}

if(strcmp(indhers,"blank")!=0&&(strcmp(weightsfile,"blank")!=0||ignoreweights!=-9999||power!=-9999||her!=-9999))
{printf("Error, when using \"--ind-hers\" you should not use \"--weights\", \"--ignore-weights\", \"--power\" or \"--her\" (\"--ind-hers %s\" is equivalent to using \"--weights %s\" with \"--power -1\" and her equal to the sum of the per-predictor heritabilities)\n\n", indhers2, indhers2);exit(1);}

if(strcmp(indhers,"blank")!=0&&mpheno==-1)
{printf("Error, you can not use \"--ind-hers\" with \"--mpheno ALL\"\n\n");exit(1);}

if(herscale!=-9999&&strcmp(indhers,"blank")==0)
{printf("Error, you can only use \"--her-scale\" with \"--ind-hers\"\n\n");exit(1);}

////////

if(loco!=-9999&&mode!=151&&mode!=152&&mode!=153&&mode!=154&&mode!=172)
{printf("Error, you can only use \"--LOCO\" with \"--ridge\", \"--bolt\", \"--bayesr\", \"--elastic\" or \"--calc-scores\"\n\n");exit(1);}

if(verbose!=-9999&&mode!=151&&mode!=152&&mode!=153&&mode!=154)
{printf("Error, you can only use \"--verbose\" with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n");exit(1);}

if(dichot!=-9999&&mode!=151&&mode!=152&&mode!=153&&mode!=154&&mode!=176)
{printf("Error, you can only use \"--binary\" with \"--ridge\", \"--bolt\", \"--bayesr\", \"--elastic\" or \"--jackknife\"\n\n");exit(1);}

if(dichot==1&&strcmp(sampwfile,"blank")!=0)
{printf("Error, you can not use \"--binary YES\" with \"--sample-weights\"\n\n");exit(1);}

if(dichot==1&&loco==0)
{printf("Error, you can not use \"--binary YES\" with \"--LOCO NO\"\n\n");exit(1);}

if(multi!=-9999&&mode!=131&&mode!=134&&mode!=151&&mode!=152&&mode!=153&&mode!=154)
{printf("Error, you can only use \"--multivariate\" with \"--linear\", \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n");exit(1);}

if(multi==1&&mpheno!=-1)
{printf("Error, you can only use \"--multivariate YES\" with \"--mpheno ALL\"\n\n");exit(1);}

if(multi==1&&mode==131&&strcmp(prsfile,"blank")==0)
{printf("Error, you can only use \"--multivariate YES\" with \"--PRS\"\n\n");exit(1);}

if(multi==1&&strcmp(sampwfile,"blank")!=0)
{printf("Error, you can not use \"--multivariate YES\" with \"--sample-weights\"\n\n");exit(1);}

if(multi==1&&dichot==1)
{printf("Error, you can not use \"--multivariate YES\" with \"--binary YES\"\n\n");exit(1);}

if(multi==1&&gctastep==1)
{printf("Error, you can not use \"--multivariate YES\" with \"--GCTA-LOCO-step1\"\n\n");exit(1);}

if(multi==1&&faststep==1)
{printf("Error, you can not use \"--multivariate YES\" with \"--fastGWA-step1\"\n\n");exit(1);}

if(fprs!=-9999&&mode!=151&&mode!=152&&mode!=153&&mode!=154)
{printf("Error, you can only use \"--force-prs\" with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n");exit(1);}

if(fprs==1&&loco!=1)
{printf("Error, you can only use \"--force-prs\" with \"--LOCO YES\"\n\n");exit(1);}

if(fastgwa!=-9999&&mode!=151)
{printf("Error, you can only use \"--fastGWA\" with \"--ridge\"\n\n");exit(1);}

if(fastgwa==1&&multi==1)
{printf("Error, you can not use \"--fastGWA YES\" with \"--multivariate YES\"\n\n");exit(1);}

if(fastgwa==1&&loco==0)
{printf("Error, you can only use \"--fastGWA YES\" with \"--LOCO YES\"\n\n");exit(1);}

////////

if(skipcv!=-9999&&mode!=152&&mode!=153&&mode!=154&&mode!=159)
{printf("Error, you can only use \"--skip-cv\" with \"--bolt\", \"--bayesr\", \"--elastic\", \"--mega-prs\" or \"--giga-prs\"\n\n");exit(1);}

if(cvprop!=-9999&&mode!=151&&mode!=152&&mode!=153&&mode!=154&&mode!=159&&mode!=176)
{printf("Error, you can only use \"--cv-proportion\" with \"--ridge\", \"--bolt\" \"--bayesr\", \"--elastic\", \"--mega-prs\", \"--giga-prs\" or \"--jackknife\"\n\n");exit(1);}

if(strcmp(bvsfile,"blank")!=0&&mode!=151&&mode!=152&&mode!=153&&mode!=154&&mode!=176)
{printf("Error, you can only use \"--cv-samples\" with \"--ridge\", \"--bolt\", \"--bayesr\", \"--elastic\" or \"--jackknife\"\n\n");exit(1);}

if(cvprop!=-9999&&skipcv==1)
{printf("Error, it is contradictory to use both \"--cv-proportion\" and \"--skip-cv YES\"\n\n");exit(1);}

if(strcmp(bvsfile,"blank")!=0&&skipcv==1)
{printf("Error, it is contradictory to use both \"--cv-samples\" and \"--skip-cv YES\"\n\n");exit(1);}

if(cvprop!=-9999&&strcmp(bvsfile,"blank")!=0)
{printf("Error, you can not use both \"--cv-proportion\" and \"--cv-samples\"\n\n");exit(1);}

if(skipcv==1&&loco==1)
{printf("Error, you can not use \"--skip-cv YES\" with \"--LOCO YES\"\n\n");exit(1);}

if(skipcv==1&&multi==1)
{printf("Error, you can not use \"--skip-cv YES\" with \"--multivariate YES\"\n\n");exit(1);}

////////

if(ldpred!=-9999&&mode!=152)
{printf("Error, you can only use \"--LDpred\" with \"--bolt\"\n\n");exit(1);}

if(pointmass!=-9999&&mode!=153)
{printf("Error, you can only use \"--point-mass\" with \"--bayesr\"\n\n");exit(1);}

////////

if(ndivs!=-9999&&mode!=151&&mode!=152&&mode!=153&&mode!=154)
{printf("Error, you can only use \"--num-divides\" with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n");exit(1);}

if(nmcmc!=-9999&&drm==-9999&&mode!=151&&mode!=152&&mode!=153&&mode!=154)
{printf("Error, you can only use \"--num-random-vectors\" with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\", or with \"--linear\" and \"--DRM\"\n\n");exit(1);}

if(maxher!=-9999&&mode!=151&&mode!=152&&mode!=153&&mode!=154&&mode!=159)
{printf("Error, you can only use \"--max-her\" with \"--ridge\", \"--bolt\", \"--bayesr\", \"--elastic\", \"--mega-prs\" or \"--giga-prs\"\n\n");exit(1);}

////////

if(checkped!=-9999&&mode!=151&&mode!=152&&mode!=153&&mode!=154)
{printf("Error, you can only use \"--check-pedigree\" with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n");exit(1);}

if(checkped==1&&loco==0)
{printf("Error, you can only use \"--check-pedigree YES\" with \"--LOCO YES\"\n\n");exit(1);}

if(checkped==0&&fastgwa==1)
{printf("Error, you can not use \"--check-pedigree NO\" with \"--fastGWA YES\"\n\n");exit(1);}

if(nped!=-9999&&mode!=151&&mode!=152&&mode!=153&&mode!=154)
{printf("Error, you can only use \"--num-pedigree-predictors\" with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n");exit(1);}

if(nped!=-9999&&loco==0)
{printf("Error, you can only use \"--num-pedigree-predictors\" with \"--LOCO YES\"\n\n");exit(1);}

if(nped!=-9999&&checkped==0)
{printf("Error, it is contradictory to use \"--num-pedigree-predictors\" with \"--check-pedigree NO\"\n\n");exit(1);}

if(maithresh!=-9999&&mode!=151&&mode!=152&&mode!=153&&mode!=154)
{printf("Error, you can only use \"--MAI-threshold\" with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n");exit(1);}

if(maithresh!=-9999&&loco==0)
{printf("Error, you can only use \"--MAI-threshold\" with \"--LOCO YES\"\n\n");exit(1);}

if(maithresh!=-9999&&checkped==0)
{printf("Error, it is contradictory to use \"--MAI-threshold\" with \"--check-pedigree NO\"\n\n");exit(1);}

if(maisig!=-9999&&mode!=151&&mode!=152&&mode!=153&&mode!=154)
{printf("Error, you can only use \"--MAI-significance\" with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n");exit(1);}

if(maisig!=-9999&&loco==0)
{printf("Error, you can only use \"--MAI-significance\" with \"--LOCO YES\"\n\n");exit(1);}

if(maisig!=-9999&&checkped==0)
{printf("Error, it is contradictory to use \"--MAI-significance\" with \"--check-pedigree NO\"\n\n");exit(1);}

if(ncal!=-9999&&mode!=151&&mode!=152&&mode!=153&&mode!=154)
{printf("Error, you can only use \"--num-calibration-predictors\" with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n");exit(1);}

if(ncal!=-9999&&loco==0)
{printf("Error, you can only use \"--num-calibration-predictors\" with \"--LOCO YES\"\n\n");exit(1);}

if(ncomp!=-9999&&mode!=152&&mode!=153&&mode!=154)
{printf("Error, you can only use \"--num-comparison-predictors\" with \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n");exit(1);}

if(cthresh!=-9999&&mode!=152&&mode!=153&&mode!=154)
{printf("Error, you can only use \"--comparison-threshold\" with \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n");exit(1);}

if(nscan!=-9999&&mode!=151&&mode!=152&&mode!=153&&mode!=154)
{printf("Error, you can only use \"--num-scans\" with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n");exit(1);}

if(revher!=-9999&&mode!=151&&mode!=152&&mode!=153&&mode!=154)
{printf("Error, you can only use \"--revise-heritability\" with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n");exit(1);}

////////

if(substep!=-9999&&mode!=151&&mode!=152&&mode!=153&&mode!=154)
{printf("Error, you can only use \"--sub-step\" with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n");exit(1);}

if(substep!=-9999&&multi!=1)
{printf("Error, you can only use \"--sub-step\" with \"--multivariate YES\"\n\n");exit(1);}

////////

if(strcmp(blockfile,"blank")!=0&&mode!=156)
{printf("Error, you can only use \"--break-points\" with \"--calc-cors\"\n\n");exit(1);}

if(strcmp(blockfile,"blank")!=0&&(window_kb!=-9999||window_cm!=-9999))
{printf("Error, you can not use \"--break-points\" with \"--window-kb\" or \"--window-cm\"\n\n");exit(1);}

if(strcmp(betafile,"blank")!=0&&mode!=156)
{printf("Error, you can only use \"--betas\" with \"--calc-cors\"\n\n");exit(1);}

////////

if(strcmp(corslist,"blank")!=0&&mode!=157&&mode!=159)
{printf("Error, you can only use \"--corslist\" with \"--join-cors\", \"--mega-prs\" or \"--giga-prs\"\n\n");exit(1);}

if(strcmp(corslist,"blank")!=0&&strcmp(sumsfile,"blank")!=0)
{printf("Error, when using \"--corslist\", you should use \"--sumslist\" (instead of \"--summary\")\n\n");exit(1);}

////////

if(subprop!=-9999&&mode!=158)
{printf("Error, you can only use \"--training-proportion\" with \"--pseudo-summaries\"\n\n");exit(1);}

////////

if(strcmp(corname,"blank")!=0&&mode!=159&&mode!=180)
{printf("Error, you can only use \"--cors\" with \"--mega-prs\", \"--giga-prs\" or \"--impute-summaries\"\n\n");exit(1);}

if(strcmp(corname,"blank")!=0&&strcmp(corslist,"blank")!=0)
{printf("Error, you can not use both \"--cors\" and \"--corslist\"\n\n");exit(1);}

if(strcmp(pseudostem,"blank")!=0&&mode!=159)
{printf("Error, you can only use \"--pseudos\" with \"--mega-prs\" or \"--giga-prs\"\n\n");exit(1);}

if(strcmp(pseudostem,"blank")!=0&&cvprop!=-9999)
{printf("Error, you can not use both \"--pseudos\" and \"--cv-proportion\"\n\n");exit(1);}

if(strcmp(pseudostem,"blank")!=0&&skipcv==1)
{printf("Error, it is contradictory to use both \"--pseudos\" and \"--skip-cv YES\"\n\n");exit(1);}

if(strcmp(pseudostem,"blank")!=0&&strcmp(sumslist,"blank")!=0)
{printf("Error, you can not use \"--pseudos\" with \"--sumslist\"\n\n");exit(1);}

////////

if(useann!=-9999&&gigaprs!=1)
{printf("Error, you can only use \"--use-annotations\" with \"--giga-prs\"\n\n");exit(1);}

if(condition!=-9999&&mode!=159)
{printf("Error, you can only use \"--conditional\" with \"--mega-prs\" or \"--giga-prs\"\n\n");exit(1);}

if(finemap!=-9999&&mode!=159)
{printf("Error, you can only use \"--fine-map\" with \"--mega-prs\" or \"--giga-prs\"\n\n");exit(1);}

if(finemap==1&&skipcv==0)
{printf("Error, you can not use \"--fine-map YES\" with \"--skip-cv NO\"\n\n");exit(1);}

if(Lnum!=-9999&&mode!=159)
{printf("Error, you can only use \"--num-effects\" with \"--mega-prs\" or \"--giga-prs\"\n\n");exit(1);}

if(Lnum!=-9999&&finemap!=1)
{printf("Error, you can only use \"--num-effects\" with \"--fine-map YES\"\n\n");exit(1);}

if(num_grids!=-9999&&mode!=159)
{printf("Error, you can only use \"--num-variances\" with \"--mega-prs\" or \"--giga-prs\"\n\n");exit(1);}

if(num_grids!=-9999&&finemap!=1)
{printf("Error, you can only use \"--num-variances\" with \"--fine-map YES\"\n\n");exit(1);}

if(grid_min!=-9999&&mode!=159)
{printf("Error, you can only use \"--min-variance\" with \"--mega-prs\" or \"--giga-prs\"\n\n");exit(1);}

if(grid_min!=-9999&&finemap!=1)
{printf("Error, you can only use \"--min-variance\" with \"--fine-map YES\"\n\n");exit(1);}

if(grid_max!=-9999&&mode!=159)
{printf("Error, you can only use \"--max-variance\" with \"--mega-prs\" or \"--giga-prs\"\n\n");exit(1);}

if(grid_max!=-9999&&finemap!=1)
{printf("Error, you can only use \"--max-variance\" with \"--fine-map YES\"\n\n");exit(1);}

if(grid_min!=-9999&&grid_max==-9999)
{printf("Error, when using \"--min-variance\" you must also use \"--max-variance\"\n\n");exit(1);}

if(grid_min==-9999&&grid_max!=-9999)
{printf("Error, when using \"--max-variance\" you must also use \"--min-variance\"\n\n");exit(1);}

if(grid_min!=-9999&&grid_min>=grid_max)
{printf("Error, \"--min-variance\" must be less than \"--max-variance\"\n\n");exit(1);}

if(csalpha!=-9999&&mode!=159)
{printf("Error, you can only use \"--CS-coverage\" with \"--mega-prs\" or \"--giga-prs\"\n\n");exit(1);}

if(csalpha!=-9999&&finemap!=1)
{printf("Error, you can only use \"--CS-coverage\" with \"--fine-map YES\"\n\n");exit(1);}

if(impsums!=-9999&&mode!=159)
{printf("Error, you can only use \"--summary-imputation\" with  \"--mega-prs\" or \"--giga-prs\"\n\n");exit(1);}

////////

if(ptype!=-9999&&mode!=159)
{printf("Error, you can only use \"--model\" with \"--mega-prs\" or \"--giga-prs\"\n\n");exit(1);}

if(ptype!=-9999&&finemap==1)
{printf("Error, you can not use \"--model\" with \"--fine-map YES\"\n\n");exit(1);}

if(strcmp(bestfile,"blank")!=0&&mode!=159)
{printf("Error, you can only use \"--best-model\" with \"--mega-prs\" or \"--giga-prs\"\n\n");exit(1);}

if(strcmp(bestfile,"blank")!=0&&ptype!=-9999)
{printf("Error, you can not use \"--best-model\" with \"--model\"\n\n");exit(1);}

if(strcmp(bestfile,"blank")!=0&&strcmp(indhers,"blank")==0)
{printf("Error, when using \"--best-model\", you must also use \"--ind-hers\"\n\n");exit(1);}

if(strcmp(bestfile,"blank")!=0&&cvprop!=-9999)
{printf("Error, you can not use \"--best-model\" with \"--cv-proportion\"\n\n");exit(1);}

if(strcmp(bestfile,"blank")!=0&&skipcv==0)
{printf("Error, you can not use \"--best-model\" with \"--skip-cv NO\"\n\n");exit(1);}

if(strcmp(bestfile,"blank")!=0&&strcmp(pseudostem,"blank")!=0)
{printf("Error, you can not use \"--best-model\" with \"--pseudos\"\n\n");exit(1);}

if(exparam!=-9999&&mode!=159)
{printf("Error, you can only use \"--extra-parameters\" with \"--mega-prs\" or \"--giga-prs\"\n\n");exit(1);}

if(exparam==1&&ptype!=-9999&&ptype!=5&&ptype!=6)
{printf("Error, you can only use \"--extra-parameters YES\" with \"--model bayesr\" or \"--model bayesr-shrink\"\n\n");exit(1);}

if(ensemble!=-9999&&mode!=159)
{printf("Error, you can only use \"--ensemble\" with \"--mega-prs\" or \"--giga-prs\"\n\n");exit(1);}

if(ensemble==1&&skipcv==1)
{printf("Error, you can not use \"--ensemble YES\" with \"--skip-cv YES\"\n\n");exit(1);}

if(ensemble==1&&finemap==1)
{printf("Error, you can not use \"--ensemble YES\" with \"--fine-map YES\"\n\n");exit(1);}

if(megasave!=-9999&&mode!=159)
{printf("Error, you can only use \"--save-all-models\" with \"--mega-prs\" or \"--giga-prs\"\n\n");exit(1);}

if(megasave==1&&skipcv==1)
{printf("Error, you can not use \"--save-all-models YES\" with \"--skip-cv YES\"\n\n");exit(1);}

if(megasave==1&&finemap==1)
{printf("Error, you can not use \"--save-all-models YES\" with \"--fine-map YES\"\n\n");exit(1);}

if(strcmp(fracfile,"blank")!=0&&mode!=152&&mode!=153&&mode!=154&&mode!=159)
{printf("Error, you can only use \"--parameters\" with \"--bolt\", \"--bayesr\", \"--elastic\", \"--mega-prs\" or \"--giga-prs\"\n\n");exit(1);}

if(strcmp(fracfile,"blank")!=0&&(mode==151||mode==152||mode==153||mode==154)&&loco==1)
{printf("Error, you can not use \"--parameters\" if using \"--LOCO YES\"\n\n");exit(1);}

if(strcmp(fracfile,"blank")!=0&&finemap==1)
{printf("Error, you can not use \"--parameters\" with \"--fine-map YES\"\n\n");exit(1);}

if(strcmp(fracfile,"blank")!=0&&(ptype==0||ptype==1||ptype==2||ptype==3))
{printf("Error, you can not use \"--parameters\" with \"--model lasso-sparse\", \"--model lasso\", \"--model ridge\" or \"--model mega\"\n\n");exit(1);}

if(strcmp(fracfile,"blank")!=0&&strcmp(bestfile,"blank")!=0)
{printf("Error, you can not use \"--parameters\" with \"--best-model\"\n\n");exit(1);}

////////

if(checkld!=-9999&&mode!=159)
{printf("Error, you can only use \"--check-high-LD\" with \"--mega-prs\" or \"--giga-prs\"\n\n");exit(1);}

if(strcmp(ldfile,"blank")!=0&&mode!=159)
{printf("Error, you can only use \"--high-LD\" with \"--mega-prs\" or \"--giga-prs\"\n\n");exit(1);}

if(checkld==0&&strcmp(ldfile,"blank")!=0)
{printf("Error, it is contradictory to use both \"--high-LD\" and \"--check-high-LD NO\"\n\n");exit(1);}

if(checkld!=-9999&&skipcv==1)
{printf("Error, you can not use \"--check-high-LD\" with \"--skip-cv YES\"\n\n");exit(1);}

if(strcmp(ldfile,"blank")!=0&&skipcv==1)
{printf("Error, you can not use \"--high-LD\" with \"--skip-cv YES\"\n\n");exit(1);}

////////

if(checkpops!=-9999&&mode!=159)
{printf("Error, you can only use \"--check-assignments\" with \"--mega-prs\" or \"--giga-prs\"\n\n");exit(1);}

if(checkfreq!=-9999&&mode!=159&&mode!=180)
{printf("Error, you can only use \"--check-frequencies\" with \"--mega-prs\", \"--giga-prs\" or \"--impute-summaries\"\n\n");exit(1);}

if(maxfreq!=-9999&&mode!=159)
{printf("Error, you can only use \"--max-frequency-difference\" with \"--mega-prs\" or \"--giga-prs\"\n\n");exit(1);}

if(resfix!=-9999&&mode!=159)
{printf("Error, you can only use \"--fix-residual-variance\" with \"--mega-prs\" or \"--giga-prs\"\n\n");exit(1);}

if(minher!=-9999&&mode!=159)
{printf("Error, you can only use \"--min-her\" with \"--mega-prs\" or \"--giga-prs\"\n\n");exit(1);}

if(mcmcfinal!=-9999&&mode!=159)
{printf("Error, you can only use \"--MCMC-solve\" with \"--mega-prs\" or \"--giga-prs\"\n\n");exit(1);}

if(mcmcfinal==1&&finemap==1)
{printf("Error, you can not use \"--MCMC-solve\" with \"--fine-map YES\"\n\n");exit(1);}

if(num_chains!=-9999&&mode!=159)
{printf("Error, you can only use \"--number-MCMC-chains\" with \"--mega-prs\" or \"--giga-prs\"\n\n");exit(1);}

if(num_chains!=-9999&&mcmcfinal==0)
{printf("Error, it is contradictory to use \"--number-MCMC-chains\" with \"--MCMC-solve NO\"\n\n");exit(1);}

if(num_burn!=-9999&&mode!=159)
{printf("Error, you can only use \"--number-MCMC-burn\" with \"--mega-prs\" or \"--giga-prs\"\n\n");exit(1);}

if(num_burn!=-9999&&mcmcfinal==0)
{printf("Error, it is contradictory to use \"--number-MCMC-burn\" with \"--MCMC-solve NO\"\n\n");exit(1);}

if(num_mcmc!=-9999&&mode!=159)
{printf("Error, you can only use \"--number-MCMC-iterations\" with \"--mega-prs\" or \"--giga-prs\"\n\n");exit(1);}

if(num_mcmc!=-9999&&mcmcfinal==0)
{printf("Error, it is contradictory to use \"--number-MCMC-iterations\" with \"--MCMC-solve NO\"\n\n");exit(1);}

if(prsvar!=-9999&&mode!=159&&mode!=172)
{printf("Error, you can only use \"--PRS-variance\" with \"--mega-prs\",  \"--giga-prs\" or \"--calc-scores\"\n\n");exit(1);}

if(prsvar==1&&skipcv==1)
{printf("Error, you can not use \"--PRS-variance YES\" with \"--skip-cv YES\"\n\n");exit(1);}

if(prsvar==1&&finemap==1)
{printf("Error, you can not use \"--PRS-variance YES\" with \"--fine-map YES\"\n\n");exit(1);}

if(prsvar==1&&mode==159&&mcmcfinal!=1)
{printf("Error, you can only use \"--PRS-variance YES\" with \"--MCMC-solve YES\"\n\n");exit(1);}

if(prsvar==1&&ensemble==1)
{printf("Error, you can not use \"--PRS-variance YES\" with \"--ensemble YES\"\n\n");exit(1);}

if(prsvar==1&&megasave==1)
{printf("Error, you can not use \"--PRS-variance YES\" with \"--save-all-models YES\"\n\n");exit(1);}

if(prsvar==1&&loco==1)
{printf("Error, you can not use \"--PRS-variance YES\" with \"--LOCO YES\"\n\n");exit(1);}

if(num_step!=-9999&&mode!=159)
{printf("Error, you can only use \"--MCMC-step\" with \"--mega-prs\" or \"--giga-prs\"\n\n");exit(1);}

if(num_step!=-9999&&prsvar!=1)
{printf("Error, you can only use \"--MCMC-step\" with \"--PRS-variance YES\"\n\n");exit(1);}

if(mtagcor!=-9999&&mode!=159)
{printf("Error, you can only use \"--MTAG-correlation\" with \"--mega-prs\" or  \"--giga-prs\"\n\n");exit(1);}

if(mtagpercent!=-9999&&mode!=159)
{printf("Error, you can only use \"--MTAG-percent\" with \"--mega-prs\" or  \"--giga-prs\"\n\n");exit(1);}

if(mtagmeta!=-9999&&mode!=159)
{printf("Error, you can only use \"--MTAG-meta\" with \"--mega-prs\" or  \"--giga-prs\"\n\n");exit(1);}

if(mtagreg!=-9999&&mode!=159)
{printf("Error, you can only use \"--MTAG-regional\" with \"--mega-prs\" or  \"--giga-prs\"\n\n");exit(1);}

if(mtagforce!=-9999&&mode!=159)
{printf("Error, you can only use \"--MTAG-force-all\" with \"--mega-prs\" or  \"--giga-prs\"\n\n");exit(1);}

if(mtagforce==1&&mtagcor!=-9999)
{printf("Error, you can not use \"--MTAG-force-all YES\" with \"--MTAG-correlation\"\n\n");exit(1);}

if(mtagforce==1&&mtagpercent!=-9999)
{printf("Error, you can not use \"--MTAG-force-all YES\" with \"--MTAG-percent\"\n\n");exit(1);}

if(mtagcopy!=-9999&&mode!=159)
{printf("Error, you can only use \"--MTAG-copy\" with \"--mega-prs\" or  \"--giga-prs\"\n\n");exit(1);}

if(mtagcopy==1&&strcmp(sumslist,"blank")==0)
{printf("Error, you can only use \"--MTAG-copy YES\" with \"--sumslist\"\n\n");exit(1);}

if(mtagcopy==1&&mtagmeta==1)
{printf("Error, you can not use \"--MTAG-copy YES\" with \"--MTAG-meta YES\"\n\n");exit(1);}

if(mtagcopy==1&&mtagreg==1)
{printf("Error, you can not use \"--MTAG-copy YES\" with \"--MTAG-regional YES\"\n\n");exit(1);}

///////////////////////////

//pca, decompose and adjust-grm

if(axes!=-9999&&mode!=161&&mode!=167)
{printf("Error, you can only use \"--axes\" with \"--pca\"\n\n");exit(1);}

if(strcmp(pcastem,"blank")!=0&&mode!=162)
{printf("Error, you can only use \"--pcastem\" with \"--calc-pca-loads\"\n\n");exit(1);}

////////

if(eigenraw!=-9999&&mode!=163)
{printf("Error, you can only use \"--eigen-raw\" with \"--decompose\"\n\n");exit(1);}

if(strcmp(eigenfile,"blank")!=0&&mode!=121&&mode!=131&&mode!=138&&mode!=161)
{printf("Error, you can only use \"--eigen\" with \"--reml\", \"--linear\", \"--calc-genes-reml\" or \"--pca\"\n\n");exit(1);}

if(strcmp(eigenfile,"blank")!=0&&strcmp(kinname,"blank")==0&&strcmp(kinlist,"blank")==0)
{printf("Error, you must use \"--grm\" to provide the kinship matrix to which the eigen-decomposition corresponds\n\n");exit(1);}

////////

if(noise!=-9999&&mode!=169)
{printf("Error, you can only use \"--noise\" with \"--gxemm-iid\"\n\n");exit(1);}

///////////////////////////

//stats, scores, making phenotypes, snps, jackknifing, folds, find gaussian

if(strcmp(scorefile,"blank")!=0&&mode!=108&&mode!=160&&mode!=172)
{printf("Error, you can only use \"--scorefile\" with \"--validate\", \"--find-tags\" or \"--calc-scores\"\n\n");exit(1);}

if(strcmp(morefile,"blank")!=0&&mode!=172)
{printf("Error, you can only use \"--multi-scorefile\" with \"--calc-scores\"\n\n");exit(1);}

if(strcmp(scorefile,"blank")!=0&&strcmp(morefile,"blank")!=0)
{printf("Error, you can not use both \"--scorefile\" and \"--multi-scorefile\"\n\n");exit(1);}

if(strcmp(morefile,"blank")!=0&&strcmp(sumsfile,"blank")!=0)
{printf("Error, you can only use \"--summary\" with \"--scorefile\" (not with \"--multi-scorefile\")\n\n");exit(1);}

if(strcmp(morefile,"blank")!=0&&prsvar==1)
{printf("Error, you can only use \"--PRS-variance YES\" with \"--scorefile\" (not with \"--multi-scorefile\")\n\n");exit(1);}

if(strcmp(cofile,"blank")!=0&&mode==122)
{printf("Error, you can not use \"--coeffsfile\" with \"--calc-blups\" (covariates details are linked to in the remlfile)\n\n");exit(1);}

if(strcmp(cofile,"blank")!=0&&mode!=172)
{printf("Error, you can only use \"--coeffsfile\" with \"--calc-scores\"\n\n");exit(1);}

if(savecounts!=-9999&&mode!=172)
{printf("Error, you can only use \"--save-counts\" with \"--calc-scores\"\n\n");exit(1);}

if(savecounts==1&&prsvar==1)
{printf("Error, you can not use \"--save-counts YES\" with \"--PRS-variance YES\"\n\n");exit(1);}

if(savecounts==1&&strcmp(morefile,"blank")!=0)
{printf("Error, you can not use \"--save-counts YES\" with \"--multi-scorefile\"\n\n");exit(1);}

if(strcmp(finalfile,"blank")!=0&&mode!=172)
{printf("Error, you can only use \"--final-effects\" with \"--calc-scores\"\n\n");exit(1);}

////////

if((num_phenos!=-9999||num_causals!=-9999)&&mode!=173)
{printf("Error, you can only use \"--num-phenos\" or \"--num-causals\" with \"--make-phenos\"\n\n");exit(1);}

if(her>1&&mode!=159)
{printf("Error, you can only use \"--her-big\" with \"--mega-prs\" or \"--giga-prs\"\n\n");exit(1);}

if(her!=-9999&&mode!=151&&mode!=152&&mode!=153&&mode!=154&&mode!=159&&mode!=173)
{printf("Error, you can only use \"--her\" with \"--ridge\", \"--bolt\", \"--bayesr\", \"--elastic\", \"--mega-prs\" or \"--make-phenos\"\n\n");exit(1);}

if(her!=-9999&&(mode==151||mode==152||mode==153||mode==154)&&power==-9999)
{printf("Error, you can only use \"--her\" if also using \"--power\"\n\n");exit(1);}

if(her!=-9999&&gigaprs==1)
{printf("Error, you can not use \"--her\" with \"--giga-prs\"\n\n");exit(1);}

if(her==0&&mode!=173)
{printf("Error, \"--her\" should be followed by a positive float (not 0)\n\n");exit(1);}

if(her!=-9999&&strcmp(partpref,"blank")!=0)
{printf("Error, you can not use \"--her\" when using partitions\n\n");exit(1);}

if(her!=-9999&&strcmp(annpref,"blank")!=0)
{printf("Error, you can not use \"--her\" when using annotations\n\n");exit(1);}

if(her!=-9999&&strcmp(taufile,"blank")!=0)
{printf("Error, you can not use both \"--her\" and \"--taus\"\n\n");exit(1);}

if(her!=-9999&&multi==1)
{printf("Error, you can not use \"--her\" with \"--multivariate YES\"\n\n");exit(1);}

if(her!=-9999&&revher==1)
{printf("Error, you can not use \"--her\" with \"--revise-heritability YES\"\n\n");exit(1);}

if(cher!=-9999&&mode!=173)
{printf("Error, you can only use \"--covar-her\" with \"--make-phenos\"\n\n");exit(1);}

if(cher!=-9999&&strcmp(covarfile,"blank")==0)
{printf("Error, you can only use \"--covar-her\" with \"--covar\"\n\n");exit(1);}

if(gxe!=-9999&&mode!=173)
{printf("Error, you can only use \"--GxE\" with \"--make-phenos\"\n\n");exit(1);}

if(gxe!=-9999&&strcmp(envfile,"blank")==0)
{printf("Error, when using \"--GxE\", you must also use \"--enviro\"\n\n");exit(1);}

if(oddhalf!=-9999&&mode!=173)
{printf("Error, you can only use \"--odd-half-her\" with \"--make-phenos\"\n\n");exit(1);}

if(bivar!=-9999&&mode!=173)
{printf("Error, you can only use \"--bivar\" with \"--make-phenos\"\n\n");exit(1);}

if(bivar2!=-9999&&mode!=173)
{printf("Error, you can only use \"--bivar-env\" with \"--make-phenos\"\n\n");exit(1);}

if(bivar2!=-9999&&strcmp(sampwfile,"blank")!=0)
{printf("Error, you can not use \"--bivar-env\" with \"--sample-weights\"\n\n");exit(1);}

if(bivar3!=-9999&&mode!=173)
{printf("Error, you can only use \"--bivar-proportion\" with \"--make-phenos\"\n\n");exit(1);}

if(bivar3!=-9999&&bivar!=-9999)
{printf("Error, you can not use \"--bivar-proportion\" and \"--bivar\"\n\n");exit(1);}

if(strcmp(probsfile,"blank")!=0&&mode!=173)
{printf("Error, you can only use \"--probabilties\" with \"--make-phenos\"\n\n");exit(1);}

if(strcmp(probsfile,"blank")!=0&&(strcmp(weightsfile,"blank")!=0||ignoreweights!=-9999||power!=-9999))
{printf("Error, when using \"--probabilties\" you should not use \"--weights\", \"--ignore-weights\" or \"--power\"\n\n");exit(1);}

if(strcmp(probsfile,"blank")!=0&&num_causals==-1)
{printf("Error, you can not use \"--probabilties\" with \"--num-causals -1\"\n\n");exit(1);}

if((strcmp(causalsfile,"blank")!=0||strcmp(effectsfile,"blank")!=0)&&mode!=173)
{printf("Error, you can only use \"--causals\" or \"--effects\" with \"--make-phenos\"\n\n");exit(1);}

if(strcmp(causalsfile,"blank")!=0&&extract==1)
{printf("Error, when using \"--causals\" you can not use \"--extract\", \"--exclude\", \"--chr\" or \"--snp\"\n\n");exit(1);}

if(strcmp(causalsfile,"blank")!=0&&strcmp(probsfile,"blank")!=0)
{printf("Error, you can not use \"--causals\" with \"--probabilities\"\n\n");exit(1);}

if(strcmp(causalsfile,"blank")!=0&&bivar!=-9999)
{printf("Error, you can not use \"--causals\" with \"--bivar\"\n\n");exit(1);}

if(strcmp(causalsfile,"blank")!=0&&bivar2!=-9999)
{printf("Error, you can not use \"--causals\" with \"--bivar-proportion\"\n\n");exit(1);}

if(strcmp(effectsfile,"blank")!=0&&bivar!=-9999)
{printf("Error, you can not use \"--effects\" with \"--bivar\"\n\n");exit(1);}

if(strcmp(effectsfile,"blank")!=0&&bivar2!=-9999)
{printf("Error, you can not use \"--effects\" with \"--bivar-proportion\"\n\n");exit(1);}

////////

if(num_inds!=-9999&&mode!=156&&mode!=174)
{printf("Error, you can only use \"--num-samples\" with \"--calc-cors\" or \"--make-snps\"\n\n");exit(1);}

if(num_snps!=-9999&&mode!=174)
{printf("Error, you can only use \"--num-snps\" with \"--make-snps\"\n\n");exit(1);}

if((maf1!=-9999||maf2!=-9999)&&mode!=174)
{printf("Error, you can only use \"--maf-low\" or \"--maf-high\" with \"--make-snps\"\n\n");exit(1);}

if(maf1!=-9999&&maf2!=-9999&&maf1>maf2)
{printf("Error, maf-low (%.6f) can not be higher than maf-high (%.6f)\n\n", maf1, maf2);exit(1);}

if(nchrom!=-9999&&mode!=174)
{printf("Error, you can only use \"--num-chr\" with \"--make-snps\"\n\n");exit(1);}

if(famsize!=-9999&&mode!=174)
{printf("Error, you can only use \"--family-size\" with \"--make-snps\"\n\n");exit(1);}

if(famsize!=-9999&&trios==1)
{printf("Error, you can not use \"--family-size\" with \"--trios YES\"\n\n");exit(1);}

if(quads!=-9999&&mode!=174)
{printf("Error, you can only use \"--quads\" with \"--make-snps\"\n\n");exit(1);}

if(quads==1&&trios==1)
{printf("Error, you can not use both \"--quads YES\" and \"--trios YES\"\n\n");exit(1);}

if(quads==1&&famsize!=-9999)
{printf("Error, you can not use \"--family-size\" with \"--quads YES\"\n\n");exit(1);}

if(closeness!=-9999&&mode!=174)
{printf("Error, you can only use \"--relatedness\" with \"--make-snps\"\n\n");exit(1);}

if(closeness!=-9999&&trios==1)
{printf("Error, you can not use \"--relatedness\" with \"--trios YES\\n\n");exit(1);}

if(closeness!=-9999&&quads==1)
{printf("Error, you can not use \"--relatedness\" with \"--quads YES\\n\n");exit(1);}

if(closeness!=-9999&&famsize==-9999)
{printf("Error, you can only use \"--relatedness\" with \"--family-size\"\n\n");exit(1);}

if(pops!=-9999&&mode!=174)
{printf("Error, you can only use \"--populations\" with \"--make-snps\"\n\n");exit(1);}

if(pops!=-9999&&trios==1)
{printf("Error, you can not use \"--populations\" with \"--trios YES\"\n\n");exit(1);}

if(pops!=-9999&&quads==1)
{printf("Error, you can not use \"--populations\" with \"--quads YES\"\n\n");exit(1);}

if(missrate!=-9999&&mode!=174)
{printf("Error, you can only use \"--missing-rate\" with \"--make-snps\"\n\n");exit(1);}

if(missrate!=-9999&&trios==1)
{printf("Error, you can not use \"--missing-rate\" with \"--trios YES\\n\n");exit(1);}

if(missrate!=-9999&&quads==1)
{printf("Error, you can not use \"--missing-rate\" with \"--quads YES\\n\n");exit(1);}

////////

if((strcmp(predlista,"blank")!=0||strcmp(predlistb,"blank")!=0)&&mode!=175)
{printf("Error, you can only use \"--lista\" or \"--listb\" with \"--calc-inflation\"\n\n");exit(1);}

if(savepairs!=-9999&&mode!=175)
{printf("Error, you can only use \"--save-pairs\" with \"--calc-inflation\"\n\n");exit(1);}

////////

if(strcmp(jackfile,"blank")!=0&&mode!=176)
{printf("Error, you can only use \"--data-pairs\" with \"--jackknife\"\n\n");exit(1);}

if(strcmp(proffile,"blank")!=0&&mode!=176)
{printf("Error, you can only use \"--profile\" with \"--jackknife\"\n\n");exit(1);}

if(strcmp(jackfile,"blank")!=0&&strcmp(proffile,"blank")!=0)
{printf("Error, you can not use both \"--data-pairs\" and \"--profile\"\n\n");exit(1);}

if(cvprop!=-9999&&strcmp(jackfile,"blank")!=0)
{printf("Error, you can only use \"--cv-proportion\" with \"--profile\" (not \"--data-pairs\")\n\n");exit(1);}

if(strcmp(bvsfile,"blank")!=0&&strcmp(jackfile,"blank")!=0)
{printf("Error, you can only use \"--cv-samples\" with \"--profile\" (not \"--data-pairs\")\n\n");exit(1);}

////////

if(num_folds!=-9999&&mode!=177)
{printf("Error, you can only use \"--num-folds\" with \"--cut-folds\"\n\n");exit(1);}

////////

if(strcmp(likefile,"blank")!=0&&mode!=178)
{printf("Error, you can only use \"--likelihoods\" with \"--find-gaussian\"\n\n");exit(1);}

if(num_means!=-9999&&mode!=178)
{printf("Error, you can only use \"--num-means\" with \"--find-gaussian\"\n\n");exit(1);}

if(num_sds!=-9999&&mode!=178)
{printf("Error, you can only use \"--num-sds\" with \"--find-gaussian\"\n\n");exit(1);}

if(minmean!=-9999||maxmean!=-9999)
{
if(mode!=178)
{printf("Error, you can only use \"--min-mean\" or \"--max-mean\" with \"--find-gaussian\"\n\n");exit(1);}

if(maxmean==-9999)
{printf("Error, when using \"--min-mean\" you must also use \"--max-mean\"\n\n");exit(1);}

if(minmean==-9999)
{printf("Error, when using \"--max-mean\" you must also use \"--min-mean\"\n\n");exit(1);}

if(maxmean<=minmean)
{printf("Error, max-mean (%.6f) must be higher than min-mean (%.6f)\n\n", maxmean, minmean);exit(1);}
}

if(maxsd!=-9999&&mode!=178)
{printf("Error, you can only use \"--max-sd\" with \"--find-gaussian\"\n\n");exit(1);}

if(omitone!=-9999&&mode!=178)
{printf("Error, you can only use \"--omit-one\" with \"--find-gaussian\"\n\n");exit(1);}

///////////////////////////

//making and condensing data

if(exsame==0&&mode!=131&&mode!=132&&mode!=134)
{printf("Error, you can only use \"--exclude-same-names NO\" with \"--linear\" or \"--logistic\"\n\n");exit(1);}

if(exlong==0&&mode!=106&&mode!=107&&mode!=131&&mode!=132&&mode!=134&&mode!=171&&mode!=181&&mode!=182&&mode!=183&&mode!=184)
{printf("Error, you can only use \"--exclude-long-alleles NO\" with \"--linear\", \"--logistic\", \"--thin\", \"--thin-tops\", \"--calc-stats\" or when making data\n\n");exit(1);}

if(exincons!=-9999&&mode!=181&&mode!=182&&mode!=183&&mode!=184)
{printf("Error, you can only use \"--exclude-inconsistent\" when making data\n\n");exit(1);}

if(speedlong!=-9999&&mode!=184&&mode!=189)
{printf("Error, you can only use \"--speed-long\" with \"--make-speed\" or \"--condense-speed\"\n\n");exit(1);}

////////

if(useminor!=-9999&&mode!=186&&mode!=187&&mode!=188&&mode!=189)
{printf("Error, you can only use \"--count-minor\" when condensing data\n\n");exit(1);}

if(useminor!=-9999&&nonsnp==1)
{printf("Error, you can not use \"--count-minor\" with \"--SNP-data NO\"\n\n");exit(1);}

if(resample!=-9999&&mode!=181&&mode!=182&&mode!=183&&mode!=184&&mode!=185)
{printf("Error, you can only use \"--resample\" when making data\n\n");exit(1);}

if(resample!=-9999&&strcmp(datalist,"blank")!=0)
{printf("Error, you can not use \"--resample\" with  \"--mbfile\", \"--msp\", \"--msped\" or \"--mspeed\"\n\n");exit(1);}

if(recomb!=-9999&&resample==-9999)
{printf("Error, you can only use \"--recombination-rate\" with \"--resample\"\n\n");exit(1);}

//////////////////////////

//gre options

if(sinv!=-9999&&mode!=193)
{printf("Error, you can only use \"--save-inverse\" with \"--join-gre\"\n\n");exit(1);}

if(strcmp(greout,"blank")!=0&&mode!=194)
{printf("Error, you can only use \"--gre-output\" with \"--solve-gre\"\n\n");exit(1);}

//////////////////////////

//common options

if(checkroot!=-9999&&mode!=122&&mode!=123&&mode!=124&&mode!=125&&mode!=162&&strcmp(eigenfile,"blank")==0&&strcmp(prsfile,"blank")==0)
{printf("Error, you can only use \"--check-root\" with \"--calc-blups\", \"--he\", \"--pcgc\", \"--reml-predict\" or \"--calc-pca-loads\", or when using \"--eigen\" or \"--PRS\"\n\n");exit(1);}

if(mincor!=-9999&&mode!=102&&mode!=104&&mode!=108&&mode!=109&&mode!=141&&mode!=151&&mode!=152&&mode!=153&&mode!=154&&mode!=156)
{printf("Error, you can only use \"--min-cor\" with \"--calc-weights\", \"--calc-weights-all\", \"--calc-tagging\", \"--ridge\", \"--bolt\", \"--bayesr\", \"--elastic\", \"--calc-cors\", \"--find-tags\" or \"--remove-tags\"\n\n");exit(1);}

if(maxcor!=-9999&&mode!=151&&mode!=152&&mode!=153&&mode!=154&&mode!=156&&mode!=193)
{printf("Error, you can only use \"--max-cor\" with \"--ridge\", \"--bolt\", \"--bayesr\", \"--elastic\", \"--calc-cors\" or \"--join-gre\"\n\n");exit(1);}

if(cutoff!=-9999&&mode!=107&&mode!=110&&mode!=146&&mode!=147&&mode!=159&&mode!=166&&mode!=179)
{printf("Error, you can only use \"--cutoff\" with \"--sum-hers\", \"--sum-cors\", \"--mega-prs\", \"--giga-prs\", \"--truncate-grm\", \"--thin-tops\" or \"--thin-common\"\n\n");exit(1);}

if(mode==110&&nonsnp==0&&cutoff>=0.5)
{printf("Error, \"--cutoff\" (%.4f) must be less than 0.5 (the default value is 0.01)\n\n", cutoff);exit(1);}

if(cutoff2!=-9999&&mode!=146&&mode!=147)
{printf("Error, you can only use \"--truncate\" with \"--sum-hers\" or \"--sum-cors\"\n\n");exit(1);}

if(cutoff!=-9999&&cutoff2!=-9999)
{printf("Error, you can not use both \"--cutoff\" and \"--truncate\"\n\n");exit(1);}

if(constrain!=-9999&&mode!=121&&mode!=126&&mode!=129&&mode!=130&&mode!=229&&mode!=230&&mode!=131&&mode!=133)
{printf("Error, you can only use \"--constrain\" with \"--reml\", \"--fast-reml\", \"--quant-her\", \"--quant-her\", \"--linear\" or \"--solve-null\"\n\n");exit(1);}

if(num_blocks!=-9999&&mode!=123&&mode!=124&&mode!=127&&mode!=128&&mode!=129&&mode!=130&&mode!=229&&mode!=230&&mode!=146&&mode!=147&&mode!=176)
{printf("Error, you can only use \"--num-blocks\" with \"--he\", \"--pcgc\", \"--fast-he\", \"--fast-pcgc\", \"--quant-her\", \"--quant-her\", \"--sum-hers\", \"--sum-cors\" or \"--jackknife\"\n\n");exit(1);}

if(num_blocks==-1&&mode!=123&&mode!=124&&mode!=129&&mode!=130&&mode!=229&&mode!=230&&mode!=176)
{printf("Error, you can only use \"--num-blocks -1\" with \"--he\", \"--pcgc\", \"--quant-her\", \"--quant-her\" or \"--jackknife\"\n\n");exit(1);}

if(permute!=-9999&&mode==133)
{printf("Error, you can not use \"--permute\" with \"--solve-null\"; you should instead permute the phenotypes manually\n\n");exit(1);}

if(permute!=-9999&&mode!=121&&mode!=123&&mode!=124&&mode!=126&&mode!=127&&mode!=128&&mode!=129&&mode!=130&&mode!=229&&mode!=230&&mode!=131&&mode!=132&&mode!=134&&mode!=138&&mode!=181&&mode!=182&&mode!=183&&mode!=184&&mode!=185)
{printf("Error, you can only use \"--permute\" with \"--reml\", \"--he\", \"--pcgc\", \"--fast-reml\", \"--fast-he\", \"--fast-pcgc\", \"--quant-her\", \"--quant-her\", \"--linear\", \"--logistic\", \"--calc-genes-reml\" or when making data\n\n");exit(1);}

if(permute==1&&strcmp(datalist,"blank")!=0)
{printf("Error, you can not use \"--permute YES\" with  \"--mbfile\", \"--msp\", \"--msped\" or \"--mspeed\"\n\n");exit(1);}

if(shrink!=-9999&&mode!=121&&mode!=138&&mode!=140&&mode!=159)
{printf("Error, you can only use \"--shrink\" with \"--reml\", \"--calc-genes-reml\", \"--mega-prs\" or \"--giga-prs\"\n\n");exit(1);}

if(shrink!=-9999&&mode!=121&&mode!=138&&mode!=140&&mode!=159)
{printf("Error, you can only use \"--shrink\" with \"--reml\", \"--calc-genes-reml\", \"--mega-prs\" or \"--giga-prs\"\n\n");exit(1);}

if(strip!=-9999&&mode!=121&&mode!=156)
{printf("Error, you can only use \"--strip\" with \"--reml\" or \"--calc-cors\"\n\n");exit(1);}

if(strip!=-9999&&mode==121&&(strcmp(kinname,"blank")!=0||strcmp(kinlist,"blank")!=0))
{printf("Error, you can not use \"--strip\" if providing a kinship matrix\n\n");exit(1);}

if(strip!=-9999&&mode==156&&num_subs!=-9999)
{printf("Error, you can not use \"--strip\" if using sample subsets\n\n");exit(1);}

if(strip!=-9999&&shrink!=-9999)
{printf("Error, you can not use both \"--strip\" and \"--shrink\"\n\n");exit(1);}

////////

if(tol!=-9999&&mode!=121&&mode!=126&&mode!=129&&mode!=130&&mode!=229&&mode!=230&&mode!=131&&mode!=133&&mode!=138&&mode!=140&&mode!=146&&mode!=147&&mode!=151&&mode!=152&&mode!=153&&mode!=154&&mode!=159&&mode!=161&&mode!=167&&mode!=179)
{printf("Error, you can only use \"--tolerance\" with \"--reml\", \"--fast-reml\", \"--quant-her\", \"--tetra-her\", \"--linear\", \"--solve-null\", \"--calc-genes-reml\", \"--sum-hers\", \"--sum-cors\", \"--ridge\", \"--bolt\", \"--bayesr\", \"--elastic\", \"--mega-prs\", \"--giga-prs\" or \"--pca\"\n\n");exit(1);}

if(maxiter!=-9999&&mode!=102&&mode!=104&&mode!=121&&mode!=126&&mode!=129&&mode!=130&&mode!=229&&mode!=230&&mode!=131&&mode!=133&&mode!=138&&mode!=140&&mode!=146&&mode!=147&&mode!=151&&mode!=152&&mode!=153&&mode!=154&&mode!=159&&mode!=179)
{printf("Error, you can only use \"--max-iter\" with \"--calc-weights\", \"--calc-weights-all\", \"--reml\", \"--fast-reml\", \"--quant-her\", \"--tetra-her\", \"--linear\", \"--solve-null\", \"--calc-genes-reml\", \"--sum-hers\", \"--sum-cors\", \"--ridge\", \"--bolt\", \"--bayesr\", \"--elastic\", \"--mega-prs\" or \"--giga-prs\"\n\n");exit(1);}

if((mode==102||mode==104)&&maxiter!=-9999&&simplex!=1)
{printf("Error, you can only use \"--max-iter\" with \"--simplex YES\"\n\n");exit(1);}

if(memsave!=-9999&&mode!=121&&mode!=123&&mode!=124&&mode!=126&&mode!=133)
{printf("Error, you can only use \"--memory-save\" with \"--reml\", \"--he\", \"--pcgc\", \"--fast-reml\" or \"--solve-null\"\n\n");exit(1);}

if(memsave==0&&diagonal==1)
{printf("Error, you can not use \"--memory-save NO\" with \"--diagonal YES\"\n\n");exit(1);}

if(memsave==0&&shrink!=-9999)
{printf("Error, you can not use \"--memory-save NO\" with \"--shrink\"\n\n");exit(1);}

////////

if(manypreds!=-9999&&mode!=151&&mode!=152&&mode!=153&&mode!=154)
{printf("Error, you can only use \"--allow-many-predictors\" with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n");exit(1);}

if(manysamples!=-9999&&mode!=102&&mode!=104&&mode!=106&&mode!=107&&mode!=110&&mode!=138&&mode!=140&&mode!=141&&mode!=156)
{printf("Error, you can only use \"--allow-many-samples\" with \"--calc-weights\", \"--calc-weights-all\", \"--calc-genes-reml\", \"--calc-tagging\", \"--calc-cors\", \"--thin\", \"--thin-tops\" or \"--thin-common\"\n\n");exit(1);}
	
//////////////////////////

