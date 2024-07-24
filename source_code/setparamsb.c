/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//set parameters for prediction models (megaprs)

///////////////////////////

//set num_try
if(strcmp(bestfile,"blank")!=0)	//considering a single model
{
num_try=1;
}
else
{
if(strcmp(fracfile,"blank")==0)	//defaults depend on mode and whether her is specified
{
if(ptype==0)	//mega - do lasso, ridge, bayesr (minus ridge) and elastic
{
if(her!=-9999){num_try=45;}
else{num_try=143;}
}
if(ptype==1)	//lasso-sparse
{
num_try=80;
}
if(ptype==2)	//lasso
{
if(her!=-9999){num_try=1;}
else{num_try=7;}
}
if(ptype==3)	//ridge
{
if(her!=-9999){num_try=1;}
else{num_try=7;}
}
if(ptype==4)	//bolt
{
if(her!=-9999){num_try=18;}
else{num_try=54;}
}
if(ptype==5)	//bayesr
{
if(her!=-9999){num_try=35;}
else{num_try=105;}
}
if(ptype==6)	//bayesr-shrink
{
if(her!=-9999){num_try=35;}
else{num_try=105;}
}
if(ptype==7)	//elastic
{
if(her!=-9999){num_try=9;}
else{num_try=27;}
}
}
else	//will set models based on values in fracfile
{
num_try=countrows(fracfile);
}
}

//allocate parameters
trytypes=malloc(sizeof(int)*num_try);
tryhers=malloc(sizeof(double)*num_try);
trylams=malloc(sizeof(double)*num_try);
tryscales=malloc(sizeof(double)*num_try);
tryps=malloc(sizeof(double)*num_try);
tryp2s=malloc(sizeof(double)*num_try);
tryp3s=malloc(sizeof(double)*num_try);
tryp4s=malloc(sizeof(double)*num_try);
tryf2s=malloc(sizeof(double)*num_try);

//set scale to 1, then remaining values to NA (then will change when required)
for(p=0;p<num_try;p++)
{
trytypes[p]=-9999;tryhers[p]=-9999;trylams[p]=-9999;tryscales[p]=1;
tryps[p]=-9999;tryp2s[p]=-9999;tryp3s[p]=-9999;tryp4s[p]=-9999;tryf2s[p]=-9999;
}

//set num_try
if(strcmp(bestfile,"blank")!=0)	//considering a single model
{
//open bestfile, skip the header row and get model
if((input=fopen(bestfile,"r"))==NULL)
{printf("Error opening %s\n\n",bestfile);exit(1);}
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}
if(fscanf(input, "%d %s ", &readint, readstring)!=2)
{printf("Error reading Row 2 of %s, suggesting the file has been changed since creation with \"--mega-prs\"\n\n", bestfile);exit(1);}
if(strcmp(readstring,"lasso-sparse")!=0&&strcmp(readstring,"lasso")!=0&&strcmp(readstring,"ridge")!=0&&strcmp(readstring,"bolt")!=0 &&strcmp(readstring,"bayesr")!=0&&strcmp(readstring,"bayesr-shrink")!=0&&strcmp(readstring,"elastic")!=0)
{printf("Error, Element 2 of Row 2 of %s should be  \"lasso-sparse\", \"lasso\", \"ridge\", \"bolt\", \"bayesr\", \"bayesr-shrink\" or \"--elastic\" (not %s), suggesting the file has been changed since creation with \"--mega-prs\"\n\n", bestfile, readstring);exit(1);}

if(strcmp(readstring,"lasso-sparse")==0)	//lasso-sparse - read lambda and s
{
trytypes[0]=1;
if(fscanf(input, "NA %lf %lf ", trylams, tryscales)!=2){printf("Error reading lambda and s from Row 2 of %s, suggesting the file has been changed since creation with \"--mega-prs\"\n\n", bestfile);exit(1);}
}

if(strcmp(readstring,"lasso")==0)	//lasso-shrink - read heritability
{
trytypes[0]=2;
if(fscanf(input, "%lf ", tryhers)!=1){printf("Error reading heritability from Row 2 of %s, suggesting the file has been changed since creation with \"--mega-prs\"\n\n", bestfile);exit(1);}
}

if(strcmp(readstring,"ridge")==0)	//ridge - read heritability
{
trytypes[0]=3;
if(fscanf(input, "%lf ", tryhers)!=1){printf("Error reading heritability from Row 2 of %s, suggesting the file has been changed since creation with \"--mega-prs\"\n\n", bestfile);exit(1);}
}

if(strcmp(readstring,"bolt")==0)	//bolt - read heritability, p and f2 (and set p2/p3)
{
trytypes[0]=4;
if(fscanf(input, "%lf NA NA %lf %lf ", tryhers, tryps, tryf2s)!=3){printf("Error reading heritability, p and f2 from Row 2 of %s, suggesting the file has been changed since creation with \"--mega-prs\"\n\n", bestfile);exit(1);}
tryp2s[0]=1-tryps[0];
}

if(strcmp(readstring,"bayesr")==0)	//bayesr-sparse - read heritability, p, p2, p3 and p4
{
trytypes[0]=5;
if(fscanf(input, "%lf NA NA NA NA %lf %lf %lf %lf ", tryhers, tryps, tryp2s, tryp3s, tryp4s)!=5){printf("Error reading heritability, p, p2, p3 and p4 from Row 2 of %s, suggesting the file has been changed since creation with \"--mega-prs\"\n\n", bestfile);exit(1);}
}

if(strcmp(readstring,"bayesr-shrink")==0)	//bayesr-shrink - read heritability, p, p2, p3 and p4
{
trytypes[0]=6;
if(fscanf(input, "%lf NA NA NA NA %lf %lf %lf %lf ", tryhers, tryps, tryp2s, tryp3s, tryp4s)!=5){printf("Error reading heritability, p, p2, p3 and p4 from Row 2 of %s, suggesting the file has been changed since creation with \"--mega-prs\"\n\n", bestfile);exit(1);}
}

if(strcmp(readstring,"elastic")==0)	//elastic - read heritability, p and f2 (and set p2/p3)
{
trytypes[0]=4;
if(fscanf(input, "%lf NA NA %lf %lf ", tryhers, &readdouble, tryf2s)!=3){printf("Error reading heritability, p and f2 from Row 2 of %s, suggesting the file has been changed since creation with \"--mega-prs\"\n\n", bestfile);exit(1);}
tryps[0]=readdouble/2;
tryp2s[0]=readdouble/2;
tryp3s[0]=1-readdouble;
}

fclose(input);
}
else
{
if(strcmp(fracfile,"blank")==0)	//using defaults
{
count=0;

if(ptype==0)	//mega - using lasso, ridge, bayesr and elastic
{
//lasso and ridge
if(her!=-9999)
{
trytypes[count]=2;tryhers[count]=her;count++;
trytypes[count]=3;tryhers[count]=her;count++;
}
else
{
//lasso and ridge	
for(j=1;j<8;j++){trytypes[count]=2;tryhers[count]=0.1*j;count++;}
for(j=1;j<8;j++){trytypes[count]=3;tryhers[count]=0.1*j;count++;}
}

//bayesr (minus ridge)
loads=malloc(sizeof(double)*5);
loads[0]=0;loads[1]=.01;loads[2]=.05;loads[3]=.1;loads[4]=.2;
if(her!=-9999)
{
for(j=0;j<5;j++){
for(j2=j;j2<5;j2++){
for(j3=j2;j3<5;j3++){
if(j+j2+j3>0)
{
trytypes[count]=5;tryhers[count]=her;tryp4s[count]=loads[j];tryp3s[count]=loads[j2];
tryp2s[count]=loads[j3];tryps[count]=1-tryp2s[count]-tryp3s[count]-tryp4s[count];count++;
}
}}}
}
else
{
for(j=0;j<5;j++){
for(j2=j;j2<5;j2++){
for(j3=j2;j3<5;j3++){
if(j+j2+j3>0)
{
trytypes[count]=5;tryhers[count]=0.1;tryp4s[count]=loads[j];tryp3s[count]=loads[j2];
tryp2s[count]=loads[j3];tryps[count]=1-tryp2s[count]-tryp3s[count]-tryp4s[count];count++;
trytypes[count]=5;tryhers[count]=0.3;tryp4s[count]=loads[j];tryp3s[count]=loads[j2];
tryp2s[count]=loads[j3];tryps[count]=1-tryp2s[count]-tryp3s[count]-tryp4s[count];count++;
trytypes[count]=5;tryhers[count]=0.5;tryp4s[count]=loads[j];tryp3s[count]=loads[j2];
tryp2s[count]=loads[j3];tryps[count]=1-tryp2s[count]-tryp3s[count]-tryp4s[count];count++;
}
}}}
}
free(loads);

//elastic
loads=malloc(sizeof(double)*3);
loads2=malloc(sizeof(double)*3);
loads[0]=.5;loads[1]=.1;loads[2]=.01;
loads2[0]=.5;loads2[1]=.3;loads2[2]=.1;
if(her!=-9999)
{
for(j=0;j<3;j++)
{
for(j2=0;j2<3;j2++)
{
trytypes[count]=7;tryhers[count]=her;tryps[count]=loads[j]/2;tryp2s[count]=loads[j]/2;
tryp3s[count]=1-loads[j];tryf2s[count]=loads2[j2];count++;
}}
}
else
{
for(j=0;j<3;j++)
{
for(j2=0;j2<3;j2++)
{
trytypes[count]=7;tryhers[count]=0.1;tryps[count]=loads[j]/2;tryp2s[count]=loads[j]/2;
tryp3s[count]=1-loads[j];tryf2s[count]=loads2[j2];count++;
trytypes[count]=7;tryhers[count]=0.3;tryps[count]=loads[j]/2;tryp2s[count]=loads[j]/2;
tryp3s[count]=1-loads[j];tryf2s[count]=loads2[j2];count++;
trytypes[count]=7;tryhers[count]=0.5;tryps[count]=loads[j]/2;tryp2s[count]=loads[j]/2;
tryp3s[count]=1-loads[j];tryf2s[count]=loads2[j2];count++;
}}
}
free(loads);free(loads2);
}	//end of ptype=0

if(ptype==1)	//lasso-sparse (her does not feature)
{
for(p=0;p<20;p++){trytypes[count]=1;tryscales[count]=0.9;trylams[count]=0.001*pow(100,(double)p/19);count++;}
for(p=0;p<20;p++){trytypes[count]=1;tryscales[count]=0.5;trylams[count]=0.001*pow(100,(double)p/19);count++;}
for(p=0;p<20;p++){trytypes[count]=1;tryscales[count]=0.2;trylams[count]=0.001*pow(100,(double)p/19);count++;}
for(p=0;p<20;p++){trytypes[count]=1;tryscales[count]=0.1;trylams[count]=0.001*pow(100,(double)p/19);count++;}
}

if(ptype==2)	//lasso
{
if(her!=-9999)
{
trytypes[count]=2;tryhers[count]=her;count++;
}
else
{	
for(j=1;j<8;j++){trytypes[count]=2;tryhers[count]=0.1*j;count++;}
}
}

if(ptype==3)	//ridge
{
if(her!=-9999)
{
trytypes[count]=3;tryhers[count]=her;count++;
}
else
{	
for(j=1;j<8;j++){trytypes[count]=3;tryhers[count]=0.1*j;count++;}
}
}

if(ptype==4)	//bolt
{
loads=malloc(sizeof(double)*6);
loads2=malloc(sizeof(double)*3);
loads[0]=.5;loads[1]=.2;loads[2]=.1;loads[3]=.05;loads[4]=.02;loads[5]=.01;
loads2[0]=.5;loads2[1]=.3;loads2[2]=.1;
if(her!=-9999)
{
for(j=0;j<6;j++)
{
for(j2=0;j2<3;j2++)
{
trytypes[count]=4;tryhers[count]=her;tryps[count]=loads[j];tryp2s[count]=1-loads[j];tryf2s[count]=loads2[j2];count++;
}}
}
else
{
for(j=0;j<6;j++)
{
for(j2=0;j2<3;j2++)
{
trytypes[count]=4;tryhers[count]=0.1;tryps[count]=loads[j];tryp2s[count]=1-loads[j];tryf2s[count]=loads2[j2];count++;
trytypes[count]=4;tryhers[count]=0.3;tryps[count]=loads[j];tryp2s[count]=1-loads[j];tryf2s[count]=loads2[j2];count++;
trytypes[count]=4;tryhers[count]=0.5;tryps[count]=loads[j];tryp2s[count]=1-loads[j];tryf2s[count]=loads2[j2];count++;
}}
}
free(loads);free(loads2);
}	//end of ptype=4

if(ptype==5)	//bayesr (make sure to include ridge)
{
loads=malloc(sizeof(double)*5);
loads[0]=0;loads[1]=.01;loads[2]=.05;loads[3]=.1;loads[4]=.2;
if(her!=-9999)
{
for(j=0;j<5;j++){
for(j2=j;j2<5;j2++){
for(j3=j2;j3<5;j3++){
if(j+j2+j3>0)
{
trytypes[count]=5;tryhers[count]=her;tryp4s[count]=loads[j];tryp3s[count]=loads[j2];
tryp2s[count]=loads[j3];tryps[count]=1-tryp2s[count]-tryp3s[count]-tryp4s[count];count++;
}
}}}
trytypes[count]=5;tryhers[count]=her;tryp4s[count]=1;tryp3s[count]=0;
tryp2s[count]=0;tryps[count]=0;count++;
}
else
{	
for(j=0;j<5;j++){
for(j2=j;j2<5;j2++){
for(j3=j2;j3<5;j3++){
if(j+j2+j3>0)
{
trytypes[count]=5;tryhers[count]=0.1;tryp4s[count]=loads[j];tryp3s[count]=loads[j2];
tryp2s[count]=loads[j3];tryps[count]=1-tryp2s[count]-tryp3s[count]-tryp4s[count];count++;
trytypes[count]=5;tryhers[count]=0.3;tryp4s[count]=loads[j];tryp3s[count]=loads[j2];
tryp2s[count]=loads[j3];tryps[count]=1-tryp2s[count]-tryp3s[count]-tryp4s[count];count++;
trytypes[count]=5;tryhers[count]=0.5;tryp4s[count]=loads[j];tryp3s[count]=loads[j2];
tryp2s[count]=loads[j3];tryps[count]=1-tryp2s[count]-tryp3s[count]-tryp4s[count];count++;
}
}}}
trytypes[count]=5;tryhers[count]=0.1;tryp4s[count]=1;tryp3s[count]=0;
tryp2s[count]=0;tryps[count]=0;count++;
trytypes[count]=5;tryhers[count]=0.3;tryp4s[count]=1;tryp3s[count]=0;
tryp2s[count]=0;tryps[count]=0;count++;
trytypes[count]=5;tryhers[count]=0.5;tryp4s[count]=1;tryp3s[count]=0;
tryp2s[count]=0;tryps[count]=0;count++;
}
free(loads);

//squeeze down the models (make sure tryps4 non-zero)
for(j=0;j<count;j++)
{
if(tryp4s[j]==0)
{
if(tryp3s[j]>0){tryp4s[j]=tryp3s[j];tryp3s[j]=tryp2s[j];tryp2s[j]=0;}
else	//so p4 and p3 are zero
{
if(tryp2s[j]>0){tryp4s[j]=tryp2s[j];tryp2s[j]=0;}
else	//p4, p3 and p2 are all zero
{tryp4s[j]=tryps[j];tryps[j]=0;}
}
}
}
}	//end of ptype=5

if(ptype==6)	//bayes-shrink
{
loads=malloc(sizeof(double)*5);
loads[0]=0;loads[1]=.01;loads[2]=.05;loads[3]=.1;loads[4]=.2;
if(her!=-9999)
{
for(j=0;j<5;j++){
for(j2=j;j2<5;j2++){
for(j3=j2;j3<5;j3++){
trytypes[count]=6;tryhers[count]=her;tryp4s[count]=loads[j];tryp3s[count]=loads[j2];
tryp2s[count]=loads[j3];tryps[count]=1-tryp2s[count]-tryp3s[count]-tryp4s[count];count++;
}}}
}
else
{
loads=malloc(sizeof(double)*5);
loads[0]=0;loads[1]=.01;loads[2]=.05;loads[3]=.1;loads[4]=.2;
for(j=0;j<5;j++){
for(j2=j;j2<5;j2++){
for(j3=j2;j3<5;j3++){
trytypes[count]=6;tryhers[count]=0.1;tryp4s[count]=loads[j];tryp3s[count]=loads[j2];
tryp2s[count]=loads[j3];tryps[count]=1-tryp2s[count]-tryp3s[count]-tryp4s[count];count++;
trytypes[count]=6;tryhers[count]=0.3;tryp4s[count]=loads[j];tryp3s[count]=loads[j2];
tryp2s[count]=loads[j3];tryps[count]=1-tryp2s[count]-tryp3s[count]-tryp4s[count];count++;
trytypes[count]=6;tryhers[count]=0.5;tryp4s[count]=loads[j];tryp3s[count]=loads[j2];
tryp2s[count]=loads[j3];tryps[count]=1-tryp2s[count]-tryp3s[count]-tryp4s[count];count++;
}}}
}
free(loads);

//squeeze down the models (make sure tryps4 non-zero)
for(j=0;j<count;j++)
{
if(tryp4s[j]==0)
{
if(tryp3s[j]>0){tryp4s[j]=tryp3s[j];tryp3s[j]=tryp2s[j];tryp2s[j]=0;}
else	//so p4 and p3 are zero
{
if(tryp2s[j]>0){tryp4s[j]=tryp2s[j];tryp2s[j]=0;}
else	//p4, p3 and p2 are all zero
{tryp4s[j]=tryps[j];tryps[j]=0;}
}
}
}
}	//end of ptype=6

if(ptype==7)	//elastic
{
loads=malloc(sizeof(double)*3);
loads2=malloc(sizeof(double)*3);
loads[0]=.5;loads[1]=.1;loads[2]=.01;
loads2[0]=.5;loads2[1]=.3;loads2[2]=.1;
if(her!=-9999)
{
for(j=0;j<3;j++)
{
for(j2=0;j2<3;j2++)
{
trytypes[count]=7;tryhers[count]=her;tryps[count]=loads[j]/2;tryp2s[count]=loads[j]/2;
tryp3s[count]=1-loads[j];tryf2s[count]=loads2[j2];count++;
}}
}
else
{
for(j=0;j<3;j++)
{
for(j2=0;j2<3;j2++)
{
trytypes[count]=7;tryhers[count]=0.1;tryps[count]=loads[j]/2;tryp2s[count]=loads[j]/2;
tryp3s[count]=1-loads[j];tryf2s[count]=loads2[j2];count++;
trytypes[count]=7;tryhers[count]=0.3;tryps[count]=loads[j]/2;tryp2s[count]=loads[j]/2;
tryp3s[count]=1-loads[j];tryf2s[count]=loads2[j2];count++;
trytypes[count]=7;tryhers[count]=0.5;tryps[count]=loads[j]/2;tryp2s[count]=loads[j]/2;
tryp3s[count]=1-loads[j];tryf2s[count]=loads2[j2];count++;
}}
}
free(loads);free(loads2);
}	//end of ptype=7

if(count!=num_try){printf("Doug error %d and %d\n", num_try, count);exit(1);}
}	//end of setting using defaults

else	//read fracfile
{
for(p=0;p<num_try;p++){trytypes[p]=ptype;}

if(ptype==2||ptype==3)	//lasso or ridge - read heritability
{
read_values(fracfile, tryhers, num_try, NULL, 1, 0, 0);
for(p=0;p<num_try;p++)
{
if(tryhers[p]<=0||tryhers[p]>=1){printf("Error, heritability in Row %d of %s (%.4f) is not within (0,1)\n\n", p+1, fracfile, tryhers[p]);exit(1);}
}
}

if(ptype==4)	//bolt - read heritability, p and f2 (also set p2)
{
read_values(fracfile, tryhers, num_try, NULL, 1, 0, 0);
for(p=0;p<num_try;p++)
{
if(tryhers[p]<=0||tryhers[p]>=1){printf("Error, heritability in Row %d of %s (%.4f) is not within (0,1)\n\n", p+1, fracfile, tryhers[p]);exit(1);}
}

read_values(fracfile, tryps, num_try, NULL, 2, 0, 0);
for(p=0;p<num_try;p++)
{
if(tryps[p]<=0||tryps[p]>=1){printf("Error, probability in Row %d of %s (%.4f) is not within (0,1=\n\n", p+1, fracfile, tryps[p]);exit(1);}
}
for(p=0;p<num_try;p++){tryp2s[p]=1-tryps[p];}

read_values(fracfile, tryf2s, num_try, NULL, 3, 0, 0);
for(p=0;p<num_try;p++)
{
if(tryf2s[p]<0||tryf2s[p]>1){printf("Error, fraction in Row %d of %s (%.4f) is not within [0,1]\n\n", p+1, fracfile, tryf2s[p]);exit(1);}
}
}

if(ptype==5||ptype==6)	//bayesr or bayesr-shrink - read heritability, p, p2, p3 and p4
{
read_values(fracfile, tryhers, num_try, NULL, 1, 0, 0);
for(p=0;p<num_try;p++)
{
if(tryhers[p]<=0||tryhers[p]>=1){printf("Error, heritability in Row %d of %s (%.4f) is not within (0,1)\n\n", p+1, fracfile, tryhers[p]);exit(1);}
}

if(ptype==7)	//elastic - read heritability, 2p and f2 (also set p2 and p3)
{
read_values(fracfile, tryhers, num_try, NULL, 1, 0, 0);
for(p=0;p<num_try;p++)
{
if(tryhers[p]<0||tryhers[p]>1){printf("Error, heritability in Row %d of %s (%.4f) is not within [0,1]\n\n", p+1, fracfile, tryhers[p]);exit(1);}
}

read_values(fracfile, tryps, num_try, NULL, 2, 0, 0);
for(p=0;p<num_try;p++)
{
if(tryps[p]<0||tryps[p]>1){printf("Error, probability in Row %d of %s (%.4f) is not within [0,1]\n\n", p+1, fracfile, tryps[p]);exit(1);}
}
for(p=0;p<num_try;p++){tryps[p]=tryps[p]/2;tryp2s[p]=tryps[p];tryp3s[p]=1-tryps[p]-tryp2s[p];}

read_values(fracfile, tryf2s, num_try, NULL, 3, 0, 0);
for(p=0;p<num_try;p++)
{
if(tryf2s[p]<0||tryf2s[p]>1){printf("Error, fraction in Row %d of %s (%.4f) is not within [0,1]\n\n", p+1, fracfile, tryf2s[p]);exit(1);}
if(tryps[p]==1&&tryf2s[p]>0){printf("Error in Row %d of %s; can not have p=1 and f2>0\n\n", p+1, fracfile);exit(1);}
if(tryps[p]==0&&tryf2s[p]<1){printf("Error in Row %d of %s; can not have p=0 and f2<1\n\n", p+1, fracfile);exit(1);}
}
}

read_values(fracfile, tryps, num_try, NULL, 2, 0, 0);
read_values(fracfile, tryp2s, num_try, NULL, 3, 0, 0);
read_values(fracfile, tryp3s, num_try, NULL, 4, 0, 0);
read_values(fracfile, tryp4s, num_try, NULL, 5, 0, 0);

for(p=0;p<num_try;p++)
{
sum=tryps[p]+tryp2s[p]+tryp3s[p]+tryp4s[p];
if(sum<0.99||sum>1.01)
{printf("Error, probabilities in Row %d of %s sum to %.4f (not one)\n\n", p+1, fracfile, sum);exit(1);}
tryps[p]=tryps[p]/sum;
tryp2s[p]=tryp2s[p]/sum;
tryp3s[p]=tryp3s[p]/sum;
tryp4s[p]=tryp4s[p]/sum;
}
}
}	//end of reading fracfile
}	//end of not using best model

//save parameters
sprintf(filename,"%s.parameters",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}

fprintf(output, "Model\tType\tHeritability\tLasso_lambda\tLasso_s\tBolt_p\tBolt_f2\tBayesR_p1\tBayesR_p2\tBayesR_p3\tBayesR_p4\tElastic_p\tElastic_f2\n");
for(p=0;p<num_try;p++)
{
if(trytypes[p]==1){fprintf(output, "%d\tlasso-sparse\tNA\t%.4f\t%.4f\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n", p+1, trylams[p], tryscales[p]);}
if(trytypes[p]==2){fprintf(output, "%d\tlasso\t%.4f\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n", p+1, tryhers[p]);}
if(trytypes[p]==3){fprintf(output, "%d\tridge\t%.4f\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n", p+1, tryhers[p]);}
if(trytypes[p]==4){fprintf(output, "%d\tbolt\t%.4f\tNA\tNA\t%.4f\t%.4f\tNA\tNA\tNA\tNA\tNA\tNA\n", p+1, tryhers[p], tryps[p], tryf2s[p]);}
if(trytypes[p]==5){fprintf(output, "%d\tbayesr\t%.4f\tNA\tNA\tNA\tNA\t%.4f\t%.4f\t%.4f\t%.4f\tNA\tNA\n", p+1, tryhers[p], tryps[p], tryp2s[p], tryp3s[p], tryp4s[p]);}
if(trytypes[p]==6){fprintf(output, "%d\tbayesr-shrink\t%.4f\tNA\tNA\tNA\tNA\t%.4f\t%.4f\t%.4f\t%.4f\tNA\tNA\n", p+1, tryhers[p], tryps[p], tryp2s[p], tryp3s[p], tryp4s[p]);}
if(trytypes[p]==7){fprintf(output, "%d\telastic\t%.4f\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t%.4f\t%.4f\n", p+1, tryhers[p], 1-tryp3s[p], tryf2s[p]);}
}
fclose(output);

///////////////////////////

