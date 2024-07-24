/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Code for quickly merging bim files

///////////////////////////

//work out from which file each predictor comes
indexer=malloc(sizeof(int)*num_preds_use);
indexer2=malloc(sizeof(int)*num_preds_use);
for(k=0;k<num_files;k++)
{
for(j=0;j<XNuse[k];j++){indexer[Xkp2[k][j]]=k;indexer2[Xkp2[k][j]]=Xkp[k][j];}
}

if(bitsize>num_preds_use){bitsize=num_preds_use;}

//have checked that all bim files contain num_samples_use samples
rowlength=(num_samples_use-1)/4+1;
rowchars=malloc(sizeof(unsigned char)*rowlength);

//open all bim files, check their length and first three digits, and move to first predictor
Xinput=malloc(sizeof(FILE*)*num_files);
Xcurrent=malloc(sizeof(int)*num_files);

for(k=0;k<num_files;k++)
{
if((Xinput[k]=fopen(datastems[k],"rb"))==NULL)
{printf("Error opening %s\n\n",datastems[k]);exit(1);}

fseeko(Xinput[k], 0, SEEK_END);
if(ftello(Xinput[k])!=(off_t)sizeof(unsigned char)*rowlength*XNall[k]+sizeof(unsigned char)*3)
{printf("Error reading %s; should have size %jd (%d ind x %d predictors), but instead has size %jd\n\n", datastems[k], (off_t)sizeof(unsigned char)*rowlength*XNall[k]+sizeof(unsigned char)*3, num_samples_use, XNall[k], ftello(Xinput[k]));exit(1);}

fseeko(Xinput[k], 0, SEEK_SET);
if(fread(startchars, sizeof(unsigned char), 3, Xinput[k])!=3)
{printf("Error reading first three values of %s\n\n", datastems[k]);exit(1);}
if(startchars[0]!=108||startchars[1]!=27)
{printf("Error reading %s; does not appear to be in binary PLINK format\n\n", datastems[k]);exit(1);}
if(startchars[2]!=1)
{printf("Error reading %s; can only read in SNP-major mode\n\n", datastems[k]);exit(1);}

Xcurrent[k]=0;
}

//open output files
sprintf(filename,"%s.progress", outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fclose(output);

sprintf(filename2,"%s.bed", outfile);
if((output2=fopen(filename2,"wb"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
onechar=108;fwrite(&onechar, sizeof(unsigned char), 1, output2);
onechar=27;fwrite(&onechar, sizeof(unsigned char), 1, output2);
onechar=1;fwrite(&onechar, sizeof(unsigned char), 1, output2);

////////

//ready for bit loop
bittotal=(num_preds_use-1)/bitsize+1;
for(bit=0;bit<bittotal;bit++)
{
bitstart=bit*bitsize;
bitend=(bit+1)*bitsize;
if(bitend>num_preds_use){bitend=num_preds_use;}
bitlength=bitend-bitstart;

printf("Making data for Chunk %d of %d\n", bit+1, bittotal);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output, "Making data for Chunk %d of %d\n", bit+1, bittotal);
fclose(output);

for(j=bitstart;j<bitend;j++)
{
k=indexer[j];j2=indexer2[j];
if(j2!=Xcurrent[k])
{
if(fseeko(Xinput[k], (off_t)sizeof(unsigned char)*rowlength*j2+sizeof(unsigned char)*3, SEEK_SET)!=0)
{printf("Error reading %s; unable to find Predictor %d\n\n", datastems[k], j2+1);exit(1);}
}
if(fread(rowchars, sizeof(unsigned char), rowlength, Xinput[k])!=rowlength)
{printf("Error reading values for Predictor %d from %s\n\n", j2+1, datastems[k]);exit(1);}
Xcurrent[k]=j2+1;

fwrite(rowchars, sizeof(unsigned char), rowlength, output2);
}
}
printf("\n");

for(k=0;k<num_files;k++){fclose(Xinput[k]);}
fclose(output2);

//write bim and fam files
sprintf(filename3,"%s.bim", outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
for(j=0;j<num_preds_use;j++)
{
fprintf(output3, "%d\t%s\t", chr[j], preds[j]);
if(cm[j]==0){fprintf(output3, "0\t");}
else{fprintf(output3, "%.6f\t", cm[j]);}
fprintf(output3, "%ld\t%c\t%c\n", (long int)bp[j], al1[j], al2[j]);
}
fclose(output3);

sprintf(filename4,"%s.fam", outfile);
if((output4=fopen(filename4,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename4);exit(1);}
for(i=0;i<num_samples_use;i++)
{fprintf(output4, "%s %s %s %s %s %s\n", ids1[i], ids2[i], pid[i], mid[i], schar[i], pchar[i]);}
fclose(output4);

free(indexer);free(indexer2);
free(rowchars);
free(Xinput);free(Xcurrent);

printf("Data for %d samples and %d predictors saved in %s, %s and %s\n\n", num_samples_use, num_preds_use, filename2, filename3, filename4);

///////////////////////////

