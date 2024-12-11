/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Save data for one predictor

///////////////////////////

if(mode==184)	//work out min and max and print header columns
{
if(nonsnp==0){minfloat=0;maxfloat=2;}
else
{
minfloat=missingvalue;maxfloat=missingvalue;
for(i=0;i<num_samples_use;i++)
{
value=data[(size_t)j*num_samples_use+i];
if(value!=missingvalue)
{
if(minfloat==missingvalue){minfloat=value;maxfloat=value;}
if(value<minfloat){minfloat=value;}
if(value>maxfloat){maxfloat=value;}
}
}
if(minfloat==missingvalue){minfloat=0;maxfloat=2;}
if(minfloat==maxfloat){maxfloat++;}
}
if(speedlong==0){writefloat=1;}
else{writefloat=2;}
fwrite(&writefloat, sizeof(float), 1, output3);
fwrite(&minfloat, sizeof(float), 1, output3);
fwrite(&maxfloat, sizeof(float), 1, output3);
writefloat=0;
for(k=0;k<13;k++){fwrite(&writefloat, sizeof(float), 1, output3);}
}

if(mode==185)	//print headers
{fprintf(output3, "%d %s %ld %c %c ", chr[bitstart+j], preds[bitstart+j], (long int)bp[bitstart+j], al1[bitstart+j], al2[bitstart+j]);}

for(i=0;i<num_samples_use;i++)
{
value=data[(size_t)j*num_samples_use+i];

if(mode==181)
{
if(i%4==0){onechar=0;}
if(value==0){onechar+=(3<<(2*(i%4)));}
if(value==1){onechar+=(2<<(2*(i%4)));}
//if(value==2){onechar+=(0<<(2*(i%4)));}
if(value==missingvalue){onechar+=(1<<(2*(i%4)));}
if(i%4==3||i==num_samples_use-1){fwrite(&onechar, sizeof(unsigned char), 1, output3);}
}
if(mode==182)
{
if(value!=missingvalue)
{
if(fabs(value-round(value))<0.00005){fprintf(output3,"%d ", (int)round(value));}
else{fprintf(output3,"%.4f ", value);}
}
else
{fprintf(output3,"NA ");}
}
if(mode==183)
{
writefloat=(float)value;
fwrite(&writefloat, sizeof(float), 1, output3);
}
if(mode==184)
{
if(speedlong==0)	//use 256
{
if(value==missingvalue){onechar=255;}
else{onechar=round((value-minfloat)/(maxfloat-minfloat)*254);}
fwrite(&onechar, sizeof(unsigned char), 1, output3);
}
else	//use 65536
{
if(value==missingvalue){oneshort=65535;}
else{oneshort=round((value-minfloat)/(maxfloat-minfloat)*65534);}
fwrite(&oneshort, sizeof(unsigned short), 1, output3);
}
}
if(mode==185)
{
if(value==missingvalue){fprintf(output3,"0.33 0.33 0.33 ");}
else
{
value=ps[0][(size_t)j*num_samples_use+i];
if(value==0||value==1){fprintf(output3,"%d ", (int)value);}
else{fprintf(output3,"%.3f ", value);}
value=ps[1][(size_t)j*num_samples_use+i];
if(value==0||value==1){fprintf(output3,"%d ", (int)value);}
else{fprintf(output3,"%.3f ", value);}
value=1-ps[0][(size_t)j*num_samples_use+i]-ps[1][(size_t)j*num_samples_use+i];
if(value==0||value==1){fprintf(output3,"%d ", (int)value);}
else{fprintf(output3,"%.3f ", value);}
}
}
}	//end of i loop
if(mode==182||mode==185){fprintf(output3,"\n");}

///////////////////////////


