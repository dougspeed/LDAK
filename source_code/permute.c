/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Permutation functions

///////////////////////////

void permute_int(int *array, int length)
//permutes integer array of size length
{
int i, j;
int temp;

for(i=length-1;i>0;i--)
{
j=rand()%(i+1);
temp=array[i];
array[i]=array[j];
array[j]=temp;
}
}

////////

void permute_double(double *array, int length)
//permutes double array of size length
{
int i, j;
double temp;

for(i=length-1;i>0;i--)
{
j=rand()%(i+1);
temp=array[i];
array[i]=array[j];
array[j]=temp;
}
}

////////

void permute_float(float *array, int length)
//permutes float array of size length
{
int i, j;
float temp;

for(i=length-1;i>0;i--)
{
j=rand()%(i+1);
temp=array[i];
array[i]=array[j];
array[j]=temp;
}
}

///////////////////////////

