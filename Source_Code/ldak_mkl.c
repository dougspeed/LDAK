/*
Copyright 2026 Doug Speed.
 
    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.
*/

///////////////////////////

//Wrapper file - this sets the variables MKL to 0 and MET to 0, then switches to ldak.c

///////////////////////////

#define MKL 1	//1 to compile with mkl, 0 to compile without mkl, 2 to compile with AOCL
#define MET 0	//0 to compile with qsopt (required if you wish to calculate weightings)
#define DRM 0	//1 to allow DRM

#include "ldak.c"

///////////////////////////

