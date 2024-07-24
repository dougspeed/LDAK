/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//List of possible modes (D/d indicates data required/possible; K/k, R/r, T/t, S/s, same for kinships, regions, tops, summaries; P/p phenotypes indicates phenotypes required/possible - ignores their use for filtering)

///////////////////////////

//101 - cut-weights D
//102 - calc-weights D
//103 - join-weights D
//104 - calc-weights-all D
//105 - adjust-weights D

//106 - thin D
//107 - thin-tops D
//108 - find-tags D
//109 - remove-tags D

//111 - cut-kins D
//112 - calc-kins D
//113 - join-kins K
//114 - calc-kins-direct D

//115 - filter K p
//116 - add-grm K
//117 - sub-grm K d
//118 - convert-gz K
//119 - convert-raw K
//120 - calc-sim-grms K

//121 - reml krt sp
//122 - calc-blups D k p
//123 - he krt P
//124 - pcgc krt P
//125 - reml-pred K p

//126 - fast-reml Krt P
//127 - fast-he D P
//128 - fast-pcgc D P
//129 - quant-her P
//130 - tetra-her P

//131 - linear D kt P
//132 - logistic D t P
//133 - solve-null K rt P
//134 - kvik-step2 D P (temporary)

//136 - cut-genes D
//137 - calc-genes-kins D
//138 - calc-genes-reml D kt sp
//139 - join-genes-reml d
//140 - kvik-step3 D

//141 - calc-tagging D
//142 - join-tagging
//143 - merge-tagging
//144 - reduce-tagging
//145 - calc-overlaps

//146 - sum-hers S
//147 - sum-cors S

//149 - calc-exps
//150 - calc-posts S

//151 - ridge D P
//152 - bolt D P
//153 - bayesr D P
//154 - elastic D P

//156 - calc-cors D
//157 - join-cors
//158 - pseudo-summaries d S
//159 - mega-prs S
//160 - validate D P

//161 - pca K
//162 - calc-pca-loads D K
//163 - decompose K
//164 - adjust-grm K

//166 - truncate-grm K
//167 - pca-grm K
//168 - square-grm K
//169 - gxemm-iid K
//170 - gxemm-free K

//171 - calc-stats D
//172 - calc-scores D sp
//173 - make-phenos D
//174 - make-snps
//175 - calc-inflation D

//176 - jackknife
//177 - cut-folds dk
//178 - find-gaussian
//179 - winners-curse S

//181 - make-bed D+
//182 - make-sp D+
//183 - make-sped D+
//184 - make-speed D+
//185 - make-gen D+

//186 - condense-bed D
//187 - condense-sp D
//188 - condense-sped D
//189 - condense-speed D
//190 - calc-sim-data DD

//191 - cut-gre D
//192 - calc-gre D
//193 - join-gre D
//194 - solve-gre DP

//229 - quant-bivar P
//230 - tetra-bivar P

//201 - speed-tests
//202 - speed-tests2
//203 - speed-tests3

///////////////////////////

