*code originally written by Jennifer Morse and editted for use with EIFB data by Kimberly Nolf.  The code will recode variable type from string to numeric, 
*recode reverse items, and create scores for Anxiety and Avoidance scales.

*changes string variables to numeric variables.
ALTER TYPE ecr1 to ecr36 (F2.0) .
execute .

ALTER TYPE Subject (F4.0) .
execute .

MISSING VALUES
 ecr1,ecr2,ecr3,ecr4,ecr5,ecr6,ecr7,ecr8,ecr9,ecr10,
 ecr11,ecr12,ecr13,ecr14,ecr15,ecr16,ecr17,ecr18,ecr19,ecr20,
 ecr21,ecr22,ecr23,ecr24,ecr25,ecr26,ecr27,ecr28,ecr29,ecr30,
 ecr31,ecr32,ecr33,ecr34,ecr35,ecr36 (9) .
EXECUTE .

*recode the reverse items for Anxiety (2, 22) and Avoidance (3, 5, 11,15, 17, 19, 25, 27, 29, 31, 33, 35) scales.
RECODE 
ecr2,ecr3,ecr5,ecr11,ecr15,ecr17,ecr19,ecr22,ecr25,ecr27,ecr29,ecr31,ecr33,ecr35
(1=7) (2=6) (3=5) (4=4) (5=3) (6=2) (7=1) INTO ecr2_R,ecr3_R,ecr5_R,ecr11_R,ecr15_R,ecr17_R,ecr19_R,ecr22_R,
ecr25_R,ecr27_R,ecr29_R,ecr31_R,ecr33_R,ecr35_R.
EXECUTE .

COMPUTE ECRanxiety=mean.16 (ecr2_R, ecr4, ecr6, ecr8, ecr10, ecr12, 
   ecr14,  ecr16, ecr18, ecr20, ecr22_R, ecr24, ecr26, ecr28, 
   ecr30, ecr32, ecr34, ecr36) .
COMPUTE ECRavoidance=mean.16 (ecr1, ecr3_R, ecr5_R, ecr7, ecr9, 
   ecr11_R, ecr13, ecr15_R, ecr17_R, ecr19_R, ecr21, ecr23, 
   ecr25_R, ecr27_R, ecr29_R, ecr31_R, ecr33_R, ecr35_R) .
EXECUTE.

COMPUTE SumValues = SUM.33(ecr1,ecr2_R,ecr3_R,ecr4,ecr5_R,ecr6,ecr7,ecr8,ecr9,ecr10,ecr11_R,ecr12,ecr13,ecr14,
ecr15_R,ecr16,ecr17_R,ecr18,ecr19_R,ecr20,ecr21,ecr22_R,ecr23,ecr24,ecr25_R,ecr26,ecr27_R,ecr28,ecr29_R,ecr30,
ecr31_R,ecr32,ecr33_R,ecr34,ecr35_R,ecr36) .
EXECUTE .

COMPUTE NumMiss = NMISS(ecr1,ecr2_R,ecr3_R,ecr4,ecr5_R,ecr6,ecr7,ecr8,ecr9,ecr10,ecr11_R,ecr12,ecr13,ecr14,
ecr15_R,ecr16,ecr17_R,ecr18,ecr19_R,ecr20,ecr21,ecr22_R,ecr23,ecr24,ecr25_R,ecr26,ecr27_R,ecr28,ecr29_R,ecr30,
ecr31_R,ecr33_R,ecr34,ecr35_R,ecr36) .
EXECUTE .
