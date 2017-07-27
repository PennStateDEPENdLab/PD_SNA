*code written by Kimberly Nolf.  creates global summary scores for each prototype and calculates theta scores based on the raw score.
* see Pilkonis, Kim,Yu, Morse manuscript, table 2, conversion of raw scores to IRT theta scores.

COMPUTE DEPEND=SUM(proto1_1,proto1_2,proto1_3,proto1_4,proto1_5) .
COMPUTE AMBIV=SUM(proto2_1,proto2_2,proto2_3,proto2_4,proto2_5) .
COMPUTE CAREGIV=SUM(proto3_1,proto3_2,proto3_3,proto3_4,proto3_5) .
COMPUTE RIGIDSC=SUM(proto4_1,proto4_2,proto4_3,proto4_4,proto4_5) .
COMPUTE DEFSEP=SUM(proto5_1,proto5_2,proto5_3,proto5_4,proto5_5) .
COMPUTE EMODETACH=SUM(proto6_1,proto6_2,proto6_3,proto6_4,proto6_5) .
COMPUTE SECURE=SUM(proto7_1,proto7_2,proto7_3,proto7_4,proto7_5) .
EXECUTE .

COMPUTE DEPtheta=0 .
IF (DEPEND=0) DEPtheta= -1.15 .
IF (DEPEND=1) DEPtheta= -0.52 .
IF (DEPEND=2) DEPtheta= -0.19 .
IF (DEPEND=3) DEPtheta= 0.10 .
IF (DEPEND=4) DEPtheta= 0.34 .
IF (DEPEND=5) DEPtheta= 0.58 .
IF (DEPEND=6) DEPtheta= 0.80 .
IF (DEPEND=7) DEPtheta= 1.04 .
IF (DEPEND=8) DEPtheta= 1.31 .
IF (DEPEND=9) DEPtheta= 1.63 .
IF (DEPEND=10) DEPtheta= 2.03 .
EXECUTE .

COMPUTE AMBIVtheta=0 .
IF (AMBIV=0) AMBIVtheta= -1.13 .
IF (AMBIV=1) AMBIVtheta= -0.50 .
IF (AMBIV=2) AMBIVtheta= -0.18 .
IF (AMBIV=3) AMBIVtheta= 0.09 .
IF (AMBIV=4) AMBIVtheta= 0.32 .
IF (AMBIV=5) AMBIVtheta= 0.53 .
IF (AMBIV=6) AMBIVtheta= 0.74 .
IF (AMBIV=7) AMBIVtheta= 0.96 .
IF (AMBIV=8) AMBIVtheta= 1.20 .
IF (AMBIV=9) AMBIVtheta= 1.48 .
IF (AMBIV=10) AMBIVtheta= 1.95 .
EXECUTE .

COMPUTE CAREGIVtheta=0 .
IF (CAREGIV=0) CAREGIVtheta= -0.67 .
IF (CAREGIV=1) CAREGIVtheta= 0.03 .
IF (CAREGIV=2) CAREGIVtheta= 0.44 .
IF (CAREGIV=3) CAREGIVtheta= 0.85 .
IF (CAREGIV=4) CAREGIVtheta= 1.16 .
IF (CAREGIV=5) CAREGIVtheta= 1.44 .
IF (CAREGIV=6) CAREGIVtheta= 1.72 .
IF (CAREGIV=7) CAREGIVtheta= 1.99 .
IF (CAREGIV=8) CAREGIVtheta= 2.25 .
IF (CAREGIV=9) CAREGIVtheta= 2.51 .
IF (CAREGIV=10) CAREGIVtheta= 2.82 .
EXECUTE .

COMPUTE RIGIDSCtheta=0 .
IF (RIGIDSC=0) RIGIDSCtheta= -0.68 .
IF (RIGIDSC=1) RIGIDSCtheta= 0.14 .
IF (RIGIDSC=2) RIGIDSCtheta= 0.53 .
IF (RIGIDSC=3) RIGIDSCtheta= 0.84 .
IF (RIGIDSC=4) RIGIDSCtheta= 1.09 .
IF (RIGIDSC=5) RIGIDSCtheta= 1.33 .
IF (RIGIDSC=6) RIGIDSCtheta= 1.56 .
IF (RIGIDSC=7) RIGIDSCtheta= 1.81 .
IF (RIGIDSC=8) RIGIDSCtheta= 2.08 .
IF (RIGIDSC=9) RIGIDSCtheta= 2.41 .
IF (RIGIDSC=10) RIGIDSCtheta= 2.80 .
EXECUTE .

COMPUTE DEFSEPtheta=0 .
IF (DEFSEP=0) DEFSEPtheta= -0.92 .
IF (DEFSEP=1) DEFSEPtheta= -0.23 .
IF (DEFSEP=2) DEFSEPtheta= 0.16 .
IF (DEFSEP=3) DEFSEPtheta= 0.49 .
IF (DEFSEP=4) DEFSEPtheta= 0.77 .
IF (DEFSEP=5) DEFSEPtheta= 1.04 .
IF (DEFSEP=6) DEFSEPtheta= 1.30 .
IF (DEFSEP=7) DEFSEPtheta= 1.56 .
IF (DEFSEP=8) DEFSEPtheta= 1.83 .
IF (DEFSEP=9) DEFSEPtheta= 2.13 .
IF (DEFSEP=10) DEFSEPtheta= 2.51 .
EXECUTE .

COMPUTE EMODETtheta=0 .
IF (EMODETACH=0) EMODETtheta= -0.95 .
IF (EMODETACH=1) EMODETtheta= -0.27 .
IF (EMODETACH=2) EMODETtheta= 0.10 .
IF (EMODETACH=3) EMODETtheta= 0.44 .
IF (EMODETACH=4) EMODETtheta= 0.75 .
IF (EMODETACH=5) EMODETtheta= 1.01 .
IF (EMODETACH=6) EMODETtheta= 1.24 .
IF (EMODETACH=7) EMODETtheta= 1.52 .
IF (EMODETACH=8) EMODETtheta= 1.84 .
IF (EMODETACH=9) EMODETtheta= 2.13 .
IF (EMODETACH=10) EMODETtheta= 2.47 .
EXECUTE .

COMPUTE SECUREtheta=0 .
IF (SECURE=0) SECUREtheta= -0.99 .
IF (SECURE=1) SECUREtheta= -0.27 .
IF (SECURE=2) SECUREtheta= 0.15 .
IF (SECURE=3) SECUREtheta= 0.48 .
IF (SECURE=4) SECUREtheta= 0.80 .
IF (SECURE=5) SECUREtheta= 1.08 .
IF (SECURE=6) SECUREtheta= 1.33 .
IF (SECURE=7) SECUREtheta= 1.60 .
IF (SECURE=8) SECUREtheta= 1.87 .
IF (SECURE=9) SECUREtheta= 2.15 .
IF (SECURE=10) SECURETtheta= 2.55 .
EXECUTE .

