*code written by Kimberly Nolf.  creates global summary scores for each prototype.

COMPUTE DEPEND=SUM(proto1_1,proto1_2,proto1_3,proto1_4,proto1_5) .
COMPUTE AMBIV=SUM(proto2_1,proto2_2,proto2_3,proto2_4,proto2_5) .
COMPUTE CAREGIV=SUM(proto3_1,proto3_2,proto3_3,proto3_4,proto3_5) .
COMPUTE RIGIDSC=SUM(proto4_1,proto4_2,proto4_3,proto4_4,proto4_5) .
COMPUTE DEFSEP=SUM(proto5_1,proto5_2,proto5_3,proto5_4,proto5_5) .
COMPUTE EMODETACH=SUM(proto6_1,proto6_2,proto6_3,proto6_4,proto6_5) .
COMPUTE SECURE=SUM(proto7_1,proto7_2,proto7_3,proto7_4,proto7_5) .
EXECUTE .
