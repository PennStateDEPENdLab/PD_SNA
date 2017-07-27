*code created by kimberly nolf on 1/9/15 to score the emotion ifb PANAS questionnaire.  This questionnaire does not 
*contain enough items to generate scores for the more specific scales, e.g., fear, hostility, guilt, sadness, joviality, self-
*assurance or attentiveness. 
*General Dimension Scales

COMPUTE
PANASNegAff=sum.9(PANAS2,PANAS4,PANAS6,PANAS7,PANAS8,PANAS11,PANAS13,PANAS15,PANAS18,PANAS20) .
EXECUTE .

COMPUTE 
PANASPosAff=sum.9(PANAS1,PANAS3,PANAS5,PANAS9,PANAS10,PANAS12,PANAS14,PANAS16,PANAS17,PANAS19) .
EXECUTE .

VARIABLE LABELS PANASNegAff 'PANAS Negative Affect' PANASPosAff 'PANAS Positive Affect' .
