* SBQ-R scoring written by Kimberly Nolf on 3/24/14 for scoring of measure in the couples study

ALTER TYPE sbq1 TO SBQ4_3MTH (f8.0).
execute.

RECODE SBQ1, SBQ1_3MTH (1=1) (2=2) (3=3) (4=3) (5=4) (6=4) INTO SBQ1r, SBQ1_3MTHr .
RECODE SBQ3, SBQ3_3MTH (1=1) (2=2) (3=2) (4=3) (5=3) INTO SBQ3r, SBQ3_3MTHr .
RECODE SBQ4, SBQ4_3MTH (1=0) (2=1) (3=2) (4=3) (5=4) (6=5) (7=6) INTO SBQ4r, SBQ4_3MTHr.
EXECUTE .

COMPUTE SBQLife=SUM(SBQ1r, SBQ2, SBQ3r, SBQ4r) .
EXECUTE .
COMPUTE SBQ3mo=SUM(SBQ1_3MTHr, SBQ2_3MTH, SBQ3_3MTHr, SBQ4_3MTHr) .
EXECUTE .

VARIABLE LABELS
SBQLife 'SBQ Lifetime score' SBQ3mo 'SBQ past 3 months score' .
