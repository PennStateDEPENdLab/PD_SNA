*Code written by kimberly nolf for the couples study on 3.17.14

COMPUTE FStotal=SUM.7(FS1,FS2,FS3,FS4,FS5,FS6,FS7,FS8) .
EXECUTE .

RECODE FStotal (SYSMIS=9999) .
EXECUTE .

VARIABLE LABELS FStotal 'Flourishing Scale Total' .

alter type fstotal (f8.0).
Execute.
