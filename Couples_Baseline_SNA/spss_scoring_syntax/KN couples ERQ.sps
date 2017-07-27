*ERQ code written by kimberly nolf on 3.17.14

COMPUTE ERQreapp = sum.5(ERQ1,ERQ3,ERQ5,ERQ7,ERQ8,ERQ10) .
EXECUTE .

COMPUTE ERQsup = sum.3(ERQ2,ERQ4,ERQ6,ERQ9) .
EXECUTE .

RECODE ERQreapp ERQSup (SYSMIS=9999) .
EXECUTE .
VARIABLE LABELS ERQreapp 'ERQ Reappraisal' ERQsup 'ERQ Suppression' .
EXECUTE .
