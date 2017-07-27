*KN code to create variable for underrepresented groups for the NIH targetted enrollment report

RECODE race ('White'=0) ('Asian;White'=0) ('Asian'=0) (SYSMIS=0) (ELSE=1) INTO under_represented. 
VARIABLE LABELS  under_rep 'under_rep'. 
EXECUTE.