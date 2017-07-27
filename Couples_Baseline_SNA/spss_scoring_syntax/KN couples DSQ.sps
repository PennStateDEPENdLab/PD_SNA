alter type dsq1 to dsq27 (f8.0).
execute.

COMPUTE DSQAlcProb=DSQ15a + DSQ15b + DSQ15c + DSQ15d + DSQ15e + DSQ15f + DSQ15g + DSQ15h + DSQ15i + DSQ15j. 
EXECUTE.

VARIABLE LABELS DSQAlcProb 'DSQ Alcohol Problems: total score for 15a-15j' .
Alter type dsqalcprob (F8.0).

RECODE DSQ1 (SYSMIS=9999).
EXECUTE.

DO IF (DSQ1 = 9999).
RECODE DSQ2a DSQ2b DSQ2c DSQ3 DSQ4 DSQ8 DSQ9 DSQ10 DSQ11 DSQ12 DSQ13 DSQ14 DSQ15a 
    DSQ15b DSQ15c DSQ15d DSQ15e DSQ15f DSQ15g DSQ15h DSQ15i DSQ15j DSQAlcProb DSQ16 DSQ17 DSQ18 DSQ19 
    DSQ20 DSQ21 DSQ22 DSQ23 DSQ24 DSQ25 DSQ26 DSQ27 (MISSING=9999).
END IF.
EXECUTE.

DO IF (DSQ4=0).
RECODE DSQ5 DSQ6 DSQ7 (MISSING=8888).
END IF.
EXECUTE.

DO IF (DSQ1 =0).
RECODE DSQ2a DSQ2b DSQ2c DSQ3 DSQ4 DSQ5 DSQ6 DSQ7 DSQ8 DSQ9 DSQ10 DSQ11 DSQ12 DSQ13 DSQ14 DSQ15a 
    DSQ15b DSQ15c DSQ15d DSQ15e DSQ15f DSQ15g DSQ15h DSQ15i DSQ15j DSQAlcProb DSQ16 DSQ17 DSQ18 DSQ19 
    DSQ20 DSQ21 DSQ22 DSQ23 DSQ24 DSQ25 DSQ26 DSQ27 (MISSING=8888).
END IF.
EXECUTE.

RECODE DSQ2a TO DSQ27 (SYSMIS=9999).
EXECUTE.

VARIABLE LABELS
dsq1 'ever had a drink of alcohol?'
dsq2a 'at what age did you have your first drink?'
dsq2b 'if you never had a drink, do you think that you will ever drink alcohol?'
dsq2c 'if you never had a drink, at what age do you think you will take your first drink?'
dsq3 'during your life, how many times have you had at least one drink of alcohol?'
dsq4 'in the last 30 days, how many times have you had at least one drink?'
dsq5 'in the last 30 days, on the days that you drank, on average, how many drinks did you have?'
dsq6 'in the last 30 days, how many times did you have 5 or more drinks of alcohol at one time?'
dsq7 'what is the largest number of drinks you had on any day in the last 30 days?'
dsq8 'which of the following best describes how often you drink alcohol?'
dsq9 'which of the following best describes how much alcohol you usually drink at one time?'
dsq10 'which of the following is true for you?'
dsq11 'which of the following is true for you?'
dsq12 'who do you usually drink with?'
dsq13 'where do you usually drink alcohol?'
dsq14 'when do you usually drink alcohol?'
dsq15a 'I have gotten a hangover from drinking alcohol'
dsq15b 'I have gotten nauseous and/or vomited from drinking alcohol'
dsq15c 'I have had a blackout while drinking alcohol'
dsq15d 'There have been times when I could not recall what I did while drinking alcohol'
dsq15e 'I have gotten in trouble with my family or significant other for drinking alcohol'
dsq15f 'I have gotten in trouble with school or work for drinking alcohol'
dsq15g 'I have gotten in trouble with friends for drinking alcohol'
dsq15h 'I have gotten into fights while drinking alcohol'
dsq15i 'I have been stopped by police for drunk driving or for being drunk and disorderly'
dsq15j 'I have committed other illegal acts (larceny, robbery, breaking and entering, vandalism, destruction of others property) when drinking alcohol'
dsq16 'in general, from what source do you learn the most about the effects of alcohol'
dsq17 'approximately how much do you spend on alcohol in one week?'
dsq18 'which type of alcoholic drink do you prefer?'
dsq19 'what is the most alcohol you have consumed at one time?'
dsq20 'have you ever been continually drunk for 2 or more days?'
dsq21 'what percent of your friends do you think drank alcohol last month?'
dsq22 'when your friends drink alcohol, on average, how many drinks do you think they have?'
dsq23 'if you have a significant other, do you think they drank in the last month?'
dsq24 'when your significant other drinks alcohol, on average, how many drinks do you think they have?'
dsq25 'over the past year how many times have you tried to cut down/stop drinking alcohol'
dsq26 'for how long did you cut down/stop your drinking?'
dsq27 'how likely is it that you will try to cut down or stop drinking alcohol in the next year?'.
Execute.

VALUE LABELS
dsq1 dsq2b dsq15a TO dsq15j 0 'no' 1 'yes'.
Execute.
VALUE LABELS
Dsq8 1 'never had a drink of alcohol' 2 'have only had 1, 2, 3, or 4 drinks of alcohol in my life'
3	'only drink alcohol 3 or 4 times a year' 4 'I drink alcohol about once a month' 5 'I drink alcohol once or twice a week' 6 'I drink alcohol almost daily'.
Execute.
VALUE LABELS
Dsq9 1 'I do not drink alcohol at all' 2 'small amounts of alcohol (one beer or one drink or less)'
3 'moderate amounts of alcohol (between 2-3 beers or drinks)' 4 'quite a bit of alcohol (between 4-8 beers or drinks)' 5 'a lot of alcohol (more than 9 beers or drinks)'.
Execute.
VALUE LABELS
dsq3 0 'never' 1 '1-2 times' 2 '3-5 times' 3 '6-10 times' 4 '11-50 times' 5 '51-100 times' 6 'over 100 times'.
Execute.
VALUE LABELS
dsq10 1 'I have never been drunk' 2 'I have been drunk once or twice in my life' 3 'I get drunk 2, 3, or 4 times a year'
4 'I get drunk about once a month' 5 'I get drunk about once a week' 6 'I get drunk more than once a week'.
Execute.
VALUE LABELS
DSQ11 1 'I dont drink alcohol' 2 'When I drink alcohol, I always stop before I get drunk' 3 'When I drink alcohol, I almost always stop before I get drink'
4 'When I drink alcohol, I stop before I get drunk more than one-half of the time' 5 'When I drink alcohol, I get drunk more than one-half the time'
6 'When I drink alcohol, I almost always get drunk'.
Execute.
VALUE LABELS
DSQ12 1 'I dont drink alcohol' 2 'Im usually with my family when I drink alcohol' 3 'Im usually with a group of friends when I drink alcohol'
4 'Im usually alone when I drink alcohol' 5 'Im usually alone with my boyfriend/girlfriend when I drink alcohol'.
EXECUTE.
VALUE LABELS
DSQ13 1 'I dont drink alcohol' 2 'I usually drink alcohol at home' 3 'I usually drink alcohol at a friends home' 4 'I usually drink alcohol just before, at, or after a sporting event'
5 'I usually drink alcohol just before, at, or after a party' 6 'I usually drink alcohol at work/school' 7 'I usually drink alcohol in a car' 8 'I usually drink alcohol at a religious service or activity'.
Execute.
VALUE LABELS
DSQ14 1 'I dont drink alcohol' 2 'I usually drink alcohol in the morning, before work/school' 3 'I usually drink alcohol during work/school hours' 
4 'I usually drink alcohol during the day on Saturday or Sunday' 5 'I usually drink alcohol during the week nights (Sun-Thurs)' 6 'I usually drink alcohol during the weekend nights (Fri or Sat)'
7 'I usually drink alcohol every day'.
EXECUTE.
VALUE LABELS
dsq16 1 'my parents' 2 'my peers/friends' 3 'my church' 4 'the mass media (TV, radio, ads, books, magazines, etc.)' 5 'my school' 6 'my own experience with alcohol'
7 'other'.
EXECUTE.
VALUE LABELS
dsq17 1 'nothing, I dont drink alcohol' 2 '1-5 dollars' 3 '5.01-10 dollars' 4 '10.01-15 dollars' 5 '15.01-20 dollars' 6 '20.01-25 dollars' 7 'more than 25 dollars' 8 'I drink but do not pay for it'.
EXECUTE.
VALUE LABELS
dsq18 1 'I dont drink' 2 'beer' 3 'wine' 4 'liquor (including mixed drinks)'.
EXECUTE.
VALUE LABELS
DSQ19 1 'I dont drink' 2 '1-2 drinks/beers' 3 '3-5 drinks/beers' 4 '6-11 drinks/beers' 5 '1 pint of liquor or 12 beers' 6 'between a pint and fifth of liquor or 12-23 beers'
7 'over a fifth of liquor or a case or more of beer'.
EXECUTE.
VALUE LABELS
dsq20 1 'no' 2 'yes, once or twice' 3 'yes, three or more times'.
EXECUTE.
VALUE LABELS
dsq23 0 'no' 1 'yes' 2 'I currently do not have a romantic partner'.
EXECUTE.
VALUE LABELS
DSQ26 1 'never/not appropriate' 2 'less than 1 day' 3 '1 day to 1 week' 4 '1 week to 1 month' 5 '1 to 3 months' 6 'more than 3 months'.
EXECUTE.
VALUE LABELS
dsq27 0 'NA' 1 'definitely wont try' 2 'probably wont try' 3 'not sure' 4 'probably will try' 5 'definitely will try'.
EXECUTE.

