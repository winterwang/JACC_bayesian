 data d1 ; missing X ; infile "D:\work\20170228\dataset\C**.DATA" ; input 			 


 /*  DATA SET FROM DR LIN*/

/* data d1 ; missing X ; infile "D:\work\20170228\data_copy20190228\C**.DATA" ; input 	 /*DATA SET COPIED FROM CD ROM ON 20190228*/

/*data d1 ; missing X ; infile "D:\work\jacc2009_121220\C**.DATA" ; input */			 /*  DATA SET FROM MANY YEARS AGO VERY OLD DATA SET USED WHEN IN NAGOYA UNIV.*/

 
#1     areano    2-4
		p_HT			36
		p_MI   		37
		p_APO       35
		p_KID        38 
		p_LIV         39
		p_DM			41
		p_ULC		42
          p_can1c   44-46
       p_can1    47
       p_can2c   48-50
       p_can2    51
       p_oth1c   52-54
       p_oth1    55
       p_oth2c   56-58
       p_oth2    59
     
#2     tr_sex    11
       tr_age    78-80
       tr_birth  03-08
	   F_APO   32
	   F_HT     33
	   F_HD     34
	   F_DM    35
	   M_APO  43
	   M_HT    44
	   M_HD    45
	   M_DM    46

      
#3      regy          3-6 
            regm          7-8 
            regd          9-10  
			sport        30
			record  		3-10
			SLEEP   	15-17
			SAMPO     31
			BF_JP  		59
			BF_WST     60
			BF_KAYU    61
			BF_OTH     62
			BF_NON     63
            
#4       T_DX              2 
           J_YEARX       3-6 
           J_MONTHX    7-8 
           J_AGEX          9-11 
           CNT_OTH      12
			MILK 			19
			YOGURT		20
			CHEESE			21
			BUTTER			22
			MARGARIN 	23
			FISH				26
			SPI				29
			FRU				44
			CAKE			45
			COFE			46
			COFEQ			47-48
			COFEM			51
			TEA1				53
			TEA2				54-55
			TEAM			58
			GreTEA1	 	61
			GreTEA2		62-63
			OLOTEA1		64
			OLOTEA2	 	65-66
           
#5	 ICD9  02-06		
ICD10  $   07-11	
DR1 12
DR1AGE 13-14
DR1F 15
ALCOHOL 27-29
SM1 43
SM1AGE 44-45
TOBACCO 46-47
SUIKOMI 48
KINEN 49
OC1 		68
OCSIT	73


#6				    
ht10      12-15
wt10      16-19
wt20      20-23
HOME1 		30
HOME2		31
HOMEO		32
SCHOOL	33-34
MARRY 		43
MAR_AGE	44-45
MARRIED	46
PREGNANC 51-52
DELIVERY 53-54
MEN_AGE 57-58
MENO_AGE 59-60
SEXHORM1 62

#7 	
            actual 76-79
			CHUKAN		2
			DIS1		22
			DIS2		23
			DIS3		24
			CANC	25

	
#8
			FOOD7	14
			FOODC7	15
			FOOD8		16
			FOOD8C	17
			FOOD9		18
			FOOD9C	19
			FOOD10		20
			FOOD10C	21
			FOOD11		22
			FOOD11C	23
 
#9 
			DR1 			2
			DR1CH		6
			SM1			16
			TOBACCO	19-20
			TOBCHAN		21
			DR1_Y			59-61
			SM1_Y			62-64
			SM2_QY			65-67
	
#10    
#11 	
#12  ICD10R1s $ 3-6
s_ymd1s 8-15
ICD10R2s $ 17-20
s_ymd2s  22-29
ICD10R3s $ 31-34
s_ymd3s 36-43
ICD10R4s $ 45-48
s_ymd4s 50-57

#13  ICD10i_1 $ 13-16
ICD2nd1 $ 17-25
s_ymd1 26-33
r_end1 $  64
endpointR1 66-73
actual_R1 75-79

#14 	    ICD10i_2 $ 13-16
ICD2nd2 17-25
s_ymd2 26-33
r_end2 $ 64
endpointR2 66-73
actual_R2 75-79

#15   ICD10i_3 $ 13-16
ICD2nd3 17-25
s_ymd3 26-33
r_end3 $ 64
endpointR3 66-73
actual_R3 75-79

#16   r_age4x    4-6
      ICD10R4x $ 8-11
      ICD10i_4 $ 13-16
      ICD2nd4    17-25
      s_ymd4     26-33
      r_end4  $   64
      endpointR4 66-73
      actual_R4  75-79

#17   
#38
		ENERGY			3-15
		PROTN			55-67
		FAT				68-80
		CHO				3-15
		CALC			55-67

	
#55   
       seizon 72
       endpointD 75-80

#56
#57
#58
#59
#60
#61
#62
#63
;



proc export data = d1
/*					outfile = 'D:\work\20170228\csv\Pool_LowDenSMK.csv'*/
					outfile = 'D:\work\20170228\csv\StrokeMilk.csv'
/*					outfile = 'D:\work\20170228\data_copy20190228\Pool_Breast_Alc_20190228.csv'*/
/*					outfile = "D:\work\jacc2009_121220\Pool_Breast_Alc_20190228.csv"*/
					dbms = csv
					replace; 
run;


DATA D1;
	SET	D1;
	IF tr_age >= 40 & tr_age < 80 THEN OUTPUT; 
RUN; 

PROC FREQ DATA = D1;
	TABLE TR_SEX;
RUN; 
	
PROC MEANS DATA = D1; 
		VAR ACTUAL; 
RUN;

PROC UNIVARIATE DATA = D1 PLOT NORMAL; 
VAR ACTUAL; 
RUN;
