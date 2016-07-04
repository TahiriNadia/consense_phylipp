#include "phylip.h"
#include "cons.h"


int main(int argc, Char *argv[])
{ 
	char * treesNewick [] ={"(4:0.6660,((1:1.2234,(10:0.8667,(13:0.7613,(2:0.6189,8:0.8176):1.4461):0.6259):0.5966):0.8599,((3:1.6857,7:0.6000):0.7439,(6:1.9101,15:1.5985):1.2384):1.3035):0.7627,(11:1.1495,((5:0.7647,9:0.6086):1.3284,(12:0.8095,(14:0.6788,16:0.8620):1.3349):0.8865):0.8509):1.4008);",
	 "(4:0.6660,((2:1.2234,(10:0.8667,(13:0.7613,(1:0.6189,8:0.8176):1.4461):0.6259):0.5966):0.8599,((3:1.6857,7:0.6000):0.7439,(6:1.9101,15:1.5985):1.2384):1.3035):0.7627,(11:1.1495,((5:0.7647,9:0.6086):1.3284,(12:0.8095,(14:0.6788,16:0.8620):1.3349):0.8865):0.8509):1.4008);",
	 "(13:0.8785,(4:1.4093,((10:0.9596,(1:1.3741,(6:1.0037,15:1.1345):1.6311):0.6104):0.7605,((7:0.7134,14:0.8217):0.8077,(9:0.7243,16:2.2580):0.8752):0.7479):0.7527):0.6755,(11:1.0347,(3:1.0933,(2:0.7588,(8:0.9821,(5:0.6977,12:0.6440):1.3816):1.3430):1.1131):0.8067):1.0069);",
	 "(13:0.8785,(4:1.4093,((10:0.9596,(1:1.3741,(6:1.0037,15:1.1345):1.6311):0.6104):0.7605,((7:0.7134,14:0.8217):0.8077,(9:0.7243,16:2.2580):0.8752):0.7479):0.7527):0.6755,(11:1.0347,(12:1.0933,(2:0.7588,(8:0.9821,(5:0.6977,3:0.6440):1.3816):1.3430):1.1131):0.8067):1.0069);",
	 "((10:0.8447,(6:3.1321,(1:0.9276,(16:1.2873,(11:0.5728,(2:0.9222,15:0.8928):1.8756):0.5969):0.7027):0.8222):0.5830):0.9111,(4:0.5815,((9:0.6427,13:0.9100):0.6754,(14:0.8118,(5:0.8935,7:1.0046):1.6148):1.1266):0.9732):1.1149,(12:1.1573,(3:1.0087,8:0.7018):0.7171):0.9952);",
	 "((10:0.8447,(6:3.1321,(2:0.9276,(16:1.2873,(11:0.5728,(1:0.9222,15:0.8928):1.8756):0.5969):0.7027):0.8222):0.5830):0.9111,(4:0.5815,((9:0.6427,13:0.9100):0.6754,(14:0.8118,(5:0.8935,7:1.0046):1.6148):1.1266):0.9732):1.1149,(12:1.1573,(3:1.0087,8:0.7018):0.7171):0.9952);",
	 "((1:0.8447,(6:3.1321,(10:0.9276,(16:1.2873,(11:0.5728,(2:0.9222,15:0.8928):1.8756):0.5969):0.7027):0.8222):0.5830):0.9111,(4:0.5815,((9:0.6427,13:0.9100):0.6754,(14:0.8118,(5:0.8935,7:1.0046):1.6148):1.1266):0.9732):1.1149,(12:1.1573,(3:1.0087,8:0.7018):0.7171):0.9952);",
	 "(10:0.8316,((2:0.6145,(5:0.6205,6:0.9002):1.2925):1.0065,(16:1.6438,((7:1.2011,(15:0.8131,(1:0.6308,11:0.6147):0.9774):1.0724):0.6522,(12:1.7871,14:1.3660):1.6932):1.4504):1.0451):0.6352,(4:0.7807,(3:1.4190,(9:0.7055,(8:0.7526,13:0.6769):1.0528):0.7920):1.3507):0.6212);",
	 "(10:0.8316,((2:0.6145,(5:0.6205,6:0.9002):1.2925):1.0065,(16:1.6438,((1:1.2011,(15:0.8131,(7:0.6308,11:0.6147):0.9774):1.0724):0.6522,(12:1.7871,14:1.3660):1.6932):1.4504):1.0451):0.6352,(4:0.7807,(3:1.4190,(9:0.7055,(8:0.7526,13:0.6769):1.0528):0.7920):1.3507):0.6212);",
	 "(10:0.8316,((4:0.6145,(5:0.6205,6:0.9002):1.2925):1.0065,(16:1.6438,((7:1.2011,(15:0.8131,(1:0.6308,11:0.6147):0.9774):1.0724):0.6522,(12:1.7871,14:1.3660):1.6932):1.4504):1.0451):0.6352,(2:0.7807,(3:1.4190,(9:0.7055,(8:0.7526,13:0.6769):1.0528):0.7920):1.3507):0.6212);"};
	 
	 char treeConsense[100000];
	 main_consense(treeConsense,treesNewick,10,16);
		
	 printf("DANS LE MAIN TREE CONSENUSE %s\n",treeConsense);
	 return 1;
}
