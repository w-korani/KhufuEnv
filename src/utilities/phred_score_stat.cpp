#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdexcept>
#include <stdlib.h>
#include <math.h>

int main(int argc, const char* argv[])
{
	try
	{
		std::string phred = "";
		int Phred;
		std::string y;
		for (int x=0;x<argc;x++)
		{
			y=argv[x];
			if(y=="-p")
			{
				phred = argv[x+1];
				Phred = atoi(phred.c_str());
			}
		}
		//#################################
		if(phred == ""){Phred = 33;}
		//#################################
		int count1 = 1;
		int sum = 0;
		int score=0;
		float acculatescore=0;
		float Pscore=0;
		float acculatePscore=0;
		float avrQ=0;
		float avrP=0;
		float stdQ=0;
		float stdP=0;
		std::string l1,l2,l3,l4;
		std::cout << "id\tQaverage\tQstd" << std::endl;
		for (std::string line; std::getline(std::cin, line);)
		{
			if(count1==1){l1=line;count1+=1;}
			else if (count1==2){l2=line;count1+=1;}
			else if (count1==3){l3=line;count1+=1;}
			else if (count1==4)
			{
				l4=line;
				count1=1;
				for (int i=0;i<l2.length();++i)
				{
					score+=((int)l4[i]-Phred);
				//	Pscore+=(pow(10, -(float(score)/10) ));
					sum+=1;
				}
				avrQ=float(score)/float(sum);
				// avrP=float(Pscore)/float(sum);
				score=0;Pscore=0;
				///////////////
				for (int i=0;i<l2.length();++i)
				{
					score=((int)l4[i]-Phred);
				//	Pscore=(pow(10, -(float(score)/10) ));
					acculatescore+=( score - avrQ );
				//	acculatePscore+=( Pscore - avrP );
				}
				stdQ = sqrt(abs(acculatescore/sum));
				// stdP = sqrt(abs(acculatePscore/sum));
				std::cout << l1.erase(0, 1) << "\t" << avrQ << "\t" << stdQ << std::endl; // << "\t" << avrP << "\t" << stdP << std::endl;
				sum=0;score=0;Pscore=0;acculatescore=0;acculatePscore=0;
			}
		}
		//#################################
	}
	catch (std::exception e)
	{
		std::cout<<e.what()<<std::endl;
	}
	return 0;
};
// rm phred_score_stat.cpp ; nano phred_score_stat.cpp; g++ phred_score_stat.cpp -o phred_score_stat
