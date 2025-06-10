#ifndef runinfo_H
#define runinfo_H
#include <string>
namespace info{
	class runinfo{
		public:
		int runnb,A,Z,gasel;
		int gasA[4]={0};
		int gasZ[4]={0};
		int gasS[4]={0};
		float Ein, density;
		char* lab; 
		void setinfo(int nb);
	};
}

#endif /* runinfo_H */
