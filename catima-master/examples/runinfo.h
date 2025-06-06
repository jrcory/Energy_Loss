#ifndef runinfo_H
#define runinfo_H
namespace info{
	class runinfo{
		public:
		int runnb,A,Z,gasel;
		int gasA[4]={0};
		int gasZ[4]={0};
		int gasS[4]={0};
		float Ein, density;
		void setinfo(int nb);
	};
}

#endif /* runinfo_H */
