#ifndef _BOUNDARIES_H
#define _BOUNDARIES_H
void boundaries_check();

struct ghostbox{
	double shiftx;
	double shifty;
	double shiftz;
	double shiftvx;
	double shiftvy;
	double shiftvz;
};

struct ghostbox get_ghostbox(int i, int j, int k);

extern int nghostx;
extern int nghosty;
extern int nghostz;

#endif
