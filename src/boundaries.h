#ifndef _BOUNDARIES_H
#define _BOUNDARIES_H
void check_boundaries();

struct ghostbox{
	double shiftx;
	double shifty;
	double shiftz;
} ghostbox;

struct ghostbox get_ghostbox(int i, int j, int k);

extern int nghostx;
extern int nghosty;
extern int nghostz;

#endif
