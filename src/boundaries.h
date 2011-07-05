#ifndef _BOUNDARIES_H
#define _BOUNDARIES_H
void check_boundaries();

struct ghostbox{
	double shiftx;
	double shifty;
	double shiftz;
} ghostbox;

struct ghostbox get_ghostbox(int i, int j, int k);

extern const int nghostx;
extern const int nghosty;
extern const int nghostz;

#endif
