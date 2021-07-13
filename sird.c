#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#define MAP_SIZE (100)
#define POS(x, y) ((x * MAP_SIZE) + y) 

struct System { 
	double *S, *I, *R;
	double tS, tI, tR;
	double *lS, *lI, *lR;
	double b, g, w;
	double t;
	double *DS, *DI, *DR;
};

struct Vec2D {
	double x, y;
};


int laplacianmap(double *M, double *oup) {
	struct Vec2D xa, ya;
	int x, y;
	double d;
	double dx2, dy2;
	int rx,ry,ox,oy;
	for (x = 0; x < MAP_SIZE; x++) {
		for (y = 0; y < MAP_SIZE; y++) {
			rx = x-1; ry = y-1;
			ox = x+1; oy = y+1;
			if (rx < 0) rx = MAP_SIZE - 1;
			if (ry < 0) ry = MAP_SIZE - 1;
			if (ox > MAP_SIZE - 1) ox = 0;
			if (oy > MAP_SIZE - 1) oy = 0;
			
			dx2 = (M[POS(ox, y)] - M[POS(x, y)]) + (M[POS(rx, y)] - M[POS(x, y)]);
			dy2 = (M[POS(x, oy)] - M[POS(x, y)]) + (M[POS(x, ry)] - M[POS(x, y)]);
		
			oup[POS(x, y)] = dx2 + dy2;
		
		//	nabla(inp_x, x, y, &xa);
		//	nabla(inp_y, x, y, &ya);
		//	oup[POS(x,y)] = inp_x[POS(x,y)] + inp_y[POS(x, y)];
		//	oup[POS(x, y)] = xa.x + ya.y;
		}
	}
	return 1;
}
			

double iterate(struct System *s, double dt) {
	double dS, dI, dR;
	double N;
	
	laplacianmap(s->S, s->lS);
	laplacianmap(s->I, s->lI);
	laplacianmap(s->R, s->lR);
	
	int x, y, cp;
	for (x = 0; x < MAP_SIZE; x++) {
		for (y = 0; y < MAP_SIZE; y++) {
			cp = POS(x,y);
			N = s->S[cp] + s->I[cp] + s->R[cp];
			dS = s->DS[cp] * s->lS[cp];
			dI = s->DI[cp] * s->lI[cp];
			dR = s->DR[cp] * s->lR[cp];
			dS += - (s->b * s->I[cp] * s->S[cp]) / N;
			dI += ((s->b * s->I[cp] * s->S[cp]) / N) - (s->g * s->I[cp]);
			dR += s->g * s->I[cp];
			dS += s->w * s->R[cp];
			dR += -(s->w * s->R[cp]);
			s->S[cp] += dS * dt;
			s->I[cp] += dI * dt;
			s->R[cp] += dR * dt;
			s->tS += s->S[cp];
			s->tI += s->I[cp];
			s->tR += s->R[cp];
		}
	}
	s->t += dt;
	return s->t;
}

int setup(struct System *s) {
	if (s == NULL)
		return 1;
	s->S = (double *) malloc(sizeof(double) * MAP_SIZE * MAP_SIZE);
	s->I = (double *) malloc(sizeof(double) * MAP_SIZE * MAP_SIZE);
	s->R = (double *) malloc(sizeof(double) * MAP_SIZE * MAP_SIZE);
	s->lS = (double *) malloc(sizeof(double) * MAP_SIZE * MAP_SIZE);
	s->lI = (double *) malloc(sizeof(double) * MAP_SIZE * MAP_SIZE);
	s->lR = (double *) malloc(sizeof(double) * MAP_SIZE * MAP_SIZE);
	s->DS = (double *) malloc(sizeof(double) * MAP_SIZE * MAP_SIZE);
	s->DI = (double *) malloc(sizeof(double) * MAP_SIZE * MAP_SIZE);
	s->DR = (double *) malloc(sizeof(double) * MAP_SIZE * MAP_SIZE);
	int x, y;
	for (x = 0; x < MAP_SIZE * MAP_SIZE; x++) {
		s->S[x] = 500;
		s->I[x] = 0;
		s->R[x] = 0;
		s->lS[x] = 0;
		s->lI[x] = 0;
		s->lR[x] = 0;
		s->DS[x] = 3;
		s->DI[x] = 3;
		s->DR[x] = 3;
	}
	for (x = 0; x < 3; x++) 
		s->I[rand() % (MAP_SIZE * MAP_SIZE)] = 1;
	s->b = 5;
	s->g = 2;
	s->w = 0.2;
	s->t = 0;
	for (x = 10; x < 50; x++) {
		for (y = 40; y < 80; y++) {
			s->DS[POS(x, y)] = 0.2;
			s->DI[POS(x, y)] = 0.2;
			s->DR[POS(x, y)] = 0.2;
		}
	}
	
	for (x = 40; x < 80; x++) {
		for (y = 10; y < 30; y++) {
			s->DS[POS(x, y)] = 6;
			s->DI[POS(x, y)] = 6;
			s->DR[POS(x, y)] = 6;
		}
	}
			
	return 0;
}

int freemap(struct System *s) {
	if (s == NULL) return 1;
	free(s->S);
	free(s->I);
	free(s->R);
	free(s->lS);
	free(s->lI);
	free(s->lR);
	free(s->DS);
	free(s->DI);
	free(s->DR);
	return 0;
}

int printmap(double *M, FILE *fp) {
	int x, y, cp;
	for (x = 0; x < MAP_SIZE; x++) {
		for (y = 0; y < MAP_SIZE; y++) {
			cp = POS(x, y);
			fprintf(fp, "%f\t", M[cp]);
		}
		fprintf(fp, "\n");
	}
	return 0;
}
		
int printsys(struct System *s, int i) {
	char fn[255];
	FILE *fp = NULL;
	sprintf(fn, "maps/S%05d.dat", i);
	fp = fopen(fn, "w");
	if (fp == NULL) return 1;
	printmap(s->S, fp);
	fclose(fp);
	fp = NULL;
	
	sprintf(fn, "maps/I%05d.dat", i);
	fp = fopen(fn, "w");
	if (fp == NULL) return 1;
	printmap(s->I, fp);
	fclose(fp);
	fp = NULL;
	
	sprintf(fn, "maps/R%05d.dat", i);
	fp = fopen(fn, "w");
	if (fp == NULL) return 1;
	printmap(s->R, fp);
	fclose(fp);
	fp = NULL;
	
	return 0;
}
	

int main(int argc, char * argv[]) {
	srand(time(NULL));
	struct System s;
	if (setup(&s) == 1) { fprintf(stderr, "Fatal error\n"); return 1; }

	FILE *fp = NULL;
	fp = fopen("out2.csv", "w");
	if (fp == NULL) { fprintf(stderr, "Fatal error\n"); return 1; }
	int i;
	double dt = 0.01;
	int cp = POS(3, 4);
	char filename[255];
	FILE *lp = fopen("laplacian.csv", "w");
	for (i = 0; i < 5000; i++) {
		s.tS = 0;
		s.tI = 0;
		s.tR = 0;
		iterate(&s, dt);
		if (i%10 == 0) printsys(&s, i);
		//printmap(s.lS, lp);
		fprintf(fp, "%d, %f, %f, %f, %f\n", i, s.t, s.tS, s.tI, s.tR);
	}
	freemap(&s);
	fclose(fp);
	fclose(lp);
	return 0;
}
