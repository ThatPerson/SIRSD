#include <stdio.h>
#include <math.h>
#include <stdlib.h>

struct System { 
	double S, I, R;
	double b, g;
	double t;
};

double iterate(struct System *s, double dt) {
	double dS, dI, dR;
	double N = s->S + s->I + s->R;
	dS = - (s->b * s->I * s->S) / N;
	dI = ((s->b * s->I * s->S) / N) - (s->g * s->I);
	dR = s->g * s->I;
	s->S += dS * dt;
	s->I += dI * dt;
	s->R += dR * dt;
	s->t += dt;
	return s->t;
}

int setup(struct System *s) {
	if (s == NULL)
		return 1;
	s->S = 500;
	s->I = 1;
	s->R = 0;
	s->b = 5;
	s->g = 2;
	s->t = 0;
	return 0;
}

int main(int argc, char * argv[]) {
	struct System s;
	if (setup(&s) == 1) { fprintf(stderr, "Fatal error\n"); return 1; }
	
	FILE *fp = NULL;
	fp = fopen("out.csv", "w");
	if (fp == NULL) { fprintf(stderr, "Fatal error\n"); return 1; }
	int i;
	double dt = 0.001;
	for (i = 0; i < 5000; i++) {
		iterate(&s, dt);
		fprintf(fp, "%d, %f, %f, %f, %f\n", i, s.t, s.S, s.I, s.R);
	}
	fclose(fp);
	return 0;
}
