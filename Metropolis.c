#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define M 2147483697
#define E 2.718281828459

float metropolis(float x_0, float x_1);
float gaussiana(float x_0, float x_1);
float dist_gauss(float x_0, float paso);

int main(int argc, char *argv[]){

	int item;
	float paso;
	sscanf(argv[1],"%i",&item);
	sscanf(argv[2],"%f",&paso);
	
	if(item){
		srand(time(NULL));
		FILE *file;
		float x_0=0.0;
		char filename[64];
		sprintf(filename, "distribucion_gaussiana_delta_%.2f.txt", paso);
		file=fopen(filename,"a");
	
		for(int i=0; i<100000; i++){
			x_0=dist_gauss(x_0,paso);
			fprintf(file,"%.3f\n", x_0);
		}
		fclose(file);
	}
}

float metropolis(float x_0, float x_1){

	if(rand()*1.0/M<gaussiana(x_0,x_1)){
		x_0=x_1;
	}
		
	return x_0;
}

float gaussiana(float x_0, float x_1){
	
	float a=powl(E,-x_1*x_1/2)/powl(E,-x_0*x_0/2);
	return a;	
}

float dist_gauss(float x_0, float paso){

	float x_1;
	x_1=x_0+paso*2*(rand()*1.0/M-0.5);
	x_0=metropolis(x_0,x_1);
	
	return x_0;
}
