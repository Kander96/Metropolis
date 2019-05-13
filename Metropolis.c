#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define M 2147483697
#define E 2.718281828459
#define K 1.38064852*powl(10,-23)

float metropolis(float x_0, float x_1);
float gaussiana(float x_0, float x_1);
float dist_gauss(float x_0, float paso);

int main(int argc, char *argv[]){

	int item;
	sscanf(argv[1],"%i",&item);
	
	if(item==1){
		srand(time(NULL));
		FILE *file;
		float paso;
		float x_0=0.0;
		char filename[64];
		sscanf(argv[2],"%f",&paso);
		sprintf(filename, "distribucion_gaussiana_delta_%.2f.txt",paso);
		file=fopen(filename,"a");
	
		for(int i=0; i<100000; i++){
			x_0=dist_gauss(x_0,paso);
			fprintf(file,"%.3f\n", x_0);
		}
		fclose(file);
	}
	
	if(item==2){
		srand(time(NULL));
		
		FILE *file;
		float B=-100.0;
		float T=6000.0,p;
		int dim=32,s;
		float beta=1.0/T;
		int m=0;
		
		int *red;
		red = (int*)malloc(dim*dim*sizeof(int));
		
		//for(int i=0; i<dim*dim; i++)
		//	*(red+i)=1;
		
		file=fopen("Ejercicio_2_a.txt","w");
		for(int j=0; j<1000; j++){
			m+=dim*dim;
			for(int i=0; i<10000; i++){
				for(int i=0; i<dim*dim; i++)
					*(red+i)=1;
				
				s=rand()%(dim*dim);	
				p=powl(E,-2*beta*B*(*(red+s)));

				if(rand()*1.0/M<p){
					*(red+s)=-*(red+s);
					m+=*(red+s)*2;
					}
				
				fprintf(file,"%.4f\t", m*1.0/(dim*dim));
				//fprintf(file,"%.4f\n", p);
			}
			
			fprintf(file,"\n");
			
		}	
		fclose(file);
		free(red);
	}
	
	if(item==3){
		srand(time(NULL));
		int dim=32;
		
		int *red;
		red = (int*)malloc(dim*dim*sizeof(int));
		
		for(int i=0; i<dim*dim; i++)
			*(red+i)=1;
		
		
	
	}
	
	if(item==100){
		printf("%i\n",563%100);
	
	
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

/*int flipear_spin(int *red, int dim){
	
	int s_u,s_d,s_l,s_r;
	int a=rand()%(dim*dim);
	
	if(a/dim==0)
		s_u=*(red+a+dim*dim-dim);
	else
		s_u=*(red+a-dim);
	
	if(a/dim==dim-1)
		s_d=*(red+a-dim*dim+dim);
	else
		s_d=*(red+a+dim);
		
	if(a%dim==0)
		s_l=*(red+a-1+dim);
	else
		s_l=*(red+a-1);
		
	if(a%dim==dim-1)
		s_r=*(red+a+1-dim);
	else
		s_r=*(red+a+1);
		
	return 0;		
}*/
