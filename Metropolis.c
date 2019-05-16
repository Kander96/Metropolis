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
int flipear_spin(int *red, int dim, int s);

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
		float T,p;
		int dim=32,s;
		
		sscanf(argv[2],"%f",&T);
		
		float beta=1.0/T;
		int m=0;
		char filename[64];
		
		sprintf(filename, "ejercicio_2_a_T_%.1f.txt",T);
		
		int *red;
		red = (int*)malloc(dim*dim*sizeof(int));
		
		//for(int i=0; i<dim*dim; i++)
		//	*(red+i)=1;
		
		file=fopen(filename,"w");
		//for(int j=0; j<1000; j++){
			m=dim*dim;
			for(int i=0; i<dim*dim; i++)
				*(red+i)=1;
				//*(red+i)=(rand()%2)*2-1;
			for(int i=0; i<100*dim*dim; i++){
				
				s=rand()%(dim*dim);	
				p=powl(E,-2*beta*B*(*(red+s)));

				if(rand()*1.0/M<p){
					*(red+s)=-*(red+s);
					m+=*(red+s)*2;
					}
				
				if(i%(dim*dim)==0)
					fprintf(file,"%.4f\n", m*1.0/(dim*dim));
				//fprintf(file,"%.4f\n", p);
			//}
			
			//fprintf(file,"\n");
			
		}	
		fclose(file);
		free(red);
	}
	
	if(item==3){
		srand(time(NULL));
		FILE *file;
		int dim=32,a,m,s;
		float J=0.1;
		float T=100.0;
		float beta=1.0/T;
		char filename[64];
			
		int *red;
		red = (int*)malloc(dim*dim*sizeof(int));

		float *p;
		p = (float*)malloc(dim*dim*sizeof(float));
		for(int i=0; i<5; i++)
			*(p+i)=powl(E,-beta*J*4*(i-2));
		
		sprintf(filename, "ejercicio_2_b.txt");
		file=fopen(filename,"w");
		
		for(int i=0; i<dim*dim; i++)
			*(red+i)=1;
			//*(red+i)=(rand()%2)*2-1;
		m=dim*dim;
		
		for(int i=0; i<5000; i++){
			
			s=rand()%(dim*dim);
			a=flipear_spin(red,dim,s);
			if(rand()*1.0/M<*(p+a)){
				*(red+s)=-*(red+s);
				m+=*(red+s)*2;
			}
		}
		
		
		for(int i=0; i<1000*dim*dim; i++){
			
			s=rand()%(dim*dim);
			a=flipear_spin(red,dim,s);
			if(rand()*1.0/M<*(p+a)){
				*(red+s)=-*(red+s);
				m+=*(red+s)*2;
			}
			if(i%(dim*dim)==0)
			fprintf(file,"%.3f\n",m*1.0/(dim*dim));
		}
		
		free(red);
		free(p);
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

int flipear_spin(int *red, int dim, int a){
	
	int s_u,s_d,s_l,s_r;
	int b;
	int i=a/dim;
	int j=a%dim;

	s_d=*(red+((i+1+dim)%dim)*dim+j);
	s_u=*(red+((i-1+dim)%dim)*dim+j);
	s_r=*(red+i*dim+(j+1+dim)%dim);
	s_l=*(red+i*dim+(j-1+dim)%dim);

/*	
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
*/	
	b=*(red+a)*(s_u+s_d+s_l+s_r);
	b=b/2+2;
	
	return b;		
}
