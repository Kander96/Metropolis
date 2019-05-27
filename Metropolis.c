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
float medir_a(int *red, int dim, float B_T, FILE *file, float m, int l);
int flipear_spin(int *red, int dim, int s);
float medir_b(int *red, int dim, float J_T, FILE *file, float m, int l);
int print_file(FILE *file);
int variar_dimension(int dim);
int flippear_spin_antiferro(int *red, int dim, int a, float B_T);
int definir_valores_de_p(float *p, float B, float J);
float termalizacion_antiferro(int *red, float *p, int dim, float B_T, float J_T, FILE *file, float m, int l);
float medir_d(int *red, float *p, int dim, float B_T, float J_T, FILE *file, float m, int l);
int flippear_spin_frustacion(int *red, int dim, int a, float B);
int definir_valores_de_p_2(float *p, float B, float J);
float medir_e(int *red, float *p, int dim, float B_T, float J_T, FILE *file, float m, int l);


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
		float B_T=-2.0; //esto es B*=B/T
		int dim=32;
		int l=1000;
		
		//sscanf(argv[2],"%f",&T);
		
		//float beta=1.0/T;
		int m=0;
		char filename[64];
		
		sprintf(filename, "ejercicio_2_a.txt");
		
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
		
		
		for(B_T=-4; B_T<-2; B_T+=0.5)
		m=medir_a(red,dim,B_T,file,m,l);
		
		for(B_T=-2; B_T<-1; B_T+=0.25)
		m=medir_a(red,dim,B_T,file,m,l);
		
		for(B_T=-1; B_T<1; B_T+=0.1)
		m=medir_a(red,dim,B_T,file,m,l);
		
		for(B_T=1; B_T<2; B_T+=0.25)
		m=medir_a(red,dim,B_T,file,m,l);
		
		for(B_T=2; B_T<=4; B_T+=0.5)
		m=medir_a(red,dim,B_T,file,m,l);
		//for(B_T=-4; B_T<-0.05; B_T=B_T/2)
		//m=medir_a(red,dim,B_T,file,m,l);
		
			
	
		fclose(file);
		free(red);
	}
	
	if(item==3){
		srand(time(NULL));
		FILE *file;
		int dim=32,m;
			
		int *red;
		red = (int*)malloc(dim*dim*sizeof(int));
		
		file=fopen("ejercicio_2_b.txt","w");
		
		for(int i=0; i<dim*dim; i++)
			*(red+i)=1;
			//*(red+i)=(rand()%2)*2-1;
		m=dim*dim;
		
		//m=medir_b(red,dim,0.1,file,m,l);				
				
		for(float i=0.1; i<2.0; i+=0.1)
			m=medir_b(red,dim,1.0/i,file,m,1000);		
		
		for(float i=2.0; i<2.8; i+=0.01)
			m=medir_b(red,dim,1.0/i,file,m,1000);
		
		for(float i=2.8; i<4.0; i+=0.05)
			m=medir_b(red,dim,1.0/i,file,m,1000);
		
		for(float i=4.0; i<=6.0; i+=0.1)
			m=medir_b(red,dim,1.0/i,file,m,1000);
		
		fclose(file);		
		free(red);
		
	}
	
	if(item==4){
		srand(time(NULL));
		int dim;
		
		for(dim=4; dim<=64; dim+=dim)
			variar_dimension(dim);
		
		
	}
	
	if(item==5){
		srand(time(NULL));
		FILE *file;
		int dim=32,m;
		float J_T,B_T;
		char filename[64];
		
		sprintf(filename,"ejercicio_2_d.txt");
		file=fopen(filename,"w");
		
		int *red;
		float *p;
		red = (int*)malloc(dim*dim*sizeof(int));
		p = (float*)malloc(10*sizeof(float));
		
		for(int i=0; i<dim*dim; i++)
			*(red+i)=1;
			//*(red+i)=(rand()%2)*2-1;
		m=dim*dim;
		
		for(B_T=-3.0; B_T<=3.0; B_T+=0.5){
			for(float i=0.1; i<=6.1; i+=0.5){
				J_T=1.0/i;
				definir_valores_de_p(p,fabs(B_T),J_T);
				m=termalizacion_antiferro(red,p,dim,B_T,J_T,file,m,100);
			}
		}
		
		//printf("%i\n",m);
		
		fclose(file);
		free(red);
		free(p);
	}
	
	if(item==6){
		srand(time(NULL));
		int dim=32,m;
		float J_T,B_T;
		char filename[64];
		int *red;
		float *p;
		red = (int*)malloc(dim*dim*sizeof(int));
		p = (float*)malloc(10*sizeof(float));
		
		for(int i=0; i<dim*dim; i++)
			*(red+i)=1;
			//*(red+i)=(rand()%2)*2-1;
		m=dim*dim;
		
		for(B_T=-9.0; B_T<=-1.0; B_T=B_T/3.0){
			FILE *file;
			sprintf(filename,"ejercicio_2_d_B_%.2f.txt", B_T);
			file=fopen(filename,"w");
		
			for(float i=0.1; i<=6.0; i+=0.1){
				J_T=1.0/i;
				m=medir_d(red,p,dim,B_T,J_T,file,m,1000);
			}
			
			fclose(file);
		}
		
		B_T=0.0;
		FILE *file;
		sprintf(filename,"ejercicio_2_d_B_%.2f.txt", B_T);
		file=fopen(filename,"w");
	
		for(float i=0.1; i<=6.0; i+=0.1){
			J_T=1.0/i;
			m=medir_d(red,p,dim,B_T,J_T,file,m,1000);
		}
		
		fclose(file);
		
		for(B_T=9.0; B_T>=1.0; B_T=B_T/3.0){
			FILE *file;
			sprintf(filename,"ejercicio_2_d_B_%.2f.txt", B_T);
			file=fopen(filename,"w");
		
			for(float i=0.1; i<=6.0; i+=0.1){
				J_T=1.0/i;
				m=medir_d(red,p,dim,B_T,J_T,file,m,1000);
			}
			
			fclose(file);
		}
		
				
		free(red);
		free(p);
	}
	
	if(item==7){
		srand(time(NULL));
		int dim=32,m;
		float J_T,B_T;
		char filename[64];
		int *red;
		float *p;
		red = (int*)malloc(dim*dim*sizeof(int));
		p = (float*)malloc(18*sizeof(float));
		
		for(int i=0; i<dim*dim; i++)
			*(red+i)=1;
			//*(red+i)=(rand()%2)*2-1;
		m=dim*dim;
		
		for(B_T=1.0; B_T>=-1.0; B_T-=0.5){
			FILE *file;
			sprintf(filename,"ejercicio_2_e_B_%.2f.txt", B_T);
			file=fopen(filename,"w");
		
			for(float i=0.1; i<10.0; i+=0.1){
				J_T=1.0/i;
				m=medir_e(red,p,dim,B_T,J_T,file,m,1000);
			}
			for(int i=10; i<=20; i++){
				J_T=1.0/i;
				m=medir_e(red,p,dim,B_T,J_T,file,m,1000);
			}
			
			fclose(file);
		}
	
	
		free(red);
		free(p);
	}
	
	if(item==100){
		/*FILE *file;
		file=fopen("Hello World","w");
		for(float i=2.0; i>0.01; i=i/2.0)
			print_file(file);
		fclose(file);
		*/
		float B=-97;
		int b;
		int s=-1;
		int s_u=1;
		int s_d=1;
		int s_l=1;
		int s_r=1;
		
		if(B)
			b=s*(s_u+s_d+s_l+s_r)+4+(B*s*1.0/fabs(B)+1)/2;
		else
			b=s*(s_u+s_d+s_l+s_r)+4;
			
		printf("%d\n",b);
	
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

float medir_a(int *red, int dim, float B_T, FILE *file, float m, int l){
	
	int s;
	float p;
	
	fprintf(file,"%.3f\t",B_T);
	for(int i=0; i<5000; i++){
		s=rand()%(dim*dim);	
		p=powl(E,-2*B_T*(*(red+s)));

		if(rand()*1.0/M<p){
			*(red+s)=-*(red+s);
			m+=*(red+s)*2;
			}			
	}

	for(int i=0; i<l*dim*dim; i++){
	
		s=rand()%(dim*dim);	
		p=powl(E,-2*B_T*(*(red+s)));

		if(rand()*1.0/M<p){
			*(red+s)=-*(red+s);
			m+=*(red+s)*2;
			}
	
		if(i%(dim*dim)==0)
			fprintf(file,"%.4f\t", m*1.0/(dim*dim));
	
	}
	fprintf(file,"\n");
	
	return m;
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

int print_file(FILE *file){
	fprintf(file,"Hello World!\n");

	return 0;
}

float medir_b(int *red, int dim, float J_T, FILE *file, float m, int l){
	
	int s,a;
	
	float *p;
	p = (float*)malloc(dim*dim*sizeof(float));
	for(int i=0; i<5; i++)
		*(p+i)=powl(E,-J_T*4*(i-2));		
	
	fprintf(file,"%f\t",J_T);
	
	for(int i=0; i<5000; i++){
		
		s=rand()%(dim*dim);
		a=flipear_spin(red,dim,s);
		if(rand()*1.0/M<*(p+a)){
			*(red+s)=-*(red+s);
			m+=*(red+s)*2;
		}
	}
	
	for(int j=0; j<l; j++){
	
		fprintf(file,"%.3f\t",m*1.0/(dim*dim));
	
		for(int i=0; i<dim*dim; i++){
		
			s=rand()%(dim*dim);
			a=flipear_spin(red,dim,s);
			if(rand()*1.0/M<*(p+a)){
				*(red+s)=-*(red+s);
				m+=*(red+s)*2;
			}
			//if(i%(dim*dim)==0)
			//fprintf(file,"%.3f\n",m*1.0/(dim*dim));
		}
	}
	fprintf(file,"\n");
	
	free(p);
	
	return m;
}

int variar_dimension(int dim){

	FILE *file;
	int m;
	int *red;
	red = (int*)malloc(dim*dim*sizeof(int));
	
	char filename[64];
	sprintf(filename, "ejercicio_2_c_dim_%i.txt",dim);
	file=fopen(filename,"w");
	
	for(int i=0; i<dim*dim; i++)
		*(red+i)=1;
		//*(red+i)=(rand()%2)*2-1;
	m=dim*dim;
	
	//m=medir_b(red,dim,0.1,file,m,l);				
			
	for(float i=0.1; i<2.0; i+=0.1)
		m=medir_b(red,dim,1.0/i,file,m,1000);		
	
	for(float i=2.0; i<2.8; i+=0.01)
		m=medir_b(red,dim,1.0/i,file,m,1000);
	
	for(float i=2.8; i<4.0; i+=0.05)
		m=medir_b(red,dim,1.0/i,file,m,1000);
	
	for(float i=4.0; i<=10.0; i+=0.1)
		m=medir_b(red,dim,1.0/i,file,m,1000);
	
	fclose(file);		
	free(red);
	
	return 0;
}

int flippear_spin_antiferro(int *red, int dim, int a, float B){
	
	int s_u,s_d,s_l,s_r,s=-*(red+a);
	int b;
	int i=a/dim;
	int j=a%dim;

	s_d=*(red+((i+1+dim)%dim)*dim+j);
	s_u=*(red+((i-1+dim)%dim)*dim+j);
	s_r=*(red+i*dim+(j+1+dim)%dim);
	s_l=*(red+i*dim+(j-1+dim)%dim);
	
	if(B)
		b=s*(s_u+s_d+s_l+s_r)+4+(B/fabs(B)*s+1)/2;
	else
		b=s*(s_u+s_d+s_l+s_r)+4;
	
	return b;
}

int definir_valores_de_p(float *p, float B, float J){
	
	*p=powl(E,-(-8*J+2*B));
	*(p+1)=powl(E,-(-8*J-2*B));
	*(p+2)=powl(E,-(-4*J+2*B));
	*(p+3)=powl(E,-(-4*J-2*B));
	*(p+4)=powl(E,-(2*B));
	*(p+5)=powl(E,-(-2*B));
	*(p+6)=powl(E,-(4*J+2*B));
	*(p+7)=powl(E,-(4*J-2*B));
	*(p+8)=powl(E,-(8*J+2*B));
	*(p+9)=powl(E,-(8*J-2*B));
	
	return 0;
}

float termalizacion_antiferro(int *red, float *p, int dim, float B_T, float J_T, FILE *file, float m, int l){
	
	int s;
	int a;
	
	fprintf(file,"%.3f\t",B_T);
	fprintf(file,"%.3f\t",J_T);
	
	for(int i=0; i<l*dim*dim; i++){
		s=rand()%(dim*dim);	
		a=flippear_spin_antiferro(red,dim,s,B_T);
		
		if(rand()*1.0/M<*(p+a)){
			*(red+s)=-*(red+s);
			m+=*(red+s)*2;
			}
		
		fprintf(file,"%.4f\t", m*1.0/(dim*dim));			
	}

	fprintf(file,"\n");
	
	return m;
}

float medir_d(int *red, float *p, int dim, float B_T, float J_T, FILE *file, float m, int l){

	int s;
	int a;
	
	definir_valores_de_p(p,fabs(B_T),J_T);
	
	fprintf(file,"%.3f\t",B_T);
	fprintf(file,"%.3f\t",J_T);
	
	for(int i=0; i<7000; i++){
		
		s=rand()%(dim*dim);	
		a=flippear_spin_antiferro(red,dim,s,B_T);
		
		if(rand()*1.0/M<*(p+a)){
			*(red+s)=-*(red+s);
			m+=*(red+s)*2;
		}
	
	}
	
	for(int j=0; j<l; j++){
		
		fprintf(file,"%.4f\t", m*1.0/(dim*dim));
		
		for(int i=0; i<dim*dim; i++){
			s=rand()%(dim*dim);	
			a=flippear_spin_antiferro(red,dim,s,B_T);
			
			if(rand()*1.0/M<*(p+a)){
				*(red+s)=-*(red+s);
				m+=*(red+s)*2;
			}
						
		}
	}
	fprintf(file,"\n");	

	return m;
}

int flippear_spin_frustacion(int *red, int dim, int a, float B){
	
	int s_u,s_d,s_l,s_r,s_ul,s_ur,s_dl,s_dr,s=-*(red+a);
	int b;
	int i=a/dim;
	int j=a%dim;

	s_d=*(red+((i+1+dim)%dim)*dim+j);
	s_u=*(red+((i-1+dim)%dim)*dim+j);
	s_r=*(red+i*dim+(j+1+dim)%dim);
	s_l=*(red+i*dim+(j-1+dim)%dim);
	s_ul=*(red+((i-1+dim)%dim)*dim+(j-1+dim)%dim);
	s_ur=*(red+((i-1+dim)%dim)*dim+(j+1+dim)%dim);
	s_dl=*(red+((i+1+dim)%dim)*dim+(j-1+dim)%dim);
	s_dr=*(red+((i+1+dim)%dim)*dim+(j+1+dim)%dim);
	
	if(B)
		b=s*(s_u+s_d+s_l+s_r-s_ul-s_ur-s_dl-s_dr)+8+(B/fabs(B)*s+1)/2;
	else
		b=s*(s_u+s_d+s_l+s_r-s_ul-s_ur-s_dl-s_dr)+8;
	
	return b;
}

int definir_valores_de_p_2(float *p, float B, float J){
	
	*p=powl(E,-(-16*J+2*B));
	*(p+1)=powl(E,-(-16*J-2*B));
	*(p+2)=powl(E,-(-12*J+2*B));
	*(p+3)=powl(E,-(-12*J-2*B));
	*(p+4)=powl(E,-(-8*J+2*B));
	*(p+5)=powl(E,-(-8*J-2*B));
	*(p+6)=powl(E,-(-4*J+2*B));
	*(p+7)=powl(E,-(-4*J-2*B));
	*(p+8)=powl(E,-(2*B));
	*(p+9)=powl(E,-(-2*B));
	*(p+10)=powl(E,-(4*J+2*B));
	*(p+11)=powl(E,-(4*J-2*B));
	*(p+12)=powl(E,-(8*J+2*B));
	*(p+13)=powl(E,-(8*J-2*B));
	*(p+14)=powl(E,-(12*J+2*B));
	*(p+15)=powl(E,-(12*J-2*B));
	*(p+16)=powl(E,-(16*J+2*B));
	*(p+17)=powl(E,-(16*J-2*B));
	
	return 0;
}

float medir_e(int *red, float *p, int dim, float B_T, float J_T, FILE *file, float m, int l){

	int s;
	int a;
	
	definir_valores_de_p_2(p,fabs(B_T),J_T);
	
	fprintf(file,"%.3f\t",B_T);
	fprintf(file,"%.3f\t",J_T);
	
	for(int i=0; i<7000; i++){
		
		s=rand()%(dim*dim);	
		a=flippear_spin_frustacion(red,dim,s,B_T);
		
		if(rand()*1.0/M<*(p+a)){
			*(red+s)=-*(red+s);
			m+=*(red+s)*2;
		}
	
	}
	
	for(int j=0; j<l; j++){
		
		fprintf(file,"%.4f\t", m*1.0/(dim*dim));
		
		for(int i=0; i<dim*dim; i++){
			s=rand()%(dim*dim);	
			a=flippear_spin_frustacion(red,dim,s,B_T);
			
			if(rand()*1.0/M<*(p+a)){
				*(red+s)=-*(red+s);
				m+=*(red+s)*2;
			}
						
		}
	}
	fprintf(file,"\n");	

	return m;
}
