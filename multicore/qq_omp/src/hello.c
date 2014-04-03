#include<stdio.h>
#include<omp.h>
#include <time.h>
static long num_steps = 1000000;

double seq()
{
	
	int i;
	double x, pi, step, sum = 0.0;
	step = 1.0/(double) num_steps;
	for (i=0;i< num_steps; i++){
		x = (i+0.5)*step;
		sum = sum + 4.0/(1.0+x*x);
	}
	pi = step * sum;
	return pi;
}

double par()
{
	double sum[4]={0.0};
	int i;
	double pi = 0;
	double step;
	step = 1.0/(double) num_steps;

	#pragma omp parallel num_threads(4)
	{
		double x;
		int id=omp_get_thread_num();
		for(i=id*(num_steps/4);i< (id+1)*num_steps/4;i++){
				x = (i+0.5)*step;
				sum[id] = sum[id] + 4.0/(1.0+x*x);
		}
		
	}	
	for(i=0; i<4;i++)
		pi += sum[i];
	
	pi = step * pi;

	return pi;


}
void main()
{
	//omp_set_num_threads(5);
	int i;
	clock_t  start, end;
	double time;
	//double out[10000];
	
	start = clock();
		par();
	end = clock();
	time = ((double)(end - start))/ CLOCKS_PER_SEC;
	printf("par time=%f\n",(double)time);

	start = clock();
		seq();
		//out[i] = seq();
	//printf("out=%f\n",out[100]);
	end = clock();
	time = ((double)(end - start))/ CLOCKS_PER_SEC;
	printf("seq time=%f\n",(double)time);
	
	printf("seq pi=%f\n",seq());
	printf("par pi=%f\n",par());

}


