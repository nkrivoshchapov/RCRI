#include <iostream>
#include <stdio.h>
#include <math.h>
using namespace std;
extern "C" void __main_pep_MOD_runtlc(double[6],double[7],double[2],double[][3],double[][3],double[][3],double[][3][3],double[][3][3],double[][3][3], int*);

double pi = 3.141592653589793238462643383279502884197;

extern "C" struct xyz_arr /* Nx3 */
{
    int nconf; /* number of elements */
    double arr[16][9][3]; /* ptr to vector elements*/
};

extern "C" struct str_data
{
    double leng[9];
	double vang[9];
	double tang[3];
};

extern "C" void genconf(str_data *inp, xyz_arr *xyz){
	double a1 = inp->vang[7]/180*pi;
	double a2 = inp->vang[8]/180*pi;
	double t = inp->tang[2]/180*pi;
	double b6 = inp->leng[6];
	double b7 = inp->leng[7];
	double b8 = inp->leng[8];

	double at4[] = {0,0,0};
	double at3[] = {b6,0,0};
	double at2[] = {b7*cos(a1)-b8*(cos(a1)*cos(a2)-cos(t)*sin(a1)*sin(a2)),b7*sin(a1)+b8*(-cos(t)*cos(a1)*sin(a2)-cos(a2)*sin(a1)),-sin(t)*sin(a2)*b8};
	double at1[] = {b7*cos(a1),b7*sin(a1),0};
	double r_n[3][3]={{at1[0],at1[1],at1[2]},{0,0,0},{0,0,0}};
	double r_a[3][3]={{at2[0],at2[1],at2[2]},{0,0,0},{at3[0],at3[1],at3[2]}};
	double r_c[3][3]={{0,0,0},{0,0,0},{at4[0],at4[1],at4[2]}};
	double sr_n[16][3][3];
	double sr_a[16][3][3];
	double sr_c[16][3][3];
	int nsol;

	double b_len[6];
	double b_ang[7];
	double t_ang[2];
	//cout << "\nLengths:";
	for(int i=0;i<6;++i){
		b_len[i] = inp->leng[i];
		//cout << b_len[i] << " ";
	}
	//cout << "\nVA:";
	for(int i=0;i<7;++i){
		b_ang[i] = inp->vang[i]/180*pi;
		//cout << inp->vang[i]<< " ";
	}
	//cout << "\nTA:";
	for(int i=0;i<2;++i){
		t_ang[i] = inp->tang[i]/180*pi;
		//cout << t_ang[i]<< " ";
	}
	__main_pep_MOD_runtlc(b_len,b_ang,t_ang,r_n,r_a,r_c,sr_n,sr_a,sr_c,&nsol);
	
	const int ind[][2] = {{1,0},{2,0},{0,1},{1,1},{2,1},{0,2},{1,2},{2,2},{0,0}};
	xyz->nconf = nsol;
	for(int c=0; c<nsol;++c){
		for(int i=0;i<9;++i){
			if(ind[i][0] == 0){
				xyz->arr[c][i][0] = sr_n[c][ ind[i][1] ][0];
				xyz->arr[c][i][1] = sr_n[c][ ind[i][1] ][1];
				xyz->arr[c][i][2] = sr_n[c][ ind[i][1] ][2];
			}
			
			if(ind[i][0] == 1){
				xyz->arr[c][i][0] = sr_a[c][ ind[i][1] ][0];
				xyz->arr[c][i][1] = sr_a[c][ ind[i][1] ][1];
				xyz->arr[c][i][2] = sr_a[c][ ind[i][1] ][2];
			}
			
			if(ind[i][0] == 2){
				xyz->arr[c][i][0] = sr_c[c][ ind[i][1] ][0];
				xyz->arr[c][i][1] = sr_c[c][ ind[i][1] ][1];
				xyz->arr[c][i][2] = sr_c[c][ ind[i][1] ][2];
			}
		}
		
	}
}
