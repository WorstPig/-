#ifndef COMMON_H_INCLUDED
#define COMMON_H_INCLUDED

#include<bits/stdc++.h>

#define MAXM 1000				// M������
#define MAXN 100				// N * N ������
#define MAXH 100				// ��ת���ܵ�����
#define MAXC 100				// ��ת���ܵĲ���
#define MAXQ 200				// Q����Ʒ
#define MAXW 100				// ÿ���������W����Ʒ 
#define MAXNN 200				// ��ʼ��NN����Ⱥ
#define debug true
#define INF 100000			
#define NPOP 0.2				// ��ʼ����Ⱥ��%NPOP lionsΪnomad
#define S 0.8					// %S females in each pride
#define MAXT 1000				// ��������
#define R 0.7					// %R females in each pride roaming
#define PI 3.1415926			// ��
#define MAXNN 5000				// ������Ⱥ����
#define MAX_PP 100				// ����pride����
#define Ma 0.2					// select %Ma=0.2 females mating with a male
#define Mu 0.2					// %Mu ���̵Ļ�����б���
#define Mi 0.3					// %Mi of the maximum number of females in pride to migrate
using namespace std;

const double EPS = 1E-6;

/*
N * N ������, N <= MAXN;
M������, M <= MAXM;
H����ת����, H <= MAXH;
ÿ����ת����C����, C <= MAXC;
Q����Ʒ, Q <= MAXQ;
ÿ���������W����Ʒ, W <= MAXW;
��ʼ��NN����Ⱥ, NN <= MAXNN;
����T��, T <= MAXT
*/
extern int N, M, H, C, Q, W, NN, T;


// [x,y]
int random(int x, int y);

int random(int x, int y, int except);

// (x,y)
int randomo(int x, int y);

// ��������p�е���Сֵ�Լ�����ֵ
pair<int, int> getMinPair(int *p, int len);

// ��ʨ�ӵ�������о���
void convert(int x[MAXM + 1]);

// ����size��[x,y]���������
vector<int> randList(int x, int y, int size);

// ��������
void escape(int prey[MAXM + 1], int lion[MAXM + 1], int preFitness, int newFitness);

// �������λ��
void updateBestPosition(int best[MAXM + 1], int lion[MAXM + 1]);

// ����ʨ���������Ӧ��
int fitness(int lion[MAXM + 1]);

// ���������ľ���
double dist(int source[MAXM + 1], int dest[MAXM + 1]);

// ��ȡ����������
vector<int> getVector(int source[MAXM + 1], int dest[MAXM + 1]);

// ����0~1�������
double random();

// �������һֻʨ��
vector<int> generateOneLion(int p);

bool isMale(int sexi);

bool isFemale(int sexi);

// ��ȡ�ö�����������
int getBatchNumber(int x[MAXM + 1]);

#endif // COMMON_H_INCLUDED
