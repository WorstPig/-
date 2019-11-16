#ifndef COMMON_H_INCLUDED
#define COMMON_H_INCLUDED

#include<bits/stdc++.h>

#define MAXM 1000				// M个订单
#define MAXN 100				// N * N 个货架
#define MAXH 100				// 周转货架的数量
#define MAXC 100				// 周转货架的槽数
#define MAXQ 200				// Q种商品
#define MAXW 100				// 每个订单最多W种商品 
#define MAXNN 200				// 初始化NN个种群
#define debug true
#define INF 100000			
#define NPOP 0.2				// 初始化种群的%NPOP lions为nomad
#define S 0.8					// %S females in each pride
#define MAXT 1000				// 迭代次数
#define R 0.7					// %R females in each pride roaming
#define PI 3.1415926			// π
#define MAXNN 5000				// 最大的种群数量
#define MAX_PP 100				// 最大的pride数量
#define Ma 0.2					// select %Ma=0.2 females mating with a male
#define Mu 0.2					// %Mu 幼崽的基因进行变异
#define Mi 0.3					// %Mi of the maximum number of females in pride to migrate
using namespace std;

const double EPS = 1E-6;

/*
N * N 个货架, N <= MAXN;
M个订单, M <= MAXM;
H个周转货架, H <= MAXH;
每个周转货架C个槽, C <= MAXC;
Q种商品, Q <= MAXQ;
每个订单最多W种商品, W <= MAXW;
初始化NN个种群, NN <= MAXNN;
迭代T次, T <= MAXT
*/
extern int N, M, H, C, Q, W, NN, T;


// [x,y]
int random(int x, int y);

int random(int x, int y, int except);

// (x,y)
int randomo(int x, int y);

// 返回数组p中的最小值以及索引值
pair<int, int> getMinPair(int *p, int len);

// 对狮子的坐标进行纠正
void convert(int x[MAXM + 1]);

// 返回size个[x,y]的随机整数
vector<int> randList(int x, int y, int size);

// 猎物逃跑
void escape(int prey[MAXM + 1], int lion[MAXM + 1], int preFitness, int newFitness);

// 更新最佳位置
void updateBestPosition(int best[MAXM + 1], int lion[MAXM + 1]);

// 计算狮子坐标的适应度
int fitness(int lion[MAXM + 1]);

// 计算两点间的距离
double dist(int source[MAXM + 1], int dest[MAXM + 1]);

// 获取两点间的向量
vector<int> getVector(int source[MAXM + 1], int dest[MAXM + 1]);

// 生成0~1的随机数
double random();

// 随机生成一只狮子
vector<int> generateOneLion(int p);

bool isMale(int sexi);

bool isFemale(int sexi);

// 获取该订单的批次数
int getBatchNumber(int x[MAXM + 1]);

#endif // COMMON_H_INCLUDED
