#include<bits/stdc++.h>
#include "common.h"
using namespace std;

int random(int x, int y) {
	if (x > y) swap(x, y);
	if (x == y) return x;
	return rand() % (y - x + 1) + x;
}

int random(int x, int y, int except) {
	if (x > y) swap(x, y);
	if (x == y) return x;
	int rest;
	do {
		rest = rand() % (y - x + 1) + x;
	} while (rest == except);
	return rest;
}

int randomo(int x, int y) {
	if (x > y) swap(x, y);
	if (x == y) return x;
	if (x + 1 == y) return y;
	return rand() % (y - x - 1) + x + 1;
}

pair<int, int> getMinPair(int *p, int len) {
	int min_value = INT_MAX;
	int index = -1;
	for (int i = 0; i < len; ++i) {
		if (*(p + i) < min_value) {
			min_value = *(p + i);
			index = i;
		}
	}
	return make_pair(min_value, index);
}


/*
转换第x只lion的编码为可行编码，X(x1,x2,x3,...,xn)
1.从小到大，去重排序
2.对于xi等于yj，则另xi=j
3.遍历所有批次，判断批次的订单数是否合法
对于不合法的批次，取出该批次的第1个订单放到订单数最少的批次，如果不合法，则放到一个新建的批次中
*/
void convert(int x[MAXM + 1]) {
	int tmp[MAXM + 1];
	for (int i = 1; i <= M; ++i) {
		tmp[i] = x[i];
	}
	sort(tmp + 1, tmp + 1 + M);
	int n = unique(tmp + 1, tmp + 1 + M) - tmp - 1;
	map<int, int> indexMap;
	for (int i = 1; i <= n; ++i) {
		indexMap.insert(make_pair(tmp[i], i));
	}
	for (int i = 1; i <= M; ++i) {
		x[i] = indexMap.find(x[i])->second;
	}
	for (int i = 1; i <= n; ++i) {
		int p[MAXM + 1] = { 0 };
		for (int j = 1; j <= M; ++j) {
			p[x[j]]++;
		}
		//printf("p[%d]=%d\n", i, p[i]);
		while (p[i] > H * C) {
			pair<int, int> minPair = getMinPair(p + 1, n);
			//cout<<"first "<<minPair.first<<" second: "<<minPair.second + 1<<endl;
			if (minPair.first < H * C) {
				for (int j = 1; j <= M; ++j) {
					if (x[j] == i) {
						p[x[j]]--;
						x[j] = minPair.second + 1;
						p[x[j]]++;
						break;
					}
				}
			}
			else {
				for (int j = 1; j <= M; ++j) {
					if (x[j] == i) {
						p[x[j]]--;
						x[j] = ++n;
						p[x[j]]++;
						break;
					}
				}
			}
		}
	}

	//    int p[M + 1] = {0};
	//    int max_value = INT_MIN;
	//    for (int i = 1; i <= M; ++i) {
	//        p[x[i]]++;
	//        max_value = max(max_value, x[i]);
	//    }
	//    cout<<"max_p:"<<max_value<<endl;
	//    for (int i = 1; i <= max_value; ++i) {
	//        cout<<p[i]<<",";
	//    }
	//    cout<<endl;
}

vector<int> randList(int x, int y, int size) {
	vector<int> res;
	res.clear();
	int* flag = new int[y + 1];
	for (int i = 0; i < size; ++i) {
		int r = random(x, y);
		while (flag[r] == 1) {
			r = random(x, y);
		}
		flag[r] = 1;
		res.push_back(r);
	}
	sort(res.begin(), res.end());
	return res;
}


void escape(int prey[MAXM + 1], int lion[MAXM + 1], int preFitness, int newFitness) {
	double pi = (preFitness - newFitness) * 1.0 / preFitness;
	double r = random();
	for (int m = 1; m <= M; ++m) {
		prey[m] = prey[m] * 1.0 + r * pi * (prey[m] - lion[m]);
	}
	convert(prey);
}

void updateBestPosition(int best[MAXM + 1], int lion[MAXM + 1]) {
	for (int m = 1; m <= M; ++m) {
		best[m] = lion[m];
	}
}

double dist(int source[MAXM + 1], int dest[MAXM + 1]) {
	long sum = 0;
	for (int m = 1; m <= M; ++m) {
		sum = sum + (dest[m] - source[m]) * (dest[m] - source[m]);
	}
	return sqrt(sum);
}

vector<int> getVector(int source[MAXM + 1], int dest[MAXM + 1]) {
	vector<int> res;
	for (int m = 1; m <= M; ++m) {
		res.push_back(dest[m] - source[m]);
	}
	return res;
}

double random() {
	return (rand() % 100 / (double)101);
}

vector<int> generateOneLion(int p) {
	int res[MAXM + 1] = {0};
	for (int j = 1; j <= M; ++j) {
		res[j] = random(1, p);
	}
	convert(res);
	return vector<int>(res, res + M + 1);
}

bool isMale(int sexi) {
	return sexi != 0;
}


bool isFemale(int sexi) {
	return sexi == 0;
}

int getBatchNumber(int x[MAXM + 1]) {
	int max_value = INT_MIN;
	for (int m = 1; m <= M; ++m) {
		max_value = max(max_value, x[m]);
	}
	return max_value;
}
