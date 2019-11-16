#pragma warning(disable:4996)
#include<bits/stdc++.h>
#include "common.h"
using namespace std;

int P;							// 最少的批次数
int PP = 4;						// pride的数量
int MN;							// 迭代过程中总的种群数量
int N=8, M=500, H=3, C=20, Q=100, W=8, NN=200, T=100;


int goods[MAXQ + 1][2];			// 第i种商品在第goods[0]行第goods[1]列货架
int cnt[MAXM + 1];				// 第i个订单有cnt[i]个商品（假设每个订单中每种商品1件）
int order[MAXM + 1][MAXW + 1];	// 第i个订单的第j个商品编号为order[i][j]
int lions[MAXNN][MAXM + 1];		// 种群矩阵
int bests[MAXNN][MAXM + 1];		// 每个种群在t次迭代中的最佳位置
int successes[MAXNN][MAXT];		//
int femaleNumber[MAX_PP];		// 第i个pride的female的数量 
int nomadFemaleNumber;
int nomadMaleNumber;


set<int> sparePos;
set<int> pride[100];
set<int> nomad;
int sex[MAXNN] = { 0 };      // 0 为female，1为male



void readInput() {
	scanf("%d %d %d %d %d %d %d %d", &M, &N, &H, &C, &Q, &W, &NN, &T);
	P = ceil(M * 1.0 / (H * C));
	MN = NN;
	// 读取Q种商品所在的货架坐标
	for (int i = 1; i <= Q; ++i) {
		int x, y;
		scanf("%d %d", &x, &y);
		goods[i][0] = x;
		goods[i][1] = y;
	}

	// 读取M个订单
	for (int i = 1; i <= M; ++i) {
		int w;			// 第i个订单包含w个商品
		scanf("%d", &w);
		cnt[i] = w;
		// 读取每个商品的编号
		for (int j = 1; j <= w; ++j) {
			int id_no;
			scanf("%d,", &id_no);
			order[i][j] = id_no;
		}

	}
}


int getSparePos() {
	if (sparePos.size() == 0) {
		MN = MN + 1;
		cout << "MN:" << MN << endl;
		return MN;
	}
	int index = *(sparePos.begin());
	sparePos.erase(index);
	return index;
}


int fitness(int lion[MAXM + 1]) {
	vector<set<pair<int, int> > > vspii(M + 1);
	int max_p = INT_MIN;
	for (int m = 1; m <= M; ++m) {
		max_p = max(max_p, lion[m]);
		for (int k = 1; k <= cnt[m]; ++k) {
			for (int i = 1; i <= k; ++i) {
				pair<int, int> position = make_pair(goods[order[m][i]][0], goods[order[m][i]][1]);
				vspii[lion[m]].insert(position);
			}
		}
	}
	int sum = 0;
	for (int i = 1; i <= max_p; ++i) {
		sum += vspii[i].size();
	}
	return sum;
}

/*
初始化每种商品所属的货架，每个订单的商品
*/
void init() {
	//N1= 8, M1 = 200, H1 = 3, C1 = 20, Q1 = 50, W = 8,
	printf("商品所在货架信息：\n");
	printf("商品编号\t横坐标\t纵坐标\n");
	for (int i = 1; i <= Q; ++i) {
		goods[i][0] = random(1, N);
		goods[i][1] = random(1, N);
		printf("%d\t%d\t%d\n", i, goods[i][0], goods[i][1]);
	}

	for (int i = 1; i <= M; ++i) {
		cnt[i] = random(1, W);

	//	printf("第%d个订单有%d种商品\n", i, cnt[i]);
	}
	printf("订单序列\t品项数\t商品编号\n");
	for (int i = 1; i <= M; ++i) {

			printf("%d\t%d\t", i,cnt[i]);
		for (int j = 1; j <= cnt[i]; ++j) {
			order[i][j] = random(1, Q);
		printf("%d\t", order[i][j]);
		}
		printf("\n");
	}

}

/*
初始化种群
*/
void initPopulation() {
	memset(successes, 0, sizeof(successes));
	for (int i = 1; i <= NN; ++i) {
		successes[i][0] = 1;
		int p = random(P, M);   // 随机第i个解的订单批次为p
								// lion[i] 为1-p的连续整数
		for (int j = 1; j <= M; ++j) {
			lions[i][j] = random(1, p);
		}
		convert(lions[i]);
		updateBestPosition(bests[i], lions[i]);
	}
}

void dividePopulation() {
	nomad.clear();
	sparePos.clear();
	for (int i = 0; i < PP + 1; ++i) {
		pride[i].clear();
	}
	int m = NN * NPOP;
	int flag[MAXNN + 1] = { 0 };
	for (int i = 1; i <= m; ++i) {
		int num = random(1, NN);
		while (flag[num]) {
			num = random(1, NN);
		}
		flag[num] = 1;
		nomad.insert(num);
	}

	int rest = NN - m;
	for (int i = 0, p = rest / PP; i < PP; ++i) {
		for (int j = 1; j <= p; ++j) {
			int num = random(1, NN);
			while (flag[num]) {
				num = random(1, NN);
			}
			flag[num] = 1;
			pride[i].insert(num);
		}
	}
	if (rest % PP != 0) {
		for (int z = 1; z <= NN; ++z) {
			if (!flag[z]) {
				pride[PP].insert(z);
			}
		}
		PP = PP + 1;
	}

	for (int i = 0; i < PP; ++i) {
		int maleSize = pride[i].size() * (1.0 - S);
		femaleNumber[i] = pride[i].size() - maleSize;
		int j = 0;
		for (set<int>::iterator iter = pride[i].begin(); iter != pride[i].end(); ++iter) {
			if (j < maleSize) {
				sex[*(iter)] = 1;
				j++;
			}
		}
	}

	int maleSize = nomad.size() * S;
	nomadMaleNumber = maleSize;
	nomadFemaleNumber = nomad.size() - maleSize;
	int j = 0;
	for (set<int>::iterator iter = nomad.begin(); iter != nomad.end(); ++iter) {
		if (j < maleSize) {
			sex[*(iter)] = 1;
			j++;
		}
	}

	for (int i = 1; i <= 100; ++i) {
		sparePos.insert(++MN);
	}
}

vector<int> hunting(int i, int t) {
	// 随机选择%60的female
	int femaleSize = 0;
	for (set<int>::iterator it = pride[i].begin(); it != pride[i].end(); ++it) {
		if (sex[*it] == 0) {
			femaleSize++;
		}
	}
	int huntingFemaleSize = femaleSize * 0.6;

	vector<int> hunter;
	vector<int> randFemaleList = randList(1, femaleSize, huntingFemaleSize);

	int j = 0, k = 0;
	for (set<int>::iterator it = pride[i].begin(); it != pride[i].end(); ++it) {
		if (sex[*it] == 0) {
			if (j < randFemaleList.size() && ++k == randFemaleList[j]) {
				hunter.push_back(*it);
				++j;
			}
		}
	}

	// 1. generate dummy prey point
	int prey[MAXM + 1] = { 0 };
	for (int h = 0; h < hunter.size(); ++h) {
		int pos = hunter[h];
		for (int m = 1; m <= M; ++m) {
			prey[m] += lions[pos][m];
		}
	}

	for (int m = 1; m <= M; ++m) {
		prey[m] /= hunter.size();        // 取整
	}

	// 2.divide hunters into three sub groups
	// 随机选择%30的 为center
	int centerSize = hunter.size() * 0.3;
	vector<int> centerNo = randList(0, hunter.size() - 1, centerSize);
	set<int> centerNoSet(centerNo.begin(), centerNo.end());
	//3.
	for (int h = 0; h < hunter.size(); ++h) {
		int pos = hunter[h];
		int preFitness = fitness(lions[pos]);
		int bestFitness = fitness(bests[pos]);
		// center move
		if (centerNoSet.find(h) != centerNoSet.end()) {

			for (int m = 1; m <= M; ++m) {
				if (lions[pos][m] < prey[m]) {
					lions[pos][m] = randomo(lions[pos][m], prey[m]);
				}
				else {
					lions[pos][m] = randomo(prey[m], lions[pos][m]);
				}
				convert(lions[pos]);
			}
			int newFitness = fitness(lions[pos]);
			if (newFitness < preFitness) {
				cout << "-------------------------OK Hunting improve--------------------" << endl;
				successes[pos][t] = 1;
				escape(prey, lions[pos], preFitness, newFitness);
			}
			if (newFitness < bestFitness) {
				cout << "-------------------------OK Hunting--------------------" << endl;
				updateBestPosition(bests[pos], lions[pos]);
			}
		}
		else { // wing move

			for (int m = 1; m <= M; ++m) {
				if (2 * prey[m] - lions[pos][m] < prey[m]) {
					lions[pos][m] = randomo(2 * prey[m] - lions[pos][m], prey[m]);
				}
				else {
					lions[pos][m] = randomo(prey[m], 2 * prey[m] - lions[pos][m]);
				}
				convert(lions[pos]);
			}
			int newFitness = fitness(lions[pos]);
			if (newFitness < preFitness) {
				successes[pos][t] = 1;
				escape(prey, lions[pos], preFitness, newFitness);
			}
			if (newFitness < bestFitness) {
				cout << "-------------------------OK Hunting--------------------" << endl;
				updateBestPosition(bests[pos], lions[pos]);
			}
		}

	}
	return hunter;
}

void moving(vector<int> hunters, int p, int t) {
	set<int> remainFemales;
	vector<int> improveFemales;
	for (set<int>::iterator iter = pride[p].begin(); iter != pride[p].end(); ++iter) {
		int pos = *iter;
		if (sex[pos] == 0) {
			remainFemales.insert(pos);
		}
	}
	for (int i = 0; i < hunters.size(); ++i) {
		remainFemales.erase(hunters[i]);
	}

	int ks = 0;
	for (set<int>::iterator iter = remainFemales.begin(); iter != remainFemales.end(); ++iter) {
		ks += successes[*iter][t - 1];
		if (successes[*iter][t - 1] == 1) {
			improveFemales.push_back(*iter);
		}
	}
	int tSize = max(2, (int)(ceil(ks / 2.0)));   // 似乎没有用到

	// for each remain female, select a place to move
	for (set<int>::iterator iter = remainFemales.begin(); iter != remainFemales.end(); ++iter) {
		if (improveFemales.size() == 0) {
			// TODO
			//cout << "moving----------" << endl;
		}
		else {
			// 随机选择一个improve的female为destination，
			//if (improveFemales.size() == 1 && improveFemales[0] == *iter) {
			//	continue;
			//}
			int r = random(0, improveFemales.size() - 1);
			int destPos = improveFemales[r];
			//while (destPos == *iter) {
			//	r = random(0, improveFemales.size() - 1);
			//	destPos = improveFemales[r];
			//}
			double distance = dist(lions[*iter], bests[destPos]);
			vector<int> viR = getVector(lions[*iter], bests[destPos]);
			
			int preFitness = fitness(lions[*iter]);
			int bestFitness = fitness(bests[*iter]);
			
			for (int m = 1; m <= M; ++m) {
				lions[*iter][m] = lions[*iter][m] + 2 * distance * random() * viR[m - 1] + distance;
			}
			convert(lions[*iter]);
			int newFitness = fitness(lions[*iter]);
			if (newFitness < preFitness) {
				cout << "-------------------------OK Moving Improve--------------------" << endl;
				successes[*iter][t] = 1;
			}
			if (newFitness < bestFitness) {
				cout << "-------------------------OK Moving--------------------" << endl;
				updateBestPosition(bests[*iter], lions[*iter]);
			}
		}
	}
}

void prideRoaming(int p, int t) {
	vector<int> male;

	// 1.get resident males
	for (set<int>::iterator it = pride[p].begin(); it != pride[p].end(); ++it) {
		if (isMale(sex[*it])) {
			male.push_back(*it);
		}
	}

	// 2.select %R places of territory randomly,
	int s = ceil(pride[p].size() * R);
	vector<int> randomList = randList(0, pride[p].size() - 1, s);
	vector<int> place;
	int i = 0, j = 0;
	for (set<int>::iterator it = pride[p].begin(); it != pride[p].end(); ++it) {
		if (j < randomList.size() && i == randomList[j]) {
			place.push_back(*it);
			j++;
		}
		++i;
	}

	for (i = 0; i < male.size(); ++i) {
		int preFitness = fitness(bests[male[i]]);
		int tmpBest[MAXM + 1];
		int bestFitness = preFitness;
		for (j = 0; j < place.size(); ++j) {
			// 3.go toward place
			double distance = dist(lions[male[i]], bests[place[j]]);
			int xunit = randomo(0, ceil(2 * distance));
			//  [-PI/6,PI/6]
			double sita = random() * (PI / 3.0) - (PI / 6.0);
			double offset = tan(sita);
			for (int m = 1; m <= M; ++m) {
				lions[male[i]][m] = lions[male[i]][m]  + xunit + lions[male[i]][m] * offset; // TODO
			}
			convert(lions[male[i]]);
			int newFitness = fitness(lions[male[i]]);
			if (newFitness < bestFitness) {
				bestFitness = newFitness;
				for (int m = 1; m <= M; ++m) {
					tmpBest[m] = lions[male[i]][m];
				}
			}
		}
		if (bestFitness < preFitness) {
			cout << "-------------------------OK Pride Roaming--------------------" << endl;
			updateBestPosition(bests[male[i]], tmpBest);
			updateBestPosition(lions[male[i]], tmpBest);
			successes[male[i]][t] = 1;
		}
	}
}

void nomadRoaming(int t) {
	for (set<int>::iterator it = nomad.begin(); it != nomad.end(); ++it) {
		int pos = *it;
		int p = random(P, P);
		int nomad_i = fitness(lions[pos]);
		int best_nomad = fitness(bests[pos]);
		double pri = 0.1 + min(0.5, (nomad_i - best_nomad) * 1.0 / best_nomad);
		vector<int> lion = generateOneLion(p);
		for (int i = 1; i < lion.size(); ++i) {
			double rand_j = random();
			if (rand_j > pri && fabs(rand_j - pri) > EPS) {
				lions[pos][i] = lions[pos][i];
			}
			else {
				lions[pos][i] = lion[i];
			}
		}
		convert(lions[pos]);
		int newFitness = fitness(lions[pos]);
		if (newFitness < nomad_i) {
			cout << "-------------------------OK Nomad Roaming Improve--------------------" << endl;
			successes[pos][t] = 1;
		}
		if (newFitness < best_nomad) {
			cout << "-------------------------OK Nomad Roaming--------------------" << endl;
			updateBestPosition(bests[pos], lions[pos]);
		}
	}
}


void residentDefense(int p) {
	vector<int> males;
	for (set<int>::iterator it = pride[p].begin(); it != pride[p].end(); ++it) {
		if (isMale(sex[*it])) {
			males.push_back(*it);
		}
	}
	int max_fitness = INT_MIN;
	int drive_out = -1;
	for (int i = 0; i < males.size(); ++i) {
		int fit = fitness(bests[males[i]]);
		if (fit > max_fitness) {
			max_fitness = fit;
			drive_out = males[i];
		}
	}
	pride[p].erase(drive_out);
	nomad.insert(drive_out);
}

void residentMating(int p, int t) {
	// select %Ma females mating with a male
	vector<int> females;
	vector<int> males;
	for (set<int>::iterator it = pride[p].begin(); it != pride[p].end(); ++it) {
		if (sex[*it] == 0) {
			females.push_back(*it);
		}
		else {
			males.push_back(*it);
		}
	}
	int ma = ceil(females.size() * Ma);
	vector<int> list = randList(0, females.size(), ma);
	vector<int> mateFemales;
	for (int i = 0, j = 0; i < females.size(); ++i) {
		if (j < ma && i == list[j]) {
			mateFemales.push_back(females[i]);
			j++;
		}
	}
	int nr = males.size();
	if (nr == 0) return;
	// mating
	for (int i = 0; i < mateFemales.size(); ++i) {
		int posi = mateFemales[i];
		int r = random(0, nr - 1);
		int posj = males[r];
		// 生成2只幼崽
		for (int k = 1; k <= 2; ++k) {
			int newPos = getSparePos();
			for (int m = 1; m <= M; ++m) {
				double b = random() * (0.6 - 0.4) + 0.4;
				if (k == 2) {
					b = 1 - b;
				}
				lions[newPos][m] = b * lions[posi][m] + (1 - b) * lions[posj][m];	// 选择最后一次迭代的解
				// %Mu 变异
				double mu_r = random();
				if (mu_r < Mu && fabs(mu_r - Mu) > EPS) {
					lions[newPos][m] = random(1, P);
				}
			}
			convert(lions[newPos]);
			updateBestPosition(bests[newPos], lions[newPos]);
			successes[newPos][t] = 1;
			pride[p].insert(newPos);			// mature
			sex[newPos] = k%2;					// one female, the other male
			if (isMale(sex[newPos])) {
				residentDefense(p);
			}
		}
	}
}

void nomadMating(int t) {
	vector<int> nomadFemales;
	vector<int> nomadMales;
	for (set<int>::iterator it = nomad.begin(); it != nomad.end(); ++it) {
		if (isMale(sex[*it])) {
			nomadMales.push_back(*it);
		}
		else {
			nomadFemales.push_back(*it);
		}
	}

	int mateFemaleSize = max(int(floor(nomadFemales.size() * 1.0 * Ma)), 1);
	vector<int> list = randList(0, nomadFemales.size(), mateFemaleSize);
	for (int i = 0, j = 0; i < nomadFemales.size(); ++i) {
		if (j < mateFemaleSize && i == list[j]) {
			int posi = nomadFemales[i];
			int r = random(0, nomadMales.size() - 1);
			int posj = nomadMales[r];
			for (int k = 1; k <= 2; ++k) {
				int newPos = getSparePos();
				for (int m = 1; m <= M; ++m) {
					double b = random() * (0.6 - 0.4) + 0.4;
					if (k == 2) {
						b = 1 - b;
					}
					lions[newPos][m] = b * lions[posi][m] + (1 - b) * lions[posj][m];	// 选择最后一次迭代的解
																						// %Mu=0.2 变异
					double mu_r = random();
					if (mu_r < Mu && fabs(mu_r - Mu) > EPS) {
						lions[newPos][m] = random(1, P);
					}
				}
				convert(lions[newPos]);
				updateBestPosition(bests[newPos], lions[newPos]);
				successes[newPos][t] = 1;
				nomad.insert(newPos);				// mature
				sex[newPos] = 0;
				if (k == 2) {						// one female, one male
					sex[newPos] = 1;
				}
			}
		}
	}
}

void nomadDefense(int t) {

	vector<int> nomadMales;
	for (set<int>::iterator it = nomad.begin(); it != nomad.end(); ++it) {
		if (isMale(sex[*it])) {
			nomadMales.push_back(*it);
		}
	}

	int bits[MAX_PP] = { 0 };
	for (int i = 0; i < PP; ++i) {
		bits[i] = random(0, 1);
	}

	for (int i = 0; i < nomadMales.size(); ++i) {
		bool flag = 0;
		for (int j = 0; j < PP && !flag; ++j) {
			if (bits[j] == 1) {
				vector<int> residentMales;
				for (set<int>::iterator it = pride[j].begin(); it != pride[j].end(); ++it) {
					if (isMale(sex[*it])) {
						residentMales.push_back(*it);
					}
				}
				int nr = residentMales.size();
				int nomadFitness = fitness(lions[nomadMales[i]]);
				for (int z = 0; z < nr; ++z) {
					int residentFitness = fitness(lions[residentMales[z]]);
					if (nomadFitness < residentFitness) {
						flag = 1;
						nomad.erase(nomadMales[i]);
						pride[j].erase(residentMales[z]);
						pride[j].insert(nomadMales[i]);
						cout << "--------------nomad defense----------" << endl;
					}
				}
			}
		}
	}

}

// 从大到小排序
bool cmp(pair<int,int> pair1, pair<int,int> pair2) {
	return pair1.second > pair2.second;
}


void migrating(int p, int t) {
	vector<int> females;
	for (set<int>::iterator it = pride[p].begin(); it != pride[p].end(); ++it) {
		if(isFemale(sex[*it])) {
			females.push_back(*it);
		}
	}
	int sum = females.size() - femaleNumber[p];
	sum = sum + floor(femaleNumber[p] * 1.0 * Mi);

	vector<pair<int, int> > list;
	for (int i = 0; i < females.size(); ++i) {
		list.push_back(make_pair(females[i], fitness(bests[females[i]])));
	}
	sort(list.begin(), list.end(), cmp);

	for (int i = 0; i < sum; ++i) {
		pride[p].erase(list[i].first);
		nomad.insert(list[i].first);
	}
}

void fillEmptyPlace() {
	vector<int> females;
	vector<pair<int,int> > femaleList;
	vector<pair<int, int> > maleList;
	for (set<int>::iterator it = nomad.begin(); it != nomad.end(); ++it) {
		if (isFemale(sex[*it])) {
			females.push_back(*it);
			femaleList.push_back(make_pair(*it, fitness(bests[*it])));
		}
		else {
			maleList.push_back(make_pair(*it, fitness(bests[*it])));
		}
	}
	sort(femaleList.begin(), femaleList.end(), cmp);
	int z = femaleList.size() - 1;
	for (int p = 0; p < PP; ++p) {
		int prideFemaleSize = 0;
		for (set<int>::iterator it = pride[p].begin(); it != pride[p].end(); ++it) {
			if (isFemale(sex[*it])) {
				prideFemaleSize++;
			}
		}
		int k = femaleNumber[p] - prideFemaleSize;
		for (int j = 1; j <= k; ++j) {
			pride[p].insert(femaleList[z].first);
			nomad.erase(femaleList[z].first);
			z--;
		}
	}

	//  remove nomad females
	int removeFemaleNum = z + 1 - nomadFemaleNumber;
	if (removeFemaleNum < 0) {
		cout << "error!!! remove Nomad Female error" << endl;
		return;
	}
	for (int i = 0; i < removeFemaleNum; ++i) {
		nomad.erase(femaleList[i].first);
		sparePos.insert(femaleList[i].first);
	}

	// remove nomad males
	sort(maleList.begin(), maleList.end(), cmp);
	int removeMaleNum = maleList.size() - nomadMaleNumber;
	if (removeMaleNum < 0) {
		cout << "error!!! remove Nomad Male error" << endl;
		return;
	}
	for (int i = 0; i < removeMaleNum; ++i) {
		nomad.erase(maleList[i].first);
		sparePos.insert(maleList[i].first);
	}
}

// 判断种群的数量是否保持不变,如果有打印出error，说明算法的逻辑有问题
void validPopulation() {
	// 1. nomad
	int nomadFemales = 0;
	int nomadMales = 0;
	for (set<int>::iterator it = nomad.begin(); it != nomad.end(); ++it) {
		if (isFemale(sex[*it])) {
			nomadFemales++;
		}
		else {
			nomadMales++;
		}
	}
	if (nomadFemales != nomadFemaleNumber) {
		cout << "error!! nomad female size"<<endl;
	}
	if (nomadMales != nomadMaleNumber) {
		cout << "error !! nomad Male size" << endl;
	}

	for (int p = 0; p < PP; ++p) {
		int prideFemales = 0;
		int prideMales = 0;
		for (set<int>::iterator it = pride[p].begin(); it != pride[p].end(); ++it) {
			if (isFemale(sex[*it])) {
				prideFemales++;
			}
			else {
				prideMales++;
			}
		}

		if (prideFemales != femaleNumber[p]) {
			cout << "error!! nomad female size" << endl;
		}
	}
	
}

void printSolution() {
	int minFitness = INT_MAX;
	int resPop = -1;
	for (set<int>::iterator it = nomad.begin(); it != nomad.end(); ++it) {
		int fit = fitness(bests[*it]);
		if (fit < minFitness) {
			minFitness = fit;
			resPop = *it;
		}
	}
	for (int p = 0; p < PP; ++p) {
		for (set<int>::iterator it = pride[p].begin(); it != pride[p].end(); ++it) {
			int fit = fitness(bests[*it]);
			if (fit < minFitness) {
				minFitness = fit;
				resPop = *it;
			}
		}
	}

	printf("%d 个订单, 每批最多%d个订单的情况，将订单分成 %d 批, 最少需要搬运货架 %d 次\n", M, H * C, getBatchNumber(bests[resPop]),minFitness);
	printf("(");
	for (int m = 1; m <= M; ++m) {
		printf("%d", bests[resPop][m]);
		if (m != M) {
			printf(",");
		}
	}
	printf(")\n");

	printf("每个订单的商品种类以及商品数量信息如下：\n");
	for (int i = 1; i <= M; ++i) {
		printf("第%d个订单有%d个商品: ", i, cnt[i]);
		for (int j = 1; j <= cnt[i]; ++j) {
			printf("(%d,[%d,%d])", order[i][j], goods[order[i][j]][0], goods[order[i][j]][1]);
			if (j != cnt[i]) {
				printf(",");
			}
		}
		printf("\n");
	}
}

int main()
{
#ifndef ONLINE_JUDGE
	freopen("in.txt", "r", stdin);
	freopen("out.txt", "w", stdout);
#endif // ONLINE_JUDGE
	srand(time(NULL));
	init();
	//initPopulation();
	/*	readInput();
	dividePopulation();

	for (int t = 1; t < T; ++t) {
		cout << "当前迭代次数："<< t << endl;
		for (int p = 0; p < PP; ++p) {
			vector<int> hunters = hunting(p, t);
			moving(hunters, p, t);
			prideRoaming(p, t);
			residentMating(p, t);
		}
		nomadRoaming(t);
		nomadMating(t);
		nomadDefense(t);

		for (int p = 0; p < PP; ++p) {
			migrating(p, t);
		}

		fillEmptyPlace();
		validPopulation();
	}
	printSolution();*/
	return 0;
}
