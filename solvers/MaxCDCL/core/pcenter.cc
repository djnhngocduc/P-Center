/************************************************************************************[SimpSolver.h]
MiniSat -- Copyright (c) 2006,      Niklas Een, Niklas Sorensson
           Copyright (c) 2007-2010, Niklas Sorensson

 Copyright (c) 2017, Mao Luo, Chu-Min LI, Fan Xiao: implementing a learnt clause minimisation approach
 Reference: M. Luo, C.-M. Li, F. Xiao, F. Manya, and Z. L. , “An effective learnt clause minimization approach for cdcl sat solvers,” in IJCAI-2017, 2017, pp.703-711.
 
Maple_CM, Based on Maple_LCM --Copyright (c) 2018, Chu-Min LI, Mao Luo, Fan Xiao: implementing a clause minimisation approach.

Copyright (c) 2021,Chu-Min Li (chu-min.li@u-picardie.fr)

A MaxSAT solver combining branch-and-bound and clause learning, implemented by Chu-Min Li with the help of Jordi Coll
Based on Maple_LCM

Reference: Chu-Min Li, Zhenxing Xu, Jordi Coll, Felip Manyà, Djamal Habet, Kun He, Combining Clause Learning and Branch and Bound for MaxSAT, in CP-2021, to appear
 
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
**************************************************************************************************/

#include <iostream>
#include <fstream>
#include <sstream>
#include "../mtl/Sort.h"
#include "../simp/SimpSolver.h"
#include "../utils/System.h"


using namespace Minisat;
using namespace std;

#define NONE -1

int SolverPcenter::edge_redundant(int fnode, int snode, double dist, int nb_edges) {
	for (unsigned i = 0; i < nb_edges; i++) {
		if (fnode == edges[i]->first && snode == edges[i]->second) {
			if (dist < edges[i]->distance)
				edges[i]->distance = dist;
			return i;
		}
	}
	return NONE;
}
bool SolverPcenter::build_simple_graph_instance(const char* input_file) {
	FILE* fp_in = fopen(input_file, "r");
	if (fp_in == NULL) return false;

	const char* ssc = strstr(input_file, "/");
	int l = 0;
	do {
		l = strlen(ssc) + 1;
		input_file = &input_file[strlen(input_file) - l + 2];
		ssc = strstr(input_file, "/");
	} while (ssc);
	printf("c | %s\n", input_file);
	pcenterCheckFilename.append("solutionCheckOutput_");
	pcenterCheckFilename.append(input_file);

	// //依次读入节点数，边的数目，center数，开始计算时的半径
	// fscanf(fp_in, "%d%d%d%lf", &nbNode, &nbEdge, &nbCenter, &initRadius);
	// currentRadius = 0;
	//initRadius = 0; //没给定初始半径时用这个
	fscanf(fp_in, "%d%d%d", &nbNode, &nbEdge, &nbCenter);

	Edge* spaceForEdges = (Edge*)malloc(nbEdge * (sizeof(Edge)));//存放全部的边，两个端点和距离
	edges = (Edge**)malloc(nbEdge * sizeof(Edge*));//指向每一条边的指针
	visited = (int*)malloc(nbNode * sizeof(int));
	totalDist = 0;//所有边长度的和

	isInitCenter = (int*)malloc(nbNode * sizeof(int));
  isInitNonCenter = (int*)malloc(nbNode * sizeof(int));
  inInitPair = (int*)malloc(nbNode * sizeof(int));

	for (int i = 0; i < nbEdge; i++)
		edges[i] = &spaceForEdges[i];

	//一维数组存放全部的距离
	matriceElements = (double*)malloc(nbNode * nbNode * sizeof(double));
	//二维数组存放指向这些距离数据的指针
	distMatrice = (double**)malloc(nbNode * sizeof(double*));
	for (int i = 0; i < nbNode; i++)
		distMatrice[i] = &matriceElements[i * nbNode];

	int fnode = 0, snode = 0;
	double dist = 0;
	for (int i = 0; i < nbEdge; i++) {
		//读入边，放到edge中
		fscanf(fp_in, "%d%d%lf", &fnode, &snode, &dist);
		//in the input file, nodes are numbered from 1, here they have to be numbered from 0
		fnode--; snode--; totalDist += dist;
		if (fnode == snode)//如果两个端点相同，这条边的数据无意义，删除
			printf("auto edge %d over %d\n", i--, nbEdge--);
		else {
			//把小一点的顶点放在前面
			if (fnode > snode) {
				int node = fnode;
				fnode = snode;
				snode = node;
			}
			int redn;
			//检查这条边是否已经存在，如果存在且当前长度大于已存在的那条，把当前边删除
			if ((redn = edge_redundant(fnode, snode, dist, i)) != NONE)
			{
				printf("c | edge redundant %d over %d with %d\n", i--, nbEdge--, redn);
			}
			else {//如果是一条新边，把两个端点对应的邻居个数增加
			  //Edge* newedge = new(edges[i]) Edge(fnode, snode, dist);
			  edges[i]->first = fnode;  edges[i]->second = snode;  edges[i]->distance = dist;
			}
		}
	}
	//printf("c | total distance %lf \n", totalDist);
	fclose(fp_in);

	for (int i = 0; i < nbNode; i++)
		for (int j = 0; j < nbNode; j++)
			distMatrice[i][j] = totalDist + 1;//初始化边距的矩阵

	//可以把检查边的重复放在这儿，应该会更省事？不对，那就还得更新totalDist
	for (int i = 0; i < nbEdge; i++) {//对每一条存在的边，更新这个边距矩阵
		distMatrice[edges[i]->first][edges[i]->second] = edges[i]->distance;
		distMatrice[edges[i]->second][edges[i]->first] = edges[i]->distance;
	}
	//每个点到自身的距离是0
	for (int i = 0; i < nbNode; i++)
		distMatrice[i][i] = 0;

	free(edges);
	free(spaceForEdges);

	nodeNeibors.growTo(nbNode);
	printf("c | Graph correctly build\n");
	return true;
}

bool SolverPcenter::build_simple_graph_instance_tsp(const char* input_file)
{
	//打开算例文件
	FILE* fp_in = fopen(input_file, "r");
	if (fp_in == NULL) return false;

	//创建一个输出文件
	const char* ssc = strstr(input_file, "/");
	int l = 0;
	do {
		l = strlen(ssc) + 1;
		input_file = &input_file[strlen(input_file) - l + 2];
		ssc = strstr(input_file, "/");
	} while (ssc);
	printf("%s\n", input_file);
	pcenterCheckFilename.append("solutionCheckOutput_");
	pcenterCheckFilename.append(input_file);

	//依次读入节点数，center数
	for(int i = 0; i < 3; i++){
      fscanf(fp_in, "%*[^\n]\n");
  }
  
  fscanf(fp_in, "DIMENSION : %d\n", &nbNode);
  fscanf(fp_in, "%*[^\n]\n");
  fscanf(fp_in, "%*[^\n]\n");

  isInitCenter = (int*)malloc(nbNode * sizeof(int));
  isInitNonCenter = (int*)malloc(nbNode * sizeof(int));
  inInitPair = (int*)malloc(nbNode * sizeof(int));

	visited = (int*)malloc(nbNode * sizeof(int));
	nodes_x = (double*)malloc(nbNode * sizeof(double));
	nodes_y = (double*)malloc(nbNode * sizeof(double));
	totalDist = 0;//所有边长度的和


  double node_x = 0, node_y = 0;
  int num = 0;
	for (int i = 0; i < nbNode; i++) {
		fscanf(fp_in, "%d %lf %lf", &num, &node_x, &node_y);
		nodes_x[i] = node_x;
		nodes_y[i] = node_y;
	}

	currentRadius = 0;

	//读入每个节点的坐标
	
	fclose(fp_in);

	if (compute_dist_matrix) 
		nodeNeibors.growTo(nbNode);
	else
		nodeNeibors.growTo(2);

	printf("c | Graph correctly build\n");

	return true;
}

// compute the shortest distance from start to every other node
void SolverPcenter::Dijkstra(double** distMatrice, int n, int start) {
	double* distance;
	int count, nextnode;

	// if a node is visited, then its shortest distance from start is known
	for (int i = 0; i < n; i++) {
		visited[i] = 0;
	}
	distance = distMatrice[start];
	distance[start] = 0;//start到自身的距离是0
	visited[start] = 1;
	count = 1;

	while (count < n - 1) {
		int mindistance = totalDist+1;

		//找到当前离start距离最小的点及其距离
		for (int i = 0; i < n; i++)
			if (distance[i] < mindistance && !visited[i]) {
				mindistance = distance[i];
				nextnode = i;
			}

		visited[nextnode] = 1;
		//如果start到i的距离大于start到nextnode再到点i的距离，将距离更新为后者
		for (int i = 0; i < n; i++)
			if (!visited[i])
				if (mindistance + distMatrice[nextnode][i] < distance[i]) {
					//if (i == 66)
					//	printf("find %d \n", i);
					distance[i] = mindistance + distMatrice[nextnode][i];
				}
		count++;
	}
}

void SolverPcenter::getNeibors(int n) {

	if (pmed_instance){
		//先清空，否则会继承上一次半径的邻居
		for (int i = 0; i < n; i++)
			nodeNeibors[i].clear();

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (i != j && distMatrice[i][j] <= currentRadius) {
					nodeNeibors[i].push(Neibor(j, distMatrice[i][j]));
				}
			}
		}

		//Print neibors 
		/*printf("_________________________\n");
		printf("NEIBORS\n");
		for (int currNode = 0; currNode < n; currNode++) {
			vec<Neibor>& mainNeibors = nodeNeibors[currNode];
			printf("Node : %d  ->  ", currNode);
			for (int j = 0; j < mainNeibors.size(); j++)
				printf("(%d, %5.2f) ", mainNeibors[j].neibor, mainNeibors[j].dist);
			printf("\n");
		}
		printf("_________________________\n");*/

	}else{

		if (compute_dist_matrix){
			for (int i = 0; i < n; i++)
				nodeNeibors[i].clear();

			double dist_x = 0, dist_y = 0, dist = 0;
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					dist_x = nodes_x[i] - nodes_x[j];
					dist_y = nodes_y[i] - nodes_y[j];
					dist = sqrt(pow(dist_x, 2) + pow(dist_y, 2));
					//printf("dist before %lf ", dist);
					dist = (int)(100.0 * dist + 0.5) / 100.0;
					if (i != j && dist <= currentRadius) {
						nodeNeibors[i].push(Neibor(j, dist));
					}
				}
			}
		}
	}
	
}


void SolverPcenter::getNeiborsNode(int i, int index){
	nodeNeibors[index].clear();
	double dist_x = 0, dist_y = 0, dist = 0;
	for (int j = 0; j < nbNode; j++) {
		dist_x = nodes_x[i] - nodes_x[j];
		dist_y = nodes_y[i] - nodes_y[j];
		dist = sqrt(pow(dist_x, 2) + pow(dist_y, 2));
		//printf("dist before %lf ", dist);
		dist = (int)(100.0 * dist + 0.5) / 100.0;
		if (i != j && dist <= currentRadius) {
			nodeNeibors[index].push(Neibor(j, dist));
		}
	}
}

int SolverPcenter::getDist(int i, int j){
	double dist_x = 0, dist_y = 0, dist = 0;
	dist_x = nodes_x[i] - nodes_x[j];
	dist_y = nodes_y[i] - nodes_y[j];
	dist = sqrt(pow(dist_x, 2) + pow(dist_y, 2));
	//printf("dist before %lf ", dist);
	dist = (int)(100.0 * dist + 0.5) / 100.0;
	return dist;
}

int SolverPcenter::cost(){

    int cost = 0;
    for (int i = 0; i < nbNode; i++){
        int min = INT_MAX;
        for (int j = 0; j < initCenters.size(); j++){
            if (getDist(i, initCenters[j]) < min){
                min = getDist(i, initCenters[j]);
            }
        }
            
        // printf("min : %d\n", min);
        if (min > cost){
            cost = min;
        }
        
    }
    initCenters.clear();
    return cost;
}

//如果一个点已经被选作非中心点，则它不再活跃
bool SolverPcenter::activeNode(int node) {
  return !isInitNonCenter[node];
}

void SolverPcenter::internalNodesRule1(int currNode, vec<Neibor>& mainNeibors){
	//determine N3, and it it is not vide, remove them
	int nbInternalNodes = 0;
	for (int j = 0; j < mainNeibors.size(); j++) {
		int node = mainNeibors[j].neibor;
		if (visited[node] == 1 && activeNode(node)) {//除了N1和N2，剩下的都是N3
			initNonCenters.push(node);
			nbInternalNodes++;
			isInitNonCenter[node] = 1;
			//printf("initNonCenters by R1 N3: %d\n", node + 1);
		}
	}

	if (nbInternalNodes > 0) {//如果N3不为空，则N2中的点也不作为center
		for (int j = 0; j < mainNeibors.size(); j++) {
			int node = mainNeibors[j].neibor;
			if (visited[node] == 3 && activeNode(node)) {
				initNonCenters.push(node);
				isInitNonCenter[node] = 1;
				//printf("initNonCenters by R1 N2: %d\n", node + 1);
			}
		}
		if (!isInitCenter[currNode]) {//currNode必须得是center
			initCenters.push(currNode);
			isInitCenter[currNode] = 1;
			//printf("initCenters: %d\n", currNode + 1);
		}
	}
}

//识别出center和noncenter的点，将他们分别加入对应数组
void SolverPcenter::rule1(int n) {
	for (int currNode = 0; currNode < n; currNode++) {
		if (!activeNode(currNode))
			continue;
		for (int node = 0; node < n; node++)
			visited[node] = 0;
		//标记currNode和它的所有邻居，相当于标记 N[currNode]
		visited[currNode] = 1;

		int index = currNode;
		if (!compute_dist_matrix){
			getNeiborsNode(currNode, 0);
			index = 0;
		}

		vec<Neibor>& mainNeibors = nodeNeibors[index];

		for (int j = 0; j < mainNeibors.size(); j++)
		visited[mainNeibors[j].neibor] = 1;

		//For each free neibor of currNode, see if it is adjacent to a node not adjacent to currNode
		// determine N1 (whose visited will be 2)
		//从currNode的邻居中找出N1集合
		for (int j = 0; j < mainNeibors.size(); j++) {
			int node = mainNeibors[j].neibor;
			if (isInitCenter[node])//如果已经是center了则跳过
				visited[node] = 2;
			else if (activeNode(node)) {//如果node依旧活跃，则检查node的邻居
					index = node;
					if (!compute_dist_matrix){
						getNeiborsNode(node, 1);
						index = 1;
					}
					vec<Neibor>& neibors = nodeNeibors[index];
				
				for (int k = 0; k < neibors.size(); k++)
					if (visited[neibors[k].neibor] == 0 && activeNode(neibors[k].neibor)) {//node的邻居有N[v]以外的点
						//node has a neibor not adjacent to currNode
						visited[node] = 2;
						break;
					}
			}
		}

		// determine N2
		for (int j = 0; j < mainNeibors.size(); j++) {
			int node = mainNeibors[j].neibor;
			if (visited[node] == 1 && activeNode(node)) {//检查N1集合之外的所有活跃点
				index = node;
				if (!compute_dist_matrix){
					getNeiborsNode(node, 1);
					index = 1;
				}
				vec<Neibor>& neibors = nodeNeibors[index];
				for (int k = 0; k < neibors.size(); k++)
					if (visited[neibors[k].neibor] == 2 && activeNode(neibors[k].neibor)) {//node与N1中的点相邻，则属于N2
						visited[node] = 3;
						break;
					}
			}
		}
		internalNodesRule1(currNode, mainNeibors);
		
	}
	//printf("\n\n");
	for (int node = 0; node < n; node++)
		visited[node] = 0;
}


void SolverPcenter::rule2(int n) {
	vec<int> allNeibors;
	for (int currNode1 = 0; currNode1 < n; currNode1++) {
		if (!activeNode(currNode1))
			continue;
		for (int currNode2 = currNode1 + 1; currNode2 < n; currNode2++) {
			if (!activeNode(currNode2))
				continue;
			
			for (int node = 0; node < n; node++)
				visited[node] = 0;
			//确定了两个currNode，标为5
			visited[currNode1] = 5;  visited[currNode2] = 5;
			//determine all free neibors of currNode1 and currNode2, put them into allNeibors
			allNeibors.clear();
			//把currNode1/2的所有邻居标记为1
			vec<Neibor>& mainNeibors1 = nodeNeibors[currNode1];
			for (int j = 0; j < mainNeibors1.size(); j++)
				if (activeNode(mainNeibors1[j].neibor) && visited[mainNeibors1[j].neibor] == 0) {
					visited[mainNeibors1[j].neibor] = 1;
					allNeibors.push(mainNeibors1[j].neibor);
				}
			vec<Neibor>& mainNeibors2 = nodeNeibors[currNode2];
			for (int j = 0; j < mainNeibors2.size(); j++)
				if (activeNode(mainNeibors2[j].neibor) && visited[mainNeibors2[j].neibor] == 0) {
					visited[mainNeibors2[j].neibor] = 1;
					allNeibors.push(mainNeibors2[j].neibor);
				}

			//Determine N1: For each active neibor of currNode1 or currNode2, see if it is adjacent
			//to an active node not adjacent to currNode1 or currNode2, its visited will be 2
			//把属于N1的所有点的visited标记为2
			for (int j = 0; j < allNeibors.size(); j++) {
				int node = allNeibors[j];
				if (isInitCenter[node] || inInitPair[node])
					visited[node] = 2;
				else {
					vec<Neibor>& neibors = nodeNeibors[node];
					for (int k = 0; k < neibors.size(); k++)
						if (visited[neibors[k].neibor] == 0 && activeNode(neibors[k].neibor)) {
							//node has a free neibor not adjacent to currNode
							visited[node] = 2;
							break;
						}
				}
			}

			//determine N2: if node is not N1 and has at least one N1 neibor (its visited will be 3)
			//把属于N2的所有点的visited标记为3
			for (int j = 0; j < allNeibors.size(); j++) {
				int node = allNeibors[j];
				if (visited[node] == 1) {
					vec<Neibor>& neibors = nodeNeibors[node];
					for (int k = 0; k < neibors.size(); k++)
						if (visited[neibors[k].neibor] == 2) { //&& freeNode(neibors[k].neibor)) {
							visited[node] = 3;
							break;
						}
				}
			}

			//determine N3: if node is not N1, nor N2 (its visited remains 1, and should be active and
			// is not init center)
			int nbInternalNodes = 0;
			for (int j = 0; j < allNeibors.size(); j++)
				if (visited[allNeibors[j]] == 1)//剩下的visited没变的点都属于N3
					nbInternalNodes++;

			//如果N3不为空
			if (nbInternalNodes > 0) {
				//N3 is not empty
				int bySingle = 0;
				for (int j = 0; j < allNeibors.size(); j++) {
					int node = allNeibors[j];
					//如果N2或者N3中的某个点邻接N3中的全部点，则跳过？
					if (visited[node] == 1 || visited[node] == 3) {
						// check if any single node in N2 or N3 dominate all N3
						int nbN3 = (visited[node] == 1) ? 1 : 0;
						vec<Neibor>& neibors = nodeNeibors[node];
						for (int k = 0; k < neibors.size(); k++)
							if (visited[neibors[k].neibor] == 1) //&& freeNode(neibors[k].neibor))
								nbN3++;
						if (nbN3 == nbInternalNodes) {
							bySingle = 1;
							break;
						}
					}
				}
				//这儿的代码也没看懂
				if (!bySingle) {
					int nbN31 = 0, nbN32 = 0;//分别检查currNode1/2的邻居中有没有N3的点
					for (int j = 0; j < mainNeibors1.size(); j++)
						if (visited[mainNeibors1[j].neibor] == 1)
							nbN31++;
					for (int j = 0; j < mainNeibors2.size(); j++)
						if (visited[mainNeibors2[j].neibor] == 1)
							nbN32++;
					for (int j = 0; j < allNeibors.size(); j++)
						if (visited[allNeibors[j]] == 1) {
							initNonCenters.push(allNeibors[j]);
							isInitNonCenter[allNeibors[j]] = 1;
							//printf("initNonCenters by R2 for N3: %d\n", allNeibors[j] + 1);
						}
					if (nbN31 == nbInternalNodes && nbN32 == nbInternalNodes) {
						//case 2.1
						if (isInitCenter[currNode1]) {
							//printf("case 2.1, %d v %d, already initCenter %d\n",
							//	currNode1 + 1, currNode2 + 1, currNode1 + 1);
							for (int j = 0; j < mainNeibors1.size(); j++)
								if (visited[mainNeibors1[j].neibor] == 3) {
									initNonCenters.push(mainNeibors1[j].neibor);
									isInitNonCenter[mainNeibors1[j].neibor] = 1;
									//printf("initNonCenters by R2.1+N2^Nv^Nw(-Nw): %d\n", mainNeibors1[j].neibor + 1);
								}
						}
						if (isInitCenter[currNode2]) {
							//printf("case 2.1, %d v %d, already initCenter %d\n",
							//	currNode1 + 1, currNode2 + 1, currNode2 + 1);
							for (int j = 0; j < mainNeibors2.size(); j++)
								if (visited[mainNeibors2[j].neibor] == 3) {
									initNonCenters.push(mainNeibors2[j].neibor);
									isInitNonCenter[mainNeibors2[j].neibor] = 1;
									//printf("initNonCenters by R2.1+N2^Nv^Nw(-Nv): %d\n", mainNeibors2[j].neibor + 1);
								}
						}
						if (!isInitCenter[currNode1] && !isInitCenter[currNode2]) {
							for (int j = 0; j < mainNeibors1.size(); j++)
								if (visited[mainNeibors1[j].neibor] == 3)
									visited[mainNeibors1[j].neibor] = 4;
							for (int j = 0; j < mainNeibors2.size(); j++)
								if (visited[mainNeibors2[j].neibor] == 4) {
									initNonCenters.push(mainNeibors2[j].neibor);
									isInitNonCenter[mainNeibors2[j].neibor] = 1;
									//printf("initNonCenters by R2.1+N2^Nv^Nw: %d\n", mainNeibors2[j].neibor + 1);
								}
							initCenterPairs.push(currNode1); initCenterPairs.push(currNode2);
							inInitPair[currNode1] = 1; 	    inInitPair[currNode2] = 1;
							//printf("init pair by R2.1: %d v %d\n", currNode1 + 1, currNode2 + 1);
						}
					}
					else if (nbN31 == nbInternalNodes) {
						//case 2.2
						for (int j = 0; j < mainNeibors1.size(); j++)
							if (visited[mainNeibors1[j].neibor] == 3) {
								initNonCenters.push(mainNeibors1[j].neibor);
								isInitNonCenter[mainNeibors1[j].neibor] = 1;
								//printf("initNonCenters by R2.2+N2^Nv: %d\n", mainNeibors1[j].neibor + 1);
							}
						if (!isInitCenter[currNode1]) {
							initCenters.push(currNode1);
							isInitCenter[currNode1] = 1;
							//printf("initCenters by R2.2v: %d(w=%d)\n", currNode1 + 1, currNode2 + 1);
						}
					}
					else if (nbN32 == nbInternalNodes) {
						//case 2.3
						for (int j = 0; j < mainNeibors2.size(); j++)
							if (visited[mainNeibors2[j].neibor] == 3) {
								initNonCenters.push(mainNeibors2[j].neibor);
								isInitNonCenter[mainNeibors2[j].neibor] = 1;
								//printf("initNonCenters by R2.2+N2^Nw: %d\n", mainNeibors2[j].neibor + 1);
							}
						if (!isInitCenter[currNode2]) {
							initCenters.push(currNode2);
							isInitCenter[currNode2] = 1;
							//printf("initCenters by R2.2w: %d(v=%d)\n", currNode2 + 1, currNode1 + 1);
						}
					}
					else {
						//case 2.4
						for (int j = 0; j < allNeibors.size(); j++) {
							int node = allNeibors[j];
							if (visited[node] == 3) {
								initNonCenters.push(node);
								isInitNonCenter[node] = 1;
								//printf("initNonCenters by R2.2+N2: %d\n", node + 1);
							}
						}
						if (!isInitCenter[currNode1]) {
							initCenters.push(currNode1);
							isInitCenter[currNode1] = 1;
						}
						if (!isInitCenter[currNode2]) {
							initCenters.push(currNode2);
							isInitCenter[currNode2] = 1;
						}
						//printf("initCenters by R2.4: %d, %d\n", currNode1 + 1, currNode2 + 1);
					}
				}
			}
		}
	}
	//printf("\n\n");
	for (int node = 0; node < n; node++)
		visited[node] = 0;
}

//根据给定的dist划定邻居，再化简这个图
void SolverPcenter::reduceGraph(int n) {
	//距离小于dist的两点互为邻居
	getNeibors(n);

	//全部都置0
	initCenters.clear();
	initNonCenters.clear();
	initCenterPairs.clear();
	for (int currNode = 0; currNode < nbNode; currNode++) {
		isInitCenter[currNode] = 0;
		isInitNonCenter[currNode] = 0;
		inInitPair[currNode] = 0;
	}

	if (verbosity > 1)
		printf("\nc Graph reduction for radius %5.2lf\n", currentRadius);

	//rule2没看太懂，先注释一下
	// rule2(n);
	// rule1(n);

	/*if (initCenters.size() > 0)
		for (int i = 0; i < initCenters.size(); i++)
			printf("initCenters: %d\n", initCenters[i] + 1);
	else printf("no init center\n");

	if (initNonCenters.size() > 0)
		for (int i = 0; i < initNonCenters.size(); i++)
			printf("initNonCenters: %d\n", initNonCenters[i] + 1);
	else printf("no init non center\n");

	if (initCenterPairs.size() > 0)
		for (int i = 0; i < initCenterPairs.size(); i += 2)
			printf("init pair by R2.1: %d, %d\n", initCenterPairs[i] + 1, initCenterPairs[i + 1] + 1);
	else printf("no init pair\n");*/

	if (verbosity > 1){
		if (initCenters.size() > 0)
			printf("c nb initCenter %d\n", initCenters.size());
		else printf("c no init center\n");

		if (initNonCenters.size() > 0)
			printf("c nb initNonCenter %d\n", initNonCenters.size());
		else printf("c no init non center\n");

		if (initCenterPairs.size() > 0)
			printf("c initCenterPairs %d\n", initCenterPairs.size());
		else printf("c no init pair\n");
	}
}


//在读入算例时变量个数就已经确定了，因此可以初始化与变量相关的部分，以及生成记录了所有半径的数组并排序
void SolverPcenter::initializeAndSortDistance()
{
	//每个点至少被一个集合覆盖，一共nbNode个硬子句
	hardWeight = nbNode + 1;
	nbOrignalVars = nbNode;
	totalWeight = 0;

	//把图中所有距离从小到大排个序
	for (int i = 0; i < nbNode - 1; i++)
		for (int j = i + 1; j < nbNode; j++)
			allDistances.push(distMatrice[i][j]);
	sort(allDistances);

	//如果有初始半径，则可以把初始半径往后的所有距离都删除
	if (initRadius)
	{
		int i;
		for (i=0; i < allDistances.size(); i++)
			if (allDistances[i] == initRadius)
				break;
		i++;
		allDistances.shrink(allDistances.size() - i);
	}


	//删除重复的半径, when allDistances is sorted in increasing order
	int i, j;
	for(i=1, j=0; i<allDistances.size(); i++)
	  if (allDistances[j] != allDistances[i])
	    allDistances[++j] = allDistances[i];
	allDistances.shrink(i-(j+1));
	

	// //只剩一个元素的数组不需要删除重复了
	// if (allDistances.size() <= 1)
	//   return;
	
	// //删除重复的半径
	// int p = 0, q = 1;
	// while (q < allDistances.size())
	// {
	// 	if (allDistances[p] != allDistances[q])
	// 	{
	// 		p++;
	// 		allDistances[p] = allDistances[q];
	// 	}
	// 	q++;
	// }
	// p++;
	// allDistances.shrink(allDistances.size() - p);
}

void SolverPcenter::initializeAndGetExtremeDistance()
{
	//每个点至少被一个集合覆盖，一共nbNode个硬子句
	hardWeight = nbNode + 1;
	nbOrignalVars = nbNode;
	totalWeight = 0;
	lowestDistance = 100000;
	double dist_x = 0, dist_y = 0, dist = 0;
	highestDistance = 0;

	//把图中所有距离从小到大排个序
	for (int i = 0; i < nbNode - 1; i++)
		for (int j = i + 1; j < nbNode; j++){
			dist_x = nodes_x[i] - nodes_x[j];
			dist_y = nodes_y[i] - nodes_y[j];
			dist = sqrt(pow(dist_x, 2) + pow(dist_y, 2));
			//printf("dist before %lf ", dist);
			dist = (int)(100.0 * dist + 0.5) / 100.0;
			if (dist < lowestDistance)
				lowestDistance = dist;
			if (dist > highestDistance)
				highestDistance = dist;
		}
}

void SolverPcenter::getDistanceInBetween(double lowRadius, double highRadius)
{
	double dist_x = 0, dist_y = 0, dist = 0;

	for (int i = 0; i < nbNode - 1; i++)
		for (int j = i + 1; j < nbNode; j++){
			dist_x = nodes_x[i] - nodes_x[j];
			dist_y = nodes_y[i] - nodes_y[j];
			dist = sqrt(pow(dist_x, 2) + pow(dist_y, 2));
			//printf("dist before %lf ", dist);
			//dist = (int)(100.0 * dist + 0.5) / 100.0;
			if (dist <= highRadius && dist >= lowRadius)
				candidates.push(dist);
		}

	sort(candidates);
	for (int i = 1, j = 1; i < candidates.size(); i++){
		if (candidates[i-1] != candidates[i])
			candidates[j++] = candidates[i];
	}

}

int SolverPcenter::isAllAssignsUndef(){
	int nbAssigned = 0;
	for (int i = 0; i < assigns.size(); i++){
		Var x = i;
		if (value(x) != l_Undef){
			nbAssigned++;
		}
	}
	return nbAssigned;
}

void SolverPcenter::resetAssigns(){
	for (int i = 0; i < assigns.size(); i++){
		Var x = i;
		assigns[x] = l_Undef;
	}
}

//将集覆盖问题转为maxSAT算例，并将那些子句加入子句库
void SolverPcenter::createMaxsatInstances()
{

	resetAssigns();

	vec<Lit> lits;

	std::ifstream inputFile(filename); // Open the file

    if (!inputFile.is_open()) {
        std::cerr << "Failed to open the file." << std::endl;
        return;
    }

    std::string line;
    int numClauses = 0;
    int numVariables = 0;
    int maxWeight = 0;

    while (std::getline(inputFile, line)) {
        std::istringstream iss(line);

        // Check if the line starts with 'p wcnf'
        if (line.compare(0, 6, "p wcnf") == 0) {
            // Parse the header line
            iss.ignore(6); // Ignore 'p wcnf'
            iss >> numVariables >> numClauses >> maxWeight;

            initLB = 0;
		    initUB = nbCenter + 1;
		    nbNode = numVariables;
			totalWeight = 0;
			for (int i = 0; i < numVariables; i++){
		        newVar();
		    }
		    nbOrignalVars = nVars();
            continue; // Skip the header line
        }

        hardWeight = nbNode + 1;
        int weight;
        int node;
        lits.clear();
        iss >> weight;

        while (iss >> node) {
        	if (node < 0)
        		lits.push(~mkLit(-node-1));
        	else if (node > 0)
           		lits.push(mkLit(node-1));
        }
        addClause_(lits, weight);
    }

    inputFile.close();


}

//输出半径，每个点的邻居点和MaxCDCL给出的解，用来检查求解是否正确。即在当前半径下，是否每个点都被某个中心点所覆盖
void SolverPcenter::printfSolutionCheckOutput(double radius)
{
	FILE* fpSolutionOut = fopen(pcenterCheckFilename.c_str(), "w");
	//文件的第一行输出了半径
	fprintf(fpSolutionOut, "%lf\n", radius);

	/*
	//文件的第一行输出了半径和变量个数
	fprintf(fpSolutionOut, "%lf %d %d\n", radius, nbOrignalVars, nbCenter);

	//文件的第二部分是每个点的邻居点
	for (int i = 0; i < nbNode; i++)
	{
		for (int j = 0; j < nodeNeibors[i].size(); j++)
		{
			int node = nodeNeibors[i][j].neibor;
			fprintf(fpSolutionOut, "%d ", node + 1);
		}
		fprintf(fpSolutionOut, "\n");
	}*/

	fclose(fpSolutionOut);
}

int SolverPcenter::readGraphInitUBSortDist() {
	if (pmed_instance){

		/*if (!build_simple_graph_instance_tsp_distMatrice(filename))
			return NONE;*/
		if (!build_simple_graph_instance(filename))
			return NONE;

		//把读入的东西再打印出来看看对不对
		/*FILE* fpOut = fopen("output", "w");
		fprintf(fpOut, " %d %d %d \n", nbNode, nbEdge, nbCenter);
		for (int i = 0; i < nbEdge; i++)
			fprintf(fpOut, " %d %d %lf \n", edges[i]->first + 1, edges[i]->second + 1, edges[i]->distance);
		fclose(fpOut);*/

		//不是任意两点之间都有路径，两点间的直接路径也许不是最短，因此要计算每两点之间的最短距离
		//dijkstra = false;//在tsp下两点之间的路径已经是最短了
		//printf("c | dist option %d\n", dijkstra);
		if (dijkstra)
			for (int i = 0; i < nbNode; i++)
				Dijkstra(distMatrice, nbNode, i);

		//检查是否距离计算出错了
		for (int i = 0; i < nbNode; i++)
			for (int j = 0; j < nbNode; j++)
				if (distMatrice[i][j] != distMatrice[j][i])
					printf("error dist %d->%d %5.2lf, %d->%d %5.2lf\n", i, j, distMatrice[i][j], j, i, distMatrice[j][i]);


		initializeAndSortDistance();	

	}else{
		if (!build_simple_graph_instance_tsp(filename))
			return NONE;

		initializeAndGetExtremeDistance();
	}

	initUB = nbCenter + 1;//initUB是算例中给出的p值+1，因为UB是不能被达到的
	initLB = 0;
	return EXIT_SUCCESS; // return 0;
	
}


void SolverPcenter::removeclauses(vec<CRef>& cs) {
  for(int i=0; i<cs.size(); i++)
    removeClause(cs[i]);
  cs.clear();
}

void SolverPcenter::cancelAll(lbool prevStatus) {
  for (int c = trail.size()-1; c >= 0; c--){
    Var      x  = var(trail[c]);
    if (!VSIDS){
      uint32_t age = conflicts - picked[x];
      if (age > 0){
	double adjusted_reward = ((double) (conflicted[x] + almost_conflicted[x])) / ((double) age);
	double old_activity = activity_CHB[x];
	activity_CHB[x] = step_size * adjusted_reward + ((1 - step_size) * old_activity);
	if (order_heap_CHB.inHeap(x)){
	  if (activity_CHB[x] > old_activity)
	    order_heap_CHB.decrease(x);
	  else
	    order_heap_CHB.increase(x);
	}
      }
#ifdef ANTI_EXPLORATION
      canceled[x] = conflicts;
#endif
    }
    assigns [x] = l_Undef;
    if (phase_saving > 1 || (phase_saving == 1) && c > trail_lim.last())
      polarity[x] = sign(trail[c]);
    insertAuxiVarOrder(x);
    insertVarOrder(x);
    seen[x] = 0;
  }
  qhead = 0;
  trail.shrink(trail.size());
  trail_lim.shrink(trail_lim.size());
  falseLits.shrink(falseLits.size());
  falseLits_lim.shrink(falseLits_lim.size());

  // clear all var data except activities
  for(Var v=0; v<nVars(); v++) {
    vardata[v].reason = CRef_Undef;
    Lit p=mkLit(v);
    picked[v]=0; conflicted[v]=0; almost_conflicted[v]=0;
#ifdef ANTI_EXPLORATION
    canceled[v]=0;
#endif
    seen[v]=0; seen2[toInt(p)]=0; seen2[toInt(~p)]=0;
    decision[v]=true; imply[toInt(p)]=lit_Undef; imply[toInt(~p)]=lit_Undef;
    softVarLocked[v]=0; unlockReason[v]=var_Undef; inConflicts[v]=NON; inConflict[v]=NON;
    involved[v]=0; eliminated[v]=0; litTrail[toInt(p)]=0; litTrail[toInt(~p)]=0;
    softLits[v] = lit_Undef;
    touched[v]=0;
    watches[p].clear();  watches[~p].clear();
    watches_bin[p].clear();  watches_bin[~p].clear();
  }

  vardata.shrink(vardata.size() - nbNode);
  
  eliminatedVars.clear();
  elimclauses.clear();

  equivLits.clear();
  feasibleNbEq=0; prevEquivLitsNb=0; nbSoftEq=0; prevNbSoftEq=0;
  dynVars.clear();
  nbCompSoftLitPairs=0;
  
  involvedLits.clear();

  if (keepLearnClauses && prevStatus == l_True){
	  printf("Sizes : %d %d %d %d\n", clauses.size(), learnts_local.size(), learnts_tier2.size(), learnts_core.size());
	  //learnts_core_storage.clear();

	  totalLearntsCore+=learnts_core.size();
	  learnts_core_storage.growTo(totalLearntsCore);
	  for(int i = 0; i < learnts_core.size(); i++){
	  	Clause& c = ca[learnts_core[i]];
	  	for(int j = 0; j < c.size(); j++){
	  		learnts_core_storage[totalLearntsCore - learnts_core.size() + i].push(c[j]);
	  	}
	  }

	  totalLearntsTier2+=learnts_tier2.size();
	  learnts_tier2_storage.growTo(totalLearntsTier2);
	  for(int i = 0; i < learnts_tier2.size(); i++){
	  	Clause& c = ca[learnts_tier2[i]];
	  	for(int j = 0; j < c.size(); j++){
	  		learnts_tier2_storage[totalLearntsTier2 - learnts_tier2.size() + i].push(c[j]);
	  	}
	  }

	  /*totalLearntsLocal+=learnts_local.size();
	  learnts_local_storage.growTo(totalLearntsLocal);
	  for(int i = 0; i < learnts_local.size(); i++){
	  	Clause& c = ca[learnts_local[i]];
	  	for(int j = 0; j < c.size(); j++){
	  		learnts_local_storage[totalLearntsLocal - learnts_local.size() + i].push(c[j]);
	  	}
	  }*/
	}

  removeclauses(clauses);
  removeclauses(learnts_local);
  removeclauses(learnts_tier2);
  removeclauses(learnts_core);
 

  removeclauses(hardens);
  removeclauses(isetClauses);
  removeclauses(cardinalityC);
  checkGarbage();

  softClauses.clear();
  solutionCost = 0;  totalWeight=0; clauses_literals=0; elim_heap.clear();
}

static double luby(double y, int x){
    
    // Find the finite subsequence that contains index 'x', and the
    // size of that subsequence:
    int size, seq;
    for (size = 1, seq = 0; size < x+1; seq++, size = 2*size+1);
    
    while (size-1 != x){
        size = (size-1)>>1;
        seq--;
        x = x % size;
    }
    
    return pow(y, seq);
}

lbool SolverPcenter::solvePcenter() {
    model.clear(); usedClauses.clear();
    conflict.clear();
    if (!ok) return l_False;
   	lbool status;
  	lbool prevStatus = status;
    status = l_Undef; fixedCostBySearch = 0;
    createMaxsatInstances();
    if (solutionCost >= initUB)
      status = l_False;
    
    max_learnts               = nClauses() * learntsize_factor;
    learntsize_adjust_confl   = learntsize_adjust_start_confl;
    learntsize_adjust_cnt     = (int)learntsize_adjust_confl;
    
    if (status == l_Undef)
      UB=initUB-solutionCost;
    else UB=initUB;

    if (verbosity > 1) printf("c start initialization\n");

    if (status == l_Undef && !initialization())
      status = l_False;

  	if (verbosity > 1)
      printf("c Initialization done !\n");
  	// exit(0);


  	// if (keepLearnClauses){

	//   	for(int i = 0; i < learnts_core_storage.size(); i++){
		  	
	// 	  	CRef cr = ca.alloc(learnts_core_storage[i], true);
	// 	  	claBumpActivity(ca[cr]);
	//         attachClause(cr);
	//         //ajouter dans learnts_core
	//         learnts_core.push(cr);
	// 	}

	// 	for(int i = 0; i < learnts_tier2_storage.size(); i++){
		  	
	// 	  	CRef cr = ca.alloc(learnts_tier2_storage[i], true);
	// 	  	claBumpActivity(ca[cr]);
	//         attachClause(cr);
	//         learnts_tier2.push(cr);
	// 	}

	// 	/*for(int i = 0; i < learnts_local_storage.size(); i++){
		  	
	// 	  	CRef cr = ca.alloc(learnts_local_storage[i], true);
	// 	  	claBumpActivity(ca[cr]);
	//         attachClause(cr);
	// 	}*/
  	// }


    if (initUB<=fixedCostBySearch+solutionCost+derivedCost)
    	status = l_False;
    
    rebuildOrderHeap();

    if (status == l_Undef)
      UB=initUB-fixedCostBySearch-solutionCost-derivedCost;
    else UB=initUB;
    uint64_t prevUB=UB;
    int nbVSIDSphase=0, nbLRBphase=0;
      add_tmp.clear(); softConflictFlag=false;
      UBconflictFlag=false; softConflictFlag=false; falseVar = var_Undef;
      fflush(stdout);

      if (UB == 1)
	harden();

      VSIDS = true;
      int init = 10000;
      while (status == l_Undef && init > 0 && prevUB==UB)
        status = search(init);
      if (verbosity > 1){
      	printf("c ends of initiationization by VSIDS at %llu conflicts with init %d\n\n", 
	     		conflicts, init);
	    } 
      //  if (!switch_mode)
	VSIDS = false;
      
      // Search:
      uint64_t phase_allotment=20000000;
      int curr_restarts = 0;
      for(; status == l_Undef && prevUB==UB ;) {

	//	uint64_t budget = phase_allotment;
	uint64_t savedUP = propagations;
	uint64_t savedConfl = conflicts;
	uint64_t savedRestarts = starts;
	
        fflush(stdout);

	bool toDo=true;

        while (status == l_Undef && propagations - savedUP < phase_allotment && prevUB==UB) {
	  if (VSIDS) {
	    int weighted = INT32_MAX;
	    status = search(weighted);
	  }
	  else{
	    int nof_conflicts = luby(restart_inc, curr_restarts) * restart_first;
	    curr_restarts++;
	    status = search(nof_conflicts);
	  }
	  if (toDo && status==l_Undef && equivLits.size()>prevEquivLitsNb &&
	      cardinalityC.size()==0 && propagations - savedUP >= phase_allotment/2) {
	    toDo=false;
	    if (!moreEliminateEqLits())
	      status=l_False;
	  }
	}
	if (VSIDS) {
	  nbVSIDSphase++;
	  if (verbosity > 1){
	  printf("c VSIDS phase %d: conflicts %llu, phase %llu, starts %llu, UP %llu\n",
		 nbVSIDSphase, conflicts-savedConfl, phase_allotment,
		 starts-savedRestarts, propagations-savedUP);
	 	}
	}
	else  {
	  nbLRBphase++;
	  if (verbosity > 1){
	  printf("c LRB phase %d: conflicts %llu, phase %llu, starts %llu, UP %llu\n",
		 nbLRBphase, conflicts-savedConfl, phase_allotment,
		 starts-savedRestarts, propagations-savedUP);
		}
	}
	
	VSIDS = !VSIDS;
        if (!VSIDS)
            phase_allotment *= 2;
	if (status==l_Undef && equivLits.size()>prevEquivLitsNb &&
	    cardinalityC.size()==0 && !moreEliminateEqLits())
	  status=l_False;
      }
      assert(status != l_Undef || prevUB>UB);
      float meanLB=0, dev=0, succRate=0;
      if (nbLKsuccess>savednbLKsuccess) {
	meanLB= (float)totalPrunedLB/(nbLKsuccess-savednbLKsuccess);
	dev = sqrt((float)totalPrunedLB2/(nbLKsuccess-savednbLKsuccess) - meanLB*meanLB);
      }
      if (LOOKAHEAD > savedLOOKAHEAD) 
	succRate = (float) (nbLKsuccess-savednbLKsuccess)/(LOOKAHEAD-savedLOOKAHEAD);


      if (status == l_False) {
      	
		cancelAll(status);
		//	printf("c derivedCost %llu, softLits %d, sup %llu\n", derivedCost, allSoftLits.size(), sup);
		next_C_reduce = 0;
		next_L_reduce = 0; next_T2_reduce=0; subconflicts = 0; curSimplify = 1; nbconfbeforesimplify=1000;
		totalPrunedLB=0; totalPrunedLB2=0; savedLOOKAHEAD = LOOKAHEAD; savednbLKsuccess=nbLKsuccess;
				
				      }
      else if (status == l_True || prevUB > UB) {
      			status = l_True;
				int openedCenters = UB+fixedCostBySearch+solutionCost+derivedCost;
				cancelAll(status);
				next_C_reduce = 0;
				next_L_reduce = 0; next_T2_reduce=0; subconflicts = 0; curSimplify = 1; nbconfbeforesimplify=1000;
				
				totalPrunedLB=0; totalPrunedLB2=0; WithNewUB=true; //curSimplify = 1; nbconfbeforesimplify=1000;
				savedLOOKAHEAD = LOOKAHEAD; savednbLKsuccess=nbLKsuccess;
      }
      else {
      	printf("Error\n");
      }

      double cpu_time = cpuTime();
      printf("c CPU time  : %g s\n", cpu_time);
    
    cancelAll(status);

    return status;
}



