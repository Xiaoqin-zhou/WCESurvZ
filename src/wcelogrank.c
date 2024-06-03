#include <math.h>     // 导入数学库
#include "survS.h"    // 导入survS头文件
#include "survproto.h"// 导入survproto头文件

// 定义函数wcelogrank
void wcelogrank(int *nn, int *nngroup,    //总样本量和总分组数，因为R语言数据结构和C语言不同，所以需要使用指针，将数据导入C语言
               double *time, int *status, int *group, double *weight,    //time，status，weight都是长度为*nn的数组，*group相当于分组标签，也是数组
               double *nrisk1, double *nrisk2,double *nriskall,
               double *obs, double *exp, double *tempt, double *var, double *risk   //观测值，方差，risk都是长度为ngroup的数组
               )
{
    // 注册变量，快速访问
    register int i, j, k; // 用于循环迭代的变量，register关键字表示将它们存储在CPU寄存器中以提高访问速度
    int kk;               // 索引变量，用于二维数组var的访问
    int n;                // 当前分层中的观测数
    int ngroup;           // 组的数量
    int ntot;             // 总观测数
    // int istart;           // 当前分层的起始索引
    // int koff;             // 当前分层组的偏移量，用于obs和exp数组的访问
    double km;            // Kaplan-Meier估计值
    double nrisk;         // 当前风险集中的样本数
    double wt;            // 权重值，根据Kaplan-Meier估计和rho计算得出
    double tmp;           // 临时变量，用于方差计算
    double deaths;        // 当前时间点的死亡数
    // double temp_i;          // 观测值-期望值


    // 赋值变量，使用指针，将R语言的数据格式转换成C语言的数据格式，对于数组指针，不需要赋值
    n = *nn;  //n = 总样本数
    ngroup = *nngroup;      // 组数
    var[0] = 0;   //初始化

    for (i = 0; i < ngroup; i++) {

        obs[i] = 0;    //初始化观测值为0
        exp[i] = 0;   //初始化期望值为0
        // tempt[i] = 0;
        risk[i] = 0;   //初始化风险数组为0
    }

    // 进行实际检验
    for (i = n - 1; i >= 0; i--) {  
        // 计算权重
        // wt = 1;
        k = group[i] - 1;
        obs[k] += status[i] * weight[i];

        exp[0] +=  (weight[i]* nrisk1[i])/(nriskall[i]);
        exp[1] +=  (weight[i]* nrisk2[i])/(nriskall[i]);

        tmp =  (weight[i]* nrisk1[i]*nrisk2[i]*(nriskall[i]-1))/ ((nriskall[i]*nriskall[i]) * (nriskall[i]-weight[i]));
        var[0] += tmp;
    }
}


