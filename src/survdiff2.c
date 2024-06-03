#include <math.h>     // 导入数学库
#include "survS.h"    // 导入survS头文件
#include "survproto.h"// 导入survproto头文件

// 定义函数survdiff2
void survdiff2(int *nn, int *nngroup, int *nstrat, 
               double *rho, double *time, int *status, 
               int *group, int *strata, double *obs, 
               double *exp, double *tempt, double *var, double *risk, 
               double *kaplan)
{
    // 注册变量，快速访问
    register int i, j, k; // 用于循环迭代的变量，register关键字表示将它们存储在CPU寄存器中以提高访问速度
    int kk;               // 索引变量，用于二维数组var的访问
    int n;                // 当前分层中的观测数
    int ngroup;           // 组的数量
    int ntot;             // 总观测数
    int istart;           // 当前分层的起始索引
    int koff;             // 当前分层组的偏移量，用于obs和exp数组的访问
    double km;            // Kaplan-Meier估计值
    double nrisk;         // 当前风险集中的样本数
    double wt;            // 权重值，根据Kaplan-Meier估计和rho计算得出
    double tmp;           // 临时变量，用于方差计算
    double deaths;        // 当前时间点的死亡数
    double temp_i;          // 观测值-期望值


    // 初始化变量
    ntot = *nn;             // 总样本数
    ngroup = *nngroup;      // 组数
    istart = 0; koff = 0;   // 初始开始索引和偏移量
    // double temp_i = 0;



    for (i = 0; i < *nstrat * ngroup; i++) {

        obs[i] = 0;
        exp[i] = 0;
        tempt[i] = 0;
    }


    // 循环遍历每个分层
        // 初始化风险数组为0
        for (i = 0; i < ngroup; i++) risk[i] = 0;

        // 找到当前分层的最后一个观测值
        // for (i = istart; i < ntot; i++)
        //     if (strata[i] == 1) break;
        // n = i + 1;  //n = 228
        n = ntot;  //n = 228


        // 进行实际检验
        for (i = n - 1; i >= istart; i--) {  
            // 计算权重
			wt = 1;

            // 计算死亡数
            deaths = 0;
            for (j = i; j >= istart && time[j] == time[i]; j--) {
                k = group[j] - 1;
                deaths += status[j];
                risk[k] += 1;
                obs[k + koff] += status[j] * wt;
                // obs[k + koff] = status[j] * wt;
            }

            i = j + 1;    

            nrisk = n - i; 


            // 如果有死亡事件，计算期望值和方差
			// printf("deaths:%f\n",deaths);
            if (deaths > 0) {
                for (k = 0; k < ngroup; k++)
				{

                    exp[k + koff] += wt * deaths * risk[k] / nrisk;
                    // exp[k + koff] = wt * deaths * risk[k] / nrisk;
                    temp_i = wt * (obs[k + koff]-exp[k + koff]);
					
                    tempt[k + koff] += temp_i;
					// printf("exp[%d]:%f\n",k,exp[k + koff]);
					// printf("obs[%d]:%f\n",k,obs[k + koff]);
					// printf("tempt[%d]:%f\n",k,tempt[k + koff]);
					
				}


                if (nrisk == 1) continue;  // 只有一个样本，无需计算方差
                kk = 0;
                wt = wt * wt;
                for (j = 0; j < ngroup; j++) {
                    tmp = wt * deaths * risk[j] * (nrisk - deaths) / (nrisk * (nrisk - 1));
                    var[kk + j] += tmp;
                    for (k = 0; k < ngroup; k++) {
                        var[kk] -= tmp * risk[k] / nrisk;
                        kk++;
                    }
                }
            }


        }
}
