/*
** 所有生存函数的原型
** 将其包含在每个例程中有助于防止参数错误
* 这些函数声明帮助防止函数定义与调用时参数不匹配的错误
*/

// 定义一个函数，用于将一维 double 数组转换为二维数组
double **dmatrix(double *array, int nrow, int ncol);

// 定义一个函数，用于将一维 int 数组转换为二维数组
int    **imatrix(int *array, int nrow, int ncol);

// 初始化循环的函数，接收最小值和最大值作为参数
void init_doloop(int min, int max);

// 执行循环的函数，接收循环次数和索引数组作为参数
int doloop(int nloops, int *index);

// 计算没有风险的个体数量的函数
// 接收参数：个体数、开始时间、结束时间、状态、排序索引1、排序索引2、分层信息
int *norisk(int n, double *time1, double *time2, double *status, 
            int *sort1, int *sort2, int *strata);

// 计算 P 值步长的函数
// 接收参数：列数、索引数组、索引数组2、权重数组、数据数组、因子数组、维度数组、切割数组、步长、边缘
double pystep(int nc, int *index, int *index2, double *wt, 
              double *data, int *fac, int *dims, double **cuts, 
              double step, int edge);

// 定义一个结构体 snode
typedef struct snode {
    int value;              // 节点的值
    int depth;              // 向前链接的数量
    struct snode *forward[1];  // 向前链接的集合
} snode;

// 定义 survfitkm 函数，用于计算生存曲线
// 这是一个计算生存曲线的函数，使用了 R 的 SEXP 类型作为参数，用于与 R 交互
// 接收多个 SEXP 类型的参数，表示 R 传递的数据
//这个函数接收多个参数，表示生存时间、权重、排序索引、类型、ID、组数、位置、影响、反向标记和入口标记等数据
SEXP survfitkm(SEXP y2, SEXP weight2, SEXP sort12, SEXP sort22, 
               SEXP type2, SEXP id2, SEXP nid2, SEXP position2, 
               SEXP influence2, SEXP reverse2, SEXP entry2);

void survdiff2(int   *nn,     int   *nngroup,    int   *nstrat, 
	       double *rho,    double *time,       int   *status, 
	       int   *group,  int   *strata,	   double *obs, 
	       double *exp, double *tempt,    double *var,        double *risk, 
	       double *kaplan);

void wcelogrank(int   *nn,     int   *nngroup,    
	          double *time,       int   *status, int   *group, double *weight, 
              double *nrisk1, double *nrisk2,  double *nriskall,
              double *obs, double *exp, double *tempt,   double *var,   double *risk);