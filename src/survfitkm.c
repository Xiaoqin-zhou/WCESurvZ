/*
** 处理原始数据的survfit例程，不适用于多状态数据，也不适用于Cox模型
** 参考文档：survival:survfitkm
**
** 传递给R的变量名使用"zed2"，在代码中变量名为"zed"
** y: 生存时间的2或3列矩阵
** weight: 案例权重的向量
** sort1: 起始时间的排序索引向量（如果y有3列），否则为null
** sort2: 结束时间的排序索引向量
** type: 生存计算类型
** id: 用于聚类方差的组标识符。0,1,2,...
** nid: 唯一组的数量
** position: 用于标记多次观察，例如(1,10)，(10,14)，(14,20)。
**          1*(第一个观察) + 2*(最后一个观察)。
**          只有一个观察值的个体的值为3。
** influence: 返回哪些影响矩阵，1*(累积风险) + 2*(生存率)
** reverse: 0=普通KM，1=估计删失分布（未完成）
** entry: 1=也保留具有唯一进入时间的行
*/
#include <math.h>
#include "survS.h"
#include "survproto.h"

SEXP survfitkm(SEXP y2, SEXP weight2, SEXP sort12, SEXP sort22, 
               SEXP type2, SEXP id2, SEXP nid2, SEXP position2,
               SEXP influence2, SEXP reverse2, SEXP entry2) {
              
    // 声明变量
    int i, i1, i2, j, k, person1, person2; // 常规索引变量
    int nused, nid, type, influence; // 使用的观测数、组数、生存分析类型、影响力变量
    int ny, ntime; // 输入数据的列数和唯一时间点的数量
    int reverse, entry; // 是否反转和是否记录进入时间的标志
    double* time1 = 0, * time2, * status, * wt; // 时间、状态和权重指针
    double v1, v2, dtemp, haz; // 中间变量和风险变量
    double temp, dtemp2, dtemp3, frac, btemp; // 临时变量和分数变量
    double d0, d1, nrisk; // 死亡数、权重死亡数、风险数
    int* sort1 = 0, * sort2, * id = 0; // 排序索引和组标识符指针
    static const char *outnames[]={"time", "n", "estimate", "std.err",
                                     "influence1", "influence2", ""};  // 输出变量名称
    SEXP rlist; // SEXP类型的R列表
    double* gwt = 0, * inf1 = 0, * inf2 = 0, * inf3 = 0; // 工作向量和影响矩阵指针
    int* gcount = 0; // 组计数器指针
    int n1, n2, n3, n4; // 风险数、事件数、删失数、添加数
    int* position = 0; // 位置指针
    double wt1, wt2, wt3, wt4; // 风险数、事件数、删失数、添加数的权重版本
                      
    /* 输出变量 */
    double  *n[8],  *dtime,  // n数组和唯一时间指针
            *kvec, *nvec, *std[3], *imat1=0, *imat2=0; // 当前估计值、标准差和影响矩阵指针 /* =0 以消除 -Wall 警告 */
    double km, nelson; // KM估计值和Nelson-Aalen估计值 /* 当前估计值 */

    /* 映射输入数据 */
    ny = ncols(y2);     /* 获取 y2 的列数，如果是2表示普通生存率数据，如果是3表示起始和停止数据 */
    nused = nrows(y2);  /* 获取 y2 的行数，即使用的数据量 */
    
    //if (ny == 3) {
    //    time1 = REAL(y2);        /* 将 y2 映射为 C 语言中的 double 类型指针，指向起始时间列 */
    //    time2 = time1 + nused;   /* time2 指向停止时间列，位置在 time1 之后 nused 个元素 */
    //    sort1 = INTEGER(sort12); /* 将 sort12 映射为 C 语言中的整数类型指针 */
    //}
    //else {
        time2 = REAL(y2);        /* 如果 y2 只有两列，则将 y2 映射为 C 语言中的 double 类型指针，指向时间列 */
        
    //}

    // 获取 time2 的长度
    //int len = XLENGTH(y2);
    //printf("%ld\n", len);
    // 打印 time2 的每一个元素
    //for (int i = 0; i < len; i++) {
    //    printf("time2[%ld] = %f\n", i, time2[i]);
    //}
    

    status = time2 + nused;      /* status 指向状态列，位置在 time2 之后 nused 个元素 */
    // printf("%d", nused);
    // printf("%d", status);
    wt = REAL(weight2);          /* 将 weight2 映射为 C 语言中的 double 类型指针，指向权重列 */
    sort2 = INTEGER(sort22);     /* 将 sort22 映射为 C 语言中的整数类型指针 */
    nused = LENGTH(sort22);      /* 获取 sort22 的长度，即实际使用的数据量 */

    type = asInteger(type2);     /* 将 type2 转换为整数，表示生存分析的类型 */
    nid = asInteger(nid2);       /* 将 nid2 转换为整数，表示唯一组的数量 */
    // if (nid > 0) {
    //     id = INTEGER(id2);       /* 如果 nid 大于 0，则将 id2 映射为 C 语言中的整数类型指针，指向组标识符 */
    // }
    position = INTEGER(position2); /* 将 position2 映射为 C 语言中的整数类型指针，指向位置标识符 */
    influence = asInteger(influence2); /* 将 influence2 转换为整数，表示影响矩阵的类型 */
    reverse = asInteger(reverse2);     /* 将 reverse2 转换为整数，表示是否反转 */
    entry = asInteger(entry2);         /* 将 entry2 转换为整数，表示是否记录进入时间 */


    /* nused 用于两个目的。第一个是输入数据y的长度，仅用于设置 time1, time2 和 status。
       第二个是我们实际使用的观察数，即 sort2 的长度。
       本例程可以多次调用，sort1/sort2 指向数据的不同子集，而 y, wt, id 和 position 保持不变。
    */

    /* 第一步，获取唯一时间的数量，用于内存分配。
       提供了 xval 组（唯一 id 值）的数量。
       我同时计算唯一的起始时间和结束时间。
    */
    // 如果 entry 标志为 1，表示需要记录进入时间
    //if (entry == 1) {
    //    ntime = 1; // 初始化唯一时间点数量为 1
    //    temp = time1[sort1[0]];  /* 将 temp 设置为第一个排序的起始时间，比任何事件/删失时间都短 */
    //    j = 1; // 初始化第二个排序索引为 1

    //    // 遍历所有使用的数据点
    //    for (i = 0; i < nused; i++) {
    //        i2 = sort2[i]; // 获取排序后的第 i 个结束时间的索引

    //        // 查找所有起始时间小于当前结束时间的情况
    //        for (; j < nused && (time1[sort1[j]] < time2[i2]); j++) {
    //            i1 = sort1[j]; // 获取排序后的第 j 个起始时间的索引

    //            // 如果当前起始时间与 temp 不同，且位置标志表示这是一个有效的起始时间
    //            if (time1[i1] != temp && (position[i1] & 1) == 1) {
    //                ntime++; // 增加唯一时间点的计数
    //                temp = time1[i1]; // 更新 temp 为当前起始时间
    //            }
    //        }

    //        // 如果当前结束时间与 temp 不同，且位置标志表示这是一个有效的结束时间或状态非零（表示事件发生）
    //        if ((time2[i2] != temp) && (position[i2] > 1 || status[i2] > 0)) {
    //            ntime++; // 增加唯一时间点的计数
    //            temp = time2[i2]; // 更新 temp 为当前结束时间
    //        }
    //    }
    //}
    //else {
        // 计算唯一的 time2 值（即结束时间），不考虑 entry 标志
        ntime = 1; // 初始化唯一时间点数量为 1
        temp = time2[sort2[0]];  /* 将 temp 设置为第一个排序的结束时间，最短的结束时间 */

        // 遍历所有使用的数据点，从第二个点开始，这一步用于得到唯一时间点，用于构建时间间隔或者计算生存概率
        for (i = 1; i < nused; i++) {
            i2 = sort2[i]; // 获取排序后的第 i 个结束时间的索引

            // 如果当前结束时间与 temp 不同，且位置标志表示这是一个有效的结束时间或状态非零（表示事件发生）
            if ((position[i2] > 1 || status[i2] > 0) && time2[i2] != temp) {
                ntime++; // 增加唯一时间点的计数
                temp = time2[i2]; // 更新 temp 为当前结束时间
            }
        }
    //}

 
    /* 为输出分配内存
        n 有 6 列，分别表示风险数、事件数、删失数，然后是
        3 个加权版本，然后可选两列用于风险集合添加数（当 entry=1 时）
    */
    // PROTECT 将 R 对象 rlist 保护起来，避免在垃圾回收过程中被清除
    PROTECT(rlist = mkNamed(VECSXP, outnames));

    // 为 dtime 分配内存，并将其添加到 rlist 中
    dtime = REAL(SET_VECTOR_ELT(rlist, 0, allocVector(REALSXP, ntime)));
    j = 6;
    // // 根据 entry 的值确定 j 的值。如果 entry 为 0，则 j 为 6；否则 j 为 8
    // if (entry == 0) {
        
    // }
    // else {
    //     j = 8;
    // }

    // 为 n[0] 分配一个 ntime * j 的矩阵，并将其添加到 rlist 中
    n[0] = REAL(SET_VECTOR_ELT(rlist, 1, allocMatrix(REALSXP, ntime, j)));

    // 设置 n 数组中每个元素的指针
    for (i = 1; i < j; i++) {
        n[i] = n[0] + i * ntime;
    }

    // 为 kvec 分配一个 ntime * 2 的矩阵，并将其添加到 rlist 中
    kvec = REAL(SET_VECTOR_ELT(rlist, 2, allocMatrix(REALSXP, ntime, 2)));

    // nvec 指向 kvec 矩阵的第二列，表示 Nelson-Aalen 估计值
    nvec = kvec + ntime;

    // 为 std 分配一个 ntime * ny 的矩阵，并将其添加到 rlist 中
    // std[0] 表示生存率的标准差，std[1] 表示累积风险的标准差，std[2] 表示 AUC 的标准差
    std[0] = REAL(SET_VECTOR_ELT(rlist, 3, allocMatrix(REALSXP, ntime, ny)));
    std[1] = std[0] + ntime;
    // if (ny == 3) {
    //     std[2] = std[0] + ntime;
    // }

    // 如果 nid 大于 0，表示需要计算 robust 方差
    // if (nid > 0) {
    //     // 分配内存给 gcount 数组，用于存储每个组的计数
    //     gcount = (int*)R_alloc(nid, sizeof(int));

    //     // 如果 type 小于 3，表示需要分配影响的工作向量
    //     if (type < 3) {
    //         // 分配 4 * nid 大小的内存给 gwt 数组，用于存储权重
    //         gwt = (double*)R_alloc(4 * nid, sizeof(double));
    //         inf1 = gwt + nid;
    //         inf2 = inf1 + nid;
    //         inf3 = inf2 + nid;

    //         // 初始化 gwt、gcount、inf1、inf2 和 inf3 数组
    //         for (i = 0; i < nid; i++) {
    //             gwt[i] = 0.0;
    //             gcount[i] = 0;
    //             inf1[i] = 0; // 每个受试者的 IJ 生存率，累积风险，AUC
    //             inf2[i] = 0;
    //             inf3[i] = 0;
    //         }
    //     }
    //     else {
    //         // 分配 2 * nid 大小的内存给 gwt 数组，用于存储权重
    //         gwt = (double*)R_alloc(2 * nid, sizeof(double));
    //         inf2 = gwt + nid;

    //         // 初始化 gwt 和 inf2 数组
    //         for (i = 0; i < nid; i++) {
    //             gwt[i] = 0.0;
    //             gcount[i] = 0;
    //             inf2[i] = 0;
    //         }
    //     }

    //     // 如果 type 小于 3，分配影响矩阵的内存
    //     if (type < 3) {
    //         if (influence == 1 || influence == 3) {
    //             imat1 = REAL(SET_VECTOR_ELT(rlist, 4, allocMatrix(REALSXP, nid, ntime)));
    //         }
    //         if (influence == 2 || influence == 3) {
    //             imat2 = REAL(SET_VECTOR_ELT(rlist, 5, allocMatrix(REALSXP, nid, ntime)));
    //         }
    //     }
    //     else if (influence != 0) {
    //         imat2 = REAL(SET_VECTOR_ELT(rlist, 5, allocMatrix(REALSXP, nid, ntime)));
    //     }
    // }


    // R_CheckUserInterrupt();  /* 检查是否有 Ctrl-C 中断 */

    /* 
    ** 首先填写唯一时间列表，重现计算唯一时间数量的代码。
    ** 现在 dtime 已经分配。
    */
    // 如果 entry 标志为 1，表示需要记录进入时间
    //if (entry == 1) {
    //    temp = time1[sort1[0]];  /* 将 temp 设置为第一个排序的起始时间，比任何事件/删失时间都短 */
    //    dtime[0] = temp; // 将 dtime 数组的第一个元素设置为 temp
    //    k = 1; // 初始化 dtime 的索引 k 为 1
    //    j = 1; // 初始化第二个排序索引 j 为 1

    //    // 遍历所有使用的数据点
    //    for (i = 0; i < nused; i++) {
    //        i2 = sort2[i]; // 获取排序后的第 i 个结束时间的索引

    //        // 查找所有起始时间小于当前结束时间的情况
    //        for (; j < nused && (time1[sort1[j]] < time2[i2]); j++) {
    //            i1 = sort1[j]; // 获取排序后的第 j 个起始时间的索引

    //            // 如果当前起始时间与 temp 不同，且位置标志表示这是一个有效的起始时间
    //            if (time1[i1] != temp && (position[i1] & 1) == 1) {
    //                temp = time1[i1]; // 更新 temp 为当前起始时间
    //                dtime[k++] = temp; // 将当前起始时间存储到 dtime 数组中，并增加 k 的值
    //            }
    //        }

    //        // 如果当前结束时间与 temp 不同，且位置标志表示这是一个有效的结束时间或状态非零（表示事件发生）
    //        if ((time2[i2] != temp) && (position[i2] > 1 || status[i2] > 0)) {
    //            temp = time2[i2]; // 更新 temp 为当前结束时间
    //            dtime[k++] = temp; // 将当前结束时间存储到 dtime 数组中，并增加 k 的值
    //        }
    //    }
    //}
    //else {
        // 如果 entry 标志不为 1，仅计算结束时间
        temp = time2[sort2[0]];  /* 将 temp 设置为第一个排序的结束时间，最短的结束时间 */
        dtime[0] = temp; // 将 dtime 数组的第一个元素设置为 temp
        k = 1; // 初始化 dtime 的索引 k 为 1

        // 遍历所有使用的数据点，从第二个点开始,此步骤主要获取dtime(2组分别唯一的时间点)
        for (i = 1; i < nused; i++) {
            i2 = sort2[i]; // 获取排序后的第 i 个结束时间的索引

            // 如果当前结束时间与 temp 不同，且位置标志表示这是一个有效的结束时间或状态非零（表示事件发生）
            if ((position[i2] > 1 || status[i2] > 0) && time2[i2] != temp) {
                temp = time2[i2]; // 更新 temp 为当前结束时间
                dtime[k++] = temp; // 将当前结束时间存储到 dtime 数组中，并增加 k 的值
            }
        }
    //}
            // for (int i = 0; i < ntime; i++) {
            //    printf("%f\n",dtime[i]);
            // }

 
    /*
    ** 接下来计算所有计数
    ** 临时变量 n1 = 风险数，n2= 事件数，n3 = 删失数，
    **   n4 = 添加数。wt1, wt2, wt3, wt4 = n1 到 n4 的加权版本。
    **   将其存储在 n[][0-5] 中以返回给 R 例程。
    ** person1, person2 分别跟踪 sort1 和 sort2，
    **   i1 和 i2 同样。
    */
    // 初始化 person1 和 person2 为 nused-1，即最后一个数据点的索引
    person1 = nused - 1;
    person2 = nused - 1;

    // 初始化 n1（风险数）和 wt1（加权风险数）为 0
    n1 = 0;
    wt1 = 0;

    // printf("ntime:%d",ntime);   ntime是唯一时间点个数
    // printf("nused:%d",nused);   nused是组内的时间点的个数

    // 逆向遍历时间点，从最后一个时间点开始
    for (k = ntime - 1; k >= 0; k--) {
        // 初始化 n2（事件数）、n3（删失数）、wt2（加权事件数）和 wt3（加权删失数）为 0
        n2 = 0;
        n3 = 0;
        wt2 = 0;
        wt3 = 0;

        // 逆向遍历数据点，根据结束时间排序
        for (; person2 >= 0; person2--) {
            i2 = sort2[person2]; // 获取排序后的第 person2 个结束时间的索引

            // 如果当前结束时间小于当前时间点，跳出循环
            if (time2[i2] < dtime[k]) break;

            // 增加风险数
            n1++;
            wt1 += wt[i2];    // 增加加权风险数

            // 如果当前状态为事件
            if (status[i2] == 1) {
                n2++;            // 增加事件数
                wt2 += wt[i2];   // 增加加权事件数
            }
            else if (position[i2] & 2) {
                /*
                ** 如果这是字符串的最后一个(a,b](b,c](c,d].. 对于
                ** 个体(position[i2]=2 或 3)，则为“真实”删失
                */
                n3++;           // 增加删失数
                wt3 += wt[i2];  // 增加加权删失数
            }
        }

        // 如果 ny == 3，表示存在起始时间和结束时间
        // if (ny == 3) {  /* 遍历进入时间 */
        //     n4 = 0;
        //     wt4 = 0;

        //     // 逆向遍历数据点，根据起始时间排序
        //     for (; person1 >= 0; person1--) {
        //         i1 = sort1[person1]; // 获取排序后的第 person1 个起始时间的索引

        //         // 如果当前起始时间小于当前时间点，跳出循环
        //         if (time1[i1] < dtime[k]) break;

        //         // 进入时间 >= dtime，从风险集合中删除
        //         n1--;
        //         wt1 -= wt[i1]; // 减少加权风险数

        //         // 如果 entry 为 1 且位置标志表示这是一个有效的起始时间，且起始时间等于当前时间点
        //         if ((entry == 1) && (position[i1] & 1) && time1[i1] == dtime[k]) {
        //             // 增加进入时间数和加权进入时间数
        //             n4++;
        //             wt4 += wt[i1];
        //         }
        //     }

        //     // 如果 entry 为 1，存储进入时间数和加权进入时间数
        //     if (entry == 1) {
        //         n[6][k] = n4;
        //         n[7][k] = wt4;
        //     }
        // }

        // 将当前时间点的风险数、事件数、删失数和相应的加权值存储到 n 数组中
        n[0][k] = n1;
        n[1][k] = n2;
        n[2][k] = n3;
        n[3][k] = wt1;
        n[4][k] = wt2;
        n[5][k] = wt3;
    }

    // R_CheckUserInterrupt();  /* 检查是否有 Ctrl-C 中断 */


    /*
    ** 计算生存率和累积风险，简单方差
    ** 参考文档：survfitkm：估计值
    */
    // 初始化 Nelson-Aalen 风险和 KM 生存率
    nelson = 0.0;
    km = 1.0;

    // 初始化方差 v1 和 v2
    v1 = 0;
    v2 = 0;

    printf("reverse:%f",reverse);
    // 遍历所有时间点
    for (i = 0; i < ntime; i++) {
        // 如果 reverse 标志为 1，表示估计删失曲线 G
        // if (reverse == 1) {
        //     d0 = n[2][i]; /* 删失数 */
        //     d1 = n[5][i]; /* 加权删失数 */
        //     nrisk = n[3][i] - n[4][i]; /* 风险数 */
        // }
        // else { // 否则，处理死亡数
            d0 = n[1][i];    /* 死亡数，未加权 */
            d1 = n[4][i];    /* 死亡数，加权 */
            nrisk = n[3][i]; /* 风险数 */
        // }

        // 处理 Nelson-Aalen 风险（type 为 1 或 3）
        if (type == 1 || type == 3) {
            if (d0 > 0 && d1 > 0) {  /* 至少有一个事件，权重 > 0 */
                nelson += d1 / nrisk; // 更新 Nelson-Aalen 风险
                v2 += d1 / (nrisk * nrisk); // 更新方差 v2
            }
            nvec[i] = nelson; // 存储 Nelson-Aalen 估计值
            std[1][i] = sqrt(v2); // 存储标准差
        }
        //else { // 处理 Fleming 风险
        //    for (j = 0; j < d0; j++) {
        //        dtemp = nrisk - j * d1 / d0;
        //        nelson += d1 / (d0 * dtemp); // 更新 Nelson-Aalen 风险
        //        v2 += d1 / (d0 * dtemp * dtemp); // 更新方差 v2
        //    }
        //    kvec[i] = exp(-nelson); // 存储 Fleming 生存率
        //    nvec[i] = nelson; // 存储 Nelson-Aalen 估计值
        //    std[1][i] = sqrt(v2); // 存储标准差
        //}

        // 处理 KM 生存率（type 小于 3）
        if (type < 3) {
            if (d0 > 0 && d1 > 0) {  /* 至少有一个事件 */
                km *= (nrisk - d1) / nrisk; // 更新 KM 生存率
                v1 += d1 / (nrisk * (nrisk - d1)); /* Greenwood 公式计算方差 */
            }
            kvec[i] = km; // 存储 KM 生存率
            std[0][i] = sqrt(v1); // 存储方差
        }
        // printf("km:%f",km);
        // printf("v1:%f",sqrt(v1));
        //else {  /* 处理 exp 生存率 */
        //    kvec[i] = exp(-nvec[i]); // 存储 exp 生存率
        //    std[0][i] = std[1][i]; // exp 生存率的标准差与 Nelson-Aalen 风险相同
        //}
    }
// printf("nid:%f",nid);
    // 如果 nid 大于 0，表示需要计算 robust 方差
    // if (nid > 0) {
    //     /* 使用 IJ 估计替换上述方差，可选地保存影响值
    //     ** 参考文档：survfitkm：robust 方差
    //     */
    //     v1 = 0;
    //     v2 = 0;
    //     person2 = 0;

    //     if (ny == 3) {
    //         person1 = 0;
    //     }
    //     else {
    //         /* 一开始，每个人都在风险中 */
    //         for (i = 0; i < nused; i++) {
    //             i2 = id[sort2[i]]; // 获取当前受试者的 ID
    //             gcount[i2]++; // 增加该受试者的计数
    //             gwt[i2] += wt[sort2[i]]; // 增加该受试者的权重
    //         }
    //     }

    //     km = 1; /* 这会滞后 */
    //     if (type == 1) {  /* KM 生存率，NA 风险 */
    //         for (i = 0; i < ntime; i++) {
    //             /* 丢弃过时的 */
    //             for (; person2 < nused; person2++) {
    //                 i2 = sort2[person2];
    //                 if (time2[i2] >= dtime[i]) break; // 跳出循环，处理当前时间点
    //                 gcount[id[i2]]--; // 减少该受试者的计数
    //                 if (gcount[id[i2]] == 0) {
    //                     gwt[id[i2]] = 0; // 如果计数为 0，设置权重为 0
    //                 }
    //                 else {
    //                     gwt[id[i2]] -= wt[i2]; // 减少该受试者的权重
    //                 }
    //             }
    //             if (ny == 3) {
    //                 /* 添加新受试者 */
    //                 for (; person1 < nused; person1++) {
    //                     i1 = sort1[person1];
    //                     if (time1[i1] >= dtime[i]) break; // 跳出循环，处理当前时间点
    //                     gcount[id[i1]]++; // 增加该受试者的计数
    //                     gwt[id[i1]] += wt[i1]; // 增加该受试者的权重
    //                 }
    //             }

    //             if (n[1][i] > 0 && n[4][i] > 0) { /* 需要更新总数 */
    //                 haz = n[4][i] / n[3][i]; // 计算风险值
    //                 for (k = 0; k < nid; k++) {
    //                     inf1[k] = inf1[k] * (1.0 - haz) + gwt[k] * km * haz / n[3][i];
    //                     inf2[k] -= gwt[k] * haz / n[3][i];
    //                 }
    //                 for (; person2 < nused; person2++) {
    //                     i2 = sort2[person2];
    //                     if (time2[i2] > dtime[i]) break; // 跳出循环，处理当前时间点
    //                     if (status[i2] == 1) {
    //                         inf1[id[i2]] -= km * wt[i2] / n[3][i];
    //                         inf2[id[i2]] += wt[i2] / n[3][i];
    //                     }
    //                     gcount[id[i2]]--;
    //                     if (gcount[id[i2]] == 0) {
    //                         gwt[id[i2]] = 0.0;
    //                     }
    //                     else {
    //                         gwt[id[i2]] -= wt[i2];
    //                     }
    //                 }
    //                 km *= (1 - haz); // 更新 KM 生存率
    //                 nelson += haz; // 更新 Nelson-Aalen 风险值

    //                 v1 = 0;
    //                 v2 = 0;
    //                 for (k = 0; k < nid; k++) {
    //                     v1 += inf1[k] * inf1[k];
    //                     v2 += inf2[k] * inf2[k];
    //                 }
    //             }

    //             std[0][i] = sqrt(v1); // 存储标准差
    //             std[1][i] = sqrt(v2); // 存储标准差
    //             if (influence == 1 || influence == 3) {
    //                 for (k = 0; k < nid; k++) {
    //                     *imat1++ = inf1[k]; // 保存影响值
    //                 }
    //             }
    //             if (influence == 2 || influence == 3) {
    //                 for (k = 0; k < nid; k++) {
    //                     *imat2++ = inf2[k]; // 保存影响值
    //                 }
    //             }
    //         }
    //     }
    //     //else if (type == 2) {  /* KM 生存率，Fleming-Harrington 风险 */
    //     //    for (i = 0; i < ntime; i++) {
    //     //        /* 丢弃过时的 */
    //     //        for (; person2 < nused; person2++) {
    //     //            i2 = sort2[person2];
    //     //            if (time2[i2] >= dtime[i]) break; // 跳出循环，处理当前时间点
    //     //            gcount[id[i2]]--; // 减少该受试者的计数
    //     //            if (gcount[id[i2]] == 0) {
    //     //                gwt[id[i2]] = 0;
    //     //            }
    //     //            else {
    //     //                gwt[id[i2]] -= wt[i2]; // 减少该受试者的权重
    //     //            }
    //     //        }
    //     //        if (ny == 3) {
    //     //            /* 添加新受试者 */
    //     //            for (; person1 < nused; person1++) {
    //     //                i1 = sort1[person1];
    //     //                if (time1[i1] >= dtime[i]) break; // 跳出循环，处理当前时间点
    //     //                gcount[id[i1]]++;
    //     //                gwt[id[i1]] += wt[i1];
    //     //            }
    //     //        }

    //     //        if (n[1][i] > 0 && n[4][i] > 0) { /* 需要更新总数 */
    //     //            dtemp = 0;  /* 工作分母 */
    //     //            dtemp2 = 0;  /* 平方和 */
    //     //            dtemp3 = 0;
    //     //            temp = n[3][i] - n[4][i];  /* 非死亡的权重和 */
    //     //            for (k = n[1][i]; k > 0; k--) {
    //     //                frac = k / n[1][i];
    //     //                btemp = 1 / (temp + frac * n[4][i]);  /* 数学中的“b” */
    //     //                dtemp += btemp;
    //     //                dtemp2 += btemp * btemp * frac;
    //     //                dtemp3 += btemp * btemp;    /* 非死亡导数 */
    //     //            }

    //     //            dtemp /= n[1][i];        /* 平均分母 */
    //     //            if (n[4][i] != n[1][i]) { /* 权重情况 */
    //     //                dtemp2 *= n[4][i] / n[1][i];
    //     //                dtemp3 *= n[4][i] / n[1][i];
    //     //            }
    //     //            nelson += n[4][i] * dtemp;

    //     //            haz = n[4][i] / n[3][i];
    //     //            for (k = 0; k < nid; k++) {
    //     //                inf1[k] = inf1[k] * (1.0 - haz) + gwt[k] * km * haz / n[3][i];
    //     //                if (gcount[k] > 0) {
    //     //                    inf2[k] -= gwt[k] * dtemp3;
    //     //                }
    //     //            }
    //     //            for (; person2 < nused; person2++) {
    //     //                /* 跟踪到此事件时间点的终点 */
    //     //                i2 = sort2[person2];
    //     //                if (time2[i2] > dtime[i]) break; // 跳出循环，处理当前时间点
    //     //                if (status[i2] == 1) {
    //     //                    inf1[id[i2]] -= km * wt[i2] / n[3][i];
    //     //                    inf2[id[i2]] += wt[i2] * (dtemp + dtemp3 - dtemp2);
    //     //                }
    //     //                gcount[id[i2]]--;
    //     //                if (gcount[id[i2]] == 0) {
    //     //                    gwt[id[i2]] = 0.0;
    //     //                }
    //     //                else {
    //     //                    gwt[id[i2]] -= wt[i2];
    //     //                }
    //     //            }
    //     //            km *= (1 - haz);

    //     //            v1 = 0;
    //     //            v2 = 0;
    //     //            for (k = 0; k < nid; k++) {
    //     //                v1 += inf1[k] * inf1[k];
    //     //                v2 += inf2[k] * inf2[k];
    //     //            }
    //     //        }

    //     //        std[0][i] = sqrt(v1); // 存储标准差
    //     //        std[1][i] = sqrt(v2); // 存储标准差
    //     //        if (influence == 1 || influence == 3) {
    //     //            for (k = 0; k < nid; k++) {
    //     //                *imat1++ = inf1[k]; // 保存影响值
    //     //            }
    //     //        }
    //     //        if (influence == 2 || influence == 3) {
    //     //            for (k = 0; k < nid; k++) {
    //     //                *imat2++ = inf2[k]; // 保存影响值
    //     //            }
    //     //        }
    //     //    }
    //     //}

    //     //else if (type == 3) {  /* exp() 生存率，NA 风险 */
    //     //    for (i = 0; i < ntime; i++) {
    //     //        /* 丢弃过时的 */
    //     //        for (; person2 < nused; person2++) {
    //     //            i2 = sort2[person2];
    //     //            if (time2[i2] >= dtime[i]) break;
    //     //            gcount[id[i2]]--; // 减少该受试者的计数
    //     //            if (gcount[id[i2]] == 0) {
    //     //                gwt[id[i2]] = 0; // 如果计数为 0，设置权重为 0
    //     //            }
    //     //            else {
    //     //                gwt[id[i2]] -= wt[i2]; // 减少该受试者的权重
    //     //            }
    //     //        }
    //     //        if (ny == 3) {
    //     //            /* 添加新受试者 */
    //     //            for (; person1 < nused; person1++) {
    //     //                i1 = sort1[person1];
    //     //                if (time1[i1] >= dtime[i]) break;
    //     //                gcount[id[i1]]++; // 增加该受试者的计数
    //     //                gwt[id[i1]] += wt[i1]; // 增加该受试者的权重
    //     //            }
    //     //        }

    //         //    if (n[1][i] > 0 && n[4][i] > 0) { /* 需要更新总数 */
    //         //        haz = n[4][i] / n[3][i]; // 计算风险值
    //         //        for (k = 0; k < nid; k++) {
    //         //            inf2[k] -= gwt[k] * haz / n[3][i];
    //         //        }
    //         //        for (; person2 < nused; person2++) {
    //         //            /* 跟踪到此事件时间点的终点 */
    //         //            i2 = sort2[person2];
    //         //            if (time2[i2] > dtime[i]) break;
    //         //            if (status[i2] == 1) {
    //         //                inf2[id[i2]] += wt[i2] / n[3][i];
    //         //            }
    //         //            gcount[id[i2]]--;
    //         //            if (gcount[id[i2]] == 0) {
    //         //                gwt[id[i2]] = 0.0;
    //         //            }
    //         //            else {
    //         //                gwt[id[i2]] -= wt[i2];
    //         //            }
    //         //        }
    //         //        nelson += haz; // 更新 Nelson-Aalen 风险值

    //         //        v2 = 0;
    //         //        for (k = 0; k < nid; k++) {
    //         //            v2 += inf2[k] * inf2[k];
    //         //        }
    //         //    }

    //         //    std[1][i] = sqrt(v2); // 存储标准差
    //         //    std[0][i] = sqrt(v2); // 存储标准差

    //         //    if (influence > 0) {
    //         //        for (k = 0; k < nid; k++) {
    //         //            *imat2++ = inf2[k]; // 保存影响值
    //         //        }
    //         //    }
    //         //}

    //         //}
    //     //else {  /* exp() 生存率，Fleming-Harrington 风险 */
    //     //        for (i = 0; i < ntime; i++) {
    //     //            /* 丢弃过时的 */
    //     //            for (; person2 < nused; person2++) {
    //     //                i2 = sort2[person2];
    //     //                if (time2[i2] >= dtime[i]) break;
    //     //                gcount[id[i2]]--;
    //     //                if (gcount[id[i2]] == 0) {
    //     //                    gwt[id[i2]] = 0;
    //     //                }
    //     //                else {
    //     //                    gwt[id[i2]] -= wt[i2];
    //     //                }
    //     //            }
    //     //            if (ny == 3) {
    //     //                /* 添加新受试者 */
    //     //                for (; person1 < nused; person1++) {
    //     //                    i1 = sort1[person1];
    //     //                    if (time1[i1] >= dtime[i]) break;
    //     //                    gcount[id[i1]]++;
    //     //                    gwt[id[i1]] += wt[i1];
    //     //                }
    //     //            }
    //     //            if (n[1][i] > 0 && n[4][i] > 0) { /* 需要更新总数 */
    //     //                dtemp = 0;  /* 工作分母 */
    //     //                dtemp2 = 0;  /* 平方和 */
    //     //                dtemp3 = 0;
    //     //                temp = n[3][i] - n[4][i];  /* 非死亡的权重和 */
    //     //                for (k = n[1][i]; k > 0; k--) {
    //     //                    frac = k / n[1][i];
    //     //                    btemp = 1 / (temp + frac * n[4][i]);  /* 数学中的“b” */
    //     //                    dtemp += btemp;
    //     //                    dtemp2 += btemp * btemp * frac;
    //     //                    dtemp3 += btemp * btemp;    /* 非死亡导数 */
    //     //                }

    //     //                dtemp /= n[1][i];        /* 平均分母 */
    //     //                if (n[4][i] != n[1][i]) { /* 权重情况 */
    //     //                    dtemp2 *= n[4][i] / n[1][i];
    //     //                    dtemp3 *= n[4][i] / n[1][i];
    //     //                }
    //     //                nelson += n[4][i] * dtemp;

    //     //                for (k = 0; k < nid; k++) {
    //     //                    if (gcount[k] > 0) {
    //     //                        inf2[k] -= gwt[k] * dtemp3;
    //     //                    }
    //     //                }
    //     //                for (; person2 < nused; person2++) {
    //     //                    i2 = sort2[person2];
    //     //                    if (time2[i2] > dtime[i]) break;
    //     //                    if (status[i2] == 1) {
    //     //                        inf2[id[i2]] += wt[i2] * (dtemp + dtemp3 - dtemp2);
    //     //                    }
    //     //                    gcount[id[i2]]--;
    //     //                    if (gcount[id[i2]] == 0) {
    //     //                        gwt[id[i2]] = 0.0;
    //     //                    }
    //     //                    else {
    //     //                        gwt[id[i2]] -= wt[i2];
    //     //                    }
    //     //                }

    //     //                v2 = 0;
    //     //                for (k = 0; k < nid; k++) {
    //     //                    v2 += inf2[k] * inf2[k];
    //     //                }
    //     //            }

    //     //            std[1][i] = sqrt(v2); // 存储标准差
    //     //            std[0][i] = sqrt(v2); // 存储标准差

    //     //            if (influence > 0) {
    //     //                for (k = 0; k < nid; k++) {
    //     //                    *imat2++ = inf2[k]; // 保存影响值
    //     //                }
    //     //            }
    //     //        }
    //     //        }
    // }

    UNPROTECT(1); // 解除对 rlist 的保护
    return (rlist); // 返回结果列表
}

