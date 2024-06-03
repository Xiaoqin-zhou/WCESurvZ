/*
** 这个文件会预加载 .C 例程的入口点
** 根据 R-core 的要求添加
** 它通过声明参数数量增加了一层保护，
** 并且可能带来一点点速度上的提升
*/
#include "survS.h" //包含自定义的头文件
#include "survproto.h" //包含自定义的头文件
#include "R_ext/Rdynload.h" //包含 R 动态加载的头文件
#include "Rversion.h" //包含 R 版本信息的头文件


static const R_CMethodDef Centries[] = {
    {"Csurvdiff2",  (DL_FUNC) &survdiff2, 14},
    {"Cwcelogrank",  (DL_FUNC) &wcelogrank, 14},
    {NULL, NULL, 0}
};

// 静态常量 R_CallMethodDef 数组，定义了要注册的 .C 例程
//这是一个静态常量数组，类型为 R_CallMethodDef，用于定义要注册的.Call 例程
static const R_CallMethodDef Callentries[] = {
    {"Csurvfitkm",    (DL_FUNC) &survfitkm,   11},  // 定义 "Csurvfitkm" 函数，指向 survfitkm，接受 11 个参数
    {NULL, NULL, 0}  // 数组的终止标记
};

// 初始化 WCESurvZ 包时调用的函数
void R_init_WCESurvZ(DllInfo *dll){
    // 注册 R 例程，第二个参数为 NULL，表示没有 .C 例程，第三个参数为 Callentries，表示要注册的 .Call 例程
    R_registerRoutines(dll, Centries,  Callentries, NULL, NULL);

    /* 
    ** 以下代码行使得只有上面定义的例程可以被外部包使用，
    ** 即像 dmatrix() 这样的内部函数现在是不可见的。
    */
    R_useDynamicSymbols(dll, FALSE); 

    /*
    ** 这一行使得它们只能通过上面定义的符号访问，
    ** 即 .Call("tmerge", ) 不会工作，但 .Call(Ctmerge, ) 可以。
    ** 这个特性在版本 2.16 中添加
    */
#if defined(R_VERSION) && R_VERSION >= R_Version(2, 16, 0)
    R_forceSymbols(dll, TRUE);
#endif
}
