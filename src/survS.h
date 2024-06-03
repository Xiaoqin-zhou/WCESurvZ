/* 
**  这个文件最初是为了支持 Splus，然后逐渐演变成支持 R 或 Splus（基于 ifdef 行），现在仅支持 R。
*/
#include "R.h"                 // 包含 R 的头文件，定义了基本的 R 接口
#include "Rinternals.h"        // 包含 R 内部结构的头文件
#include <R_ext/Utils.h>       // 包含 R 扩展工具的头文件

/*
** 使用 ALLOC 定义的内存由 S 自动移除。
**  使用 "Calloc" 定义的内存需要我自己移除。使用后者用于需要在调用之间持久存在的对象。
*/
#define ALLOC(a,b)  R_alloc(a,b)  // 定义宏 ALLOC，将其映射到 R 的内存分配函数 R_alloc
#define CALLOC(a,b) R_Calloc(a,b) // 定义宏 CALLOC，将其映射到 R 的内存分配函数 R_Calloc
#define FREE(a)     R_Free(a)     // 定义宏 FREE，将其映射到 R 的内存释放函数 R_Free
