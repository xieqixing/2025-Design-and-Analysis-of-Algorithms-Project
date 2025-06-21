## Lab2 report

### 实验目标

​	本次 Lab 的目标是实现复杂 DNA 序列的比对，需要同时兼顾算法的时空复杂度和最终结果准确性。要求分别在 PPT 给出的两组 query 和 reference 中，实现 query 和reference 的匹配，并输出具体的匹配情况。其中，匹配过程中可能会遇到：单核苷酸突变，插入突变，删除突变，重复突变，倒位突变，片段移位等。

​	代码已提交至github:https://github.com/xieqixing/2025-Design-and-Analysis-of-Algorithms-Project

### 实验过程

#### 1. 代码框架

​	这一次的实验中，我们可以大致沿用lab1的整体思路，主体为`solve()`函数，如下所示：

```c++
void solve(std::string& reference, std::string& query){
    // allocate memory
    score_t** score_matrix = get_score_matrix(query.size() + 1, reference.size() + 1);
    dot_t** dot_matrix = get_dot_matrix(query.size() + 1, reference.size() + 1);

    // get reversed reference
    std::string reversed_reference;
    Get_reversed_reference(reference, reversed_reference);

    // setup dot matrix 
    Setup_dot_matrix(dot_matrix, reference, query, reversed_reference);
  

    // using dynamic programming to get a score matrix
    Setup_score_matrix_2(score_matrix, dot_matrix, query, reference);

    // output the answer
    Output_the_answer2(score_matrix, dot_matrix, query, reference, reversed_reference);
    
    // free the memory
    clear_matrix(score_matrix, dot_matrix, query.size() + 1);
    
}
```

​	其中只有`Setup_score_matrix_2()`以及`Output_the_answer2()`两个函数有较大的改动。其他函数的功能和lab1差不多：`Get_reversed_reference()`获取倒位的reference序列;`Setup_dot_matrix()`构建了一个`n*m`的矩阵`dot_matrix`记录reference，query之间的关系；`Setup_score_matrix_2()`则从之前的`dot_matrix`中计算出得分最高的切割方式；`Output_the_answer2()`则负责输出结果。



#### 2.数据类型

​	在这次lab中定义了三种数据类型，如下所示：

```c++
// score_matrix中每个点。score和tmp_score记录当前点得分；start_position, length记录当前切割的长度以及是否是起始点；increase_num，decrease_num记录ref和query移多少位之后匹配；parent记录前一个碱基的位置
typedef struct score{
    int score, tmp_score, start_position, length;
    int increase_num, decrease_num1, decrease_num2;
    std::pair<int, int> parent;
} score_t;

// dot_matrix中每一个点。strand为1说明对应碱基相等；为-1说明对应碱基互补；为0说明对应碱基无关
typedef struct dot{
    int strand;
} dot_t;

// 用于结果的输出
struct sequence{
    int ref_st, ref_en, query_st, query_en;
};
```

​	

#### 3.实现细节

​	代码的核心即为`Setup_score_matrix_2()`，其核心步骤为：

```c++
for(size_t i = 1; i <= query.size(); i++){
    for(size_t j = 2; j <= reference.size(); j++){
        
        // 处理正向匹配的情况
        if(dot_matrix[i][j].strand == 1){
            // 处理基因突变的情况         
            // 处理删除变异的情况
            // 处理重复片段的情况
        }
        
        // 处理反向匹配的情况
        if(dot_matrix[i][j].strand == -1){
            // 处理基因突变的情况         
            // 处理删除变异的情况
            // 处理重复片段的情况
        }
        
        // 处理不匹配的情况
        if(dot_matrix[i][j].strand == 0){
            // 一般可以被考虑为是插入突变的情况
        }
    }
}
```

​	其中处理重复片段的情况在lab1中已经实现，而处理其他变异需要定义常量：`INCREASE_NUM, DECREASE_NUM, MUTATION_NUM`。实现如下：

```c++
//基因突变的情况
for(size_t k = 1; k <= MUTATION_NUM; k++){
    if(i - k < 1 || j - k < 1) break; // 边界条件检验
    
    if(dot_matrix[i-k][j-k].strand == 1){
        // 计算tmp_score
        if(tmp_score > score_matrix[i][j].score + score_matrix[i][j].tmp_score){
            // 更新score_matrix[i][j]
        }
    }
}

// 删除变异的情况
for(size_t k = 1; k <= DECREASE_NUM; k++){
    if(j - k < 1) break; // 边界条件检测
    if(score_matrix[i-1][j-k].decrease_num1 + (k - 1) > DECREASE_NUM) continue; //删除的数量不能超过DECREASE_NUM
    
    if(dot_matrix[i-1][j-k].strand == 1){
        // 计算tmp_score
        if(tmp_score > score_matrix[i][j].score + score_matrix[i][j].tmp_score){
            // 更新score_matrix[i][j]
        }
    }
}
```

​	经过这些处理，程序就能处理多种多样的变异了，实验结果也已经提交至评分网站，第一组序列得分为29817，第二组序列得分为2103，基本已完成基线的要求。



#### 4.时间复杂度

​	构建`dot_matrix`的时间复杂度为：`O(mn)`，构建`score_matrix`的过程中由于`INCREASE_NUM, DECREASE_NUM, MUTATION_NUM`都是常数，所以时间复杂度依旧是`O(mn)`，而输出的过程中，时间复杂度为：`O(m)`。

​	所以总的时间复杂度为`O(mn)`。







