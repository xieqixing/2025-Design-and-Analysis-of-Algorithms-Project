#include<iostream>
#include<string>
#include<utility>
#include<map>
#include<vector>
#include<algorithm>

#define MAX 3000
// seq2: 6 6 3 30
#define INCREASE_NUM 6
#define DECREASE_NUM 6
#define MUTATION_NUM 3
#define MIN_SUBSEQ 30

struct duplication{
    int start_idx, size, strand;

    bool operator<(const duplication& other) const{
        return this->size < other.size;
    }
};

struct score
{
    int score, tmp_score, max_position_in_reference, start_position, length;
    int increase_num, decrease_num1, decrease_num2;
    std::pair<int, int> parent;
}score_matrix[MAX][MAX];

struct dot
{
    int strand;
    int position_in_query, position_in_reference;
} dot_matrix[MAX][MAX];

struct sequence
{
    int ref_st, ref_en, query_st, query_en;
};

std::map<duplication, int> output;
std::string query, reference, reversed_reference;
std::vector<sequence> seq;

int max_score_idx;
int query_to_reference[MAX], judge_duplication[MAX];

void Get_reversed_reference(){
    reversed_reference = reference;
    for (size_t i = 0; i < reversed_reference.size(); i++)
    {
        switch (reversed_reference[i])  
        {
        case 'A':
            reversed_reference[i] = 'T';
            break;
        case 'T':
            reversed_reference[i] = 'A';
            break;
        case 'C':
            reversed_reference[i] = 'G';
            break;
        case 'G':
            reversed_reference[i] = 'C';
            break;
        default: 
            break;
        }
    }
}

void Setup_dot_matrix(){

    for(size_t i = 0; i < query.size(); i++){
        for(size_t j = 0; j < reference.size(); j++){
            if(query[i] == reference[j]){
                dot_matrix[i+1][j+1].strand = 1;
            }else if(query[i] == reversed_reference[j]){
                dot_matrix[i+1][j+1].strand = -1;
            }

            dot_matrix[i+1][j+1].position_in_query = i+1;
            dot_matrix[i+1][j+1].position_in_reference = j+1;
        }
    }
}

// void Setup_score_matrix(){
//     score_matrix[1][1].score = 1;
//     for(size_t i = 2; i <= query.size(); i++){

//         // initialize temp score idx
//         for(size_t j = reference.size(); j >= 1; j--) temp_max_score[j] = 0, temp_max_score_idx[j] = 0;
//         //max_score_idx[reference.size()+1] = temp_max_score_idx;

//         for(size_t j = reference.size(); j >= 1; j--){
//             // dot doesn't exist
//             if(dot_matrix[i][j].strand == 0) continue;
            
//             score_matrix[i][j].max_position_in_reference = j;
//             // match successfully when strand equals to 1
//             if(dot_matrix[i][j].strand == 1 && dot_matrix[i-1][j-1].strand == 1)
//             {
//                 score_matrix[i][j].score = score_matrix[i-1][j-1].score + 1;
//                 score_matrix[i][j].parent = std::make_pair(i-1, j-1);
//                 score_matrix[i][j].max_position_in_reference = std::max(score_matrix[i][j].max_position_in_reference, score_matrix[i-1][j-1].max_position_in_reference);
//             }
//             else if(dot_matrix[i][j].strand == -1 && dot_matrix[i-1][j+1].strand == -1) // match successfully when strand equals to -1
//             {
//                 score_matrix[i][j].score = score_matrix[i-1][j+1].score + 1;
//                 score_matrix[i][j].parent = std::make_pair(i-1, j+1);
//                 score_matrix[i][j].max_position_in_reference = std::max(score_matrix[i][j].max_position_in_reference, score_matrix[i-1][j+1].max_position_in_reference);
//             }

//             // duplication occurred
//             if(score_matrix[i-1][max_score_idx[j]].score > score_matrix[i][j].score)
//             {
//                 score_matrix[i][j].score = score_matrix[i-1][max_score_idx[j]].score;
//                 score_matrix[i][j].parent = std::make_pair(i-1, max_score_idx[j]);
//                 score_matrix[i][j].max_position_in_reference = std::max(score_matrix[i][j].max_position_in_reference, score_matrix[i-1][max_score_idx[j]].max_position_in_reference);
//             }

//             // update temp max score
//             if(score_matrix[i][j].score > temp_max_score[score_matrix[i][j].max_position_in_reference]) temp_max_score[score_matrix[i][j].max_position_in_reference] = score_matrix[i][j].score, temp_max_score_idx[score_matrix[i][j].max_position_in_reference] = j;

//         }

//         // initialize max score idx
//         for(size_t j = reference.size(); j >= 1; j--){
//             if(temp_max_score[j] > score_matrix[i][max_score_idx[j+1]].score) max_score_idx[j] = temp_max_score_idx[j];
//             else max_score_idx[j] = max_score_idx[j+1];
//         }
//     }
// }

// void Output_the_answer(){
//     std::cout << "The max score is: " << score_matrix[(int)query.size()][max_score_idx[1]].score << std::endl;
//     std::cout << "The length of query is: " << query.size() << std::endl;
//     std::cout << std::endl;

//     int max_idx = max_score_idx[1];
//     for (size_t i = query.size(); i >= 1; i--)
//     {
//         query_to_reference[i] = max_idx;
//         //std::cout << max_idx << ' ';
//         max_idx = score_matrix[i][max_idx].parent.second;
//         //std::cout<<dot_matrix[i][temp_max_score_idx].strand;
//     }
    
//     int size = 0, strand = 0, start_position = 0, end_position = 0;
//     for (size_t i = 1; i <= query.size(); i++)
//     {
//         if(judge_duplication[query_to_reference[i]] == 0){
//             judge_duplication[query_to_reference[i]]++;

//             if(size == 0){
//                 continue;
//             }
//             else
//             {
//                 // join the map and initialize
//                 duplication dup;

//                 if(strand == 1) dup = {end_position, size, strand};
//                 else dup = {start_position, size, strand};

//                 if(output.find(dup) == output.end()) output[dup] = 1;
//                 else output[dup]++;

//                 start_position = end_position = size = strand = 0;
//             }
//         }else{
//             if(start_position == 0){
//                 start_position = end_position = query_to_reference[i];
//                 size++;
//                 strand = dot_matrix[i][query_to_reference[i]].strand;
//             }else{
//                 if(end_position + strand == query_to_reference[i]){
//                     size++;
//                     end_position = query_to_reference[i];
//                 }else{
//                     //join the map and initialize
//                     duplication dup;

//                     if(strand == 1) dup = {end_position, size, strand};
//                     else dup = {start_position, size, strand};
    
//                     if(output.find(dup) == output.end()) output[dup] = 1;
//                     else output[dup]++;

//                     start_position = end_position = query_to_reference[i];
//                     size = 1;
//                     strand = dot_matrix[i][query_to_reference[i]].strand;
//                 }
//             }
//         }
//     }

//     if(start_position != 0){
//         duplication dup;

//         if(strand == 1) dup = {end_position, size, strand};
//         else dup = {start_position, size, strand};

//         if(output.find(dup) == output.end()) output[dup] = 1;
//         else output[dup]++;
//     }
    
//     for (auto i = output.begin(); i != output.end(); i++)
//     {
//         std::cout << "Position in reference: " << i->first.start_idx <<std::endl;
//         std::cout << "Size: " << i->first.size <<std::endl;
//         std::cout << "Duplication times: " << i->second <<std::endl;
//         std::cout << "Strand: " << i->first.strand <<std::endl;
//         std::cout << std::endl;
//     }
// }

void Setup_score_matrix_2(){
    score_matrix[1][1].tmp_score = 1;
    score_matrix[1][1].start_position = 1;
    score_matrix[1][1].length = 1;
    max_score_idx = 1;


    for(size_t i = 1; i <= query.size(); i++){
        for(size_t j = 1; j <= reference.size(); j++){

            // 正向匹配
            if(dot_matrix[i][j].strand == 1){
                //基因突变的情况
                for(size_t k = 1; k <= MUTATION_NUM; k++){
                    if(i - k < 1 || j - k < 1) break;
                    
                    if(dot_matrix[i-k][j-k].strand == 1){
                        int tmp_score = score_matrix[i-k][j-k].score + score_matrix[i-k][j-k].tmp_score + 1 - 0.5 * (k - 1); 
                        if(tmp_score > score_matrix[i][j].score + score_matrix[i][j].tmp_score){
                            score_matrix[i][j].tmp_score = tmp_score - score_matrix[i-k][j-k].score;
                            score_matrix[i][j].score = score_matrix[i-k][j-k].score;
                            score_matrix[i][j].parent = std::make_pair(i-k, j-k);
                            score_matrix[i][j].length = score_matrix[i-k][j-k].length + 1;
                            score_matrix[i][j].decrease_num1 = score_matrix[i-k][j-k].decrease_num1;
                            score_matrix[i][j].increase_num = score_matrix[i-k][j-k].increase_num;
                        }
                    }
                }

                // 删除变异的情况
                for(size_t k = 1; k <= DECREASE_NUM; k++){
                    if(j - k < 1) break;
                    if(score_matrix[i-1][j-k].decrease_num1 + (k - 1) > DECREASE_NUM) continue; //删除的数量不能超过DECREASE_NUM
                    
                    if(dot_matrix[i-1][j-k].strand == 1){
                        int tmp_score = score_matrix[i-1][j-k].score + score_matrix[i-1][j-k].tmp_score + 1 - (k - 1); // 删除变异需要减去相应的分数
                        if(tmp_score > score_matrix[i][j].score + score_matrix[i][j].tmp_score){
                            score_matrix[i][j].tmp_score = tmp_score - score_matrix[i-1][j-k].score;
                            score_matrix[i][j].score = score_matrix[i-1][j-k].score;
                            score_matrix[i][j].parent = std::make_pair(i-1, j-k);
                            score_matrix[i][j].length = score_matrix[i-1][j-k].length + 1;
                            score_matrix[i][j].decrease_num1 = score_matrix[i-1][j-k].decrease_num1 + (k - 1);
                            score_matrix[i][j].increase_num = score_matrix[i-1][j-k].increase_num;
                        }
                    }
                }

                // 重复的情况(会不会有冲突？)
                int tmp_score = (score_matrix[i-1][max_score_idx].length > MIN_SUBSEQ) ? score_matrix[i-1][max_score_idx].score + score_matrix[i-1][max_score_idx].tmp_score : score_matrix[i-1][max_score_idx].score;
                if(tmp_score > score_matrix[i][j].score + score_matrix[i][j].tmp_score){
                    score_matrix[i][j].score = tmp_score;
                    score_matrix[i][j].tmp_score = 0;
                    score_matrix[i][j].length = 1;
                    score_matrix[i][j].parent = std::make_pair(i-1, max_score_idx);
                    score_matrix[i][j].start_position = 1;
                }
            }

            // 反向匹配
            if(dot_matrix[i][j].strand == -1){
                //基因突变的情况
                for(size_t k = 1; k <= MUTATION_NUM; k++){
                    if(i - k < 1 || j + k > reference.size()) break;
                    
                    if(dot_matrix[i-k][j+k].strand == -1){
                        int tmp_score = score_matrix[i-k][j+k].score + score_matrix[i-k][j+k].tmp_score + 1 - 0.5 * (k - 1);
                        if(tmp_score > score_matrix[i][j].score + score_matrix[i][j].tmp_score){
                            score_matrix[i][j].tmp_score = tmp_score - score_matrix[i-k][j+k].score;
                            score_matrix[i][j].score = score_matrix[i-k][j+k].score;
                            score_matrix[i][j].parent = std::make_pair(i-k, j+k);
                            score_matrix[i][j].length = score_matrix[i-k][j+k].length + 1;
                            score_matrix[i][j].decrease_num1 = score_matrix[i-k][j+k].decrease_num1;
                            score_matrix[i][j].increase_num = score_matrix[i-k][j+k].increase_num;
                        }
                    }
                }

                // 删除变异的情况
                for(size_t k = 1; k <= DECREASE_NUM; k++){
                    if(j + k > reference.size()) break;
                    if(score_matrix[i-1][j+k].decrease_num2 + (k - 1) > DECREASE_NUM) continue; //删除的数量不能超过DECREASE_NUM

                    if(dot_matrix[i-1][j+k].strand == -1){
                        int tmp_score = score_matrix[i-1][j+k].score + score_matrix[i-1][j+k].tmp_score + 1 - (k - 1); // 删除变异需要减去相应的分数
                        if(tmp_score > score_matrix[i][j].score + score_matrix[i][j].tmp_score){
                            score_matrix[i][j].tmp_score = tmp_score - score_matrix[i-1][j+k].score;
                            score_matrix[i][j].score = score_matrix[i-1][j+k].score;
                            score_matrix[i][j].parent = std::make_pair(i-1, j+k);
                            score_matrix[i][j].length = score_matrix[i-1][j+k].length + 1;
                            score_matrix[i][j].decrease_num2 = score_matrix[i-1][j+k].decrease_num2 + (k - 1);
                            score_matrix[i][j].increase_num = score_matrix[i-1][j+k].increase_num;
                        }
                    }
                }

                // 重复的情况(会不会有冲突？)
                int tmp_score = (score_matrix[i-1][max_score_idx].length > MIN_SUBSEQ) ? score_matrix[i-1][max_score_idx].score + score_matrix[i-1][max_score_idx].tmp_score : score_matrix[i-1][max_score_idx].score;
                if(tmp_score > score_matrix[i][j].score + score_matrix[i][j].tmp_score){
                    score_matrix[i][j].score = tmp_score;
                    score_matrix[i][j].tmp_score = 0;
                    score_matrix[i][j].length = 1;
                    score_matrix[i][j].parent = std::make_pair(i-1, max_score_idx);
                    score_matrix[i][j].start_position = 1;
                }
            }

            // 不匹配情况
            if(dot_matrix[i][j].strand == 0){
                for(size_t k = 1; k <= INCREASE_NUM; k++){
                    if(i - k < 1) break;
                    if(score_matrix[i-k][j].increase_num + (k - 1) > INCREASE_NUM) continue; //增加的数量不能超过INCREASE_NUM

                    if(dot_matrix[i-k][j].strand == 1){
                        int tmp_score = score_matrix[i-k][j].score + score_matrix[i-k][j].tmp_score - (k - 1); // 增加变异需要减去相应的分数
                        if(tmp_score > score_matrix[i][j].score + score_matrix[i][j].tmp_score){
                            dot_matrix[i][j].strand = 1;
                            score_matrix[i][j].tmp_score = tmp_score - score_matrix[i-k][j].score;
                            score_matrix[i][j].score = score_matrix[i-k][j].score;
                            score_matrix[i][j].parent = std::make_pair(i-k, j);
                            score_matrix[i][j].length = score_matrix[i-k][j].length + 1;
                            score_matrix[i][j].increase_num = score_matrix[i-k][j].increase_num + (k - 1);
                            score_matrix[i][j].decrease_num1 = score_matrix[i-k][j].decrease_num1;
                            score_matrix[i][j].decrease_num2 = score_matrix[i-k][j].decrease_num2;
                        }
                    }else if(dot_matrix[i-k][j].strand == -1){
                        int tmp_score = score_matrix[i-k][j].score + score_matrix[i-k][j].tmp_score - (k - 1); // 增加变异需要减去相应的分数
                        if(tmp_score > score_matrix[i][j].score + score_matrix[i][j].tmp_score){
                            dot_matrix[i][j].strand = -1;
                            score_matrix[i][j].score = score_matrix[i-k][j].score;
                            score_matrix[i][j].tmp_score = tmp_score - score_matrix[i-k][j].score;
                            score_matrix[i][j].parent = std::make_pair(i-k, j);
                            score_matrix[i][j].length = score_matrix[i-k][j].length + 1;
                            score_matrix[i][j].increase_num = score_matrix[i-k][j].increase_num + (k - 1);
                            score_matrix[i][j].decrease_num1 = score_matrix[i-k][j].decrease_num1;
                            score_matrix[i][j].decrease_num2 = score_matrix[i-k][j].decrease_num2;
                        }
                    }
                }
            }
            
        }

        // 更新max_score_idx
        int max_score = 0;
        for(size_t j = 1; j <= reference.size(); j++){
            int tmp_score = (score_matrix[i][j].length > MIN_SUBSEQ) ? score_matrix[i][j].score + score_matrix[i][j].tmp_score : score_matrix[i][j].score;
            if(tmp_score > max_score){
                max_score = tmp_score;
                max_score_idx = j;
            }
        }
    }
}

void Output_the_answer2(){
    int query_st = 0, query_en = 0, ref_st = 0, ref_en = 0;
    std::cout << "The max score is: " << score_matrix[(int)query.size()][max_score_idx].score << std::endl;

    for(size_t i = query.size(); i >= 1;){
        //std::cout << i << ' ' << max_score_idx << std::endl;
        if(query_en == 0){
            query_en = i;
            ref_en = max_score_idx;
        }

        if(score_matrix[i][max_score_idx].start_position == 1){
            query_st = i;
            ref_st = max_score_idx;
            if(dot_matrix[i][max_score_idx].strand == -1) std::swap(ref_st, ref_en);
            seq.push_back({ref_st, ref_en, query_st, query_en});
            query_en = 0;
            ref_en = 0;
        }

        int tmp1 = score_matrix[i][max_score_idx].parent.first, tmp2 = score_matrix[i][max_score_idx].parent.second;

        i = tmp1;
        max_score_idx = tmp2;
    }

    std::cout << '[';
    for(auto& i : seq){
        std::cout << '(' << i.query_st << ',' << i.query_en << ',' << i.ref_st << ',' << i.ref_en << ')' << ',';
    }
    std::cout << ']';
}

int main(){
    std::cin >> reference >> query;

    std::cout << "The length of query is: " << query.size() << std::endl;
    std::cout << "The length of reference is: " << reference.size() << std::endl;
    // get reversed reference
    Get_reversed_reference();

    std::cout << "The length of reversed reference is: " << reversed_reference.size() << std::endl;
    std::cout << "Start setting up the dot matrix..." << std::endl;
    // setup dot matrix 
    Setup_dot_matrix();
    std::cout << "Finish setting up the dot matrix..." << std::endl;

    // // using dynamic programming to get a score matrix
    // Setup_score_matrix();

    // // output the answer
    // Output_the_answer();

    std::cout << "Start calculating the score matrix..." << std::endl;
    // using dynamic programming to get a score matrix
    Setup_score_matrix_2();

    std::cout << "Finish calculating the score matrix..." << std::endl;
    std::cout << "Start outputting the answer..." << std::endl;
    Output_the_answer2();
    
    return 0;
}