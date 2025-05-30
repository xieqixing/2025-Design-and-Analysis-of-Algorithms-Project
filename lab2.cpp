#include<iostream>
#include<string>
#include<utility>
#include<map>
#include<vector>
#include<algorithm>

#define MAX 30000
// seq2: 6 6 3 30
#define INCREASE_NUM 6
#define DECREASE_NUM 6
#define MUTATION_NUM 3
#define MIN_SUBSEQ 30

typedef struct score{
    int score, tmp_score, start_position, length;
    int increase_num, decrease_num1, decrease_num2;
    std::pair<int, int> parent;
} score_t;

typedef struct dot{
    int strand;
} dot_t;

struct sequence{
    int ref_st, ref_en, query_st, query_en;
};

int max_score_idx;

void Get_reversed_reference(std::string& reference, std::string& reversed_reference){
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

void Setup_dot_matrix(dot_t** dot_matrix, std::string& reference, std::string& query, std::string& reversed_reference){

    for(size_t i = 0; i < query.size(); i++){
        for(size_t j = 0; j < reference.size(); j++){
            if(query[i] == reference[j]){
                dot_matrix[i+1][j+1].strand = 1;
            }else if(query[i] == reversed_reference[j]){
                dot_matrix[i+1][j+1].strand = -1;
            }
        }
    }
}

void Setup_score_matrix_2(score_t** score_matrix, dot_t** dot_matrix, std::string& query, std::string& reference){
    // score_matrix[1][1].tmp_score = 1;
    // score_matrix[1][1].start_position = 1;
    // score_matrix[1][1].length = 1;
    // max_score_idx = 1;

    // 初始化
    for(size_t j = 1; j <= reference.size(); j++){
        if(dot_matrix[1][j].strand != 0){
            score_matrix[1][j].tmp_score = 1;
            score_matrix[1][j].score = 0;
            score_matrix[1][j].start_position = 1;
            score_matrix[1][j].length = 1;
        }
    }


    for(size_t i = 1; i <= query.size(); i++){
        for(size_t j = 2; j <= reference.size(); j++){

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
                    //score_matrix[i][j].decrease_num1 = 0;
                    //score_matrix[i][j].increase_num = 0;
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
                    //score_matrix[i][j].decrease_num2 = 0;
                    //score_matrix[i][j].increase_num = 0;
                    
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

void Output_the_answer2(score_t** score_matrix, dot_t** dot_matrix, std::string& query, std::string& reference, std::string& reversed_reference){
    std::vector<sequence> seq;
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
            if(dot_matrix[i][max_score_idx].strand == -1){
                std::swap(ref_st, ref_en);
                seq.push_back({ref_st - 1, ref_en, query_st, query_en + 1});
            }else{
                seq.push_back({ref_st, ref_en + 1, query_st, query_en + 1});
            }
            
            query_en = 0;
            ref_en = 0;
        }

        int tmp1 = score_matrix[i][max_score_idx].parent.first, tmp2 = score_matrix[i][max_score_idx].parent.second;

        i = tmp1;
        max_score_idx = tmp2;
    }

    // little check
    std::sort(seq.begin(), seq.end(), [](const sequence& a, const sequence& b) {
        return a.query_st < b.query_st;
    });

    std::vector<sequence> seq_ans;
    for(size_t i = 1; i < seq.size() - 1; i++){
        if(seq[i].query_en == seq[i].query_st + 1 && abs(seq[i-1].ref_en - seq[i+1].ref_st) <= 2){
            int query_st = seq[i-1].query_st, query_en = seq[i+1].query_en, ref_st = seq[i-1].ref_st, ref_en = seq[i+1].ref_en;
            seq_ans.push_back({ref_st, ref_en, query_st, query_en});
            //std::cout << '(' << ref_st << ',' << ref_en << ',' << query_st << ',' << query_en << ')' << std::endl;
            i += 2;
        }else{
            seq_ans.push_back(seq[i-1]);
        }
    }
    
    seq_ans.push_back(seq[seq.size()-2]);
    seq_ans.push_back(seq[seq.size()-1]);

    // output the answer
    std::cout << '[';
    for(auto& i : seq_ans){
        std::cout << '(' << i.query_st << ',' << i.query_en << ',' << i.ref_st << ',' << i.ref_en << ')' << ',';
    }
    std::cout << ']';
}

score_t** get_score_matrix(int idx1, int idx2){
    score_t** score_matrix = new score_t*[idx1];
    for(int i = 0; i < idx1; i++){
        score_matrix[i] = new score_t[idx2];
        for(int j = 0; j < idx2; j++){
            score_matrix[i][j].score = 0;
            score_matrix[i][j].tmp_score = 0;
            score_matrix[i][j].start_position = 0;
            score_matrix[i][j].length = 0;
            score_matrix[i][j].increase_num = 0;
            score_matrix[i][j].decrease_num1 = 0;
            score_matrix[i][j].decrease_num2 = 0;
        }
    }
    return score_matrix;
}

dot_t** get_dot_matrix(int idx1, int idx2){
    dot_t** dot_matrix = new dot_t*[idx1];
    for(int i = 0; i < idx1; i++){
        dot_matrix[i] = new dot_t[idx2];
        for(int j = 0; j < idx2; j++){
            dot_matrix[i][j].strand = 0;
        }
    }
    return dot_matrix;
}

void clear_matrix(score_t** score_matrix, dot_t** dot_matrix, int idx1){
    for(int i = 0; i < idx1; i++){
        delete[] score_matrix[i];
    }
    delete[] score_matrix;

    for(int i = 0; i < idx1; i++){
        delete[] dot_matrix[i];
    }
    delete[] dot_matrix;
}

void solve(std::string& reference, std::string& query){
    score_t** score_matrix = get_score_matrix(query.size() + 1, reference.size() + 1);
    dot_t** dot_matrix = get_dot_matrix(query.size() + 1, reference.size() + 1);

    // get reversed reference
    std::string reversed_reference;
    Get_reversed_reference(reference, reversed_reference);

    // setup dot matrix 
    std::cout << "Start setting up the dot matrix..." << std::endl;
    Setup_dot_matrix(dot_matrix, reference, query, reversed_reference);
    std::cout << "Finish setting up the dot matrix..." << std::endl;

    // using dynamic programming to get a score matrix
    std::cout << "Start calculating the score matrix..." << std::endl;
    Setup_score_matrix_2(score_matrix, dot_matrix, query, reference);
    std::cout << "Finish calculating the score matrix..." << std::endl;

    // output the answer
    std::cout << "Start outputting the answer..." << std::endl;
    Output_the_answer2(score_matrix, dot_matrix, query, reference, reversed_reference);
    clear_matrix(score_matrix, dot_matrix, query.size() + 1);
    std::cout << std::endl;
    std::cout << "Finish outputting the answer..." << std::endl;
}

void Get_Seq(score_t** score_matrix, dot_t** dot_matrix, std::string& query, std::vector<sequence>& seq, int bias){
    int query_st = 0, query_en = 0, ref_st = 0, ref_en = 0;

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
            seq.push_back({ref_st, ref_en, query_st + bias, query_en + bias});
            query_en = 0;
            ref_en = 0;
        }

        int tmp1 = score_matrix[i][max_score_idx].parent.first, tmp2 = score_matrix[i][max_score_idx].parent.second;

        i = tmp1;
        max_score_idx = tmp2;
    }
}

void Output_vector(std::vector<sequence>& seq){
    std::cout << '[';
    for(auto& i : seq){
        std::cout << '(' << i.query_st << ',' << i.query_en << ',' << i.ref_st << ',' << i.ref_en << ')' << ',';
    }
    std::cout << ']' << std::endl;
}

void solve_by_part(std::string& reference, std::string& query){
    std::vector<sequence> seq;
    int query_length = query.size(), reference_length = reference.size();
    int sub_query_length = 30000;

    for(int i = 1; sub_query_length * (i - 1) < query_length; i++){
        int query_st = sub_query_length * (i - 1);
        int query_en = std::min(sub_query_length * i, query_length);
        std::string sub_query = query.substr(query_st, query_en - query_st);

        // get the answer
        std::cout << "Processing the " << query_st << " to " << query_en << " part of the query..." << std::endl;
        score_t** score_matrix = get_score_matrix(sub_query.size() + 1, reference_length + 1);
        dot_t** dot_matrix = get_dot_matrix(sub_query.size() + 1, reference_length + 1);

        std::string reversed_reference;
        Get_reversed_reference(reference, reversed_reference);

        Setup_dot_matrix(dot_matrix, reference, sub_query, reversed_reference);
        Setup_score_matrix_2(score_matrix, dot_matrix, sub_query, reference);
        Get_Seq(score_matrix, dot_matrix, sub_query, seq, sub_query_length * (i - 1));
        clear_matrix(score_matrix, dot_matrix, sub_query.size() + 1);

    }

    Output_vector(seq);

}

int main(){
    std::string reference, query;
    std::cout << "Please input the reference and query sequence: " << std::endl;
    std::cin >> reference >> query;

    int query_length = query.size(), reference_length = reference.size();
    std::cout << "The length of query is: " << query_length << std::endl;
    std::cout << "The length of reference is: " << reference_length << std::endl;


    solve(reference, query);
    
    return 0;
}