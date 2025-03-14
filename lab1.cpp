#include<iostream>
#include<string>
#include<utility>
#include<map>


#define MAX 2000

struct duplication{
    int start_idx, size, strand;

    bool operator<(const duplication& other) const{
        return this->size < other.size;
    }
};

struct score
{
    int score;
    std::pair<int, int> parent;
}score_matrix[MAX][MAX];


struct dot
{
    int strand;
    int position_in_query, position_in_reference;
} dot_matrix[MAX][MAX];

std::map<duplication, int> output;
std::string query, reference, reversed_reference;

int main(){
    std::cin >> reference >> query;

    // get reversed reference
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
    

    // setup dot matrix 
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
            
    // using dynamic programming to get a score matrix
    int max_score_idx = 1, temp_max_score, temp_max_score_idx;
    score_matrix[1][1].score = 1;

    for(size_t i = 2; i <= query.size(); i++){
        temp_max_score = 0;

        for(size_t j = reference.size(); j >= 1; j--){
            // dot doesn't exist
            if(dot_matrix[i][j].strand == 0) continue;

            // match successfully when strand equals to 1
            if(dot_matrix[i][j].strand == 1 && dot_matrix[i-1][j-1].strand == 1)
            {
                score_matrix[i][j].score = score_matrix[i-1][j-1].score + 1;
                score_matrix[i][j].parent = std::make_pair(i-1, j-1);

                if(score_matrix[i][j].score > temp_max_score) temp_max_score = score_matrix[i][j].score, temp_max_score_idx = j;
            }
            else if(dot_matrix[i][j].strand == -1 && dot_matrix[i-1][j+1].strand == -1) // match successfully when strand equals to -1
            {
                score_matrix[i][j].score = score_matrix[i-1][j+1].score + 1;
                score_matrix[i][j].parent = std::make_pair(i-1, j+1);

                if(score_matrix[i][j].score > temp_max_score) temp_max_score = score_matrix[i][j].score, temp_max_score_idx = j;
            }

            // duplication occurred
            if(score_matrix[i-1][max_score_idx].score > score_matrix[i][j].score)
            {
                score_matrix[i][j].score = score_matrix[i-1][max_score_idx].score;
                score_matrix[i][j].parent = std::make_pair(i-1, max_score_idx);

                if(score_matrix[i][j].score > temp_max_score) temp_max_score = score_matrix[i][j].score, temp_max_score_idx = j;
            }

        }

        max_score_idx = temp_max_score_idx;
    }

    // output the answer
    std::cout << "The max score is: " << temp_max_score << std::endl;
    std::cout << "The length of query is: " << query.size() << std::endl;
    
    
    int query_to_reference[MAX], judge_duplication[MAX] = {0};
    for (size_t i = query.size(); i >= 1; i--)
    {
        query_to_reference[i] = temp_max_score_idx;
        std::cout << temp_max_score_idx << ' ';
        temp_max_score_idx = score_matrix[i][temp_max_score_idx].parent.second;
        //std::cout<<dot_matrix[i][temp_max_score_idx].strand;
    }
    
    int size = 0, strand = 0, start_position = 0, end_position = 0;
    for (size_t i = 1; i <= query.size(); i++)
    {
        if(judge_duplication[query_to_reference[i]] == 0){
            judge_duplication[query_to_reference[i]]++;

            if(size == 0){
                continue;
            }
            else
            {
                // join the map and initialize
                duplication dup;

                if(strand == 1) dup = {end_position, size, strand};
                else dup = {start_position, size, strand};

                if(output.find(dup) == output.end()) output[dup] = 1;
                else output[dup]++;

                start_position = end_position = size = strand = 0;
            }
        }else{
            if(start_position == 0){
                start_position = end_position = query_to_reference[i];
                size++;
                strand = dot_matrix[i][query_to_reference[i]].strand;
            }else{
                if(end_position + strand == query_to_reference[i]){
                    size++;
                    end_position = query_to_reference[i];
                }else{
                    //join the map and initialize
                    duplication dup;

                    if(strand == 1) dup = {end_position, size, strand};
                    else dup = {start_position, size, strand};
    
                    if(output.find(dup) == output.end()) output[dup] = 1;
                    else output[dup]++;

                    start_position = end_position = query_to_reference[i];
                    size = 1;
                    strand = dot_matrix[i][query_to_reference[i]].strand;
                }
            }
        }
    }

    if(start_position != 0){
        duplication dup;

        if(strand == 1) dup = {end_position, size, strand};
        else dup = {start_position, size, strand};

        if(output.find(dup) == output.end()) output[dup] = 1;
        else output[dup]++;
    }
    
    for (auto i = output.begin(); i != output.end(); i++)
    {
        std::cout << "Position in reference: " << i->first.start_idx <<std::endl;
        std::cout << "Size: " << i->first.size <<std::endl;
        std::cout << "Duplication times: " << i->second <<std::endl;
        std::cout << "Strand: " << i->first.strand <<std::endl;
        std::cout << std::endl;
    }
    

    return 0;
}