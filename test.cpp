#include<iostream>
#include<string>

# define MAX 5010

using namespace std;

/*
求两个字符串的最长公共子串，但是两字串不能在相同的位置

特点：
1.重复子串是连续的
2.只会有一个断点吗

*/

int dp[MAX][MAX];
int ans;

string s1, s2;

int main(){

    cin>>s1>>s2;
    
    int sz1 = s1.size(), sz2 = s2.size();

    for (int i = 0; i < sz1; i++)
    {
        for (int j = i+1; j < sz2; j++)
        {
            if(s1[i] == s2[j]) {
                if(i == 0) dp[i][j] = 1;
                else dp[i][j] = min(dp[i-1][j-1]+1, j-i);
            }
            else dp[i][j] = 0;

            if(dp[i][j] == 400) cout<<i+1<<" "<<j+1<<endl;

            ans = max(ans, dp[i][j]);
        }

        // 分析一下对于每个i, 最大重复序列的数量

        /*
         for j = i+1, j < sz2 ans[dp[i][j]]++
         for j = 1, j < sz1 find MAX answer
         */
    }
    
    cout<<ans;

    return 0;
}