#include <bits/stdc++.h>

#define vvi vector<vector<int> >
using namespace std;

struct Mat{
    int R, C;
    Mat(){
        R = C = 0;
    }
};

istream& operator>> (istream &in, Mat &obj){
    cin >> obj.R >> obj.C;
    return in;
}

int N; // number of Matrices;
vvi DP, splits;
vector<Mat> Dimension;

void chainMult(){
    for(int i=N; i>=1; i--){
        for(int j=i; j<=N; j++){
            DP[i][j] = 1e6;
            for(int k=i; k<j;k++){
                int tmp = DP[i][k] + DP[k+1][j] + Dimension[i].R * Dimension[k].C * Dimension[j].C;
                if(DP[i][j]>tmp){
                    DP[i][j] = tmp; splits[i][j] = k;
                }
            }
            if(i==j)
                DP[i][j] = 0;
        }
    }
}

void PrintChain(int i, int j) {
    if(i==j) {
        cout << "A" << i;
        return;
    }
    cout << "( ";
    PrintChain(i,splits[i][j]);
    PrintChain(splits[i][j]+1, j);
    cout << " ) ";
}

int main(){

    //freopen("Input.txt", "r", stdin);
    cin >> N;

    DP.resize(N+1,vector<int>(N+1,0));
    splits.resize(N+1, vector<int>(N+1,0));

    Dimension.resize(N+1,Mat());

    for(int i=1;i<=N;i++){
        cin >> Dimension[i];
        cout << Dimension[i].R << " " << Dimension[i].C << "\n";
    }

    chainMult();
    DP.clear();
    
    PrintChain(1, 6);
    cout << "\n";

    return 0;
}
