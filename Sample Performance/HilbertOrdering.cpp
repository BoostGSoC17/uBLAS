#include <bits/stdc++.h>

using namespace std;

void Hilbert(int startx, int starty, int endx, int endy, int orientation) {
    if(startx == endx && starty == endy) {
        cout << startx << " " << starty << "\n";
        return;
    }
    int midx = (startx+endx)>>1;    // as it is for square matrix (2^k x 2^k)
    int midy = (starty+endy)>>1;
    if(orientation == 1) {
        Hilbert(midx+1, starty, endx, midy, 2);
        Hilbert(startx, starty, midx, midy, 1);
        Hilbert(startx, midy+1, midx, endy, 1);
        Hilbert(midx+1, midy+1, endx, endy, 3);
    }
    else if(orientation == 2) {
        Hilbert(midx+1, starty, endx, midy, 1);
        Hilbert(midx+1, midy+1, endx, endy, 2);
        Hilbert(startx, midy+1, midx, endy, 2);
        Hilbert(startx, starty, midx, midy, 4);
    }
    else if(orientation == 3) {
        Hilbert(startx, midy+1, midx, endy, 4);
        Hilbert(startx, starty, midx, midy, 3);
        Hilbert(midx+1, starty, endx, midy, 3);
        Hilbert(midx+1, midy+1, endx, endy, 1);
    }
    else {
        Hilbert(startx, midy+1, midx, endy, 3);
        Hilbert(midx+1, midy+1, endx, endy, 4);
        Hilbert(midx+1, starty, endx, midy, 4);
        Hilbert(startx, starty, midx, midy, 2);
    }
}

int main() {

    //freopen("A.txt", "w", stdout);
    Hilbert(1, 1, 8, 8, 1);
    return 0;
}
