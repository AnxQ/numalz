//
// Created by 安啸琪 on 2019/12/23.
//

#include <iostream>
#include "../Matrix.h"
#include "../Solver.h"

using namespace std;

int main() {
    Mat A({
                  {2, 1, 1, 0},
                  {4, 3, 3, 1},
                  {8, 7, 9, 5},
                  {6, 7, 9, 8}
          });
    Mat b ({1, 2, 2, -1});
    cout << hstack(A, b);
    auto g_s = GaussSolver(A, b);
    cout << g_s.solve();
    auto gcp_s = GaussColPivotSolver(A, b);
    cout << gcp_s.solve();
    auto j_is = JacobiIterSolver(A, b);
    cout << j_is.solve();
    auto gs_is = GaussSeidelIterSolver(A, b);
    cout << gs_is.solve();
    auto sor_is = SORIterSolver(A, b, 1.3);
    cout << sor_is.solve();
    return 0;
}