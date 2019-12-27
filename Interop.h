//
// Created by 安啸琪 on 2019/12/20.
//

#ifndef NUMALZ_INTEROP_H
#define NUMALZ_INTEROP_H

#include <functional>
#include <vector>

using namespace std;

class Callable : public function<int(int)> {
    int b;
    vector<int> X;
public:
    Callable(int b) : function([this](int a) {
        cout << this << endl;
        return a + this->b;
    }), b(b) {
        X.push_back(1);
    }
};

#endif //NUMALZ_INTEROP_H
