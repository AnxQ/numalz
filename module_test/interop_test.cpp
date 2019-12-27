//
// Created by 安啸琪 on 2019/12/20.
//

#include <iostream>
#include "../Interop.h"

using namespace std;

int main() {
    Callable obj(2);
    obj(10);
    cout << &obj;
    int size;
    cin >> size;
    int *a = new int[size];
    return 0;
}
