//
// Created by 安啸琪 on 2019/12/26.
//

#include <opencv2/opencv.hpp>
#include <iostream>

using namespace std;

int main() {
    auto img = cv::imread("../module_test/img.jpg");
    cv::imshow("Environment test", img);
    cv::waitKey();
    return 0;
}