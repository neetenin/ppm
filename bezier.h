#ifndef BEZIER
#define BEZIER

#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <vector>
#include "CGMath.hpp"

using namespace std;

class BezierCurve3D {
public:
    BezierCurve3D();
    BezierCurve3D(Vec3* p, int dg);// p,dg分别表示控制点和曲线次数
    void SetCurve(Vec3* p, int dg);
    Vec3 Point(double t);          // 用De Casteljau算法计算曲线上的点
    void gen(FILE* output, int d1, int d2);
    ~BezierCurve3D();
private:
    Vec3* cp;
    int degree; // degree of the curves
};

#endif