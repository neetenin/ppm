#include "Bezier.h"

#define PI ((double)3.14159265358979)

BezierCurve3D::BezierCurve3D()
    : degree(0)
{
}

BezierCurve3D::~BezierCurve3D()
{
}

BezierCurve3D::BezierCurve3D(Vec3* p, int dg)
{
    degree = dg;
    cp = p;
}

void BezierCurve3D::SetCurve(Vec3* p, int dg)
{
    degree = dg;
    cp = p;
}

Vec3 BezierCurve3D::Point(double t)
{
    Vec3* cq = new Vec3[degree + 1];
    for (int i = 0; i <= degree; i++)
        cq[i] = cp[i];
    for (int k = 1; k <= degree; k++) {
        for (int i = 0; i <= degree - k; i++) {
            cq[i] = (1 - t) * cq[i] + t * cq[i + 1];
        }
    }
    return cq[0];
}

void BezierCurve3D::gen(FILE* output, int d1, int d2)
{
    printf("Generate\n");
    double dt = 1.0 / d1; // 参数
    double ds = 2 * PI / d2; // 角度
    double ct = 0.0;
    for (int i = 0; i <= d1; i++)
    {
        double cs = 0.0;
        Vec3 p = Point(ct);
        for (int j = 0; j <= d2; j++)
        {
            fprintf(output, "v %f %f %f\n", p.x * cos(cs), p.y, p.x * sin(cs));
            cs += ds;
        }
        ct += dt;
    }
}

int main()
{
    int cnt = 300;
    FILE *output;
    output = fopen("curve.obj", "w");
    Vec3 cp1[] = {Vec3(1.045, 44.52, 0), Vec3(0.42, 17.64, 0), Vec3(10.92, 18.375, 0), Vec3(20.16, 17.22, 0)};
    BezierCurve3D b1(cp1, 3);
    b1.gen(output, cnt, cnt);
    Vec3 cp2[] = {Vec3(20.16, 17.22, 0), Vec3(23.52, 16.8, 0), Vec3(20.16, 16.38, 0)};
    BezierCurve3D b2(cp2, 2);
    b2.gen(output, cnt, cnt);
    Vec3 cp3[] = {Vec3(20.16, 16.38, 0), Vec3(1.68, 14.7, 0), Vec3(3.36, 1.68, 0), Vec3(0.84, 0, 0)};
    BezierCurve3D b3(cp3, 3);
    b3.gen(output, cnt, cnt);
    for (int i = 0; i < cnt * 3; i++)
    {
        for (int j = 0; j <= cnt; j++)
        {
            int x, y, z, w;
            x = i * (cnt + 1) + j + 1;
            y = i * (cnt + 1) + (j + 1) % cnt + 1;
            w = (i + 1) * (cnt + 1) + j + 1;
            z = (i + 1) * (cnt + 1) + (j + 1) % (cnt) + 1;
            fprintf(output, "f %d %d %d %d\n", x, y, z, w);
        }
    }
}