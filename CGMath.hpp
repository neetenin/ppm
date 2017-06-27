#ifndef CGMATH
#define CGMATH

#include <cmath>
#include <cstdio>
// 三维向量，点类
class Vec3
{
public:
    union {
        struct { double x, y, z; };
        struct { double s, t, p; };
        struct { double r, g, b; };
    };
    // Vec3() : x(0.0), y(0.0), z(0.0) {} // 默认构造函数
    ~Vec3() {} // 析构
    Vec3(double x = 0.000, double y = 0.000, double z = 0.000) : x(x), y(y), z(z) {}
    Vec3(const Vec3 &u) : x(u.x), y(u.y), z(u.z) {}
    void print() const {printf("(%.2f, %.2f, %.2f)\n", x, y, z);}
    // 运算符重载
    Vec3& operator = (const Vec3 &u) { x = u.x; y = u.y; z = u.z; return *this; }
    Vec3 operator - () { return Vec3(-x, -y, -z); }
    double* operator & () { return (double*)this; }
    Vec3& operator += (double num) { x += num; y += num; z += num; return *this; }
    Vec3& operator += (const Vec3 &u) { x += u.x; y += u.y; z += u.z; return *this; }
    Vec3& operator -= (double num) { x -= num; y -= num; z -= num; return *this; }
    Vec3& operator -= (const Vec3 &u) { x -= u.x; y -= u.y; z -= u.z; return *this; }
    Vec3& operator *= (double num) { x *= num; y *= num; z *= num; return *this; }
    Vec3& operator *= (const Vec3 &u) { x *= u.x; y *= u.y; z *= u.z; return *this; }
    Vec3& operator /= (double num) { x /= num; y /= num; z /= num; return *this; }
    Vec3& operator /= (const Vec3 &u) { x /= u.x; y /= u.y; z /= u.z; return *this; }
    friend Vec3 operator + (const Vec3 &u, double num) { return Vec3(u.x + num, u.y + num, u.z + num); }
    friend Vec3 operator + (double num, const Vec3 &u) { return Vec3(num + u.x, num + u.y, num + u.z); }
    friend Vec3 operator + (const Vec3 &u, const Vec3 &v) { return Vec3(u.x + v.x, u.y + v.y, u.z + v.z); }
    friend Vec3 operator - (const Vec3 &u, double num) { return Vec3(u.x - num, u.y - num, u.z - num); }
    friend Vec3 operator - (double num, const Vec3 &u) { return Vec3(num - u.x, num - u.y, num - u.z); }
    friend Vec3 operator - (const Vec3 &u, const Vec3 &v) { return Vec3(u.x - v.x, u.y - v.y, u.z - v.z); }
    friend Vec3 operator * (const Vec3 &u, double num) { return Vec3(u.x * num, u.y * num, u.z * num); }
    friend Vec3 operator * (double num, const Vec3 &u) { return Vec3(num * u.x, num * u.y, num * u.z); }
    friend Vec3 operator * (const Vec3 &u, const Vec3 &v) { return Vec3(u.x * v.x, u.y * v.y, u.z * v.z); }
    friend Vec3 operator / (const Vec3 &u, double num) { return Vec3(u.x / num, u.y / num, u.z / num); }
    friend Vec3 operator / (double num, const Vec3 &u) { return Vec3(num / u.x, num / u.y, num / u.z); }
    friend Vec3 operator / (const Vec3 &u, const Vec3 &v) { return Vec3(u.x / v.x, u.y / v.y, u.z / v.z); }
    friend bool operator != (const Vec3&u, const Vec3 &v) {return fabs(u.x - v.x) > 1e-4 || fabs(u.y - v.y) > 1e-4 || fabs(u.z - v.z) > 1e-4;}
};


Vec3 cross(const Vec3 &u, const Vec3 &v)
{
    return Vec3(u.y * v.z - u.z * v.y, u.z * v.x - u.x * v.z, u.x * v.y - u.y * v.x);
}

double dot(const Vec3 &u, const Vec3 &v)
{
    return u.x * v.x + u.y * v.y + u.z * v.z;
}

double length(const Vec3 &u)
{
    return sqrt(u.x * u.x + u.y * u.y + u.z * u.z);
}

double length2(const Vec3 &u)
{
    return u.x * u.x + u.y * u.y + u.z * u.z;
}

Vec3 normalize(const Vec3 &u)
{
    return u * (1.0 / sqrt(u.x * u.x + u.y * u.y + u.z * u.z));
}

#endif