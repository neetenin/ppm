#ifndef PPM_H
#define PPM_H

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <vector>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include "Eigen/Core"
#include "Eigen/Dense"

#include "CGMath.hpp"

using namespace Eigen;
using namespace std;

#define Color Vec3
// =========================HitPoint======================= //

struct HitPoint
{
	HitPoint() {n = 0; lc = NULL; rc = NULL;}
	Vec3 pos, norm, dir, flux; // 位置，法向，射线方向，光通量
	Vec3 weight; // 对像素点的贡献权重
	int dim; // kdTree 维度
	int pixel_id; // 对应像素点
	double R2; // 半径平方
	unsigned int n; // 光子数量
	HitPoint *lc, *rc;
	void print() {pos.print(); printf("R2: %lf\n", R2);}
};

bool smallerX(HitPoint* a, HitPoint* b) {return a->pos.x < b->pos.x;}
bool smallerY(HitPoint* a, HitPoint* b) {return a->pos.y < b->pos.y;}
bool smallerZ(HitPoint* a, HitPoint* b) {return a->pos.z < b->pos.z;}

class KdPoints
{
public:
	KdPoints(double r, double a) : radius(r), alpha(a) {}
	void build(vector<HitPoint*> hps);
	void update(const Vec3 &photon, const Color &color, const Vec3 &d);
	void print()
	{
		root->print();
		root->lc->print();
		root->lc->lc->print();
		root->lc->rc->print();
		root->rc->print();
		root->rc->lc->print();
		root->rc->rc->print();
	}
private:
	HitPoint* root;
	HitPoint* build(vector<HitPoint*> hps, int dim);
	void update(HitPoint* r, const Vec3 &photon, const Color &color, const Vec3 &d);
	Vec3 min, max; // 包围盒
	double radius;
	double alpha;
};

// ===================== Ray ======================= //

struct Ray {Vec3 o, d; Ray() {}; Ray(Vec3 o_, Vec3 d_) : o(o_), d(d_) {}}; // 射线

// ===================== Light ======================= //

class Light
{
public:
	Light(Vec3 p, Color c): pos(p), color(c) {}
	virtual void genPhoton(Ray& pr, Vec3& f, int i);
private:
	Vec3 pos; // 点光源位置
	Color color; // 光源颜色
};

class SpotLight : public Light
{
public:
	SpotLight(Vec3 p, Color c): Light(p, c), pos(p), color(c) {}
	void genPhoton(Ray & pr, Vec3 & f, int i);
private:
	Vec3 pos; // 点光源位置
	Color color; // 光源颜色
};
// ====================== Camera ====================== //

class Camera
{
public:
	Camera(Vec3 o_, Vec3 d_, double lw, double lh, int w_, int h_):
		o(o_), d(normalize(d_)), lens_W(lw), lens_H(lh), w(w_), h(h_)
	{
		cx = Vec3(w * lens_W / h), cy = normalize(cross(cx , d)) * lens_H;
	}
	Ray getEyeray(int x, int y);
private:
	Vec3 o; // 视点
	Vec3 d; // 方向
	int w, h; // 镜头长宽
	double lens_W, lens_H;	//镜头的长（或宽）与感光点到镜头距离的比值
	Vec3 cx, cy;
};

// ===================== Object ======================= //

class interInfo // 交点信息
{
public:
	interInfo() {u = 0, v = 0, dis = 1e20;}
	int id;// 待判断相交物体
	double dis;
	Vec3 n; // 法向
	Vec3 x; // 相交点
	Vec3 d; // 入射方向
	bool inside; // 判断是否内部
	double u, v;
	void print()
	{
		printf("距离: %f\n交点: ", dis);
		x.print(); printf("法向: "); n.print();
		printf("入射方向: "); d.print();
		printf("内部: %d\n", inside);
	}
};

class Material // 材质
{
public:
	Material(Vec3 d, Vec3 s, Vec3 r, double re) : diff(d), spec(s), refr(r), refr_rate(re) {}
	Vec3 diff; // 漫反射系数
	Vec3 spec; // 镜面反射系数
	Vec3 refr; // 折射系数
	double refr_rate; // 折射率
};

class Object
{
public:
	Object(Material m_) : m(m_) {}
	virtual interInfo intersect(const Ray &r) = 0; // 求交
	virtual Vec3 diff(const interInfo &info) {return m.diff;}
	virtual Vec3 spec() {return m.spec;}
	virtual Vec3 refr() {return m.refr;}
	virtual double refr_rate() {return m.refr_rate;}
	Material m;
	Vec3 u, v;
};

class Sphere: public Object
{
public:
	Sphere(double r, Vec3 o, Material m_)
		: Object(m_), radius(r), center(o) {}
	interInfo intersect(const Ray &r); // 求交
private:
	double radius; // 半径
	Vec3 center; // 球心
};

class Plane : public Object
{
public:
	Plane(Vec3 o, Vec3 n, Material m, const char* texture = ""): Object(m), o(o), norm(n) {loadTexture(texture);}
	interInfo intersect(const Ray &r); // 求交
	Vec3 diff(const interInfo &info);
	// Vec3 u, v;
private:
	void loadTexture(const char* t);
	Vec3 o;
	Vec3 norm;
	bool has_texture;
	cv::Mat texture;
};

class BezierObj : public Object
{
public:
	BezierObj(Vec3* cps, int d, Vec3 a, Vec3 start, Material m, const char* texture = "")
		: Object(m), axis(a), start(start), dim(d)
	{
		controls = new Vec3[d + 1];
		for (int i = 0; i <= d; i++)
			controls[i] = cps[i];
		loadTexture(texture);
	}
	interInfo intersect(const Ray &r); // 求交
	Vec3 getPoint(double t);
	Vec3 getPoint3D(double t, double theta);
	Vec3 getDif(double t);
	Vec3 diff(const interInfo &info);
private:
	Vector3f targetFunc(const Vector3f &arg, const Ray& r);
	Vector3f NewtonIt(Vector3f arg, const Ray& r, int maxit, double th);
	Vector3f merge(double u1, double u2, const Ray& r);
	void loadTexture(const char* t);
	int dim;
	Vec3* controls;
	Vec3 axis;
	Vec3 start;
	bool has_texture;
	cv::Mat texture;
};

// ===================== Scene ====================== //

class Scene
{
public:
	Scene();
	void Raytracing(const Ray &ray, bool mode, int depth, Vec3 flux, Vec3 weight, int index);
	void estimation(int i);
	void initHps();
	void PhotonTracing();
private:
	int w, h; // 长宽
	int epochs; // 光子迭代数
	vector<Object*> objects;
	vector<Light*> light;
	Camera* camera;
	vector<HitPoint*> hps;
	KdPoints* hpsMap;
	double radius; // 初始邻域半径
	double th; // 最低衰减阈值
	double alpha; // 论文中的alpha
	void buildTree();
	bool intersect(const Ray& ray, interInfo& info);
};

#endif