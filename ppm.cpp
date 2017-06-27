#include "ppm.h"
// #include "omp.h"
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <assert.h>

#define PI ((double)3.14159265358979)
#define INF (1e20)

// ======================== QMC =========================== //

// Halton sequence with reverse permutation
int primes[61] = {
	2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79,
	83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181,
	191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283
};
inline int rev(const int i, const int p) {
	if (i == 0) return i; else return p - i;
}
double hal(const int b, int j) {
	const int p = primes[b];
	double h = 0.0, f = 1.0 / (double)p, fct = f;
	while (j > 0) {
		h += rev(j % p, p) * fct; j /= p; fct *= f;
	}
	return h;
}

// =========================HitPoint======================= //

void KdPoints::build(vector<HitPoint*> hps)
{
	root = build(hps, 0);
	min = Vec3(INF, INF, INF);
	max = Vec3(-INF, -INF, -INF);
	vector<HitPoint*>::iterator it;
	for (it = hps.begin(); it != hps.end(); it++)
	{
		(*it)->R2 = radius * radius;
		min.x = (*it)->pos.x < min.x ? (*it)->pos.x : min.x;
		min.y = (*it)->pos.y < min.y ? (*it)->pos.y : min.y;
		min.z = (*it)->pos.z < min.z ? (*it)->pos.z : min.z;
		max.x = (*it)->pos.x > max.x ? (*it)->pos.x : max.x;
		max.y = (*it)->pos.y > max.y ? (*it)->pos.y : max.y;
		max.z = (*it)->pos.z > max.z ? (*it)->pos.z : max.z;
	}
	min -= radius;
	max += radius;
}

HitPoint* KdPoints::build(vector<HitPoint*> hps, int dim)
{
	if (hps.size() < 1) return NULL;
	if (hps.size() < 2) return hps[0];

	if (dim == 0) sort(hps.begin(), hps.end(), smallerX);
	else if (dim == 1) sort(hps.begin(), hps.end(), smallerY);
	else sort(hps.begin(), hps.end(), smallerZ);

	int mid = hps.size() >> 1;
	vector<HitPoint*> lhp, rhp;
	for (int i = 0; i < mid; i++)
		lhp.push_back(hps[i]);
	for (int i = mid + 1; i < hps.size(); i++)
		rhp.push_back(hps[i]);
	hps[mid]->lc = build(lhp, (dim + 1) % 3);
	hps[mid]->rc = build(rhp, (dim + 1) % 3);
	hps[mid]->dim = dim;
	hps[mid]->R2 = radius * radius;
	return hps[mid];
}

void KdPoints::update(const Vec3 &photon, const Color &color, const Vec3 &n)
{
	if (photon.x < min.x && photon.y < min.y && photon.z < min.z) return;
	if (photon.x > max.x && photon.y > max.y && photon.z > max.z) return;
	update(root, photon, color, n);
}

int cnt = 0;

/*
	更新采样点信息
	r: 待更新的HitPoint
	photon: 光子位置
	n: 光子碰撞处法向
	color: 光子颜色
*/
void KdPoints::update(HitPoint* r, const Vec3 &photon, const Color &color, const Vec3 &n)
{
	if (!r) return;
	if (length2(r->pos - photon) <= r->R2 && dot(r->norm, n) > 1e-3)
	{
		double g = (r->n * alpha + alpha) / (r->n * alpha + 1.0);
		r->R2 = r->R2 * g;
		r->n++;
		r->flux = (r->flux + r->weight * color * (1. / PI)) * g;
		cnt++;
	}
	if (r->dim == 0) // x轴
	{
		if (photon.x >= r->pos.x - sqrt(r->R2)) update(r->rc, photon, color, n);
		if (photon.x <= r->pos.x + sqrt(r->R2)) update(r->lc, photon, color, n);
	}
	else if (r->dim == 1) // y轴
	{
		if (photon.y >= r->pos.y - sqrt(r->R2)) update(r->rc, photon, color, n);
		if (photon.y <= r->pos.y + sqrt(r->R2)) update(r->lc, photon, color, n);
	}
	else // z轴
	{
		if (photon.z >= r->pos.z - sqrt(r->R2)) update(r->rc, photon, color, n);
		if (photon.z <= r->pos.z + sqrt(r->R2)) update(r->lc, photon, color, n);
	}
}

// ========================= Light ======================= //

void Light::genPhoton(Ray& pr, Vec3& f, int i)
{
	f = color * (PI * 4.0); // flux
	double p = 2.*PI * hal(0, i), t = 2.*acos(sqrt(1. - hal(1, i)));
	double st = sin(t);
	pr.d = Vec3(cos(p) * st, cos(t), sin(p) * st);
	pr.o = pos;
}

void SpotLight::genPhoton(Ray& pr, Vec3& f, int i)
{
	f = color * (PI * 4.0); // flux
	double p = 2.*PI * hal(0, i), t = .3 * acos(sqrt(1. - hal(1, i)));
	double st = sin(t);
	pr.d = Vec3(cos(p) * st, cos(t), sin(p) * st);
	pr.o = pos;
}

// ====================== Camera ====================== //

Ray Camera::getEyeray(int x, int y)
{
	Vec3 dd = cx * ((x + 0.5) / w - 0.5) + cy * (-(y + 0.5) / h + 0.5) + d;
	Ray res = Ray(o + dd * 140, normalize(dd));
	return res;
}

// ===================== Object ======================= //

interInfo Sphere::intersect(const Ray &r)
{
	interInfo info;

	Vec3 op = center - r.o;
	if (dot(op, op) - dot(radius, radius) < 1e-3) info.inside = true;
	else info.inside = false;
	double t, b = dot(op, r.d), det = b * b - dot(op, op) + radius * radius;
	if (det < 0) {info.dis = INF;}
	else
	{
		det = sqrt(det);
		info.dis = (t = b - det) > 1e-8 ? t : ((t = b + det) > 1e-8 ? t : INF);
		info.x = r.o + r.d * t;
		info.n = normalize(info.x - center);
		info.inside = dot(info.n, r.d) > 0;
		info.n = info.inside ? info.n * -1 : info.n;
		info.d = normalize(r.d);
	}
	return info;
}

interInfo Plane::intersect(const Ray &r)
{
	interInfo info; info.dis = INF;
	Vec3 op(r.o - o), n;
	if (fabs(dot(norm, r.d)) <= 1e-8) return info;
	info.d = r.d;
	if (dot(norm, op) > 1e-8) {info.n = norm; info.inside = false;}
	else { info.n = norm * -1.; info.inside = true;}
	double dis = dot(info.n, op);
	double project = dot(info.n, r.d);
	if (project > 1e-8) return info;
	info.dis = dis / (-project);
	info.x = r.o + r.d * info.dis;
	return info;
}

void Plane::loadTexture(const char* t)
{
	has_texture = false;
	if (t == "") return;
	printf("Loading pic\n");
	texture = cv::imread(t);
	has_texture = true;
	// TODO;
	u = Vec3(0, 0, 1);
	v = Vec3(1, 0, 0);
}

Vec3 Plane::diff(const interInfo &info)
{
	Vec3 x = info.x;
	int h = texture.rows, w = texture.cols;
	if (h == 0 || w == 0) return m.diff;
	Vec3 dis = x - o;
	int du = (int(fabs(dis.x) / 110 * h)), dv = (int(fabs(dis.z) / 210 * w));
	if (du >= h || dv >= w) return m.diff;
	cv::Vec3b v3b = texture.at<cv::Vec3b>(du, dv);
	return Vec3((double)v3b[2] / 255, (double)v3b[1] / 255, (double)v3b[0] / 255);
}

Vec3 BezierObj::getPoint(double t)
{
	if (dim == 3)
		return  Vec3((1.0 - t) * (1.0 - t) * (1.0 - t) * controls[0]  +  (1.0 - t) * (1.0 - t) * t * 3 * controls[1] +
		             (1.0 - t) * t * t * 3 * controls[2] +  t * t * t * controls[3]);
	else
		return  Vec3((1.0 - t) * (1.0 - t) * controls[0] +  (1.0 - t) * t * 2 * controls[1] + t * t * controls[2]);
}

Vec3 BezierObj::getDif(double t)
{
	if (dim == 3)
		return Vec3(-3 * (1 - t) * (1 - t) * controls[0] + (9 * t * t - 12 * t + 3) * controls[1] + (-9 * t * t + 6 * t) * controls[2] + 3 * t * t * controls[3]);
	else
		return Vec3(-2 * (1 - t) * controls[0] + (- 4 * t + 2) * controls[1] + 2 * t * controls[2]);
}

Vec3 BezierObj::getPoint3D(double t, double theta)
{
	Vec3 temp = getPoint(t);
	temp.z = temp.x * sin(theta);
	temp.x = temp.x * cos(theta);
	return temp;
}

Vector3f BezierObj::targetFunc(const Vector3f &arg, const Ray& r)
{
	Vec3 p1 = getPoint3D(arg[1], arg[2]);
	Vec3 p2 = r.o + arg[0] * r.d;
	Vec3 temp = p2 - p1;
	return Vector3f(temp.x, temp.y, temp.z);
}

Vector3f BezierObj::NewtonIt(Vector3f arg, const Ray& r, int maxit, double th)
{
	Matrix3f f;
	Vec3 P = r.o + arg[0] * r.d;
	Vector3f res = arg;
	Vector3f func = targetFunc(res, r);
	for (int i = 0; i < maxit; i++)
	{
		if (arg[1] - 1 > 0.1 || arg[1] < -0.1) return Vector3f(INF, INF, INF);
		if (func[0] * func[0] + func[1] * func[1] + func[2] * func[2] < th) return res;
		Vec3 Point = getPoint(arg[1]);
		Vec3 Dif = getDif(arg[1]);
		f << r.d.x, -cos(arg[2]) * Dif.x, sin(arg[2]) * Point.x,
		r.d.y, -Dif.y, 0,
		r.d.z, -sin(arg[2]) * Dif.x, -cos(arg[2]) * Point.x;
		if (fabs(f.determinant()) < 1e-8) return Vector3f(INF, INF, INF);
		f = f.inverse().eval();
		res = res - f * func;
		func = targetFunc(res, r);
	}
	return Vector3f(INF, INF, INF);
}

int beziercnt = 0;
int totalcnt = 0;

Vector3f BezierObj::merge(double u1, double u2, const Ray& r)
{
	// 圆柱包围盒
	Vec3 p1 = getPoint(u1), p2 = getPoint(u2);
	double r1 = p1.x, r2 = p2.x;
	vector<double> cand;
	if (fabs(r.d.y) > 1e-8) // 圆面相交
	{
		double t1 = (p1.y - r.o.y) / r.d.y;
		Vec3 tempP; double PR2;
		if (t1 > 1e-8)
		{
			tempP = r.o + t1 * r.d;
			PR2 = sqrt(tempP.x * tempP.x + tempP.z * tempP.z);
			if ((PR2 - r1) * (PR2 - r2) < 0) cand.push_back(t1);
		}
		double t2 = (p2.y - r.o.y) / r.d.y;
		if (t2 > 1e-8)
		{
			tempP = r.o + t2 * r.d;
			PR2 = sqrt(tempP.x * tempP.x + tempP.z * tempP.z);
			if ((PR2 - r1) * (PR2 - r2) < 0) cand.push_back(t2);
		}
	}
	double a = r.o.x, b = r.d.x, c = r.o.z, d = r.d.z;
	if ((d * d + b * b) > 1e-8) // 侧面相交
	{
		double temp = -a * a * d * d + 2.*a * b * c * d + b * b * (-c * c) + b * b * r1 * r1 + d * d * r1 * r1;
		if (temp > 0)
		{
			double t3 = (-sqrt(temp) - a * b - c * d) / (b * b + d * d);
			if (t3 > 1e-8 && (p1.y - r.o.y - t3 * r.d.y) * (p2.y - r.o.y - t3 * r.d.y) < 0) cand.push_back(t3);
			double t4 = (sqrt(temp) - a * b - c * d) / (b * b + d * d);
			if (t4 > 1e-8 && (p1.y - r.o.y - t4 * r.d.y) * (p2.y - r.o.y - t4 * r.d.y) < 0) cand.push_back(t4);
		}
		temp = -a * a * d * d + 2.*a * b * c * d + b * b * (-c * c) + b * b * r2 * r2 + d * d * r2 * r2;
		if (temp > 0)
		{
			double t5 = (-sqrt(temp) - a * b - c * d) / (b * b + d * d);
			if (t5 > 1e-8 && (p1.y - r.o.y - t5 * r.d.y) * (p2.y - r.o.y - t5 * r.d.y) < 0) cand.push_back(t5);
			double t6 = (sqrt(temp) - a * b - c * d) / (b * b + d * d);
			if (t6 > 1e-8 && (p1.y - r.o.y - t6 * r.d.y) * (p2.y - r.o.y - t6 * r.d.y) < 0) cand.push_back(t6);
		}
	}
	int size = cand.size();
	if (size == 0) return Vector3f(INF, INF, INF);
	if (fabs(u1 - u2) < 1e-2)
	{
		Vector3f arg(INF, INF, INF);
		for (int i = 0; i < size; i++)
		{
			Vector3f temp_arg(cand[i], u1, atan2(r.o.z + cand[i] * r.d.z, r.o.x + cand[i] * r.d.x));
			Vector3f res = NewtonIt(temp_arg, r, 20, 1e-8);
			if (res[0] > 0 && res[0] < arg[0]) arg = res;
		}
		return arg;
	}
	Vector3f res1 = merge(u1, (u1 + u2) / 2, r), res2 = merge((u1 + u2) / 2, u2, r);
	return res1[0] < res2[0] ? res1 : res2;
}

interInfo BezierObj::intersect(const Ray& r)
{
	totalcnt++;
	Vector3f arg = merge(0, 1, r);
	interInfo info; info.dis = arg[0];
	if (info.dis < INF)
	{
		// printf("bingo\n");
		beziercnt++;
		info.d = r.d;
		info.x = r.o + info.dis * r.d;
		Vec3 xs = info.x - start;
		xs.y = 0;
		Vec3 xdif = getDif(arg[1]);
		xdif.z = xdif.x * sin(arg[2]);
		xdif.x = xdif.x * cos(arg[2]);
		info.n = cross(cross(axis, xdif), xdif);
		info.n = dot(info.n, r.d) > 0 ? info.n * -1 : info.n;
		info.n = normalize(info.n);
		info.u = arg[1], info.v = arg[2];
	}
	return info;
}

void BezierObj::loadTexture(const char* t)
{
	has_texture = false;
	if (t == "") return;
	printf("Loading pic\n");
	texture = cv::imread(t);
	has_texture = true;
}

Vec3 BezierObj::diff(const interInfo &info)
{
	int h = texture.rows, w = texture.cols;
	if (h == 0 || w == 0) return m.diff;
	int du = int(h * info.u), dv = ((int)((info.v + PI) / (2 * PI) * w) + 30) % w;
	cv::Vec3b v3b = texture.at<cv::Vec3b>(du, dv);
	return Vec3((double)v3b[2] / 255, (double)v3b[1] / 255, (double)v3b[0] / 255);
}

// ===================== Scene ====================== //

Scene::Scene()
{
	// 定义参数
	radius = 10.0;
	th = 1e-3;
	alpha = 0.8;
	w = 4096, h = 3072;
	epochs = 100000;

	// 定义光源

	// light.push_back(new Light(Vec3(-15, 60, 20), Color(10000, 10000, 10000)));
	// light.push_back(new Light(Vec3(65, 60, -20), Color(10000, 10000, 10000)));
	// light.push_back(new SpotLight(Vec3(65, 81.5, 10), Color(17500, 17500, 17500)));
	light.push_back(new SpotLight(Vec3(-29.5, 81.5, 0), Color(17500, 17500, 17500)));
	light.push_back(new SpotLight(Vec3(79.5, 81.5, 0), Color(17500, 17500, 17500)));
	light.push_back(new SpotLight(Vec3(25, 81.5, 20), Color(17500, 17500, 17500)));

	// 定义相机
	camera = new Camera(Vec3(25, 58, 210.6), Vec3(0, -0.122612, -1), .5135, .5135, w, h);

	// 定义材料

	Material WHITE_DIFF(Vec3(.75, .75, .75), Vec3(), Vec3(), 0);
	Material BLACK_DIFF(Vec3(), Vec3(), Vec3(), 0);
	Material RED_DIFF(Vec3(.75, .3, .3), Vec3(), Vec3(), 0);
	Material BLUE_DIFF(Vec3(.3, .3, .75), Vec3(), Vec3(), 0);
	Material LIGHT_BLUE_DIFF(Vec3(0.45, 0.45, .75), Vec3(), Vec3(), 0);
	// Material LIGHT_BLUE_DIFF2(Vec3(0.75, 0.6, .65), Vec3(), Vec3(), 0);
	Material LIGHT_RED_DIFF(Vec3(0.705, 0.45, .45), Vec3(), Vec3(), 0);
	Material MIRROR(Vec3(), Vec3(1, 1, 1)*.999, Vec3(), 0);
	Material MIRROR2(Vec3(.25, .25, .25), Vec3(1, 1, 1)*.999, Vec3(), 0);
	Material REFR0(Vec3(), Vec3(), Vec3(.999, .999, .999),  1.5);
	Material REFR1(Vec3(), Vec3(), Vec3(.75, .999, .999),  1.5);
	Material REFR2(Vec3(), Vec3(), Vec3(.999, .75, .75),  1.5);
	Material REFR3(Vec3(), Vec3(), Vec3(.75, .999, .75), 1.5);
	Material REFR4(Vec3(), Vec3(), Vec3(.999, .75, .999),  1.5);

	// Bezier曲线物体
	// Vec3 cp1[] = {Vec3(0.975, 41.34, 0), Vec3(0.39, 16.38, 0), Vec3(10.14, 17.0625, 0), Vec3(18.72, 15.99, 0)};
	// Vec3 cp2[] = {Vec3(18.72, 15.99, 0), Vec3(21.84, 15.6, 0), Vec3(18.72, 15.21, 0)};
	// Vec3 cp3[] = {Vec3(18.72, 15.21, 0), Vec3(1.56, 13.65, 0), Vec3(3.12, 1.56, 0), Vec3(0.78, 0, 0)};

	Vec3 cp1[] = {Vec3(1.045, 44.52, 0), Vec3(0.42, 17.64, 0), Vec3(10.92, 18.375, 0), Vec3(20.16, 17.22, 0)};
	Vec3 cp2[] = {Vec3(20.16, 17.22, 0), Vec3(23.52, 16.8, 0), Vec3(20.16, 16.38, 0)};
	Vec3 cp3[] = {Vec3(20.16, 16.38, 0), Vec3(1.68, 14.7, 0), Vec3(3.36, 1.68, 0), Vec3(0.84, 0, 0)};

	Object* obj[] = {
		//WHITE_DIFF, "texture.png"
		new Plane(Vec3(-30, 0, -85), Vec3(0, 1, 0), WHITE_DIFF, "texture.png"), // 地面
		// new Sphere(1e5, Vec3(-30 + 1e5, 40, 0), RED_DIFF),
		new Plane(Vec3(-30, 0, 125), Vec3(1, 0, 0), RED_DIFF), // 左面
		new Plane(Vec3(80, 0, 0), Vec3(-1, 0, 0), BLUE_DIFF), // 右面
		// new Sphere(1e5, Vec3(0, 40, -85 - 1e5), MIRROR),
		new Plane(Vec3(-30, 0, -85), Vec3(0, 0, 1), WHITE_DIFF), // 背面
		// new Sphere(1e5, Vec3(0, 81.6 - 1e5, 0), MIRROR),
		new Plane(Vec3(-30, 81.6, -85), Vec3(0, 1, 0), WHITE_DIFF), // 顶面
		new Plane(Vec3(0, 0, 125), Vec3(0, 0, 1), WHITE_DIFF), // 正面

		new Sphere(8, Vec3(32, 8, -20), REFR2),//
		new Sphere(1, Vec3(25, 80, 20), REFR0),// 顶灯
		new Sphere(1, Vec3(-29.5, 80, 0), REFR0),// 顶灯
		new Sphere(1, Vec3(79.5, 80, 0), REFR0),// 顶灯
		new Sphere(7.5, Vec3(24, 7.5, 40), REFR3), //
		// new Sphere(6.5, Vec3(0, 13, 50), MIRROR), // 玻璃球
		new Sphere(8, Vec3(48, 8, 0), MIRROR), //
		new Sphere(8, Vec3(60, 8, 27), MIRROR), //

		new BezierObj(cp1, 3, Vec3(0, 1, 0), Vec3(0, 8.2, 0), WHITE_DIFF, "bl2.jpg"),
		new BezierObj(cp2, 2, Vec3(0, 1, 0), Vec3(0, 7.8, 0), WHITE_DIFF, "bl2.jpg"),
		new BezierObj(cp3, 3, Vec3(0, 1, 0), Vec3(0, 0, 0), WHITE_DIFF, "bl2.jpg")
	};

	for (int i = 0; i < 16; i++)
		objects.push_back(obj[i]);
}

// ========================= RayTracing ======================== //

// 交点
bool Scene::intersect(const Ray& ray, interInfo& info)
{
	int n = objects.size();
	double inf = INF;
	info.dis = inf;
	for (int i = 0; i < n; i++)
	{
		interInfo info_t = objects[i]->intersect(ray);
		if (info_t.dis < info.dis)
		{
			info = info_t;
			info.id = i;
		}
	}
	return info.dis < inf;
}

void Scene::Raytracing(const Ray &ray, bool mode, int depth, Vec3 flux, Vec3 weight, int index)
{
	interInfo info;

	if (!intersect(ray, info)) return;
	if (weight.x < th || weight.y < th || weight.z < th || depth >= 19) return;

	depth++;
	int depth3 = depth * 3;
	Object* obj = objects[info.id]; // 碰撞物体

	Vec3 diff = obj->diff(info);
	Vec3 spec = obj->spec();
	Vec3 refr = obj->refr();
	double refr_rate = obj->refr_rate();

	if (diff != Vec3()) // 漫反射
	{
		Vec3 diff_weight = weight * diff;
		if (mode) // 采样点
		{
			HitPoint* hp = new HitPoint;
			hp->weight = diff_weight;
			hp->pos = info.x;
			hp->norm = info.n;
			hp->dir = info.d;
			hp->pixel_id = index;
			hps.push_back(hp);
		}
		else // 光子
		{
			Vec3 diff_flux = flux * diff;
			hpsMap->update(info.x, diff_flux, info.n);

			double r1 = 2.*PI * hal(depth3 - 1, index);
			double r2 = hal(depth3 + 0, index);
			double r2s = sqrt(r2);
			Vec3 w = info.n;
			Vec3 u = normalize(cross((fabs(w.x) > .1 ? Vec3(0, 1) : Vec3(1)), w));
			Vec3 v = cross(w, u);
			Vec3 new_dir = normalize(u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2));
			double p = diff.x > diff.y && diff.x > diff.z ? diff.x : diff.y > diff.z ? diff.y : diff.z; // max(x,y,z)
			if (hal(depth3 + 1, index) < p) Raytracing(Ray(info.x, new_dir), mode, depth, diff_flux * (1. / p), diff_weight, index);
		}
	}

	Vec3 reflect = normalize(info.d - info.n * 2.0 * dot(info.n, info.d));
	Ray refl_ray = Ray(info.x, reflect);
	if (spec != Vec3()) // 镜面反射
	{
		Vec3 spec_weight = weight * spec;
		Vec3 spec_flux = flux * spec;
		Raytracing(refl_ray, mode, depth, spec_flux, spec_weight, index);
	}
	if (refr != Vec3()) // 折射
	{
		// 折射方向
		double cos1 = dot(info.n, info.d);
		double nn = info.inside ? refr_rate : 1.0 / refr_rate; // 相对折射系数
		double cos2 = sqrt(1.0 -  (1.0 - cos1 * cos1) * nn * nn);

		if (cos2 < 0) return Raytracing(refl_ray, mode, depth, flux, weight, index); // 非全反射

		Vec3 T = normalize(info.d * nn - (cos2 + cos1 * nn) * info.n); // 折射方向
		Vec3 refr_weight = weight * refr;
		double a = refr_rate - 1.0, b =  refr_rate + 1.0;
		double R0 = a * a / (b * b), c = 1 - (info.inside ? fabs(dot(T, info.n)) : -cos1);
		double Re = R0 + (1 - R0) * c * c * c;
		Ray refr_ray = Ray(info.x, T);
		if (mode)
		{
			Raytracing(refl_ray, mode, depth, flux, refr_weight * Re, index);
			Raytracing(refr_ray, mode, depth, flux, refr_weight * (1.0 - Re), index);
		}
		else
		{
			(hal(depth3 - 1, index) < Re) ? Raytracing(refl_ray, mode, depth, flux, refr_weight, index) : Raytracing(refr_ray, mode, depth, flux, refr_weight, index);
		}
	}
}

void Scene::buildTree()
{
	hpsMap = new KdPoints(radius, alpha);
	hpsMap->build(hps);
}

void Scene::initHps()
{
	for (int y = 0; y < h; y++)
	{
		fprintf(stderr, "\rHitPointPass %5.2f%%", 100.0 * y / (h - 1));
		for (int x = 0; x < w; x++)
			Raytracing(camera->getEyeray(x, y), true, 0, Vec3(), Vec3(1, 1, 1), x + y * w);
	}
	fprintf(stderr, "\n");
	printf("\n%lu", hps.size());
	buildTree();
}

// 进行伽马校正
int toInt(double x) {
	return int(pow(1 - exp(-x), 1 / 2.2) * 255 + .5);
}

// 发射光子
void Scene::PhotonTracing()
{
	#pragma omp parallel for schedule(dynamic, 1)
	for (int i = 0; i < epochs; i++)	// 发射epochs轮光子
	{
		double p = 100. * (i + 1) / epochs;
		fprintf(stderr, "\rPhotonPass %5.3f%%", p);
		int m = 1000 * i;	// 总光子数
		Ray r;				// 光子方向
		Vec3 flux;			// 光子颜色
		for (int j = 0; j < 1000; j++)
		{
			for (int k = 0; k < light.size(); k++)
			{
				light[k]->genPhoton(r, flux, m + j);					// 每个光源产生光子
				Raytracing(r, false, 0, flux, Vec3(1, 1, 1), m + j);	// 对光子进行追踪
			}
		}
		if (i % 100 == 0)
		{
			estimation(i);
		}
	}
	estimation(epochs);
	printf("PPP: %d\n", cnt);
}

void Scene::estimation(int ix)
{
	Vec3 *c = new Vec3[w * h];
	vector<HitPoint*>::iterator it;
	for (it = hps.begin(); it != hps.end(); it++)
	{
		HitPoint* hp = *it;
		int i = hp->pixel_id;
		c[i] = c[i] + hp->flux * (1.0 / (PI * hp->R2 * ix * 1000.0));
	}

	// 保存图片
	char s[30]; sprintf(s, "result/image%d.ppm", ix);
	FILE* f = fopen(s, "w"); fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
	for (int i = 0; i < w * h; i++) {
		fprintf(f, "%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
	}
}

// ========================= Main ========================== //

int main()
{
	Scene s;
	s.initHps();
	cout << "Bezier" << beziercnt << endl;
	// // cout << "totalcnt" << totalcnt << endl;
	s.PhotonTracing();
	cout << "totalcnt" << totalcnt << endl;
}
