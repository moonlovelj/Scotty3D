#include "environment_light.h"

namespace CMU462 { namespace StaticScene {

EnvironmentLight::EnvironmentLight(const HDRImageBuffer* envMap)
    : envMap(envMap) {
  // TODO: initialize things here as needed
  double w = envMap->w;
  double h = envMap->h;
  vTheta.resize(h, 0);
  vThetaPhi.resize(h, std::vector<double>(w));
  vPhiGivenTheta.resize(h, std::vector<double>(w));

  double C = 0;
  for (int y = 0; y < h; y++) {
    double theta = (y + 0.5) / h * M_PI;
    double sinTheta = std::sin(theta);

    for (int x = 0; x < w; x++) {
      vThetaPhi[y][x] = envMap->data[x + y * w].illum() * sinTheta;
      C += vThetaPhi[y][x];
    }
  }

  for (int y = 0; y < h; y++) {
    for (int x = 0; x < w; x++) {
      vThetaPhi[y][x] /= C;
      vTheta[y] += vThetaPhi[y][x];
    }

    if (vTheta[y] != 0) {
      for (int x = 0; x < w; x++) {
        vPhiGivenTheta[y][x] = vThetaPhi[y][x] / vTheta[y];
      }
    }
  }

  for (int y = 0; y < h; y++) {
    if (y > 0) {
      vTheta[y] += vTheta[y - 1];
    }

    for (int x = 0; x < w; x++) {
      if (x > 0) {
        vPhiGivenTheta[y][x] += vPhiGivenTheta[y][x - 1];
      }
    }
  }
}

void EnvironmentLight::importance_sample(Vector3D *wi, float *pdf) const {
  double r1 = (double)(std::rand()) / RAND_MAX;
  double r2 = (double)(std::rand()) / RAND_MAX;
  double theta, phi;
  double x, y;
  int t, q;

  // Calculate theta.
  r1 *= vTheta.back();
  std::vector<double>::const_iterator itr = std::lower_bound(vTheta.begin(), vTheta.end(), r1);
  t = itr - vTheta.begin();
  double prev = t > 0 ? *(itr - 1) : 0;
  y = t + (r1 - prev) / (*itr - prev);
  theta = std::min(y / envMap->h, 1.0) * PI;

  // Calculate phi.
  r2 *= vPhiGivenTheta[t].back();
  itr = std::lower_bound(vPhiGivenTheta[t].begin(), vPhiGivenTheta[t].end(), r2);
  q = itr - vPhiGivenTheta[t].begin();
  prev = q > 0 ? *(itr - 1) : 0;
  x = q + (r2 - prev) / (*itr - prev);
  phi = std::min(x / envMap->w, 1.0) * 2 * PI;

  double sinTheta = std::sin(theta);
  double cosTheta = std::cos(theta);

  *pdf = vThetaPhi[t][q];
  *pdf /= (sinTheta * (2 * PI / envMap->w) *(PI / envMap->h));
  *wi = Vector3D(sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta);
}

Spectrum EnvironmentLight::sample_L(const Vector3D& p, Vector3D* wi,
                                    float* distToLight,
                                    float* pdf) const {
  // TODO: Implement
  // 	*pdf = 1.0 / (4 * M_PI);
  // 
  // 	double Xi1 = (double)(std::rand()) / RAND_MAX * 2.0 - 1.0;
  // 	double Xi2 = (double)(std::rand()) / RAND_MAX;
  //  double cosTheta = Xi1;
  //  double sinTheta = sqrt(1 - cosTheta*cosTheta);
  // 	double phi = 2.0 * PI * Xi2;
  // 	double xs = sinTheta * cosf(phi);
  // 	double ys = sinTheta * sinf(phi);
  // 	double zs = cosTheta;
  // 
  // 	*wi = Vector3D(xs, ys, zs);
  importance_sample(wi, pdf);
  *distToLight = INF_D;
	Ray r(p, *wi);
	return sample_dir(r);
}

Spectrum EnvironmentLight::sample_dir(const Ray& r) const {
  // TODO: Implement

	double theta = std::acos(r.d[2]);
	double phi = std::acos(cos_phi(r.d));
	if (r.d[1] < 0.0)
		phi = 2 * M_PI - phi;

	double u = phi / (2.0 * M_PI);
	double v = theta / M_PI;
	double tu = u * envMap->w - 0.5;
	double tv = v * envMap->h - 0.5;

	int su = (int)tu;
	int sv = (int)tv;

	double a, b;
	int px1, px2, py1, py2;

	if (tu < 0.0)
	{
		a = tu + 1;
		px1 = envMap->w - 1;
		px2 = 0;
	}
	else if (tu >= envMap->w - 1)
	{
		a = tu - envMap->w + 1;
		px1 = envMap->w - 1;
		px2 = 0;
	}
	else
	{
		a = tu - su;
		px1 = su;
		px2 = su + 1;
	}

	if (tv < 0.0)
	{
		b = tv + 1;
		py1 = envMap->h - 1;
		py2 = 0;
	}
	else if (tv >= envMap->h - 1)
	{
		b = tv - envMap->h + 1;
		py1 = envMap->h - 1;
		py2 = 0;
	}
	else
	{
		b = tv - sv;
		py1 = sv;
		py2 = sv + 1;
	}

	Spectrum z11 = envMap->data[px1 + envMap->w*py1];
	Spectrum z21 = envMap->data[px2 + envMap->w*py1];
	Spectrum z12 = envMap->data[px1 + envMap->w*py2];
	Spectrum z22 = envMap->data[px2 + envMap->w*py2];

	Spectrum zy1 = z11 * (1 - a) + z21 * a;
	Spectrum zy2 = z12 * (1 - a) + z22 * a;

	return zy1 * (1 - b) + zy2 * b;
}

} // namespace StaticScene
} // namespace CMU462
