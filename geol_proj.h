#ifndef GEOL_PROJ_H
#define GEOL_PROJ_H

void WGS84toECEF_deg(
	double* x,   double* y,   double* z,
	double  lat, double  lon, double  h
);

void WGS84toECEF_rad(
	double* x,   double* y,   double* z,
	double  lat, double  lon, double  h
);

void ECEFtoWGS84_deg(
	double* lat, double* lon, double* h,
	double  x,   double  y,   double  z
);

void ECEFtoWGS84_rad(
	double* lat, double* lon, double* h,
	double  x,   double  y,   double  z
);

#ifdef GEOL_PROJ_IMPLEMENTATION

#include <math.h>

#define _geol_pi 3.1415926535898  // pi
#define _geol_a  6378137.0        // constant - earth major axis
#define _geol_b  6356752.3142     // constant - earth minor axis
#define _geol_f  0.00335281066475 // (a - b) / a - ellipsoid flatness
#define _geol_E  0.00669437999014 // f * (2 - f) - square of eccentricity
#define _geol_Es 0.00673949674228 // E / (1 - E)

#define _geol_degtorad(deg) ((deg) * _geol_pi / 180.0)
#define _geol_radtodeg(rad) ((rad) * 180.0 / _geol_pi)

void WGS84toECEF_deg(
	double* x,   double* y,   double* z,
	double  lat, double  lon, double  h
) {
	WGS84toECEF_rad(
		x, y, z,
		_geol_degtorad(lat),
		_geol_degtorad(lon),
		h
	);
}

void WGS84toECEF_rad(
	double* x,   double* y,   double* z,
	double  lat, double  lon, double  h
) {
	const double N = _geol_a / sqrt(1 - _geol_E * sin(lat) * sin(lat));

	*x = (h + N) * cos(lat) * cos(lon);
	*y = (h + N) * cos(lat) * sin(lon);
	*z = (h + N * (1 - _geol_E)) * sin(lat);
}

void ECEFtoWGS84_deg(
	double* lat, double* lon, double* h,
	double  x,   double  y,   double  z
) {
	ECEFtoWGS84_rad(lat, lon, h, x, y, z);
	*lat = _geol_radtodeg(*lat);
	*lon = _geol_radtodeg(*lon);
}

void ECEFtoWGS84_rad(
	double* lat, double* lon, double* h,
	double  x,   double  y,   double  z
) {
	const double p     = sqrt(x * x + y * y);
	const double q     = atan2((z * _geol_a), (p * _geol_b));
	const double qs    = sin(q);
	const double qc    = cos(q);
	const double qs3   = qs * qs * qs;
	const double qc3   = qc * qc * qc;
	const double phi_a = z + _geol_Es * _geol_b * qs3;
	const double phi_b = p - _geol_E * _geol_a * qc3;
	const double phi   = atan2(phi_a, phi_b);
	const double lam   = atan2(y, x);
	const double v     = _geol_a / sqrt(1 - _geol_E * sin(phi) * sin(phi));
	*lat = phi;
	*lon = lam;
	*h = (p / cos(phi)) - v;
}

#endif
#endif
