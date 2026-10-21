#pragma once
#include <math.h>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <string>
#include <io.h>
#include <fcntl.h>

template <typename T, unsigned int D>
struct vec {
	T components[D];

	T& operator [](unsigned int i) {
		return components[i];
	}
	vec operator+ (vec v) {
		for (int i = 0; i < D; i++) {
			v[i] += components[i];
		}
		return v;
	}
	void operator+= (vec v) {
		for (int i = 0; i < D; i++) {
			components[i] += v[i];
		}
	}
	vec operator- (vec v) {
		for (int i = 0; i < D; i++) {
			v[i] = components[i] - v[i];
		}
		return v;
	}
	void operator-= (vec v) {
		for (int i = 0; i < D; i++) {
			components[i] -= v[i];
		}
	}
	float operator* (vec v) {
		float val = 0.0;
		for (int i = 0; i < D; i++) {
			val += (*this)[i] * v[i];
		}
		return val;
	}
	template <typename U>
	vec operator* (U s) {
		vec out;
		for (int i = 0; i < D; i++) {
			out[i] = components[i] * s;
		}
		return out;
	}
	template <typename U>
	vec operator*= (U s) {
		for (int i = 0; i < D; i++) {
			components[i] *= s;
		}
	}
};

using fvec2 = vec<float, 2>;
using fvec3 = vec<float, 3>;
using fvec4 = vec<float, 4>;
using dvec2 = vec<double, 2>;
using dvec3 = vec<double, 3>;
using dvec4 = vec<double, 4>;
using ivec2 = vec<int, 2>;
using ivec3 = vec<int, 3>;
using ivec4 = vec<int, 4>;
using uvec2 = vec<unsigned int, 2>;
using uvec3 = vec<unsigned int, 3>;
using uvec4 = vec<unsigned int, 4>;

template <typename T, unsigned int R, unsigned int C>
struct mat {
	vec<T, C> rows[R];

	vec<T, C>& operator [](unsigned int i) {
		return rows[i];
	}
	mat<T, R, C> operator +(mat<T, R, C> m) {
		for (int i = 0; i < R; i++) {
				m[i] += rows[i];
		}
		return m;
	}
	void operator +=(mat<T, R, C> m) {
		for (int i = 0; i < R; i++) {
				rows[i] += m[i];
		}
	}
	vec<T, R> operator *(vec<T, C> v) {
		vec<T, R> out;
		for (int i = 0; i < R; i++) {
			float sum = 0.0;
			for (int j = 0; j < C; j++) {
				sum += rows[i][j] * v[j];
			}
			out[i] = sum;
		}
		return out;
	}
	mat<T, R, C> operator *(T s) {
		vec<T, R> out;
		for (int i = 0; i < R; i++) {
			out[i] = rows[i] * s;
		}
		return out;
	}
	void operator *=(T s) {
		for (int i = 0; i < R; i++) {
			rows[i] *= s;
		}
	}
};

template <typename T, unsigned int D>
vec<T, D> init_vec() {
	vec<T, D> out;
	for (int i = 0; i < D; i++) {
		out[i] = T(0);
	}
	return out;
}
template <typename T, unsigned int D>
vec<T, D> init_vec(T s) {
	vec<T, D> out;
	for (int i = 0; i < D; i++) {
		out[i] = s;
	}
	return out;
}
template <typename T, unsigned int D>
vec<T, D> init_vec(const T (&rows)[D]) {
	vec<T, D> out;
	for (int i = 0; i < D; i++) {
		out[i] = rows[i];
	}
	return out;
}
template <typename T, unsigned int D>
float vecs_dist_sq(vec<T, D> a, vec<T, D> b) {
	float dist = 0.0;
	for (int i = 0; i < D; i++) {
		float delta = (float)a[i] - (float)b[i];
		dist += delta * delta;
	}
	return dist;
}

template <typename T, typename U, unsigned int D>
vec<T, D> scale_vec(vec<T, D> v, U scale) {
	for (int i = 0; i < D; i++) {
		v[i] *= scale;
	}
	return v;
}
template <typename T, unsigned int D>
float dot_vec(vec<T, D> a, vec<T, D> b) {
	float product = 0.0f;
	for (int i = 0; i < D; i++) {
		product += a[i] * b[i];
	}
	return product;
}
template <typename T, unsigned int D>
vec<T, D> mod_vec(vec<T, D> v, float period) {
	for (int i = 0; i < D; i++) {
		v[i] = fmodf((float)v[i], period);
	}
	return v;
}

template <typename T, unsigned int D>
vec<T, D> scale_vec(vec<T, D> v, vec<T, D> scale) {
	for (int i = 0; i < D; i++) {
		v[i] *= scale[i];
	}
	return v;
}

template <typename T, unsigned int D>
vec<T, D> sub_vecs(vec<T, D> a, vec<T, D> b) {

	vec<T, D> out;
	for (int i = 0; i < D; i++) {
		out[i] = a[i] - b[i];
	}
	return out;
}

template <typename T>
vec<T, 3> cross_vecs(vec<T, 3> a, vec<T, 3> b) {
	fvec3 v;
	v[0] = a[1] * b[2] - a[2] * b[1];
	v[1] = a[2] * b[0] - a[0] * b[2];
	v[2] = a[0] * b[1] - a[1] * b[0];
	return v;
}
template <typename T, unsigned int D>
vec<T, D> scale_vec(vec<T, D> v, T s) {
	vec<T, D> out;
	for (int i = 0; i < D; i++) {
		out[i] = v[i] * s;
	}
	return out;
}

template <typename T, unsigned int D>
float vec_length(vec<T, D> v) {
	return sqrt(dot_vec(v, v));
}
template <typename T, unsigned int D>
vec<T, D> normalize(vec<T, D> v) {
	float vec_length = sqrt(v * v);
	return scale_vec(v, 1.0 / vec_length);
}

template <typename T>
mat<T, 3, 3> rotation_mat(vec<T, 3> pivot, T angle) {
	float s = sin(angle);
	float c = cos(angle);
	return mat<T, 3, 3>({ {
		{(pivot[0] * pivot[0]) * (1 - c) + c,							(pivot[0] * pivot[1]) * (1 - c) - pivot[2] * s,	(pivot[0] * pivot[2]) * (1 - c) + pivot[1] * s},
		{(pivot[0] * pivot[1]) * (1 - c) + pivot[2] * s,	(pivot[1] * pivot[1]) * (1 - c) + c,							(pivot[1] * pivot[2]) * (1 - c) - pivot[0] * s},
		{(pivot[0] * pivot[2]) * (1 - c) - pivot[1] * s,	(pivot[1] * pivot[2]) * (1 - c) + pivot[0] * s,	(pivot[2] * pivot[2]) * (1 - c) + c}
	} });
}
template <typename T, unsigned int R, unsigned int C>
mat<T, R, C> zero_mat() {
	vec<T, R> out;
	for (int i = 0; i < R; i++) {
		for (int j = 0; j < C; j++) {
			out[i][j] = T(0);
		}
	}
	return out;
}
template <typename T, unsigned int R, unsigned int C>
vec<T, R> mat_mult_vec(mat<T, R, C> m, vec<T, C> v) {
	vec<T, R> out;
	for (int i = 0; i < R; i++) {
		float sum = 0.0;
		for (int j = 0; j < C; j++) {
			sum += m[i][j] * v[j];
		}
		out[i] = sum;
	}
	return out;
}

template <typename T, unsigned int D>
float unit_vec_arc_length(vec<T, D> a, vec<T, D> b) {
	//return the arc length between two unit vectors
	return acos(a * b);
}
template <typename T, unsigned int D>
float vec_arc_length(vec<T, D> a, vec<T, D> b) {
	//return the arc length between two unit vectors
	return acos(a * b / (vec_length(a) * vec_length(b)));
}

template <typename T, unsigned int D>
vec<T, D> solve_system_of_equations(mat<T, D, D> M, vec<T, D> augment) {
	vec<T, D> solution;

	//forward elimination
	for (int i = 0; i < D; i++) {
		for (int j = i + 1; j < D; j++) {
			float scaling_factor = -1.0 * M[j][i] / M[i][i];
			for (int k = i + 1; k < D; k++) {
				M[j][k] += scaling_factor * M[i][k];
			}
			augment[j] += scaling_factor * augment[i];
		}
	}

	//back substitution:
	solution[D-1] = augment[D - 1] / M[D-1][D-1];
	for (int i = D - 2; i >= 0; i--) {
		float sum = augment[i];
		for (int j = i + 1; j < D; j++) {
			sum -= M[i][j] * solution[j];
		}
		solution[i] = sum / M[i][i];
	}

	return solution;
}


template <typename T, unsigned int R, unsigned int C>
mat<T, C, R> transpose(mat<T, R, C> m) {
	mat<T, C, R> out;
	for (int i = 0; i < C; i++) {
		for (int j = 0; j < R; j++) {
			out[j][i] = m.rows[i][j];
		}
	}
	return out;
}

template <typename T, unsigned int D>
struct span {
	vec<T, D> min;
	vec<T, D> max;

	span<T, D> operator+(vec<T, D> v) {
		span<T, D> s;
		s.min = min + v;
		s.max = max + v;
		return s;
	}
	void operator+=(vec<T, D> v) {
		min += v;
		max += v;
	}
	span<T, D> operator-(vec<T, D> v) {
		span<T, D> s;
		s.min = min - v;
		s.max = max - v;
		return s;
	}
	void operator-=(vec<T, D> v) {
		min -= v;
		max -= v;
	}
	template <typename U>
	span<T, D> operator*(U x) {
		span<T, D> s;
		s.min = min * x;
		s.max = max * x;
		if (x < 0) {
			return {s.max, s.min};
		}
		return s;
	}
};

template <typename T, unsigned int D>
span<T, D> scale_span(span<T, D> s, vec<T, D> scale) {
	s.min = scale_vec(s.min, scale);
	s.max = scale_vec(s.max, scale);
	return s;
}

template <typename T, unsigned int D>
span<T, D> span_overlap(span<T, D> a, span<T, D> b) {
	span<T, D> overlap;
	for (int i = 0; i < D; i++) {
		if (a.max[i] < b.min[i] || b.max[i] < a.min[i]) {
			return { {{0.0, 0.0}},{{0.0, 0.0}} };
		}
		overlap.min[i] = std::clamp(b.min[i], a.min[i], a.max[i]);
		overlap.max[i] = std::clamp(b.max[i], a.min[i], a.max[i]);
	}
	return overlap;
}

template <typename T, unsigned int D>
bool vec_in_span(vec<T, D> v, span<T, D> s) {
	for (int i = 0; i < D; i++) {
		if (v[i] < s.min[i] || v[i] > s.max[i]) {
			return false;
		}
	}
	return true;
}

template <typename T, unsigned int D>
span<T, D> expand_span_to_grid(span<T, D> s, vec<T, D> offset, vec<T, D> increment) {
	span<T, D> scaled = scale_span(s, increment);
	span<T, D> out;
	for (int i = 0; i < D; i++) {
		out.min[i] = offset[i] + increment[i] * std::floor(scaled.min[i]);
		out.max[i] = offset[i] + increment[i] * std::ceil(scaled.max[i]);
	}
	return out;
}

using frect = span<float, 2>;
using irect = span<int, 2>;
using urect = span<unsigned int, 2>;

template <typename T, unsigned int D>
void print_vector(vec<T, D> v) {
	std::cout << "[";
	for (int c = 0; c < D; c++) {
		std::cout << v[c];
		if (c < D - 1) {
			std::cout << ", ";
		}
	}
	std::cout << "]";
}

template <typename T, unsigned int R, unsigned int C>
void print_matrix(mat<T, R, C> M) {
	std::wstring output = L"";
	for (int i = 0; i < R; i++) {
		//determine opening and closing characters:
		std::wstring open = L"\u2502";
		std::wstring close = L"\u2502";
		if (i == 0) {
			open = L'\u250C';
			close = L'\u2510';
		}
		else if (i == R - 1) {
			open = L'\u2514';
			close = L'\u2518';  
		}
		output += open;

		for (int j = 0; j < C; j++) {
			output += std::to_wstring(M[i][j]);
			if (j != C - 1) {
				output += L"\t";
			}
		}
		output += close + L'\n';
	}
	_setmode(_fileno(stdout), _O_U16TEXT);
	std::wcout << output;
	//set back to normal so that printing works correctly in other function
	_setmode(_fileno(stdout), _O_TEXT);
}