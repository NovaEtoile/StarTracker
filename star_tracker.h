// Steps: find the stars (adaptive threshold, brightness-aware filter, optional hot-pixel
// map and second-picture check) -> optionally estimate the real field of view from the
// stars -> identify each star (log-likelihood on neighbour angles,
// -> reject stars that disagree with each other -> solve the rotation (Wahba's problem)
// 
// Camera frame of the answer: x = right, y = DOWN, z = forward (out of the lens) 
// Sky directions are unit vectors in the catalogue's frame (x, y, z from RA and Dec)
//   camera direction = attitude * sky direction
// The direction the camera points at, in the sky, is the third row of the attitude
//
// Everything here was tested on simulated pictures only
#pragma once

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <functional>
#include "linear_algebra.h"

// reading the star catalogue and the rotation maths

struct hip_star {
    double x, y, z;   //direction (unit vector)
    float mag;        //brightness (smaller number = brighter)
};

// skips ahead to the next number in the text and reads it
static bool read_number(const char*& p, double& value) {
    while (*p && !((*p >= '0' && *p <= '9') || *p == '.')) {
        p++;
    }
    if (!*p) {
        return false;
    }
    char* end = NULL;
    value = std::strtod(p, &end);
    p = end;
    return true;
}

// text file. One line looks like (Is this even important):
//   3: R.A: 00h\t00m\t01s DEC: +38°\t51m\t33s brightness: 6.61
static bool parse_logan_line(const std::string& line, hip_star* out) {
    const char* text = line.c_str();
    const char* ra = std::strstr(text, "R.A:");
    const char* dec = std::strstr(text, "DEC:");
    const char* br = std::strstr(text, "brightness:");

    if (!ra || !dec || !br) {
        return false;
    }

    double h = 0, m = 0, s = 0;
    const char* q = ra + 4;
    if (!read_number(q, h) || !read_number(q, m) || !read_number(q, s)) {
        return false;
    }

    q = dec + 4;
    while (*q && *q != '+' && *q != '-') {
        q++;
    }
    if (!*q) {
        return false;
    }
    double sign = (*q == '-') ? -1.0 : 1.0;
    q++;

    double dd = 0, dm = 0, ds = 0;
    if (!read_number(q, dd) || !read_number(q, dm) || !read_number(q, ds)) {
        return false;
    }

    q = br + 11;
    while (*q == ' ') {
        q++;
    }
    double mag_sign = 1.0;
    if (*q == '-') {
        mag_sign = -1.0;
        q++;
    }
    double mag = 0;
    if (!read_number(q, mag)) {
        return false;
    }
    mag *= mag_sign;

    const double two_pi = 6.283185307179586;

    double ra_angle = two_pi * (h / 24.0 + m / 1440.0 + s / 86400.0);
    double dec_angle = sign * two_pi * (dd / 360.0 + dm / 21600.0 + ds / 1296000.0);

    out->x = cos(ra_angle) * cos(dec_angle);
    out->y = sin(ra_angle) * cos(dec_angle);
    out->z = sin(dec_angle);
    out->mag = (float)mag;

    return true;
}

// reads one number out of a piece of text. Returns false if it is blank
static bool field_to_number(const std::string& text, double* value) {
    const char* start = text.c_str();
    while (*start == ' ') {
        start++;
    }
    if (*start == 0) {
        return false;
    }
    char* end = NULL;
    double v = std::strtod(start, &end);
    if (end == start) {
        return false;
    }
    *value = v;
    return true;
}

// the official ESA/CDS file (hip_main.dat). Fields are separated by '|'
// field 5 = Vmag (brightness), field 8 = RA in degrees, field 9 = Dec in degrees
static bool parse_official_line(const std::string& line, hip_star* out) {
    std::vector<std::string> f;

    size_t begin = 0;
    while (true) {
        size_t bar = line.find('|', begin);
        if (bar == std::string::npos) {
            f.push_back(line.substr(begin));
            break;
        }
        f.push_back(line.substr(begin, bar - begin));
        begin = bar + 1;
    }

    if (f.size() < 10) {
        return false;
    }

    double vmag = 0, ra_deg = 0, dec_deg = 0;

    if (!field_to_number(f[5], &vmag) ||
        !field_to_number(f[8], &ra_deg) ||
        !field_to_number(f[9], &dec_deg)) {
        return false;
    }

    if (ra_deg < 0.0 || ra_deg > 360.0 || dec_deg < -90.0 || dec_deg > 90.0) {
        return false;
    }

    const double radians_per_degree = 0.017453292519943295;
    double ra_angle = ra_deg * radians_per_degree;
    double dec_angle = dec_deg * radians_per_degree;

    out->x = cos(ra_angle) * cos(dec_angle);
    out->y = sin(ra_angle) * cos(dec_angle);
    out->z = sin(dec_angle);
    out->mag = (float)vmag;

    return true;
}

static bool parse_hipparcos_line(const std::string& line, hip_star* out) {
    if (line.find("R.A:") != std::string::npos) {
        return parse_logan_line(line, out);
    }
    return parse_official_line(line, out);
}

// returns -1 if the file cannot be opened, otherwise the number of stars read
// 'first_line' gets the first line of the file (to help if nothing could be read)
static int load_hipparcos_file(
    const char* path,
    std::vector<hip_star>& stars,
    std::string* first_line
) {
    std::ifstream file(path, std::ios::binary);
    if (!file) {
        return -1;
    }

    std::string line;
    bool first = true;

    while (std::getline(file, line)) {
        if (first) {
            *first_line = line.substr(0, 150);
            first = false;
        }

        hip_star s;
        if (parse_hipparcos_line(line, &s)) {
            stars.push_back(s);
        }
    }

    return (int)stars.size();
}

// Working out which way the camera points (Wahba's problem)

struct vec3d {
    double v[3];
};

struct mat3d {
    double m[3][3];
};

static vec3d make_vec3d(double x, double y, double z) {
    vec3d r;
    r.v[0] = x;
    r.v[1] = y;
    r.v[2] = z;
    return r;
}

static double dot3(const vec3d& a, const vec3d& b) {
    return a.v[0] * b.v[0] + a.v[1] * b.v[1] + a.v[2] * b.v[2];
}

static vec3d cross3(const vec3d& a, const vec3d& b) {
    return make_vec3d(
        a.v[1] * b.v[2] - a.v[2] * b.v[1],
        a.v[2] * b.v[0] - a.v[0] * b.v[2],
        a.v[0] * b.v[1] - a.v[1] * b.v[0]
    );
}

static vec3d normalize3(const vec3d& a) {
    double length = sqrt(dot3(a, a));
    return make_vec3d(a.v[0] / length, a.v[1] / length, a.v[2] / length);
}

static vec3d mat_vec(const mat3d& R, const vec3d& x) {
    return make_vec3d(
        R.m[0][0] * x.v[0] + R.m[0][1] * x.v[1] + R.m[0][2] * x.v[2],
        R.m[1][0] * x.v[0] + R.m[1][1] * x.v[1] + R.m[1][2] * x.v[2],
        R.m[2][0] * x.v[0] + R.m[2][1] * x.v[1] + R.m[2][2] * x.v[2]
    );
}

// angle between two vectors (accurate even when the angle is tiny)
static double angle_between3(const vec3d& a, const vec3d& b) {
    vec3d c = cross3(a, b);
    double sine = sqrt(dot3(c, c));
    double cosine = dot3(a, b);
    return atan2(sine, cosine);
}

// finds eigenvalues and eigenvectors of a symmetric 4x4 matrix
// eigenvectors come back as the COLUMNS of 'vectors'
static void jacobi_eigen_4x4(double A[4][4], double values[4], double vectors[4][4]) {
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            vectors[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }

    for (int sweep = 0; sweep < 60; sweep++) {
        double off = 0.0;
        for (int p = 0; p < 3; p++) {
            for (int q = p + 1; q < 4; q++) {
                off += A[p][q] * A[p][q];
            }
        }
        if (off < 1.0e-30) {
            break;
        }

        for (int p = 0; p < 3; p++) {
            for (int q = p + 1; q < 4; q++) {
                if (fabs(A[p][q]) < 1.0e-300) {
                    continue;
                }

                double theta = (A[q][q] - A[p][p]) / (2.0 * A[p][q]);
                double t = (theta >= 0.0 ? 1.0 : -1.0) /
                    (fabs(theta) + sqrt(theta * theta + 1.0));
                double c = 1.0 / sqrt(t * t + 1.0);
                double s = t * c;

                for (int k = 0; k < 4; k++) {
                    double akp = A[k][p];
                    double akq = A[k][q];
                    A[k][p] = c * akp - s * akq;
                    A[k][q] = s * akp + c * akq;
                }
                for (int k = 0; k < 4; k++) {
                    double apk = A[p][k];
                    double aqk = A[q][k];
                    A[p][k] = c * apk - s * aqk;
                    A[q][k] = s * apk + c * aqk;
                }
                for (int k = 0; k < 4; k++) {
                    double vkp = vectors[k][p];
                    double vkq = vectors[k][q];
                    vectors[k][p] = c * vkp - s * vkq;
                    vectors[k][q] = s * vkp + c * vkq;
                }
            }
        }
    }

    for (int i = 0; i < 4; i++) {
        values[i] = A[i][i];
    }
}

// best rotation R such that  body[i] = R * ref[i]  for the chosen stars
// 'body' = direction seen by the camera, 'ref' = direction in the sky catalogue
// This is Wahba's problem 
// solved with the quaternion method, which gives the same answer as the SVD method and also
// works when only 2 or 3 stars are available
static bool wahba_solve(
    const std::vector<vec3d>& body,
    const std::vector<vec3d>& ref,
    const std::vector<int>& active,
    mat3d* R
) {
    if (active.size() < 2) {
        return false;
    }

    double B[3][3] = { {0,0,0},{0,0,0},{0,0,0} };

    for (size_t k = 0; k < active.size(); k++) {
        const vec3d& b = body[active[k]];
        const vec3d& r = ref[active[k]];
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                B[i][j] += b.v[i] * r.v[j];
            }
        }
    }

    double sigma = B[0][0] + B[1][1] + B[2][2];
    double Z[3] = {
        B[1][2] - B[2][1],
        B[2][0] - B[0][2],
        B[0][1] - B[1][0]
    };

    double K[4][4];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            K[i][j] = B[i][j] + B[j][i] - ((i == j) ? sigma : 0.0);
        }
        K[i][3] = Z[i];
        K[3][i] = Z[i];
    }
    K[3][3] = sigma;

    double values[4];
    double vectors[4][4];
    jacobi_eigen_4x4(K, values, vectors);

    int best = 0;
    for (int i = 1; i < 4; i++) {
        if (values[i] > values[best]) {
            best = i;
        }
    }

    double q[4];
    double length = 0.0;
    for (int i = 0; i < 4; i++) {
        q[i] = vectors[i][best];
        length += q[i] * q[i];
    }
    length = sqrt(length);
    if (length == 0.0) {
        return false;
    }
    for (int i = 0; i < 4; i++) {
        q[i] /= length;
    }

    // q[0..2] = vector part, q[3] = scalar part
    double qv2 = q[0] * q[0] + q[1] * q[1] + q[2] * q[2];
    double skew[3][3] = {
        { 0.0,   -q[2],  q[1] },
        { q[2],   0.0,  -q[0] },
        { -q[1],  q[0],  0.0  }
    };

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            double value = 2.0 * q[i] * q[j] - 2.0 * q[3] * skew[i][j];
            if (i == j) {
                value += q[3] * q[3] - qv2;
            }
            R->m[i][j] = value;
        }
    }

    return true;
}

// solves for the camera rotation in two stages, so one or two wrongly
// identified stars cannot spoil the answer:
//   1. Agreement check: the distance between two stars in the picture must
//      match the distance between their catalogue stars. A wrongly identified
//      star disagrees with almost everyone, so it is thrown out
//   2. Fit, then throw out the star that fits worst and try again, until
//      everything left agrees
// Returns false (no answer) rather than risk giving a wrong answer
static bool estimate_attitude(
    const std::vector<vec3d>& body,
    const std::vector<vec3d>& ref,
    double tolerance_radians,
    int min_stars,
    mat3d* R,
    int* stars_used
) {
    int n = (int)body.size();

    // Stage 1: agreement check
    double pair_tolerance = 2.0 * tolerance_radians;

    std::vector<int> agree(n, 0);

    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            double in_picture = angle_between3(body[i], body[j]);
            double in_catalogue = angle_between3(ref[i], ref[j]);

            if (fabs(in_picture - in_catalogue) < pair_tolerance) {
                agree[i]++;
                agree[j]++;
            }
        }
    }

    int max_agree = 0;
    for (int i = 0; i < n; i++) {
        if (agree[i] > max_agree) {
            max_agree = agree[i];
        }
    }

    if (max_agree < min_stars - 1) {
        return false;
    }

    std::vector<int> active;
    for (int i = 0; i < n; i++) {
        if (agree[i] * 2 >= max_agree && agree[i] >= 1) {
            active.push_back(i);
        }
    }

    // Stage 2: fit and trim
    while ((int)active.size() >= min_stars) {
        mat3d fit;
        if (!wahba_solve(body, ref, active, &fit)) {
            return false;
        }

        double worst = -1.0;
        int worst_position = -1;

        for (int k = 0; k < (int)active.size(); k++) {
            vec3d predicted = mat_vec(fit, ref[active[k]]);
            double e = angle_between3(predicted, body[active[k]]);
            if (e > worst) {
                worst = e;
                worst_position = k;
            }
        }

        if (worst <= tolerance_radians) {
            *R = fit;
            *stars_used = (int)active.size();
            return true;
        }

        active.erase(active.begin() + worst_position);
    }

    return false;
}

// angle (radians) of the rotation that takes one orientation to the other
static double rotation_error(const mat3d& a, const mat3d& b) {
    double trace = 0.0;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            trace += a.m[i][j] * b.m[i][j];
        }
    }
    double c = (trace - 1.0) / 2.0;
    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;
    return acos(c);
}

// angle (radians) between where the two orientations point the camera
static double boresight_error(const mat3d& a, const mat3d& b) {
    vec3d za = make_vec3d(a.m[2][0], a.m[2][1], a.m[2][2]);
    vec3d zb = make_vec3d(b.m[2][0], b.m[2][1], b.m[2][2]);
    return angle_between3(za, zb);
}

// Matching a star to the catalogue

struct star_candidate {
    unsigned int id;
    float error;   // negative log-likelihood: LOWER is a better match
};

static float clamp_dot(fvec3 a, fvec3 b) {
    float value = a * b;

    if (value > 1.0f) value = 1.0f;
    if (value < -1.0f) value = -1.0f;

    return value;
}

static float angular_distance(fvec3 a, fvec3 b) {
    return acos(clamp_dot(a, b));
}

static fvec3 cross_product(fvec3 a, fvec3 b) {
    return fvec3({
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]
    });
}

static fvec3 camera_ray_from_pixel(
    float x,
    float y,
    float width,
    float height,
    float fov,
    float aspect_ratio
) {
    float normalized_x = (2.0f * x / width - 1.0f) * aspect_ratio;
    float normalized_y = 1.0f - 2.0f * y / height;
    float scale = tan(fov / 2.0f);

    fvec3 ray = fvec3({
        normalized_x * scale,
        normalized_y * scale,
        1.0f
    });

    return normalize(ray);
}

// returns the 'count' smallest angles from star 'index' to every other star,
// in order from smallest to largest
static std::vector<float> nearest_angles(
    const std::vector<fvec3>& stars,
    int index,
    int count
) {
    std::vector<float> dots;
    dots.reserve(stars.size());

    fvec3 centre = stars[index];

    for (int i = 0; i < (int)stars.size(); i++) {
        if (i == index) {
            continue;
        }

        fvec3 other = stars[i];
        dots.push_back(centre * other);
    }

    size_t keep = std::min((size_t)count, dots.size());

    std::partial_sort(
        dots.begin(),
        dots.begin() + keep,
        dots.end(),
        std::greater<float>()
    );

    std::vector<float> angles;
    angles.reserve(keep);

    for (size_t k = 0; k < keep; k++) {
        float d = dots[k];
        if (d > 1.0f) d = 1.0f;
        if (d < -1.0f) d = -1.0f;
        angles.push_back(acos(d));
    }

    return angles;
}

static float log_likelihood(
    const std::vector<float>& observed,
    const std::vector<float>& catalogue,
    float sigma
) {
    const float window = 7.0f * sigma;
    const float two_sigma_sq = 2.0f * sigma * sigma;

    float total = 0.0f;

    for (size_t j = 0; j < observed.size(); j++) {
        float phi = observed[j];

        std::vector<float>::const_iterator it =
            std::lower_bound(catalogue.begin(), catalogue.end(), phi - window);

        float sum = 0.0f;

        for (; it != catalogue.end() && *it <= phi + window; ++it) {
            float d = phi - *it;
            sum += exp(-(d * d) / two_sigma_sq);
        }

        if (sum < 1.0e-12f) {
            sum = 1.0e-12f;
        }

        total += log(sum);
    }

    return total;
}

static std::vector<star_candidate> find_candidates(
    const std::vector<fvec3>& observed_rays,
    int observed_index,
    const std::vector<std::vector<float>>& catalogue_features,
    int neighbour_count,
    float sigma,
    int candidate_count
) {
    std::vector<float> observed_features =
        nearest_angles(observed_rays, observed_index, neighbour_count);

    std::vector<star_candidate> candidates;
    candidates.reserve(catalogue_features.size());

    for (int i = 0; i < (int)catalogue_features.size(); i++) {
        star_candidate candidate;
        candidate.id = (unsigned int)i;
        candidate.error = -log_likelihood(observed_features, catalogue_features[i], sigma);
        candidates.push_back(candidate);
    }

    size_t keep = std::min((size_t)candidate_count, candidates.size());

    std::partial_sort(
        candidates.begin(),
        candidates.begin() + keep,
        candidates.end(),
        [](const star_candidate& a, const star_candidate& b) {
            return a.error < b.error;
        }
    );

    candidates.resize(keep);

    return candidates;
}

// Finding stars in a picture

inline int pixel_index(int x, int y, int width) {
    return y * width + x;
}

// a bright patch found in the image, with its size and shape
struct blob {
    fvec2 centroid;     // brightness-weighted centre (same method as before)
    int area;           // number of pixels in the patch
    int peak;           // brightest pixel
    double elongation;  // 1 = round, bigger = longer and thinner
};

// same flood fill as detect_stars, but it also measures each patch's size and shape
static std::vector<blob> detect_blobs(
    const unsigned char* image_data,
    int width,
    int height,
    int threshold
) {
    std::vector<unsigned char> visited(width * height, 0);
    std::vector<blob> blobs;

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {

            int start = pixel_index(x, y, width);

            if (visited[start] || image_data[start] < threshold) {
                continue;
            }

            std::vector<int> queue;
            queue.push_back(start);
            visited[start] = 1;

            double weighted_x = 0.0, weighted_y = 0.0, total_weight = 0.0;
            double sum_x = 0.0, sum_y = 0.0, sum_xx = 0.0, sum_yy = 0.0, sum_xy = 0.0;
            int peak = 0;

            for (size_t q = 0; q < queue.size(); q++) {

                int current = queue[q];
                int current_x = current % width;
                int current_y = current / width;
                int value = image_data[current];

                weighted_x += current_x * (double)value;
                weighted_y += current_y * (double)value;
                total_weight += value;

                sum_x += current_x;
                sum_y += current_y;
                sum_xx += (double)current_x * current_x;
                sum_yy += (double)current_y * current_y;
                sum_xy += (double)current_x * current_y;
                if (value > peak) peak = value;

                for (int dy = -1; dy <= 1; dy++) {
                    for (int dx = -1; dx <= 1; dx++) {

                        if (dx == 0 && dy == 0) continue;

                        int nx = current_x + dx;
                        int ny = current_y + dy;

                        if (nx < 0 || nx >= width || ny < 0 || ny >= height) continue;

                        int neighbour = pixel_index(nx, ny, width);

                        if (visited[neighbour] || image_data[neighbour] < threshold) continue;

                        visited[neighbour] = 1;
                        queue.push_back(neighbour);
                    }
                }
            }

            if (total_weight > 0.0) {
                double n = (double)queue.size();
                double mean_x = sum_x / n;
                double mean_y = sum_y / n;

                // spread of the patch, with one pixel's own size (1/12) added so a
                // single pixel counts as perfectly round
                double cxx = sum_xx / n - mean_x * mean_x + 1.0 / 12.0;
                double cyy = sum_yy / n - mean_y * mean_y + 1.0 / 12.0;
                double cxy = sum_xy / n - mean_x * mean_y;

                double half_sum = (cxx + cyy) / 2.0;
                double half_diff = sqrt(((cxx - cyy) / 2.0) * ((cxx - cyy) / 2.0) + cxy * cxy);
                double long_axis = half_sum + half_diff;
                double short_axis = half_sum - half_diff;
                if (short_axis < 1.0e-9) short_axis = 1.0e-9;

                blob b;
                b.centroid = fvec2({ (float)(weighted_x / total_weight), (float)(weighted_y / total_weight) });
                b.area = (int)queue.size();
                b.peak = peak;
                b.elongation = sqrt(long_axis / short_axis);
                blobs.push_back(b);
            }
        }
    }

    return blobs;
}

static std::vector<fvec2> blob_centroids(const std::vector<blob>& blobs) {
    std::vector<fvec2> out;
    for (size_t i = 0; i < blobs.size(); i++) {
        out.push_back(blobs[i].centroid);
    }
    return out;
}

// estimates the background level and noise from the whole image 
// stars cover well under 1% of the pixels, so the middle of the distribution is background
static void estimate_background(const unsigned char* image_data, int width, int height,
    double* background, double* noise_sigma) {

    long total = (long)width * height;

    std::vector<long> histogram(256, 0);
    for (long i = 0; i < total; i++) {
        histogram[image_data[i]]++;
    }

    long running = 0;
    int median = 0;
    for (int v = 0; v < 256; v++) {
        running += histogram[v];
        if (running * 2 >= total) {
            median = v;
            break;
        }
    }

    std::vector<long> deviation(256, 0);
    for (int v = 0; v < 256; v++) {
        deviation[abs(v - median)] += histogram[v];
    }

    running = 0;
    int mad = 0;
    for (int d = 0; d < 256; d++) {
        running += deviation[d];
        if (running * 2 >= total) {
            mad = d;
            break;
        }
    }

    *background = median;
    *noise_sigma = 1.4826 * mad;
}

// brightness-aware filter. A real star of a given brightness has a predictable size:
// the brighter it is, the more pixels it covers. A blob far SMALLER than that (a bright
// single pixel, a small bright square) is not a star. A faint blob may be just one pixel
static std::vector<fvec2> brightness_filter(
    const std::vector<blob>& blobs,
    double background,
    double threshold,
    double star_sigma_pixels,
    double fraction,
    double max_elongation
) {
    const double pi = 3.14159265358979;
    double threshold_above = threshold - background;
    if (threshold_above < 1.0) threshold_above = 1.0;

    std::vector<fvec2> out;

    for (size_t i = 0; i < blobs.size(); i++) {

        double peak_above = (double)blobs[i].peak - background;
        double expected_area = 1.0;
        if (peak_above > threshold_above) {
            expected_area = pi * 2.0 * star_sigma_pixels * star_sigma_pixels * log(peak_above / threshold_above);
        }

        int needed = (int)std::floor(fraction * expected_area);
        if (needed < 1) needed = 1;

        if (blobs[i].area >= needed && blobs[i].elongation <= max_elongation) {
            out.push_back(blobs[i].centroid);
        }
    }

    return out;
}

// hot-pixel map: take a picture with no stars (a "dark frame") and note every pixel
// bright enough to be mistaken for a star. These get replaced by their neighbours' value
static std::vector<char> make_hot_pixel_map(const unsigned char* dark_data, int width, int height, int threshold) {
    long total = (long)width * height;
    std::vector<char> mask(total, 0);
    for (long i = 0; i < total; i++) {
        if (dark_data[i] >= threshold) {
            mask[i] = 1;
        }
    }
    return mask;
}

static void fill_masked_pixels(unsigned char* image_data, int width, int height,
    const std::vector<char>& mask, const std::vector<long>& masked_list, double background) {

    for (size_t m = 0; m < masked_list.size(); m++) {
        long i = masked_list[m];
        int x = (int)(i % width);
        int y = (int)(i / width);

        std::vector<int> values;
        for (int dy = -1; dy <= 1; dy++) {
            for (int dx = -1; dx <= 1; dx++) {
                if (dx == 0 && dy == 0) continue;
                int nx = x + dx;
                int ny = y + dy;
                if (nx < 0 || nx >= width || ny < 0 || ny >= height) continue;
                long j = (long)ny * width + nx;
                if (mask[j]) continue;
                values.push_back(image_data[j]);
            }
        }

        if (values.empty()) {
            image_data[i] = (unsigned char)background;
        }
        else {
            std::sort(values.begin(), values.end());
            image_data[i] = (unsigned char)values[values.size() / 2];
        }
    }
}

// two-picture check: keep a detection only if the second picture has one at the same place
static std::vector<fvec2> confirm_in_both(std::vector<fvec2> first, std::vector<fvec2> second, double tolerance_pixels) {
    std::vector<fvec2> out;
    for (size_t i = 0; i < first.size(); i++) {
        fvec2 a = first[i];
        for (size_t j = 0; j < second.size(); j++) {
            fvec2 b = second[j];
            double dx = a[0] - b[0];
            double dy = a[1] - b[1];
            if (dx * dx + dy * dy <= tolerance_pixels * tolerance_pixels) {
                out.push_back(a);
                break;
            }
        }
    }
    return out;
}

// Estimating the field of view from the stars

// tries many field-of-view guesses and keeps the one that makes the most stars agree
// For each guess, the picture's star directions are recomputed, a handful of stars are
// matched to their best catalogue star, and we count how many of those matches are
// consistent with each other (the distance between two stars in the picture must equal
// the distance between their catalogue stars) At the right field of view many matches
// agree. At a wrong one they do not, because the guesses are scattered at random.
// returns the estimated error as a fraction (0.01 = real field of view is 1% bigger)
static double estimate_fov_error(
    const std::vector<fvec2>& detected,
    const std::vector<std::vector<float>>& features,
    const std::vector<fvec3>& catalogue_dirs,
    float fov,
    int width,
    int height,
    float aspect_ratio,
    int neighbours
) {
    int n = (int)detected.size();
    if (n < 5) {
        return 0.0;
    }

    // use up to 20 stars, spread through the list
    std::vector<int> subset;
    int use = n < 20 ? n : 20;
    for (int k = 0; k < use; k++) {
        subset.push_back((int)((long)k * n / use));
    }

    const float radians_per_pixel = fov / (float)height;

    struct guess_result {
        int agreeing;        // size of the biggest group of mutually consistent matches
        double mismatch;     // total matching error (used to break ties)
        double refinement;   // small extra correction (fraction) suggested by the agreeing pairs
    };

    auto try_guess = [&](double error, float sigma_pixels, double abs_tolerance_pixels, double relative_tolerance) -> guess_result {

        float fov_guess = (float)(fov * (1.0 + error));

        std::vector<fvec3> rays;
        for (int i = 0; i < n; i++) {
            fvec2 p = detected[i];
            rays.push_back(camera_ray_from_pixel(p[0], p[1], (float)width, (float)height, fov_guess, aspect_ratio));
        }

        int m = (int)subset.size();
        std::vector<int> ids(m, -1);
        std::vector<fvec3> observed(m);
        double mismatch = 0.0;

        for (int k = 0; k < m; k++) {
            std::vector<star_candidate> c =
                find_candidates(rays, subset[k], features, neighbours, sigma_pixels * radians_per_pixel, 1);
            if (!c.empty()) {
                ids[k] = (int)c[0].id;
                mismatch += c[0].error;
            }
            observed[k] = rays[subset[k]];
        }

        std::vector<int> agree(m, 0);
        std::vector<double> pair_observed;
        std::vector<double> pair_catalogue;
        std::vector<int> pair_first;
        std::vector<int> pair_second;

        for (int i = 0; i < m; i++) {
            for (int j = i + 1; j < m; j++) {
                if (ids[i] < 0 || ids[j] < 0 || ids[i] == ids[j]) continue;

                fvec3 ci = catalogue_dirs[ids[i]];
                fvec3 cj = catalogue_dirs[ids[j]];
                double a_observed = angular_distance(observed[i], observed[j]);
                double a_catalogue = angular_distance(ci, cj);

                double tolerance = abs_tolerance_pixels * radians_per_pixel + relative_tolerance * a_observed;

                if (fabs(a_observed - a_catalogue) < tolerance) {
                    agree[i]++;
                    agree[j]++;
                    pair_observed.push_back(a_observed);
                    pair_catalogue.push_back(a_catalogue);
                    pair_first.push_back(i);
                    pair_second.push_back(j);
                }
            }
        }

        int best_agree = 0;
        for (int i = 0; i < m; i++) {
            if (agree[i] > best_agree) best_agree = agree[i];
        }

        // extra correction from the agreeing pairs: how much to stretch the picture's angles
        double refinement = 0.0;
        double sum_cross = 0.0;
        double sum_square = 0.0;
        int pairs_used = 0;
        for (size_t q = 0; q < pair_observed.size(); q++) {
            if (agree[pair_first[q]] * 2 >= best_agree && agree[pair_second[q]] * 2 >= best_agree) {
                sum_cross += pair_observed[q] * pair_catalogue[q];
                sum_square += pair_observed[q] * pair_observed[q];
                pairs_used++;
            }
        }
        if (pairs_used >= 3 && sum_square > 0.0) {
            refinement = sum_cross / sum_square - 1.0;
        }

        guess_result r;
        r.agreeing = best_agree;
        r.mismatch = mismatch;
        r.refinement = refinement;
        return r;
    };

    auto better = [](const guess_result& a, const guess_result& b) {
        if (a.agreeing != b.agreeing) return a.agreeing > b.agreeing;
        return a.mismatch < b.mismatch;
    };

    // coarse search: -6% to +6% in steps of 0.25%, forgiving matching
    double best_error = 0.0;
    guess_result best = try_guess(0.0, 2.0f, 2.0, 0.0015);    // start from "no error" so ties keep the assumed value

    for (int step = -24; step <= 24; step++) {
        double e = step * 0.0025;
        guess_result r = try_guess(e, 2.0f, 2.0, 0.0015);
        if (better(r, best)) {
            best = r;
            best_error = e;
        }
    }

    if (best.agreeing < 2) {
        return 0.0;      // no group of stars agreed anywhere: keep the assumed field of view
    }

    // fine search: around the best coarse guess in steps of 0.1%, stricter matching
    double centre = best_error;
    guess_result fine_best = try_guess(centre, 1.0f, 1.0, 0.001);
    best_error = centre;

    for (int step = -6; step <= 6; step++) {
        double e = centre + step * 0.001;
        guess_result r = try_guess(e, 1.0f, 1.0, 0.001);
        if (better(r, fine_best)) {
            fine_best = r;
            best_error = e;
        }
    }

    // apply the small extra correction measured from the agreeing pairs
    return (1.0 + best_error) * (1.0 + fine_best.refinement) - 1.0;
}

// The star tracker: settings, catalogue, and one call that does everything

struct st_settings {
    // the picture and the assumed camera
    int width = 1920;
    int height = 1080;
    float fov = 0.13f;                     // vertical field of view in radians (what you believe it is)

    // finding stars
    double threshold_sigmas = 6.0;         // a pixel must be this many noise-widths above the background
    int min_threshold = 20;                // out of 255
    int max_threshold = 200;
    double star_sigma_pixels = 1.5;        // width of a star's blob (set this from your optics)
    double size_fraction = 0.5;            // a blob must cover this fraction of the size a star of that brightness should
    double max_elongation = 2.0;           // 1 = round
    int max_detections = 300;              // more than this means the picture is too noisy to use
    double two_picture_tolerance_pixels = 1.0;

    // identifying stars
    int neighbours = 6;
    float sigma_pixels = 0.3f;             // how far a measured angle may be off, in pixels

    // working out the direction
    double trim_tolerance_pixels = 2.0;    // a star farther than this from its fitted place is thrown out
    int min_stars = 3;

    // calibration
    bool self_calibrate_fov = true;        // search for the real field of view (-6% to +6%)
};

struct st_catalogue {
    std::vector<fvec3> directions;
    std::vector<float> magnitudes;
    std::vector<std::vector<float>> features;      // each star's nearest neighbour angles, smallest first
};

struct st_hot_pixel_map {
    std::vector<char> flag;
    std::vector<long> list;
};

struct st_result {
    bool solved = false;
    mat3d attitude = { { {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0} } };
    int stars_found = 0;           // stars left after the filters
    int stars_used = 0;            // stars that agreed and went into the answer
    double fov_error = 0.0;        // estimated error of the assumed field of view (0.01 = 1% bigger)
    int threshold = 0;             // brightness threshold used (out of 255)
    std::string message;
};

// keeps stars at least as bright as 'magnitude_cutoff'. Pick a cutoff that matches what the camera
// can see: not shallower than the camera, and at most about one magnitude deeper
inline st_catalogue st_build_catalogue(const std::vector<hip_star>& stars, float magnitude_cutoff, int catalogue_neighbours = 20) {
    st_catalogue c;
    for (size_t i = 0; i < stars.size(); i++) {
        if (stars[i].mag <= magnitude_cutoff) {
            c.directions.push_back(fvec3({ (float)stars[i].x, (float)stars[i].y, (float)stars[i].z }));
            c.magnitudes.push_back(stars[i].mag);
        }
    }
    for (int i = 0; i < (int)c.directions.size(); i++) {
        c.features.push_back(nearest_angles(c.directions, i, catalogue_neighbours));
    }
    return c;
}

inline int st_threshold(double background, double noise_sigma, const st_settings& s) {
    double t = background + s.threshold_sigmas * noise_sigma;
    if (t < s.min_threshold) t = s.min_threshold;
    if (t > s.max_threshold) t = s.max_threshold;
    return (int)std::ceil(t);
}

// map of hot pixels, from a picture taken with no stars (a dark frame). Stuck-bright pixels are
// repaired from their neighbours before stars are looked for. Remake it now and then: sensors
// in orbit develop new hot pixels
inline st_hot_pixel_map st_make_hot_pixel_map(const unsigned char* dark_picture, const st_settings& s) {
    double background = 0.0, noise_sigma = 0.0;
    estimate_background(dark_picture, s.width, s.height, &background, &noise_sigma);

    st_hot_pixel_map m;
    m.flag = make_hot_pixel_map(dark_picture, s.width, s.height, st_threshold(background, noise_sigma, s));
    for (long i = 0; i < (long)m.flag.size(); i++) {
        if (m.flag[i]) m.list.push_back(i);
    }
    return m;
}

// star centres (pixel coordinates, x right, y down) found in one picture
inline std::vector<fvec2> st_find_stars(const unsigned char* picture, const st_settings& s,
    const st_hot_pixel_map* hot_map, int* threshold_used) {

    long total = (long)s.width * s.height;

    double background = 0.0, noise_sigma = 0.0;
    estimate_background(picture, s.width, s.height, &background, &noise_sigma);
    int threshold = st_threshold(background, noise_sigma, s);
    if (threshold_used) *threshold_used = threshold;

    std::vector<unsigned char> repaired;
    const unsigned char* data = picture;
    if (hot_map && !hot_map->list.empty()) {
        repaired.assign(picture, picture + total);
        fill_masked_pixels(repaired.data(), s.width, s.height, hot_map->flag, hot_map->list, background);
        data = repaired.data();
    }

    std::vector<blob> blobs = detect_blobs(data, s.width, s.height, threshold);
    return brightness_filter(blobs, background, threshold, s.star_sigma_pixels, s.size_fraction, s.max_elongation);
}

// Which way is the camera pointing?
//   picture         the picture (8-bit, width*height bytes, row by row from the top)
//   second_picture  optional second picture of the same sky; only stars seen in both are used
//   hot_map         optional hot-pixel map. Pass nullptr to skip
inline st_result st_solve(
    const st_settings& s,
    const st_catalogue& catalogue,
    const unsigned char* picture,
    const unsigned char* second_picture = nullptr,
    const st_hot_pixel_map* hot_map = nullptr
) {
    st_result r;
    float aspect_ratio = (float)s.width / (float)s.height;

    std::vector<fvec2> stars = st_find_stars(picture, s, hot_map, &r.threshold);

    if (second_picture) {
        std::vector<fvec2> second = st_find_stars(second_picture, s, hot_map, nullptr);
        stars = confirm_in_both(stars, second, s.two_picture_tolerance_pixels);
    }

    r.stars_found = (int)stars.size();

    if (r.stars_found > s.max_detections) {
        r.message = "too many detections: the picture is probably too noisy";
        return r;
    }
    if (r.stars_found < s.min_stars) {
        r.message = "too few stars found";
        return r;
    }

    float fov = s.fov;
    if (s.self_calibrate_fov) {
        r.fov_error = estimate_fov_error(stars, catalogue.features, catalogue.directions,
            s.fov, s.width, s.height, aspect_ratio, s.neighbours);
        fov = (float)(s.fov * (1.0 + r.fov_error));
    }
    float radians_per_pixel = fov / (float)s.height;

    std::vector<fvec3> rays;
    for (size_t i = 0; i < stars.size(); i++) {
        fvec2 p = stars[i];
        rays.push_back(camera_ray_from_pixel(p[0], p[1], (float)s.width, (float)s.height, fov, aspect_ratio));
    }

    // identify every star (best guess), then let the agreement check sort out the wrong ones
    std::vector<vec3d> body;
    std::vector<vec3d> ref;

    for (int i = 0; i < (int)rays.size(); i++) {
        std::vector<star_candidate> candidates =
            find_candidates(rays, i, catalogue.features, s.neighbours, s.sigma_pixels * radians_per_pixel, 1);

        if (candidates.empty()) continue;

        int guess = (int)candidates[0].id;
        fvec3 ray = rays[i];
        fvec3 sky = catalogue.directions[guess];

        body.push_back(make_vec3d(ray[0], -ray[1], ray[2]));     // camera frame: x right, y down, z forward
        ref.push_back(make_vec3d(sky[0], sky[1], sky[2]));
    }

    if ((int)body.size() < s.min_stars) {
        r.message = "too few stars could be matched";
        return r;
    }

    if (!estimate_attitude(body, ref, s.trim_tolerance_pixels * radians_per_pixel, s.min_stars, &r.attitude, &r.stars_used)) {
        r.message = "the identified stars did not agree with each other";
        return r;
    }

    r.solved = true;
    r.message = "ok";
    return r;
}
