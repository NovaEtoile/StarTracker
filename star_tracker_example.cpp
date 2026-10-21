// Reads a star catalogue and a picture, and prints which way the camera is pointing.
//
// Usage:
//   star_tracker_example <catalogue file> <picture.pgm> [options]
//
// The picture must be an 8-bit greyscale PGM file (binary "P5" format). Most image programs
// can save one, for example:  magick photo.png -colorspace Gray -depth 8 photo.pgm
//
// Options (write them as name=value):
//   fov=0.13        vertical field of view of the camera, in radians
//   cutoff=7        faintest catalogue star to use (magnitude); match this to what the camera can see
//   star_sigma=1.5  width of a star's blob in pixels
//   selfcal=1       1 = search for the real field of view (-6% to +6%), 0 = trust fov exactly
//   second=<file>   a second picture of the same sky (keeps only stars seen in both; removes radiation hits)
//   dark=<file>     a picture with no stars (makes a map of hot pixels, which are then ignored)
//
// Catalogue: either star_data_8.0_hipparcos.txt or the official hip_main.dat

#include <iostream>
#include <iomanip>
#include <sstream>
#include <cctype>
#include "star_tracker.h"

static bool read_pgm(const std::string& path, std::vector<unsigned char>* data, int* width, int* height) {
    std::ifstream f(path.c_str(), std::ios::binary);
    if (!f) return false;

    std::string magic;
    f >> magic;
    if (magic != "P5") return false;

    auto next_number = [&]() -> int {
        char c;
        while (f.get(c)) {
            if (c == '#') { while (f.get(c) && c != '\n') {} }
            else if (!std::isspace((unsigned char)c)) { f.unget(); break; }
        }
        int v = 0;
        f >> v;
        return v;
    };

    *width = next_number();
    *height = next_number();
    int maximum = next_number();
    if (*width <= 0 || *height <= 0 || maximum != 255) return false;

    f.get();    // the single space after the header

    data->resize((size_t)(*width) * (*height));
    f.read((char*)data->data(), (std::streamsize)data->size());
    return (size_t)f.gcount() == data->size();
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cout << "Usage: star_tracker_example <catalogue file> <picture.pgm> [fov=0.13] [cutoff=7] [star_sigma=1.5]\n"
                  << "                           [selfcal=1] [second=<file.pgm>] [dark=<file.pgm>]\n";
        return 1;
    }

    st_settings settings;
    float cutoff = 7.0f;
    std::string second_path, dark_path;

    for (int i = 3; i < argc; i++) {
        std::string arg = argv[i];
        size_t eq = arg.find('=');
        if (eq == std::string::npos) { std::cout << "Could not understand: " << arg << "\n"; return 1; }
        std::string name = arg.substr(0, eq);
        std::string value = arg.substr(eq + 1);

        if (name == "fov") settings.fov = (float)std::atof(value.c_str());
        else if (name == "cutoff") cutoff = (float)std::atof(value.c_str());
        else if (name == "star_sigma") settings.star_sigma_pixels = std::atof(value.c_str());
        else if (name == "selfcal") settings.self_calibrate_fov = std::atoi(value.c_str()) != 0;
        else if (name == "second") second_path = value;
        else if (name == "dark") dark_path = value;
        else { std::cout << "Unknown option: " << name << "\n"; return 1; }
    }

    std::vector<hip_star> stars;
    std::string first_line;
    if (load_hipparcos_file(argv[1], stars, &first_line) <= 0) {
        std::cout << "Could not read any stars from " << argv[1] << "\n";
        if (!first_line.empty()) std::cout << "First line of that file: " << first_line << "\n";
        return 1;
    }

    std::vector<unsigned char> picture, second, dark;
    int w = 0, h = 0, w2 = 0, h2 = 0;

    if (!read_pgm(argv[2], &picture, &w, &h)) {
        std::cout << "Could not read " << argv[2] << " as an 8-bit binary PGM picture.\n";
        return 1;
    }
    settings.width = w;
    settings.height = h;

    if (!second_path.empty() && (!read_pgm(second_path, &second, &w2, &h2) || w2 != w || h2 != h)) {
        std::cout << "Could not read the second picture, or it is a different size.\n";
        return 1;
    }

    st_hot_pixel_map hot_map;
    bool have_hot_map = false;
    if (!dark_path.empty()) {
        if (!read_pgm(dark_path, &dark, &w2, &h2) || w2 != w || h2 != h) {
            std::cout << "Could not read the dark picture, or it is a different size.\n";
            return 1;
        }
        hot_map = st_make_hot_pixel_map(dark.data(), settings);
        have_hot_map = true;
        std::cout << "Hot pixels mapped: " << hot_map.list.size() << "\n";
    }

    std::cout << "Catalogue: " << stars.size() << " stars read; using those up to magnitude " << cutoff << "\n";
    st_catalogue catalogue = st_build_catalogue(stars, cutoff);
    std::cout << "           " << catalogue.directions.size() << " stars in use\n";

    st_result r = st_solve(settings, catalogue, picture.data(),
        second_path.empty() ? nullptr : second.data(), have_hot_map ? &hot_map : nullptr);

    std::cout << "Stars found: " << r.stars_found << "   brightness threshold: " << r.threshold << " / 255\n";

    if (!r.solved) {
        std::cout << "NO ANSWER: " << r.message << "\n";
        return 2;
    }

    std::cout << "Stars used: " << r.stars_used << "\n";
    if (settings.self_calibrate_fov) {
        std::cout << "Field of view correction found: " << std::showpos << std::fixed << std::setprecision(2)
                  << 100.0 * r.fov_error << "%" << std::noshowpos << "  (field of view used: "
                  << std::setprecision(5) << settings.fov * (1.0 + r.fov_error) << " rad)\n";
    }

    // the direction the camera points at, in the sky = third row of the attitude
    const double two_pi = 6.283185307179586;
    const double pi = two_pi / 2.0;
    vec3d bore = make_vec3d(r.attitude.m[2][0], r.attitude.m[2][1], r.attitude.m[2][2]);
    double ra = atan2(bore.v[1], bore.v[0]);
    if (ra < 0.0) ra += two_pi;
    double dec = asin(bore.v[2]);

    double hours = ra / two_pi * 24.0;
    int rh = (int)hours;
    int rm = (int)((hours - rh) * 60.0);
    double rs = ((hours - rh) * 60.0 - rm) * 60.0;

    double dec_deg = dec * 180.0 / pi;
    double dec_abs = std::fabs(dec_deg);
    int dd = (int)dec_abs;
    int dm = (int)((dec_abs - dd) * 60.0);
    double ds = ((dec_abs - dd) * 60.0 - dm) * 60.0;

    std::cout << std::fixed << std::setprecision(1);
    std::cout << "Camera points at:  RA " << rh << "h " << rm << "m " << rs << "s   Dec "
              << (dec_deg < 0 ? "-" : "+") << dd << "d " << dm << "' " << ds << "\"\n";
    std::cout << std::setprecision(5);
    std::cout << "                   RA " << ra * 180.0 / pi << " deg   Dec " << dec_deg << " deg\n";

    // roll: position angle of the picture's "up" direction, measured from north towards east
    vec3d up_dir = make_vec3d(-r.attitude.m[1][0], -r.attitude.m[1][1], -r.attitude.m[1][2]);   // y points down
    double east_norm = sqrt(bore.v[0] * bore.v[0] + bore.v[1] * bore.v[1]);
    if (east_norm > 1.0e-9) {
        vec3d east = make_vec3d(-bore.v[1] / east_norm, bore.v[0] / east_norm, 0.0);
        vec3d north = cross3(bore, east);
        double pa = atan2(dot3(up_dir, east), dot3(up_dir, north)) * 180.0 / pi;
        if (pa < 0.0) pa += 360.0;
        std::cout << std::setprecision(3) << "Roll (position angle of picture up): " << pa << " deg\n";
    }

    std::cout << std::setprecision(6) << "Attitude (camera direction = attitude * sky direction; camera x right, y down, z forward):\n";
    for (int i = 0; i < 3; i++) {
        std::cout << "  " << std::setw(10) << r.attitude.m[i][0] << " " << std::setw(10) << r.attitude.m[i][1]
                  << " " << std::setw(10) << r.attitude.m[i][2] << "\n";
    }
    return 0;
}
