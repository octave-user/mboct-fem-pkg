#include "config.h"

#define NDEBUG

#include <cassert>
#include <cstring>
#include <cstddef>
#include <cstdlib>
#include <cmath>
#include <array>
#include <algorithm>
#ifdef HAVE_WIN32_IO
#include <io.h>
#include <fcntl.h>
#elif defined(HAVE_UNIX_IO)
#include <unistd.h>
#endif

using namespace std;

template <typename T, size_t N>
T interp1(const array<T, N>& xi, const array<T, N>& yi, T x)
{
        constexpr ssize_t Ns = N;
        
        static_assert(N >= 2, "N must be at least two");

        if (x < xi[0]) {
                return yi[0];
        }

        if (x > xi[N - 1]) {
                return yi[N - 1];
        }
        
        ssize_t i2 = 0;

        for (i2 = 1; i2 < Ns; ++i2) {
                assert(xi[i2 - 1] <= xi[i2]);
        }
        
        for (i2 = 1; i2 < Ns; ++i2) {
                if (x <= xi[i2]) {
                        break;
                }
        }

        ssize_t i1 = i2 - 1;

        return (yi[i2] - yi[i1]) / (xi[i2] - xi[i1]) * (x - xi[i1]) + yi[i1];
}

int diaphragm(int argc, char* argv[])
{
        constexpr int nvars = 11;

        if (argc != nvars + 1) {
                fprintf(stderr, "%s: invalid number of arguments\n", argv[0]);
                return 3;
        }
        
        union {
                struct {
                        double Wo, Di, Do, d1, d2, d3, d4, d7, h1, h2, ht;
                };

                double vars[nvars];
        } par;

        memset(&par, 0, sizeof(par));

        static const char format[nvars][6] = {"Wo=%g",
                                              "Di=%g",
                                              "Do=%g",
                                              "d1=%g",
                                              "d2=%g",
                                              "d3=%g",
                                              "d4=%g",
                                              "d7=%g",
                                              "h1=%g",
                                              "h2=%g",
                                              "ht=%g"};
        
        for (int i = 1; i < argc; ++i) {
                int j = 0;
                while (argv[i][j] != '\0') {
                        if (argv[i][j++] == '=') {
                                break;
                        }
                }

                int k = 0;

                while (format[i - 1][k] != '\0') {
                        if (format[i - 1][k++] == '=') {
                                break;
                        }
                }

                if (0 != strncmp(argv[i], format[i - 1], k)) {
                        fprintf(stderr, "invalid argument \"%s\"\n", argv[i]);
                        return 1;
                }

                char* endptr = nullptr;
                
                par.vars[i - 1] = strtod(&argv[i][j], &endptr);

                if (endptr && *endptr) {
                        fprintf(stderr, "invalid argument \"%s\"\n", argv[i]);
                        return 1;
                }
                
                fprintf(stderr, format[i - 1], i - 1, par.vars[i - 1]);
                fputc('\n', stderr);
        }

        constexpr int N = 3;
        array<double, N> X;

        array<double, 4> x1 = {0.5 * par.Di,
                               0.5 * par.Di + par.d1,
                               0.5 * par.Di + par.d1 + par.ht,
                               0.5 * par.Do};

        array<double, 4> y1 = {0, 0, 1, 1};

        array<double, 5> x2 = {0,
                               par.d3 + 0.5 * par.d7,
                               par.d3 + 0.5 * par.d7 + par.ht,
                               par.Di * M_PI - 0.5 * par.d7 - par.ht,
                               par.Di * M_PI};

        array<double, 5> y2 = {0, 0, 1, 1, 0};
        array<double, 6> x3;
        array<double, 6> y3 = {1, 1, 0, 0, 1, 1};
        array<double, 6> x4 = {0,
                               0.5 * par.d7,
                               0.5 * par.d7 + par.d3,
                               0.5 * par.d7 + par.d3 + par.ht,
                               par.Di * M_PI - par.ht,
                               par.Di * M_PI};
        array<double, 6> y4 = {par.d2,
                               par.d2,
                               par.d4,
                               0,
                               0,
                               par.d2};
        array<double, 2> xh = {0, 1};
        array<double, 2> yh = {par.h1, par.h2};
        array<double, 3> f;
        
        while (read(fileno(stdin), &X[0], sizeof(X)) == sizeof(X)) {                
#ifndef NDEBUG
                for (int i = 0; i < N; ++i) {
                        fprintf(stderr, "X[%d]=%g\n", i, X[i]);
                }
#endif                
                for (int i = 0; i < N; ++i) {
                        if (!isfinite(X[i])) {
                                return 0;
                        }
                }

                double r = sqrt(X[0] * X[0] + X[1] * X[1]);
                
#ifndef NDEBUG
                fprintf(stderr, "r=%g\n", r);
#endif
                
                double Phi = atan2(X[1], X[0]);

                if (Phi < 0) {
                        Phi += 2 * M_PI;
                }

#ifndef NDEBUG
                fprintf(stderr, "Phi=%g\n", Phi * 180 / M_PI);
#endif
                
                double x = 0.5 * par.Di * Phi;

                assert(x >= 0);
                assert(x <= par.Di * M_PI);
#ifndef NDEBUG
                fprintf(stderr, "x=%g\n", x);
#endif  
                double z = X[2];
                
                double w = interp1(x4, y4, x);

                assert(w >= 0);
#ifndef NDEBUG
                fprintf(stderr, "w=%g\n", w);
#endif
                
                f[0] = interp1(x1, y1, r);
                f[1] = interp1(x2, y2, x);

                x3 = {-0.5 * par.Wo,
                      -w - par.ht,
                      -w,
                      w,
                      w + par.ht,
                      0.5 * par.Wo};

                f[2] = interp1(x3, y3, z);

#ifndef NDEBUG
                for (size_t i = 0; i < f.size(); ++i) {
                        fprintf(stderr, "f[%ld]=%g\n", i, f[i]);
                }                                
#endif
                
                double h = interp1(xh, yh, *max_element(f.begin(), f.end()));

#ifndef NDEBUG
                fprintf(stderr, "h=%g\n",h);
#endif
                
                if (sizeof(h) != write(fileno(stdout), &h, sizeof(h))) {
                        return 3;
                }
        }
        
        return !feof(stdin) ? 0 : 2;
}

int main(int argc, char* argv[])
{
#ifdef HAVE_WIN32_IO
        setmode(fileno(stdin), _O_BINARY);
        setmode(fileno(stdout), _O_BINARY);
#endif

        if (argc < 2) {
                fprintf(stderr, "usage: fem_pre_mesh_size <geometry> <args>\n");
                return 1;
        }

        static const struct {
                int (*pfn)(int, char**);
                char name[10];
        } funcs[] = {
                {&diaphragm, "diaphragm"}
        };
        
        for (const auto& f:funcs) {
                if (!strcmp(argv[1], f.name)) {
                        return f.pfn(argc - 1, argv + 1);
                }
        }

        fprintf(stderr, "function \"%s\" not found\n", argv[1]);
        
        return 1;
}
