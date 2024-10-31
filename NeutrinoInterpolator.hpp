#ifndef __NeutrinoInterpolator_hpp__
#define __NeutrinoInterpolator_hpp__

#include <gsl/gsl_spline.h>
#include <gsl/gsl_spline2d.h>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>

class NeutrinoInterpolator
{
private:
    gsl_interp_accel *acc;
    gsl_spline *spline;

public:
    NeutrinoInterpolator(const std::string& file)
    {
        // Check if file exists
        std::ifstream f(file);
        if (!f.is_open())
        {
            std::cerr << "Error: Data file '" << file << "' does not exist.\n";
            throw std::runtime_error("Data file not found.");
        }

        // Read file contents
        std::vector<std::pair<double, double> > array;
        std::string line;
        while (getline(f, line))
        {
            if (line.empty() || line[0] == '#') continue;

            std::stringstream iss(line);
            std::pair<double, double> point;
            iss >> point.first;
            iss.ignore();
            iss >> point.second;
            array.push_back(point);
        }
        f.close();
        unsigned int npoints = array.size();

        // Fill axes
        double x[npoints], y[npoints];
        for (unsigned int i = 0; i < npoints; i++)
        {
            x[i] = array[i].first;
            y[i] = array[i].second;
        }

        // Initialize spline
        acc = gsl_interp_accel_alloc();
        spline = gsl_spline_alloc(gsl_interp_cspline, npoints);
        gsl_spline_init(spline, x, y, npoints);
    }

    // Disable copy constructor and assignment operator
    NeutrinoInterpolator(const NeutrinoInterpolator&) = delete;
    NeutrinoInterpolator& operator=(const NeutrinoInterpolator&) = delete;

    ~NeutrinoInterpolator()
    {
        gsl_spline_free(spline);
        gsl_interp_accel_free(acc);
    }

    double eval(double x)
    {
        return gsl_spline_eval(spline, x, acc);
    }
};

class NeutrinoInterpolator2D
{
private:
    gsl_interp_accel *xacc;
    gsl_interp_accel *yacc;
    gsl_spline2d *spline2d;

public:
    NeutrinoInterpolator2D(const std::string& file)
    {
        // Check if file exists
        std::ifstream f(file);
        if (!f.is_open())
        {
            std::cerr << "Error: Data file '" << file << "' does not exist.\n";
            throw std::runtime_error("Data file not found.");
        }

        std::vector<double> xvals, yvals;
        std::map<std::pair<double, double>, double> zvals;
        double maxzval = 0;
        std::string line;

        while (getline(f, line))
        {
            if (line.empty() || line[0] == '#') continue;

            std::stringstream iss(line);
            double xval, yval, zval;

            iss >> xval;
            if (std::find(xvals.begin(), xvals.end(), xval) == xvals.end())
                xvals.push_back(xval);
            iss.ignore();

            iss >> yval;
            if (std::find(yvals.begin(), yvals.end(), yval) == yvals.end())
                yvals.push_back(yval);
            iss.ignore();

            iss >> zval;
            zvals[std::make_pair(xval, yval)] = zval;
            if (zval > maxzval) maxzval = zval;
        }
        f.close();
        size_t xsize = xvals.size();
        size_t ysize = yvals.size();

        // Fill axes
        double x[xsize], y[ysize], z[xsize * ysize];
        for (size_t j = 0; j < ysize; j++)
        {
            y[j] = yvals[j];
            for (size_t i = 0; i < xsize; i++)
            {
                x[i] = xvals[i];
                auto index = std::make_pair(x[i], y[j]);
                z[j * xsize + i] = zvals.count(index) ? zvals[index] : maxzval;
            }
        }

        // Initialize spline
        xacc = gsl_interp_accel_alloc();
        yacc = gsl_interp_accel_alloc();
        spline2d = gsl_spline2d_alloc(gsl_interp2d_bicubic, xsize, ysize);
        gsl_spline2d_init(spline2d, x, y, z, xsize, ysize);
    }

    // Disable copy constructor and assignment operator
    NeutrinoInterpolator2D(const NeutrinoInterpolator2D&) = delete;
    NeutrinoInterpolator2D& operator=(const NeutrinoInterpolator2D&) = delete;

    ~NeutrinoInterpolator2D()
    {
        gsl_spline2d_free(spline2d);
        gsl_interp_accel_free(xacc);
        gsl_interp_accel_free(yacc);
    }

    double eval(double x, double y)
    {
        return gsl_spline2d_eval(spline2d, x, y, xacc, yacc);
    }
};

#endif // __NeutrinoInterpolator_hpp__