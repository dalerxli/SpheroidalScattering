
#include <iostream>
#include <ctime>

#include <splinter/datatable.h>
#include <splinter/bspline.h>
#include <splinter/bsplinebuilder.h>

using std::cout;
using std::endl;

using namespace SPLINTER;

// Six-hump camelback function
double f(DenseVector x)
{
    assert(x.rows() == 2);
    return (4 - 2.1*x(0)*x(0)
            + (1/3.)*x(0)*x(0)*x(0)*x(0))*x(0)*x(0)
           + x(0)*x(1)
           + (-4 + 4*x(1)*x(1))*x(1)*x(1);
}

void testInterpolation()
{
    // Create new DataTable to manage samples
    DataTable samples;

    // Sample the function
    DenseVector x(2);
    double y;
    for(int i = 0; i < 20; i++)
    {
        for(int j = 0; j < 20; j++)
        {
            // Sample function at x
            x(0) = i*0.1;
            x(1) = j*0.1;
            y = f(x);

            // Store sample
            samples.addSample(x,y);
        }
    }

    // Build B-splines that interpolate the samples
    BSpline bspline1 = BSpline::Builder(samples).degree(1).build();
    BSpline bspline3 = BSpline::Builder(samples).degree(3).build();

    // Build penalized B-spline (P-spline) that smooths the samples
    BSpline pspline = BSpline::Builder(samples)
            .degree(3)
            .smoothing(BSpline::Smoothing::PSPLINE)
            .alpha(0.03)
            .build();

    /* Evaluate the approximants at x = (1,1)
     * Note that the error will be 0 at that point (except for the P-spline, which may introduce an error
     * in favor of a smooth approximation) because it is a point we sampled at.
     */
    x(0) = 1; x(1) = 1;
    cout << "-----------------------------------------------------" << endl;
    cout << "Function at x:                 " << f(x)               << endl;
    cout << "Linear B-spline at x:          " << bspline1.eval(x)   << endl;
    cout << "Cubic B-spline at x:           " << bspline3.eval(x)   << endl;
    cout << "P-spline at x:                 " << pspline.eval(x)    << endl;
    cout << "-----------------------------------------------------" << endl;
}


double f3D(DenseVector x)
{
    assert(x.rows() == 3);
    return std::sin(x(0) + x(1))*x(2) + std::cos(x(2));
}



void testInterpolation3D()
{
    // Create new DataTable to manage samples
    DataTable samples;

    // Sample the function
    DenseVector x(3);
    double x0_0 = -1.0;
    double x0_1 = 1.0;
    double x1_0 = -1.0;
    double x1_1 = 1.0;
    double x2_0 = -1.0;
    double x2_1 = 1.0;
    double y;

    int nx0 = 7;
    int nx1 = 7;
    int nx2 = 7;

    for(int i = 0; i < nx0; i++)
    {
        for(int j = 0; j < nx1; j++)
        {
            for(int k = 0; k < nx2; k++)
            {
                // Sample function at x
                x(0) = x0_0 + (double)i/(nx0 - 1) * (x0_1 - x0_0);
                x(1) = x1_0 + (double)j/(nx1 - 1) * (x1_1 - x1_0);
                x(2) = x2_0 + (double)k/(nx2 - 1) * (x2_1 - x2_0);
                y = f3D(x);

                // Store sample
                samples.addSample(x,y);
            }
        }
    }

    std::srand(std::time(nullptr));

    // Build B-splines that interpolate the samples
    BSpline bspline = BSpline::Builder(samples).degree(4).build();

    /* Evaluate the approximants at x = (1,1)
     * Note that the error will be 0 at that point (except for the P-spline, which may introduce an error
     * in favor of a smooth approximation) because it is a point we sampled at.
     */
    cout << "-----------------------------------------------------" << endl;
    for(int i = 0; i < 100; ++i) {
        x(0) = (double)std::rand()/RAND_MAX;
        x(1) = (double)std::rand()/RAND_MAX;
        x(2) = (double)std::rand()/RAND_MAX;

        std::cout << "(" << x(0) << "," << x(1) << "," << x(2) << "): " << f3D(x) << "  " << bspline.eval(x) << std::endl;
    }
    cout << "-----------------------------------------------------" << endl;
}

