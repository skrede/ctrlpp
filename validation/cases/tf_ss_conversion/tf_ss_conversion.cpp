// tf_ss_conversion.cpp -- TF<->SS round-trip matching tf_ss_conversion.m
// Compares poles and round-trip TF coefficients (canonical form invariant).
// Usage: ./tf_ss_conversion > tf_ss_conversion_cpp.csv

#include "ctrlpp/model/conversion.h"
#include "ctrlpp/model/transfer_function.h"

#include <Eigen/Eigenvalues>

#include <algorithm>
#include <complex>
#include <cstdio>

int main()
{
    using Scalar = double;

    ctrlpp::transfer_function<Scalar, 1, 2> tf{};
    tf.numerator = {2.0, 3.0};
    tf.denominator = {1.0, 0.5, 1.0};

    auto sys = ctrlpp::tf2ss(tf);

    // Poles from A
    Eigen::EigenSolver<Eigen::Matrix<Scalar, 2, 2>> es(sys.A);
    auto eigs = es.eigenvalues();
    std::array<std::complex<Scalar>, 2> poles = {eigs(0), eigs(1)};
    std::sort(poles.begin(), poles.end(), [](auto a, auto b) {
        return a.real() < b.real() || (a.real() == b.real() && a.imag() < b.imag());
    });

    // Round-trip
    auto tf_rt = ctrlpp::ss2tf<Scalar, 2>(sys);

    Scalar d0 = tf_rt.denominator[0];
    Scalar num0 = tf_rt.numerator[0] / d0;
    Scalar num1 = tf_rt.numerator[1] / d0;
    Scalar num2 = tf_rt.numerator[2] / d0;
    Scalar den0 = 1.0;
    Scalar den1 = tf_rt.denominator[1] / d0;
    Scalar den2 = tf_rt.denominator[2] / d0;

    std::printf("pole0_re,pole0_im,pole1_re,pole1_im,D_00,num_rt_0,num_rt_1,den_rt_0,den_rt_1,den_rt_2\n");
    std::printf("%.15e,%.15e,%.15e,%.15e,%.15e,%.15e,%.15e,%.15e,%.15e,%.15e\n",
                poles[0].real(), poles[0].imag(),
                poles[1].real(), poles[1].imag(),
                sys.D(0, 0),
                num0, num1,
                den0, den1, den2);

    (void)num2;
}
