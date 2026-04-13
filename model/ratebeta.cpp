//
// Created by chordata on 8/27/25.
//

#include "ratebeta.h"
#include <boost/math/distributions/beta.hpp>
#include <iostream>

#include "ncl/nxsstring.h"

RateBeta::RateBeta(){
}

double* RateBeta::SampleOmegas(int n,double alpha, double beta) {
    boost::math::beta_distribution<double> distribution(alpha,beta);
    double *omega = new double[n];
    double midpoint = (1.0/n)/2;
    for (int i = 0; i < n; i++) {
        omega[i] = (boost::math::quantile(distribution,(double(i)/n)+ midpoint));
    }
    return omega;
}

