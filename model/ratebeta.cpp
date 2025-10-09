//
// Created by chordata on 8/27/25.
//

#include "ratebeta.h"
#include <boost/math/distributions/beta.hpp>
#include <iostream>

#include "ncl/nxsstring.h"

RateBeta::RateBeta(){
}

double* RateBeta::SampleOmegas(double p, double q) {
    boost::math::beta_distribution distribution(p,q);
    double *omega = new double[10];
    //omega[0] = (boost::math::quantile(distribution,0.1)+boost::math::quantile(distribution,0.005))/2.0;
    omega[0] = (boost::math::quantile(distribution,0.05));
    for (int i = 1; i < 10; i++) {
        //omega[i] = (boost::math::quantile(distribution,(double(i)/10)+boost::math::quantile(distribution,double(i+1)/10))/2.0;
        omega[i] = (boost::math::quantile(distribution,(double(i)/10)+0.05));
        //if (omega[i] == 1.0) {
            //std::cout << std::setprecision(10) << "quantile: " << (double(i)/10)+0.05 << " P: " << p << " Q: " << q << endl;
        //}

    }
    return omega;
}


