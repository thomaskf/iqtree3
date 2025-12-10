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
    boost::math::beta_distribution distribution(alpha,beta);
    double *omega = new double[n];
    //omega[0] = (boost::math::quantile(distribution,0.1)+boost::math::quantile(distribution,0.005))/2.0;
    //omega[0] = (boost::math::quantile(distribution,0.05));
    double midpoint = (1.0/n)/2;
    for (int i = 0; i < n; i++) {
        //omega[i] = (boost::math::quantile(distribution,(double(i)/10)+boost::math::quantile(distribution,double(i+1)/10))/2.0;
        omega[i] = (boost::math::quantile(distribution,(double(i)/n)+ midpoint));
        //if (omega[i] == 1.0) {
            //std::cout << std::setprecision(10) << "quantile: " << (double(i)/10)+0.05 << " P: " << p << " Q: " << q << endl;
        //}

    }
    return omega;
}


