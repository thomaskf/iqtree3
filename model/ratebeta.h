//
// Created by chordata on 8/27/25.
//

#ifndef IQTREE_RATEBETA_H
#define IQTREE_RATEBETA_H
#include <boost/math/distributions/beta.hpp>


class RateBeta {
public:
    RateBeta();

    static double* SampleOmegas(double alpha, double beta);
};

#endif //IQTREE_RATEBETA_H