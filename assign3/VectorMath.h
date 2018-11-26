/*
 * VectorMath.h
 *
 *  Created on: Nov 24, 2018
 *      Author: Matt
 */

#ifndef _VECTOR_MATH_H_
#define _VECTOR_MATH_H_

#include <stdlib.h>
#include <math.h>

class VectorMath
{


public:
    /**
     * point - represents a point in OpenGL Space
     */
    struct point
    {
        double x;
        double y;
        double z;
    };

    /**
     * Constructor
     */
    VectorMath();

    /**
     * Destructor
     */
    ~VectorMath();

    /**
     * addVectors - Does Vector Addition
     *              on Two Vectors
     */
    static VectorMath::point addVectors(point v1, point v2);

    /**
     * subtractVectors - Does Vector Subtraction
     *                   on two Vectors
     */
    static VectorMath::point subtractVectors(point v1, point v2);

    /**
     * dotProduct - Performs the Dot Product on Vector A and Vector B
     *
     * param vectorA
     * param vectorB
     */
    static double dotProduct(point vectorA, point vectorB);

    /**
     * quadratic - Performs the Quadratic Equation
     *             on the provided values
     */
    static bool quadratic(double a, double b, double c, double &root1, double &root2);
};

#endif /* _VECTOR_MATH_H_ */
