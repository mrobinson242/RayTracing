/*
 * VectorMath.cpp
 *
 *  Created on: Nov 24, 2018
 *      Author: Matt
 */

#include "VectorMath.h"

/**
 * Constructor
 */
VectorMath::VectorMath()
{
    //N/A
}

/**
 * Destructor
 */
VectorMath::~VectorMath()
{
    // N/A
}

/**
 * addVectors - Does Vector Addition
 *              on Two Vectors
 */
VectorMath::point VectorMath::addVectors(point v1, point v2)
{
    // Initialize Added Vector
    point returnVector;

    returnVector.x = v1.x + v2.x;
    returnVector.y = v1.y + v2.y;
    returnVector.z = v1.z + v2.z;

    return returnVector;
}

/**
 * subtractVectors - Does Vector Subtraction
 *                   on two Vectors
 */
VectorMath::point VectorMath::subtractVectors(VectorMath::point v1, VectorMath::point v2)
{
    // Initialize Subtracted Vector
    VectorMath::point returnVector;

    returnVector.x = v1.x - v2.x;
    returnVector.y = v1.y - v2.y;
    returnVector.z = v1.z - v2.z;

    return returnVector;
}

/**
 * dotProduct - Performs the Dot Product on Vector A and Vector B
 *              (Projection of Vector A onto Vector B)
 *
 * param v1
 * param v2
 */
double VectorMath::dotProduct(VectorMath::point v1, VectorMath::point v2)
{
    // Initialize the Return Vector
    double returnVal;

    // Calculate the Dot Product between the two Vectors
    returnVal = (v1.x * v2.x) + (v1.y * v2.y) + (v1.z * v2.z);

    return returnVal;
}

/**
 * quadratic - Performs the Quadratic Equation
 *             on the provided values
 */
bool VectorMath::quadratic(double a, double b, double c, double &root1, double &root2)
{
    // Initialize Return Indicator
    bool isValid = true;

    // Calculate the Descriminant
    double discriminant = (b*b) - (4*a*c);

    // Check if Determinant is Negative
    if(discriminant < 0)
    {
        // Equation is Invalid (Divide by 0)
        // There will be an imaginary root
        isValid = false;
    }
    // Check if Determinant is 0
    else if(discriminant == 0)
    {
        // Update Roots
        root1 = -0.5 * (b/a);
        root2 = root1;
    }
    // Check if Determinant is Positive
    else
    {
        // Update Real Roots
        root1 = (-b + sqrt(discriminant)) / (2*a);
        root2 = (-b - sqrt(discriminant)) / (2*a);
    }

    return isValid;
}
