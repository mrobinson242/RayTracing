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
     * scalarMultiply - Does Scalar Multiplication
     *                  of a Vector v by Scalar s
     */
    static VectorMath::point scalarAddition(double s, point v);

    /**
     * scalarMultiply - Does Scalar Multiplication
     *                  of a Vector v by Scalar s
     */
    static VectorMath::point scalarMultiply(double s, point v);

    /**
     * scalarMultiply - Does Scalar Division
     *                  of a Vector v by Scalar s
     */
    static VectorMath::point scalarDivision(double s, point v);

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
     * normalize - Normalizes a Vector
     *             to get the Unit Vector
     */
    static VectorMath::point normalize(point v);

    /**
     * magnitude - Gets the Magnitude
     *             of a Vector
     */
    static double magnitude(point v);

    /**
     * dotProduct - Performs the Dot Product on Vector A and Vector B
     *
     * param vectorA
     * param vectorB
     */
    static double dotProduct(point vectorA, point vectorB);

    /**
     * crossProduct - Calcuates the Cross Product
     *                betweens Vectors A and B, returning
     *                an Orthogonal ("Normal") Vector C
     *
     * @param a - Vector A
     * @param b - Vector B
     */
    static VectorMath::point crossProduct(point a, point b);

    /**
     * quadratic - Performs the Quadratic Equation
     *             on the provided values
     */
    static bool quadratic(double a, double b, double c, double &root1, double &root2);

private:

};

#endif /* _VECTOR_MATH_H_ */
