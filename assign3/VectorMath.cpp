#include "VectorMath.h"

/**
 * Constructor
 */
VectorMath::VectorMath()
{
    // N/As
}

/**
 * Destructor
 */
VectorMath::~VectorMath()
{
    // N/A
}

/**
 * scalarMultiply - Does Scalar Multiplication
 *                  of a Vector v by Scalar s
 */
VectorMath::point VectorMath::scalarAddition(double s, VectorMath::point v)
{
    // Initialize Scaled Vector
    VectorMath::point sVector;

    sVector.x = s + v.x;
    sVector.y = s + v.y;
    sVector.z = s + v.z;

    return sVector;
}

/**
 * scalarMultiply - Does Scalar Multiplication
 *                  of a Vector v by Scalar s
 */
VectorMath::point VectorMath::scalarMultiply(double s, VectorMath::point v)
{
    // Initialize Scaled Vector
    VectorMath::point sVector;

    sVector.x = s * v.x;
    sVector.y = s * v.y;
    sVector.z = s * v.z;

    return sVector;
}

/**
 * scalarMultiply - Does Scalar Division
 *                  of a Vector v by Scalar s
 */
VectorMath::point VectorMath::scalarDivision(double s, point v)
{
    // Initialize Scaled Vector
    VectorMath::point sVector;

    sVector.x = v.x / s;
    sVector.y = v.y / s;
    sVector.z = v.z / s;

    return sVector;
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
 * normalize - Normalizes a Vector
 *             to get the Unit Vector
 */
VectorMath::point VectorMath::normalize(VectorMath::point v)
{
    // Create Unit Vector
    VectorMath::point unitVector;

    // Calculate Length (Magnitude) of the Vector
    double length = magnitude(v);

    // Check Divide by zero case
    if(length > 0)
    {
        unitVector.x = v.x/length;
        unitVector.y = v.y/length;
        unitVector.z = v.z/length;
    }
    else
    {
        unitVector.x = v.x;
        unitVector.y = v.y;
        unitVector.z = v.z;
    }

    return unitVector;
}

/**
 * magnitude - Gets the Magnitude
 *             of a Vector
 */
double VectorMath::magnitude(VectorMath::point v)
{
    // Calculate Length (Magnitude) of the Vector
    double magnitude = sqrt(pow(v.x,2) + pow(v.y,2) + pow(v.z,2));

    return magnitude;
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

/**
 * crossProduct - Calcuates the Cross Product
 *                betweens Vectors A and B, returning
 *                an Orthogonal ("Normal") Vector C
 *
 * @param a - Vector A
 * @param b - Vector B
 */
VectorMath::point VectorMath::crossProduct(VectorMath::point a, VectorMath::point b)
{
    // Initialize the Return Vector
    VectorMath::point v;

    // Calculate the Cross Product
    v.x = (a.y * b.z) - (a.z * b.y);
    v.y = (a.z * b.x) - (a.x * b.z);
    v.z = (a.x * b.y) - (a.y * b.x);

    return v;
}
