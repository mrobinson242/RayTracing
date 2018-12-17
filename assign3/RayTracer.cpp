/*                                            */
/*  CSCI 420 Computer Graphics                */
/*  Assignment 3: RayTracer                   */
/*  Author: Matt Robinson                     */
/*  Student Id: 9801107811                    */
/*                                            */

/** Open GL Libraries */
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>

/** Standard Libraries */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sstream>
#include <iostream>
#include <vector>
#include <string.h>

/** Custom Libraries */
#include <pic.h>
#include "VectorMath.h"

/** Constants */
#define MAX_TRIANGLES 2000
#define MAX_SPHERES 10
#define MAX_LIGHTS 10
#define epsilon 0.000000001
#define MAX_DISTANCE -1000000.0;

/** Filename */
char *filename=0;

/** Different display modes */
#define MODE_DISPLAY 1
#define MODE_JPEG 2
int mode = MODE_DISPLAY;

/** Width/Height of the Image */
#define WIDTH 640
#define HEIGHT 480

/** Field of View of the camera */
#define fov 60.0

/** Memory Buffer for the Image */
unsigned char buffer[HEIGHT][WIDTH][3];

/** Sphere Structure */
typedef struct _Sphere
{
    double position[3];
    double color_diffuse[3];
    double color_specular[3];
    double shininess;
    double radius;
} Sphere;

/** Light Structure */
typedef struct _Light
{
    double position[3];
    double color[3];
} Light;

/** Vertex Structure */
struct Vertex
{
    double position[3];
    double color_diffuse[3];
    double color_specular[3];
    double normal[3];
    double shininess;
};

/** Triangle Structure */
typedef struct _Triangle
{
    struct Vertex v[3];
} Triangle;

/** Color (char) */
struct Color
{
    char red;
    char green;
    char blue;
};

/** Ray Structure */
struct Ray
{
    VectorMath::point origin;
    VectorMath::point direction;
};

/** Color (double) */
struct RayColor
{
    double r;
    double g;
    double b;

    RayColor& operator += (const RayColor& other)
    {
        // Add Red Intensity
        r += other.r;

        // Bound the Red Intensity
        if(r > 1.0)
        {
            r = 1.0;
        }

        // Add Green Intensity
        g += other.g;

        // Bound the Green Intensity
        if(g > 1.0)
        {
            g = 1.0;
        }

        // Add Blue Intensity
        b += other.b;

        // Bound the Blue Intensity
        if(b > 1.0)
        {
            b = 1.0;
        }

        return *this;
    }
};

// Initialize Triangle Data Struct
Triangle _triangles[MAX_TRIANGLES];

// Initialize Sphere Data Struct
Sphere _spheres[MAX_SPHERES];

// Initialize Lights Data Struct
Light _lights[MAX_LIGHTS];

// Ambient Light
double _ambientLight[3];

// Object Count Indicators
int _numTriangles = 0;
int _numSpheres = 0;
int _numLights = 0;

//hack to make it only draw once
static int once = 0;

// Initialize Camera Origin Point for each Rays
VectorMath::point _cameraOrigin = {0.0, 0.0, 0.0};

/** Define Functions */
void plotPixelDisplay(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plotPixelJpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plotPixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);

/**
 * pointInTriangle - Checks if intersection point is in triangle
 *
 * param triangle - The Triangle
 * param p - The Intersection Point
 * param n - The Normal Vector of the Triangle
 */
bool pointInTriangle(Triangle triangle, VectorMath::point p, VectorMath::point n)
{
    // Initialize Inside Triangle Indicator
    bool isInTriangle = true;

    // Get Triangle Vertices
    Vertex a = triangle.v[0];
    Vertex b = triangle.v[1];
    Vertex c = triangle.v[2];

    // Calculate Edges of Triangle
    VectorMath::point edgeCB = {c.position[0] - b.position[0],
                                c.position[1] - b.position[1],
                                c.position[2] - b.position[2]};

    VectorMath::point edgeBA = {b.position[0] - a.position[0],
                                b.position[1] - a.position[1],
                                b.position[2] - a.position[2]};

    VectorMath::point edgeAC = {a.position[0] - c.position[0],
                                a.position[1] - c.position[1],
                                a.position[2] - c.position[2]};

    // Calculate Vector between Intersection Point and each Vertex
    VectorMath::point vInterA = {p.x - a.position[0],
                                 p.y - a.position[1],
                                 p.z - a.position[2]};

    VectorMath::point vInterB = {p.x - b.position[0],
                                 p.y - b.position[1],
                                 p.z - b.position[2]};

    VectorMath::point vInterC = {p.x - c.position[0],
                                 p.y - c.position[1],
                                 p.z - c.position[2]};

    // Calculate Cross Product between each Edge and the Intersection/Vertex Vector
    VectorMath::point crossA = VectorMath::crossProduct(edgeBA, vInterA);
    VectorMath::point crossB = VectorMath::crossProduct(edgeCB, vInterB);
    VectorMath::point crossC = VectorMath::crossProduct(edgeAC, vInterC);

    // Check if point is inside Triangle
    if( (VectorMath::dotProduct(n, crossA) < 0) ||
        (VectorMath::dotProduct(n, crossB) < 0) ||
        (VectorMath::dotProduct(n, crossC) < 0))
    {
        // No Intersection
        isInTriangle = false;
    }

    return isInTriangle;
}

/**
 * intersectSphere - Checks if the Ray intersects with the Sphere
 *
 * @param ray - The Ray
 * @param center - The center of the sphere
 * @param radius - The radius of the sphere
 * @param interPoint - The Intersection Point
 */
bool intersectSphere(Ray ray, VectorMath::point center, double radius, VectorMath::point &interPoint)
{
    // Initialize Intersection Indicator
    bool intersectsSphere = true;

    // Create Intersection Values
    double t0 = -1.0;
    double t1 = -1.0;
    double t = -1.0;

    VectorMath::point l = VectorMath::subtractVectors(ray.origin, center);

    // Obtain Scalar Values from Quadratic Equation
    double a = VectorMath::dotProduct(ray.direction, ray.direction);
    double b = 2 * (VectorMath::dotProduct(ray.direction, l));
    double c = VectorMath::dotProduct(l, l) - pow(radius, 2);

    // Check if valid intersection points
    if(VectorMath::quadratic(a, b, c, t0, t1))
    {
        // Check if at least one intersection
        if(t0 >= epsilon || t1 >= epsilon)
        {
            if(t0 >= 0.0 && t1 >= 0.0)
            {
                t = fmin(t0, t1);
            }
            else if(t0 >= 0.0 && t1 <= 0.0)
            {
                t = t0;
            }
            else if(t1 >= 0.0 && t0 <= 0.0)
            {
                t = t1;
            }
        }
        else
        {
            // Does not Intersect Sphere
            intersectsSphere = false;
        }
    }
    else
    {
        // Does not Intersect Sphere
        intersectsSphere = false;
    }

    // Calculate Intersection Point [p = v0 + (vd*t)]
    VectorMath::point sphereIntersection = VectorMath::addVectors(ray.origin, VectorMath::scalarMultiply(t, ray.direction));
    interPoint.x = sphereIntersection.x;
    interPoint.y = sphereIntersection.y;
    interPoint.z = sphereIntersection.z;

    return intersectsSphere;
}

/**
 * intersectTriangle - Checks if the ray intersects with the Triangle
 *
 * @param triangle - Triangle to Check
 * @param ray - The Ray
 * @param interPoint - The Intersection Point
 */
bool intersectTriangle(Triangle triangle, Ray ray, VectorMath::point &interPoint)
{
    // Initialize Intersection Indicator
    bool intersectsTriangle = true;

    // Get Origin Vector
    VectorMath::point origin = ray.origin;

    // Get Direction Vector
    VectorMath::point direction = ray.direction;

    // Get Triangle Vertices
    Vertex a = triangle.v[0];
    Vertex b = triangle.v[1];
    Vertex c = triangle.v[2];

    // Calculate Edge BA
    VectorMath::point ba = {b.position[0] - a.position[0],
                            b.position[1] - a.position[1],
                            b.position[2] - a.position[2]};

    // Calculate Edge CA
    VectorMath::point ca = {c.position[0] - a.position[0],
                            c.position[1] - a.position[1],
                            c.position[2] - a.position[2]};

    // Calculated Normalized Vector for Ray-Plane Intersection
    VectorMath::point n = VectorMath::normalize(VectorMath::crossProduct(ba, ca));

    // Step 2: Check if Ray Parallel to Plane
    double nDotD = VectorMath::dotProduct(n, direction);

    // Check if [n dot d] is close to 0
    if(abs(nDotD) < epsilon)
    {
        // Ray is Parallel to Plane
        intersectsTriangle = false;
    }
    else
    {
        // Step 3: Calculate d Vector
        VectorMath::point vertexA = {a.position[0], a.position[1], a.position[2]};
        VectorMath::point dist =  VectorMath::subtractVectors(vertexA, origin);
        double d = VectorMath::dotProduct(dist, n);

        // Step 4: Calculate Intersection Point t [t = -(n dot p0 + d) / (n dot d)]
        double t =  (d / nDotD);

        // Step 5: Check if Intersection is behind Ray Origin
        if(t < 0)
        {
            // Intersection is behind Ray Origin
            intersectsTriangle = false;
        }
        else
        {
            // Compute the Intersection Point [P = origin + (t * Direction)]
            VectorMath::point intersection = VectorMath::addVectors(origin, VectorMath::scalarMultiply(t, direction));
            interPoint.x = intersection.x;
            interPoint.y = intersection.y;
            interPoint.z = intersection.z;

            // Check if point not in triangle
            if(!pointInTriangle(triangle, intersection, n))
            {
                // No Triangle Intersection
                intersectsTriangle = false;
            }
        }
    }

    return intersectsTriangle;
}

/**
 * computeDiffuse - Calculates the Diffuse Component for Illumination
 *
 * @param kD   - Diffuse Material Coefficient
 * @param lVec - Unit Vector to Light Source
 * @param nVec - Unit Surface Normal
 */
double computeDiffuse(double kD, VectorMath::point lVec, VectorMath::point nVec)
{
    // Calculate [l dot n]
    double lDotN = VectorMath::dotProduct(lVec, nVec);

    // If L dot N is Less than or close to 0
    if(lDotN < epsilon)
    {
        // Clamp to 0
        lDotN = 0.0;
    }

    // Calculate the Diffuse Component
    double diffuse =  kD * lDotN;

    return diffuse;
}

/**
 * computeSpecular - Calculates the Specular Component for Illumination
 *
 * @param kS    - Specular Material Coefficient
 * @param alpha - The Shininess Coefficient
 * @param rVec  - Unit Vector of Reflected Light
 * @param vVec  - Unit Vector to Camera
 */
double computeSpecular(double kS, double alpha, VectorMath::point rVec, VectorMath::point vVec)
{
    // Calculate [r dot v]
    double rDotV = VectorMath::dotProduct(rVec, vVec);

    if(rDotV < epsilon)
    {
        rDotV = 0.0;
    }

    // Calculate the Specular Component
    double specular = kS * pow(rDotV, alpha);

    return specular;
}

/**
 * computeRVec - Computes the Unit Vector of the Reflected Light
 *
 * @param lVec - Unit Vector to Light Source
 * @parma nVec - Unit Surface Normal
 */
VectorMath::point computeRVec(VectorMath::point lVec, VectorMath::point nVec)
{
    // Calculate [l dot N]
    double lDotN = VectorMath::dotProduct(lVec, nVec);

    // If [l dot N]is Less than or close to 0
    if(lDotN < epsilon)
    {
        // Clamp to 0
        lDotN = 0.0;
    }

    // Calculate [2(l dot N)]
    lDotN *= 2.0;

    // Calculate [2(l dot N)N]
    VectorMath::point tempVec = VectorMath::scalarMultiply(lDotN, nVec);

    // Calculate R Vector [2(l dot N)N - l]
    VectorMath::point rVec = VectorMath::subtractVectors(tempVec, lVec);

    return rVec;
}

/**
 * computeSphereIllumination - Calculates the Colors to display at the Sphere
 *                             based on the Phong Illumination Model
 *
 * param s          - The Sphere
 * param light      - The Light
 * param interPoint - The Intersection Point of the Viewing Vector with the Sphere
 */
RayColor computeSphereIllumination(Sphere s, Light light, VectorMath::point interPoint)
{
    // Get the Diffuse for each Color Channel
    double kDRed = s.color_diffuse[0];
    double kDGreen = s.color_diffuse[1];
    double kDBlue = s.color_diffuse[2];

    // Get the Specular for each Color Channel
    double kSRed = s.color_specular[0];
    double kSGreen = s.color_specular[1];
    double kSBlue = s.color_specular[2];

    // Get the Shininess Coefficient
    double alpha = s.shininess;

    // Calculate V Vector (Vector from Intersection Point to Camera)
    VectorMath::point vVec = VectorMath::subtractVectors(_cameraOrigin, interPoint);
    vVec = VectorMath::normalize(vVec);

    // Calculate N Vector
    VectorMath::point center = {s.position[0], s.position[1], s.position[2]};
    VectorMath::point nVec = VectorMath::scalarDivision(s.radius, VectorMath::subtractVectors(interPoint, center));

    // Get Light Position
    VectorMath::point lightPos = {light.position[0],
                                  light.position[1],
                                  light.position[2]};

    // Calculate L Vector (Vector from Intersection Point to Light)
    VectorMath::point lVec = VectorMath::normalize(VectorMath::subtractVectors(lightPos, interPoint));

    // Calculate R Vector
    VectorMath::point rVec = computeRVec(lVec, nVec);

    // Calculate the Red Intensity
    double diffuseRed = computeDiffuse(kDRed, lVec, nVec);
    double specularRed = computeSpecular(kSRed, alpha, rVec, vVec);
    double iRed = light.color[0] * (diffuseRed + specularRed);

    // Calculate the Green Intensity
    double diffuseGreen = computeDiffuse(kDGreen, lVec, nVec);
    double specularGreen = computeSpecular(kSGreen, alpha, rVec, vVec);
    double iGreen = light.color[1] * (diffuseGreen + specularGreen);

    // Calculate the Blue Intensity
    double diffuseBlue = computeDiffuse(kDBlue, lVec, nVec);
    double specularBlue = computeSpecular(kSBlue, alpha, rVec, vVec);
    double iBlue = light.color[2] * (diffuseBlue + specularBlue);

    // Set Color
    RayColor c = {iRed, iGreen, iBlue};

    return c;
}

/**
 * computeBarycentricCoordinates - Calculates the Barycentric Coordinates used
 *                                 to assist with the Triangle Interpolation
 *
 * @param t - The Triangle
 * @param interPoint - Intersection Point with the Triangle
 * @param alpha - Barycentric Coordinate Alpha
 * @param beta  - Barycentric Coordinate Beta
 * @param gamma - Barycentric Coordinate Gamma
 */
void computeBarycentricCoordinates(Triangle t, VectorMath::point interPoint, double &alpha, double &beta, double &gamma)
{
    // Get Triangle Vertices
    Vertex vA = t.v[0];
    Vertex vB = t.v[1];
    Vertex vC = t.v[2];

    // Convert Vertices into Vector Points
    VectorMath::point p1 = {vA.position[0], vA.position[1], vA.position[2]};
    VectorMath::point p2 = {vB.position[0], vB.position[1], vB.position[2]};
    VectorMath::point p3 = {vC.position[0], vC.position[1], vC.position[2]};

    // Calculate Edge BA
    VectorMath::point ba = {vB.position[0] - vA.position[0],
                            vB.position[1] - vA.position[1],
                            vB.position[2] - vA.position[2]};

    // Calculate Edge CA
    VectorMath::point ca = {vC.position[0] - vA.position[0],
                            vC.position[1] - vA.position[1],
                            vC.position[2] - vA.position[2]};

    // Calculate Edges of Triangle
    VectorMath::point edgeCB = {vC.position[0] - vB.position[0],
                                vC.position[1] - vB.position[1],
                                vC.position[2] - vB.position[2]};

    VectorMath::point edgeAC = {vA.position[0] - vC.position[0],
                                vA.position[1] - vC.position[1],
                                vA.position[2] - vC.position[2]};


    VectorMath::point vInterB = {interPoint.x - vB.position[0],
                                 interPoint.y - vB.position[1],
                                 interPoint.z - vB.position[2]};

    VectorMath::point vInterC = {interPoint.x - vC.position[0],
                                 interPoint.y - vC.position[1],
                                 interPoint.z - vC.position[2]};

    // Calculate Cross Product between each Edge and the Intersection/Vertex Vector
    VectorMath::point crossB = VectorMath::crossProduct(edgeCB, vInterB);
    VectorMath::point crossC = VectorMath::crossProduct(edgeAC, vInterC);

    VectorMath::point planar = VectorMath::crossProduct(ba, ca);
    double denominator = VectorMath::dotProduct(planar, planar);

    // Compute Barycentric coordinates
    alpha = VectorMath::dotProduct(planar, crossB) / denominator;
    beta = VectorMath::dotProduct(planar, crossC) / denominator;
    gamma = 1.0 - alpha - beta;
}

/**
 * computeTriangleIllumination - Calculates the Colors to display at the Triangle
 *                                based on the Phong Illumination Model
 *
 * param t          - The Triangle
 * param light      - The Light Source
 * param interPoint - The Intersection Point of the Viewing Vector with the Triangle
 */
RayColor computeTriangleIllumination(Triangle t, Light light, VectorMath::point interPoint)
{
    // Initialize Barycentric Coordinates
    double alpha;
    double beta;
    double gamma;

    // Calculate Barycentric Coordinates
    computeBarycentricCoordinates(t, interPoint, alpha, beta, gamma);

    // Calculate V Vector (Vector from Intersection Point to Camera)
    VectorMath::point vVec = VectorMath::subtractVectors(_cameraOrigin, interPoint);
    vVec = VectorMath::normalize(vVec);

    // Calculate N Vector
    VectorMath::point nVec;
    nVec.x = (alpha * t.v[0].normal[0]) + (beta * t.v[1].normal[0]) + (gamma * t.v[2].normal[0]);
    nVec.y = (alpha * t.v[0].normal[1]) + (beta * t.v[1].normal[1]) + (gamma * t.v[2].normal[1]);
    nVec.z = (alpha * t.v[0].normal[2]) + (beta * t.v[1].normal[2]) + (gamma * t.v[2].normal[2]);
    nVec = VectorMath::normalize(nVec);

    // Calculate Diffuse
    double kDRed =   (alpha * t.v[0].color_diffuse[0]) + (beta * t.v[1].color_diffuse[0]) + (gamma * t.v[2].color_diffuse[0]);
    double kDGreen = (alpha * t.v[0].color_diffuse[1]) + (beta * t.v[1].color_diffuse[1]) + (gamma * t.v[2].color_diffuse[1]);
    double kDBlue =  (alpha * t.v[0].color_diffuse[2]) + (beta * t.v[1].color_diffuse[2]) + (gamma * t.v[2].color_diffuse[2]);

    // Calculate Specular
    double kSRed =   (alpha * t.v[0].color_specular[0]) + (beta * t.v[1].color_specular[0]) + (gamma * t.v[2].color_specular[0]);
    double kSGreen = (alpha * t.v[0].color_specular[1]) + (beta * t.v[1].color_specular[1]) + (gamma * t.v[2].color_specular[1]);
    double kSBlue =  (alpha * t.v[0].color_specular[2]) + (beta * t.v[1].color_specular[2]) + (gamma * t.v[2].color_specular[2]);

    // Calculate Shininess Value
    double shininess = (alpha * t.v[0].shininess) + (beta * t.v[1].shininess) + (gamma * t.v[2].shininess);

    // Get Light Position
    VectorMath::point lightPos = {light.position[0],
                                  light.position[1],
                                  light.position[2]};

    // Calculate L Vector (Vector from Intersection Point to Light)
    VectorMath::point lVec = VectorMath::normalize(VectorMath::subtractVectors(lightPos, interPoint));

    // Calculate R Vector
    VectorMath::point rVec = VectorMath::normalize(computeRVec(lVec, nVec));

    // Calculate the Red Intensity
    double diffuseRed = computeDiffuse(kDRed, lVec, nVec);
    double specularRed = computeSpecular(kSRed, shininess, rVec, vVec);
    double iRed = light.color[0] * (diffuseRed + specularRed);

    // Calculate the Green Intensity
    double diffuseGreen = computeDiffuse(kDGreen, lVec, nVec);
    double specularGreen = computeSpecular(kSGreen, shininess, rVec, vVec);
    double iGreen = light.color[1] * (diffuseGreen + specularGreen);

    // Calculate the Blue Intensity
    double diffuseBlue = computeDiffuse(kDBlue, lVec, nVec);
    double specularBlue = computeSpecular(kSBlue, shininess, rVec, vVec);
    double iBlue = light.color[2] * (diffuseBlue + specularBlue);

    // Set Color
    RayColor c = {iRed, iGreen, iBlue};

    return c;
}

/**
 * calcSphereShadowRay - Calculate the Shadow Rays hitting
 *                       all the Spheres
 *
 * @param interPoint    - Intersection Point of the first Ray
 * @param shadowRay     - The Shadow Ray
 * @param lightPos      - Position of the Light Source
 * @param sphereIndex   - Index of Initial Sphere
 */
bool calcSphereShadowRay(VectorMath::point interPoint, Ray shadowRay, VectorMath::point lightPos, int sphereIndex)
{
    // Check if Illuminated
    bool isIlluminated = true;

    // Iterate over all the Spheres
    for(int k = 0; k < _numSpheres; ++k)
    {
        // Ignore own object
        if(k != sphereIndex)
        {
            // Get the Sphere
            Sphere s = _spheres[k];

            // Get Center Point of Sphere
            VectorMath::point center = {s.position[0], s.position[1], s.position[2]};

            // Create Intersection Point
            VectorMath::point shadowInterPoint;

            // Check if Intersection
            if(intersectSphere(shadowRay, center, s.radius, shadowInterPoint))
            {
                // Calculate Distance Vectors
                VectorMath::point dist1 = VectorMath::subtractVectors(shadowInterPoint, interPoint);
                VectorMath::point dist2 = VectorMath::subtractVectors(lightPos, interPoint);

                // Get the Magnitude of the Distance Vectors
                double mag1 = VectorMath::magnitude(dist1);
                double mag2 = VectorMath::magnitude(dist2);

                // Ensure Shadow Intersection Point is not past the light
                if(mag1 < mag2)
                {
                    isIlluminated = false;
                    break;
                }
            }
        }
    }

    return isIlluminated;
}

/**
 * calcTriangleShadowRay - Calculate the Shadow Rays hitting
 *                         all the Triangles
 *
 * interPoint    - Intersection Point
 * shadowRay     - The Shadow Ray
 * lightPos      - Position of the Light Source
 * triangleIndex - Index of Initial Triangle
 */
bool calcTriangleShadowRay(VectorMath::point interPoint, Ray shadowRay, VectorMath::point lightPos, int triangleIndex)
{
    // Check if Illuminated
    bool isIlluminated = true;

    // Iterate over all the Triangles
    for(int k = 0; k < _numTriangles; ++k)
    {
        // Ignore own object
        if(k != triangleIndex)
        {
            // Create Intersection Point
            VectorMath::point shadowInterPoint;

            // Get Current Triangle
            Triangle triangle = _triangles[k];

            // Check if Ray Intersects Triangle
            if(intersectTriangle(triangle, shadowRay, shadowInterPoint))
            {
                // Calculate Distance Vectors
                VectorMath::point dist1 = VectorMath::subtractVectors(shadowInterPoint, interPoint);
                VectorMath::point dist2 = VectorMath::subtractVectors(lightPos, interPoint);

                // Get the Magnitude of the Distance Vectors
                double mag1 = VectorMath::magnitude(dist1);
                double mag2 = VectorMath::magnitude(dist2);

                // Ensure Shadow Intersection Point is not past the light
                if(mag1 < mag2)
                {
                    isIlluminated = false;
                    break;
                }
            }
        }
    }

    return isIlluminated;
}

/**
 * processTriangles - Processes all the Ray Intersections
 *                    that hit a Triangle in the Scene
 *
 * @param ray          - The Ray shot from the Camera
 * @param rayColor     - Color of the Pixel according to the Ray
 * @param closestPoint - The Closest Point the ray interacts with
 */
RayColor processTriangles(Ray ray, RayColor& rayColor, double& closestPoint)
{
    // Initialize Color
    RayColor c = rayColor;

    // Iterate over all the Triangles
    for(int i = 0; i < _numTriangles; ++i)
    {
        // Get the Triangle
        Triangle triangle = _triangles[i];

        // Create Intersection Point
        VectorMath::point intersection;

        // Check if Ray Intersects Triangle
        if(intersectTriangle(triangle, ray, intersection) && (intersection.z > closestPoint))
        {
            // If Intersection, default Color to black
            c.r = 0.0;
            c.g = 0.0;
            c.b = 0.0;

            // Iterate over the Lights in the Scene
            for(int j = 0; j < _numLights; ++j)
            {
                // Get Light Position
                VectorMath::point light = {_lights[j].position[0],
                                           _lights[j].position[1],
                                           _lights[j].position[2]};

                // Calculate the Shadow Ray
                VectorMath::point shadowDirection = VectorMath::subtractVectors(light, intersection);
                Ray shadowRay = {intersection, VectorMath::normalize(shadowDirection)};

                // Check if Shadow Ray hits a Sphere
                bool isSphereIlluminated = calcSphereShadowRay(intersection, shadowRay, light, -1);

                // Check if Shadow Ray hits a Triangle
                bool isTriangleIlluminated = calcTriangleShadowRay(intersection, shadowRay, light, i);

                // Check if Illuminated
                if(isSphereIlluminated && isTriangleIlluminated)
                {
                    // Calculate Phong Illumination for the Triangle
                    RayColor triangleColor = computeTriangleIllumination(triangle, _lights[j], intersection);

                    // Update Color
                    c.r += triangleColor.r;
                    c.g += triangleColor.g;
                    c.b += triangleColor.b;
                }
            }

            // Update Closest Point
            closestPoint = intersection.z;
        }
    }

    return c;
}

/**
 * processSpheres - Processes all the Ray Intersections
 *                  that hit a Sphere in the Scene
 *
 * @param ray          - The Ray shot from the Camera
 * @param rayColor     - Color of the Pixel according to the Ray
 * @param closestPoint - The Closest Point the ray interacts with
 */
RayColor processSpheres(Ray ray, RayColor& rayColor, double& closestPoint)
{
    // Initialize Color
    RayColor c = rayColor;

    // Iterate over all the Spheres
    for(int i = 0; i < _numSpheres; ++i)
    {
        // Get the Sphere
        Sphere s = _spheres[i];

        // Get Center Point of Sphere
        VectorMath::point center = {s.position[0], s.position[1], s.position[2]};

        // Create Intersection Point
        VectorMath::point intersection;

        // Check if Ray Intersects Sphere
        if(intersectSphere(ray, center, s.radius, intersection) && (intersection.z > closestPoint))
        {
            // If Intersection, default Color to black
            c.r = 0.0;
            c.g = 0.0;
            c.b = 0.0;

            // Iterate over the Lights in the Scene
            for(int j = 0; j < _numLights; ++j)
            {
                // Get Light Position
                VectorMath::point light = {_lights[i].position[0],
                                           _lights[i].position[1],
                                           _lights[i].position[2]};

                // Calculate the Shadow Ray
                VectorMath::point shadowDirection = VectorMath::subtractVectors(light, intersection);
                Ray shadowRay = {intersection, VectorMath::normalize(shadowDirection)};

                // Check if Shadow Ray hits a Sphere
                bool isSphereIlluminated = calcSphereShadowRay(intersection, shadowRay, light, i);

                // Check if Shadow Ray hits a Triangle
                bool isTriangleIlluminated = calcTriangleShadowRay(intersection, shadowRay, light, -1);

                // Check if Illuminated
                if(isSphereIlluminated && isTriangleIlluminated)
                {
                    // Calculate Phong Illumination for the Sphere
                    RayColor sphereColor = computeSphereIllumination(s, _lights[j], intersection);

                    // Update Color
                    c.r += sphereColor.r;
                    c.g += sphereColor.g;
                    c.b += sphereColor.b;
                }
            }

            // Update Closest Intersection Point
            closestPoint = intersection.z;
        }
    }

    return c;
}

/**
 * drawScene - Draws the Ray-Traced Scene
 */
void drawScene()
{
    // Initialize Aspect
    double aspect = (double)WIDTH/(double)HEIGHT;

    // Initialize Degree to Radian Coversion
    double degToRad = M_PI/180.0;

    // Initial Ray Direction Position
    double angle = tan((fov/2.0) * degToRad);

    // Start Drawing
    glPointSize(2.0);
    glBegin(GL_POINTS);

    // Iterate over Height of Image
    for(int y = 0; y < HEIGHT; y++)
    {
        // Iterate over Width of Image
        for(int x = 0; x < WIDTH; x++)
        {
            // Compute Delta X/Y
            double dx = (x + 0.5) / (double)WIDTH;
            double dy = (y + 0.5) / (double)HEIGHT;

            // Compute Screen Coordinates
            double xScreen = (2 * dx) - 1;
            double yScreen = (2 * dy) - 1;

            // Calculate Ray Direction
            VectorMath::point direction;
            direction.x = (xScreen * angle * aspect) - _cameraOrigin.x;
            direction.y = (yScreen * angle) - _cameraOrigin.y;
            direction.z = -1 - _cameraOrigin.z;

            // Create new Ray
            Ray ray = {_cameraOrigin, VectorMath::normalize(direction)};

            // Initialize Closest Point
            double closestPoint = MAX_DISTANCE;

            // Initialize Color from Ray
            RayColor c = {1.0, 1.0, 1.0};

            // Perform Sphere Intersection Calculations
            c = processSpheres(ray, c, closestPoint);

            // Perform Triangle Intersection Calculations
            c = processTriangles(ray, c, closestPoint);

            // Add Ambient Light
            c.r += _ambientLight[0];
            c.g += _ambientLight[1];
            c.b += _ambientLight[2];

            // Bound the Red Intensity
            if(c.r > 1.0)
            {
                c.r = 1.0;
            }

            // Bound the Green Intensity
            if(c.g > 1.0)
            {
                c.g = 1.0;
            }

            // Bound the Blue Intensity
            if(c.b > 1.0)
            {
                c.b = 1.0;
            }

            // Convert to Char Values
            double r = c.r * 255.0;
            double g = c.g * 255.0;
            double b = c.b * 255.0;

            // Set Pixel Color
            Color pixColor = {(char)r, (char)g, (char)b};

            // Draw Pixel
            plotPixel(x, y, pixColor.red, pixColor.green, pixColor.blue);
        }
    }

    // End Drawing
    glEnd();
    glFlush();

    // Log Debug
    std::cout << "Trace Completed \n" << std::endl;
}

/**
 * plotPixelDisplay
 */
void plotPixelDisplay(int x,int y,unsigned char r,unsigned char g,unsigned char b)
{
    glColor3f(((double)r)/256.f,((double)g)/256.f,((double)b)/256.f);
    glVertex2i(x,y);
}

/**
 * plotPixelJpeg
 */
void plotPixelJpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b)
{
    // Store the RGB values in the Buffer
    buffer[HEIGHT-y-1][x][0]=r;
    buffer[HEIGHT-y-1][x][1]=g;
    buffer[HEIGHT-y-1][x][2]=b;
}

/**
 * plotPixel
 */
void plotPixel(int x,int y,unsigned char r,unsigned char g, unsigned char b)
{
    plotPixelDisplay(x,y,r,g,b);

    // Check if JPEG Mode
    if(mode == MODE_JPEG)
    {
        plotPixelJpeg(x,y,r,g,b);
    }
}

/**
 * saveJpeg - Saves Image to JPEG File
 */
void saveJpeg()
{
    Pic *in = NULL;

    in = pic_alloc(640, 480, 3, NULL);
    printf("Saving JPEG file: %s\n", filename);

    memcpy(in->pix,buffer,3*WIDTH*HEIGHT);
    if (jpeg_write(filename, in))
    {
        printf("File saved Successfully\n");
    }
    else
    {
        printf("Error in Saving\n");
    }

    pic_free(in);
}

/**
 * parseCheck
 */
void parseCheck(char *expected,char *found)
{
    if(strcasecmp(expected,found))
    {
        char error[100];
        printf("Expected '%s ' found '%s '\n",expected,found);
        printf("Parse error, abnormal abortion\n");
        exit(0);
    }
}

/**
 * parseDoubles
 */
void parseDoubles(FILE*file, char *check, double p[3])
{
    char str[100];
    fscanf(file,"%s",str);
    parseCheck(check,str);
    fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
    printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

/**
 * parseRad
 */
void parseRad(FILE*file, double *r)
{
    char str[100];
    fscanf(file,"%s",str);
    parseCheck("rad:",str);
    fscanf(file,"%lf",r);
    printf("rad: %f\n",*r);
}

/**
 * parseShi - Read in Shininess Value from Image File
 */
void parseShi(FILE*file, double *shi)
{
    char s[100];
    fscanf(file,"%s",s);
    parseCheck("shi:",s);
    fscanf(file,"%lf",shi);
    printf("shi: %f\n",*shi);
}

/**
 * loadScene - Loads in Values from Image File
 */
int loadScene(char *argv)
{
    // Load in the Scene File
    FILE *file = fopen(argv,"r");

    // Create Number of Objects
    int numObjects;

    char type[50];
    int i;
    Triangle t;
    Sphere s;
    Light l;

    // Read in the Number of Objects
    fscanf(file,"%i",&numObjects);

    // Debug Stmt
    printf("number of objects: %i\n",numObjects);
    char str[200];

    parseDoubles(file,"amb:",_ambientLight);

    // Iterate over the Objects
    for(i=0; i<numObjects; i++)
    {
        fscanf(file,"%s\n",type);
        printf("%s\n",type);

        // Check if Object is a Triangle
        if(strcasecmp(type,"triangle") == 0)
        {
            // Log Triangle Object
            printf("found triangle\n");

            for(int j=0; j<3; j++)
            {
                parseDoubles(file,"pos:",t.v[j].position);
                parseDoubles(file,"nor:",t.v[j].normal);
                parseDoubles(file,"dif:",t.v[j].color_diffuse);
                parseDoubles(file,"spe:",t.v[j].color_specular);
                parseShi(file,&t.v[j].shininess);
            }

            // Check if Max Triangle Threshold reached
            if(_numTriangles == MAX_TRIANGLES)
            {
                printf("too many triangles, you should increase MAX_TRIANGLES!\n");
                exit(0);
            }
            _triangles[_numTriangles++] = t;
        }
        else if(strcasecmp(type,"sphere")==0)
        {
            printf("found sphere\n");

            parseDoubles(file,"pos:",s.position);
            parseRad(file,&s.radius);
            parseDoubles(file,"dif:",s.color_diffuse);
            parseDoubles(file,"spe:",s.color_specular);
            parseShi(file,&s.shininess);

            // Check if Max Sphere Threshold reached
            if(_numSpheres == MAX_SPHERES)
            {
                // Log Error, and exit program
                printf("too many spheres, you should increase MAX_SPHERES!\n");
                exit(0);
            }

            // Store Sphere
            _spheres[_numSpheres++] = s;
        }
        else if(strcasecmp(type,"light")==0)
        {
            printf("found light\n");
            parseDoubles(file,"pos:",l.position);
            parseDoubles(file,"col:",l.color);

            if(_numLights == MAX_LIGHTS)
            {
                printf("too many lights, you should increase MAX_LIGHTS!\n");
                exit(0);
            }

            // Store Light
            _lights[_numLights++] = l;
        }
        else
        {
            printf("unknown type in scene description:\n%s\n",type);
            exit(0);
        }
    }
    return 0;
}

/**
 * display
 */
void display()
{
    // N/A
}

/**
 * idle
 */
void idle()
{
    // Ensure drawn only once
    if(!once)
    {
        // Draw the Scene
        drawScene();

        // Check if JPEG Mode
        if(mode == MODE_JPEG)
        {
            // Save the Image
            saveJpeg();
        }
    }

    once=1;
}

/**
 * init
 */
void init()
{
    glMatrixMode(GL_PROJECTION);
    glOrtho(0,WIDTH,0,HEIGHT,1,-1);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glClearColor(0,0,0,0);
    glClear(GL_COLOR_BUFFER_BIT);
}

/**
 * main
 */
int main (int argc, char ** argv)
{
    // Ensure Program has only 2-3 Arguments
    if (argc < 2 || argc > 3)
    {
        printf ("usage: %s <scenefile> [jpegname]\n", argv[0]);
        exit(0);
    }

    // Check if Program has 3 Arguments
    if(argc == 3)
    {
        // Set as JPEG Mode
        mode = MODE_JPEG;

        // Save JPEG Filename
        filename = argv[2];
    }
    else if(argc == 2)
    {
        // Set as Display Mode
        mode = MODE_DISPLAY;
    }

    // Initialize GLUT
    glutInit(&argc,argv);

    // Load Scene
    loadScene(argv[1]);

    // Request Single Buffer, and Color
    glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);

    // Set window size/position
    glutInitWindowSize(WIDTH,HEIGHT);
    glutInitWindowPosition(0,0);

    // Create Window
    glutCreateWindow("Ray Tracer");

    // GLUT Callbacks
    glutDisplayFunc(display);
    glutIdleFunc(idle);

    // Initialize States
    init();

    glutMainLoop();
    return (0);
}
