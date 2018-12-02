/*
 * CSCI 420
 * Assignment 3 Raytracer
 *
 * Name: Matt Robinson
 */

#include <stdlib.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#include <pic.h>
#include "VectorMath.h"
#include "Ray.h"
#include <string.h>
#include <vector>
#include <iostream>
#include <unistd.h>
#include <sstream>

#define MAX_TRIANGLES 2000
#define MAX_SPHERES 10
#define MAX_LIGHTS 10

#define epsilon 0.000001

char *filename=0;

//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2
int mode=MODE_DISPLAY;

//you may want to make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480

//the field of view of the camera
#define fov 70.0

unsigned char buffer[HEIGHT][WIDTH][3];

struct Vertex
{
    double position[3];
    double color_diffuse[3];
    double color_specular[3];
    double normal[3];
    double shininess;
};

struct Color
{
    char red;
    char green;
    char blue;
};

typedef struct _Triangle
{
    struct Vertex v[3];
} Triangle;

typedef struct _Sphere
{
    double position[3];
    double color_diffuse[3];
    double color_specular[3];
    double shininess;
    double radius;
} Sphere;

typedef struct _Light
{
    double position[3];
    double color[3];
} Light;

// Initialize Triangle Data Struct
Triangle _triangles[MAX_TRIANGLES];

// Initialize Sphere Data Struct
Sphere _spheres[MAX_SPHERES];

// Initialize Lights Data Struct
Light _lights[MAX_LIGHTS];

double _alpha;
double _beta;
double _gamma;

// Ambient Light
double _ambientLight[3];

// Number Indicators
int _numTriangles = 0;
int _numSpheres = 0;
int _numLights = 0;

//hack to make it only draw once
static int once = 0;

int _sphereIntersections;

// Initialize Camera Origin Point for each Rays
VectorMath::point _cameraOrigin = {0.0, 0.0, 0.0};

void plotPixelDisplay(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plotPixelJpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plotPixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);

/**
 * calcArea - Calculates the Area of a Triangle in 3-D
 */
double calcArea(VectorMath::point a, VectorMath::point b, VectorMath::point c)
{
    // Calculate the Area
    double area = 0.5 * ((b.x-a.x)*(c.y-a.y) - (c.x-a.x)*(b.y-a.y));

    return area;
}

/**
 * pointInTriangle - Checks if intersection point is in triangle
 */
bool pointInTriangle(Triangle triangle, VectorMath::point p, VectorMath::point n)
{
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
    else
    {
        // Convert Vertices into Vector Points
        VectorMath::point p1 = {a.position[0], a.position[1], a.position[2]};
        VectorMath::point p2 = {b.position[0], b.position[1], b.position[2]};
        VectorMath::point p3 = {c.position[0], c.position[1], c.position[2]};

        // Calculate Area of Triangle
        double areaAlpha = calcArea(p, p2, p3);
        double areaBeta = calcArea(p1, p, p3);
        double area = calcArea(p1, p2, p3);

        // Set Barycentric Coordinates for Triangle
        _alpha = areaAlpha / area;
        _beta = areaBeta / area;
        _gamma = 1.0 - _alpha - _beta;
    }

    return isInTriangle;
}

/**
 * intersectSphere
 */
double intersectSphere(Ray &ray, VectorMath::point center, double radius)
{
    // Create Intersection Values
    double t0 = -1.0;
    double t1 = -1.0;
    double t = -1.0;

    VectorMath::point l = VectorMath::subtractVectors(ray.getOrigin(), center);

    // Obtain Scalar Values from Quadratic Equation
    double a = VectorMath::dotProduct(ray.getDirection(), ray.getDirection());
    double b = 2 * (VectorMath::dotProduct(ray.getDirection(), l));
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
            // Set to invalid value
            t = -1.0;
        }
    }

    return t;
}

/**
 * intersectTriangle
 */
double intersectTriangle(Triangle triangle, Ray &ray)
{
    // Initialize Intersection Value
    double t = -1.0;

    // Get Origin Vector
    VectorMath::point origin = ray.getOrigin();

    // Get Direction Vector
    VectorMath::point direction = ray.getDirection();

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
        t = -1.0;
    }
    else
    {
        // Step 3: Calculate d Vector
        VectorMath::point v0 = {a.position[0], a.position[1], a.position[2]};
        double d = VectorMath::dotProduct(n, v0);

        // Step 4: Calculate Intersection Point t [t = -(n dot p0 + d) / (n dot d)]
        t = ((VectorMath::dotProduct(n, origin) + d)) / nDotD;

        // Step 5: Check if Intersection is behind Ray Origin
        if(t <= 0)
        {
            // Set to invalid value
            t = -1.0;
        }
        else
        {
            // Compute the Intersection Point [P = origin + (t * Direction)]
            VectorMath::point p = VectorMath::addVectors(origin, VectorMath::scalarMultiply(t, direction));

            // Check if point not in triangle
            if(!pointInTriangle(triangle, p, n))
            {
                // Set to invalid value
                t = -1.0;
            }
        }
    }

    return t;
}

/**
 * computeDiffuse
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
 * computeSpecular
 */
double computeSpecular(double kS, double alpha, VectorMath::point rVec, VectorMath::point vVec)
{
    // Calculate [r dot v]
    double rDotV = VectorMath::dotProduct(rVec, vVec);

    if(rDotV < epsilon)
    {
        rDotV = 0;
    }

    // Calculate the Specular Component
    double specular = kS * pow(rDotV, alpha);

    return specular;
}

/**
 * computeRVec
 */
VectorMath::point computeRVec(VectorMath::point lVec, VectorMath::point nVec)
{
    // Calculate [l dot N]
    double lDotN = VectorMath::dotProduct(lVec, nVec);

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
 * param interPoint - The Intersection Point of the Viewing Vector with the Sphere
 */
Color computeSphereIllumination(Sphere s, VectorMath::point interPoint)
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

    // Calculate N Vector (Unit Normal)
    VectorMath::point center = {s.position[0], s.position[1], s.position[2]};
    VectorMath::point sphereVector = VectorMath::subtractVectors(interPoint, center);
    double length = VectorMath::magnitude(sphereVector);

    // Calculate N Vector
    VectorMath::point nVec = VectorMath::scalarDivision(s.radius, VectorMath::subtractVectors(interPoint, center));

    // Initialize Light Intensity
    double iRed = 0.0;
    double iGreen = 0.0;
    double iBlue = 0.0;

    // Iterate over the Lights in the Scene
    for(int i = 0; i < _numLights; ++i)
    {
        // Get Light Position
        VectorMath::point light = {_lights[i].position[0],
                _lights[i].position[1],
                _lights[i].position[2]};

        // Calculate L Vector (Vector from Intersection Point to Light)
        VectorMath::point lVec = VectorMath::normalize(VectorMath::subtractVectors(light, interPoint));

        // Calculate R Vector
        VectorMath::point rVec = computeRVec(lVec, nVec);

        // Calculate the Red Intensity
        double diffuseRed = computeDiffuse(kDRed, lVec, nVec);
        double specularRed = computeSpecular(kSRed, alpha, rVec, vVec);
        iRed += _lights[i].color[0] * (diffuseRed + specularRed);

        // Calculate the Green Intensity
        double diffuseGreen = computeDiffuse(kDGreen, lVec, nVec);
        double specularGreen = computeSpecular(kSGreen, alpha, rVec, vVec);
        iGreen += _lights[i].color[1] * (diffuseGreen + specularGreen);

        // Calculate the Blue Intensity
        double diffuseBlue = computeDiffuse(kDBlue, lVec, nVec);
        double specularBlue = computeSpecular(kSBlue, alpha, rVec, vVec);
        iBlue += _lights[i].color[2] * (diffuseBlue + specularBlue);
    }

    // Avg Intensity per Light
    iRed /= _numLights;
    iGreen /= _numLights;
    iBlue /= _numLights;

    // Add Ambient Light
    iRed += _ambientLight[0];
    iGreen += _ambientLight[1];
    iBlue += _ambientLight[2];

    // Bound the Red Intensity
    if(iRed >= 1.0)
    {
        iRed = 1.0;
    }

    // Bound the Green Intensity
    if(iGreen >= 1.0)
    {
        iGreen = 1.0;
    }

    // Bound the Blue Intensity
    if(iBlue >= 1.0)
    {
        iBlue = 1.0;
    }

    // Convert to Char Values
    double r = iRed * 255.0;
    double g = iGreen * 255.0;
    double b = iBlue * 255.0;

    // Set Color
    Color c = {(char)r, (char)g, (char)b};

    return c;
}

/**
 * computeTriangleIllumination - Calculates the Colors to display at the Triangle
 *                                based on the Phong Illumination Model
 *
 * param t          - The Triangle
 * param interPoint - The Intersection Point of the Viewing Vector with the Triangle
 */
Color computeTriangleIllumination(Triangle t, VectorMath::point interPoint)
{
    // Initialize Light Intensity
    double iRed, iGreen, iBlue;

    // Calculate V Vector (Vector from Intersection Point to Camera)
    VectorMath::point vVec = VectorMath::subtractVectors(_cameraOrigin, interPoint);
    vVec = VectorMath::normalize(vVec);

    // Calculate N Vector
    VectorMath::point nVec;
    nVec.x = (_alpha * t.v[0].normal[0]) + (_beta * t.v[1].normal[0]) + (_gamma * t.v[2].normal[0]);
    nVec.y = (_alpha * t.v[0].normal[1]) + (_beta * t.v[1].normal[1]) + (_gamma * t.v[2].normal[1]);
    nVec.z = (_alpha * t.v[0].normal[2]) + (_beta * t.v[1].normal[2]) + (_gamma * t.v[2].normal[2]);
    nVec = VectorMath::normalize(nVec);

    // Calculate Diffuse
    double kDRed =   (_alpha * t.v[0].color_diffuse[0]) + (_beta * t.v[1].color_diffuse[0]) + (_gamma * t.v[2].color_diffuse[0]);
    double kDGreen = (_alpha * t.v[0].color_diffuse[1]) + (_beta * t.v[1].color_diffuse[1]) + (_gamma * t.v[2].color_diffuse[1]);
    double kDBlue =  (_alpha * t.v[0].color_diffuse[2]) + (_beta * t.v[1].color_diffuse[2]) + (_gamma * t.v[2].color_diffuse[2]);

    // Calculate Specular
    double kSRed =   (_alpha * t.v[0].color_specular[0]) + (_beta * t.v[1].color_specular[0]) + (_gamma * t.v[2].color_specular[0]);
    double kSGreen = (_alpha * t.v[0].color_specular[1]) + (_beta * t.v[1].color_specular[1]) + (_gamma * t.v[2].color_specular[1]);
    double kSBlue =  (_alpha * t.v[0].color_specular[2]) + (_beta * t.v[1].color_specular[2]) + (_gamma * t.v[2].color_specular[2]);

    // Calculate Shininess Value
    double shininess = (_alpha * t.v[0].shininess) + (_beta * t.v[1].shininess) + (_gamma * t.v[2].shininess);

    // Iterate over the Lights in the Scene
    for(int i = 0; i < _numLights; ++i)
    {
        // Get Light Position
        VectorMath::point light = {_lights[i].position[0],
                                   _lights[i].position[1],
                                   _lights[i].position[2]};

        // Calculate L Vector (Vector from Intersection Point to Light)
        VectorMath::point lVec = VectorMath::normalize(VectorMath::subtractVectors(light, interPoint));

        // Calculate R Vector
        VectorMath::point rVec = computeRVec(lVec, nVec);

        // Calculate the Red Intensity
        double diffuseRed = computeDiffuse(kDRed, lVec, nVec);
        double specularRed = computeSpecular(kSRed, shininess, rVec, vVec);
        iRed += _lights[i].color[0] * (diffuseRed + specularRed);

        // Calculate the Green Intensity
        double diffuseGreen = computeDiffuse(kDGreen, lVec, nVec);
        double specularGreen = computeSpecular(kSGreen, shininess, rVec, vVec);
        iGreen += _lights[i].color[1] * (diffuseGreen + specularGreen);

        // Calculate the Blue Intensity
        double diffuseBlue = computeDiffuse(kDBlue, lVec, nVec);
        double specularBlue = computeSpecular(kSBlue, shininess, rVec, vVec);
        iBlue += _lights[i].color[2] * (diffuseBlue + specularBlue);
    }

    // Avg Intensity per Light
    iRed /= _numLights;
    iGreen /= _numLights;
    iBlue /= _numLights;

    // Add Ambient Light
    iRed += _ambientLight[0];
    iGreen += _ambientLight[1];
    iBlue += _ambientLight[2];

    // Bound the Red Intensity
    if(iRed > 1.0)
    {
        iRed = 1.0;
    }

    // Bound the Green Intensity
    if(iGreen > 1.0)
    {
        iGreen = 1.0;
    }

    // Bound the Blue Intensity
    if(iBlue > 1.0)
    {
        iBlue = 1.0;
    }

    // Convert to Char Values
    double r = iRed * 255.0;
    double g = iGreen * 255.0;
    double b = iBlue * 255.0;

    // Set Color
    Color c = {(char)r, (char)g, (char)b};

    return c;
}

/**
 * drawBackground - Draws the Background of the Image
 */
void drawBackground()
{
    // Iterate over Width of Image
    for(int x = 0; x < WIDTH; x++)
    {
        // Iterate over Height of Image
        for(int y = 0; y < HEIGHT; y++)
        {
            // Draw Pixel
            plotPixel(x, y, 255.0, 255.0, 255.0);
        }
    }
}

/**
 * drawScene
 */
void drawScene()
{
    glPointSize(2.0);
    glBegin(GL_POINTS);

    // Draw Background
    drawBackground();

    // Initialize Aspect
    float aspect = WIDTH/HEIGHT;

    // Initialize Degree to Radian Coversion
    float degToRad = M_PI/180.0;

    // Initial Ray Direction Position
    float initX = -aspect * tan((fov/2) * degToRad);
    float initY = tan((fov/2) * degToRad);
    float initZ = -1;

    // Iterate over Height of Image
    for(int y = 0; y < HEIGHT; y++)
    {
        // Iterate over Width of Image
        for(int x = 0; x < WIDTH; x++)
        {
            // Calculate Delta X/Y
            float xLength = (-initX) - initX;
            float yLength = (-initY) - initY;
            float dx = (x * xLength) / WIDTH;
            float dy = (y * yLength) / HEIGHT;

            // Calculate Ray Direction
            VectorMath::point direction;
            direction.x = initX + dx - _cameraOrigin.x;
            direction.y = initY + dy - _cameraOrigin.y;
            direction.z = initZ - _cameraOrigin.z;

            // Create new Ray
            Ray *ray = new Ray(_cameraOrigin, direction);

            // Initialize Closest Point
            double closestPoint = 1000000;

            // Initialize Color
            Color c = {0.0, 0.0, 0.0};

            // Iterate over all the Triangles
            for(int j = 0; j < _numTriangles; ++j)
            {
                // Get the Triangle
                Triangle triangle = _triangles[j];

                // Check if there is an intersection with the Triangle
                double tTriangle = intersectTriangle(triangle, *ray);

                // Check if Ray Intersects Triangle
                if(tTriangle > 0)
                {
                    // Calculate Intersection Point [p = v0 + (vd*t)]
                    VectorMath::point triangleInterPoint = VectorMath::addVectors(ray->getOrigin(), VectorMath::scalarMultiply(tTriangle, ray->getDirection()));

                    // Check if Closest Point to Camera (Z Value)
                    if(triangleInterPoint.z < closestPoint)
                    {
                        // Update Closest Point
                        closestPoint = triangleInterPoint.z;

                        // Calculate Phong Illumination for the Sphere
                        c = computeTriangleIllumination(triangle, triangleInterPoint);

                        // Draw Pixel
                        plotPixel(x, HEIGHT - y, c.red, c.green, c.blue);
                    }
                }
            }

            // Iterate over all the Spheres
            for(int i = 0; i < _numSpheres; ++i)
            {
                // Get the Sphere
                Sphere s = _spheres[i];

                // Get Center Point of Sphere
                VectorMath::point center = {s.position[0], s.position[1], s.position[2]};

                // Check if there is an intersection with the Sphere
                double tSphere = intersectSphere(*ray, center, s.radius);

                // Check if Ray Intersects Sphere
                if(tSphere > 0)
                {
                    // Calculate Intersection Point [p = v0 + (vd*t)]
                    VectorMath::point sphereInterPoint = VectorMath::addVectors(ray->getOrigin(), VectorMath::scalarMultiply(tSphere, ray->getDirection()));

                    // Check if Closest Point to Camera (Z Value)
                    //if(sphereInterPoint.z < closestPoint)
                    //{
                        // Update Closest Point
                        //closestPoint = sphereInterPoint.z;

                        // Calculate Phong Illumination for the Sphere
                        Color c = computeSphereIllumination(s, sphereInterPoint);

                        // Draw Pixel
                        plotPixel(x, HEIGHT-y, c.red, c.green, c.blue);
                    //}
                }
            }
        }
    }

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
 * saveJpeg
 */
void saveJpeg()
{
    Pic *in = NULL;

    in = pic_alloc(640, 480, 3, NULL);
    printf("Saving JPEG file: %s\n", filename);

    memcpy(in->pix,buffer,3*WIDTH*HEIGHT);
    if (jpeg_write(filename, in))
        printf("File saved Successfully\n");
    else
        printf("Error in Saving\n");

    pic_free(in);

}

void parse_check(char *expected,char *found)
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
    parse_check(check,str);
    fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
    printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE*file,double *r)
{
    char str[100];
    fscanf(file,"%s",str);
    parse_check("rad:",str);
    fscanf(file,"%lf",r);
    printf("rad: %f\n",*r);
}

void parse_shi(FILE*file,double *shi)
{
    char s[100];
    fscanf(file,"%s",s);
    parse_check("shi:",s);
    fscanf(file,"%lf",shi);
    printf("shi: %f\n",*shi);
}

/**
 * loadScene
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
                parse_shi(file,&t.v[j].shininess);
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
            parse_rad(file,&s.radius);
            parseDoubles(file,"dif:",s.color_diffuse);
            parseDoubles(file,"spe:",s.color_specular);
            parse_shi(file,&s.shininess);

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
    return 0;
}
