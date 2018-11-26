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

#define MAX_TRIANGLES 2000
#define MAX_SPHERES 10
#define MAX_LIGHTS 10

#define DEG_TO_RAD 3.14159265/180

char *filename=0;

//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2
int mode=MODE_DISPLAY;

//you may want to make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480

//the field of view of the camera
#define fov 60.0

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
Triangle triangles[MAX_TRIANGLES];

// Initialize Sphere Data Struct
Sphere spheres[MAX_SPHERES];

// Initialize Lights Data Struct
Light lights[MAX_LIGHTS];

double ambient_light[3];

// Initialize Amount Indicators
int _numTriangles = 0;
int _numSpheres = 0;
int _numLights = 0;

int _sphereIntersections;

void plotPixelDisplay(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plotPixelJpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plotPixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);

// Initialize Camera Origin Point for each Rays
VectorMath::point _cameraOrigin = {0.0, 0.0, 0.0};

/**
 * intersectSphere
 */
double intersectSphere(Ray &ray, VectorMath::point center, double radius)
{
    // Create Intersection Values
    double t0 = -1.0;
    double t1 = -1.0;

    VectorMath::point l = VectorMath::subtractVectors(ray.getOrigin(), center);

    // Obtain Scalar Values from Quadratic Equation
    double a = VectorMath::dotProduct(ray.getDirection(), ray.getDirection());
    double b = 2 * (VectorMath::dotProduct(ray.getDirection(), l));
    double c = VectorMath::dotProduct(l, l) - pow(radius, 2);

    // Check if valid intersection points
    if(VectorMath::quadratic(a, b, c, t0, t1))
    {
        // Ensure Ascending Order
        if(t0 > t1)
        {
            // Update Intersection Values
            std::swap(t0,t1);
        }

        // Check if t0 is negative
        if (t0 < 0)
        {
            // Use t1 instead
            t0 = t1;

            // Check if t1 is negative
            if (t0 < 0)
            {
                // Set t0 to be negative
                t0 = -1.0;
            }
       }
    }

    return t0;
}

/**
 * computeSphereIllumination - Calculates the Colors to display at the Sphere
 *                             based on the Phong Illumination Model
 *
 * param s          - The Sphere
 * param interPoint - The Intersection Point of the Viewing Vector with the Sphere
 */
Color computeSphereIllumination(Sphere s, double interPoint)
{
    Color c = {0.0, 255.0, 0.0};

    // Get the Diffuse for each Color Channel
    double kDRed = s.color_diffuse[0];
    double kDGreen = s.color_diffuse[1];
    double kDBlue = s.color_diffuse[2];

    // Get the Specular for each Color Channel
    double kSRed = s.color_specular[0];
    double kSGreen = s.color_specular[1];
    double kSBlue = s.color_specular[2];

    // Calculate V Vector (Vector from Intersection Point to Camera)
    //VectorMath::point vVector = VectorMath::subtractVectors(v1, v2);

    // Iterate over the Lights in the Scene
    for(int i = 0; i < _numLights; ++i)
    {
        // Calculate L Vector (Vector from Intersection Point to Light
        // TODO

        // Calculate the Red Intensity
        // TODO

        // Calculate the Green Intensity
        // TODO

        // Calculate the Blue Intensity
        // TODO
    }

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
        glPointSize(2.0);
        glBegin(GL_POINTS);

        // Iterate over Height of Image
        for(int y = 0; y < HEIGHT; y++)
        {
            // Draw Pixel
            plotPixel(x, y, 255.0, 255.0, 255.0);
        }
        glEnd();
        glFlush();
    }
}

/**
 * drawScene
 */
void drawScene()
{
    // Draw Background
    drawBackground();

    // Initialize Aspect
    float aspect = WIDTH/HEIGHT;

    // Initial Ray Direction Position
    float initX = -aspect * tan((fov/2) * DEG_TO_RAD);
    float initY = tan((fov/2) * DEG_TO_RAD);
    float initZ = -1;

    // Iterate over Height of Image
    for(int y = 0; y < HEIGHT; y++)
    {
        // Iterate over Width of Image
        for(int x = 0; x < WIDTH; x++)
        {
            // Calculate Delta X/Y
            float xLength = -initX - initX;
            float yLength = -initY - initY;
            float dx = (x * xLength) / WIDTH;
            float dy = (y * yLength) / HEIGHT;

            // Calculate Ray Direction
            VectorMath::point direction;
            direction.x = initX + dx - _cameraOrigin.x;
            direction.y = initY + dy - _cameraOrigin.y;
            direction.z = initZ - _cameraOrigin.z;

            // Create new Ray
            Ray *ray = new Ray(_cameraOrigin, direction);

            // Iterate over all the Spheres
            for(int i = 0; i < _numSpheres; ++i)
            {
                // Get the Sphere
                Sphere s = spheres[i];

                // Get Center Point of Sphere
                VectorMath::point center = {s.position[0], s.position[1], s.position[2]};

                // Get the Intersection Point
                double interPoint = intersectSphere(*ray, center, s.radius);

                // Check if Ray Intersects Sphere
                if(interPoint > 0)
                {
                    // Calculate Phong Illumination for the Sphere
                    Color c = computeSphereIllumination(s, interPoint);

                    // Draw Pixel
                    plotPixel(x, y, c.red, c.green, c.blue);
                }
            }

            // Iterate over all the Triangles
            for(int j = 0; j < _numTriangles; ++j)
            {
                // Get the Triangle
                Triangle t = triangles[j];
            }
        }
    }

    // Log Debug
    printf("Done!\n");
    fflush(stdout);
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

    parseDoubles(file,"amb:",ambient_light);

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
            triangles[_numTriangles++] = t;
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
            spheres[_numSpheres++] = s;
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
            lights[_numLights++] = l;
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
 * intersectTriangle
 */
bool intersectTriangle(Ray &ray)
{

}


/**
 * idle
 */
void idle()
{
    //hack to make it only draw once
    static int once = 0;

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

        glFinish();
    }

    // Update Draw Once Indicator
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
    glClearColor(0.0, 0.0, 0.0, 0.0); // set background color
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
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGBA);

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
