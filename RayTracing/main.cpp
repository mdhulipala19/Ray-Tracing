/*
Mahesh Dhulipala HW3
Finished features
1. Beach
2. Beach ball
3. Parasol
4. Sand castle
5. Ocean (requires Beach)
6. Waves (requires Ocean)
7. Shadows (requires Beach)
8. Palm tree
9. Pinking (requires Palm tree)
10. Flotsam
11. This side up (requires Flotsam)
12. Disco
13. Sky
14. Animation
*/

#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#if defined(__APPLE__)
#include <GLUT/GLUT.h>
#include <OpenGL/gl3.h>
#include <OpenGL/glu.h>
#else
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__)
#include <windows.h>
#endif
#include <GL/glew.h>
#include <GL/freeglut.h>
#endif

#include <vector>


#include "vec2.h"
#include "vec3.h"
#include "vec4.h"
#include "matrix4x4.h"
#include <limits>
#include <iostream>
#include <algorithm>

const unsigned int windowWidth = 512, windowHeight = 512;
int majorVersion = 3, minorVersion = 0;

void getErrorInfo(unsigned int handle)
{
    int logLen;
    glGetShaderiv(handle, GL_INFO_LOG_LENGTH, &logLen);
    if (logLen > 0)
    {
        char * log = new char[logLen];
        int written;
        glGetShaderInfoLog(handle, logLen, &written, log);
        printf("Shader log:\n%s", log);
        delete log;
    }
}

void checkShader(unsigned int shader, char * message)
{
    int OK;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &OK);
    if (!OK)
    {
        printf("%s!\n", message);
        getErrorInfo(shader);
    }
}

void checkLinking(unsigned int program)
{
    int OK;
    glGetProgramiv(program, GL_LINK_STATUS, &OK);
    if (!OK)
    {
        printf("Failed to link shader program!\n");
        getErrorInfo(program);
    }
}

class Shader
{
    protected:
    unsigned int shaderProgram;
    
    public:
    Shader()
    {
        const char *vertexSource = "\n\
        #version 410 \n\
        precision highp float; \n\
        \n\
        in vec2 vertexPosition;    \n\
        in vec2 vertexTexCoord; \n\
        out vec2 texCoord; \n\
        \n\
        void main() \n\
        { \n\
        texCoord = vertexTexCoord; \n\
        gl_Position = vec4(vertexPosition.x, vertexPosition.y, 0, 1); \n\
        } \n\
        ";
        
        const char *fragmentSource = "\n\
        #version 410 \n\
        precision highp float; \n\
        \n\
        uniform sampler2D samplerUnit; \n\
        in vec2 texCoord;  \n\
        out vec4 fragmentColor; \n\
        \n\
        void main() { \n\
        fragmentColor = texture(samplerUnit, texCoord);  \n\
        } \n\
        ";
        
        unsigned int vertexShader = glCreateShader(GL_VERTEX_SHADER);
        if (!vertexShader) { printf("Error in vertex shader creation\n"); exit(1); }
        
        glShaderSource(vertexShader, 1, &vertexSource, NULL);
        glCompileShader(vertexShader);
        checkShader(vertexShader, "Vertex shader error");
        
        unsigned int fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
        if (!fragmentShader) { printf("Error in fragment shader creation\n"); exit(1); }
        
        glShaderSource(fragmentShader, 1, &fragmentSource, NULL);
        glCompileShader(fragmentShader);
        checkShader(fragmentShader, "Fragment shader error");
        
        shaderProgram = glCreateProgram();
        if (!shaderProgram) { printf("Error in shader program creation\n"); exit(1); }
        
        glAttachShader(shaderProgram, vertexShader);
        glAttachShader(shaderProgram, fragmentShader);
        
        glBindAttribLocation(shaderProgram, 0, "vertexPosition");
        glBindAttribLocation(shaderProgram, 1, "vertexTexCoord");
        
        glBindFragDataLocation(shaderProgram, 0, "fragmentColor");
        
        glLinkProgram(shaderProgram);
        checkLinking(shaderProgram);
    }
    
    ~Shader()
    {
        if(shaderProgram) glDeleteProgram(shaderProgram);
    }
    
    void Run()
    {
        if(shaderProgram) glUseProgram(shaderProgram);
    }
    
    void UploadSamplerID()
    {
        int samplerUnit = 0;
        int location = glGetUniformLocation(shaderProgram, "samplerUnit");
        glUniform1i(location, samplerUnit);
        glActiveTexture(GL_TEXTURE0 + samplerUnit);
    }
};

Shader *shader = 0;

// Simple material class, with object color, and headlight shading.

class Material
{
    public:
    vec3 color;
    vec3 kd=vec3(1,1,1);
    vec3 ks=vec3(1,1,1);
    vec3 ka=vec3(1,1,1);
    bool isMetal=false;
    
    public:
    Material(vec3 color2=vec3(1,0,0))
    {
        color=color2;
    }
    
    public:
    
    virtual vec3 getColor(
                          vec3 position,
                          vec3 normal,
                          vec3 viewDir)
    {
        return color*normal.dot(viewDir);
        
    }
    virtual vec3 shade(vec3 position,
                       vec3 normal,
                       vec3 viewDir,vec3 lightDir,vec3 powerDensity)=0;
    //    virtual vec3 shade(vec3 position,
    //                       vec3 normal,
    //                       vec3 viewDir,vec3 lightDir,vec3 powerDensity)
    //    {
    //        std::cout<<normal.norm()<<std::endl;
    //        vec3 answer=powerDensity*(kd*(abs(normal.dot(lightDir))));
    //        return answer;
    //    }
};
float snoise(vec3 r) {
    unsigned int x = 0x0625DF73;
    unsigned int y = 0xD1B84B45;
    unsigned int z = 0x152AD8D0;
    float f = 0;
    for(int i=0; i<32; i++) {
        vec3 s(    x/(float)0xffffffff,
               y/(float)0xffffffff,
               z/(float)0xffffffff);
        f += sin(s.dot(r));
        x = x << 1 | x >> 31;
        y = y << 1 | y >> 31;
        z = z << 1 | z >> 31;
    }
    return f / 64.0 + 0.5;
};
vec3 snoise2(vec3 r) {
    unsigned int x = 0x0625DF73;
    unsigned int y = 0xD1B84B45;
    unsigned int z = 0x152AD8D0;
    float f = 0;
    for(int i=0; i<32; i++) {
        vec3 s(    x/(float)0xffffffff,
               y/(float)0xffffffff,
               z/(float)0xffffffff);
        f += sin(s.dot(r));
        x = x << 1 | x >> 31;
        y = y << 1 | y >> 31;
        z = z << 1 | z >> 31;
    }
    return vec3(1,2,3)*(f / 64.0 + 0.5);
};


class Wood : public Material
{
    float scale;
    float turbulence;
    float period;
    float sharpness;
    public:
    Wood():
    Material(vec3(1, 1, 1))
    {
        scale = 10;
        turbulence = 500;
        period = 10;
        sharpness = 1;
    }
     vec3 getColor(
                          vec3 position,
                          vec3 normal,
                          vec3 viewDir)
    {
        //return normal;
        float w = position.x * period + pow(snoise(position * scale), sharpness)*turbulence + 10000.0;
        w -= int(w);
        return (vec3(1, 0.3, 0) * w + vec3(0.35, 0.1, 0.05) * (1-w)) * normal.dot(viewDir);
    }
    virtual vec3 shade(vec3 position,
                       vec3 normal,
                       vec3 viewDir,vec3 lightDir,vec3 powerDensity)
    {
        return getColor(position,normal,viewDir);
    }
};
class Ball : public Material
{
    float scale;
    float turbulence;
    float period;
    float sharpness;
public:
    Ball():
    Material(vec3(1, 1, 1))
    {
        scale = 0.01;
        turbulence = 500;
        period = 10;
        sharpness = 1;
    }
    virtual vec3 getColor(
                          vec3 position,
                          vec3 normal,
                          vec3 viewDir)
    {
        vec3 color3;
        float angle=atan2(normal.x, normal.z);
        float value=abs((fmod((int)(angle*180/M_PI),(int)(180))));
        if(value<22.5||(45<value&&value<67.5))
        {
            kd=vec3(1,0,0);
        }
        else
        {
            kd=vec3(1,1,1);
        }
        
        
        return kd;
    }
    virtual vec3 shade(vec3 position,
                       vec3 normal,
                       vec3 viewDir,vec3 lightDir,vec3 powerDensity)
    {
        float gamma=2;
        kd=getColor(position,normal,viewDir);
        ks=vec3(0.5,0.5,0.5);
        vec3 answer=powerDensity*(kd*(normal.dot(lightDir)))+powerDensity*((ks)*(pow(normal.dot((lightDir+viewDir).normalize()),gamma)));
        return answer;
    }
    
   
    
    
    
    
};
class Sand : public Material
{
    float scale;
    float turbulence;
    float period;
    float sharpness;
public:
    Sand():
    Material(vec3(1, 1, 1))
    {
        scale = 100;
        turbulence = 1;
        period = 5;
        sharpness = 10;
    }
     vec3 getColor(
                          vec3 position,
                          vec3 normal,
                          vec3 viewDir)
    {
        //return normal;
        float w = position.x * period + pow(snoise(position * scale), sharpness)*turbulence + 10000.0;
        w -= int(w);
        return (vec3(1, 0.3, 0) * w + vec3(0.35, 0.1, 0.05) * (1-w)) * normal.dot(viewDir);
    }
     vec3 shade(vec3 position,
                       vec3 normal,
                       vec3 viewDir,vec3 lightDir,vec3 powerDensity)
    {
        return getColor(position,normal,viewDir);
    }
};
class Marble : public Material
{
    float scale;
    float turbulence;
    float period;
    float sharpness;
    public:
    Marble():
    Material(vec3(1, 1, 1))
    {
        scale = 32;
        turbulence = 50;
        period = 32;
        sharpness = 1;
    }
     vec3 getColor(
                          vec3 position,
                          vec3 normal,
                          vec3 viewDir)
    {
        //return normal;
        float w = position.x * period + pow(snoise(position * scale), sharpness)*turbulence;
        w = pow(sin(w)*0.5+0.5, 4);
        return (vec3(0, 0, 1) * w + vec3(1, 1, 1) * (1-w)) * normal.dot(viewDir);
    }
     vec3 shade(vec3 position,
                       vec3 normal,
                       vec3 viewDir,vec3 lightDir,vec3 powerDensity)
    {
        return getColor(position,normal,viewDir);
    }
    
    
};

class DiffuseMaterial :public Material
{
    
    
    public:
    DiffuseMaterial(vec3 c=vec3(1,1,1))
    {
        //color=c;
        kd=c;
    }
     vec3 shade(vec3 position,
                       vec3 normal,
                       vec3 viewDir,vec3 lightDir,vec3 powerDensity)
    {
        ks=vec3(0,0,0);
        float gamma=1;
        vec3 answer=powerDensity*(kd*(normal.dot(lightDir)))+powerDensity*((ks)*(pow(normal.dot((lightDir+viewDir).normalize()),gamma)));
        return answer;
    }
};



class SpecularMaterial :public Material
{
    
    
    public:
    SpecularMaterial(vec3 c=vec3(1,0,1))
    {
        color=c;
        kd=c;
    }
     vec3 shade(vec3 position,
                       vec3 normal,
                       vec3 viewDir,vec3 lightDir,vec3 powerDensity)
    {
        
        float gamma=2;
        ks=vec3(1,1,1);
        vec3 answer=powerDensity*(kd*(normal.dot(lightDir)))+powerDensity*((ks)*(pow(normal.dot((lightDir+viewDir).normalize()),gamma)));
        return answer;
    }
};
class Metal:public Material
{
public:
    Metal(vec3 c=vec3(1,1,1))
    {
        color=c;
        kd=c;
    }
     virtual vec3 shade(vec3 position,
                       vec3 normal,
                       vec3 viewDir,vec3 lightDir,vec3 powerDensity)
    {
        float gamma=10;
        kd=vec3(0.0,0,1);
        ks=vec3(0,0,0);
        vec3 answer=powerDensity*(kd*(normal.dot(lightDir)))+powerDensity*((ks)*(pow(normal.dot(((lightDir+viewDir).normalize())),gamma)));
        return answer;
    }
    
};

class Water:public Metal
{
public:
    Water(vec3 c=vec3(1,1,1))
    {
        color=c;
        kd=c;
    }
    vec3 noiseGrad(vec3 r) {
        vec3 s = vec3(7502, 22777, 4767);
        vec3 f = vec3(0.0, 0.0, 0.0);
        for(int i=0; i<16; i++) {
            
            f +=  (s - vec3(32768, 32768, 32768)) * 40.0*cos((s - vec3(32768, 32768, 32768)).dot(r*40.0)/ 65536.0);
           
            
            s = s%32768.0 * 2.0 + vec3(floor(fmod(s.x,32768.0)),floor(fmod(s.x,32768.0)),floor(fmod(s.x,32768.0)));
        }
        return vec3(f.x / 65536.0,f.y / 65536.0,f.z / 65536.0);
    }
    
    vec3 shade(vec3 position,
               vec3 normal,
               vec3 viewDir,vec3 lightDir,vec3 powerDensity)
    {
       
        float gamma=10;
        kd=vec3(0.0,0,1);
        ks=vec3(0,0,0);
        normal=(normal+noiseGrad(position)).normalize();
        vec3 answer=powerDensity*(kd*(normal.dot(lightDir)))+powerDensity*((ks)*(pow(normal.dot(((lightDir+viewDir).normalize())),gamma)));
        return answer;
    }
    
};
// Camera class.

class Camera
{
    vec3 eye;        // World space camera position.
    vec3 lookAt;    // Center of window in world space.
    vec3 right;        // Vector from window center to window right-mid (in world space).
    vec3 up;        // Vector from window center to window top-mid (in world space).
    
    public:
    Camera()
    {
        eye = vec3(0, 0, 2);
        lookAt = vec3(0, 0, 1);
        right = vec3(1, 0, 0);
        up = vec3(0, 1, 0);
    }
    vec3 getEye()
    {
        return eye;
    }
    
    // Compute ray through pixel at normalized device coordinates.
    
    vec3 rayDirFromNdc(float x, float y) {
        return (lookAt - eye
                + right * x
                + up    * y
                ).normalize();
    }
};

// Ray structure.

class Ray
{
    public:
    vec3 origin;
    vec3 dir;
    Ray(vec3 o, vec3 d)
    {
        origin = o;
        dir = d;
    }
};

// Hit record structure. Contains all data that describes a ray-object intersection point.

class Hit
{
    public:
    Hit()
    {
        t = -1;
    }
    float t;                // Ray paramter at intersection. Negative means no valid intersection.
    vec3 position;            // Intersection coordinates.
    vec3 normal;            // Surface normal at intersection.
    Material* material;        // Material of intersected surface.
};

// Abstract base class.

class Intersectable
{
    protected:
    Material* material;
    public:
    Intersectable(Material* material):material(material) {}
    virtual Hit intersect(const Ray& ray)=0;
};

// Simple helper class to solve quadratic equations with the Quadratic Formula [-b +- sqrt(b^2-4ac)] / 2a, and store the results.

class QuadraticRoots
{
    public:
    float t1;
    float t2;
    
    // Solves the quadratic a*t*t + b*t + c = 0 using the Quadratic Formula [-b +- sqrt(b^2-4ac)] / 2a, and sets members t1 and t2 to store the roots.
    
    QuadraticRoots(float a, float b, float c)
    {
        float discr = b * b - 4.0 * a * c;
        if ( discr < 0 ) // no roots
        {
            t1 = -1;
            t2 = -1;
            return;
        }
        float sqrt_discr = sqrt( discr );
        t1 = (-b + sqrt_discr)/2.0/a;
        t2 = (-b - sqrt_discr)/2.0/a;
    }
    
    // Returns the lesser of the positive solutions, or a negative value if there was no positive solution.
    
    float getLesserPositive()
    {
        return (0 < t1 && (t2 < 0 || t1 < t2)) ? t1 : t2;
    }
};

// Object realization.

class Sphere : public Intersectable
{
    vec3 center;
    float radius;
    public:
    Sphere(const vec3& center, float radius, Material* material):
    Intersectable(material),
    center(center),
    radius(radius)
    {
    }
    QuadraticRoots solveQuadratic(const Ray& ray)
    {
        vec3 diff = ray.origin - center;
        float a = ray.dir.dot(ray.dir);
        float b = diff.dot(ray.dir) * 2.0;
        float c = diff.dot(diff) - radius * radius;
        return QuadraticRoots(a, b, c);
        
    }
    vec3 getNormalAt(vec3 r)
    {
        return (r - center).normalize();
    }
    Hit intersect(const Ray& ray)
    {
        // This is a generic intersect that works for any shape with a quadratic equation. solveQuadratic should solve the proper equation (+ ray equation) for the shape, and getNormalAt should return the proper normal.
        
        float t = solveQuadratic(ray).getLesserPositive();
        
        Hit hit;
        hit.t = t;
        hit.material = material;
        hit.position = ray.origin + ray.dir * t;
        hit.normal = getNormalAt(hit.position);
        
        return hit;
    }
};

// CLASS PLANE COULD COME HERE
class Plane : public Intersectable
{
    vec3 normal;
    vec3 r0;
    public:
    Plane(vec3 r1,vec3 normal2, Material* material):
    Intersectable(material)
    {
        normal=normal2;
        r0=r1;
    }
    Hit intersect(const Ray& ray)
    {
        // This is a generic intersect that works for any shape with a quadratic equation. solveQuadratic should solve the proper equation (+ ray equation) for the shape, and getNormalAt should return the proper normal.
        
        float t=(r0-ray.origin).dot(normal)/(ray.dir.dot(normal));
        
        Hit hit;
        hit.t = t;
        hit.material = material;
        hit.position = ray.origin + ray.dir * t;
        hit.normal = normal;
        
        return hit;
    }
    
};
// CLASS QUADRIC COULD COME HERE
class Quadric: public Intersectable
{
    mat4x4 coeffs;
    
    public:
    Quadric(Material* material=0,mat4x4 c= mat4x4(1,0,0,0,
                                                  0,2,0,0,
                                                  0,0,1,0,
                                                  0,0,0,-1)):
    Intersectable(material)
    
    {
        coeffs=c;
    }
    QuadraticRoots solveQuadratic(const Ray& ray)
    {
        // vec3 diff = ray.origin - center;
        vec4 d=vec4(ray.dir.x,ray.dir.y,ray.dir.z,0);
        vec4 e=vec4(ray.origin.x,ray.origin.y,ray.origin.z,1);
        
        float a = d.dot(coeffs*d);
        
        float b = d.dot(coeffs*e)+e.dot(coeffs*d);
        float c = e.dot(coeffs*e);
        return QuadraticRoots(a, b, c);
        
    }
    vec3 getNormalAt(vec3 r)
    {
        //return (r - center).normalize();
        vec4 temp=vec4( coeffs*vec4(r)+vec4(r)*coeffs);
        vec3 temp2=vec3(temp.x,temp.y,temp.z).normalize();
        return temp2;
    }
    Hit intersect(const Ray& ray)
    {
        // This is a generic intersect that works for any shape with a quadratic equation. solveQuadratic should solve the proper equation (+ ray equation) for the shape, and getNormalAt should return the proper normal.
        
        float t = solveQuadratic(ray).getLesserPositive();
        
        Hit hit;
        hit.t = t;
        hit.material = material;
        hit.position = ray.origin + ray.dir * t;
        hit.normal = getNormalAt(hit.position);
        
        return hit;
    }
    bool contains(vec3 r)
    {
        vec4 rhomo(r);
        // evaluate implicit eq
        vec4 temp=coeffs*rhomo;
        float result=rhomo.dot(temp);
        // return true if negative
        if(result<=0.0001) return true;
        // return false if positive
        else return false;
    }
    // infinite slab, ideal for clipping
    Quadric* parallelPlanes(int axis=2) {
        coeffs._00 = (axis==1);
        coeffs._11 = (axis==2);
        coeffs._22 = (axis==3);
        coeffs._33 = -1;
        return this;
    }
    Quadric* parallelPlanes2(int axis=2) {
        coeffs._00 = (axis==1);
        coeffs._11 = (axis==2);
        coeffs._22 = (axis==3);
        coeffs._33 = -2;
        return this;
    }

    
    
    Quadric* sphere()
    {
        coeffs= mat4x4(1,0,0,0,
                       0,1,0,0,
                       0,0,1,0,
                       0,0,0,-1);
        return this;
    }
    
    Quadric* plane()
    {
        coeffs= mat4x4(0,0,0,0,
                       0,1,0,0,
                       0,0,0,0,
                       0,0,0,-1);
        return this;
    }
    
    Quadric* cylinder()
    {
        coeffs= mat4x4(1,0,0,0,
                       0,0,0,0,
                       0,0,1,0,
                       0,0,0,-0.5);
        return this;
    }
    Quadric* cone()
    {
        coeffs= mat4x4(1,0,0,0,
                       0,-1,0,0,
                       0,0,1,0,
                       0,0,0,0);
        return this;
    }
    Quadric* paraboloid()
    {
        coeffs= mat4x4(1,0,0,0,
                       0,0,0,0,
                       0,0,1,0,
                       0,-1,0,0);
        return this;
        
    }
    Quadric* transform(mat4x4 t)
    {
        mat4x4 inv=t.invert();
        coeffs=inv*coeffs*inv.transpose();
        return this;
    }
    
};
// CLASS CLIPPEDQUADRIC COULD COME HERE
class ClippedQuadric:public Intersectable
{
    public:
    Quadric shape;
    Quadric clipper;
    ClippedQuadric(Material* material):
    shape(Quadric(material)),clipper(Quadric(material)),
    Intersectable(material)
    {
        
    }
    virtual Hit intersect(const Ray& ray)
    {
        // This is a generic intersect that works for any shape with a quadratic equation. solveQuadratic should solve the proper equation (+ ray equation) for the shape, and getNormalAt should return the proper normal.
        
        QuadraticRoots temp=shape.solveQuadratic(ray);
        vec3 intersection1=ray.origin + ray.dir * temp.t1;
        vec3 intersection2=ray.origin + ray.dir * temp.t2;
        
        if(!(clipper.contains(intersection1)))
        {
            temp.t1=-1;
        }
        if(!(clipper.contains(intersection2)))
        {
            temp.t2=-1;
        }
        
        float t = temp.getLesserPositive();
        
        Hit hit;
        hit.t = t;
        hit.material = material;
        hit.position = ray.origin + ray.dir * t;
        hit.normal = shape.getNormalAt(hit.position);
        
        return hit;

    }
    ClippedQuadric* transform(mat4x4 t)
    {
        shape.transform(t);
        clipper.transform(t);
        return this;
    }
    
    ClippedQuadric* cylinder(float height)
    {
        shape.cylinder();
        clipper.parallelPlanes()->transform(mat4x4::scaling(vec3(0,2/height,0)));
        return this;
    }
    
    ClippedQuadric* cylinder2(float height)
    {
        shape.cylinder();
        clipper.parallelPlanes(1)->transform(mat4x4::scaling(vec3(0,2/height,0)));
        return this;
    }
    ClippedQuadric* cone(float height)
    {
        shape.cone()->transform(mat4x4::translation(vec3(0,4*height,0)));
        clipper.parallelPlanes()->transform(mat4x4::scaling(vec3(0,1/height,0)));
        return this;
    }
    
    ClippedQuadric* sphere(float height)
    {
        shape.sphere();
        clipper.parallelPlanes()->transform(mat4x4::scaling(vec3(1,height,1))*mat4x4::translation(vec3(0,1-height,0)));
        
        return this;
    }
    
    ClippedQuadric* plane()
    {
        shape.plane();
        clipper.parallelPlanes();
        
        return this;
    }
    
    
};
class Box:public Intersectable
{
public:
    Quadric plane1;
    Quadric plane2;
    Quadric plane3;
    
    Box(Material* material):
    plane1(Quadric(material)),plane2(Quadric(material)),plane3(Quadric(material)),
    Intersectable(material)
    {
        
    }
    
    Box* box2()
    {
        plane1.parallelPlanes(1);
        plane2.parallelPlanes(2);
        plane3.parallelPlanes(3);
        return this;
    }
    Box* transform(mat4x4 t)
    {
        plane1.transform(t);
        plane2.transform(t);
        plane3.transform(t);
        return this;
    }
    Hit intersect(const Ray& ray)
    {
        // This is a generic intersect that works for any shape with a quadratic equation. solveQuadratic should solve the proper equation (+ ray equation) for the shape, and getNormalAt should return the proper normal.
        
        Hit hits[3];
        hits[0]=plane1.intersect(ray);
        hits[1]=plane2.intersect(ray);
        hits[2]=plane3.intersect(ray);
        
        float tMin=100000000;
        for(Hit hitIterator:hits)
        {
            //if hit is outside the box, ignore it by setting its t to -1
           if(!contains(hitIterator.position))
           {
               hitIterator.t=-1;
               
           }
            if(hitIterator.t<tMin&&hitIterator.t>0)
            {
                tMin=hitIterator.t;
            }
        }
        
  
        float t=tMin;
        Hit hit;
        hit.t = t;
        hit.material = material;
        hit.position = ray.origin + ray.dir * t;
        
        if(plane1.contains(hit.position)&&plane2.contains(hit.position))
        {
            hit.normal = hits[2].normal;
            hit.position=hits[2].position;
        }
        else if(plane2.contains(hit.position)&&plane3.contains(hit.position))
        {
            hit.normal = hits[0].normal;
            hit.position=hits[0].position;
        }
        
        else if(plane3.contains(hit.position)&&plane1.contains(hit.position))
        {
            hit.normal = hits[1].normal;
            hit.position=hits[1].position;
        }
        else// Doesnt matter
        {
            hit.t=-1;
            hit.normal = plane2.getNormalAt(hit.position);
        }
        
        return hit;
        
    }
    bool contains(vec3 r)
    {
        return (plane1.contains(r)&&plane2.contains(r)&&plane3.contains(r));
    }
};

class Pinking:public ClippedQuadric
{
public:
    
    Pinking(Material* material):
    ClippedQuadric(material)
    
    {
        
    }
     Hit intersect(const Ray& ray)
    {
       Hit hit= ClippedQuadric::intersect(ray);
        
        float scale = 10;
        float turbulence = 500;
        float period = 10;
        float sharpness = 1;
        float w = hit.position.x * period + pow(snoise(hit.position * scale), sharpness)*turbulence + 10000.0;
        w -= int(w);
       // std::cout<<w<<std::endl;
        if(w<0.5)
        {
            hit.t=-1;
        }
        
        return hit;
        
        
    }
};
class LightSource
{
    public:
    
    vec3 position;
    vec3 PowerDensity;
    vec3 LightDir;
    float DistanceFrom;
    
    //    LightSource(vec3 p=vec3(1,1,1),vec3 pd=vec3(1,1,1),vec3 LD=vec3(1,1,1),float DistFrom=1)
    //    {
    //        position=p;
    //        PowerDensity=pd;
    //        LightDir=LD;
    //        DistanceFrom=DistFrom;
    //    }
    public:
    virtual vec3 getPowerDensityAt ( vec3 x )=0;
    virtual vec3 getLightDirAt     ( vec3 x )=0;
    virtual float  getDistanceFrom ( vec3 x )=0;
};
class DirectionalLight : public LightSource
{
    
    public:
    DirectionalLight(vec3 p=vec3(1,1,1),vec3 pd=vec3(0.8,0.8,0.8),vec3 LD=vec3(1,1,1),float DistFrom=1)
    {
        position=p;
        PowerDensity=pd;
        LightDir=LD;
        DistanceFrom=DistFrom;
    }
    // constant for PowerDensityAt
    vec3 getPowerDensityAt ( vec3 x )
    {
        return PowerDensity;
    }
    
    //getLightDirAt -position ray
    vec3 getLightDirAt(vec3 x)
    {
        return position;
    }
    
    //getDistanceFrom returns max float
    float getDistanceFrom(vec3 x)
    {
        return std::numeric_limits<float>::max();
    }
    
};
class PointLight : public LightSource
{
    
    public:
    
    vec3 position;
    vec3 PowerDensity;

    
    PointLight(vec3 p=vec3(1,1,1),vec3 pd=vec3(1,0,0))
    {
        position=p;
        PowerDensity=pd;

    }
    // getPowerDensityAt divide by squared of distance
    
    vec3 getPowerDensityAt ( vec3 x )
    {
        float factor=(x-position).norm2();
        std::cout<<factor<<std::endl;
        
        return PowerDensity*(1/factor);
    }
    
    //getLightDirAt (position-x).normalize();
    vec3 getLightDirAt     ( vec3 x )
    {
        return (x-position).normalize();
    }
    
    //getDistanceFrom returns p-l
    float getDistanceFrom(vec3 x)
    {
        return (x-position).norm();
    }
    
};
class Scene
{
    public:
    Camera camera;
    std::vector<Intersectable*> objects;
    std::vector<Material*> materials;
    std::vector<LightSource*> lights;
    
    public:
    Scene()
    {
        materials.push_back(new DiffuseMaterial(vec3(0.7,0.7,0.8)));//0 Sand
        materials.push_back(new SpecularMaterial(vec3(1,1,0)));//1 Not Used
        materials.push_back(new SpecularMaterial(vec3(1,0,0)));//2 Not Used
        materials.push_back(new DiffuseMaterial(vec3(0,1,1)));//3 Parasol top
        materials.push_back(new Ball());//4 Ball
        materials.push_back(new SpecularMaterial(vec3(0,0.5,0)));//5 Palm Leaves
        materials.push_back(new Wood());//6 Box
        materials.push_back(new Metal());//7
         materials.push_back(new Water());//8 //Water
        materials.push_back(new SpecularMaterial(vec3(0,0,1)));
       

        
//        //Sand Dune
        objects.push_back(((new Quadric(materials[0]))->paraboloid())->transform(mat4x4::scaling(vec3(3,-1,2))*mat4x4::translation(vec3(-2.35,0,0))));

        //Beach Ball 1
        objects.push_back(((new Quadric(materials[7]))->sphere())->transform(mat4x4::scaling(vec3(0.125,0.125,0.125))*mat4x4::translation(vec3(2.5,-0.5,1-2))));
        
        //Beach Ball 2
        objects.push_back(((new Quadric(materials[4]))->sphere())->transform(mat4x4::scaling(vec3(0.125,0.125,0.125))*mat4x4::translation(vec3(0,-0.5,1.3))));
////
//        //Parasol
        //Parasol Top
           objects.push_back(((new ClippedQuadric(materials[3]))->sphere(0.5))->transform(mat4x4::scaling(vec3(0.125,0.25,0.125))*mat4x4::translation(vec3(-1.8,0.6,0)))); //

        //Parasol Pole
        objects.push_back(((new ClippedQuadric(materials[5]))->cylinder(0.25))->transform(mat4x4::scaling(vec3(0.03,1,0.03))*mat4x4::translation(vec3(-1.8,-0.2,0))));

//
//        //Sand Castle
        //Cylinder
         objects.push_back(((new ClippedQuadric(materials[0]))->cylinder(0.25))->transform(mat4x4::scaling(vec3(.0625,.25,.0625))*mat4x4::translation(vec3(-1.2,0,0.5))));
      //Cone
        objects.push_back(((new ClippedQuadric(materials[0]))->cone(0.25))->transform(mat4x4::scaling(vec3(.03,.1,.03))*mat4x4::translation(vec3(-1.2,0.3,0.5))));

        //Cylinder
        objects.push_back(((new ClippedQuadric(materials[0]))->cylinder(0.25))->transform(mat4x4::scaling(vec3(.0625,.25,.0625))*mat4x4::translation(vec3(-.9,0,0.35))));
        //Cone
        objects.push_back(((new ClippedQuadric(materials[0]))->cone(0.25))->transform(mat4x4::scaling(vec3(.03,.1,.03))*mat4x4::translation(vec3(-.9,0.3,0.35))));
        
        //Middle house
        objects.push_back(((new Box(materials[0]))->box2())->transform(mat4x4::scaling(vec3(0.0625,0.0625,0.125))*mat4x4::rotation(vec3(0,5,0), 1)*mat4x4::translation(vec3(-0.6,-0.08,1.1))));
        


//        //Palm Tree
//
        //Trunk
        objects.push_back(((new ClippedQuadric(materials[6]))->cone(0.25))->transform(mat4x4::scaling(vec3(.05,.1,.05))*mat4x4::translation(vec3(0,-0.1,1))));

        objects.push_back(((new ClippedQuadric(materials[6]))->cone(0.25))->transform(mat4x4::scaling(vec3(.05,.1,.05))*mat4x4::translation(vec3(0,-0.2,1))));
         objects.push_back(((new ClippedQuadric(materials[6]))->cone(0.25))->transform(mat4x4::scaling(vec3(.05,.1,.05))*mat4x4::translation(vec3(0,0,1))));

        //Leaves
        objects.push_back(((new Pinking(materials[5]))->sphere(0.5))->transform(mat4x4::scaling(vec3(0.125,0.155,0.125))*mat4x4::translation(vec3(-0.08,0,1))*mat4x4::rotation(vec3(1,0,10), 0)));
//
objects.push_back(((new Pinking(materials[5]))->sphere(0.5))->transform(mat4x4::scaling(vec3(0.125,0.155,0.125))*mat4x4::translation(vec3(0.1,0,1))*mat4x4::rotation(vec3(1,0,10), 0)));
objects.push_back(((new Pinking(materials[5]))->sphere(0.5))->transform(mat4x4::scaling(vec3(0.0625,0.0625,0.0625))*mat4x4::translation(vec3(0,-00.05,1.1))*mat4x4::rotation(vec3(1,0,10), 0)));
        
        //Ocean
         objects.push_back(((new Quadric(materials[8]))->parallelPlanes())->transform(mat4x4::scaling(vec3(1,1,1))*mat4x4::translation(vec3(0,-2.5,1))*mat4x4::rotation(vec3(0,0,0), 0)));
        
//        //Box
       objects.push_back(((new Box(materials[6]))->box2())->transform(mat4x4::scaling(vec3(0.0625,0.125,0.125))*mat4x4::rotation(vec3(0,10,0), 1)*mat4x4::translation(vec3(-0.6,-0.4,1.2))));
        
        LightSource* k= new DirectionalLight(vec3(1,1,1));
        
        
        
        lights.push_back(k);
        
        lights.push_back(new PointLight(vec3(-1.4,-1,0.2),vec3(0,0,1)));
        lights.push_back(new PointLight(vec3(0,-1,0),vec3(1,1,0)));
        
        
    
        
        
        
        // BUILD YOUR SCENE HERE
    }
    Hit firstIntersect(const Ray& ray)
    {
        Hit hit;
        float tMin=100000000;
        
        
        for(int i=0;i<(int)objects.size();i++)
        {
            Hit temp=objects[i]->intersect(ray);
            
            if(temp.t<tMin&&temp.t>0)
            {
                tMin=temp.t;
                hit=temp;
            }
            
            
        }
        
        
        
        
        return hit;
    }
    
    ~Scene()
    {
        // UNCOMMENT THESE WHEN APPROPRIATE
        for (std::vector<Material*>::iterator iMaterial = materials.begin(); iMaterial != materials.end(); ++iMaterial)
        delete *iMaterial;
        for (std::vector<Intersectable*>::iterator iObject = objects.begin(); iObject != objects.end(); ++iObject)
        delete *iObject;
    }
    
    public:
    Camera& getCamera()
    {
        return camera;
    }
    

    vec3 trace(const Ray& ray,int depth)
    {
        float epsilon=0.01;
        Hit hit=firstIntersect(ray);
        if(hit.t <0 || hit.t ==100000000)
        {
            float scale = 8;
            float turbulence = 300;
            float period = 1;
            float sharpness = 10;
            
            float w = ray.dir.x * period + pow(snoise(ray.dir* scale), sharpness)*turbulence + 10000.0;
            w -= int(w);
//            return ((vec3(0, 0.4, 0.7) * w + vec3(0, 0.7, 0.2) * (1-w)) * ray.dir.dot(ray.dir))*ray.dir;
            return ((vec3(1, 0.4, 0.7) * w + vec3(0, 0.2, 0.7) * (1-w)) * ray.dir.dot(ray.dir));
            //return ray.dir;
        }

        vec3 sum=vec3(0,0,0);
        for(LightSource* light:lights)
        {
            //Ray k(hit.position+hit.normal*0.01,lights[i]->getLightDirAt(hit.position+hit.normal*0.001));
            vec3 lightDirection=light->getLightDirAt(hit.position);
            
            Ray k(hit.position+hit.normal*epsilon,lightDirection);
            Hit hit2=firstIntersect(k);
            if(hit2.t<0||hit2.t==100000000)
                
            {
                sum+=hit.material->shade(hit.position,
                                         hit.normal,
                                         -ray.dir,light->getLightDirAt(hit.position),light->getPowerDensityAt(hit.position));
                if(dynamic_cast< Water*>(hit.material)||dynamic_cast< Metal*>(hit.material))
                {

                    if(depth==0)
                    {
                        return vec3(0,0,0);
                    }
                        
                        vec3 normal=hit.normal;
                        vec3 incoming= ray.dir;
                        vec3 origin=hit.position + hit.normal*1e-1;

                    
                        vec3 newDir=-((normal*2)*(normal.dot(incoming))) + incoming;
                        Ray reflectedRay=Ray(origin,newDir);
                        return sum+trace(reflectedRay,depth-1)*vec3(0.8,0.8,0.8);
                        
                        //alan's
                        // Gamma=2
                        // Epsilon=0.001
                        
                    
                    
                }

               
                
            }
            
        }
        
        return sum;
        
        
    }
};

Scene scene;




class FrameBuffer {
    unsigned int textureId;
    vec3 image[windowWidth * windowHeight];
    
    public:
    FrameBuffer() {
        for(int i = 0; i < windowWidth * windowHeight; i++) image[i] = vec3(0.0, 0.0, 0.0);
        
        glGenTextures(1, &textureId);
        glBindTexture(GL_TEXTURE_2D, textureId);
        
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, windowWidth, windowHeight, 0, GL_RGB, GL_FLOAT, image);
        
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    }
    
    void Bind(Shader* s)
    {
        s->UploadSamplerID();
        glBindTexture(GL_TEXTURE_2D, textureId);
        
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, windowWidth, windowHeight, 0, GL_RGB, GL_FLOAT, image);
    }
    
    bool ComputeImage()
    {
        static unsigned int iPart = 0;
        
        if(iPart >= 64)
        return false;
        for(int j = iPart; j < windowHeight; j+=64)
        {
            for(int i = 0; i < windowWidth; i++)
            {
                float ndcX = (2.0 * i - windowWidth) / windowWidth;
                float ndcY = (2.0 * j - windowHeight) / windowHeight;
                Camera& camera = scene.getCamera();
                Ray ray = Ray(camera.getEye(), camera.rayDirFromNdc(ndcX, ndcY));
                
                image[j*windowWidth + i] = scene.trace(ray,5);
            }
        }
        iPart++;
        return true;
    }
};

class Screen {
    FrameBuffer frameBuffer;
    unsigned int vao;
    
    public:
    Screen()
    {
        glGenVertexArrays(1, &vao);
        glBindVertexArray(vao);
        
        unsigned int vbo[2];
        glGenBuffers(2, &vbo[0]);
        
        glBindBuffer(GL_ARRAY_BUFFER, vbo[0]);
        static float vertexCoords[] = { -1, -1,        1, -1,        -1, 1,        1, 1 };
        
        glBufferData(GL_ARRAY_BUFFER, sizeof(vertexCoords), vertexCoords, GL_STATIC_DRAW);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, NULL);
        
        glBindBuffer(GL_ARRAY_BUFFER, vbo[1]);
        static float vertexTextureCoords[] = { 0, 0,    1, 0,        0, 1,        1, 1 };
        
        glBufferData(GL_ARRAY_BUFFER, sizeof(vertexTextureCoords), vertexTextureCoords, GL_STATIC_DRAW);
        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 0, NULL);
    }
    
    void Draw(Shader* s)
    {
        if(frameBuffer.ComputeImage())
        glutPostRedisplay();
        
        s->Run();
        frameBuffer.Bind(s);
        
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glBindVertexArray(vao);
        glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
        glDisable(GL_BLEND);
    }
};

Screen *screen = 0;


void onDisplay( ) {
    glClearColor(0.1f, 0.2f, 0.3f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    screen->Draw(shader);
    
    glutSwapBuffers();
}

void onInitialization()
{
    glViewport(0, 0, windowWidth, windowHeight);
    
    shader = new Shader();
    
    screen = new Screen();
}

void onExit()
{
    delete screen; screen = 0;
    delete shader; shader = 0;
    printf("exit");
}

int main(int argc, char * argv[]) {
    glutInit(&argc, argv);
#if !defined(__APPLE__)
    glutInitContextVersion(majorVersion, minorVersion);
#endif
    glutInitWindowSize(windowWidth, windowHeight);
    glutInitWindowPosition(100, 100);
#if defined(__APPLE__)
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH | GLUT_3_2_CORE_PROFILE);
#else
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
#endif
    glutCreateWindow("Ray Casting");
    
#if !defined(__APPLE__)
    glewExperimental = true;
    glewInit();
#endif
    
    printf("GL Vendor    : %s\n", glGetString(GL_VENDOR));
    printf("GL Renderer  : %s\n", glGetString(GL_RENDERER));
    printf("GL Version (string)  : %s\n", glGetString(GL_VERSION));
    glGetIntegerv(GL_MAJOR_VERSION, &majorVersion);
    glGetIntegerv(GL_MINOR_VERSION, &minorVersion);
    printf("GL Version (integer) : %d.%d\n", majorVersion, minorVersion);
    printf("GLSL Version : %s\n", glGetString(GL_SHADING_LANGUAGE_VERSION));
    
    glViewport(0, 0, windowWidth, windowHeight);
    
    onInitialization();
    
    glutDisplayFunc(onDisplay);
    
    glutMainLoop();
    
    onExit();
    
    return 1;
}

