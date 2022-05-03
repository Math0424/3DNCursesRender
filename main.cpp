#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include <string>
#include <math.h>
#include <map>
#include <time.h>
#include <sys/time.h>
#include <ncurses.h>

#define ABS(a) ((a) > 0 ? (a) : (a) * -1)
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define LERP(a, b, x) ((a) + ((b) - (a)) * (x))

const char brightness[] = "$@B%8&WM#*oahkbdpqwmZO0QLCJUYXzcvunxrjft/\\|()1{}[]?-_+~<>i!lI;:,\"^`'.";

//fast inverse sqrt
float invSqrt(float n) {
   const float threehalfs = 1.5F;
   float y = n;
   long i = * ( long * ) &y;
   i = 0x5f3759df - ( i >> 1 );
   y = * ( float * ) &i;
   y = y * ( threehalfs - ( (n * 0.5F) * y * y ) );
   return y;
}

typedef struct vector {
    float x, y, z;
    float dot(const vector a) {
        return (x * a.x) + (y * a.y) + (z * a.z);  
    }
    vector cross(const vector b) {
        return vector{y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x};  
    }
    vector operator *(float value) {
        return {x * value, y * value, z * value};
    }
    void operator *=(vector value) {
        x *= value.x;
        y *= value.y;
        z *= value.z;
    }
    void operator *=(float value) {
        x *= value;
        y *= value;
        z *= value;
    }
    void operator +=(vector value) {
        x += value.x;
        y += value.y;
        z += value.z;
    }
    void operator +=(float value) {
        x += value;
        y += value;
        z += value;
    }
    vector operator *(vector value) {
        return {x * value.x, y * value.y, z * value.z};
    }
    vector operator -(vector value) {
        return {x - value.x, y - value.y, z - value.z};
    }
    vector operator +(vector value) {
        return {x + value.x, y + value.y, z + value.z};
    }
    vector operator /(float value) {
        return {x / value, y / value, z / value};
    }
    float length() {
        return sqrt((x * x) + (y * y) + (z * z));
    }
    void normalize() {
        float inv = invSqrt((x * x) + (y * y) + (z * z));
        x *= inv;
        y *= inv;
        z *= inv;
    }
    vector normalized() {
        float inv = invSqrt((x * x) + (y * y) + (z * z));
        return {x * inv, y * inv, z * inv};
    }
} vector;

typedef struct matrix4x4 {
    //identity matrix
    float m11 = 1, m12 = 0, m13 = 0, m14 = 0;
    float m21 = 0, m22 = 1, m23 = 0, m24 = 0;
    float m31 = 0, m32 = 0, m33 = 1, m34 = 0;
    float m41 = 0, m42 = 0, m43 = 0, m44 = 1;

    //https://stackoverflow.com/questions/18404890/how-to-build-perspective-projection-matrix-no-api
    void createProjectionMatrix(float aspectRatio, float fov, float zFar, float zNear) {
        float rad = 1.0f / tanf(fov * 0.5f / 180.0f * 3.14159f);
        m11 = aspectRatio * rad;
        m22 = rad;
        m33 = zFar / (zFar - zNear);
        m43 = (-zFar * zNear) / (zFar - zNear);
        m34 = 1;
        m44 = 0;
    }

    vector transform(const vector in) {
        float v = 1 / (in.x * m14 + in.y * m24 + in.z * m34 + m44);
        return {
            v * ( (in.x * m11 + in.y * m21 + in.z * m31) + m41 ),
            v * ( (in.x * m12 + in.y * m22 + in.z * m32) + m42 ),
            v * ( (in.x * m13 + in.y * m23 + in.z * m33) + m43 ),
        };
    }

    vector transformNormal(const vector in) {
        return {
            in.x * m11 + in.y * m21 + in.z * m31,
            in.x * m12 + in.y * m22 + in.z * m32,
            in.x * m13 + in.y * m23 + in.z * m33,
        };
    }

    void createRotationX(float radians) {
        float num = cos(radians);
        float num2 = sin(radians);
        m11 = num;
        m13 = -num2;
        m31 = num2;
        m33 = num;
    }

    void createRotationY(float radians) {
        float num = cos(radians);
        float num2 = sin(radians);
        m22 = num;
        m23 = num2;
        m32 = -num2;
        m33 = num;
    }

    void createRotationZ(float radians) {
        float num = cos(radians);
        float num2 = sin(radians);
        m11 = num;
        m12 = num2;
        m21 = -num2;
        m22 = num;
    }

    matrix4x4 inverted() {
        double num1 = (double)m33 * m44 - (double)m43 * m34;
        double num2 = (double)m32 * m44 - (double)m42 * m34;
        double num3 = (double)m32 * m43 - (double)m42 * m33;
        double num4 = (double)m31 * m44 - (double)m41 * m34;
        double num5 = (double)m31 * m43 - (double)m41 * m33;
        double num6 = (double)m31 * m42 - (double)m41 * m32;
        double num7 = (double)m23 * m44 - (double)m43 * m24;
        double num8 = (double)m22 * m44 - (double)m42 * m24;
        double num9 = (double)m22 * m43 - (double)m42 * m23;
        double numa = (double)m21 * m44 - (double)m41 * m24;
        double numb = (double)m21 * m43 - (double)m41 * m23;
        double numc = (double)m21 * m42 - (double)m41 * m22;
        double numd = (double)m23 * m34 - (double)m33 * m24;
        double nume = (double)m22 * m34 - (double)m32 * m24;
        double numf = (double)m22 * m33 - (double)m32 * m23;
        double numg = (double)m21 * m34 - (double)m31 * m24;
        double numh = (double)m21 * m33 - (double)m31 * m23;
        double numi = (double)m21 * m32 - (double)m31 * m22;

        double a11 =  ( (m22 * num1) - (m23 * num2) + (m24 * num3) );
        double a12 = -( (m21 * num1) - (m23 * num4) + (m24 * num5) );
        double a13 =  ( (m21 * num2) - (m22 * num4) + (m24 * num6) );
        double a14 = -( (m21 * num3) - (m22 * num5) + (m23 * num6) );
        double a21 = -( (m12 * num1) - (m13 * num2) + (m14 * num3) );
        double a22 =  ( (m11 * num1) - (m13 * num4) + (m14 * num5) );
        double a23 = -( (m11 * num2) - (m12 * num4) + (m14 * num6) );
        double a24 =  ( (m11 * num3) - (m12 * num5) + (m13 * num6) );
        double a31 =  ( (m12 * num7) - (m13 * num8) + (m14 * num9) );
        double a32 = -( (m11 * num7) - (m13 * numa) + (m14 * numb) );
        double a33 =  ( (m11 * num8) - (m12 * numa) + (m14 * numc) );
        double a34 = -( (m11 * num9) - (m12 * numb) + (m13 * numc) );
        double a41 = -( (m12 * numd) - (m13 * nume) + (m14 * numf) );
        double a42 =  ( (m11 * numd) - (m13 * numg) + (m14 * numh) );
        double a43 = -( (m11 * nume) - (m12 * numg) + (m14 * numi) );
        double a44 =  ( (m11 * numf) - (m12 * numh) + (m13 * numi) );
        
        double det = 1 / ((m11 * a11) + (m12 * a12) + (m13 * a13) + (m14 * a14));

        return matrix4x4 {
            (float)(a11 * det), (float)(a21 * det), (float)(a31 * det), (float)(a41 * det),
            (float)(a12 * det), (float)(a22 * det), (float)(a32 * det), (float)(a42 * det),
            (float)(a13 * det), (float)(a23 * det), (float)(a33 * det), (float)(a43 * det),
            (float)(a14 * det), (float)(a24 * det), (float)(a34 * det), (float)(a44 * det),
        };
    }
    matrix4x4 transposed() {
        return matrix4x4 {
            m11, m21, m31, m41,
            m12, m22, m32, m42,
            m13, m23, m33, m43,
            m14, m24, m34, m44,
        };
    }
    vector getUp() {
        return { m11, m21, m31 };
    }
    vector getForward() {
        return { m12, m22, m32 };
    }
    vector getLeft() {
        return { m13, m23, m33 };
    }
    vector getPos() {
        return { m14, m24, m34 };
    }
} matrix4x4;


typedef struct triangle {
    vector vA, vB, vC;
    vector normal;

    void operator *=(matrix4x4 matrix) {
        vA = matrix.transform(vA);
        vB = matrix.transform(vB);
        vC = matrix.transform(vC);
        
        //todo: fix normal transfomation
        normal = matrix.inverted().transformNormal(normal);
    }

    void operator *=(float vec) {
        vA *= vec;
        vB *= vec;
        vC *= vec;
    }

    void operator +=(vector vec) {
        vA += vec;
        vB += vec;
        vC += vec;
    }

    vector calculateNormal() {
        vector r = (vA - vC).cross(vB - vA);
        r.normalize();
        return r;
    }

    void drawTriangle(matrix4x4 *projection) {

    }
} triangle;


float scale = 1;
std::string objectName;
std::map<int, vector> normals;
std::map<int, vector> vertecies;
std::vector<triangle> triangles;

std::vector<std::string> split(std::string str, char delimiter) {
    std::vector<std::string> vec;
    int curr = -1, tmp;
    do {
        tmp = curr + 1;
        curr = str.find(delimiter, tmp);
        vec.push_back(str.substr(tmp, curr - tmp));
    } while(curr != -1);
    return vec;
}

void parseLine(std::string line) {
    std::vector<std::string> args = split(line, ' ');
    if (args[0] == "#") {
        if (vertecies.size() == 0) {
            printf("%s\n", line.c_str());
        }
    } else if (args[0] == "o") {
        objectName = line.substr(2);
    } else if(args[0] == "v") {
        vertecies.insert(
            std::pair<int, vector>(
                vertecies.size() + 1, 
                {
                    std::stof(args[1]), 
                    std::stof(args[2]), 
                    std::stof(args[3])
                }
            )
        );
    } else if(args[0] == "vn") {
        normals.insert(
            std::pair<int, vector>(
                normals.size() + 1, 
                {
                    std::stof(args[1]), 
                    std::stof(args[2]), 
                    std::stof(args[3])
                }
            )
        );
    } else if(args[0] == "f") {
        std::vector<std::string> vnA = split(args[1], '/');
        std::vector<std::string> vnB = split(args[2], '/');
        std::vector<std::string> vnC = split(args[3], '/');
        triangles.push_back(
            triangle {
                .vA = vertecies[std::stoi(vnA[0])],
                .vB = vertecies[std::stoi(vnB[0])],
                .vC = vertecies[std::stoi(vnC[0])],
                
                .normal = (normals[std::stoi(vnA[2])] + normals[std::stoi(vnB[2])] + normals[std::stoi(vnC[2])]) / 3
            }
        );
    } else if(args[0] == "scale") {
        scale = std::stof(args[1]);
    } else {
        printf("Unknown delimiter '%s'\n", args[0].c_str());
    }
}

void drawLine(int x1, int y1, int x2, int y2, char c) {
    float dist = sqrt((x1 - x2) * (x1 - x2) + (y1 - x2) * (y1 - x2));
    for (int x = 0; x < dist; x++) {
        mvprintw(LERP(y1, y2, x / dist), LERP(x1, x2, x / dist), "%c", c);
    }
}


int main(int argc, char **argv) {
    //matrix4x4 matrix {
    //    1, 3, 5, 4,
    //    2, 4, 6, 2,
    //    4, 6, 8, 5, 
    //    3, 1, 7, 9,
    //};
    //printf("%f, %f, %f, %f\n", matrix.m11, matrix.m12, matrix.m13, matrix.m14);
    //printf("%f, %f, %f, %f\n", matrix.m21, matrix.m22, matrix.m23, matrix.m24);
    //printf("%f, %f, %f, %f\n", matrix.m31, matrix.m32, matrix.m33, matrix.m34);
    //printf("%f, %f, %f, %f\n", matrix.m41, matrix.m42, matrix.m43, matrix.m44);
    //printf("Inverted:\n");
    //matrix = matrix.inverted();
    //printf("%f, %f, %f, %f\n", matrix.m11, matrix.m12, matrix.m13, matrix.m14);
    //printf("%f, %f, %f, %f\n", matrix.m21, matrix.m22, matrix.m23, matrix.m24);
    //printf("%f, %f, %f, %f\n", matrix.m31, matrix.m32, matrix.m33, matrix.m34);
    //printf("%f, %f, %f, %f\n", matrix.m41, matrix.m42, matrix.m43, matrix.m44);

    std::ifstream file(argv[1]);

    std::string line;
    while (std::getline(file, line)) {
        parseLine(line);
    }
    file.close();

    printf("Loaded '%s' ", objectName.c_str());
    printf("- %d vertecies, %d normals, %d faces\n", (int)vertecies.size(), (int)normals.size(), (int)triangles.size());


    matrix4x4 projectionMatrix;
    projectionMatrix.createProjectionMatrix(80/40, 90, 0.1f, 100.0f);

    initscr();
    noecho();
    curs_set(false);
    system("resize -s 40 80");
    start_color();
    clear();

    init_pair(1, COLOR_WHITE, COLOR_BLACK);
    init_pair(2, COLOR_RED, COLOR_RED);
    attron(COLOR_PAIR(1));

    for(int i = 9; i < 256; i++) {
        init_color(i, i, i, i);
        init_pair(i, i, i);
    }

    int frame = 0;
    matrix4x4 rotationX, rotationZ;

    //render at 30 fps
    struct timeval start, current;
    while (true) {
        gettimeofday(&current, NULL);
        if ((current.tv_sec * 1000000 + current.tv_usec) - 
            (start.tv_sec * 1000000 + start.tv_usec) < 50000) { //20fps
            continue;
        }
        gettimeofday(&start, NULL);

        frame++;

        rotationX.createRotationX(frame / 40.0f);
        rotationZ.createRotationZ(frame / 80.0f);
        
        triangle translated;
        int triCount = 0;
        for(auto triangle : triangles) {
            translated = triangle;

            //rotate scale and transform
            translated *= rotationX;
            translated *= rotationZ;
            translated *= scale;
            translated += {0, -2, 5};

            vector normal = translated.calculateNormal();

            if (normal.dot(translated.vA) < 0) {
    
                //project to screen
                translated *= projectionMatrix;
                translated += {1, 1, 1};
                translated *= 40;

                float dotP = normal.dot( (vector{.5, .5, -1}).normalized() );
                attron(COLOR_PAIR( MAX(12, (int)(dotP * 255)) ));

                drawLine(translated.vA.x, translated.vA.y, translated.vB.x, translated.vB.y, ' ');
                drawLine(translated.vB.x, translated.vB.y, translated.vC.x, translated.vC.y, ' ');
                drawLine(translated.vC.x, translated.vC.y, translated.vA.x, translated.vA.y, ' ');
                
                triCount++;
            }

        }


        attron(COLOR_PAIR(1));
        mvprintw(0, 0, "frame %d", frame);
        mvprintw(1, 0, "%s", argv[1]);
        mvprintw(2, 0, "%d/%d tris", triCount, triangles.size());

        refresh();
        clear();
    }


    return 0;
}