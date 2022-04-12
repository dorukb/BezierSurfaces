#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#define _USE_MATH_DEFINES
#include <math.h>
#include <GL/glew.h>
//#include <OpenGL/gl3.h>   // The GL Header File
#include <GLFW/glfw3.h> // The GLFW header
#include <glm/glm.hpp> // GL Math library header
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp> 
#include <glm/gtc/quaternion.hpp>
#include <glm/gtx/quaternion.hpp>
#include <chrono>
#include <ctime>  

#define BUFFER_OFFSET(i) ((char*)NULL + (i))

using namespace std;

GLuint gProgram[2];
int gWidth, gHeight;

GLint modelingMatrixLoc[2];
GLint viewingMatrixLoc[2];
GLint projectionMatrixLoc[2];
GLint eyePosLoc[2];

glm::mat4 projectionMatrix;
glm::mat4 viewingMatrix;
glm::mat4 modelingMatrix;
glm::vec3 eyePos(0, 0, 0);

struct Vertex
{
    Vertex(GLfloat inX, GLfloat inY, GLfloat inZ) : x(inX), y(inY), z(inZ) { }
    GLfloat x, y, z;
};

struct Texture
{
    Texture(GLfloat inU, GLfloat inV) : u(inU), v(inV) { }
    GLfloat u, v;
};

struct Normal
{
    Normal(GLfloat inX, GLfloat inY, GLfloat inZ) : x(inX), y(inY), z(inZ) { }
    GLfloat x, y, z;
};

struct Face
{
	Face(int v[], int t[], int n[]) {
		vIndex[0] = v[0];
		vIndex[1] = v[1];
		vIndex[2] = v[2];
		tIndex[0] = t[0];
		tIndex[1] = t[1];
		tIndex[2] = t[2];
		nIndex[0] = n[0];
		nIndex[1] = n[1];
		nIndex[2] = n[2];
	}
    GLuint vIndex[3], tIndex[3], nIndex[3];
};

vector<Vertex> gVertices;
vector<Texture> gTextures;
vector<Normal> gNormals;
vector<Face> gFaces;

GLuint gVertexAttribBuffer, gIndexBuffer;
GLint gInVertexLoc, gInNormalLoc;
int gVertexDataSizeInBytes, gNormalDataSizeInBytes, gTextureDataSizeInBytes;


// effectively: use frag2.glsl and vert2.glsl
int activeProgramIndex = 1;

vector<GLfloat[16]> patchControlPointsX;
vector<GLfloat[16]> patchControlPointsY;
vector<GLfloat[16]> patchControlPointsYCopy;
vector<GLfloat[16]> patchControlPointsZ;

// i'th entry: all tri id's of which vertex is a member. (can appear in at most 6 triangles). 6*N*N memory overhead where N = #samples.
vector<vector<int>> triangleIndicesPerVertex;
vector<Normal> triNormals;
std::chrono::steady_clock::time_point starttime;

GLfloat controlPointsX[16];
GLfloat controlPointsY[16];
GLfloat controlPointsZ[16];

GLfloat controlPointsYCopy[16];
GLfloat controlPointsZCopy[16];

vector<GLfloat[16]> patchControlPoints;

// Global working variables.
int samples = 10; // Should be 10 by default.

int patchPerAxis = 1;
int totalPatchCount = patchPerAxis * patchPerAxis;

float increment;
int vertexCount = samples * samples * totalPatchCount;
int currIndices[3];
glm::vec3 a, b, normal;

// Theory
int fact(int n);
float bernstein3(int i, float t);
Vertex bezierSurface(float s, float t);
Vertex patchedBezierSurface(float s, float t, int patchIndex);

// Animation
void modifyControlPointsY(int colIndex, float deltaVal, int incr);
void modifyControlPointsZ(int colIndex, float deltaVal, int increment);
void patchModifyY(int colIndex, float deltaVal, int increment, int patchInd);
void patchModifyZ(int colIndex, float val, int increment, int patchInd);

// Graphics.
void calculateVertexNormals();
void updateTriangleNormals();
void updateVertexNormalsNoAlloc();
void createControlPoints();
void createMultiPatchedBezierSurface();
void createFaces(int patchInd);
void calculateTriNormal();

// TEST
void modifyCp(int cpCol, float val, int patchInd);

//void modifyCP5and8inZ(float val, int patchInd);
void modifyCP6and9inZ(float val, int patchInd);
void modifyCP7and10inZ(float val, int patchInd);

void anim(int colIndex, float v);
//TEST


bool ReadDataFromFile(const string& fileName,string&data)
{
    fstream myfile;

    // Open the input 
    myfile.open(fileName.c_str(), std::ios::in);

    if (myfile.is_open())
    {
        string curLine;

        while (getline(myfile, curLine))
        {
            data += curLine;
            if (!myfile.eof())
            {
                data += "\n";
            }
        }

        myfile.close();
    }
    else
    {
        return false;
    }

    return true;
}

GLuint createVS(const char* shaderName)
{
    string shaderSource;

    string filename(shaderName);
    if (!ReadDataFromFile(filename, shaderSource))
    {
        cout << "Cannot find file name: " + filename << endl;
        exit(-1);
    }

    GLint length = shaderSource.length();
    const GLchar* shader = (const GLchar*) shaderSource.c_str();

    GLuint vs = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vs, 1, &shader, &length);
    glCompileShader(vs);

    char output[1024] = {0};
    glGetShaderInfoLog(vs, 1024, &length, output);
    printf("VS compile log: %s\n", output);

	return vs;
}

GLuint createFS(const char* shaderName)
{
    string shaderSource;

    string filename(shaderName);
    if (!ReadDataFromFile(filename, shaderSource))
    {
        cout << "Cannot find file name: " + filename << endl;
        exit(-1);
    }

    GLint length = shaderSource.length();
    const GLchar* shader = (const GLchar*) shaderSource.c_str();

    GLuint fs = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fs, 1, &shader, &length);
    glCompileShader(fs);

    char output[1024] = {0};
    glGetShaderInfoLog(fs, 1024, &length, output);
    printf("FS compile log: %s\n", output);

	return fs;
}

void initShaders()
{
	// Create the programs

    gProgram[0] = glCreateProgram();
	gProgram[1] = glCreateProgram();

	// Create the shaders for both programs

    GLuint vs1 = createVS("vert.glsl");
    GLuint fs1 = createFS("frag.glsl");

	GLuint vs2 = createVS("vert2.glsl");
	GLuint fs2 = createFS("frag2.glsl");

	// Attach the shaders to the programs

	glAttachShader(gProgram[0], vs1);
	glAttachShader(gProgram[0], fs1);

	glAttachShader(gProgram[1], vs2);
	glAttachShader(gProgram[1], fs2);

	// Link the programs

    glLinkProgram(gProgram[0]);
	GLint status;
	glGetProgramiv(gProgram[0], GL_LINK_STATUS, &status);

	if (status != GL_TRUE)
	{
		cout << "Program 0 link failed" << endl;
		exit(-1);
	}

	glLinkProgram(gProgram[1]);
	glGetProgramiv(gProgram[1], GL_LINK_STATUS, &status);

	if (status != GL_TRUE)
	{
		cout << "Program 1 link failed" << endl;
		exit(-1);
	}

	// Get the locations of the uniform variables from both programs

	for (int i = 0; i < 2; ++i)
	{
		modelingMatrixLoc[i] = glGetUniformLocation(gProgram[i], "modelingMatrix");
		viewingMatrixLoc[i] = glGetUniformLocation(gProgram[i], "viewingMatrix");
		projectionMatrixLoc[i] = glGetUniformLocation(gProgram[i], "projectionMatrix");
		eyePosLoc[i] = glGetUniformLocation(gProgram[i], "eyePos");
	}
}

void initVBO()
{
    GLuint vao;
    glGenVertexArrays(1, &vao);
    assert(vao > 0);
    glBindVertexArray(vao);
    //cout << "vao = " << vao << endl;

	glEnableVertexAttribArray(0);
	glEnableVertexAttribArray(1);
    glEnableVertexAttribArray(2);
	assert(glGetError() == GL_NONE);

	glGenBuffers(1, &gVertexAttribBuffer);
	glGenBuffers(1, &gIndexBuffer);

	assert(gVertexAttribBuffer > 0 && gIndexBuffer > 0);

	glBindBuffer(GL_ARRAY_BUFFER, gVertexAttribBuffer);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, gIndexBuffer);

	gVertexDataSizeInBytes = gVertices.size() * 3 * sizeof(GLfloat);
	gNormalDataSizeInBytes = gNormals.size() * 3 * sizeof(GLfloat);
    gTextureDataSizeInBytes = gTextures.size() * 2 * sizeof(GLfloat);
    //cout << "tsize:" << gTextureDataSizeInBytes << "normal:"<<gNormalDataSizeInBytes << endl;

	int indexDataSizeInBytes = gFaces.size() * 3 * sizeof(GLuint);
	GLfloat* vertexData = new GLfloat [gVertices.size() * 3];
	GLfloat* normalData = new GLfloat [gNormals.size() * 3];
    GLfloat* textureData = new GLfloat[gTextures.size() * 2];
	GLuint* indexData = new GLuint [gFaces.size() * 3];

    //cout << "vert size: " << gVertices.size() << " norm size: " << gNormals.size() << " facesSize: " << gFaces.size() << endl;
    float minX = 1e6, maxX = -1e6;
    float minY = 1e6, maxY = -1e6;
    float minZ = 1e6, maxZ = -1e6;

	for (int i = 0; i < gVertices.size(); ++i)
	{
		vertexData[3*i] = gVertices[i].x;
		vertexData[3*i+1] = gVertices[i].y;
		vertexData[3*i+2] = gVertices[i].z;

        minX = std::min(minX, gVertices[i].x);
        maxX = std::max(maxX, gVertices[i].x);
        minY = std::min(minY, gVertices[i].y);
        maxY = std::max(maxY, gVertices[i].y);
        minZ = std::min(minZ, gVertices[i].z);
        maxZ = std::max(maxZ, gVertices[i].z);
	}

    //std::cout << "minX = " << minX << std::endl;
    //std::cout << "maxX = " << maxX << std::endl;
    //std::cout << "minY = " << minY << std::endl;
    //std::cout << "maxY = " << maxY << std::endl;
    //std::cout << "minZ = " << minZ << std::endl;
    //std::cout << "maxZ = " << maxZ << std::endl;

	for (int i = 0; i < gNormals.size(); ++i)
	{
		normalData[3*i] = gNormals[i].x;
		normalData[3*i+1] = gNormals[i].y;
		normalData[3*i+2] = gNormals[i].z;
	}

	for (int i = 0; i < gFaces.size(); ++i)
	{
		indexData[3*i] = gFaces[i].vIndex[0];
		indexData[3*i+1] = gFaces[i].vIndex[1];
		indexData[3*i+2] = gFaces[i].vIndex[2];
	}
    
    for (int i = 0; i < gTextures.size(); i++) 
    {   
        textureData[2 * i] = gTextures[i].u;
        textureData[2 * i+1] = gTextures[i].v;
    }


	glBufferData(GL_ARRAY_BUFFER, gVertexDataSizeInBytes + gNormalDataSizeInBytes + gTextureDataSizeInBytes, 0, GL_STATIC_DRAW);
	glBufferSubData(GL_ARRAY_BUFFER, 0, gVertexDataSizeInBytes, vertexData);
	glBufferSubData(GL_ARRAY_BUFFER, gVertexDataSizeInBytes, gNormalDataSizeInBytes, normalData);
    glBufferSubData(GL_ARRAY_BUFFER, gVertexDataSizeInBytes+ gNormalDataSizeInBytes, gTextureDataSizeInBytes, textureData);

	glBufferData(GL_ELEMENT_ARRAY_BUFFER, indexDataSizeInBytes, indexData, GL_STATIC_DRAW);

	// done copying; can free now
	delete[] vertexData;
	delete[] normalData;
	delete[] indexData;
    delete[] textureData;

	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(gVertexDataSizeInBytes));
    glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(gVertexDataSizeInBytes+ gNormalDataSizeInBytes));

}

void init() 
{
    int width, height, nrChannels;
    //unsigned char* data = stbi_load("metu_flag.jpg", &width, &height, &nrChannels, 0);
    //unsigned char* data = stbi_load("kurtbayrak.jpg", &width, &height, &nrChannels, 0);  

    //unsigned char* data = stbi_load("w3.jpg", &width, &height, &nrChannels, 0);
    //unsigned char* data = stbi_load("nilf.jpeg", &width, &height, &nrChannels, STBI_rgb_alpha);

    unsigned char* data = stbi_load("w3wp.jpg", &width, &height, &nrChannels, 0);

    cout << "Texture W:" << width << " H:" << height << endl;

    unsigned int flag;
    glGenTextures(1, &flag);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, flag);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, data);
    //glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, data);
    glGenerateMipmap(GL_TEXTURE_2D);
    stbi_image_free(data);

    srand(static_cast <unsigned> (time(0)));

    // create Control Points
    createControlPoints();

    // Create the surface together with vertices and uv and tri normals
    createMultiPatchedBezierSurface();

    // Calculate per vertex normals using Triangle normals and vertex data.
    calculateVertexNormals();

    starttime = std::chrono::steady_clock::now();

    glEnable(GL_DEPTH_TEST);
    initShaders();

    initVBO();

}
void createMultiPatchedBezierSurface() {

    gNormals.clear();
    gVertices.clear();
    gFaces.clear();
    gTextures.clear();

    vertexCount = samples * samples * totalPatchCount;
    triangleIndicesPerVertex = vector<vector<int>>(vertexCount, vector<int>());

    increment = 1.0f / (samples - 1);


    int vIndex = 0;
    int indexOffset = 0;
    float uvPiece = 1.0f / patchPerAxis;
    // if two patches in either direction, each patch gets 0.5, 0.5 sized pieces. + offset value\

    for (int k = 0; k < patchPerAxis; k++)
    {
        float vOffset = k * uvPiece;
        for (int p = 0; p < patchPerAxis; p++)
        {
            float uOffset = p * uvPiece;

            float sampledT = 0.f;
            float sampledS = 0.f;

            int patchInd = k * patchPerAxis + p;
            //cout << "tOffset for patch: " << patchInd << " is " << tOffset <<" sOffset: " << sOffset << endl;
            //cout << "UVs patch: " << patchInd << endl;
            for (int s = 0; s < samples; s++) {
                for (int t = 0; t < samples; t++)
                {
                    sampledT = t * increment;
                    sampledS = s * increment;

                    //auto vert = patchedBezierSurface(sampledS, sampledT, patchInd);
                    //cout <<"(t,s):"<< sampledT <<","<< sampledS << " V" <<vIndex<<": " << v.x << "," << v.y <<endl;
                    gVertices.push_back(patchedBezierSurface(sampledS, sampledT, patchInd));
                    vIndex++;

                    float u = sampledT * uvPiece + uOffset;
                    float v = sampledS * uvPiece + vOffset;
                    gTextures.push_back(Texture(u, v));
                }
            }
            indexOffset = patchInd * (samples * samples);
            createFaces(indexOffset);
        }
    }
}

void createFaces(int indexOffset)
{
    int faceIndex = 0;

    int vIndex[3], nIndex[3], tIndex[3];
    int sres = 0;
    int first, second, last;
    int res = samples;

    //cout << "patch: " << patchInd << " index offset " << indexOffset << endl;
    for (int s = 0; s < res; s++) {
        for (int t = 0; t < res - 1; t++) {
            first = t + s * res + indexOffset;
            second = first + 1;
            last = first + res;
            //cout << "Tri ind: " << first << "," << second << "," << last << endl;
            if (s < res - 1) {
                currIndices[0] = vIndex[0] = nIndex[0] = tIndex[0] = first;
                currIndices[1] = vIndex[1] = nIndex[1] = tIndex[1] = second;
                currIndices[2] = vIndex[2] = nIndex[2] = tIndex[2] = last;
                //cout << first << "," << second << "," << last << endl;
                triangleIndicesPerVertex[first].push_back(faceIndex);
                triangleIndicesPerVertex[second].push_back(faceIndex);
                triangleIndicesPerVertex[last].push_back(faceIndex);    

                gFaces.push_back(Face(vIndex, tIndex, nIndex));
                faceIndex++;
                calculateTriNormal();
            }
            if (s > 0) {
                last = second - res;
                currIndices[0] = vIndex[0] = nIndex[0] = tIndex[0] = second;
                currIndices[1] = vIndex[1] = nIndex[1] = tIndex[1] = first;
                currIndices[2] = vIndex[2] = nIndex[2] = tIndex[2] = last;
                //cout << second << "," << first << "," << last << endl;
                triangleIndicesPerVertex[first].push_back(faceIndex);
                triangleIndicesPerVertex[second].push_back(faceIndex);
                triangleIndicesPerVertex[last].push_back(faceIndex);

                gFaces.push_back(Face(vIndex, tIndex, nIndex));
                faceIndex++;
                calculateTriNormal();
            }
        }
    }
}
void calculateTriNormal() 
{
    int first = currIndices[0];
    int second = currIndices[1];
    int last = currIndices[2];

    a = glm::vec3(gVertices[first].x - gVertices[second].x,
        gVertices[first].y - gVertices[second].y,
        gVertices[first].z - gVertices[second].z);

    b = glm::vec3(gVertices[last].x - gVertices[second].x,
        gVertices[last].y - gVertices[second].y,
        gVertices[last].z - gVertices[second].z);

    // Nx = Ay * Bz - Az * By;    //Ny = Az * Bx - Ax * Bz;    //Nz = Ax * By - Ay * Bx;
    normal.x = a.y * b.z - a.z * b.y;
    normal.y = a.z * b.x - a.x * b.z;
    normal.z = a.x * b.y - a.y * b.x;
    normal = glm::normalize(-normal);
    //cout << "normal: " << normal.x << "," << normal.y << "," << normal.z << endl;
    triNormals.push_back(Normal(normal.x, normal.y, normal.z));
}

void createControlPoints() {
    // create Control Points
    patchControlPointsX = vector<GLfloat[16]>(totalPatchCount);
    patchControlPointsY = vector<GLfloat[16]>(totalPatchCount);
    //patchControlPointsYCopy = vector<GLfloat[16]>(totalPatchCount);
    patchControlPointsZ = vector<GLfloat[16]>(totalPatchCount);

    float xMax = 5.0f;
    float xDiv = xMax / patchPerAxis;

    float increment = xDiv / 3;
    for (int k = 0; k < patchPerAxis; k++)
    {
        // rows first. Y is same.
        float yOffset = k * xDiv;
        for (int l = 0; l < patchPerAxis; l++)
        {
            float xOffset = l * xDiv;
            int cpIndex = k * patchPerAxis + l;
            //cout << "CPs for patch: " << cpIndex << endl;
            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 4; j++) {

                    int ind = 4 * i + j;

                    patchControlPointsX[cpIndex][ind] = (j * increment) + xOffset;
                    patchControlPointsY[cpIndex][ind] = -((i * increment) + yOffset);
                    patchControlPointsZ[cpIndex][ind] = 0.f;
                    //cout << "Control point: " << patchControlPointsX[cpIndex][ind] << "," << patchControlPointsY[cpIndex][ind] << endl;
                }
            }
        }
    }

    patchControlPointsY[0][0] = 0.35f;
    patchControlPointsY[totalPatchCount - patchPerAxis][12] -= 0.35f;

}
void calculateVertexNormals() 
{
    int vertexCount = totalPatchCount * samples * samples;

    for (int i = 0; i < vertexCount; i++) {

        glm::vec3 sum(0);
        int count = triangleIndicesPerVertex[i].size();
        for(int triIndex : triangleIndicesPerVertex[i]) {
            // get that triangles normal
            sum.x = sum.x + triNormals[triIndex].x;
            sum.y = sum.y + triNormals[triIndex].y;
            sum.z = sum.z + triNormals[triIndex].z;
        }
        sum.x /= count;
        sum.y /= count;
        sum.z /= count;
        //cout << "vertex normal: " << sum.x << "," << sum.y << "," << sum.z << endl;
        gNormals.push_back(Normal(sum.x, sum.y, sum.z));
    }
    //gNormals.push_back(Normal(0.f, 0.f, 1.0f));
}

void updateVertexNormalsNoAlloc() 
{
    int vertexCount = totalPatchCount * samples * samples;
    glm::vec3 sum(0);
    for (int i = 0; i < vertexCount; i++) 
    {

        sum.x = sum.y = sum.z = 0.f;
        int count = triangleIndicesPerVertex[i].size();
        for (int triIndex : triangleIndicesPerVertex[i]) {
            // get that triangles normal
            sum.x += triNormals[triIndex].x;
            sum.y += triNormals[triIndex].y;
            sum.z += triNormals[triIndex].z;
        }
        sum.x /= count;
        sum.y /= count;
        sum.z /= count;
        //cout << "vertex normal: " << sum.x << "," << sum.y << "," << sum.z << endl;
        gNormals[i].x = sum.x;
        gNormals[i].y = sum.y;
        gNormals[i].z = sum.z;
    }


    // Recalculate/fix vertex normals along patch Borders.
    int currPatchOffset = 0;
    int otherPatchOffset = 0;
    int rightPatchOffset = 0;

    int lastRowStart = samples * (samples - 1);
    int vertexCountPerPatch = samples * samples;

    Normal blendedNormal(0, 0, 0);
     //there are patchPerAxis rows. skip last one.
    for (int row = 0; row < patchPerAxis; row++) {
        // for all intermediate rows.
        
        for (int col = 0; col < patchPerAxis; col++) {
            // for every patch in that row?
            int currPatchInd = row * patchPerAxis + col;
            currPatchOffset = currPatchInd * vertexCountPerPatch;

            int patchBelow = currPatchInd + patchPerAxis;
            otherPatchOffset = patchBelow * vertexCountPerPatch;

            int patchRight = currPatchInd + 1;
            rightPatchOffset = currPatchOffset + vertexCountPerPatch;

            // now find the last-row vertices of currPatch
            // and first-row vertices of patchBelow
            // blend their normals together.

            if (row != patchPerAxis - 1) {
                //cout << "CurrPatchOffset: " << currPatchOffset << " patchbelowoffset: " << otherPatchOffset << " rightpatchoffset: " << rightPatchOffset << endl;
                for (int i = 0; i < samples; i++)
                {
                    int currNormalInd = currPatchOffset + lastRowStart + i;
                    int otherNormalInd = otherPatchOffset + i;

                    //cout << "currNormalInd: " << currNormalInd << " otherNormInd: " << otherNormalInd << " vCount: " << vertexCount << endl;
                    // exactly samples many vertices in a row.
                    auto currV = gNormals[currNormalInd];
                    //cout << "duplicate vertices Curr: " << currV.x << "," << currV.y << "," << currV.z << endl;

                    auto other = gNormals[otherNormalInd];
                    //cout << " other: " << other.x << "," << other.y << "," << other.z << endl;

                    blendedNormal.x = (currV.x + other.x) / 2.f;
                    blendedNormal.y = (currV.y + other.y) / 2.f;
                    blendedNormal.z = (currV.z + other.z) / 2.f;

                    //test
                    /*blendedNormal.x = 0;
                    blendedNormal.y = 0;
                    blendedNormal.z = 1.0f;*/
                    //test

                    gNormals[currNormalInd].x = blendedNormal.x;
                    gNormals[currNormalInd].y = blendedNormal.y;
                    gNormals[currNormalInd].z = blendedNormal.z;

                    gNormals[otherNormalInd].x = blendedNormal.x;
                    gNormals[otherNormalInd].y = blendedNormal.y;
                    gNormals[otherNormalInd].z = blendedNormal.z;
                }
            }
            
            if (col != patchPerAxis - 1) {
                // except the last column, fix the column/vertical edges.

                for (int i = 0; i < samples; i++)
                {
                    int currNormalInd = currPatchOffset + ((i+1) * samples)-1;
                    int otherNormalInd = rightPatchOffset + (i * samples);
                    //cout << "COL: currNormalInd: " << currNormalInd << " otherNormInd: " << otherNormalInd << " vCount: " << vertexCount << endl;

                    // access vertices in the last COLUMN of the current
                    auto currV = gNormals[currNormalInd];
                    //cout << "Curr norm: " << currV.x << "," << currV.y << "," << currV.z << endl;

                    auto other = gNormals[otherNormalInd];
                    //cout << " others norm: " << other.x << "," << other.y << "," << other.z << endl;

                    //cout << " blended norm: " << blendedNormal.x << "," << blendedNormal.y << "," << blendedNormal.z << endl;

                    blendedNormal.x = (currV.x + other.x) / 2.f;
                    blendedNormal.y = (currV.y + other.y) / 2.f;
                    blendedNormal.z = (currV.z + other.z) / 2.f;
                    //test
                 /*   blendedNormal.x = 0;
                    blendedNormal.y = 0;
                    blendedNormal.z = 1.0f;*/
                    //test
                    gNormals[currNormalInd].x = blendedNormal.x;
                    gNormals[currNormalInd].y = blendedNormal.y;
                    gNormals[currNormalInd].z = blendedNormal.z;

                    gNormals[otherNormalInd].x = blendedNormal.x;
                    gNormals[otherNormalInd].y = blendedNormal.y;
                    gNormals[otherNormalInd].z = blendedNormal.z;
                }
            }
           
        }
    }

    if (totalPatchCount == 1) {
        // no artifact.
    }
    else if (totalPatchCount == 4) {
        int vPerPatch = samples * samples;
        int v1Ind = vPerPatch - 1;
        int v2Ind = (2 * vPerPatch) - samples;
        int v3ind = (2 * vPerPatch) - 1 + samples;
        int v4ind = (3 * vPerPatch);

        GLfloat x = (gNormals[v1Ind].x + gNormals[v2Ind].x + gNormals[v3ind].x + gNormals[v4ind].x) / 4.0f;
        GLfloat y = (gNormals[v1Ind].y + gNormals[v2Ind].y + gNormals[v3ind].y + gNormals[v4ind].y) / 4.0f;
        GLfloat z = (gNormals[v1Ind].z + gNormals[v2Ind].z + gNormals[v3ind].z + gNormals[v4ind].z) / 4.0f;

        gNormals[v1Ind].x = gNormals[v2Ind].x = gNormals[v3ind].x = gNormals[v4ind].x = x;
        gNormals[v1Ind].y = gNormals[v2Ind].y = gNormals[v3ind].y = gNormals[v4ind].y = y;
        gNormals[v1Ind].z = gNormals[v2Ind].z = gNormals[v3ind].z = gNormals[v4ind].z = z;
    }
    else 
    {
        int leftPatchOffset = 0;
        for (int row = 1; row < patchPerAxis - 1; row++) {
            for (int col = 1; col < patchPerAxis - 1; col++) {
                // for every patch in that row?
                int currPatchInd = row * patchPerAxis + col;
                currPatchOffset = currPatchInd * vertexCountPerPatch;

                int patchBelow = currPatchInd + patchPerAxis;
                if (patchBelow < totalPatchCount) {
                    // has below neighbor
                    otherPatchOffset = patchBelow * vertexCountPerPatch;

                    int patchRight = currPatchInd + 1;
                    if (patchRight < totalPatchCount) {
                        // has right neighbor
                        int patchLeft = currPatchInd - 1;
                        if (patchLeft < totalPatchCount) {

                            leftPatchOffset = currPatchOffset - vertexCountPerPatch;
                            int patchAbove = currPatchInd - patchPerAxis;

                            int patchTopLeft = patchAbove - 1;
                            int patchTopRight = patchAbove + 1;

                            int patchBottomLeft = patchBelow - 1;
                            int patchBottomRight = patchBelow + 1;

                            //cout << "pInd: " << currPatchInd << " has all 8 neighbors." << endl;

                            // Fix Top-Left vertex of current. 4 patches contribute.
                            // TopLeft, above, bottomLeft, current
                            int currNormalInd = currPatchOffset; //v0id
                            int topleft = patchTopLeft * vertexCountPerPatch + vertexCountPerPatch - 1;
                            int above = (patchAbove + 1)* vertexCountPerPatch - samples;
                            int left = patchLeft * vertexCountPerPatch + samples - 1;
                            //cout << "tl:" << topleft << " above:" << above << " bl:" << bottomleft << endl;

                            blendedNormal.x = (gNormals[currNormalInd].x + gNormals[topleft].x
                                             + gNormals[above].x + gNormals[left].x ) / 4.f;

                            blendedNormal.y = (gNormals[currNormalInd].y + gNormals[topleft].y
                                + gNormals[above].y + gNormals[left].y) / 4.f;

                            blendedNormal.z = (gNormals[currNormalInd].z + gNormals[topleft].z
                                + gNormals[above].z + gNormals[left].z) / 4.f;
                   
                            gNormals[currNormalInd].x = gNormals[topleft].x = gNormals[above].x = gNormals[left].x = blendedNormal.x;
                            gNormals[currNormalInd].y = gNormals[topleft].y = gNormals[above].y = gNormals[left].y = blendedNormal.y;
                            gNormals[currNormalInd].z = gNormals[topleft].z = gNormals[above].z = gNormals[left].z = blendedNormal.z;


                            // Fix Top-Right vertex of current patch.
                            // contributing patches: above, top right, right, current
                            currNormalInd = currPatchOffset + samples - 1; //v0id
                            int topRight = (patchTopRight + 1) * vertexCountPerPatch - samples;
                            above = patchAbove*vertexCountPerPatch + vertexCountPerPatch - 1;
                            int right = patchRight * vertexCountPerPatch;
                            //cout << "tl:" << topleft << " above:" << above << " bl:" << bottomleft << endl;

                            blendedNormal.x = (gNormals[currNormalInd].x + gNormals[topRight].x
                                + gNormals[above].x + gNormals[right].x) / 4.f;

                            blendedNormal.y = (gNormals[currNormalInd].y + gNormals[topRight].y
                                + gNormals[above].y + gNormals[right].y) / 4.f;

                            blendedNormal.z = (gNormals[currNormalInd].z + gNormals[topRight].z
                                + gNormals[above].z + gNormals[right].z) / 4.f;

                            gNormals[currNormalInd].x = gNormals[topRight].x = gNormals[above].x = gNormals[right].x = blendedNormal.x;
                            gNormals[currNormalInd].y = gNormals[topRight].y = gNormals[above].y = gNormals[right].y = blendedNormal.y;
                            gNormals[currNormalInd].z = gNormals[topRight].z = gNormals[above].z = gNormals[right].z = blendedNormal.z;

                            // Fix Bottom-Left vertex of current patch.
                            // contributing patches: below, bottom left, left, current
                            currNormalInd = patchRight*vertexCountPerPatch - samples;
                            left = patchLeft* vertexCountPerPatch + vertexCountPerPatch - 1;
                            int bottomleft = patchBottomLeft * vertexCountPerPatch + samples - 1;
                            int below = patchBelow * vertexCountPerPatch;
                            //cout << "tl:" << topleft << " above:" << above << " bl:" << bottomleft << endl;

                            blendedNormal.x = (gNormals[currNormalInd].x + gNormals[left].x
                                + gNormals[bottomleft].x + gNormals[below].x) / 4.f;

                            blendedNormal.y = (gNormals[currNormalInd].y + gNormals[left].y
                                + gNormals[bottomleft].y + gNormals[below].y) / 4.f;

                            blendedNormal.z = (gNormals[currNormalInd].z + gNormals[left].z
                                + gNormals[bottomleft].z + gNormals[below].z) / 4.f;

                            gNormals[currNormalInd].x = gNormals[left].x = gNormals[bottomleft].x = gNormals[below].x = blendedNormal.x;
                            gNormals[currNormalInd].y = gNormals[left].y = gNormals[bottomleft].y = gNormals[below].y = blendedNormal.y;
                            gNormals[currNormalInd].z = gNormals[left].z = gNormals[bottomleft].z = gNormals[below].z = blendedNormal.z;

                            // Fix Bottom-Right vertex of current patch.
                            // contributing patches: below, bottom left, left, current
                            currNormalInd = patchRight * vertexCountPerPatch -1;
                            right = (patchRight + 1) * vertexCountPerPatch - samples;
                            int bottomright = patchBottomRight * vertexCountPerPatch;
                            below = patchBelow * vertexCountPerPatch + samples - 1;
                            //cout << "tl:" << topleft << " above:" << above << " bl:" << bottomleft << endl;

                            blendedNormal.x = (gNormals[currNormalInd].x + gNormals[right].x
                                + gNormals[bottomright].x + gNormals[below].x) / 4.f;

                            blendedNormal.y = (gNormals[currNormalInd].y + gNormals[right].y
                                + gNormals[bottomright].y + gNormals[below].y) / 4.f;

                            blendedNormal.z = (gNormals[currNormalInd].z + gNormals[right].z
                                + gNormals[bottomright].z + gNormals[below].z) / 4.f;

                            gNormals[currNormalInd].x = gNormals[right].x = gNormals[bottomright].x = gNormals[below].x = blendedNormal.x;
                            gNormals[currNormalInd].y = gNormals[right].y = gNormals[bottomright].y = gNormals[below].y = blendedNormal.y;
                            gNormals[currNormalInd].z = gNormals[right].z = gNormals[bottomright].z = gNormals[below].z = blendedNormal.z;

                            // testt
                            /*gNormals[currNormalInd].x = gNormals[right].x = gNormals[bottomright].x = gNormals[below].x = 0.0f;
                            gNormals[currNormalInd].y = gNormals[right].y = gNormals[bottomright].y = gNormals[below].y = 0.0f;
                            gNormals[currNormalInd].z = gNormals[right].z = gNormals[bottomright].z = gNormals[below].z = 1.0f;*/
                        }
                    }
                }
            }
        }
    }
    
}
void updateTriangleNormals() 
{
    glm::vec3 a, b, normal;

    int faceCount = gFaces.size();
    for (int i = 0; i < faceCount; i++) 
    {
        auto inds = gFaces[i].vIndex;
        a = glm::vec3(
            gVertices[inds[0]].x - gVertices[inds[1]].x,
            gVertices[inds[0]].y - gVertices[inds[1]].y,
            gVertices[inds[0]].z - gVertices[inds[1]].z);

        b = glm::vec3(
            gVertices[inds[2]].x - gVertices[inds[1]].x,
            gVertices[inds[2]].y - gVertices[inds[1]].y,
            gVertices[inds[2]].z - gVertices[inds[1]].z);

        normal.x = a.z * b.y - a.y * b.z;
        normal.y = a.x * b.z - a.z * b.x;
        normal.z = a.y * b.x - a.x * b.y;
        normal = glm::normalize(-normal);
        //cout << "New normal: " << normal.x << "," << normal.y << "," << normal.z << endl;
        triNormals[i].x = normal.x;
        triNormals[i].y = normal.y;
        triNormals[i].z = normal.z;
    }
    
}
void drawModel()
{
	glBindBuffer(GL_ARRAY_BUFFER, gVertexAttribBuffer);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, gIndexBuffer);

	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(gVertexDataSizeInBytes));

	glDrawElements(GL_TRIANGLES, gFaces.size() * 3, GL_UNSIGNED_INT, 0);
    //glDrawElements(GL_LINES, gFaces.size() * 3, GL_UNiSIGNED_INT, 0);
}

void display()
{
    glClearColor(0, 0, 0, 1);
    glClearDepth(1.0f);
    glClearStencil(0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

	static float angle = 0;

	float angleRad = (float) (angle / 180.0) * M_PI;
	
	// Compute the modeling matrix
    modelingMatrix = glm::mat4(1.0);

	glm::mat4 matT = glm::translate(glm::mat4(1.0), glm::vec3(-2.0f, 2.0f, -15.0f));
    glm::mat4 matRx = glm::rotate<float>(glm::mat4(1.0), (-15. / 180.) * M_PI, glm::vec3(1.0, 0.0, 0.0));
    modelingMatrix = matT* matRx * modelingMatrix;

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    float change1 = 0.12f;
    float change2 = 0.37f;

    float sinFactor = 0.25f;
    long long msCount = std::chrono::duration_cast<std::chrono::milliseconds>(end - starttime).count();

    angle = msCount * sinFactor;

    angleRad = (float)(angle / 180.0) * M_PI;
    float sinVal = glm::sin(angleRad);
    float sinVal2 = glm::sin((1-sinVal) *1.1f);
    float change1Animated = change1 * sinVal2 * 2.f;

    float sinVal4 = glm::sin(sinVal2 * 3.f) * 0.55f;
    float change2Animated = change2 * sinVal4;
    float sinVal3 = glm::sin((sinVal+sinVal2) * 0.6f);

    float fastChange = change2 * glm::sin((0.25f*sinVal-sinVal2) * 4.f);    
    //float fastChange = change2 * glm::sin((1- sinVal2) * 8.f);

    // edit
    fastChange = change2Animated * fastChange;
    //edit

    float try0 = change1 * sinVal * 4.f;

    change1Animated *= 3.f;

    float val = try0 * 3.0f;
    anim(1, val);

    //float endchange = fastChange * 15.f;
    float endchange = fastChange * 15.f;

    float smallChange = sinVal3 *0.5f + 4.5;


    // Have first CP of first patch and (12th CP of last Row's first patch) slightly left positioned starting points!
    patchControlPointsX[0][0] = -0.25f;

    patchControlPointsX[totalPatchCount-patchPerAxis][12] = -0.25f;
    //patchControlPointsY[totalPatchCount - patchPerAxis][12] = 0.25f;


    for (int i = 0; i < patchPerAxis; i++) {
        int pIndex = (i + 1) * patchPerAxis -1; 
        int start = 2;
        patchControlPointsZ[pIndex][start] = endchange;
        patchControlPointsZ[pIndex][start+4] = endchange;
        patchControlPointsZ[pIndex][start+8] = endchange;
        patchControlPointsZ[pIndex][start+12] = endchange;

        //edit

        start = 3;
        patchControlPointsZ[pIndex][start] = change1Animated;
        patchControlPointsZ[pIndex][start + 4] = change1Animated;
        patchControlPointsZ[pIndex][start + 8] = change1Animated;
        patchControlPointsZ[pIndex][start + 12] = change1Animated;
        //edit
        start = 3;
        patchControlPointsX[pIndex][start] = smallChange;
        patchControlPointsX[pIndex][start + 4] = smallChange;
        patchControlPointsX[pIndex][start + 8] = smallChange;
        patchControlPointsX[pIndex][start + 12] = smallChange;

        if (i == 0) {
            patchControlPointsX[pIndex][3] = smallChange - (0.25f * (1.0f/patchPerAxis)) ;

        }
        else if (i == patchPerAxis - 1) {
             patchControlPointsX[pIndex][15] = smallChange - (0.6f * (1.0f/patchPerAxis));
        }
    }


    increment = 1.0f / (samples-1);
    float sampledT = 0.f;
    float sampledS = 0.f;
    for (int k = 0; k < patchPerAxis; k++)
    {
        for (int p = 0; p < patchPerAxis; p++)
        {
            int patchInd = k * patchPerAxis + p;
            //cout << "tOffset for patch: " << patchInd << " is " << tOffset << " sOffset: " << sOffset << endl;

            for (int s = 0; s < samples; s++) {
                for (int t = 0; t < samples; t++) {
                    //glm::vec3 vertex = bezierSurface(s, t);

                    sampledT = t * increment;
                    sampledS = s * increment;
                    //cout << "t:" << sampledT << "s:" << sampledS << endl;

                    auto v = patchedBezierSurface(sampledS, sampledT, patchInd);
                    //cout << "(t,s):" << sampledT << "," << sampledS << " V:" << v.x << "," << v.y << endl;
                    gVertices[patchInd * samples * samples + (s * samples + t)] = v;
                }
            }

        }
    }

    updateTriangleNormals();
    updateVertexNormalsNoAlloc(); // update vertex normals.

    initVBO();

	// Set the active program and the values of its uniform variables
	glUseProgram(gProgram[activeProgramIndex]);
	glUniformMatrix4fv(projectionMatrixLoc[activeProgramIndex], 1, GL_FALSE, glm::value_ptr(projectionMatrix));
	glUniformMatrix4fv(viewingMatrixLoc[activeProgramIndex], 1, GL_FALSE, glm::value_ptr(viewingMatrix));
	glUniformMatrix4fv(modelingMatrixLoc[activeProgramIndex], 1, GL_FALSE, glm::value_ptr(modelingMatrix));
	glUniform3fv(eyePosLoc[activeProgramIndex], 1, glm::value_ptr(eyePos));

	// Draw the scene
    drawModel();
}

void modifyControlPointsY(int colIndex, float deltaVal, int increment) 
{
    if (increment == 1) {
        for (int i = colIndex; i <= colIndex+3; i += increment) {
            controlPointsY[i] = controlPointsYCopy[i] + deltaVal;
        }
    }
    else {
        for (int i = colIndex; i <= 15; i += increment) {
            controlPointsY[i] = controlPointsYCopy[i] + deltaVal;
        }
    }
}

void patchModifyY(int colIndex, float deltaVal, int increment, int patchInd)
{
    if (increment == 1) {
        for (int i = colIndex; i <= colIndex + 3; i += increment) {
            patchControlPointsY[patchInd][i] = patchControlPointsYCopy[patchInd][i] + deltaVal;
        }
    }
    else {
        for (int i = colIndex; i <= 15; i += increment) {
            patchControlPointsY[patchInd][i] = patchControlPointsYCopy[patchInd][i] + deltaVal;
        }
    }
}
void patchModifyZ(int colIndex, float val, int increment, int patchInd)
{
    if (increment == 1) {
        for (int i = colIndex; i <= colIndex + 3; i += increment) {
            patchControlPointsZ[patchInd][i] = val;
        }
    }
    else {
        for (int i = colIndex; i <= 15; i += increment) {
            patchControlPointsZ[patchInd][i] = val;
        }
    }
}

void anim(int colIndex, float val) {
    //int colIndex = 1;

    int col1, col2;
    if (colIndex == 0) {
        col1 = 0;
        col2 = 3;
    }
    else if (colIndex == 1) {
        col1 = 1;
        col2 = 2;
    }
    else if (colIndex == 2) {
        col1 = 2;
        col2 = 1;
    }
    else if (colIndex == 3) {
        col1 = 3;
        col2 = 0;
    }

    for (int i = 0; i < patchPerAxis; i++)
    {
        for (int j = 0; j < patchPerAxis; j++)
        {
            int patchIndex = i * patchPerAxis + j;
            if (col1 == 0 && patchIndex % patchPerAxis == 0) continue;

            if (j % 2 == 0) {
                for (int j = col1; j < 16; j += 4) {

                    //patchControlPointsZ[patchIndex][j] = val;
                    patchControlPointsZ[patchIndex][j] = val;
                }
            }
            else {
                if (col1 == 0 || col1 == 3) 
                {
                    for (int j = col2; j < 16; j += 4) {

                        patchControlPointsZ[patchIndex][j] = -val;
                    }
                }
                else {
                    for (int j = col2; j < 16; j += 4) {

                        patchControlPointsZ[patchIndex][j] = -val;
                    }
                }
                
            }
        }
    }
}
void modifyCp(int cpCol, float val, int patchInd) 
{
    int rightPatch = patchInd + 1;
    int leftPatch = patchInd - 1;
    int patchBelow = patchInd + patchPerAxis;
    int patchAbove = patchInd - patchPerAxis;

    if (cpCol == 1) 
    {
        patchControlPointsZ[patchInd][1] = val;
        patchControlPointsZ[patchInd][5] = val;
        patchControlPointsZ[patchInd][9] = val;
        patchControlPointsZ[patchInd][13] = val;
        if (patchBelow < totalPatchCount) {
            // one below patch exists.
            patchControlPointsZ[patchBelow][1] = val;
            patchControlPointsZ[patchBelow][5] = -val;
            patchControlPointsZ[patchBelow][9] = -val;
            patchControlPointsZ[patchBelow][13] = -val;
        }
        // also check for left & top patches?
    }
    else if (cpCol == 2) 
    {
        patchControlPointsZ[patchInd][2] = val;
        patchControlPointsZ[patchInd][6] = val;
        patchControlPointsZ[patchInd][10] = val;
        patchControlPointsZ[patchInd][14] = val;

        if (rightPatch < totalPatchCount && ((patchInd+1)%patchPerAxis != 0))
        {
            // if not the last patch, and not a patch in the last column.

            // fix right neighbor patch cp's
            patchControlPointsZ[rightPatch][1] = -val;
            patchControlPointsZ[rightPatch][5] = -val;
            patchControlPointsZ[rightPatch][9] = -val;
            patchControlPointsZ[rightPatch][13] = -val;
        }
        if (patchBelow < totalPatchCount) {
            // one below patch exists.
            patchControlPointsZ[patchBelow][2] = val;
            patchControlPointsZ[patchBelow][6] = -val;
            patchControlPointsZ[patchBelow][10] = -val;
            patchControlPointsZ[patchBelow][14] = -val;

        }
    }
}
void modifyCP6and9inZ(float val, int patchInd) {
    patchControlPointsZ[patchInd][6] = val;
    patchControlPointsZ[patchInd][9] = val;
}
void modifyCP7and10inZ(float val, int patchInd) {
    patchControlPointsZ[patchInd][7] = val;
    patchControlPointsZ[patchInd][10] = val;
}
void modifyControlPointsZ(int colIndex, float deltaVal, int increment)
{
    if (increment == 1) {
        for (int i = colIndex; i <= colIndex + 3; i += increment) {
            controlPointsZ[i] = controlPointsZCopy[i] + deltaVal;
        }
    }
    else {
        for (int i = colIndex; i <= 15; i += increment) {
            controlPointsZ[i] = controlPointsZCopy[i] + deltaVal;
        }

    }
}

Vertex bezierSurface(float s, float t)
{
    float resX = 0;
    float resY = 0;
    float resZ = 0;

    for (int i = 0; i <= 3; i++) {
        for (int j = 0; j <= 3; j++) {
            float bernsteinTerm = bernstein3(i, s) * bernstein3(j, t);
            resX += bernsteinTerm * controlPointsX[4*j+i];
            resY += bernsteinTerm * controlPointsY[4*j+i];
            resZ += bernsteinTerm * controlPointsZ[4*j+i];
        }
    }
    return Vertex(resX, resY, resZ);
}

Vertex patchedBezierSurface(float s, float t, int patchIndex)
{
    float resX = 0;
    float resY = 0;
    float resZ = 0;
    int ind = 0;
    for (int i = 0; i <= 3; i++) {
        for (int j = 0; j <= 3; j++) {
            float bernsteinTerm = bernstein3(i, s) * bernstein3(j, t);
            ind = 4 * i + j;
            resX += bernsteinTerm * patchControlPointsX[patchIndex][ind];
            resY += bernsteinTerm * patchControlPointsY[patchIndex][ind];
            resZ += bernsteinTerm * patchControlPointsZ[patchIndex][ind];
        }
    }
    return Vertex(resX, resY, resZ);
}

float bernstein3(int i, float t) {
    float factTerm = 6 / (fact(i) * fact(3 - i));
    return factTerm * pow(t, i) * pow(1 - t, 3 - i);
}

int fact(int n) {
    int res = 1;
    for (int i = 1; i <= n; i++) {
        res *= i;
    }
    return res;
}
void reshape(GLFWwindow* window, int w, int h)
{
    w = w < 1 ? 1 : w;
    h = h < 1 ? 1 : h;

    gWidth = w;
    gHeight = h;

    glViewport(0, 0, w, h);

    //glMatrixMode(GL_PROJECTION);
    //glLoadIdentity();
    //glOrtho(-10, 10, -10, 10, -10, 10);
    //gluPerspective(45, 1, 1, 100);

	// Use perspective projection

	float fovyRad = (float) (45.0 / 180.0) * M_PI;
	projectionMatrix = glm::perspective(fovyRad, 1.0f, 1.0f, 100.0f);

	// Assume default camera position and orientation (camera is at
	// (0, 0, 0) with looking at -z direction and its up vector pointing
	// at +y direction)

	viewingMatrix = glm::mat4(1);

    //glMatrixMode(GL_MODELVIEW);
    //glLoadIdentity();
}
void keyboard(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if (key == GLFW_KEY_Q && action == GLFW_PRESS)
    {
        glfwSetWindowShouldClose(window, GLFW_TRUE);
    }
    else if (key == GLFW_KEY_W && action == GLFW_PRESS)
    {
        // Tapping W button should increase SAMPLES value by 1. There is no upper value limit.
        samples++;
        vertexCount = samples * samples * totalPatchCount;
        createMultiPatchedBezierSurface();
        calculateVertexNormals();

        initVBO();
    }
    else if (key == GLFW_KEY_S && action == GLFW_PRESS)
    {
        //• Tapping S button should decrease SAMPLES value by 1, down to a minimum value of 2
        if (samples > 3) {
            samples--;
            // what should change now that samples has changed?
            // Everything except Control Point calculation.
            vertexCount = samples * samples * totalPatchCount;

            createMultiPatchedBezierSurface();
            calculateVertexNormals();

            initVBO();
        }

    }
    else if (key == GLFW_KEY_E && action == GLFW_PRESS)
    {
        //Tapping E button should increase PATCH PER AXIS value by 1. There is no upper value limit.
        patchPerAxis++;
        totalPatchCount = patchPerAxis * patchPerAxis;

        // create Control Points
        createControlPoints();

        // Create the surface together with vertices and uv and tri normals
        createMultiPatchedBezierSurface();

        // Calculate per vertex normals using Triangle normals and vertex data.
        calculateVertexNormals();

        initVBO();
    }
    else if (key == GLFW_KEY_D && action == GLFW_PRESS)
    {
        //Tapping D button should decrease PATCH PER AXIS value by 1, down to a minimum value of 1
        if (patchPerAxis > 1) 
        {
            patchPerAxis--;
            totalPatchCount = patchPerAxis * patchPerAxis;

            // create Control Points
            createControlPoints();

            // Create the surface together with vertices and uv and tri normals
            createMultiPatchedBezierSurface();

            // Calculate per vertex normals using Triangle normals and vertex data.
            calculateVertexNormals();

            initVBO();
        }
    }
}

void mainLoop(GLFWwindow* window)
{
    while (!glfwWindowShouldClose(window))
    {
        display();
        glfwSwapBuffers(window);
        glfwPollEvents();
    }
}

int main(int argc, char** argv)   // Create Main Function For Bringing It All Together
{
    GLFWwindow* window;
    if (!glfwInit())
    {
        exit(-1);
    }

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 6);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    //glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_COMPAT_PROFILE);
    //glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

    int width = 640, height = 480;
    window = glfwCreateWindow(width, height, "Simple Example", NULL, NULL);

    if (!window)
    {
        glfwTerminate();
        exit(-1);
    }

    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    // Initialize GLEW to setup the OpenGL Function pointers
    if (GLEW_OK != glewInit())
    {
        std::cout << "Failed to initialize GLEW" << std::endl;
        return EXIT_FAILURE;
    }

    char rendererInfo[512] = {0};
    strcpy_s(rendererInfo, (const char*) glGetString(GL_RENDERER));
    strcat_s(rendererInfo, " - ");
    strcat_s(rendererInfo, (const char*) glGetString(GL_VERSION));
    glfwSetWindowTitle(window, rendererInfo);

    init();

    glfwSetKeyCallback(window, keyboard);
    glfwSetWindowSizeCallback(window, reshape);

    reshape(window, width, height); // need to call this once ourselves
    mainLoop(window); // this does not return unless the window is closed

    glfwDestroyWindow(window);
    glfwTerminate();

    return 0;
}
