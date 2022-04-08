
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

#define BUFFER_OFFSET(i) ((char*)NULL + (i))

using namespace std;

GLuint gProgram[2];
int gWidth, gHeight;

GLint modelingMatrixLoc[2];
GLint viewingMatrixLoc[2];
GLint projectionMatrixLoc[2];
GLint eyePosLoc[2];

GLfloat controlPointsX[16];
GLfloat controlPointsY[16];
GLfloat controlPointsZ[16];

glm::mat4 projectionMatrix;
glm::mat4 viewingMatrix;
glm::mat4 modelingMatrix;
glm::vec3 eyePos(0, 0, 0);

int activeProgramIndex = 0;

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
int gVertexDataSizeInBytes, gNormalDataSizeInBytes;

int fact(int n);
float bernstein3(int i, float t);
Vertex bezierSurface(float s, float t);

bool ParseObj(const string& fileName)
{
    fstream myfile;

    // Open the input 
    myfile.open(fileName.c_str(), std::ios::in);

    if (myfile.is_open())
    {
        string curLine;

        while (getline(myfile, curLine))
        {
            stringstream str(curLine);
            GLfloat c1, c2, c3;
            GLuint index[9];
            string tmp;

            if (curLine.length() >= 2)
            {
                if (curLine[0] == 'v')
                {
                    if (curLine[1] == 't') // texture
                    {
                        str >> tmp; // consume "vt"
                        str >> c1 >> c2;
                        gTextures.push_back(Texture(c1, c2));
                    }
                    else if (curLine[1] == 'n') // normal
                    {
                        str >> tmp; // consume "vn"
                        str >> c1 >> c2 >> c3;
                        gNormals.push_back(Normal(c1, c2, c3));
                    }
                    else // vertex
                    {
                        str >> tmp; // consume "v"
                        str >> c1 >> c2 >> c3;
                        gVertices.push_back(Vertex(c1, c2, c3));
                    }
                }
                else if (curLine[0] == 'f') // face
                {
                    str >> tmp; // consume "f"
					char c;
					int vIndex[3],  nIndex[3], tIndex[3];
					str >> vIndex[0]; str >> c >> c; // consume "//"
					str >> nIndex[0]; 
					str >> vIndex[1]; str >> c >> c; // consume "//"
					str >> nIndex[1]; 
					str >> vIndex[2]; str >> c >> c; // consume "//"
					str >> nIndex[2]; 

					assert(vIndex[0] == nIndex[0] &&
						   vIndex[1] == nIndex[1] &&
						   vIndex[2] == nIndex[2]); // a limitation for now

					// make indices start from 0
					for (int c = 0; c < 3; ++c)
					{
						vIndex[c] -= 1;
						nIndex[c] -= 1;
						tIndex[c] -= 1;
					}

                    gFaces.push_back(Face(vIndex, tIndex, nIndex));
                }
                else
                {
                    cout << "Ignoring unidentified line in obj file: " << curLine << endl;
                }
            }

            //data += curLine;
            if (!myfile.eof())
            {
                //data += "\n";
            }
        }

        myfile.close();
    }
    else
    {
        return false;
    }

	/*
	for (int i = 0; i < gVertices.size(); ++i)
	{
		Vector3 n;

		for (int j = 0; j < gFaces.size(); ++j)
		{
			for (int k = 0; k < 3; ++k)
			{
				if (gFaces[j].vIndex[k] == i)
				{
					// face j contains vertex i
					Vector3 a(gVertices[gFaces[j].vIndex[0]].x, 
							  gVertices[gFaces[j].vIndex[0]].y,
							  gVertices[gFaces[j].vIndex[0]].z);

					Vector3 b(gVertices[gFaces[j].vIndex[1]].x, 
							  gVertices[gFaces[j].vIndex[1]].y,
							  gVertices[gFaces[j].vIndex[1]].z);

					Vector3 c(gVertices[gFaces[j].vIndex[2]].x, 
							  gVertices[gFaces[j].vIndex[2]].y,
							  gVertices[gFaces[j].vIndex[2]].z);

					Vector3 ab = b - a;
					Vector3 ac = c - a;
					Vector3 normalFromThisFace = (ab.cross(ac)).getNormalized();
					n += normalFromThisFace;
				}

			}
		}

		n.normalize();

		gNormals.push_back(Normal(n.x, n.y, n.z));
	}
	*/

	assert(gVertices.size() == gNormals.size());

    return true;
}

bool ReadDataFromFile(
    const string& fileName, ///< [in]  Name of the shader file
    string&       data)     ///< [out] The contents of the file
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
		cout << "Program link failed" << endl;
		exit(-1);
	}

	glLinkProgram(gProgram[1]);
	glGetProgramiv(gProgram[1], GL_LINK_STATUS, &status);

	if (status != GL_TRUE)
	{
		cout << "Program link failed" << endl;
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
    cout << "vao = " << vao << endl;

	glEnableVertexAttribArray(0);
	glEnableVertexAttribArray(1);
	assert(glGetError() == GL_NONE);

	glGenBuffers(1, &gVertexAttribBuffer);
	glGenBuffers(1, &gIndexBuffer);

	assert(gVertexAttribBuffer > 0 && gIndexBuffer > 0);

	glBindBuffer(GL_ARRAY_BUFFER, gVertexAttribBuffer);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, gIndexBuffer);

	gVertexDataSizeInBytes = gVertices.size() * 3 * sizeof(GLfloat);
	gNormalDataSizeInBytes = gNormals.size() * 3 * sizeof(GLfloat);
	int indexDataSizeInBytes = gFaces.size() * 3 * sizeof(GLuint);
	GLfloat* vertexData = new GLfloat [gVertices.size() * 3];
	GLfloat* normalData = new GLfloat [gNormals.size() * 3];
	GLuint* indexData = new GLuint [gFaces.size() * 3];

    cout << "vert size: " << gVertices.size() << " norm size: " << gNormals.size() << " facesSize: " << gFaces.size() << endl;
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

    std::cout << "minX = " << minX << std::endl;
    std::cout << "maxX = " << maxX << std::endl;
    std::cout << "minY = " << minY << std::endl;
    std::cout << "maxY = " << maxY << std::endl;
    std::cout << "minZ = " << minZ << std::endl;
    std::cout << "maxZ = " << maxZ << std::endl;

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


	glBufferData(GL_ARRAY_BUFFER, gVertexDataSizeInBytes + gNormalDataSizeInBytes, 0, GL_STATIC_DRAW);
	glBufferSubData(GL_ARRAY_BUFFER, 0, gVertexDataSizeInBytes, vertexData);
	glBufferSubData(GL_ARRAY_BUFFER, gVertexDataSizeInBytes, gNormalDataSizeInBytes, normalData);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, indexDataSizeInBytes, indexData, GL_STATIC_DRAW);

	// done copying; can free now
	delete[] vertexData;
	delete[] normalData;
	delete[] indexData;

	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(gVertexDataSizeInBytes));
}

void init() 
{
    // create Control Points
   /* GLfloat controlPointsX[16];
    GLfloat controlPointsY[16];
    GLfloat controlPointsZ[16];*/
    float x, y, z;
    x = y = z = 0.f;
    float yDiff = 0.75f;
    float zDiff = 0.0f;

    srand(static_cast <unsigned> (time(0)));

    for (int i = 0; i < 4; i++) {
        x = 0.f;
        for (int j = 0; j < 4; j++) {
            controlPointsX[4 * i + j] = x + j;
            controlPointsY[4 * i + j] = y - i;

            float randZOffset = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / 2.0f));
            controlPointsZ[4 * i + j] = 0.f;
            //controlPointsZ[4 * i + j] = z - randZOffset;

            cout << "Control point: " << x + j << "," << controlPointsY[4* i + j] << "," << z << endl;
        }
    }

    //controlPointsZ[0] += 1.5f;

   /* controlPointsZ[6] += 0.25f;
    controlPointsZ[7] += 0.3f;*/
    //controlPointsZ[10] += 0.f;
    cout << "BEFORE: " << controlPointsX[15] << endl;
    //controlPointsX[3] -= 0.05f;
    cout << "AFTER: " << controlPointsX[15] << endl;

    float change = -0.35f;
    float change2 = -0.3f;
    controlPointsY[1] += change;
    controlPointsY[5] += change;
    controlPointsY[9] += change;
    controlPointsY[13] += change;

    //controlPointsY[3] -= change2;/*
    //controlPointsY[7] -= change2;
    //controlPointsY[11] -= change2;
    //controlPointsY[15] -= change2;*/

    controlPointsZ[15] += 0.5f;


    gNormals.clear();
    gVertices.clear();
    gFaces.clear();
    gTextures.clear();

    int sMax = 3;
    int tMax = 3;


    int samples = 10;
    float increment = (float)sMax / samples;
    cout << "Increment: " << increment << endl;

    vector<glm::vec3> bezierVertices;
    // create bezier surface, with 4 vertices? a rectangle really.

    float sampledT = 0.f;
    float sampledS = 0.f;

    for (int s = 0; s < samples; s++) {
        for (int t = 0; t < samples; t++) {
            //glm::vec3 vertex = bezierSurface(s, t);

            sampledT = t * increment;
            sampledS = s * increment;
            //cout << "t:" << sampledT << "s:" << sampledS << endl;

            gVertices.push_back(bezierSurface(sampledT, sampledS));
            gNormals.push_back(Normal(0.f, 0.f, 1.0f));
            gTextures.push_back(Texture((float)t, (float)s));
           
        }
    }

    int vIndex[3], nIndex[3], tIndex[3];
    int sres = 0;
    int first, second, last;
    int res = samples;

    for (int s = 0; s < res; s++) {
        for (int t = 0; t < res-1; t++) {

            //cout << "s: " << s << " t: " << t << endl;
            sres = s * res;
            first = t + sres;
            second = first + 1;
            last = first + res;

            if (s < res - 1) {
                vIndex[0] = nIndex[0] = tIndex[0] = first;
                vIndex[1] = nIndex[1] = tIndex[1] = second;
                vIndex[2] = nIndex[2] = tIndex[2] = last;
                //cout << first << "," << second << "," << last << endl;
                gFaces.push_back(Face(vIndex, tIndex, nIndex));
            }
            if (s > 0) {
                last = second - res;

                vIndex[0] = nIndex[0] = tIndex[0] = second;
                vIndex[1] = nIndex[1] = tIndex[1] = first;
                vIndex[2] = nIndex[2] = tIndex[2] = last;
                //cout << second << "," << first << "," << last << endl;
                gFaces.push_back(Face(vIndex, tIndex, nIndex));
            }
        }
    }
    // top left vertex
    //gTextures.push_back(Texture(0.f,0.f));
    //gNormals.push_back(Normal(0.f, 0.f, -1.0f));

  /*  int vIndex[3], nIndex[3], tIndex[3];
    for (int k = 0; k < 3; k++) {
        vIndex[k] = k;
        nIndex[k] = k;
        tIndex[k] = k;
    }
    gFaces.push_back(Face(vIndex, tIndex, nIndex));

    for (int k = 0; k < 3; k++) {
        vIndex[k] = 3-k;
        nIndex[k] = 3-k;
        tIndex[k] = 3-k;
    }
    gFaces.push_back(Face(vIndex, tIndex, nIndex));*/

    //for (int i = 0; i < gFaces.size(); i++) {
    //    cout << "Face:" << gFaces[i].vIndex[0] << "," << gFaces[i].vIndex[1] << "," << gFaces[i].vIndex[2] << endl;
    //}


    glEnable(GL_DEPTH_TEST);
    initShaders();

    initVBO();
}

void drawModel()
{
	glBindBuffer(GL_ARRAY_BUFFER, gVertexAttribBuffer);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, gIndexBuffer);

	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(gVertexDataSizeInBytes));

	glDrawElements(GL_TRIANGLES, gFaces.size() * 3, GL_UNSIGNED_INT, 0);
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

	modelingMatrix = glm::translate(glm::mat4(1.0), glm::vec3(-4.0f, 4.0f, -30.0f));
	//modelingMatrix = glm::rotate(modelingMatrix, angleRad, glm::vec3(0.0, 1.0, 0.0));
    //glm::mat4 matT = glm::translate(glm::mat4(1.0), glm::vec3(-0.5f, -0.4f, -5.0f));   // same as above but more clear
    ////glm::mat4 matR = glm::rotate(glm::mat4(1.0), angleRad, glm::vec3(0.0, 1.0, 0.0));
    //glm::mat4 matRx = glm::rotate<float>(glm::mat4(1.0), (-90. / 180.) * M_PI, glm::vec3(1.0, 0.0, 0.0));
    //glm::mat4 matRy = glm::rotate<float>(glm::mat4(1.0), (-90. / 180.) * M_PI, glm::vec3(0.0, 1.0, 0.0));
    //glm::mat4 matRz = glm::rotate<float>(glm::mat4(1.0), angleRad, glm::vec3(0.0, 0.0, 1.0));
    //modelingMatrix = matRy * matRx;

    // Let's make some alternating roll rotation
    static float rollDeg = 0.f;
    //static float changeRoll = 2.5;
    static float changeRoll = 0.5f;



    float startRoll = 0;
    float endRoll = 90.f;

    float rollRad = (float) (rollDeg / 180.f) * M_PI;
    rollDeg += changeRoll;
    if (rollDeg >= endRoll || rollDeg <= startRoll) // was 10.
    {
        changeRoll *= -1.f;
    }
    //glm::mat4 matRoll = glm::rotate<float>(glm::mat4(1.0), rollRad, glm::vec3(1.0, 0.0, 0.0));

    // Let's make some pitch rotation
    static float pitchDeg = 0;
    //static float changePitch = 0.1;
    static float changePitch = 0.08;

  /*  float startPitch = 0;
    float endPitch = 90;
    float pitchRad = (float) (pitchDeg / 180.f) * M_PI;
    pitchDeg += changePitch;
    if (pitchDeg >= endPitch)
    {
        changePitch = 0;
    }*/

    float startPitch = 0;
    float endPitch = 90;
    float pitchRad = (float)(pitchDeg / 180.f) * M_PI;
    pitchDeg += changePitch;
    if (pitchDeg >= endPitch)// || pitchDeg <= startPitch)
    {
        changePitch = 0.f;
    }
    //glm::mat4 matPitch = glm::rotate<float>(glm::mat4(1.0), pitchRad, glm::vec3(0.0, 0.0, 1.0));
    //modelingMatrix = matRoll * matPitch * modelingMatrix; // gimbal lock
    //modelingMatrix = matPitch * matRoll * modelingMatrix;   // no gimbal lock

     
    //glm::quat q0(0, 1, 0, 0); // along x
    //glm::quat q1(0, 0, 1, 0); // along y
    //glm::quat q = glm::mix(q0, q1, (pitchDeg - startPitch) / (endPitch - startPitch));

    //float sint = sin(rollRad / 2);
    //glm::quat rollQuat(cos(rollRad/2), sint * q.x, sint * q.y, sint * q.z);
    //glm::quat pitchQuat(cos(pitchRad/2), 0, 0, 1 * sin(pitchRad/2));
    //modelingMatrix = matT * glm::toMat4(rollQuat)  * glm::toMat4(pitchQuat) * modelingMatrix;
    // first pitch, then roll.

   
    //cout << rollQuat.w << " " << rollQuat.x << " " << rollQuat.y << " " << rollQuat.z << endl;

	// Set the active program and the values of its uniform variables


    //glm::quat q3(0, 0, 0, 1); // along +z
    //glm::quat q4(0, 0, -1, 0); // along -z?

    //float initRad = (float)(20.f / 180.f) * M_PI;
    //glm::quat rot(cos(initRad / 2), 1 * sin(initRad / 2), 0, 0);

    //glm::vec4 newAxis = glm::toMat4(rot)* glm::vec4(0, 0, 1, 1);

    //glm::quat q4(0, newAxis.x, newAxis.y, newAxis.z);

    // q4 now should be along -z.

    //glm::quat q = glm::mix(q3, q4, (rollDeg - startRoll) / (endRoll - startRoll));

    //glm::quat rollQuat(cos(rollRad / 2), 1 * sin(rollRad/2),0, 0);

    //float sint = sin(pitchRad / 2);
    //glm::quat pitchQuat(cos(pitchRad / 2), sint * q.x, sint * q.y, sint * q.z);

    //modelingMatrix = matT * glm::toMat4(rollQuat) * glm::toMat4(pitchQuat) * modelingMatrix;
    // first roll, then pitch!

    //modelingMatrix = matT * modelingMatrix;

    //modelingMatrix = glm::mat4(1);


	glUseProgram(gProgram[activeProgramIndex]);
	glUniformMatrix4fv(projectionMatrixLoc[activeProgramIndex], 1, GL_FALSE, glm::value_ptr(projectionMatrix));
	glUniformMatrix4fv(viewingMatrixLoc[activeProgramIndex], 1, GL_FALSE, glm::value_ptr(viewingMatrix));
	glUniformMatrix4fv(modelingMatrixLoc[activeProgramIndex], 1, GL_FALSE, glm::value_ptr(modelingMatrix));
	glUniform3fv(eyePosLoc[activeProgramIndex], 1, glm::value_ptr(eyePos));

	// Draw the scene
    drawModel();

	angle += 0.5;
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
    else if (key == GLFW_KEY_G && action == GLFW_PRESS)
    {
        //glShadeModel(GL_SMOOTH);
        activeProgramIndex = 0;
    }
    else if (key == GLFW_KEY_P && action == GLFW_PRESS)
    {
        //glShadeModel(GL_SMOOTH);
        activeProgramIndex = 1;
    }
    else if (key == GLFW_KEY_F && action == GLFW_PRESS)
    {
        //glShadeModel(GL_FLAT);
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

    //glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
    //glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
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
