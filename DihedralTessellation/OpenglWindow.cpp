#include"tilingOpt.h"

using namespace std;
namespace Tiling_tiles{

	GLdouble cont[300][3];
	GLdouble cont2[300][3];
	int Size = 0;
	int Size2 = 0;
	int ww = 1200;
	int wh = 800;
	GLdouble quad[12][3] = {
		{ 300, 330, 0 }, { 309, 300, 0 }, { 280, 320, 0 }, { 2, 3, 0 },
		{ -1, 2, 0 }, { 0, 0.5, 0 }, { 1, 1, 0 }, { 1, 2, 0 },
		{ -0.5, 1, 0 }, { -0.5, 2, 0 }, { 0.5, 2, 0 }, { 0.5, 1, 0 } };
	GLfloat size = 3.0;
	//int n;
	//
	
	int points_cv2gl(vector<Point2f> input, double z, GLdouble contour[300][3])
	{		
		int insize = input.size();
		for (int i = 0; i < insize; i++)
		{
			contour[i][0] = input[i].x;
			contour[i][1] = input[i].y;
			contour[i][2] = z;
		}	
		cout << "tttt" << endl;
		//cout << "t : " << coords[0] << " " << coords[1] << " "<<coords[2] << endl;
		return insize;
	}

	void myIdle()
	{
		glutPostRedisplay();
	}
	////------------------------------------------------------------  OnDraw()  
	////  
	void CALLBACK vertexCallback(GLvoid *vertex)
	{
		const GLdouble *pointer = (GLdouble *)vertex;
		//glColor3dv(pointer + 3);//在此设置颜色  
		glColor3d(0, 0, 0);
		glVertex3dv(pointer);
	}
	void CALLBACK beginCallback(GLenum which)
	{
		glBegin(which);
	}
	void CALLBACK endCallback()
	{
		glEnd();
	}
	void CALLBACK errorCallback(GLenum errorCode)
	{
		const GLubyte *estring;
		estring = gluErrorString(errorCode);
		fprintf(stderr, "Tessellation Error: %s\n", estring);
		exit(0);
	}
	void CALLBACK combineCallback(GLdouble coords[3],
		GLdouble *vertex_data[4],
		GLfloat weight[4], GLdouble **dataOut)
	{
		GLdouble *vertex;
		int i;
		vertex = (GLdouble *)malloc(6 * sizeof(GLdouble));
		vertex[0] = coords[0];
		vertex[1] = coords[1];
		vertex[2] = coords[2];
		for (i = 3; i < 7; i++)
		{
			vertex[i] = weight[0] * vertex_data[0][i]
				+ weight[1] * vertex_data[1][i]
				+ weight[2] * vertex_data[2][i]
				+ weight[3] * vertex_data[3][i];
		}
		*dataOut = vertex;
	}
	void OnDraw()
	{
		// clear the screen & depth buffer  
		glClear(GL_COLOR_BUFFER_BIT);
		// clear the previous transform  
		glLoadIdentity();

		GLUtesselator *tobj = gluNewTess();
		if (!tobj) { return; }

		gluTessCallback(tobj, GLU_TESS_VERTEX, (void (CALLBACK *)())vertexCallback);
		gluTessCallback(tobj, GLU_TESS_BEGIN, (void (CALLBACK *)())beginCallback);
		gluTessCallback(tobj, GLU_TESS_END, (void (CALLBACK *)())endCallback);
		gluTessCallback(tobj, GLU_TESS_ERROR, (void (CALLBACK *)())errorCallback);
		gluTessCallback(tobj, GLU_TESS_COMBINE, (void (CALLBACK *)())combineCallback);


		GLUtesselator *tobj2 = gluNewTess();
		if (!tobj2) { return; }

		gluTessCallback(tobj2, GLU_TESS_VERTEX, (void (CALLBACK *)())vertexCallback);
		gluTessCallback(tobj2, GLU_TESS_BEGIN, (void (CALLBACK *)())beginCallback);
		gluTessCallback(tobj2, GLU_TESS_END, (void (CALLBACK *)())endCallback);
		gluTessCallback(tobj2, GLU_TESS_ERROR, (void (CALLBACK *)())errorCallback);
		gluTessCallback(tobj2, GLU_TESS_COMBINE, (void (CALLBACK *)())combineCallback);

		// glShadeModel(GL_FLAT);  

		// gluTessProperty(tobj,GLU_TESS_WINDING_RULE,GLU_TESS_WINDING_POSITIVE); //GLU_TESS_WINDING_ODD  
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glViewport(0, 0, ww/3, wh/2);
		// set up the projection matrix  


		// just use a perspective projection  
		//gluPerspective(45,(float)w/h,0.1,100);  
		glOrtho(0, 600, 0, 600, 0.0, -10.0);

		gluTessBeginPolygon(tobj, NULL);
		gluTessBeginContour(tobj);	
		for (int i = 0; i < Size; i++)
		{
			gluTessVertex(tobj, cont[i], cont[i]);
			//gluTessVertex(tobj, contour[1], contour[1]);
		}
		gluTessEndContour(tobj);
		gluTessEndPolygon(tobj);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();

		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glViewport(0, wh / 2, ww / 3, wh / 2);
		// set up the projection matrix  


		// just use a perspective projection  
		//gluPerspective(45,(float)w/h,0.1,100);  
		glOrtho(0, 600, 0, 600, -10.0, -20.0);
		gluTessBeginPolygon(tobj2, NULL);
		gluTessBeginContour(tobj2);
		for (int i = 0; i < Size2; i++)
		{
			gluTessVertex(tobj2, cont2[i], cont2[i]);
			//gluTessVertex(tobj, contour[1], contour[1]);
		}
		gluTessEndContour(tobj2);

		//gluTessBeginContour(tobj);                      // inner quad (hole)  
		//gluTessVertex(tobj, quad[0], quad[0]);
		//gluTessVertex(tobj, quad[1], quad[1]);
		//gluTessVertex(tobj, quad[2], quad[2]);
		////gluTessVertex(tobj, quad[7], quad[7]);
		//gluTessEndContour(tobj);

		gluTessEndPolygon(tobj2);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();

		gluDeleteTess(tobj);
		gluDeleteTess(tobj2);
		glutSwapBuffers();

	}
	//------------------------------------------------------------  OnInit()  
	//  
	void OnInit()
	{
		glClearColor(1,1,1,0); 
	}
	//------------------------------------------------------------  OnExit()  
	//  
	void OnExit()
	{
	}
	//------------------------------------------------------------  OnReshape()    

	void OnReshape(int w, int h)
	{
		//// prevents division by zero when minimising window 
		if (h == 0)
		{
			h = 1;
		}
		//cout << w << "  " << h << endl;
		ww = w;
		wh = h;
	}

	void myMouse(int button, int state, int x, int y)
	{
		if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN);

		if (button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN)
			exit(0);
	}

	void OpenWindow(int w, int h, vector<Point2f> contour1, vector<Point2f> contour2)
	{
		glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
		// set the initial window size  
		glutInitWindowSize(w, h);
		// create the window  
		glutCreateWindow("Dual Tessellation");
		// run our custom initialisation  
		OnInit();
		// set the function to use to draw our scene  
		//Point2f t(-200, 300);
		//Point2f t1(-200, 0);
		//Point2f t2(200, 0);
		//vector<Point2f> contour1;
		//contour1.push_back(t);
		//contour1.push_back(t1);
		//contour1.push_back(t2);
		Size = points_cv2gl(contour1,0.0, cont);
		Size2 = points_cv2gl(contour2, 15.0,cont2);
		cout << "Size1 : " << Size << "  Size2: "<<Size2 << endl;

		glutDisplayFunc(OnDraw);
		
		// set the function to handle changes in screen size  
		glutReshapeFunc(OnReshape);
		
		//  glutIdleFunc(&myIdle);  
		// set the function to be called when we exit  
		atexit(OnExit);

		// this function runs a while loop to keep the program running.  
		glutMainLoop();
	}

	
}