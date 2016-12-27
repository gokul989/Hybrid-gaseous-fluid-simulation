#include <stdio.h>
#include <stdlib.h>
#include "glut.h"
#include <GL/gl.h>
#include <GL/glu.h>
#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>
#include "imageloader.h"

#define N 32
#define asize (N+2) * (N+2) *(N+2)
#define swap(x1,x2) {float *temp=x1;x1=x2;x2=temp;}
inline int grid(int i, int j, int k){
  return (i+j*(N+2)+k*(N+2)*(N+2));
}
using namespace std;
float rotateangle = 0.0f;
const float PI = 3.14f;
float gravity = 3.0f;
const int particle_count = 250;
bool visualizeTemp = false;
float timestep = 0.01f;
float aaa = 0.0f;
float particleSize = 2.5*(1.0f/N);
float dt=0.05f,density=3.0f, force = 3.30f;
int winWidth, winHeight;

float u[asize], v[asize],w[asize], u_prev[asize], v_prev[asize],w_prev[asize];
float densityField[asize], dens_prev[asize];

GLuint textureId;
float threshold = 50*(particleSize);

struct vec3 {
	float x;
	float y;
	float z;
};
struct vec2 {
	float x;
	float y;
};

float randomFloat() {
	return (float)rand() / ((float)RAND_MAX + 1);
}

struct particle {
	vec3 pos;
	vec3 velocity;
	vec3 color;
	vec2 acceleration;
	vec2 pressure;
	vec2 force;
	float density;
	float mass;
	float timeAlive;
	float lifespan;
	
};
//TBD
vector<int> findNeighbors(particle *prt,particle* pp);


class ParticleSystem {
public:	
		GLuint textureId;
		particle particles[particle_count];
		
		float timebwSteps;
		float angle;			
		
		void createParticle(particle* p) {
			p->pos.x = 0.30f+0.1*randomFloat();
			p->pos.y = 0.0f;
			p->pos.z = 0.5f+0.1*randomFloat();//0.30f+0.1*randomFloat()+timestep;
			float q= 0.5f * randomFloat() - 0.25f;
			p->density = 0.58f;
			p->mass = 0.005f;
			p->velocity.x =  1.0f*cos(angle) + q;
			p->velocity.y = 2.0f + q;
			p->velocity.z = 1.0f*sin(angle) + q;
			p->color.x = p->color.y = p->color.z = 0.5f;
			p->timeAlive = 0;
			p->lifespan = randomFloat()+1 ;
		}
		
		
		void step() {			
			
			angle +=  aaa + timestep;
			while (angle > 2 * PI) {
				angle = angle -(2 * PI);
			}
			                                                                      
			for(int i = 0; i < particle_count; i++) {
				particle* p = &particles[i];
				//vector<int> neigh = findNeighbors(particles,p);
				//if(neigh.size()<20)
				//sph() //WIP-TBD
				//cout<<"neighbor size...."<<neigh.size()<<endl;
				p->pos.x += p->velocity.x * timestep;
				p->pos.y += p->velocity.y * timestep;
				p->pos.z += p->velocity.z * timestep;
				p->velocity.y += (gravity * timestep);
				p->timeAlive += timestep;
				
				if ( p->pos.y>(0.25f- 0.2*randomFloat())) {
					
					int a,b,c;
					bool skip = false;
					if(p->pos.x >= 1.0f)
						skip = true;
					if(p->pos.x <=0.0f)
						skip = true;
				    if(p->pos.y >= 1.0f)
						skip = true;
					if(p->pos.y <=0.0f)
						skip = true;
					if(p->pos.z >= 1.0f)
						skip = true;
					if(p->pos.z <=0.0f)
						skip = true;
					if(skip == false)
					{
					a = (int)(p->pos.x*(N+2));
					 b = (int)((p->pos.y)*(N+2));
					 c = (int)((p->pos.z)*(N+2));
					densityField[grid(a,b,c)]=density;					
					v[grid(a,b,c)] = force;
					}
					createParticle(p);
					


				}
			}
		}
	
		ParticleSystem(GLuint textureId1) {
			textureId = textureId1;
			timebwSteps = 0;
			
			angle = 0;
			for(int i = 0; i < particle_count; i++) {
				createParticle(&particles[i]);
			}
			for(int i = 0; i < 5 / timestep; i++) {
				step();
			}
		}
		
		
		void nextFrame(float dt) {
			while (dt > 0) {
				if (timebwSteps < dt) {
					dt -= timebwSteps;
					step();
					timebwSteps = timestep;
				}
				else {
					timebwSteps -= dt;
					dt = 0;
				}
			}
		}		
		
	
};

ParticleSystem* particlesystem;
float distance1(particle* p1,particle* p2)
{
	return abs(p1->pos.x - p2->pos.x)+abs(p1->pos.y - p2->pos.y);
}


vector<int> findNeighbors(particle *prt,particle* pp)
{
	int i,j;
	vector<int> list;
	for(i=0;i<particle_count;i++)
	{
		//float d = distance1(prt,pp);
			if(distance1(prt,pp) <=threshold)
				list.push_back(i);

	}
	return list;
}


void drawDensity()
{
	int i, j;
	float x, y,z, c1, c2, c3, c4,c5, c6, c7, c8;
	float h = 1.0f/N;

	glBegin(GL_QUADS);
	for(i = 0; i <= N; i++)
	{
		x = (i-0.5f)*h;
		for(j = 0; j <= N; j++)
		{
			y = (j-0.5f)*h;
			for(int k=0;k<=N;k++)
			{
			z = (k-0.5f)*h;
			c1 = densityField[grid(i,j,k)]*0.25;
			c2 = densityField[grid(i+1,j,k)]*0.25;
			c3 = densityField[grid(i+1,j+1,k)]*0.25;
			c4 = densityField[grid(i,j+1,k)]*0.25;
			if(c1==0 || c2 ==0 || c3 ==0 || c4 == 0)
				continue;
			//cout<<c1<<endl<<c2<<endl;
			
				c1 = densityField[grid(i,j,k)]*0.25f;
				c2 = densityField[grid(i,j+1,k)]*0.25f;
				c3 = densityField[grid(i+1,j,k)]*0.25f;
				c4 = densityField[grid(i+1,j+1,k)]*0.25f;

				c5 = densityField[grid(i,j,k+1)]*0.25f;
				c6 = densityField[grid(i,j+1,k+1)]*0.25f;
				c7 = densityField[grid(i+1,j,k+1)]*0.25f;
				c8 = densityField[grid(i+1,j+1,k+1)]*0.25f;				
                
				
				float alpha = 0.04f;
				glColor4f ( c8, c8, c8, alpha );
				glVertex3f ( x+h,y+h,z+h );
				glColor4f ( c6, c6, c6, alpha );
				glVertex3f ( x, y+h, z+h);
				glColor4f ( c5, c5, c5, alpha );
				glVertex3f ( x, y, z+h );
				glColor4f ( c7, c7, c7, alpha );
				glVertex3f ( x+h, y, z+h );

				glColor4f ( c4, c4, c4, alpha );
				glVertex3f ( x+h, y+h, z );
				glColor4f ( c8, c8, c8, alpha );
				glVertex3f ( x+h,y+h,z+h );
				glColor4f ( c7, c7, c7, alpha );
				glVertex3f ( x+h, y, z+h );
				glColor4f ( c3, c3, c3, alpha );
				glVertex3f ( x+h, y, z );

				glColor4f ( c2, c2, c2, alpha );
				glVertex3f ( x, y+h, z );
				glColor4f ( c4, c4, c4, alpha );
				glVertex3f ( x+h, y+h, z );
				glColor4f ( c3, c3, c3, alpha );
				glVertex3f ( x+h, y, z );
				glColor4f ( c1, c1, c1, alpha );
				glVertex3f ( x, y, z );
                
				glColor4f ( c6, c6, c6, alpha );
				glVertex3f ( x, y+h, z+h);
				glColor4f ( c2, c2, c2, alpha );
				glVertex3f ( x, y+h, z );
				glColor4f ( c1, c1, c1, alpha );
				glVertex3f ( x, y, z );
				glColor4f ( c5, c5, c5, alpha );
				glVertex3f ( x, y, z+h );

				glColor4f ( c3, c3, c3, alpha );
				glVertex3f ( x+h, y, z );
				glColor4f ( c1, c1, c1, alpha );
				glVertex3f ( x, y, z );
				glColor4f ( c5, c5, c5, alpha );
				glVertex3f ( x, y, z+h );
				glColor4f ( c7, c7, c7, alpha );
				glVertex3f ( x+h, y, z+h );

				glColor4f ( c4, c4, c4, alpha );
				glVertex3f ( x+h, y+h, z );
				glColor4f ( c2, c2, c2, alpha );
				glVertex3f ( x, y+h, z );
				glColor4f ( c6, c6, c6, alpha );
				glVertex3f ( x, y+h, z+h);
				glColor4f ( c8, c8, c8, alpha );
				glVertex3f ( x+h, y+h, z+h );	
			}
		}
	}
	glEnd();
}
bool showVelocity = false;


/* void displayLines()
{
	
	int i, j;
	float x, y, c1, c2, c3, c4;
	float h = 1.0f/N;


	glBegin(GL_LINES);
	for(i = 0; i <= N; i+=2)
	{
		x = (i)*h;
		for(j = 0; j <= N; j+=2)
		{
			y = (j)*h;
			c1 = densityField[grid(i,j)]*0.25;
			c2 = densityField[grid(i+1,j)]*0.25;
			c3 = densityField[grid(i+1,j+1)]*0.25;
			c4 = densityField[grid(i,j+1)]*0.25;
			if(c1==0 || c2 ==0 || c3 ==0 || c4 == 0)
				continue;
			float vx = u[grid(i,j)];
			float vy = v[grid(i,j)];
			float color = 100*(abs(vy)-y/2);
			float green = y;
			//cout<<color<<endl;
			if(vx>0)
			{
				if(vy>0)
				{
				glColor4f(color, green, 1-color,1.0f); 
				glVertex2f(x, y);

				//glColor4f(c2, c2, c2,1.0f); 
				glVertex2f(x+h, y+h);
				}

				else
				{
				glColor4f(color, green, 1-color,1.0f); 
				glVertex2f(x, y+h);

				//glColor4f(c2, c2, c2,1.0f); 
				glVertex2f(x+h, y);

				}


			}
			else
			{
				if(vy>0)
				{
					glColor4f(color, green, 1-color,1.0f);  
					glVertex2f(x+h, y);
				

				//glColor4f(c2, c2, c2,1.0f); 
				glVertex2f(x, y+h);


				}

				else
				{
				glColor4f(color, green, 1-color,1.0f); 
				glVertex2f(x+h, y+h);

				//glColor4f(c2, c2, c2,1.0f); 
				glVertex2f(x, y);
				
				}
			}
				
			
		}
	}
	glEnd();


} */

void displaybox()
{
	int i, j;
	float x, y, c1, c2, c3, c4;
	float h = 1.0f/N;


	glBegin(GL_LINES);
	
				glColor4f(1.0, 1.0, 0.0,1.0f); 
				glVertex3f(0.0, 0.0,0.0);

				glColor4f(1.0, 1.0, 0.0,1.0f); 
				glVertex3f(1.0, 0.0,0.0);

				glColor4f(1.0, 1.0, 0.0,1.0f);  
				glVertex3f(1.0, 0.0,0.0);

				glColor4f(1.0, 1.0, 0.0,1.0f); 
				glVertex3f(1.0, 1.0,0.0);

				glColor4f(1.0, 1.0, 0.0,1.0f);  
				glVertex3f(1.0, 1.0,0.0);

				glColor4f(1.0, 1.0, 0.0,1.0f);  
				glVertex3f(0.0, 1.0,0.0);

				glColor4f(1.0, 1.0, 0.0,1.0f);  
				glVertex3f(0.0, 1.0,0.0);

				glColor4f(1.0, 1.0, 0.0,1.0f);  
				glVertex3f(0.0, 0.0,0.0);



				//
				glColor4f(1.0, 1.0, 0.0,1.0f); 
				glVertex3f(0.0, 0.0,1.0);

				glColor4f(1.0, 1.0, 0.0,1.0f);  
				glVertex3f(1.0, 0.0,1.0);

				glColor4f(1.0, 1.0, 0.0,1.0f); 
				glVertex3f(1.0, 0.0,1.0);

				glColor4f(1.0, 1.0, 0.0,1.0f); 
				glVertex3f(1.0, 1.0,1.0);

				glColor4f(1.0, 1.0, 0.0,1.0f); 
				glVertex3f(1.0, 1.0,1.0);

				glColor4f(1.0, 1.0, 0.0,1.0f);  
				glVertex3f(0.0, 1.0,1.0);

				glColor4f(1.0, 1.0, 0.0,1.0f); 
				glVertex3f(0.0, 1.0,1.0);

				glColor4f(1.0, 1.0, 0.0,1.0f); 
				glVertex3f(0.0, 0.0,1.0);
				
				//
				glColor4f(1.0, 1.0, 0.0,1.0f);  
				glVertex3f(0.0, 0.0,0.0);

				glColor4f(1.0, 1.0, 0.0,1.0f);  
				glVertex3f(1.0, 0.0,0.0);

				glColor4f(1.0, 1.0, 0.0,1.0f);  
				glVertex3f(1.0, 0.0,0.0);

				glColor4f(1.0, 1.0, 0.0,1.0f);  
				glVertex3f(1.0, 0.0,1.0);

				glColor4f(1.0, 1.0, 0.0,1.0f); 
				glVertex3f(1.0, 0.0,1.0);

				glColor4f(1.0, 1.0, 0.0,1.0f); 
				glVertex3f(0.0, 0.0,1.0);

				glColor4f(1.0, 1.0, 0.0,1.0f);  
				glVertex3f(0.0, 0.0,1.0);

				glColor4f(1.0, 1.0, 0.0,1.0f);  
				glVertex3f(0.0, 0.0,0.0);

				//
				glColor4f(1.0, 1.0, 0.0,1.0f); 
				glVertex3f(0.0, 1.0,0.0);

				glColor4f(1.0, 1.0, 0.0,1.0f); 
				glVertex3f(1.0, 1.0,0.0);

				glColor4f(1.0, 1.0, 0.0,1.0f); 
				glVertex3f(1.0, 1.0,0.0);

				glColor4f(1.0, 1.0, 0.0,1.0f); 
				glVertex3f(1.0, 1.0,1.0);

				glColor4f(1.0, 1.0, 0.0,1.0f); 
				glVertex3f(1.0, 1.0,1.0);

				glColor4f(1.0, 1.0, 0.0,1.0f); 
				glVertex3f(0.0, 1.0,1.0);

				glColor4f(1.0, 1.0, 0.0,1.0f);  
				glVertex3f(0.0, 1.0,1.0);

				glColor4f(1.0, 1.0, 0.0,1.0f); 
				glVertex3f(0.0, 1.0,0.0);


	glEnd();



}

	void draw(ParticleSystem* pps) {
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE);
    
		glEnable(GL_ALPHA_TEST);
		glAlphaFunc(GL_GREATER, 0);
		drawDensity();
		
		glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
			vector<particle*> ps;
			for(int i = 0; i < particle_count; i++) {
				ps.push_back(&pps->particles[i]);
			}
			//cout<<pps->textureId<<endl;
			glEnable(GL_TEXTURE_2D);
			glBindTexture(GL_TEXTURE_2D, pps->textureId);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
			glBegin(GL_QUADS);
			for(unsigned int i = 0; i < ps.size(); i++) {
				particle* p = ps[i];
				/*if(visualizeTemp)
				glColor4f(1 - p->timeAlive / p->lifespan, 0.0f, p->timeAlive / p->lifespan,
						  (1 - p->timeAlive / p->lifespan));
				else*/
					glColor4f(0.50f, 0.50f, 0.50f,
						  (1 -((p->pos.y)/0.45f)));//* p->timeAlive / p->lifespan));
				float psize = particleSize / 2;
				vec3 pos = p->pos;
				//cout<<"position is ..."<<pos.x<<" "<<pos.y<<"  "<<pos.z<<endl;
				glTexCoord2f(0, 0);
				//glColor4f(0.1f,0.1f,0.1f,0.8f);
				glVertex3f(pos.x - psize, pos.y - psize,pos.z - psize);
				glTexCoord2f(0, 1);
				//glColor4f(0.9f,0.9f,0.9f,0.8f);
				glVertex3f(pos.x - psize, pos.y + psize,pos.z - psize);
				glTexCoord2f(1, 1);
				//glColor4f(0.8f,0.8f,0.8f,0.3f);
				glVertex3f(pos.x + psize, pos.y + psize,pos.z - psize);
				glTexCoord2f(1, 0);
				//glColor4f(0.2f,0.2f,0.2f,0.2f);
				glVertex3f(pos.x + psize, pos.y - psize,pos.z - psize);
			}
			glEnd(); 
			glDisable(GL_TEXTURE_2D);
			glDisable(GL_BLEND);

			displaybox();
		}



// standard texture loading routines from openGL tutorial
char* addalphaImage(Image* image, Image* alphaImage) {
	char* pixels = new char[image->width * image->height * 4];
	for(int y = 0; y < image->height; y++) {
		for(int x = 0; x < image->width; x++) {
			for(int j = 0; j < 3; j++) {
				pixels[4 * (y * image->width + x) + j] =
					image->pixels[3 * (y * image->width + x) + j];
			}
			pixels[4 * (y * image->width + x) + 3] =
				alphaImage->pixels[3 * (y * image->width + x)];
		}
	}
	
	return pixels;
}


GLuint loadAlphaTexture(Image* image, Image* alphaImage) {
	char* pixels = addalphaImage(image, alphaImage);
	
	GLuint textureId;
	glGenTextures(1, &textureId);
	glBindTexture(GL_TEXTURE_2D, textureId);
	glTexImage2D(GL_TEXTURE_2D,
				 0,
				 GL_RGBA,
				 image->width, image->height,
				 0,
				 GL_RGBA,
				 GL_UNSIGNED_BYTE,
				 pixels);
	
	delete pixels;
	return textureId;
}


void	display(void)
{
	glViewport ( 0, 0, winWidth, winHeight);
	glMatrixMode ( GL_PROJECTION );
	glLoadIdentity ();	
	//gluOrtho2D ( 0.0, 1.0, 0.0, 1.0 );
	gluPerspective(30,winWidth/winHeight,1,1000);
	gluLookAt(-0.5,3.5,4,0,1.0,0.8,0,1,0);
	//glOrtho(0.0,1.0,0.0,1.0,1.0,1000.0);
	//gluLookAt(0,1,0,0,1,0,0,1,0);
	float tx=0.220f,ty=0.0f,tz=0.42f;
	//glTranslatef(tx*-1,0.00f,tz*-1);
	glRotatef(rotateangle,0.0f,1.0f,0.0f);
	//glTranslatef(tx,-0.00f,tz);
	glClearColor( 0.0f, 0.0f, 0.00f, 1.0f );
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	draw(particlesystem);
	glutSwapBuffers();
	particlesystem->nextFrame(10/3000.0f);
	glutPostRedisplay();

    //glutSwapBuffers();  // 
}

void	resize(int w,int h)
{
    glViewport(0,0,w,h);
    winWidth = w;
    winHeight = h; 
}


void	keyboard(unsigned char key, int x, int y)
{
	int i = (int)(N/2);
	int j = (int)(N/2);
	int top = N-1;
	int bottom = 2;
    switch(key) {
	case 'z':
			aaa+=0.001f;
			visualizeTemp = true;
			break;
	case 'x':
			aaa-=0.001f;
			visualizeTemp = false;
			break;
	case 'r':
			rotateangle+=1.0f;
			break;
	case 't':
			rotateangle-=1.0f;
			break;
   
	
	
    default:
		break;
    }
}


//grid equations and routines based on Stam's stable fluids paper and notes.
void addDensity(float * x, float * s, float dt)
{
	for(int i = 0; i < asize; i++)
		x[i] += dt*s[i];
}

void computeDiffuse(int b, float * x, float * x0,  float dt)
{
	int i, j, k;
		
	for(k = 0; k < 30; k++) 
	{
		for(i = 1; i <= N; i++)
		{
			for(j = 1; j <= N; j++)
			{
				for(k=1;k<=N;k++)
				{
				// *((x+i*(N+2)) + j) =  *((x0+i*(N+2)) + j);
				x[grid(i,j,k)] = x0[grid(i,j,k)];
				}
			}
		}
		
	}
}
 

void computeAdvection(int b, float * d, float * d0, float * u, float * v,float *w, float dt)
{
	int i, j,k, i0, j0,k0, i1, j1,k1;
	float x, y,z, s0, t0,r0, s1, t1,r1, dt0;

	dt0 = dt*N;
	for(i = 1; i <= N; i++)
	{
		for(j = 1; j <= N; j++)
		{
			for(k=1;k<=N;k++)
			{
			//x = i - dt0 * *((u+i*(N+2)) + j);
			//y = j - dt0 * *((v+i*(N+2)) + j);	
			x = i - dt0 * u[grid(i,j,k)];
			y = j - dt0 * v[grid(i,j,k)];
			z = k - dt0 * w[grid(i,j,k)];
			
			i0 = (int)x; i1 = i0+1;
			j0 = (int)y; j1 = j0+1;
			k0 = (int)z; k1 = k0+1;

			s1 = x - i0; s0 = 1 - s1;
			t1 = y - j0; t0 = 1 - t1;
			r1 = z - k0; r0 = 1 - r1;

			float q1,q2;
			//*((d+i*(N+2)) + j) = s0 * (t0 * *((d0+i0*(N+2)) + j0) + t1 * *((d0+i0*(N+2)) + j1)) + s1 * (t0 * *((d0+i1*(N+2)) + j0) + t1 * *((d0+i1*(N+2)) + j1));
			q1 = s0*(t0*d0[grid(i0,j0,k0)] + t1*d0[grid(i0,j1,k0)]) + s1*(t0*d0[grid(i1,j0,k0)] + t1*(d0[grid(i1,j1,k0)]));
			q2 = s0*(t0*d0[grid(i0,j0,k1)] + t1*d0[grid(i0,j1,k1)]) + s1*(t0*d0[grid(i1,j0,k1)] + t1*(d0[grid(i1,j1,k1)]));
			d[grid(i,j,k)] = r0*q1 + r1*q2;
		}
		}
	}
	
}


void computeProject(float * u, float * v,float *w, float * p, float * div)
{
	int i, j, k,f;
	float h = 1.0/N;

	for(i = 1; i <= N; i++)
	{
		for(j = 1; j <= N; j++)
		{
			for(k=0;k<=N;k++)
			{
			//*((div+i*(N+2)) + j) = -0.5 * h * (*((u+i*(N+2)) + j) - *((u+(i-1)*(N+2)) + j) + *((v+i*(N+2)) + j+1) - *((v+i*(N+2)) + j-1));
			//*((p+i*(N+2)) + j) = 0;
			div[grid(i,j,k)] = -0.33 * h * (u[grid(i+1,j,k)] - u[grid(i-1,j,k)] + v[grid(i,j+1,k)] - v[grid(i,j-1,k)]+w[grid(i,j,k+1)] - w[grid(i,k,k-1)]);
			p[grid(i,j,k)] = 0;
			}
		}
	}	

	for(f = 0; f < 20; f++)
	{
		for(i = 1; i <= N; i++)
		{
			for(j = 1; j <= N; j++)
			{
				for(k=0;k<=N;k++)
			{
				//*((p+i*(N+2)) + j) = (*((div+i*(N+2)) + j) + *((p+(i-1)*(N+2)) + j) + *((p+(i+1)*(N+2)) + j) + *((p+i*(N+2)) + j-1) + *((p+i*(N+2)) + j+1))/4;
				p[grid(i,j,k)] = (div[grid(i,j,k)] + p[grid(i-1,j,k)] + p[grid(i+1,j,k)] + p[grid(i,j-1,k)] + p[grid(i,j+1,k)]+p[grid(i,j,k-1)]+p[grid(i,j,k+1)])/6;
				}
			}
		}
	}

	for(i = 1; i <= N; i++)
	{
		for(j = 1; j <= N; j++)
		{
			for(k=0;k<=N;k++)
			{
			//*((u+i*(N+2)) + j) -= 0.5 * (*((p+(i+1)*(N+2)) + j) - *((p+(i-1)*(N+2)) + j))/h;
			//*((v+i*(N+2)) + j) -= 0.5 * (*((p+(i)*(N+2)) + j+1) - *((p+(i)*(N+2)) + j))/h;
			u[grid(i,j,k)] -= 0.33 * (p[grid(i+1,j,k)] - p[grid(i-1,j,k)])/h;
			v[grid(i,j,k)] -= 0.33 * (p[grid(i,j+1,k)] - p[grid(i,j-1,k)])/h;
			w[grid(i,j,k)] -= 0.33 * (p[grid(i,j,k+1)] - p[grid(i,j,k-1)])/h;
			}
		}
	}
	
}


void velocityUpdate(float * u, float * v,float *w, float * u0, float * v0,float *w0,  float dt)
{
	addDensity( u, u0, dt); 
	addDensity( v, v0, dt);
	addDensity( w, w0, dt);
	swap(u0, u); 
	computeDiffuse(1, u, u0, dt);
	swap(v0, v); 
	computeDiffuse(2, v, v0, dt);
	swap(w0, w); 
	computeDiffuse(3, w, w0, dt);

	computeProject(u, v,w, u0, v0);
	swap(u0, u);
	swap(v0, v);
	swap(w0,w);
	computeAdvection(1, u, u0, u0, v0,w0, dt);
	computeAdvection(2, v, v0, u0, v0,w0, dt);
	computeAdvection(3, w, w0, u0, v0,w0, dt);
	computeProject(u, v,w, u0, v0);			
}


void densityUpdate( float * x, float * x0, float * u, float * v,float *w,  float dt)
{
	addDensity(x, x0, dt);
	swap (x0, x); computeDiffuse (0, x, x0, dt);
	swap (x0, x); computeAdvection (0, x, x0, u, v,w, dt);
}


void	gridUpdate()
{
	//draw(particlesystem);
	//drawDensity();
	//glutSwapBuffers();
	//particlesystem->nextFrame(10/3000.0f);
	for ( int i=0 ; i<N+2 ; i++ ) 
		for(int j=0;j<N+2;j++)
			for(int k=0;k<N+2;k++)
		u_prev[grid(i,j,k)] = v_prev[grid(i,j,k)]=w_prev[grid(i,j,k)] = dens_prev[grid(i,j,k)] = 0.0f;
	velocityUpdate((float *) u,(float *) v,(float *) w,(float *) u_prev,(float *) v_prev,(float *) w_prev, dt);
	densityUpdate((float *) densityField,(float *) dens_prev,(float *) u,(float *) v, (float *) w,dt);
	glutPostRedisplay ();
}


int main(int argc, char* argv[])
{    
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
  	glEnable(GL_DEPTH_TEST);
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
	glutInitWindowSize(600,600);											
    glutCreateWindow("Graphicsfinaldemo");	
	Image* image = loadBMP("circle.bmp");
	Image* alphaImage = loadBMP("circlealpha.bmp");
	textureId = loadAlphaTexture(image, alphaImage);
	cout<<textureId<<endl;
    glutDisplayFunc(display);
    glutReshapeFunc(resize);
    glutKeyboardFunc(keyboard);
	glutIdleFunc(gridUpdate);
	particlesystem = new ParticleSystem(textureId);
	for ( int i=0 ; i<N+2 ; i++ ) 
		for(int j=0;j<N+2;j++)
			for(int k=0;k<N+2;k++)
		u_prev[grid(i,j,k)] = v_prev[grid(i,j,k)] = dens_prev[grid(i,j,k)]  =u[grid(i,j,k)]=v[grid(i,j,k)]=densityField[grid(i,j,k)]= 0.0f;
	
    glutMainLoop();
    return 0;        
}