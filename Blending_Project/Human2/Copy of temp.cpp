/*
cc -o fdm_surface11 fdm_surface11.c -lglut -lGLU -lGL -lXmu -lXext -lX11 -lm
*/

#include <gl/glut.h>
//#include "gl\gl.h"
//#include "gl\glu.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#define jb  10

#define nu 15
#define nv 15
#define nuv 289
#define nuvm 277

#define nn 3
#define icurve 8



int ciggj(double **,int,double [jb]);
int Bsurface(double *, double *, double *, double **,double **);      
int Gellipse(int,int,double,double,double,double *,double *,double *);
int cons(double,double,double,double,double,double,double [6]);
int sincos(double,double,double,double,double,double,double,
	   double,double,double,double,double,double,double csc[11]);

FILE *in1,*in2,*out,*out2,*faceout1;

#ifndef CALLBACK
#define CALLBACK
#endif
#define pi 3.1415926

GLuint startList;

void prod1(int *n,double *uu,double *un) 
{
  int i;
  *un = 1.0;
  for (i = 1; i <= *n; i++)
  *un = *un * *uu;
}


void CALLBACK errorCallback(GLenum errorCode)
{
   const GLubyte *estring;

   estring = gluErrorString(errorCode);
   fprintf(stderr, "Quadric Error: %s\n", estring);
   exit(0);
}

void init(void) 
{
    int i,j,n,n0,l,m,im[100],iu,iv,jj,iuu; 
    double xm[50][200],ym[50][200],zm[50][200],x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
		aa[nn-1][nn-1],bb[nn-1],anx[100][10],any[100][10],anz[100][10],dv,xi,ff,
		aac[nn-1][nn-1],bbc[nn-1],vv0,vv1,nxi[100][10],nyi[100][10],nzi[100][10],
		x0,x00,y0,y00,al,all,bs,bss,fmin,dx0,dy0,ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,
		h0,h01,h02,h1,h11,h12,w1,w2,cx0[6],cy0[6],cz0[6],ai0,ai1,ai2,ai3,ai4,ai5,
		csx[12],csy[12],csz[12],ccx[12],ccy[12],ccz[12],di2,di4,di6,du,u01,u02,u03,
		u04,u05,u11,u12,u13,u14,u15,un,nx1,nx2,nx3,nx4,ny1,ny2,ny3,ny4,nz1,nz2,
		nz3,nz4,nxyz,txi[100],tyi[100],tzi[100],tx1,ty1,tz1,xcentre[50],ycentre[50],
		zcentre[50],angle[50][200],anglem,a0[50],b0[50],a0u,a0b,b0u,b0b,dpi,
		ai[50][50],bi[50][50];

 //  clock_t start, finish, duration;
    GLUquadricObj *qobj;
    GLfloat mat_ambient[] = { 0.8, 0.8, 0.8, 1.0 };
    GLfloat mat_specular[] = { 0.3, 0.3, 0.3, 1.0 };
    GLfloat mat_shininess[] = { 100.0 };
    GLfloat light_position[] = { 1000.0, 70000.0, -100000.0, 0.0 };
    GLfloat model_ambient[] = { 0.9, 0.9, 0.9, 1.0 };
    glClearColor(0.382, 0.382, 0.382, 1.0);

    glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
    glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
    glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
    glLightfv(GL_LIGHT0, GL_POSITION, light_position);
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, model_ambient);

    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_DEPTH_TEST);


    startList = glGenLists(1);
    qobj = gluNewQuadric();
    gluQuadricDrawStyle(qobj, GLU_FILL); /* flat shaded */
    gluQuadricNormals(qobj, GLU_SMOOTH);
    glNewList(startList, GL_COMPILE);

	n=2*jb+1;
	n0=2*jb;
	double *bbx = new double [n];
	double **aax = new double* [n];
	for(i=0; i<n; i++)
    {
		aax[i] = new double [n];
	}
	double *bbx0 = new double [n0];
	double **aax0 = new double* [n0];
	for(i=0; i<n0; i++)
    {
		aax0[i] = new double [n0];
	}

    in1 = fopen("in.dat","rt"); 
	in2 = fopen("in1.dat","rt");
    out = fopen("result.out","wt");
    out2 = fopen("out2.mel","wt");
    faceout1 = fopen("PDE_Surface.mel","wt");

	for(i=0; i <= icurve-1; i++)
	{
		fscanf(in1,"%d\n",&im[i]);
 		fprintf(out2,"curve -d 3\n");
       for(j=0; j <= im[i]-1; j++)
		{
            fscanf(in1,"%le %le %le\n",&xm[i][j],&ym[i][j],&zm[i][j]);
		    fprintf(out2,"-p   %le %le %le\n",xm[i][j],ym[i][j],zm[i][j]);
		}
		fprintf(out2,";\n");
	}

	for(i=0; i <= icurve-1; i++)
	{
		xcentre[i]=0.0;
		ycentre[i]=0.0;
		zcentre[i]=0.0;
        for(j=0; j <= im[i]-1; j++)
		{
            xcentre[i]=xcentre[i]+xm[i][j];
            ycentre[i]=ycentre[i]+ym[i][j];
            zcentre[i]=zcentre[i]+zm[i][j];
		}
        xcentre[i]=xcentre[i]/im[i];
        ycentre[i]=ycentre[i]/im[i];
        zcentre[i]=zcentre[i]/im[i];
 		fprintf(out2,"curve -d 3\n");
        for(j=0; j <= 1; j++)
		{
		    fprintf(out2,"-p   %le %le %le\n",xcentre[i],ycentre[i],zcentre[i]);
		}
		fprintf(out2,";\n");
	}

	for(i=0; i <= icurve-1; i++)
	{
        for(j=0; j <= im[i]-1; j++)
		{
            angle[i][j]=0.0;
			anglem=(zm[i][j]-zcentre[i])/sqrt((xm[i][j]-xcentre[i])*(xm[i][j]-xcentre[i])
				+(ym[i][j]-ycentre[i])*(ym[i][j]-ycentre[i])+(zm[i][j]-zcentre[i])*(zm[i][j]-zcentre[i]));
			anglem=asin(anglem);
			if((zm[i][j]-zcentre[i])>=0.0 && (xm[i][j]-xcentre[i])>=0.0)angle[i][j]=anglem; 
			if((zm[i][j]-zcentre[i])>=0.0 && (xm[i][j]-xcentre[i])< 0.0)angle[i][j]=pi-anglem; 
			if((zm[i][j]-zcentre[i])< 0.0 && (xm[i][j]-xcentre[i])< 0.0)angle[i][j]=pi-anglem; 
			if((zm[i][j]-zcentre[i])< 0.0 && (xm[i][j]-xcentre[i])> 0.0)angle[i][j]=2.0*pi+anglem; 
		}
	}

	for(i=0; i <= icurve-1; i++)
	{
		a0u=0.0;
		a0b=0.0;
		b0u=0.0;
		b0b=0.0;
        for(j=0; j <= im[i]-1; j++)
		{
            a0u=a0u+(xm[i][j]-xcentre[i])*cos(angle[i][j]);
            a0b=a0b+cos(angle[i][j])*cos(angle[i][j]);
            b0u=b0u+(zm[i][j]-zcentre[i])*sin(angle[i][j]);
            b0b=b0b+sin(angle[i][j])*sin(angle[i][j]);
		}
		a0[i]=a0u/a0b;
		b0[i]=b0u/b0b;
	}


	dpi=pi/100.0;
//  generate ellipses with x=xc+a*cos(angle) and z=zc+b*sin(angle)
	for(i=0; i <= icurve-1; i++)
	{
 		fprintf(out2,"curve -d 3\n");
        for(j=0; j <= im[i]-1; j++)
		{
            x1=xcentre[i]+a0[i]*cos(angle[i][j]);
		    y1=ycentre[i];
            z1=zcentre[i]+b0[i]*sin(angle[i][j]);
		    fprintf(out2,"-p    %le %le %le\n",x1,y1,z1);
		}
		fprintf(out2,";\n");
	}


/*
	for(ic=0; ic <= icurve-1; ic++)
	{
// 		fprintf(out2,"curve -d 3\n");
 		fprintf(faceout1,"curve -d 3\n");
		for(i=0; i <= im[ic]-1; i++)
		{
			x1=xcentre[ic]+ai[ic][0];
			y1=ycentre[ic];
			z1=zcentre[ic]+bi[ic][0];
			for(int j=1;j<=jb;j++)
			{        

                x1=x1+ai[ic][2*j-1]*cos(j*angle[ic][i])
					+ai[ic][2*j]*sin(j*angle[ic][i]);
				y1=y1+0.0;
				z1=z1+bi[ic][2*j-1]*sin(j*angle[ic][i])
					+bi[ic][2*j]*cos(j*angle[ic][i]);
			}
//		    fprintf(out2,"-p    %le %le %le\n",x1,y1,z1);
		    fprintf(faceout1,"-p    %le %le %le\n",x1,y1,z1);
		}
//		fprintf(out2,";\n");
		fprintf(faceout1,";\n");
	}

    int ns=50;
	double anglei;
	du=1.0/ns;
	dv=2.0*pi/ns;
	for(ic=0; ic <= icurve-1; ic++)
	{
 		fprintf(faceout1,"curve -d 3\n");
		for(i=0; i <= ns; i++)
		{
			anglei=i*dv;
			x1=xcentre[ic]+ai[ic][0];
			y1=ycentre[ic];
			z1=zcentre[ic]+bi[ic][0];
			for(int j=1;j<=jb;j++)
			{        

                x1=x1+ai[ic][2*j-1]*cos(j*anglei)
					+ai[ic][2*j]*sin(j*anglei);
				y1=y1+0.0;
				z1=z1+bi[ic][2*j-1]*sin(j*anglei)
					+bi[ic][2*j]*cos(j*anglei);
			}
		    fprintf(faceout1,"-p    %le %le %le\n",x1,y1,z1);
		}
		fprintf(faceout1,";\n");
			if(ic==0)printf("%d %le\n",ic,y1);
	}
*/
//	double xcentrep[4],ycentrep[4],zcentrep[4],aip[4][50],bip[4][50];
	int ic;
	double *xcentrep = new double [4];
	double *ycentrep = new double [4];
	double *zcentrep = new double [4];
	double **aip = new double* [4];
	double **bip = new double* [4];
	for(i=0; i<4; i++)
    {
		aip[i] = new double [50];
		bip[i] = new double [50];
	}

	for(ic=0;ic<=3;ic++)
	{
		xcentrep[ic]=xcentre[ic];
		ycentrep[ic]=ycentre[ic];
		zcentrep[ic]=zcentre[ic];
		for(j=0;j<n;j++)
		{
			aip[ic][j]=ai[ic][j];
			bip[ic][j]=bi[ic][j];
		}
	}
	Bsurface(xcentrep, ycentrep, zcentrep, aip,bip); 

	int ic0=2;
	for(ic=ic0;ic<=ic0+3;ic++)
	{
		xcentrep[ic-ic0]=xcentre[ic];
		ycentrep[ic-ic0]=ycentre[ic];
		zcentrep[ic-ic0]=zcentre[ic];
		for(j=0;j<n;j++)
		{
			aip[ic-ic0][j]=ai[ic][j];
			bip[ic-ic0][j]=bi[ic][j];
		}
	}
	Bsurface(xcentrep, ycentrep, zcentrep, aip,bip);      

/*
	for(ic=0; ic <= icurve-1; ic++)
	{
		for( j=0;j<n;j++)
		{ 
			bbx[j]=0.0;
			for( k=0;k<n;k++)
			{ 
				aax[j][k]=0.0;
			}
		}
		for(i=0; i <= im[ic]-1; i++)
		{
			aax[0][0]=aax[0][0]+1.0;
			aax[0][1]=aax[0][1]+cos(angle[ic][i]);
			aax[0][2]=aax[0][2]+sin(angle[ic][i]);
			bbx[0]=bbx[0]+(xm[ic][i]-xcentre[ic]);
			aax[1][0]=aax[1][0]+cos(angle[ic][i]);
			aax[1][1]=aax[1][1]+cos(angle[ic][i])*cos(angle[ic][i]);
			aax[1][2]=aax[1][2]+sin(angle[ic][i])*cos(angle[ic][i]);
			bbx[1]=bbx[1]+(xm[ic][i]-xcentre[ic])*cos(angle[ic][i]);
			aax[2][0]=aax[2][0]+sin(angle[ic][i]);
			aax[2][1]=aax[2][1]+cos(angle[ic][i])*sin(angle[ic][i]);
			aax[2][2]=aax[2][2]+sin(angle[ic][i])*sin(angle[ic][i]);
			bbx[2]=bbx[2]+(xm[ic][i]-xcentre[ic])*sin(angle[ic][i]);
//				fprintf(out,"ic,k,i,bbx[2]= %d %d %d %le\n",ic,k,i,bbx[2]);
		}
		for(i=0;i<n;i++)
		{
			if(ic==0)fprintf(out,"x--hand %d %le \n",i,bbx[i]);
			for(j=0;j<n;j++)
			{
				if(ic==0)fprintf(out,"x--hand      %d %d %le \n",i,j,aax[i][j]);
			}
		}
		ciggj(aax,n,bbx);
		for( k=0;k<=n-1;k++)ai[ic][k]=bbx[k];
	}

	for(ic=0; ic <= icurve-1; ic++)
	{
		for( j=0;j<n;j++)
		{ 
			bbx[j]=0.0;
			for(k=0;k<n;k++)
			{ 
				aax[j][k]=0.0;
			}
		}
		for(i=0; i <= im[ic]-1; i++)
		{
			aax[0][0]=aax[0][0]+1.0;
			aax[0][1]=aax[0][1]+sin(angle[ic][i]);
			aax[0][2]=aax[0][2]+cos(angle[ic][i]);
			bbx[0]=bbx[0]+(zm[ic][i]-zcentre[ic]);
			aax[1][0]=aax[1][0]+sin(angle[ic][i]);
			aax[1][1]=aax[1][1]+sin(angle[ic][i])*sin(angle[ic][i]);
			aax[1][2]=aax[1][2]+cos(angle[ic][i])*sin(angle[ic][i]);
			bbx[1]=bbx[1]+(zm[ic][i]-zcentre[ic])*sin(angle[ic][i]);
			aax[2][0]=aax[2][0]+cos(angle[ic][i]);
			aax[2][1]=aax[2][1]+cos(angle[ic][i])*sin(angle[ic][i]);
			aax[2][2]=aax[2][2]+cos(angle[ic][i])*cos(angle[ic][i]);
			bbx[2]=bbx[2]+(zm[ic][i]-zcentre[ic])*cos(angle[ic][i]);
		}
		for(i=0;i<n;i++)
		{
			if(ic==0)fprintf(out,"z--hand %d %le \n",i,bbx[i]);
			for(j=0;j<n;j++)
			{
				if(ic==0)fprintf(out,"z--hand      %d %d %le \n",i,j,aax[i][j]);
			}
		}
		ciggj(aax,n,bbx);
		for( k=0;k<=n-1;k++)bi[ic][k]=bbx[k];
	}
	for(ic=0; ic <= icurve-1; ic++)
	{
 		fprintf(out2,"curve -d 3\n");
		for(i=0; i <= im[ic]-1; i++)
		{
			x1=xcentre[ic]+ai[ic][0]+ai[ic][1]*cos(angle[ic][i])+ai[ic][2]*sin(angle[ic][i]);
			y1=ycentre[ic];
			z1=zcentre[ic]+bi[ic][0]+bi[ic][1]*sin(angle[ic][i])+bi[ic][2]*cos(angle[ic][i]);
			fprintf(out2,"-p    %le %le %le\n",x1,y1,z1);
		}
		fprintf(out2,";\n");
	}

	for(ic=0; ic <= icurve-1; ic++)
	{
		fprintf(out,"hand %d %le %le\n", ic,ai[ic][0],bi[ic][0]);
		for(int j=1;j<n;j++)
		{ 
			fprintf(out,"hand      %d %d %le %le\n", ic,j,ai[ic][j],bi[ic][j]);
		}
	}
*/


/*
	int k,ic;
	for(ic=0; ic <= icurve-1; ic++)
	{
		for(j=0;j<n;j++)
		{
			bbx[j]=0.0;
			for(k=0;k<n;k++)aax[j][k]=0.0;
		}
		for(k=1;k<=jb;k++)
		{
			for(i=0; i <= im[ic]-1; i++)
			{
				bbx[k-1]=bbx[k-1]+(xm[ic][i]-xcentre[ic])*cos(k*angle[ic][i]);
				bbx[jb+k-1]=bbx[jb+k-1]+(xm[ic][i]-xcentre[ic])*sin(k*angle[ic][i]);
				aax[k-1][0]=aax[k-1][0]+cos(k*angle[ic][i]);
				aax[jb+k-1][0]=aax[jb+k-1][0]+sin(k*angle[ic][i]);
			}
			for(j=1;j<=jb;j++)
			{
				for(i=0; i <= im[ic]-1; i++)
				{				
					aax[k-1][2*j-2]=aax[k-1][2*j-2]+cos(j*angle[ic][i])*
						cos(k*angle[ic][i]);			
					aax[k-1][2*j-1]=aax[k-1][2*j-1]+sin(j*angle[ic][i])*
						cos(k*angle[ic][i]);
					aax[jb+k-1][2*j-2]=aax[jb+k-1][2*j-2]+cos(j*angle[ic][i])*
						sin(k*angle[ic][i]);			
					aax[jb+k-1][2*j-1]=aax[jb+k-1][2*j-1]+sin(j*angle[ic][i])*
						sin(k*angle[ic][i]);
				}
			}
		}
		for(i=0;i<n;i++)
		{
			if(ic==0)fprintf(out,"%d %le \n",i,bbx[i]);
			for(j=0;j<n;j++)
			{
				if(ic==0)fprintf(out,"%d %d %le \n",i,j,aax[i][j]);
			}
		}
	    ciggj(aax,n,bbx);
		for( k=0;k<=n-1;k++)ai[ic][k+1]=bbx[k];
	}

	for(ic=0; ic <= icurve-1; ic++)
	{
		for(j=0;j<n;j++)
		{
			bbx[j]=0.0;
			for(k=0;k<n;k++)aax[j][k]=0.0;
		}
		for(k=1;k<=jb;k++)
		{
			for(i=0; i <= im[ic]-1; i++)
			{
				bbx[k-1]=bbx[k-1]+(zm[ic][i]-zcentre[ic])*cos(k*angle[ic][i]);
				bbx[jb+k-1]=bbx[jb+k-1]+(zm[ic][i]-zcentre[ic])*sin(k*angle[ic][i]);
				aax[k-1][0]=aax[k-1][0]+cos(k*angle[ic][i]);
				aax[jb+k-1][0]=aax[jb+k-1][0]+sin(k*angle[ic][i]);
			}
			for(j=1;j<=jb;j++)
			{
				for(i=0; i <= im[ic]-1; i++)
				{				
					aax[k-1][2*j-2]=aax[k-1][2*j-2]+cos(j*angle[ic][i])*
						cos(k*angle[ic][i]);			
					aax[k-1][2*j-1]=aax[k-1][2*j-1]+sin(j*angle[ic][i])*
						cos(k*angle[ic][i]);
					aax[jb+k-1][2*j-2]=aax[jb+k-1][2*j-2]+cos(j*angle[ic][i])*
						sin(k*angle[ic][i]);			
					aax[jb+k-1][2*j-1]=aax[jb+k-1][2*j-1]+sin(j*angle[ic][i])*
						sin(k*angle[ic][i]);
				}
			}
		}
		for(i=0;i<n;i++)
		{
			if(ic==0)fprintf(out,"%d %le \n",i,bbx[i]);
			for(j=0;j<n;j++)
			{
				if(ic==0)fprintf(out,"%d %d %le \n",i,j,aax[i][j]);
			}
		}
	    ciggj(aax,n,bbx);
		for( k=0;k<=n-1;k++)bi[ic][k+1]=bbx[k];
	}

	for(ic=0; ic <= icurve-1; ic++)
	{
 		fprintf(out2,"curve -d 3\n");
		x1=xcentre[ic];
		y1=ycentre[ic];
		z1=zcentre[ic];
		for(int j=1;j<=jb;j++)
		{        
			for(i=0; i <= im[i]-1; i++)
			{
                x1=x1+ai[ic][2*j-1]*cos(j*angle[ic][i])
					+ai[ic][2*j]*sin(j*angle[ic][i]);
				y1=y1+0.0;
				z1=z1+bi[ic][2*j-1]*cos(j*angle[ic][i])
					+bi[ic][2*j]*sin(j*angle[ic][i]);
			}
		    if(i==0 && j<5)printf("\n");
		    fprintf(out2,"-p    %le %le %le\n",x1,y1,z1);
		}
		fprintf(out2,";\n");
	}

	for(ic=0; ic <= icurve-1; ic++)
	{
		for(int j=1;j<n;j++)
		{ 
			fprintf(out,"%d %d %le %le\n", ic,j,ai[ic][j],bi[ic][j]);
		}
	}
*/
/*
	for(i=0; i <= icurve-1; i++)
	{
		for(int k=1;k<=jb;k++)
		{
            bbx[k-1]=0.0;
			for(j=0; j <= im[i]-1; j++)bbx[k-1]=bbx[k-1]+(xm[i][j]-xcentre[i])*cos(k*angle[i][j]);
			for(int l=1; l<=jb;l++)
			{
				aax[k-1][l-1]=0;
				for(j=0; j <= im[i]-1; j++)
				{
					aax[k-1][l-1]=aax[k-1][l-1]+cos(l*angle[i][j])*cos(k*angle[i][j]);
				}
			}
		}
	    ciggj(aax,jb,bbx);
		for( k=1;k<=jb;k++)ai[i][k]=bbx[k-1];
	}

	for(i=0; i <= icurve-1; i++)
	{
		for(int k=1;k<=jb;k++)
		{
            bbx[k-1]=0.0;
			for(j=0; j <= im[i]-1; j++)bbx[k-1]=bbx[k-1]+(zm[i][j]-zcentre[i])*sin(k*angle[i][j]);
			for( l=1; l<=jb;l++)
			{
				aax[k-1][l-1]=0;
				for(j=0; j <= im[i]-1; j++)
				{
					aax[k-1][l-1]=aax[k-1][l-1]+sin(l*angle[i][j])*sin(k*angle[i][j]);
				}
			}
		}
	    ciggj(aax,jb,bbx);
		for(k=1;k<=jb;k++)bi[i][k]=bbx[k-1];
	}

	for(i=0; i <= icurve-1; i++)
	{
 		fprintf(out2,"curve -d 3\n");
        for(j=0; j <= im[i]-1; j++)
		{
			x1=xcentre[i];
			y1=ycentre[i];
			z1=zcentre[i];
			for(int k=1;k<=jb;k++)
			{
                x1=x1+ai[i][k]*cos(k*angle[i][j]);
				y1=y1+0.0;
				z1=z1+bi[i][k]*sin(k*angle[i][j]);
//		    if(i==icurve-1 && j>im[i]-5)printf("%d %d %d %le %le %le %le\n",i,j,k,z1,zcentre[i],bi[i][k],angle[i][j]);
			}
		    if(i==0 && j<5)printf("\n");
		    fprintf(out2,"-p    %le %le %le\n",x1,y1,z1);
		}
		fprintf(out2,";\n");
	}

	for(i=0; i <= icurve-1; i++)
	{
        for(j=1; j <= 2; j++)
		{
			fprintf(out,"%d %d %le %le\n",i,j,ai[i][j],bi[i][j]);
		}
	}

	for(i=0; i <= icurve-1; i++)
	{
            bbx[0]=0.0;
            bbx[1]=0.0;
			aax[0][0]=0.0;
			aax[0][1]=0.0;
			aax[1][0]=0.0;
			aax[1][1]=0.0;
			for(j=0; j <= im[i]-1; j++)
			{
				bbx[0]=bbx[0]+(xm[i][j]+xcentre[i])*cos(angle[i][j]);
				bbx[1]=bbx[1]+(xm[i][j]+xcentre[i])*cos(2.0*angle[i][j]);
				aax[0][0]=aax[0][0]+cos(angle[i][j])*cos(angle[i][j]);
				aax[0][1]=aax[0][1]+cos(2.0*angle[i][j])*cos(angle[i][j]);
				aax[1][0]=aax[1][0]+cos(angle[i][j])*cos(2.0*angle[i][j]);
				aax[1][1]=aax[1][1]+cos(2.0*angle[i][j])*cos(2.0*angle[i][j]);
			}
			ai[i][1]=(bbx[0]*aax[1][1]-bbx[1]*aax[0][1])/(aax[0][0]*aax[1][1]-aax[0][1]*aax[1][0]);
			ai[i][2]=(bbx[1]*aax[0][0]-bbx[0]*aax[1][0])/(aax[0][0]*aax[1][1]-aax[0][1]*aax[1][0]);
	}
*/
/*
	for(i=0; i <= icurve-1; i++)
	{
            bbx[0]=0.0;
            bbx[1]=0.0;
			aax[0][0]=0.0;
			aax[0][1]=0.0;
			aax[1][0]=0.0;
			aax[1][1]=0.0;
			for(j=0; j <= im[i]-1; j++)
			{
				bbx[0]=bbx[0]+(zm[i][j]+zcentre[i])*sin(angle[i][j]);
				bbx[1]=bbx[1]+(zm[i][j]+zcentre[i])*sin(2.0*angle[i][j]);
				aax[0][0]=aax[0][0]+sin(angle[i][j])*sin(angle[i][j]);
				aax[0][1]=aax[0][1]+sin(2.0*angle[i][j])*sin(angle[i][j]);
				aax[1][0]=aax[1][0]+sin(angle[i][j])*sin(2.0*angle[i][j]);
				aax[1][1]=aax[1][1]+sin(2.0*angle[i][j])*sin(2.0*angle[i][j]);
			}
			bi[i][1]=(bbx[0]*aax[1][1]-bbx[1]*aax[0][1])/(aax[0][0]*aax[1][1]-aax[0][1]*aax[1][0]);
			bi[i][2]=(bbx[1]*aax[0][0]-bbx[0]*aax[1][0])/(aax[0][0]*aax[1][1]-aax[0][1]*aax[1][0]);
	}

	for(i=0; i <= icurve-1; i++)
	{
        for(j=1; j <= 2; j++)
		{
			fprintf(out,"jb   %d %d %le %le\n",i,j,ai[i][j],bi[i][j]);
		}
	}
*/

/*
	im[0]=200;
	i=0;
		fprintf(out2,"curve -d 3\n");
        for(j=0; j <= im[0]-1; j++)
		{
            fscanf(in1,"%le %le %le\n",&xm[i][j],&ym[i][j],&zm[i][j]);
            fprintf(out,"%le %le %le\n",xm[i][j],ym[i][j],zm[i][j]);
		    fprintf(out2,"-p   %le %le %le\n",xm[i][j],ym[i][j],zm[i][j]);
		}
		fprintf(out2,";\n");
		i=0;
		j=0;
		fprintf(out2,"curve -d 3\n");
        do
		{
            fscanf(in1,"%le %le %le\n",&xm[i][j],&ym[i][j],&zm[i][j]);
            fprintf(out,"%le %le %le\n",xm[i][j],ym[i][j],zm[i][j]);
		    fprintf(out2,"-p   %le %le %le\n",xm[i][j],ym[i][j],zm[i][j]);
			j++;
		}while (xm[i][j-1]!=0.0 && ym[i][j-1] != 0.0 && zm[i][j-1] != 0.0);
		fprintf(out2,";\n");

	for(i=0; i <= icurve-1; i++)
	{
//		fscanf(in1,"%d\n",&im[i]);
//		fprintf(out,"%d\n",im[i]);
		fprintf(out2,"curve -d 3\n");
        for(j=0; j <= im[i]-1; j++)
		{
            fscanf(in1,"%le %le %le\n",&xm[i][j],&ym[i][j],&zm[i][j]);
            fprintf(out,"%le %le %le\n",xm[i][j],ym[i][j],zm[i][j]);
		    fprintf(out2,"-p   %le %le %le\n",xm[i][j],ym[i][j],zm[i][j]);
		}
		fprintf(out2,";\n");
	}
    printf("pass 0\n");
*/
/*
	for(i=0; i <= icurve-1; i++)
	{
        fscanf(in2,"%le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le\n",&anx[i][0],
		   &anx[i][1],&any[i][0],&any[i][1],&anz[i][0],&txi[i],&tyi[i],&tzi[i],&nxi[i][0],&nyi[i][0],
		   &nzi[i][0],&nxi[i][1],&nyi[i][1],&nzi[i][1],&nxi[i][2],&nyi[i][2],&nzi[i][2]);
        fprintf(out,"%le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le\n",anx[i][0],
		   anx[i][1],any[i][0],any[i][1],anz[i][0],txi[i],tyi[i],tzi[i],nxi[i][0],nyi[i][0],
		   nzi[i][0],nxi[i][1],nyi[i][1],nzi[i][1],nxi[i][2],nyi[i][2],nzi[i][2]);
	}
    printf("pass 1\n");
*/

   for(j=0; j <= icurve-1; j++)
   for (i = 0; i <= im[j]-2; i++)
   {
	 glColor3f(0.0, 0.0, 1.0); 
     glBegin(GL_LINES);
      glVertex3f (xm[j][i]/10.0,ym[j][i]/10.0,zm[j][i]/10.0);
      glVertex3f (xm[j][i+1]/10.0,ym[j][i+1]/10.0,zm[j][i+1]/10.0);
     glEnd(); 
   }
   for(j=0; j <= icurve-1; j++)
   {
	 glColor3f(1.0, 0.0, 0.0); 
     glBegin(GL_POINTS);
      glVertex3f (xm[j][i]/10.0,ym[j][i]/10.0,zm[j][i]/10.0);
     glEnd(); 
   }


/*
   dv=2.0*pi/100.0;
   for(j=0; j <= icurve-1; j++)
   {
	   nx1=nxi[j][0];
	   ny1=nyi[j][0];
	   nz1=nzi[j][0];
	   nx2=nxi[j][1];
	   ny2=nyi[j][1];
	   nz2=nzi[j][1];
	   nx3=nxi[j][2];
	   ny3=nyi[j][2];
	   nz3=nzi[j][2];
	   tx1=txi[j];
	   ty1=tyi[j];
	   tz1=tzi[j];
	   fprintf(out2,"curve -d 3\n");
   for (i = 0; i <= 99; i++)
   {
     vv0 = i * dv;
     vv1 = (i + 1) * dv;
 	 z1=anz[j][0]-tz1;
	 z2=z1;
	 ff=nx1*ny2*nz3+ny1*nz2*nx3+nz1*nx2*ny3-nz1*ny2*nx3-nx1*nz2*ny3-ny1*nx2*nz3;
	 x3=(ny2*nz3*(anx[j][0]-tx1+anx[j][1]*sin(vv0))+ny1*nz2*z1+nz1*ny3*(any[j][0]-ty1+any[j][1]
	     *cos(vv0))-nz1*ny2*z1-nz2*ny3*(anx[j][0]-tx1+anx[j][1]*sin(vv0))-ny1*nz3*(any[j][0]-ty1
		 +any[j][1]*cos(vv0)))/ff;
	 x4=(ny2*nz3*(anx[j][0]-tx1+anx[j][1]*sin(vv1))+ny1*nz2*z2+nz1*ny3*(any[j][0]-ty1+any[j][1]
	     *cos(vv1))-nz1*ny2*z2-nz2*ny3*(anx[j][0]-tx1+anx[j][1]*sin(vv1))-ny1*nz3*(any[j][0]-ty1
		 +any[j][1]*cos(vv1)))/ff;
	 y3=(nx1*nz3*(any[j][0]-ty1+any[j][1]*cos(vv0))+nz2*nx3*(anx[j][0]-tx1+anx[j][1]*sin(vv0))
		 +nz1*nx2*z1-nz1*nx3*(any[j][0]-ty1+any[j][1]*cos(vv0))-nx1*nz2*z1-nx2*nz3*(anx[j][0]-tx1
		 +anx[j][1]*sin(vv0)))/ff;
	 y4=(nx1*nz3*(any[j][0]-ty1+any[j][1]*cos(vv1))+nz2*nx3*(anx[j][0]-tx1+anx[j][1]*sin(vv1))
		 +nz1*nx2*z2-nz1*nx3*(any[j][0]-ty1+any[j][1]*cos(vv1))-nx1*nz2*z2-nx2*nz3*(anx[j][0]-tx1
		 +anx[j][1]*sin(vv1)))/ff;
	 z3=(nx1*ny2*z1+ny1*nx3*(any[j][0]-ty1+any[j][1]*cos(vv0))+nx2*ny3*(anx[j][0]-tx1+anx[j][1]
	     *sin(vv0))-ny2*nx3*(anx[j][0]-tx1+anx[j][1]*sin(vv0))-nx1*ny3*(any[j][0]-ty1+any[j][1]
		 *cos(vv0))-ny1*nx2*z1)/ff;
	 z4=(nx1*ny2*z2+ny1*nx3*(any[j][0]-ty1+any[j][1]*cos(vv1))+nx2*ny3*(anx[j][0]-tx1+anx[j][1]
	     *sin(vv1))-ny2*nx3*(anx[j][0]-tx1+anx[j][1]*sin(vv1))-nx1*ny3*(any[j][0]-ty1+any[j][1]
		 *cos(vv1))-ny1*nx2*z2)/ff;
 
 	 glColor3f(1.0,1.0,0.0); 
     glBegin(GL_LINES);
      glVertex3f (x3,y3,z3);
      glVertex3f (x4,y4,z4);
     glEnd(); 
	 fprintf(out2,"-p %le %le %le\n",x3,y3,z3);
 
   }
     fprintf(out2,";\n");  
   }

   jj=0;
L20: printf("");
//  constant term of x
   j=jj;
   w1=0.0;
   w2=0.0;
	   nx1=nxi[j][0];
	   ny1=nyi[j][0];
	   nz1=nzi[j][0];
	   nx2=nxi[j][1];
	   ny2=nyi[j][1];
	   nz2=nzi[j][1];
	   nx3=nxi[j][2];
	   ny3=nyi[j][2];
	   nz3=nzi[j][2];
	   tx1=txi[j];
	   ty1=tyi[j];
	   tz1=tzi[j];
 	 z1=anz[j][0]-tz1;
	 z2=z1;
	 ff=nx1*ny2*nz3+ny1*nz2*nx3+nz1*nx2*ny3-nz1*ny2*nx3-nx1*nz2*ny3-ny1*nx2*nz3;
	 h0=(ny2*nz3*(anx[j][0]-tx1)+ny1*nz2*z1+nz1*ny3*(any[j][0]-ty1)-nz1*ny2*z1-nz2*ny3
		 *(anx[j][0]-tx1)-ny1*nz3*(any[j][0]-ty1))/ff;
   h01=w1*h0;
   h02=w2*h0;

   j=jj+1;
   w1=0.0;
   w2=0.0;
	   nx1=nxi[j][0];
	   ny1=nyi[j][0];
	   nz1=nzi[j][0];
	   nx2=nxi[j][1];
	   ny2=nyi[j][1];
	   nz2=nzi[j][1];
	   nx3=nxi[j][2];
	   ny3=nyi[j][2];
	   nz3=nzi[j][2];
	   tx1=txi[j];
	   ty1=tyi[j];
	   tz1=tzi[j];
 	 z1=anz[j][0]-tz1;
	 z2=z1;
	 ff=nx1*ny2*nz3+ny1*nz2*nx3+nz1*nx2*ny3-nz1*ny2*nx3-nx1*nz2*ny3-ny1*nx2*nz3;
	 h1=(ny2*nz3*(anx[j][0]-tx1)+ny1*nz2*z1+nz1*ny3*(any[j][0]-ty1)-nz1*ny2*z1-nz2*ny3
		 *(anx[j][0]-tx1)-ny1*nz3*(any[j][0]-ty1))/ff;
   h11=w1*h1;
   h12=w2*h1;
   cons(h0,h01,h02,h1,h11,h12,cx0);
   printf("cx0= %le %le %le %le %le %le\n",cx0[0],cx0[1],cx0[2],cx0[3],cx0[4],cx0[5]);

//  constant term of y
   j=jj;
   w1=0.0;
   w2=0.0;
	   nx1=nxi[j][0];
	   ny1=nyi[j][0];
	   nz1=nzi[j][0];
	   nx2=nxi[j][1];
	   ny2=nyi[j][1];
	   nz2=nzi[j][1];
	   nx3=nxi[j][2];
	   ny3=nyi[j][2];
	   nz3=nzi[j][2];
	   tx1=txi[j];
	   ty1=tyi[j];
	   tz1=tzi[j];
 	 z1=anz[j][0]-tz1;
	 z2=z1;
	 ff=nx1*ny2*nz3+ny1*nz2*nx3+nz1*nx2*ny3-nz1*ny2*nx3-nx1*nz2*ny3-ny1*nx2*nz3;
	 h0=(nx1*nz3*(any[j][0]-ty1)+nz2*nx3*(anx[j][0]-tx1)+nz1*nx2*z1-nz1*nx3*(any[j][0]-ty1)
		 -nx1*nz2*z1-nx2*nz3*(anx[j][0]-tx1))/ff;
   h01=w1*h0;
   h02=w2*h0;

   j=jj+1;
   w1=1.0;
   w2=1.0;
	   nx1=nxi[j][0];
	   ny1=nyi[j][0];
	   nz1=nzi[j][0];
	   nx2=nxi[j][1];
	   ny2=nyi[j][1];
	   nz2=nzi[j][1];
	   nx3=nxi[j][2];
	   ny3=nyi[j][2];
	   nz3=nzi[j][2];
	   tx1=txi[j];
	   ty1=tyi[j];
	   tz1=tzi[j];
 	 z1=anz[j][0]-tz1;
	 z2=z1;
	 ff=nx1*ny2*nz3+ny1*nz2*nx3+nz1*nx2*ny3-nz1*ny2*nx3-nx1*nz2*ny3-ny1*nx2*nz3;
	 h1=(nx1*nz3*(any[j][0]-ty1)+nz2*nx3*(anx[j][0]-tx1)+nz1*nx2*z1-nz1*nx3*(any[j][0]-ty1)
		 -nx1*nz2*z1-nx2*nz3*(anx[j][0]-tx1))/ff;
   h11=w1*h1;
   h12=w2*h1;
   cons(h0,h01,h02,h1,h11,h12,cy0);
   printf("cy0= %le %le %le %le %le %le\n",cy0[0],cy0[1],cy0[2],cy0[3],cy0[4],cy0[5]);

//  constant term of z
   j=jj;
   w1=0.0;
   w2=0.0;
	   nx1=nxi[j][0];
	   ny1=nyi[j][0];
	   nz1=nzi[j][0];
	   nx2=nxi[j][1];
	   ny2=nyi[j][1];
	   nz2=nzi[j][1];
	   nx3=nxi[j][2];
	   ny3=nyi[j][2];
	   nz3=nzi[j][2];
	   tx1=txi[j];
	   ty1=tyi[j];
	   tz1=tzi[j];
 	 z1=anz[j][0]-tz1;
	 z2=z1;
	 ff=nx1*ny2*nz3+ny1*nz2*nx3+nz1*nx2*ny3-nz1*ny2*nx3-nx1*nz2*ny3-ny1*nx2*nz3;
	 h0=(nx1*ny2*z1+ny1*nx3*(any[j][0]-ty1)+nx2*ny3*(anx[j][0]-tx1)-ny2*nx3*(anx[j][0]-tx1)
		 -nx1*ny3*(any[j][0]-ty1)-ny1*nx2*z1)/ff;
   h01=w1*h0;
   h02=w2*h0;

   j=jj+1;
   w1=0.0;
   w2=0.0;
	   nx1=nxi[j][0];
	   ny1=nyi[j][0];
	   nz1=nzi[j][0];
	   nx2=nxi[j][1];
	   ny2=nyi[j][1];
	   nz2=nzi[j][1];
	   nx3=nxi[j][2];
	   ny3=nyi[j][2];
	   nz3=nzi[j][2];
	   tx1=txi[j];
	   ty1=tyi[j];
	   tz1=tzi[j];
 	 z1=anz[j][0]-tz1;
	 z2=z1;
	 ff=nx1*ny2*nz3+ny1*nz2*nx3+nz1*nx2*ny3-nz1*ny2*nx3-nx1*nz2*ny3-ny1*nx2*nz3;
	 h1=(nx1*ny2*z1+ny1*nx3*(any[j][0]-ty1)+nx2*ny3*(anx[j][0]-tx1)-ny2*nx3*(anx[j][0]-tx1)
		 -nx1*ny3*(any[j][0]-ty1)-ny1*nx2*z1)/ff;
   h11=w1*h1;
   h12=w2*h1;
   cons(h0,h01,h02,h1,h11,h12,cz0);
   printf("cz0= %le %le %le %le %le %le\n",cz0[0],cz0[1],cz0[2],cz0[3],cz0[4],cz0[5]);

//   goto L10;

   ax = 1.0;
   bx = 1.0;
   cx = 1.0;
   dx = 1.0;
   ay = 1.0;
   by = 1.0;
   cy = 1.0;
   dy = 1.0;
   az = 1.0;
   bz = 1.0;
   cz = 1.0;
   dz = 1.0;
   di2 = - 1.0;
   di4 = 1.0;
   di6 = - 1.0;
//  constant sinv term of x
   j=jj;
   w1=0.0;
   w2=0.0;
	   nx1=nxi[j][0];
	   ny1=nyi[j][0];
	   nz1=nzi[j][0];
	   nx2=nxi[j][1];
	   ny2=nyi[j][1];
	   nz2=nzi[j][1];
	   nx3=nxi[j][2];
	   ny3=nyi[j][2];
	   nz3=nzi[j][2];
	   tx1=txi[j];
	   ty1=tyi[j];
	   tz1=tzi[j];
 	 z1=anz[j][0]-tz1;
	 z2=z1;
	 ff=nx1*ny2*nz3+ny1*nz2*nx3+nz1*nx2*ny3-nz1*ny2*nx3-nx1*nz2*ny3-ny1*nx2*nz3;
	 ai0=(ny2*nz3*(anx[j][1])-nz2*ny3*(anx[j][1]))/ff;
   ai1=w1*ai0;
   ai2=w2*ai0;

   j=jj+1;
   w1=0.0;
   w2=0.0;
	   nx1=nxi[j][0];
	   ny1=nyi[j][0];
	   nz1=nzi[j][0];
	   nx2=nxi[j][1];
	   ny2=nyi[j][1];
	   nz2=nzi[j][1];
	   nx3=nxi[j][2];
	   ny3=nyi[j][2];
	   nz3=nzi[j][2];
	   tx1=txi[j];
	   ty1=tyi[j];
	   tz1=tzi[j];
 	 z1=anz[j][0]-tz1;
	 z2=z1;
	 ff=nx1*ny2*nz3+ny1*nz2*nx3+nz1*nx2*ny3-nz1*ny2*nx3-nx1*nz2*ny3-ny1*nx2*nz3;
	 ai3=(ny2*nz3*(anx[j][1])-nz2*ny3*(anx[j][1]))/ff;
   ai4=w1*ai3;
   ai5=w2*ai3;
   sincos(ax,bx,cx,dx,ai0,ai1,ai2,ai3,ai4,ai5,di2,di4,di6,csx);
   printf("csx[0-11]= %le %le %le %le %le %le %le %le %le %le %le %le\n",csx[0],csx[1],csx[2],
	   csx[3],csx[4],csx[5],csx[6],csx[7],csx[8],csx[9],csx[10],csx[11]);

//  constant sinv term of y
   j=jj;
   w1=0.0;
   w2=0.0;
	   nx1=nxi[j][0];
	   ny1=nyi[j][0];
	   nz1=nzi[j][0];
	   nx2=nxi[j][1];
	   ny2=nyi[j][1];
	   nz2=nzi[j][1];
	   nx3=nxi[j][2];
	   ny3=nyi[j][2];
	   nz3=nzi[j][2];
	   tx1=txi[j];
	   ty1=tyi[j];
	   tz1=tzi[j];
 	 z1=anz[j][0]-tz1;
	 z2=z1;
	 ff=nx1*ny2*nz3+ny1*nz2*nx3+nz1*nx2*ny3-nz1*ny2*nx3-nx1*nz2*ny3-ny1*nx2*nz3;
	 ai0=(nz2*nx3*(anx[j][1])-nx2*nz3*(+anx[j][1]))/ff;
   ai1=w1*ai0;
   ai2=w2*ai0;

	 y3=(nx1*nz3*(any[j][0]-ty1+any[j][1]*cos(vv0))+nz2*nx3*(anx[j][0]-tx1+anx[j][1]*sin(vv0))
		 +nz1*nx2*z1-nz1*nx3*(any[j][0]-ty1+any[j][1]*cos(vv0))-nx1*nz2*z1-nx2*nz3*(anx[j][0]-tx1
		 +anx[j][1]*sin(vv0)))/ff;

   j=jj+1;
   w1=0.0;
   w2=0.0;
	   nx1=nxi[j][0];
	   ny1=nyi[j][0];
	   nz1=nzi[j][0];
	   nx2=nxi[j][1];
	   ny2=nyi[j][1];
	   nz2=nzi[j][1];
	   nx3=nxi[j][2];
	   ny3=nyi[j][2];
	   nz3=nzi[j][2];
	   tx1=txi[j];
	   ty1=tyi[j];
	   tz1=tzi[j];
 	 z1=anz[j][0]-tz1;
	 z2=z1;
	 ff=nx1*ny2*nz3+ny1*nz2*nx3+nz1*nx2*ny3-nz1*ny2*nx3-nx1*nz2*ny3-ny1*nx2*nz3;
	 ai3=(nz2*nx3*(anx[j][1])-nx2*nz3*(+anx[j][1]))/ff;
	 ai4=w1*ai3;
   ai5=w2*ai3;
   sincos(ax,bx,cx,dx,ai0,ai1,ai2,ai3,ai4,ai5,di2,di4,di6,csy);
   printf("csy[0-11]= %le %le %le %le %le %le %le %le %le %le %le %le\n",csy[0],csy[1],csy[2],
	   csy[3],csy[4],csy[5],csy[6],csy[7],csy[8],csy[9],csy[10],csy[11]);

//  constant sinv term of z
   j=jj;
   w1=0.0;
   w2=0.0;
	   nx1=nxi[j][0];
	   ny1=nyi[j][0];
	   nz1=nzi[j][0];
	   nx2=nxi[j][1];
	   ny2=nyi[j][1];
	   nz2=nzi[j][1];
	   nx3=nxi[j][2];
	   ny3=nyi[j][2];
	   nz3=nzi[j][2];
	   tx1=txi[j];
	   ty1=tyi[j];
	   tz1=tzi[j];
 	 z1=anz[j][0]-tz1;
	 z2=z1;
	 ff=nx1*ny2*nz3+ny1*nz2*nx3+nz1*nx2*ny3-nz1*ny2*nx3-nx1*nz2*ny3-ny1*nx2*nz3;
	 ai0=(nx2*ny3*(anx[j][1])-ny2*nx3*(anx[j][1]))/ff;
   ai1=w1*ai0;
   ai2=w2*ai0;

   j=jj+1;
   w1=0.0;
   w2=0.0;
	   nx1=nxi[j][0];
	   ny1=nyi[j][0];
	   nz1=nzi[j][0];
	   nx2=nxi[j][1];
	   ny2=nyi[j][1];
	   nz2=nzi[j][1];
	   nx3=nxi[j][2];
	   ny3=nyi[j][2];
	   nz3=nzi[j][2];
	   tx1=txi[j];
	   ty1=tyi[j];
	   tz1=tzi[j];
 	 z1=anz[j][0]-tz1;
	 z2=z1;
	 ff=nx1*ny2*nz3+ny1*nz2*nx3+nz1*nx2*ny3-nz1*ny2*nx3-nx1*nz2*ny3-ny1*nx2*nz3;
	 ai3=(nx2*ny3*(anx[j][1])-ny2*nx3*(anx[j][1]))/ff;
	 ai4=w1*ai3;
   ai5=w2*ai3;
   sincos(ax,bx,cx,dx,ai0,ai1,ai2,ai3,ai4,ai5,di2,di4,di6,csz);
   printf("csz[0-11]= %le %le %le %le %le %le %le %le %le %le %le %le\n",csz[0],csz[1],csz[2],
	   csz[3],csz[4],csz[5],csz[6],csz[7],csz[8],csz[9],csz[10],csz[11]);

//   goto L10;

//  constant cosv term of x
   j=jj;
   w1=0.0;
   w2=0.0;
	   nx1=nxi[j][0];
	   ny1=nyi[j][0];
	   nz1=nzi[j][0];
	   nx2=nxi[j][1];
	   ny2=nyi[j][1];
	   nz2=nzi[j][1];
	   nx3=nxi[j][2];
	   ny3=nyi[j][2];
	   nz3=nzi[j][2];
	   tx1=txi[j];
	   ty1=tyi[j];
	   tz1=tzi[j];
 	 z1=anz[j][0]-tz1;
	 z2=z1;
	 ff=nx1*ny2*nz3+ny1*nz2*nx3+nz1*nx2*ny3-nz1*ny2*nx3-nx1*nz2*ny3-ny1*nx2*nz3;
	 ai0=(nz1*ny3*(any[j][1])-ny1*nz3*(any[j][1]))/ff;
   ai1=w1*ai0;
   ai2=w2*ai0;

   j=jj+1;
   w1=0.0;
   w2=0.0;
	   nx1=nxi[j][0];
	   ny1=nyi[j][0];
	   nz1=nzi[j][0];
	   nx2=nxi[j][1];
	   ny2=nyi[j][1];
	   nz2=nzi[j][1];
	   nx3=nxi[j][2];
	   ny3=nyi[j][2];
	   nz3=nzi[j][2];
	   tx1=txi[j];
	   ty1=tyi[j];
	   tz1=tzi[j];
 	 z1=anz[j][0]-tz1;
	 z2=z1;
	 ff=nx1*ny2*nz3+ny1*nz2*nx3+nz1*nx2*ny3-nz1*ny2*nx3-nx1*nz2*ny3-ny1*nx2*nz3;
	 ai3=(nz1*ny3*(any[j][1])-ny1*nz3*(any[j][1]))/ff;
	 ai4=w1*ai3;
   ai5=w2*ai3;
   sincos(ax,bx,cx,dx,ai0,ai1,ai2,ai3,ai4,ai5,di2,di4,di6,ccx);
   printf("ccx[0-11]= %le %le %le %le %le %le %le %le %le %le %le %le\n",ccx[0],ccx[1],ccx[2],
	   ccx[3],ccx[4],ccx[5],ccx[6],ccx[7],ccx[8],ccx[9],ccx[10],ccx[11]);

//  constant cosv term of y
   j=jj;
   w1=0.0;
   w2=0.0;
	   nx1=nxi[j][0];
	   ny1=nyi[j][0];
	   nz1=nzi[j][0];
	   nx2=nxi[j][1];
	   ny2=nyi[j][1];
	   nz2=nzi[j][1];
	   nx3=nxi[j][2];
	   ny3=nyi[j][2];
	   nz3=nzi[j][2];
	   tx1=txi[j];
	   ty1=tyi[j];
	   tz1=tzi[j];
 	 z1=anz[j][0]-tz1;
	 z2=z1;
	 ff=nx1*ny2*nz3+ny1*nz2*nx3+nz1*nx2*ny3-nz1*ny2*nx3-nx1*nz2*ny3-ny1*nx2*nz3;
	 ai0=(nx1*nz3*(any[j][1])-nz1*nx3*(any[j][1]))/ff;
   ai1=w1*ai0;
   ai2=w2*ai0;

   j=jj+1;
   w1=0.0;
   w2=0.0;
	   nx1=nxi[j][0];
	   ny1=nyi[j][0];
	   nz1=nzi[j][0];
	   nx2=nxi[j][1];
	   ny2=nyi[j][1];
	   nz2=nzi[j][1];
	   nx3=nxi[j][2];
	   ny3=nyi[j][2];
	   nz3=nzi[j][2];
	   tx1=txi[j];
	   ty1=tyi[j];
	   tz1=tzi[j];
 	 z1=anz[j][0]-tz1;
	 z2=z1;
	 ff=nx1*ny2*nz3+ny1*nz2*nx3+nz1*nx2*ny3-nz1*ny2*nx3-nx1*nz2*ny3-ny1*nx2*nz3;
	 ai3=(nx1*nz3*(any[j][1])-nz1*nx3*(any[j][1]))/ff;
   ai4=w1*ai3;
   ai5=w2*ai3;
   sincos(ax,bx,cx,dx,ai0,ai1,ai2,ai3,ai4,ai5,di2,di4,di6,ccy);
   printf("ccy[0-11]= %le %le %le %le %le %le %le %le %le %le %le %le\n",ccy[0],ccy[1],ccy[2],
	   ccy[3],ccy[4],ccy[5],ccy[6],ccy[7],ccy[8],ccy[9],ccy[10],ccy[11]);

//  constant cosv term of z
   j=jj;
   w1=0.0;
   w2=0.0;
	   nx1=nxi[j][0];
	   ny1=nyi[j][0];
	   nz1=nzi[j][0];
	   nx2=nxi[j][1];
	   ny2=nyi[j][1];
	   nz2=nzi[j][1];
	   nx3=nxi[j][2];
	   ny3=nyi[j][2];
	   nz3=nzi[j][2];
	   tx1=txi[j];
	   ty1=tyi[j];
	   tz1=tzi[j];
 	 z1=anz[j][0]-tz1;
	 z2=z1;
	 ff=nx1*ny2*nz3+ny1*nz2*nx3+nz1*nx2*ny3-nz1*ny2*nx3-nx1*nz2*ny3-ny1*nx2*nz3;
	 ai0=(ny1*nx3*(any[j][1])-nx1*ny3*(any[j][1]))/ff;
   ai1=w1*ai0;
   ai2=w2*ai0;

   j=jj+1;
   w1=0.0;
   w2=0.0;
	   nx1=nxi[j][0];
	   ny1=nyi[j][0];
	   nz1=nzi[j][0];
	   nx2=nxi[j][1];
	   ny2=nyi[j][1];
	   nz2=nzi[j][1];
	   nx3=nxi[j][2];
	   ny3=nyi[j][2];
	   nz3=nzi[j][2];
	   tx1=txi[j];
	   ty1=tyi[j];
	   tz1=tzi[j];
 	 z1=anz[j][0]-tz1;
	 z2=z1;
	 ff=nx1*ny2*nz3+ny1*nz2*nx3+nz1*nx2*ny3-nz1*ny2*nx3-nx1*nz2*ny3-ny1*nx2*nz3;
	 ai3=(ny1*nx3*(any[j][1])-nx1*ny3*(any[j][1]))/ff;
   ai4=w1*ai3;
   ai5=w2*ai3;
   sincos(ax,bx,cx,dx,ai0,ai1,ai2,ai3,ai4,ai5,di2,di4,di6,ccz);
   printf("ccy[0-11]= %le %le %le %le %le %le %le %le %le %le %le %le\n",ccz[0],ccz[1],ccz[2],
	   ccz[3],ccz[4],ccz[5],ccz[6],ccz[7],ccz[8],ccz[9],ccz[10],ccz[11]);

//  generate PDE surface 
   du = 1.0 / 1.0;
   dv = 2.0 * pi / 100.0;
   for (iv = 0; iv <= 100; iv++)
   {
     vv0 = iv * dv;
     vv1 = (iv + 1) * dv;
     for (iu = 0; iu <= 1; iu ++)
     {
       u01 = iu * du;
       u11 = (iu + 1) * du; 
	   u02=u01*u01;
	   u03=u02*u01;
	   u04=u02*u02;
	   u05=u03*u02;
	   u12=u11*u11;
	   u13=u12*u11;
	   u14=u12*u12;
	   u15=u13*u12;
       x1=cx0[0]+cx0[1]*u01+cx0[2]*u02+cx0[3]*u03+cx0[4]*u04+cx0[5]*u05+csx[0]*sin(vv0)+ccx[0]*cos(vv0);
       x2=cx0[0]+cx0[1]*u01+cx0[2]*u02+cx0[3]*u03+cx0[4]*u04+cx0[5]*u05+csx[0]*sin(vv1)+ccx[0]*cos(vv1);
       x3=cx0[0]+cx0[1]*u11+cx0[2]*u12+cx0[3]*u13+cx0[4]*u14+cx0[5]*u15+csx[0]*sin(vv1)+ccx[0]*cos(vv1);
       x4=cx0[0]+cx0[1]*u11+cx0[2]*u12+cx0[3]*u13+cx0[4]*u14+cx0[5]*u15+csx[0]*sin(vv0)+ccx[0]*cos(vv0);
       y1=cy0[0]+cy0[1]*u01+cy0[2]*u02+cy0[3]*u03+cy0[4]*u04+cy0[5]*u05+csy[0]*sin(vv0)+ccy[0]*cos(vv0);
       y2=cy0[0]+cy0[1]*u01+cy0[2]*u02+cy0[3]*u03+cy0[4]*u04+cy0[5]*u05+csy[0]*sin(vv1)+ccy[0]*cos(vv1); 
       y3=cy0[0]+cy0[1]*u11+cy0[2]*u12+cy0[3]*u13+cy0[4]*u14+cy0[5]*u15+csy[0]*sin(vv1)+ccy[0]*cos(vv1); 
       y4=cy0[0]+cy0[1]*u11+cy0[2]*u12+cy0[3]*u13+cy0[4]*u14+cy0[5]*u15+csy[0]*sin(vv0)+ccy[0]*cos(vv0);     
       z1=cz0[0]+cz0[1]*u01+cz0[2]*u02+cz0[3]*u03+cz0[4]*u04+cz0[5]*u05+csz[0]*sin(vv0)+ccz[0]*cos(vv0);
       z2=cz0[0]+cz0[1]*u01+cz0[2]*u02+cz0[3]*u03+cz0[4]*u04+cz0[5]*u05+csz[0]*sin(vv1)+ccz[0]*cos(vv1);
       z3=cz0[0]+cz0[1]*u11+cz0[2]*u12+cz0[3]*u13+cz0[4]*u14+cz0[5]*u15+csz[0]*sin(vv1)+ccz[0]*cos(vv1);
       z4=cz0[0]+cz0[1]*u11+cz0[2]*u12+cz0[3]*u13+cz0[4]*u14+cz0[5]*u15+csz[0]*sin(vv0)+ccz[0]*cos(vv0);
       for (i=1; i <= 11; i++)
       {
         prod1(&i,&u01,&un);
         x1=x1+csx[i]*un*sin(vv0)+ccx[i]*un*cos(vv0);
         x2=x2+csx[i]*un*sin(vv1)+ccx[i]*un*cos(vv1);
         y1=y1+csy[i]*un*sin(vv0)+ccy[i]*un*cos(vv0);
         y2=y2+csy[i]*un*sin(vv1)+ccy[i]*un*cos(vv1);
         z1=z1+csz[i]*un*sin(vv0)+ccz[i]*un*cos(vv0);
         z2=z2+csz[i]*un*sin(vv1)+ccz[i]*un*cos(vv1);
         prod1(&i,&u11,&un);
         x3=x3+csx[i]*un*sin(vv1)+ccx[i]*un*cos(vv1);
         x4=x4+csx[i]*un*sin(vv0)+ccx[i]*un*cos(vv0);
         y3=y3+csy[i]*un*sin(vv1)+ccy[i]*un*cos(vv1);
         y4=y4+csy[i]*un*sin(vv0)+ccy[i]*un*cos(vv0);
         z3=z3+csz[i]*un*sin(vv1)+ccz[i]*un*cos(vv1);
         z4=z4+csz[i]*un*sin(vv0)+ccz[i]*un*cos(vv0);
       }
       nx1 = (y2 - y1) * (z4 - z1) - (z2 - z1) * (y4 - y1);
       ny1 = (z2 - z1) * (x4 - x1) - (x2 - x1) * (z4 - z1);
       nz1 = (x2 - x1) * (y4 - y1) - (y2 - y1) * (x4 - x1);
       nxyz = sqrt(nx1 * nx1 + ny1 * ny1 + nz1 * nz1);
       nx1 = nx1 / nxyz;
       ny1 = ny1 / nxyz;
       nz1 = nz1 / nxyz;
       nx2 = (y3 - y2) * (z1 - z2) - (z3 - z2) * (y1 - y2);
       ny2 = (z3 - z2) * (x1 - x2) - (x3 - x2) * (z1 - z2);
       nz2 = (x3 - x2) * (y1 - y2) - (y3 - y2) * (x1 - x2);
       nxyz = sqrt(nx2 * nx2 + ny2 * ny2 + nz2 * nz2);
       nx2 = nx2 / nxyz;
       ny2 = ny2 / nxyz;
       nz2 = nz2 / nxyz;
       nx3 = (y4 - y3) * (z2 - z3) - (z4 - z3) * (y2 - y3);
       ny3 = (z4 - z3) * (x2 - x3) - (x4 - x3) * (z2 - z3);
       nz3 = (x4 - x3) * (y2 - y3) - (y4 - y3) * (x2 - x3);
       nxyz = sqrt(nx3 * nx3 + ny3 * ny3 + nz3 * nz3);
       nx3 = nx3 / nxyz;
       ny3 = ny3 / nxyz;
       nz3 = nz3 / nxyz;
       nx4 = (y1 - y4) * (z3 - z4) - (z1 - z4) * (y3 - y4);
       ny4 = (z1 - z4) * (x3 - x4) - (x1 - x4) * (z3 - z4);
       nz4 = (x1 - x4) * (y3 - y4) - (y1 - y4) * (x3 - x4);
       nxyz = sqrt(nx4 * nx4 + ny4 * ny4 + nz4 * nz4);
       nx4 = nx4 / nxyz;
       ny4 = ny4 / nxyz;
       nz4 = nz4 / nxyz;

if ((iu == 0 || iu == 1) && iv == 10) printf("3 jj iu,iv,x1,y1,z1= %d %d %d %le %le %le\n",jj,iu,iv,x1,y1,z1);

if (iu == 0 || iu == 1)
{
   glBegin(GL_LINES);
      glColor3f(1.0, 0.0, 0.0 ); 
      glVertex3f (x1, y1, z1);
      glVertex3f (x2, y2, z2);
     glEnd(); 
}
*/
/*

        glColor3f(0.5 * (iu - 0) / 100, 0.4 * (1.0 - iu / 100), 0.7 );
   glBegin(GL_POLYGON);
      glNormal3f (nx1,ny1,nz1);
      glVertex3f (x1, y1, z1);
      glNormal3f (nx2,ny2,nz2);
      glVertex3f (x2, y2, z2);
      glNormal3f (nx3,ny3,nz3);
      glVertex3f (x3, y3, z3);
      glNormal3f (nx4,ny4,nz4);
      glVertex3f (x4, y4, z4);
     glEnd(); 

     }     
   }   
*/
//goto L10;
/*
   iuu = 0;
   if(jj == icurve-2) iuu=1;
   du = 1.0 / 1.0;
   dv = 2.0 * pi / 20.0;
   for (iu = 0; iu <= iuu; iu ++)
   {
       u01 = iu * du;
	   u02=u01*u01;
	   u03=u02*u01;
	   u04=u02*u02;
	   u05=u03*u02;
       fprintf(faceout1,"curve -d 3\n");
	   for (iv = 0; iv <= 20; iv++)
       {
       vv0 = iv * dv;
       x1=cx0[0]+cx0[1]*u01+cx0[2]*u02+cx0[3]*u03+cx0[4]*u04+cx0[5]*u05+csx[0]*sin(vv0)+ccx[0]*cos(vv0);
       y1=cy0[0]+cy0[1]*u01+cy0[2]*u02+cy0[3]*u03+cy0[4]*u04+cy0[5]*u05+csy[0]*sin(vv0)+ccy[0]*cos(vv0);
       z1=cz0[0]+cz0[1]*u01+cz0[2]*u02+cz0[3]*u03+cz0[4]*u04+cz0[5]*u05+csz[0]*sin(vv0)+ccz[0]*cos(vv0);
       for (i=1; i <= 11; i++)
       {
         prod1(&i,&u01,&un);
         x1=x1+csx[i]*un*sin(vv0)+ccx[i]*un*cos(vv0);
         y1=y1+csy[i]*un*sin(vv0)+ccy[i]*un*cos(vv0);
         z1=z1+csz[i]*un*sin(vv0)+ccz[i]*un*cos(vv0);
       }
	   fprintf(faceout1,"-p %le %le %le\n",x1,y1,z1);
     } 
     fprintf(faceout1,";\n");  
   }  


   if(jj < icurve-2)
   {
	   jj=jj+1;
	   goto L20;
   }


L10: printf("");
*/
   fclose(in1);
   fclose(in2);
   fclose(out);
   fclose(out2);
   fclose(faceout1);
    
    delete [] bbx;
	for(i=0; i<n; i++)
    {
		delete aax[i];
	}
    delete [] aax;
    delete [] bbx0;
	for(i=0; i<n0; i++)
    {
		delete aax0[i];
	}

    delete [] aax0;   glEndList();

}


void display(void)
{

   glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
   glPushMatrix();

   glEnable(GL_COLOR_MATERIAL);
   glEnable (GL_LIGHTING);
   glShadeModel (GL_SMOOTH);

   glColor3f(0.5, 0.0, 1.0);
   glTranslatef(0.0, 0.0, 0.0);
   glPushMatrix();
   glScalef(0.1, 0.1, 0.1);
/* p3.gif  */
//   glRotatef(-70.0, 0.0, 0.0, 1.0);

/* p3x.gif	*/  
 
   glRotatef(45.0,45.0, 45.0, 1.0);


   glCallList(startList);
   glPopMatrix();

   glFlush();
}

void reshape (int w, int h)
{
   glViewport(0, 0, (GLsizei) w, (GLsizei) h);
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   if (w <= h)
      glOrtho(-2.5, 2.5, -2.5*(GLfloat)h/(GLfloat)w,
         2.5*(GLfloat)h/(GLfloat)w, -10.0, 10.0);
   else
      glOrtho(-2.5*(GLfloat)w/(GLfloat)h,
         2.5*(GLfloat)w/(GLfloat)h, -2.5, 2.5, -10.0, 10.0);
   glMatrixMode(GL_MODELVIEW);
   glLoadIdentity();
}

void keyboard(unsigned char key, int x, int y)
{
   switch (key) {
      case 27:
         exit(0);
         break;
   }
}

int main(int argc, char** argv)
{
   glutInit(&argc, argv);
   glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB | GLUT_DEPTH);
   glutInitWindowSize(500, 500); 
   glutInitWindowPosition(100, 100);
   glutCreateWindow(argv[0]);
   init();
   glutDisplayFunc(display); 
   glutReshapeFunc(reshape);
   glutKeyboardFunc(keyboard);
   glutMainLoop();
   return 0;
}




/*
   int cons(double h0,double h01,double h02,double h1,double h11,double h12,double cxyz[6])
//   double cxyz[6];
   {
	   double matrixa[6][6],columnb[6],f1,f2,f3,f4,f5,f6;
       matrixa[0][0] = 1.0;
       matrixa[0][1] = 0.0;
       matrixa[0][2] = 0.0;
       matrixa[0][3] = 0.0;
       matrixa[0][4] = 0.0;
       matrixa[0][5] = 0.0;
       matrixa[1][0] = 0.0;
       matrixa[1][1] = 1.0;
       matrixa[1][2] = 0.0;
       matrixa[1][3] = 0.0;
       matrixa[1][4] = 0.0;
       matrixa[1][5] = 0.0;
       matrixa[2][0] = 0.0;
       matrixa[2][1] = 0.0;
       matrixa[2][2] = 2.0;
       matrixa[2][3] = 0.0;
       matrixa[2][4] = 0.0;
       matrixa[2][5] = 0.0;
       matrixa[3][0] = 1.0;
       matrixa[3][1] = 1.0;
       matrixa[3][2] = 1.0;
       matrixa[3][3] = 1.0;
       matrixa[3][4] = 1.0;
       matrixa[3][5] = 1.0;
       matrixa[4][0] = 0.0;
       matrixa[4][1] = 1.0;
       matrixa[4][2] = 2.0;
       matrixa[4][3] = 3.0;
       matrixa[4][4] = 4.0;
       matrixa[4][5] = 5.0;
       matrixa[5][0] = 0.0;
       matrixa[5][1] = 0.0;
       matrixa[5][2] = 2.0;
       matrixa[5][3] = 6.0;
       matrixa[5][4] = 12.0;
       matrixa[5][5] = 20.0;
       columnb[0] = h0;
       columnb[1] = h01; 
       columnb[2] = h02;
       columnb[3] = h1;
       columnb[4] = h11;
       columnb[5] = h12;
       ciggj(matrixa,6,columnb);
       cxyz[0] = columnb[0];
       cxyz[1] = columnb[1];
       cxyz[2] = columnb[2];
       cxyz[3] = columnb[3];
       cxyz[4] = columnb[4];
       cxyz[5] = columnb[5];
       f1 = cxyz[0] - h0;
       f2 = cxyz[1] - h01;
       f3 = 2.0 * cxyz[2] - h02; 
       f4 = cxyz[0] + cxyz[1] + cxyz[2] + cxyz[3] + cxyz[4] + cxyz[5] - h1;
       f5 =  cxyz[1] + 2.0 * cxyz[2] + 3.0 * cxyz[3] + 4.0 * cxyz[4] + 5.0 * cxyz[5] - h11;
       f6 = 2.0 * cxyz[2] + 6.0 * cxyz[3] + 12.0 * cxyz[4] + 20.0 * cxyz[5] - h12; 
	   printf("f1-6= %le %le %le %le %le %le\n",f1,f2,f3,f4,f5,f6);
	   return(1);
   }


   int sincos(double ax,double bx,double cx,double dx,double ai0,double ai1,double ai2,
	   double ai3,double ai4,double ai5,double di2, double di4,double di6,double csc[11])
   {
	   int i;
	   double d1,d2,d3,d4,d5,d6,matrixa[6][6],columnb[6];
     d1 = - (24.0 * bx * (15.0 * ai0 + 8.0 * ai1 + 1.5 * ai2 - 15.0 * ai3 
          + 7.0 * ai4 - ai5) * di2 + cx * ai2 * di4 + dx * ai0 * di6);
     d2 = - (120.0 * bx * (- 6.0 * ai0 - 3.0 * ai1 - 0.5 * ai2 + 6.0 * ai3
          - 3.0 * ai4 + 0.5 * ai5) * di2 + 6.0 * cx * (- 10.0 * ai0 - 6.0 
          * ai1 - 1.5 * ai2 + 10.0 * ai3 - 4.0 * ai4 + 0.5 * ai5) * di4
          + dx * ai1 * di6);
     d3 = - (12.0 * cx * (15.0 * ai0 + 8.0 * ai1 + 1.5 * ai2 - 15.0 * ai3
          + 7.0 * ai4 - ai5) * di4 + 0.5 * dx * ai2 * di6);
     d4 = - (20.0 * cx * (- 6.0 * ai0 - 3.0 * ai1 - 0.5 * ai2 + 6.0 * ai3
          - 3.0 * ai4 + 0.5 * ai5) * di4 + dx * (- 10.0 * ai0 - 6.0 * ai1 
          - 1.5 * ai2 + 10.0 * ai3 - 4.0 * ai4 + 0.5 * ai5) * di6);
     d5 = - dx *(15.0 * ai0 + 8.0 * ai1 + 1.5 * ai2 - 15.0 * ai3 + 7.0 * ai4 
          - ai5) * di6;
     d6 = - dx *(- 6.0 * ai0 - 3.0 * ai1 - 0.5 * ai2 + 6.0 * ai3 - 3.0 * ai4 
          + 0.5 * ai5) * di6;
     matrixa[0][0] = 72.0 * (10.0 * ax + bx * di2);
     matrixa[0][1] = 24.0 * 8.0 * bx * di2;
     matrixa[0][2] = 24.0 * 15.0 * bx * di2;
     matrixa[0][3] = 24.0 * 24.0 * bx * di2;
     matrixa[0][4] = 24.0 * 35.0 * bx * di2;
     matrixa[0][5] = 24.0 * 48.0 * bx * di2;
     matrixa[1][0] = - 6.0 * (60.0 * bx * di2 + cx * di4);
     matrixa[1][1] = 18.0 * (280.0 * ax - 40.0 * bx * di2 - cx * di4);
     matrixa[1][2] = - 12.0 * (100.0 * bx * di2 + 3.0 * cx * di4);
     matrixa[1][3] = - 60.0 * (30.0 * bx * di2 + cx * di4);
     matrixa[1][4] = - 90.0 * (28.0 * bx * di2 + cx * di4);
     matrixa[1][5] = - 42.0 * (80.0 * bx * di2 + 3.0 * cx * di4); 
     matrixa[2][0] = 36.0 * (10.0 * bx * di2 + cx * di4);
     matrixa[2][1] = 96.0 * cx * di4;
     matrixa[2][2] = 180.0 * (112.0 * ax + cx * di4);
     matrixa[2][3] = 12.0 * 24.0 * cx * di4;
     matrixa[2][4] = 12.0 * 35.0 * cx * di4;
     matrixa[2][5] = 12.0 * 48.0 * cx * di4;
     matrixa[3][0] = - (60.0 * cx * di4 + dx * di6);
     matrixa[3][1] = 3.0 * (280.0 * bx * di2 - 40.0 * cx * di4 - dx * di6);
     matrixa[3][2] = - 2.0 * (100.0 * cx * di4 + 3.0 * dx * di6);
     matrixa[3][3] = 10.0 * (6048.0 * ax - 30.0 * cx * di4 - dx * di6);
     matrixa[3][4] = - 15.0 * (28.0 * cx * di4 + dx * di6);
     matrixa[3][5] = - 7.0 * (80.0 * cx * di4 + 3.0 * dx * di6);
     matrixa[4][0] = 3.0 * (10.0 * cx * di4 + dx * di6);
     matrixa[4][1] = 8.0 * dx * di6;
     matrixa[4][2] = 15.0 * (112.0 * bx * di2 + dx * di6);
     matrixa[4][3] = 24.0 * dx * di6;
     matrixa[4][4] = 35.0 * (4320.0 * ax + dx * di6);
     matrixa[4][5] = 48.0 * dx * di6;
     matrixa[5][0] = - 3.0 * dx * di6;
     matrixa[5][1] = 6.0 * (7.0 * cx * di4 - dx * di6);
     matrixa[5][2] = - 10.0 * dx * di6;
     matrixa[5][3] = 3.0 * (1008.0 * bx * di2 - 5.0 * dx * di6);
     matrixa[5][4] = - 21.0 * dx * di6;
     matrixa[5][5] = 28.0 * (11880.0 * ax - dx * di6);
     columnb[0] = d1;
     columnb[1] = d2; 
     columnb[2] = d3;
     columnb[3] = d4;
     columnb[4] = d5;
     columnb[5] = d6;
     ciggj(matrixa,6,columnb);
     csc[6] = columnb[0];
     csc[7] = columnb[1];
     csc[8] = columnb[2];
     csc[9] = columnb[3];
     csc[10] = columnb[4];
     csc[11] = columnb[5];

	 csc[0] = ai0;
     csc[1] = ai1;
     csc[2] = 0.5 * ai2;
     csc[3] = - 10.0 * ai0 - 6.0 * ai1 - 1.5 * ai2 + 10.0 * ai3 - 4.0 * ai4 + 0.5
              * ai5;
     csc[4] = 15.0 * ai0 + 8.0 * ai1 + 1.5 * ai2 - 15.0 * ai3 + 7.0 * ai4 - ai5;
     csc[5] = - 6.0 * ai0 - 3.0 * ai1 - 0.5 * ai2 + 6.0 * ai3 - 3.0 * ai4 + 0.5 
              * ai5;
     for (i = 6; i <= 11; i++)
     {
       csc[3] = csc[3] - (0.5 * i * i - 4.5 * i + 10.0) * csc[i];
       csc[4] = csc[4] + (i * i - 8.0 * i + 15.0) * csc[i];
       csc[5] = csc[5] - 0.5 * (i * i - 7.0 * i + 12.0) * csc[i]; 
     }
	 return(1);
   }

*/
int ciggj(double **a,int n,double b[jb])      
{
	int i, j, k, is, u, v;
	int *js = new int [n];
	double d, t;
	for (k=0;k<=n-1;k++)
	{
		d=0.0;
		for (i=k;i<=n-1;i++)
		for (j=k;j<=n-1;j++)
			{t=fabs(a[i][j]);			
			if(t>d) {d=t;js[k]=j;is=i;}
			}
			if(d+1.0==1.0)
				{delete js;printf("fail\n");return(0);}
			if (is!=k)
				{for(j=k;j<=n-1;j++)
					{u=k*n+j;v=is*n+j;
					 t=a[k][j];a[k][j]=a[is][j];a[is][j]=t;
					 }
				  t=b[k];b[k]=b[is];b[is]=t;
				 }
			if (js[k]!=k)
			for (i=0;i<=n-1;i++)
				{u=i*n+k;v=i*n+js[k];
				 t=a[i][k];a[i][k]=a[i][js[k]];a[i][js[k]]=t;
				}
			t=a[k][k];
			for (j=k+1;j<=n-1;j++)
				{u=k*n+j;
				 if (a[k][j]!=0.0)a[k][j]=a[k][j]/t;
				}
			b[k]=b[k]/t;
			for (j=k+1;j<=n-1;j++)
				{u=k*n+j;
				 if(a[k][j]!=0.0)
				 	{for (i=0;i<=n-1;i++)
				 		{v=i*n+k;
				 		 if((i!=k) && (a[i][k]!=0.0))
				 		 	{is=i*n+j;
				 		 	 a[i][j]=a[i][j]-a[i][k]*a[k][j];
				 		 	}
				 		 }
				 	}
				 }
				 for(i=0;i<=n-1;i++)
				 	{u=i*n+k;
				 	 if((i!=k) && (a[i][k]!=0.0))
				 	 	b[i]=b[i]-a[i][k]*b[k];
				 	 }
			}
			for (k=n-1;k>=0;k--)
				if(k!=js[k])
					{t=b[k];b[k]=b[js[k]];b[js[k]]=t;}
				delete js;
				return(1);
			}

int Bsurface(double *xcentre, double *ycentre, double *zcentre, double **ai,double **bi)      
{
	int i,j,k,ns;
	double u1,u12,u13,u2,u22,u23,u3,u32,u33,D,D10,D11,D12,D13,D20,D21,D22,D23,D30,
		D31,D32,D33,u0,v0,s0vx,s1vx,s2vx,s3vx,s0vy,s1vy,s2vy,s3vy,s0vz,
		s1vz,s2vz,s3vz,f0vx,f1vx,f2vx,f3vx,f0vy,f1vy,f2vy,f3vy,f0vz,f1vz,
		f2vz,f3vz,ddu,du,dv,x1,y1,z1;
	ns=50;
	du=1.0/ns;
	dv=2.0*pi/ns;
	ddu=1.0/3.0;
	u1=ddu;
	u12=u1*u1;
	u13=u1*u12;
	u2=2.0*ddu;
	u22=u2*u2;
	u23=u2*u22;
	u3=1.0;
	u32=u3*u3;
	u33=u3*u32;
    D=u1*u22*u33+u12*u23*u3+u13*u2*u32-u3*u22*u13-u32*u23*u1-u33*u2*u12;
	D10=u22*u13+u32*u23+u33*u12-u22*u33-u12*u23-u13*u32;
	D11=u22*u33-u32*u23;
	D12=u13*u32-u33*u12;
	D13=u12*u23-u22*u13;
	D20=u3*u13+u23*u1+u33*u2-u23*u3-u1*u33-u13*u2;
	D21=u23*u3-u33*u2;
	D22=u1*u33-u3*u13;
	D23=u13*u2-u23*u1;
	D30=u3*u22+u1*u32+u12*u2-u1*u22-u12*u3-u2*u32;
	D31=u2*u32-u3*u22;
	D32=u12*u3-u1*u32;
	D33=u1*u22-u12*u2;
	for(i=0;i<=ns;i++)
	{
		u0=i*du;
 		fprintf(faceout1,"curve -d 3\n");
		for(k=0;k<=ns;k++)
		{
			v0=k*dv;
			s0vx=xcentre[0]+ai[0][0];
			s0vy=ycentre[0];
			s0vz=zcentre[0]+bi[0][0];
			s1vx=xcentre[1]+ai[1][0];
			s1vy=ycentre[1];
			s1vz=zcentre[1]+bi[1][0];
			s2vx=xcentre[2]+ai[2][0];
			s2vy=ycentre[2];
			s2vz=zcentre[2]+bi[2][0];
			s3vx=xcentre[3]+ai[3][0];
			s3vy=ycentre[3];
			s3vz=zcentre[3]+bi[3][0];
			for(j=1;j<=jb;j++)
			{
				s0vx=s0vx+ai[0][2*j-1]*cos(j*v0)
					+ai[0][2*j]*sin(j*v0);
				s0vy=s0vy+0.0;
				s0vz=s0vz+bi[0][2*j-1]*sin(j*v0)
					+bi[0][2*j]*cos(j*v0);
				s1vx=s1vx+ai[1][2*j-1]*cos(j*v0)
					+ai[1][2*j]*sin(j*v0);
				s1vy=s1vy+0.0;
				s1vz=s1vz+bi[1][2*j-1]*sin(j*v0)
					+bi[1][2*j]*cos(j*v0);
				s2vx=s2vx+ai[2][2*j-1]*cos(j*v0)
					+ai[2][2*j]*sin(j*v0);
				s2vy=s2vy+0.0;
				s2vz=s2vz+bi[2][2*j-1]*sin(j*v0)
					+bi[2][2*j]*cos(j*v0);
				s3vx=s3vx+ai[3][2*j-1]*cos(j*v0)
					+ai[3][2*j]*sin(j*v0);
				s3vy=s3vy+0.0;
				s3vz=s3vz+bi[3][2*j-1]*sin(j*v0)
					+bi[3][2*j]*cos(j*v0);
			}
			f0vx=s0vx;
			f0vy=s0vy;
			f0vz=s0vz;
			f1vx=(D10*s0vx+D11*s1vx+D12*s2vx+D13*s3vx)/D;
			f1vy=(D10*s0vy+D11*s1vy+D12*s2vy+D13*s3vy)/D;
			f1vz=(D10*s0vz+D11*s1vz+D12*s2vz+D13*s3vz)/D;
			f2vx=(D20*s0vx+D21*s1vx+D22*s2vx+D23*s3vx)/D;
			f2vy=(D20*s0vy+D21*s1vy+D22*s2vy+D23*s3vy)/D;
			f2vz=(D20*s0vz+D21*s1vz+D22*s2vz+D23*s3vz)/D;
			f3vx=(D30*s0vx+D31*s1vx+D32*s2vx+D33*s3vx)/D;
			f3vy=(D30*s0vy+D31*s1vy+D32*s2vy+D33*s3vy)/D;
			f3vz=(D30*s0vz+D31*s1vz+D32*s2vz+D33*s3vz)/D;
			x1=f0vx+u0*f1vx+u0*u0*f2vx+u0*u0*u0*f3vx;
			y1=f0vy+u0*f1vy+u0*u0*f2vy+u0*u0*u0*f3vy;
			z1=f0vz+u0*f1vz+u0*u0*f2vz+u0*u0*u0*f3vz;
		    fprintf(faceout1,"-p    %le %le %le\n",x1,y1,z1);
		}
		if(i==0)printf("%d %le\n",i,f0vy);
		fprintf(faceout1,";\n");
	}
	return(1);
}
/*
	int k,ic;
	for(ic=0; ic <= icurve-1; ic++)
	{
	}
*/
int Gellipse(int jb,int imic,double xcentreic,double ycentreic,double zcentreic,
			 double *angleic,double *xmic,double *ymic)
{
//  fitting of generalized ellipe using x=xc+a0+(j=1 to J)(a2j-1*cos(j*angle)+a2j*sin(j*angle) and z=zc+b0+(j=1 to J)(b2j-1)*sin(j*angle)+b2j*cos(j*angle)
//  The curves can be well approximated with the constant term a0.Therefore,the following algorithm is very good.
/*
	int i,n;
	double *bbx = new double [2*jb+1];
	double **aax = new double* [2*jb+1];
	for(i=0; i<2*jb+1; i++)
    {
		aax[i] = new double [2*jb+1];
	}
	n=2*jb+1;
		bbx[0]=0.0;
		aax[0][0]=0.0;
		for(i=0; i <= imic-1; i++)
		{
			bbx[0]=bbx[0]+(xmic[i]-xcentreic);
			aax[0][0]=aax[0][0]+1.0;
		}
		for(j=1;j<=jb;j++)
		{
			aax[0][2*j-1]=0.0;
			aax[0][2*j]=0.0;
			for(i=0; i <= imic-1; i++)
			{
				aax[0][2*j-1]=aax[0][2*j-1]+cos(j*angleic[i]);
				aax[0][2*j]=aax[0][2*j]+sin(j*angleic[i]);
			}
		}
		for(j=1;j<n;j++)
		{
			bbx[j]=0.0;
			for(k=0;k<n;k++)aax[j][k]=0.0;
		}
		for(k=1;k<=jb;k++)
		{
			for(i=0; i <= imic-1; i++)
			{
				bbx[k]=bbx[k]+(xmic[i]-xcentreic)*cos(k*angleic[i]);
				bbx[jb+k]=bbx[jb+k]+(xmic[i]-xcentreic)*sin(k*angleic[i]);
				aax[k][0]=aax[k][0]+cos(k*angleic[i]);
				aax[jb+k][0]=aax[jb+k][0]+sin(k*angleic[i]);
			}
			for(j=1;j<=jb;j++)
			{
				for(i=0; i <= imic-1; i++)
				{				
					aax[k][2*j-1]=aax[k][2*j-1]+cos(j*angleic[i])*
						cos(k*angleic[i]);			
					aax[k][2*j]=aax[k][2*j]+sin(j*angleic[i])*
						cos(k*angleic[i]);
					aax[jb+k][2*j-1]=aax[jb+k][2*j-1]+cos(j*angleic[i])*
						sin(k*angleic[i]);			
					aax[jb+k][2*j]=aax[jb+k][2*j]+sin(j*angleic[i])*
						sin(k*angleic[i]);
				}
			}
		}
		for(i=0;i<n;i++)
		{
			for(j=0;j<n;j++)
			{
			}
		}
	    ciggj(aax,n,bbx);
		for( k=0;k<=n-1;k++)aiic[k]=bbx[k];
	}

		bbx[0]=0.0;
		aax[0][0]=0.0;
		for(i=0; i <= imic-1; i++)
		{
			bbx[0]=bbx[0]+(zmic[i]-zcentreic);
			aax[0][0]=aax[0][0]+1.0;
		}
		for(j=1;j<=jb;j++)
		{
			aax[0][2*j-1]=0.0;
			aax[0][2*j]=0.0;
			for(i=0; i <= im[ic]-1; i++)
			{
				aax[0][2*j-1]=aax[0][2*j-1]+sin(j*angleic[i]);
				aax[0][2*j]=aax[0][2*j]+cos(j*angleic[i]);
			}
		}
		for(j=1;j<n;j++)
		{
			bbx[j]=0.0;
			for(k=0;k<n;k++)aax[j][k]=0.0;
		}
		for(k=1;k<=jb;k++)
		{
			for(i=0; i <= imic-1; i++)
			{
				bbx[k]=bbx[k]+(zmic[i]-zcentreic)*sin(k*angleic[i]);
				bbx[jb+k]=bbx[jb+k]+(zmic[i]-zcentreic)*cos(k*angleic[i]);
				aax[k][0]=aax[k][0]+sin(k*angleic[i]);
				aax[jb+k][0]=aax[jb+k][0]+cos(k*angleic[i]);
			}
			for(j=1;j<=jb;j++)
			{
				for(i=0; i <= imic-1; i++)
				{				
					aax[k][2*j-1]=aax[k][2*j-1]+sin(j*angleic[i])*
						sin(k*angleic[i]);			
					aax[k][2*j]=aax[k][2*j]+cos(j*angleic[i])*
						sin(k*angleic[i]);
					aax[jb+k][2*j-1]=aax[jb+k][2*j-1]+sin(j*angleic[i])*
						cos(k*angleic[i]);			
					aax[jb+k][2*j]=aax[jb+k][2*j]+cos(j*angleic[i])*
						cos(k*angleic[i]);
				}
			}
		}
		for(i=0;i<n;i++)
		{
			for(j=0;j<n;j++)
			{
			}
		}
	    ciggj(aax,n,bbx);
		for( k=0;k<=n-1;k++)bi[ic][k]=bbx[k];
*/		return(1);
}
/*

	for(ic=0; ic <= icurve-1; ic++)
	{
		for( j=0;j<n;j++)
		{ 
			bbx[j]=0.0;
			for( k=0;k<n;k++)
			{ 
				aax[j][k]=0.0;
			}
		}
		for(i=0; i <= im[ic]-1; i++)
		{
			aax[0][0]=aax[0][0]+cos(angle[ic][i])*cos(angle[ic][i]);
			aax[0][1]=aax[0][1]+sin(angle[ic][i])*cos(angle[ic][i]);
			bbx[0]=bbx[0]+(xm[ic][i]-xcentre[ic])*cos(angle[ic][i]);
			aax[1][0]=aax[1][0]+cos(angle[ic][i])*sin(angle[ic][i]);
			aax[1][1]=aax[1][1]+sin(angle[ic][i])*sin(angle[ic][i]);
			bbx[1]=bbx[1]+(xm[ic][i]-xcentre[ic])*sin(angle[ic][i]);
		}
	    ciggj(aax,n,bbx);
		for( k=0;k<=n-1;k++)ai[ic][k]=bbx[k];
	}

	for(ic=0; ic <= icurve-1; ic++)
	{
		for( j=0;j<n;j++)
		{ 
			bbx[j]=0.0;
			for(k=0;k<n;k++)
			{ 
				aax[j][k]=0.0;
			}
		}
		for(i=0; i <= im[ic]-1; i++)
		{
			aax[0][0]=aax[0][0]+sin(angle[ic][i])*sin(angle[ic][i]);
			aax[0][1]=aax[0][1]+sin(angle[ic][i])*cos(angle[ic][i]);
			bbx[0]=bbx[0]+(zm[ic][i]-zcentre[ic])*sin(angle[ic][i]);
			aax[1][0]=aax[1][0]+cos(angle[ic][i])*sin(angle[ic][i]);
			aax[1][1]=aax[1][1]+cos(angle[ic][i])*cos(angle[ic][i]);
			bbx[1]=bbx[1]+(zm[ic][i]-zcentre[ic])*cos(angle[ic][i]);
		}
	    ciggj(aax,n,bbx);
		for( k=0;k<=n-1;k++)bi[ic][k]=bbx[k];
	}
	for(ic=0; ic <= icurve-1; ic++)
	{
 		fprintf(out2,"curve -d 3\n");
		for(i=0; i <= im[ic]-1; i++)
		{
			x1=xcentre[ic]+ai[ic][0]*cos(angle[ic][i])+ai[ic][1]*sin(angle[ic][i]);
			y1=ycentre[ic];
			z1=zcentre[ic]+bi[ic][0]*sin(angle[ic][i])+bi[ic][1]*cos(angle[ic][i]);
			fprintf(out2,"-p    %le %le %le\n",x1,y1,z1);
		}
		fprintf(out2,";\n");
	}
*/
//  fitting of generalized ellipe using x=xc+(j=1 to J)(a2j-2)*cos(j*angle)+(a2j-1)*sin(j*angle) and z=zc+(j=1 to J)(b2j-2)*sin(j*angle)+(b2j-1)*cos(j*angle)
//  The curves cannot be well approximated without the constant term a0.Therefore, the following algorithm is not good.
/*
	for(ic=0; ic <= icurve-1; ic++)
	{
		for(j=0;j<n0;j++)
		{
			bbx0[j]=0.0;
			for(k=0;k<n0;k++)aax0[j][k]=0.0;
		}
		for(k=1;k<=jb;k++)
		{
			for(i=0; i <= im[ic]-1; i++)
			{
				bbx0[k-1]=bbx0[k-1]+(xm[ic][i]-xcentre[ic])*cos(k*angle[ic][i]);
				bbx0[jb+k-1]=bbx0[jb+k-1]+(xm[ic][i]-xcentre[ic])*sin(k*angle[ic][i]);
			}
			for(j=1;j<=jb;j++)
			{
				for(i=0; i <= im[ic]-1; i++)
				{				
					aax0[k-1][2*j-2]=aax0[k-1][2*j-2]+cos(j*angle[ic][i])*
						cos(k*angle[ic][i]);			
					aax0[k-1][2*j-1]=aax0[k-1][2*j-1]+sin(j*angle[ic][i])*
						cos(k*angle[ic][i]);
					aax0[jb+k-1][2*j-2]=aax0[jb+k-1][2*j-2]+cos(j*angle[ic][i])*
						sin(k*angle[ic][i]);			
					aax0[jb+k-1][2*j-1]=aax0[jb+k-1][2*j-1]+sin(j*angle[ic][i])*
						sin(k*angle[ic][i]);
				}
			}
		}
	    ciggj(aax0,n0,bbx0);
		for( k=0;k<=n0-1;k++)ai[ic][k]=bbx0[k];
	}

	for(ic=0; ic <= icurve-1; ic++)
	{
		for(j=0;j<n0;j++)
		{
			bbx0[j]=0.0;
			for(k=0;k<n0;k++)aax0[j][k]=0.0;
		}
		for(k=1;k<=jb;k++)
		{
			for(i=0; i <= im[ic]-1; i++)
			{
				bbx0[k-1]=bbx0[k-1]+(zm[ic][i]-zcentre[ic])*sin(k*angle[ic][i]);
				bbx0[jb+k-1]=bbx0[jb+k-1]+(zm[ic][i]-zcentre[ic])*cos(k*angle[ic][i]);
			}
			for(j=1;j<=jb;j++)
			{
				for(i=0; i <= im[ic]-1; i++)
				{				
					aax0[k-1][2*j-2]=aax0[k-1][2*j-2]+sin(j*angle[ic][i])*
						sin(k*angle[ic][i]);			
					aax0[k-1][2*j-1]=aax0[k-1][2*j-1]+cos(j*angle[ic][i])*
						sin(k*angle[ic][i]);
					aax0[jb+k-1][2*j-2]=aax0[jb+k-1][2*j-2]+sin(j*angle[ic][i])*
						cos(k*angle[ic][i]);			
					aax0[jb+k-1][2*j-1]=aax0[jb+k-1][2*j-1]+cos(j*angle[ic][i])*
						cos(k*angle[ic][i]);
				}
			}
		}
	    ciggj(aax0,n0,bbx0);
		for( k=0;k<=n0-1;k++)bi[ic][k]=bbx0[k];
	}

	for(ic=0; ic <= icurve-1; ic++)
	{
 		fprintf(out2,"curve -d 3\n");
		for(i=0; i <= im[ic]-1; i++)
		{
			x1=xcentre[ic];
			y1=ycentre[ic];
			z1=zcentre[ic];
			for(int j=1;j<=jb;j++)
			{        

                x1=x1+ai[ic][2*j-2]*cos(j*angle[ic][i])
					+ai[ic][2*j-1]*sin(j*angle[ic][i]);
				y1=y1+0.0;
				z1=z1+bi[ic][2*j-2]*sin(j*angle[ic][i])
					+bi[ic][2*j-1]*cos(j*angle[ic][i]);
			}
		    fprintf(out2,"-p    %le %le %le\n",x1,y1,z1);
		}
		fprintf(out2,";\n");
	}
*/