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

#define jb  20
#define icurve 39



int ciggj(double **,int,double [jb]);
int Bsurface(int, double *, double *, double *, double **,double **,int &, double *, 
			 double *, double *, double *, double *, double *);      
int Gellipse(int,int,double,double,double,double *,double *,double *,double *,double *);
void CentreAngle(int,double &,double &,double &,double *,double *,double *, double *,
				 double *,double *);
void u0if(int &, double &,double &);
int ciggj4(double [][4],int,double [4]);

FILE *in1,*out,*out2,*faceout1,*out3;

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
		anx[100][10],any[100][10],anz[100][10],dv,xi,ff,
		vv0,vv1,nxi[100][10],nyi[100][10],nzi[100][10],
		x0,x00,y0,y00,al,all,bs,bss,fmin,dx0,dy0,ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,
		h0,h01,h02,h1,h11,h12,w1,w2,cx0[6],cy0[6],cz0[6],ai0,ai1,ai2,ai3,ai4,ai5,
		csx[12],csy[12],csz[12],ccx[12],ccy[12],ccz[12],di2,di4,di6,du,u01,u02,u03,
		u04,u05,u11,u12,u13,u14,u15,un,nx1,nx2,nx3,nx4,ny1,ny2,ny3,ny4,nz1,nz2,
		nz3,nz4,nxyz,txi[100],tyi[100],tzi[100],tx1,ty1,tz1,xcentre[50],ycentre[50],
		zcentre[50],angle[50][200],anglem,a0[70],b0[70],a0u,a0b,b0u,b0b,dpi,
		ai[50][70],bi[50][70];
	char szLine[256];
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
    out = fopen("result.out","wt");
    out2 = fopen("out2.mel","wt");
    out3 = fopen("out3.mel","wt");
    faceout1 = fopen("PDE_Surface.mel","wt");

	for(i=0; i <= icurve-1; i++)
	{
		fscanf(in1,"%d\n",&im[i]);
// 		fprintf(out2,"curve -d 3\n");
       for(j=0; j <= im[i]-1; j++)
		{
            fscanf(in1,"%le %le %le\n",&xm[i][j],&ym[i][j],&zm[i][j]);
//		    fprintf(out2,"-p   %le %le %le\n",xm[i][j],ym[i][j],zm[i][j]);
		}
		fgets(szLine,256,in1);
//		fprintf(out2,"-p   %le %le %le\n",xm[i][0],ym[i][0],zm[i][0]);
//		fprintf(out2,";\n");
	}

	int k,ic,imic;
	double xcentreic,ycentreic,zcentreic,angleic[200],xmic[200],ymic[200],zmic[200],aiic[70],
		biic[70];

//  Calculate the center of the ith curve and determine the angle of each curve point.
	for(ic=0; ic<= icurve-1; ic++)
	{
		imic=im[ic];
		for(j=0;j<=imic-1;j++)
		{
			xmic[j]=xm[ic][j];
			ymic[j]=ym[ic][j];
			zmic[j]=zm[ic][j];
		}
		CentreAngle(imic,xcentreic,ycentreic,zcentreic,angleic,xmic,ymic,zmic,aiic,biic);
		xcentre[ic]=xcentreic;
		ycentre[ic]=ycentreic;
		zcentre[ic]=zcentreic;
		for(j=0;j<=imic-1;j++)
		{
		angle[ic][j]=angleic[j];
		}
	}

//  Order the data according to the angle from 0 to 2*pi.
	int j0;
	double anglemin, anglemax,angleg[50][200],xmg[50][200],ymg[50][200],zmg[50][200];
	for(ic=0; ic<= icurve-1; ic++)
	{
		anglemin=0.0;
		anglemax=1.e10;
		for(k=0;k<=im[ic]-1;k++)
		{
			for(j=0;j<=im[ic]-1;j++)
			{
				if(angle[ic][j]>anglemin && angle[ic][j]<anglemax)
				{
					j0=j;
					anglemax=angle[ic][j0];
				}
			}
			angleg[ic][k]=anglemax;
			xmg[ic][k]=xm[ic][j0];
			ymg[ic][k]=ym[ic][j0];
			zmg[ic][k]=zm[ic][j0];
			anglemin=anglemax;
			anglemax=1.e10;
		}
	}
	for(ic=0; ic<= icurve-1; ic++)
	for(j=0;j<=im[ic]-1;j++)fprintf(out,"ic,j,angleg[ic][j]= %d %d %le\n",ic,j,angleg[ic][j]);

//   draw the original curves with ordered data

	for(i=0; i <= icurve-1; i++)
	{
 		fprintf(out2,"curve -d 3\n");
       for(j=0; j <= im[i]-1; j++)
		{
			fprintf(out2,"-p   %le %le %le\n",xmg[i][j],ymg[i][j],zmg[i][j]);
		}
		    fprintf(out2,"-p   %le %le %le\n",xmg[i][0],ymg[i][0],zmg[i][0]);
		fprintf(out2,";\n");
	}

//	icsize=4*(icurve1-icurve0+1);
	double v0,vj,cca[4][4],ccb[4],ccab[4][4],cp[3][4],xc[4],xmgb[300],ymgb[300],zmgb[300];
	int icurve0,icurve1,curvei,ip0,ip1,ip2,ip3,nv;
	curvei=30;
	icurve0=30;
	icurve1=35;
	nv=5;
		ip0=3;
		ip1=42;
		ip2=45;
		ip3=81;
		dv=1.0/3.0;
		for(i=0;i<=3;i++)
		{
			v0=i*dv;
			for(j=0;j<=3;j++)
			{
				u0if(j,v0,vj);
                cca[i][j]=vj;
				ccab[i][j]=cca[i][j];
			}
		}
		for(i=0;i<=3;i++)
		{
			ccb[i]=xmg[curvei][ip1-i];
		}
		ciggj4(cca,4,ccb);
		xc[0]=xmg[curvei][ip1];
		xc[1]=ccb[1]+2.0*ccb[2]+3.0*ccb[3];
		for(i=0;i<=3;i++)
		{
			for(j=0;j<=3;j++)
			{
 				cca[i][j]=ccab[i][j];
			}
		}
		for(i=0;i<=3;i++)
		{
			ccb[i]=xmg[curvei][ip2+i];
		}
		ciggj4(cca,4,ccb);
		xc[2]=xmg[curvei][ip2];
		xc[3]=ccb[1];
		cp[0][0]=xc[0];
		cp[0][1]=xc[1];
		cp[0][2]=3.0*(xc[2]-xc[1]-xc[0])-(xc[3]-xc[1]);
		cp[0][3]=xc[3]-xc[1]-2.0*(xc[2]-xc[1]-xc[0]);

		for(i=0;i<=3;i++)
		{
			ccb[i]=ymg[curvei][ip0+i];
		}
		ciggj4(cca,4,ccb);
		xc[2]=ymg[curvei][ip0];
		xc[3]=ccb[1];
		for(i=0;i<=3;i++)
		{
			for(j=0;j<=3;j++)
			{
 				cca[i][j]=ccab[i][j];
			}
		}
		for(i=0;i<=3;i++)
		{
			ccb[i]=ymg[curvei][ip3-i];
		}
		ciggj4(cca,4,ccb);
		xc[0]=ymg[curvei][ip3];
		xc[1]=ccb[1]+2.0*ccb[2]+3.0*ccb[3];
		cp[1][0]=xc[0];
		cp[1][1]=xc[1];
		cp[1][2]=3.0*(xc[2]-xc[1]-xc[0])-(xc[3]-xc[1]);
		cp[1][3]=xc[3]-xc[1]-2.0*(xc[2]-xc[1]-xc[0]);

		for(i=0;i<=3;i++)
		{
			ccb[i]=zmg[curvei][ip0+i];
		}
		ciggj4(cca,4,ccb);
		xc[2]=zmg[curvei][ip0];
		xc[3]=ccb[1];
		for(i=0;i<=3;i++)
		{
			for(j=0;j<=3;j++)
			{
 				cca[i][j]=ccab[i][j];
			}
		}
		for(i=0;i<=3;i++)
		{
			ccb[i]=zmg[curvei][ip3-i];
		}
		ciggj4(cca,4,ccb);
		xc[0]=zmg[curvei][ip3];
		xc[1]=ccb[1]+2.0*ccb[2]+3.0*ccb[3];
		cp[2][0]=xc[0];
		cp[2][1]=xc[1];
		cp[2][2]=3.0*(xc[2]-xc[1]-xc[0])-(xc[3]-xc[1]);
		cp[2][3]=xc[3]-xc[1]-2.0*(xc[2]-xc[1]-xc[0]);

		for(i=ip0;i<=ip1;i++)
		{
            xmgb[i-ip0]=xmg[curvei][i];
            ymgb[i-ip0]=ymg[curvei][i];
            zmgb[i-ip0]=zmg[curvei][i];
		}
		dv=1.0/nv;
		for(i=1;i<nv;i++)
		{
			xmgb[ip1-ip0+i]=0.0;
			ymgb[ip1-ip0+i]=0.0;
			zmgb[ip1-ip0+i]=0.0;
			v0=i*dv;
			for(j=0;j<=3;j++)
			{
				u0if(j,v0,vj);
                xmgb[ip1-ip0+i]=xmgb[ip1-ip0+i]+cp[0][j]*vj;
                ymgb[ip1-ip0+i]=ymgb[ip1-ip0+i]+cp[1][j]*vj;
                zmgb[ip1-ip0+i]=zmgb[ip1-ip0+i]+cp[2][j]*vj;
			}
		}
		for(i=ip2;i<=ip3;i++)
		{
            xmgb[ip1-ip0+nv+i-ip2+1]=xmg[curvei][i];
            ymgb[ip1-ip0+nv+i-ip2+1]=ymg[curvei][i];
            zmgb[ip1-ip0+nv+i-ip2+1]=zmg[curvei][i];
		}

		for(i=0;i<=3;i++)
		{
			ccb[i]=xmg[curvei][ip0+i];
		}
		ciggj4(cca,4,ccb);
		xc[2]=xmg[curvei][ip0];
		xc[3]=ccb[1];
		for(i=0;i<=3;i++)
		{
			for(j=0;j<=3;j++)
			{
 				cca[i][j]=ccab[i][j];
			}
		}
		for(i=0;i<=3;i++)
		{
			ccb[i]=xmg[curvei][ip3-i];
		}
		ciggj4(cca,4,ccb);
		xc[0]=xmg[curvei][ip3];
		xc[1]=ccb[1]+2.0*ccb[2]+3.0*ccb[3];
		cp[0][0]=xc[0];
		cp[0][1]=xc[1];
		cp[0][2]=3.0*(xc[2]-xc[1]-xc[0])-(xc[3]-xc[1]);
		cp[0][3]=xc[3]-xc[1]-2.0*(xc[2]-xc[1]-xc[0]);

		for(i=0;i<=3;i++)
		{
			ccb[i]=ymg[curvei][ip0+i];
		}
		ciggj4(cca,4,ccb);
		xc[2]=ymg[curvei][ip0];
		xc[3]=ccb[1];
		for(i=0;i<=3;i++)
		{
			for(j=0;j<=3;j++)
			{
 				cca[i][j]=ccab[i][j];
			}
		}
		for(i=0;i<=3;i++)
		{
			ccb[i]=ymg[curvei][ip3-i];
		}
		ciggj4(cca,4,ccb);
		xc[0]=ymg[curvei][ip3];
		xc[1]=ccb[1]+2.0*ccb[2]+3.0*ccb[3];
		cp[1][0]=xc[0];
		cp[1][1]=xc[1];
		cp[1][2]=3.0*(xc[2]-xc[1]-xc[0])-(xc[3]-xc[1]);
		cp[1][3]=xc[3]-xc[1]-2.0*(xc[2]-xc[1]-xc[0]);

		for(i=0;i<=3;i++)
		{
			ccb[i]=zmg[curvei][ip0+i];
		}
		ciggj4(cca,4,ccb);
		xc[2]=zmg[curvei][ip0];
		xc[3]=ccb[1];
		for(i=0;i<=3;i++)
		{
			for(j=0;j<=3;j++)
			{
 				cca[i][j]=ccab[i][j];
			}
		}
		for(i=0;i<=3;i++)
		{
			ccb[i]=zmg[curvei][ip3-i];
		}
		ciggj4(cca,4,ccb);
		xc[0]=zmg[curvei][ip3];
		xc[1]=ccb[1]+2.0*ccb[2]+3.0*ccb[3];
		cp[2][0]=xc[0];
		cp[2][1]=xc[1];
		cp[2][2]=3.0*(xc[2]-xc[1]-xc[0])-(xc[3]-xc[1]);
		cp[2][3]=xc[3]-xc[1]-2.0*(xc[2]-xc[1]-xc[0]);
		for(i=1;i<nv;i++)
		{
			k=ip1-ip0+nv+ip3-ip2+1+i;
			xmgb[k]=0.0;
			ymgb[k]=0.0;
			zmgb[k]=0.0;
			v0=i*dv;
			for(j=0;j<=3;j++)
			{
				u0if(j,v0,vj);
                xmgb[k]=xmgb[k]+cp[0][j]*vj;
                ymgb[k]=ymgb[k]+cp[1][j]*vj;
                zmgb[k]=zmgb[k]+cp[2][j]*vj;
			}
		}
/*
 		k=ip1-ip0+nv+ip3-ip2+nv;
		fprintf(out2,"curve -d 3\n");
        for(j=0; j <= k; j++)
		{
			fprintf(out2,"-p   %le %le %le\n",xmgb[j],ymgb[j],zmgb[j]);
		}
		fprintf(out2,";\n");
*/
exit(1);
//  Calculate the main axis a0[i] and minus axis b0[i] of the ith ellipse.
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
//  draw ellipses with x=xc+a*cos(angle) and z=zc+b*sin(angle)
/*
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
        x1=xcentre[i]+a0[i]*cos(angle[i][0]);
		y1=ycentre[i];
        z1=zcentre[i]+b0[i]*sin(angle[i][0]);
		fprintf(out2,"-p    %le %le %le\n",x1,y1,z1);
		fprintf(out2,";\n");
	}
*/
	for(i=0; i <= icurve-1; i++)
	{
// 		fprintf(out2,"curve -d 3\n");
        for(j=0; j <= im[i]-1; j++)
		{
            x1=xcentre[i]+a0[i]*cos(angleg[i][j]);
		    y1=ycentre[i];
            z1=zcentre[i]+b0[i]*sin(angleg[i][j]);
//		    fprintf(out2,"-p    %le %le %le\n",x1,y1,z1);
		}
        x1=xcentre[i]+a0[i]*cos(angleg[i][0]);
		y1=ycentre[i];
        z1=zcentre[i]+b0[i]*sin(angleg[i][0]);
//		fprintf(out2,"-p    %le %le %le\n",x1,y1,z1);
//		fprintf(out2,";\n");
	}
//  fitting of generalized ellipe using x=xc+a0+(j=1 to J)(a2j-1*cos(j*angle)+a2j*sin(j*angle) and z=zc+b0+(j=1 to J)(b2j-1)*sin(j*angle)+b2j*cos(j*angle)
//  The curves can be well approximated with the constant term a0.Therefore,the following algorithm is very good.
	for(ic=0; ic<=icurve-1;ic++)
	{
		imic=im[ic];
		xcentreic=xcentre[ic];
		ycentreic=ycentre[ic];
		zcentreic=zcentre[ic];
		for(j=0;j<=im[ic]-1;j++)
		{
			angleic[j]=angle[ic][j];
			xmic[j]=xm[ic][j];
			zmic[j]=zm[ic][j];
		}
		Gellipse(jb,imic,xcentreic,ycentreic,zcentreic,angleic,xmic,zmic,aiic,biic);
		for(j=0;j<=im[ic]-1;j++)
		{
			ai[ic][j]=aiic[j];
			bi[ic][j]=biic[j];
		}
	}
//  draw generalized ellipses with x=xc+a0+(j=1 to J)(a2j-1*cos(j*angle)+a2j*sin(j*angle) and z=zc+b0+(j=1 to J)(b2j-1)*sin(j*angle)+b2j*cos(j*angle)
    int ns=50;
	double anglei;
	du=1.0/ns;
	dv=2.0*pi/ns;
	for(ic=0; ic <= icurve-1; ic++)
	{
 		fprintf(faceout1,"curve -d 3\n");
 		fprintf(out2,"curve -d 3\n");
		for(i=0; i <= ns-1; i++)
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
		    fprintf(out2,"-p    %le %le %le\n",x1,y1,z1);
		}
		fprintf(faceout1,";\n");
		fprintf(out2,";\n");
			if(ic==0)printf("%d %le\n",ic,y1);
	}

	double *xcentrep = new double [4];
	double *ycentrep = new double [4];
	double *zcentrep = new double [4];
	double **aip = new double* [4];
	double **bip = new double* [4];
	double *xm0 = new double [71];
	double *ym0 = new double [71];
	double *zm0 = new double [71];
	double *xm1 = new double [71];
	double *ym1 = new double [71];
	double *zm1 = new double [71];
	int ic0,id;  
// id=2 stands for 4 curves without boundary tangents; 
// id=0 stand for the boundary tangent at ic0 represented 
//      with aip[ic0][] and the boundary curve at ic0 represented with aip[ic0][].
// id=1 stand for the boundary tangent at ic0+2 represented 
//      with aip[ic0+2][] and the boundary curve at ic0+2 represented with aip[ic0+2][].
	id=2;
	ic0=0;
	for(i=0; i<4; i++)
    {
		aip[i] = new double [70];
		bip[i] = new double [70];
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

	Bsurface(id,xcentrep, ycentrep, zcentrep, aip, bip, imic, xm0, ym0, zm0, xm1, ym1, zm1); 


	id=0;
	ic0=3;
	xcentrep[0]=xcentre[ic0];
	ycentrep[0]=ycentre[ic0];
	zcentrep[0]=zcentre[ic0];
	for(j=0;j<n;j++)
	{
		aip[0][j]=ai[ic0][j];
		bip[0][j]=bi[ic0][j];
	}
	for(ic=ic0+1;ic<=ic0+2;ic++)
	{
		xcentrep[ic-ic0+1]=xcentre[ic];
		ycentrep[ic-ic0+1]=ycentre[ic];
		zcentrep[ic-ic0+1]=zcentre[ic];
		for(j=0;j<n;j++)
		{
			aip[ic-ic0+1][j]=ai[ic][j];
			bip[ic-ic0+1][j]=bi[ic][j];
		}
	}

	Bsurface(id,xcentrep, ycentrep, zcentrep, aip, bip, imic, xm0, ym0, zm0, xm1, ym1, zm1); 

	id=0;
	ic0=5;
	xcentrep[0]=xcentre[ic0];
	ycentrep[0]=ycentre[ic0];
	zcentrep[0]=zcentre[ic0];
	for(j=0;j<n;j++)
	{
		aip[0][j]=ai[ic0][j];
		bip[0][j]=bi[ic0][j];
	}
	for(ic=ic0+1;ic<=ic0+2;ic++)
	{
		xcentrep[ic-ic0+1]=xcentre[ic];
		ycentrep[ic-ic0+1]=ycentre[ic];
		zcentrep[ic-ic0+1]=zcentre[ic];
		for(j=0;j<n;j++)
		{
			aip[ic-ic0+1][j]=ai[ic][j];
			bip[ic-ic0+1][j]=bi[ic][j];
		}
	}
	Bsurface(id,xcentrep, ycentrep, zcentrep, aip, bip, imic, xm0, ym0, zm0, xm1, ym1, zm1); 
    
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
    delete [] aax0;   

    delete [] xcentrep;
    delete [] ycentrep;
    delete [] zcentrep;
    delete [] xm0;
    delete [] ym0;
    delete [] zm0;
    delete [] xm1;
    delete [] ym1;
    delete [] zm1;
	for(i=0; i<4; i++)
    {
		delete aip[i];
		delete bip[i];
	}
    delete [] aip;
    delete [] bip;
	printf("Pass 1\n");

   fclose(in1);
   fclose(out);
   fclose(out2);
   fclose(out3);
   fclose(faceout1);
   
	printf("Pass 0\n");

	glEndList();

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

int ciggj4(double a[][4],int n,double b[4])      
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

int Bsurface(int id, double *xcentre, double *ycentre, double *zcentre, double **ai,double **bi,
			 int &imic, double *xm0, double *ym0, double *zm0, double *xm1, double *ym1, 
			 double *zm1)      
{
	int i,j,k,ns,i0;
	double u1,u12,u13,u2,u22,u23,u3,u32,u33,D,D10,D11,D12,D13,D20,D21,D22,D23,D30,
		D31,D32,D33,s3v,u0,v0,s0vx,s1vx,s2vx,s3vx,s0vy,s1vy,s2vy,s3vy,s0vz,
		s1vz,s2vz,s3vz,f0vx,f1vx,f2vx,f3vx,f0vy,f1vy,f2vy,f3vy,f0vz,f1vz,
		f2vz,f3vz,ddu,du,dv,x1,y1,z1;
	if(imic == 0)i0=0;
	if(imic != 0)i0=1;
	ns=50;
	du=1.0/ns;
	dv=2.0*pi/ns;
	if(id==2)
	{
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
		for(i=i0;i<=ns;i++)
		{
			u0=i*du;
			fprintf(out3,"curve -d 3\n");
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
				fprintf(out3,"-p    %le %le %le\n",x1,y1,z1);
				if(i==0)
				{
					xm0[k]=f1vx;
					ym0[k]=f1vy;
					zm0[k]=f1vz;
				}
				if(i==ns)
				{
					xm1[k]=f1vx+2.0*f2vx+3.0*f3vx;
					ym1[k]=f1vy+2.0*f2vy+3.0*f3vy;
					zm1[k]=f1vz+2.0*f2vz+3.0*f3vz;
					fprintf(out,"1st tan=   %d %le %le %le\n",k,xm1[k],ym1[k],zm1[k]);
				}
			}
			if(i==0)printf("%d %le\n",i,f0vy);
			fprintf(faceout1,";\n");
			fprintf(out3,";\n");
		}
		imic=ns;
		for(k=0;k<=ns;k++)
		{
			fprintf(out,"1st k,x-y-zm1= %d %le %le %le\n",k,xm1[k],ym1[k],zm1[k]);
		}
	}

	if(id==0)
	{
		ddu=1.0/2.0;
		u1=ddu;	
		u12=u1*u1;
		u13=u1*u12;
		u2=1.0;
		u22=u2*u2;
		u23=u2*u22;
		D=u12*u23-u13*u22;
		D10=u13-u23;
		D11=u2*u13-u1*u23;
		D12=u23;
		D13=-u13;
		D20=u22-u12;
		D21=u1*u22-u2*u12;
		D22=-u22;
		D23=u12;
		for(j=0; j<=3;j++)printf("in j,ycentre[j]= %d %le\n",j,ycentre[j]);
		for(i=i0;i<=ns;i++)
		{
			u0=i*du;
//			fprintf(out2,"curve -d 3\n");
			if(i==0)fprintf(out3,"curve -d 3\n");
 			fprintf(faceout1,"curve -d 3\n");
			for(k=0;k<=ns;k++)
			{
				v0=k*dv;
				s0vx=xcentre[0]+ai[0][0];
				s0vy=ycentre[0];
				s0vz=zcentre[0]+bi[0][0];
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
				s1vx=xm1[k];
				s1vy=ym1[k];
				s1vz=zm1[k];
				f0vx=s0vx;
				f0vy=s0vy;
				f0vz=s0vz;
			if(i==0)fprintf(out,"2nd tan=   %d %le %le %le\n",k,xm1[k],ym1[k],zm1[k]);
				f1vx=s1vx;
				f1vy=s1vy;
				f1vz=s1vz;
				f2vx=(D10*s0vx+D11*s1vx+D12*s2vx+D13*s3vx)/D;
				f2vy=(D10*s0vy+D11*s1vy+D12*s2vy+D13*s3vy)/D;
				f2vz=(D10*s0vz+D11*s1vz+D12*s2vz+D13*s3vz)/D;
				f3vx=(D20*s0vx+D21*s1vx+D22*s2vx+D23*s3vx)/D;
				f3vy=(D20*s0vy+D21*s1vy+D22*s2vy+D23*s3vy)/D;
				f3vz=(D20*s0vz+D21*s1vz+D22*s2vz+D23*s3vz)/D;
				x1=f0vx+u0*f1vx+u0*u0*f2vx+u0*u0*u0*f3vx;
				y1=f0vy+u0*f1vy+u0*u0*f2vy+u0*u0*u0*f3vy;
				z1=f0vz+u0*f1vz+u0*u0*f2vz+u0*u0*u0*f3vz;
//				if(i==ns)fprintf(out,"i,x1,y1,z1= %d %le %le %le
				if(i==0)
				{
					fprintf(out3,"-p    %le %le %le\n",x1,y1,z1);
					fprintf(out,"1st tan=   %d %le %le %le\n",k,xm1[k],ym1[k],zm1[k]);
				}
//				fprintf(out2,"-p    %le %le %le\n",x1,y1,z1);
				fprintf(faceout1,"-p    %le %le %le\n",x1,y1,z1);
				if(i==0)
				{
					xm0[k]=f1vx;
					ym0[k]=f1vy;
					zm0[k]=f1vz;
				}
				if(i==ns)
				{
					xm1[k]=f1vx+2.0*f2vx+3.0*f3vx;
					ym1[k]=f1vy+2.0*f2vy+3.0*f3vy;
					zm1[k]=f1vz+2.0*f2vz+3.0*f3vz;
				}
			}
			if(i==0)printf("%d %le\n",i,f0vy);
//			fprintf(out2,";\n");
			if(i==ns)fprintf(out3,";\n");
			fprintf(faceout1,";\n");
		}
		imic=ns;
	}
	return(1);
}


int Gellipse(int jb1,int imic,double xcentreic,double ycentreic,double zcentreic,
			 double *angleic,double *xmic,double *zmic,double *aiic,double *biic)
{
//  fitting of generalized ellipe using x=xc+a0+(j=1 to J)(a2j-1*cos(j*angle)+a2j*sin(j*angle) and z=zc+b0+(j=1 to J)(b2j-1)*sin(j*angle)+b2j*cos(j*angle)
//  The curves can be well approximated with the constant term a0.Therefore,the following algorithm is very good.

	int i,j,k,n;
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
	    ciggj(aax,n,bbx);
		for( k=0;k<=n-1;k++)aiic[k]=bbx[k];
	

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
			for(i=0; i <= imic-1; i++)
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
		for( k=0;k<=n-1;k++)biic[k]=bbx[k];
//		printf("in imic= %d\n",imic);
		int ic=0;
		double x1,y1,z1;
	if(imic==50)
	{
// 		fprintf(out2,"curve -d 3\n");
		for(i=0; i <= imic; i++)
		{
			x1=xcentreic+aiic[0];
			y1=ycentreic;
			z1=zcentreic+biic[0];
			for(int j=1;j<=jb;j++)
			{        

                x1=x1+aiic[2*j-1]*cos(j*angleic[i])
					+aiic[2*j]*sin(j*angleic[i]);
				y1=y1+0.0;
				z1=z1+biic[2*j-1]*sin(j*angleic[i])
					+biic[2*j]*cos(j*angleic[i]);
			}
//		    fprintf(out,"i,angle,x1,y1,z1= %d %le %le %le %le\n",i,angleic[i],x1,y1,z1);
//		fprintf(out2,"-p   %le %le %le\n",x1,y1,z1);
		}
//	fprintf(out2,";\n");
	}

		return(1);

}

void CentreAngle(int imic,double &xcentreic,double &ycentreic,double &zcentreic,
				 double *angleic,double *xmic,double *ymic, double *zmic,double *aiic,
				 double *biic)
{
	int j;
	double anglem;
	xcentreic=0.0;
	ycentreic=0.0;
	zcentreic=0.0;
    for(j=0; j <= imic-1; j++)
	{
		xcentreic=xcentreic+xmic[j];
        ycentreic=ycentreic+ymic[j];
        zcentreic=zcentreic+zmic[j];
	}
    xcentreic=xcentreic/imic;
    ycentreic=ycentreic/imic;
    zcentreic=zcentreic/imic;
	
//	printf("in  x-y-zcentre[i]= %le %le %le\n",xcentreic,ycentreic,zcentreic);

// 	fprintf(out2,"curve -d 3\n");
    for(j=0; j <= 1; j++)
	{
//		fprintf(out2,"-p   %le %le %le\n",xcentreic,ycentreic,zcentreic);
	}
//	fprintf(out2,";\n");

	for(j=0; j <= imic-1; j++)
	{
		angleic[j]=0.0;
		anglem=(zmic[j]-zcentreic)/sqrt((xmic[j]-xcentreic)*(xmic[j]-xcentreic)
			+(ymic[j]-ycentreic)*(ymic[j]-ycentreic)+(zmic[j]-zcentreic)*(zmic[j]
			-zcentreic));
		anglem=asin(anglem);
		if((zmic[j]-zcentreic)>=0.0 && (xmic[j]-xcentreic)>=0.0)angleic[j]=anglem; 
		if((zmic[j]-zcentreic)>=0.0 && (xmic[j]-xcentreic)< 0.0)angleic[j]=pi-anglem; 
		if((zmic[j]-zcentreic)< 0.0 && (xmic[j]-xcentreic)< 0.0)angleic[j]=pi-anglem; 
		if((zmic[j]-zcentreic)< 0.0 && (xmic[j]-xcentreic)> 0.0)angleic[j]=2.0*pi+anglem; 
	}
}

void u0if(int &i, double &u0,double &u0i)
  { 
    int j;
    u0i = 1.0;
    if(i < 1)
    {
      u0i = 1.0;
      goto L1;
    }
    for (j=1; j <= i; j++)
    {
      u0i = u0i * u0;
    }
  L1: printf("");
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


*/