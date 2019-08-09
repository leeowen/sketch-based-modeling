/*
cc -o fdm_surface11 fdm_surface11.c -lglut -lGLU -lGL -lXmu -lXext -lX11 -lm
*/

#include <GL/glut.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>

FILE *mel1;

void main()
{
	int i;
	mel1=fopen("mel1.mel","wt");
	fprintf(mel1,"select -r polySurface1.vtx[0] \n");
//	for(i=1; i<=607; i++)fprintf(mel1,"polySurface1.f[" d% "] \n",i);
	for(i=1; i<=607; i++)fprintf(mel1,"polySurface1.vtx[%d] \n",i);
	fprintf(mel1,";\n");

}

 