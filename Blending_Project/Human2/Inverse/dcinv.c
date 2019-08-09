#include <stdlib.h>
#include <math.h>
#include <stdio.h>

     
int  dcinv();
void damul();

main()
{
int i,j,k;
double a[4][4],b[4][4],c[4][4],ff;

a[0][0]=0.2368;
a[0][1]=0.2471;
a[0][2]=0.2568;
a[0][3]=1.2671;
a[1][0]=1.1161;
a[1][1]=0.1254;
a[1][2]=0.1397;
a[1][3]=0.149;
a[2][0]=0.1582;
a[2][1]=1.1675;
a[2][2]=0.1768;
a[2][3]=0.1871;
a[3][0]=0.1968;
a[3][1]=0.2071;
a[3][2]=1.2168;
a[3][3]=0.2271;
for(i=0;i<=3;i++)
for(j=0;j<=3;j++)b[i][j]=a[i][j];
dcinv(a,4);
for(i=0;i<=3;i++)
for(j=0;j<=3;j++)
printf("i,j,a[i][j]= %d %d %le\n",i,j,a[i][j]);
damul(b,a,4,4,4,c);
for(i=0;i<=3;i++)
for(j=0;j<=3;j++)
printf("i,j,c[i][j]= %d %d %le\n",i,j,c[i][j]);
}

int dcinv(a,n)
int n;
double a[];
{
  int *is,*js,i,j,k,l,u,v;
  double d,p;
  is=malloc(n* sizeof(int));
  js=malloc(n* sizeof(int));
  for(k=0;k <= n-1;k++)
  {
    d=0.0;
    for(i=k;i <= n-1; i++)
    for(j=k;j <= n-1; j++)
    {
      l=i*n+j;
      p=fabs(a[l]);
      if(p>d)
      {  
        d=p;
        is[k]=i;
        js[k]=j;
      }
    }
    if(d+1.0 == 1.0)
    {
      free(is);
      free(js);
      printf("err * * not inv\n");
      return(0);
    }
    if(is[k] != k)
    for (j=0;j <= n-1; j++)
    {
      u=k*n+j;
      v=is[k]*n+j;
      p=a[u];
      a[u]=a[v];
      a[v]=p;
    }
    if(js[k] != k)
    for (i=0; i <= n-1; i++)
    {
      u=i*n+k;
      v=i*n+js[k];
      p=a[u];
      a[u]=a[v];
      a[v]=p;
    }
    l=k*n+k;
    a[l]=1.0/a[l];
    for (j=0; j <= n-1;j++)
    if(j != k)
    {
      u=k*n+j;
      a[u]=a[u]*a[l];
    }
    for (i=0; i <= n-1; i++)
    if( i != k)
    for (j=0; j <= n-1; j++)
    if ( j != k)
    {
      u=i*n+j;
      a[u]=a[u]-a[i*n+k]*a[k*n+j];
    }
    for (i=0; i <= n-1; i++)
    if(i != k)
    {
      u=i*n+k;
      a[u]=-a[u]*a[l];
    }
    }
    for (k=n-1; k >= 0; k--)
    {
      if(js[k] != k)
      for (j=0; j <= n-1; j++)
      {
        u=k*n+j;
        v=js[k]*n+j;
        p=a[u];
        a[u]=a[v];
        a[v]=p;
      }
      if(is[k] != k)
      for (i=0; i <= n-1; i++)
      {
        u=i*n+k;
        v=i*n+is[k];
        p=a[u];
        a[u]=a[v];
        a[v]=p; 
      }
    }
    free(is);
    free(js);
    return(1);
}

void damul(a,b,m,n,k,c)
int m,n,k;
double a[],b[],c[];
{
  int i,j,l,u;
  for (i=0;i <= m-1; i++)
  for (j=0;j <= k-1; j++)
  {
    u=i*k+j;
    c[u]=0.0;
    for (l=0;l <= n-1; l++)
    c[u]=c[u]+a[i*n+l]*b[l*k+j];
  }
  return;
}
