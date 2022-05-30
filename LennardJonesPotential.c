#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

/*
double x[16] = { 1.09,3.12,0.08,0.54,2.52,3.03,4.25,0.89,2.76,3.14,0.23,1.91,4.77,5.10,4.97,3.90 };
double y[16]={0.98,5.25,2.38,4.08,4.39,2.94,3.01,3.11,0.31,1.91,5.71,2.46,0.96,4.63,5.88,0.20};
double vx[16]={-0.33,0.12,-0.08,-1.94,0.75,1.70,0.84,-1.04,1.64,0.38,-1.58,-1.55,-0.23,-0.31,1.18,0.46};
double vy[16]={0.78,-1.19,-0.10,-0.56,0.34,-1.08,0.47,0.06,1.36,-1.24,0.55,-0.16,-0.83,0.65,1.48,-0.51};
double ax[16],ay[16];
*/
double *x, *y, *vx, *vy, *ax, *ay;

FILE *md3;
FILE *data;
FILE *xball3;

int nflag=0;
int nprint,ff=1,datatype = 0;
double N,Lx,Ly,dt,dt2,flag;
double ds,L,ke,pe,kecum,pecum,vcum,area,t;
double dx,dy,fx,fy,pot,virial,rcm;

void ini(double *t,double *ke,double *area);
void red_cuadrada();
void red_aleatoria();
void red_triangular();
void force(double dx, double dy, double *fx, double *fy, double *pot);
void accel(double *pe, double *virial);
double fseparation(double ds,double L);
double fpbc(double pos,double L);
void verlet(double *t, double *ke, double *pe, double *virial);



//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void main()
{ 
  rcm = 0.5*pow(2,1./6);
  md3 = fopen("md3.dat","w");
  xball3 = fopen("xball3.dat","w");

  srand(time(0));
  ini(&t,&ke,&area);

  fprintf(xball3,"%lf %lf %lf %lf %lf %lf\n",N,t,0.,0.,Lx,Ly);
    for(int i = 0;i<N;i++)
    {
      fprintf(xball3,"%lf %lf %lf %lf %lf\n",x[i],y[i],vx[i],vy[i],rcm);
    }
  accel(&pe,&virial);

  double E = ke + pe;
  double ncum = 0;
 
do
  {
    verlet(&t,&ke,&pe,&virial);

    if(ff==0)
      {
	if(nflag%nprint==0)
	  {
	    fprintf(xball3,"%lf %lf %lf %lf %lf %lf\n",N,t,0.,0.,Lx,Ly);
	    for(int i = 0;i<N;i++)
	      {
		fprintf(xball3,"%lf %lf %lf %lf %lf\n",x[i],y[i],vx[i],vy[i],rcm);
	      }
	  }
      }
    else
      {
	fprintf(xball3,"%lf %lf %lf %lf %lf %lf\n",N,t,0.,0.,Lx,Ly);
	for(int i = 0;i<N;i++)
	  {
	    fprintf(xball3,"%lf %lf %lf %lf %lf\n",x[i],y[i],vx[i],vy[i],rcm);
	  }
      }
    E = ke + pe;
    fprintf(md3,"%lf %lf %lf %lf\n",t,ke,pe,E);
  }while(t < flag);
  
 fclose(md3);
 fclose(xball3); 
}
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void ini(double *t,double *ke,double *area)
{
//--------------------------------------------
  printf("Indique el numero de particulas. N = ");
  scanf("%lf",&N);
  printf("Indique el ancho de la caja Lx = ");
  scanf("%lf",&Lx);
  printf("Indique el largo de la caja Ly = ");
  scanf("%lf",&Ly);
  printf("Indique el paso del tiempo dt = ");
  scanf("%lf",&dt);
  printf("Indique las unidades de tiempo totales = ");
  scanf("%lf",&flag);
  printf("\n Escriba 0 si quiere correr el programa con datos aleatorios\n Escriba 1 si quiere leer del un archivo externo\n Escriba 2 si quiere generar una red cuadrada\n Escriba 3 si quiere genera una red triangular\n Escriba 4 para problema de gould\n  ---> ");
  scanf("%d",&datatype);
   printf("\nEscriba 0 si quiere imprimir intervalos de tiempo mas grandes para xball o cualquier otro numero entero para imprimir todos los pasos de tiempo\n --->");
   scanf("%d",&ff);
  
  if(ff==0)
    {
      printf("Indique el intervalo ENTERO entre los pasos de tiempo a imprimir\n --->");
      scanf("%d",&nprint);
    }
//--------------------------------------------
  x = (double*)calloc(N, sizeof(double));
  y = (double*)calloc(N, sizeof(double));
  vx = (double*)calloc(N, sizeof(double));
  vy = (double*)calloc(N, sizeof(double));
  ax = (double*)calloc(N, sizeof(double));
  ay = (double*)calloc(N, sizeof(double));
//--------------------------------------------
  if(datatype == 0)
    {
      red_aleatoria();
    }
  else if(datatype == 1)
    {
      data = fopen("data.dat","r");
      for(int i = 0;i < N;i++)
	{
	  fscanf(data,"%lf %lf %lf %lf\n",&x[i],&y[i],&vx[i],&vy[i]);
	  printf("%lf \t %lf \t %lf \t %lf\n",x[i],y[i],vx[i],vy[i]);
	}
      fclose(data);
    }
  else if(datatype == 2)
    {
      red_cuadrada();
    }
  else if(datatype == 3)
    {
      red_triangular();
    }
  else if(datatype == 4)
    {
      for(int i=0;i<N;i++)
	{
	  if(i==15)
	    {
	      vx[i]=0.99;
	      vy[i]=0.01; 
	    }
	  else
	    {
	      vx[i]=1;
	      vy[i]=0;
	    }
	  x[i]=Lx/2;
	  y[i]=(i-0.5)*Ly/N;
	 
	  printf("%lf \t %lf \t %lf \t %lf\n",x[i],y[i],vx[i],vy[i]);
	}
    }
  else
    {
      printf("No ingreso un valor valido para el programa. Adios, saboteador!\n");
      exit(0);
    }
//-------------------------------------------------
  dt2 = dt*dt;

  *ke = 0;
  for(int i=0;i < N;i++)
    {
      *ke = *ke + vx[i]*vx[i] + vy[i]*vy[i];
    }
  *ke = 0.5*(*ke);
  *area = Lx*Ly;
  *t = 0;
  kecum = 0;
  pecum = 0;
  vcum = 0;
}
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
double fseparation(double ds,double L) //criterio de imagen minima debido al potencial de corto alcance y la condicion de borde periodica.
{
  if(ds > 0.5*L) ds = ds - L;
  else if(ds < -0.5*L) ds = ds + L;
  else ds = ds;
  return ds;
}
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void force(double dx, double dy, double *fx, double *fy, double *pot)
{
  double r2 = dx*dx + dy*dy;
  double rm2 = 1/r2;
  double rm6 = rm2*rm2*rm2;
  double f_over_r = 24*rm6*(2*rm6-1)*rm2;
  *fx = f_over_r*dx;
  *fy = f_over_r*dy;
  *pot = 4*(rm6*rm6 - rm6);
}
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void accel(double *pe, double *virial) //cuidado con el pase de parametro por valor o por referencia.
{
  for(int k = 0;k<N;k++) //inicializando las matrices en 0.
    {
      ax[k]=0;
      ay[k]=0;
    }

  *virial = 0;
  *pe = 0;

  for(int i=0;i<N-1;i++)
    {
      for(int j=i+1;j<N;j++)
	{
	  double a = x[i]-x[j];
	  double b = y[i]-y[j];
	  dx = fseparation(a,Lx);
	  dy = fseparation(b,Ly);
	  //la aceleracion es la fuerza porque m = 1 en unid. reducidas
	  force(dx,dy,&fx,&fy,&pot);

	  ax[i] = ax[i] + fx;
	  ay[i] = ay[i] + fy;
          ax[j] = ax[j] - fx;
	  ay[j] = ay[j] - fy;
	  *pe = *pe + pot;
	  *virial = *virial +dx*ax[i] + dy*ay[i];
	}
    } 
}
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
double fpbc(double pos,double L) //condicion de borde periodica.
{
  double pbc;
  if(pos < 0) pbc = pos + L;
  else if(pos > L) pbc = pos - L;
  else pbc = pos;
  return pbc;
}
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void verlet(double *t, double *ke, double *pe, double *virial)
{
  for(int i = 0;i < N;i++)
    {
      double xnew = x[i] + vx[i]*dt + 0.5*ax[i]*dt2;
      double ynew = y[i] + vy[i]*dt + 0.5*ay[i]*dt2;

      vx[i] = vx[i] + 0.5*ax[i]*dt;
      vy[i] = vy[i] + 0.5*ay[i]*dt;

      x[i] = fpbc(xnew,Lx);
      y[i] = fpbc(ynew,Ly);
    }

  accel(pe,virial);

   *ke = 0;

  for(int i = 0; i < N; i++)
    {
      vx[i] = vx[i] + 0.5*ax[i]*dt;
      vy[i] = vy[i] + 0.5*ay[i]*dt;
      *ke = *ke + vx[i]*vx[i] + vy[i]*vy[i];
    }
  *ke = 0.5*(*ke);
  *t = *t + dt;
   nflag = nflag + 1;
}
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void red_aleatoria()
{
  double rc = pow(2,(1./6));
  x[0] = (rand()/(double)RAND_MAX)*(Lx+1);
  y[0] = (rand()/(double)RAND_MAX)*(Ly+1);
  vx[0]=0.;
  vy[0]=0.;
  
  for(int i = 1;i<=N;i++)
    {
      //--------POSICIONES Y VELOCIDADES----------//
      x[i] = (rand()/(double)RAND_MAX)*Lx+1;
      y[i] = (rand()/(double)RAND_MAX)*Ly+1;
      
      vx[i]=0.;
      vy[i]=0.;
      //------------------------------------------//
      for(int j = 0;j<i;j++)
	{
	  //printf("%lf \t %lf \t %lf \t %lf\n",x[i],y[i],vx[i],vy[i]);
	  double d1 = x[j]-x[i];
	  double d2 = y[j]-y[i];
	  double d12 = pow(d1*d1+d2*d2,(1./2));
	  
	  if(d12<rc)
	    {
	      x[i] = (rand()/(double)RAND_MAX)*Lx+1;
	      y[i] = (rand()/(double)RAND_MAX)*Ly+1;
	      j=0;
	    }
	}
    }
  for(int i = 0; i < N; i++)
    {
      printf("%lf \t %lf \t %lf \t %lf\n",x[i],y[i],vx[i],vy[i]);
    }
}
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void red_cuadrada()
{
  int rf;
  double cx = 0,cy = 0,dlx,dly;
  int n = floor(pow(N,1/2.));
  double ni = pow(N,1/2.);
  /*
  printf("Indique la separacion entre las particulas para la red cuadrada en X= ");
  scanf("%lf",&dlx);
  printf("Indique la separacion entre las particulas para la red cuadrada en Y= ");
  scanf("%lf",&dly);
  */

  printf("\nEscriba 0 si desea empezar la red en una ubicacion distinta al (0,0) o cualquier otro entero si no lo desea\n--->");
  scanf("%d",&rf);

  if(rf == 0)
  {
  printf("Indique la distancia de inicio en X = ");
  scanf("%lf",&cx);
  printf("Indique la distancia de inicio en Y = ");
  scanf("%lf",&cy);
  }
  
  printf("n = %d\n",n);

 if(ni-n>0)
    {
      for(int i=n;i<=(n+1);i++)
	{
	  for(int j=0;j<n;j++)
	    {
	      x[i*n+j]=cx+(2*j+1)*rcm;
	      y[i*n+j]=cy+(2*i+1)*rcm;
	    }
	}
    }  

  
  for(int i=0;i<n;i++)
    {
      for(int j=0;j<n;j++)
	{
	  //-------------------PARTICULAS JUNTAS-------------------------
	  x[i*n+j]=cx+(2*j+1)*rcm;
	  y[i*n+j]=cy+(2*i+1)*rcm;
	  //-------------------VARIANDO DISTANCIAS-----------------------
	  //x[i*n+j]= cx + rcm + 2*j*rcm + dlx; //REVISAR
	  //y[i*n+j]= cy + rcm + 2*i*rcm + dly; //REVISAR
	  //-------------------------------------------------------------
	}
    }

  for(int k = 0;k<N;k++)
    {
      vx[k]=0;
      vy[k]=0;
      printf("%lf \t %lf \t %lf \t %lf\n",x[k],y[k],vx[k],vy[k]);
    }
}
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void red_triangular()
{
   int rf;
  double cx = 0,cy = 0,dlx,dly;
  int n = floor(pow(N,1/2.));
  double ni = pow(N,1/2.);
  /*
  printf("Indique la separacion entre las particulas para la red cuadrada en X= ");
  scanf("%lf",&dlx);
  printf("Indique la separacion entre las particulas para la red cuadrada en Y= ");
  scanf("%lf",&dly);
  */

  printf("\nEscriba 0 si desea empezar la red en una ubicacion distinta al (0,0) o cualquier otro entero si no lo desea\n--->");
  scanf("%d",&rf);

  if(rf == 0)
  {
  printf("Indique la distancia de inicio en X = ");
  scanf("%lf",&cx);
  printf("Indique la distancia de inicio en Y = ");
  scanf("%lf",&cy);
  }
  
  printf("n = %d\n",n);

 if(ni-n>0)
    {
      for(int i=n;i<=(n+1);i++)
	{
	  for(int j=0;j<n;j++)
	    {
	      if(i % 2 == 0)
		{
		  x[i*n+j]=cx+(2*j+1)*rcm;
		  y[i*n+j]=cy+(2*i+1)*rcm;
		}
	      else
		{
		  x[i*n+j]=rcm+cx+(2*j+1)*rcm;
		  y[i*n+j]=cy+(2*i+1)*rcm;
		}
	    }
	}
    }  

  
  for(int i=0;i<n;i++)
    {
      for(int j=0;j<n;j++)
	{
	  if(i % 2 == 0)
	    {
	      //-------------------PARTICULAS JUNTAS-------------------------
	      x[i*n+j]=cx+(2*j+1)*rcm;
	      y[i*n+j]=cy+(2*i+1)*rcm;
	      //-------------------VARIANDO DISTANCIAS-----------------------
	      //x[i*n+j]= cx + rcm + 2*j*rcm + dlx; //REVISAR
	      //y[i*n+j]= cy + rcm + 2*i*rcm + dly; //REVISAR
	      //-------------------------------------------------------------
	    }
	  else
	    {
	      x[i*n+j]=rcm+cx+(2*j+1)*rcm;
	      y[i*n+j]=cy+(2*i+1)*rcm;
	    }
	}
    }

  for(int k = 0;k<N;k++)
    {
      vx[k]=0;
      vy[k]=0;
      printf("%lf \t %lf \t %lf \t %lf\n",x[k],y[k],vx[k],vy[k]);
    }
}
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
