#include <math.h> 

/*--------------------------------------------------------------*/
double evalfpoly(int order, double x, double *c)
/*--------------------------------------------------------------*/
{
  int i,j;
  double y=c[0];
  x=(x-(0.000000000000000))/(5.729577951308232);
  int cond = order+1;
  for(i=1,j=1;j<cond;j+=2,i++)
    y += c[j]*cos(i*x)+c[j+1]*sin(i*x);
  return y;
}
 
/*--------------------------------------------------------------*/
double SFR_10Gyr(double x)
/*--------------------------------------------------------------*
   TableCurve Function: F:\\BACKUP\\Documents\\!SCHOOL\\Masters\\SFH- inside out formation\\10Gyr.c Dec 7, 2010 9:58:52 PM 
   C:\\Program Files\\TableCurve2Dv5.01\\CLIPBRD.PRN 
   X= Radius 
   Y= SFR 
   Eqn# 6844  Fourier Series Polynomial 4x2 
   r2=0.9997121162120392 
   r2adj=0.9974090459083524 
   StdErr=0.3788328357774642 
   Fstat=434.0779847717518 
   a= -19.61853168008389 
   b= -2.104522027080407 
   c= 43.32387183696036 
   d= 21.62121906454373 
   e= 12.73212733295415 
   f= 4.002422643119834 
   g= -2.622259504214302 
   h= 0.1125790218136833 
   i= 0.8656379107122378 
 *--------------------------------------------------------------*/
{
  double y;
  static double c[]= { 
    -19.61853168008389, 
    -2.104522027080407, 
    43.32387183696036, 
    21.62121906454373, 
    12.73212733295415, 
    4.002422643119834, 
    -2.622259504214302, 
    0.1125790218136833, 
    0.8656379107122378, 
    }; 
  y=evalfpoly(8,x,c); 
  return(y);
}
 
