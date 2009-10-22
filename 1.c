#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <fftc8.h>
#include <fftfreq.h>

enum { P=10,
       N=1024 // =2^P, size of fourier transform
 };

typedef complex double c;

double
chop(double a)
{
  if(fabsl(a)<1e-12)
    return 0;
  return a;
}

void
fft(c*in)
{
  fftc8_1024((complex8*)in);
  fftc8_scale1024((complex8*)in);  
}

unsigned int ind[N];
void
shuffle_init()
{
  fftfreq_ctable(ind,N);
}

void
shuffle(c*in,c*out)
{
  int i;
  for(i=0;i<N;i++)
    out[ind[i]]=in[i];
}

c in[N*N],out[N*N],out_row[N*N],in_col[N],out_col[N];

void
fft2(c*in,c*out)
{
  int j;
  for(j=0;j<N;j++){  
    fft(in+N*j);   // transform j-th row 
    shuffle(in+N*j,out_row+N*j); // unshuffle result into out_row
  }
  int i;
  for(i=0;i<N;i++){  // do ft column wise
    for(j=0;j<N;j++)
      in_col[j]=out_row[i+N*j];
    fft(in_col);
    shuffle(in_col,out_col);
    for(j=0;j<N;j++)
      out[i+N*j]=out_col[j];
  }
}

double oabs[N*N];
unsigned char oabsc[N*N];

double
maxim(double*a,int n)
{
  double m=a[0];
  int i;
  for(i=1;i<N;i++)
    if(a[i]>m)
      m=a[i];
  return m;
}

void
writepgm(char*name,unsigned char*buf,int w,int h)
{
  FILE*f=fopen(name,"w");
  fprintf(f,"P5\n%d %d\n255\n",w,h);
  fwrite(buf,w,h,f);
  fclose(f);
}

int
main()
{
  shuffle_init();

  in[1024/2+1024*512]=1;
  in[1024/2+1024*514]=1;
  
  fft2(in,out);
  /*
  int i;
  for(i=0;i<N*N;i++)
    oabs[i]=cabs(out[i]);
  double m=maxim(oabs,N*N);
  for(i=0;i<N*N;i++)
    oabsc[i]=(unsigned char)(oabs[i]*255./m);
  
  writepgm("/dev/shm/o.pgm",oabsc,N,N);
  */
  return 0;
}
