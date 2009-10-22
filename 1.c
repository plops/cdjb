#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <fftc8.h>
#include <fftfreq.h>

#include <malloc.h>

typedef complex double c;


// pointers to all functions, index with p where number of points is n=2^p
void (*fft_a[])(complex8*in)={
  fftc8_2,
  fftc8_4,
  fftc8_8,
  fftc8_16,
  fftc8_32,
  fftc8_64,
  fftc8_128,
  fftc8_256,
  fftc8_512,
  fftc8_1024,
  fftc8_2048,
  fftc8_4096,
  fftc8_8192};

void (*fft_scale_a[])(complex8*in)={
  fftc8_scale2,
  fftc8_scale4,
  fftc8_scale8,
  fftc8_scale16,
  fftc8_scale32,
  fftc8_scale64,
  fftc8_scale128,
  fftc8_scale256,
  fftc8_scale512,
  fftc8_scale1024,
  fftc8_scale2048,
  fftc8_scale4096,
  fftc8_scale8192};


void
fft(int p,c*in)
{
  fft_a[p]((complex8*)in);
  fft_scale_a[p]((complex8*)in);  
}

unsigned int*ind=0;
int ind_len=0;
void
shuffle_init(int p)
{
  int n=1<<p;
  if(ind==0||n!=ind_len){
    free(ind);
    ind=(unsigned int*)malloc(n*sizeof(ind));
    ind_len=n;
  }
  fftfreq_ctable(ind,n);
}

void
shuffle(int p,c*in,c*out)
{
  int i,n=1<<p;
  if(ind_len!=n)
    shuffle_init(p);
  for(i=0;i<n;i++)
    out[ind[i]]=in[i];
}

// for saving first fouriertransforms
c*out_row=0;
int out_row_w=0,out_row_h=0;

c*in_col=0;
int in_col_len=0;

c*out_col=0;
int out_col_len=0;

// transform image with cols=2^p columns and rows=2^q rows
void
fft2(int p,int q,c*in,c*out)
{
  int j,cols=1<<p,rows=1<<q;
  if(!out_row || out_row_w!=cols || out_row_h!=rows){
    free(out_row);
    out_row=(c*)malloc(sizeof(c)*rows*cols);
    out_row_w=cols;
    out_row_h=rows;
  }
  for(j=0;j<rows;j++){// rows  
    fft(p,in+cols*j);   // transform j-th row 
    shuffle(p,in+cols*j,out_row+cols*j); // unshuffle result into out_row
  }

  if(!in_col || in_col_len!=rows){
    free(in_col);
    in_col=(c*)malloc(sizeof(c)*rows);
    in_col_len=rows;
  }

  if(!out_col || out_col_len!=rows){
    free(out_col);
    out_col=(c*)malloc(sizeof(c)*rows);
    out_col_len=rows;
  }

  int i;
  for(i=0;i<cols;i++){  // do ft column wise
    for(j=0;j<rows;j++)
      in_col[j]=out_row[i+cols*j];
    fft(q,in_col);
    shuffle(q,in_col,out_col);
    for(j=0;j<rows;j++)
      out[i+cols*j]=out_col[j];
  }
}

double
maximum(int n,double*a)
{
  double m=a[0];
  int i;
  for(i=1;i<n;i++)
    if(a[i]>m)
      m=a[i];
  return m;
}

double *oabs=0;
unsigned char *oabsc=0;
int oabs_w=0,oabs_h=0;

void
writepgm(char*name,c*out,int w,int h,int type)
{
  if(!oabs||!oabsc||oabs_w!=w||oabs_h!=h){
    free(oabs);
    free(oabsc);
    oabs=(double*)malloc(sizeof(double)*w*h);
    oabsc=(unsigned char*)malloc(w*h);
    oabs_w=w;
    oabs_h=h;
  }

  int i;
  switch(type){
  case 0: // absolute value
    for(i=0;i<w*h;i++)
      oabs[i]=cabs(out[i]);
    double m=maximum(w*h,oabs);
    for(i=0;i<w*h;i++)
      oabsc[i]=(unsigned char)(oabs[i]*255./m);
    break;
  case 1: // log of absolute value
    for(i=0;i<w*h;i++){
      double v=cabs(out[i]);
      oabs[i]=v?log(v):0;
    }
    double s=255./maximum(w*h,oabs);
    for(i=0;i<w*h;i++)
      oabsc[i]=(unsigned char)(oabs[i]*s);
    break;
  case 2:{ // argument (-pi .. pi) scaled to 0..255
    double s=.5/M_PI*255.;
    for(i=0;i<w*h;i++)
      oabsc[i]=(unsigned char)((carg(out[i])+M_PI)*s);
    break;}
  }
  
  FILE*f=fopen(name,"w");
  fprintf(f,"P5\n%d %d\n255\n",w,h);
  fwrite(oabsc,w,h,f);
  fclose(f);
}


enum {P=8,Q=9,W=1<<P,H=1<<Q};
c in[W*H],out[W*H];

int
main()
{
  in[H/2+W*(H/2+3)]=1;
  in[H/2+W*(H/2+12)]=1;
  
  fft2(P,Q,in,out);
  
  writepgm("/dev/shm/o.pgm",out,W,H,0);
  
  return 0;
}
