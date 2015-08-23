#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define MAT1 3
#define TINY 1e-20
#define a(i,j) a[(i)*MAT1+(j)]

void h_pivot_decomp(float *a, int *p, int *q){
    int i,j,k;
    int n=MAT1;
    int pi,pj,tmp;
    float max;
    float ftmp;
    for (k=0;k<n;k++){
        pi=-1,pj=-1,max=0.0;
        //find pivot in submatrix a(k:n,k:n)
        for (i=k;i<n;i++) {
            for (j=k;j<n;j++) {
                if (fabs(a(i,j))>max){
                    max = fabs(a(i,j));
                    pi=i;
                    pj=j;
                }
            }
        }
        //Swap Row
        tmp=p[k];
        p[k]=p[pi];
        p[pi]=tmp;
        for (j=0;j<n;j++){
            ftmp=a(k,j);
            a(k,j)=a(pi,j);
            a(pi,j)=ftmp;
        }
        //Swap Col
        tmp=q[k];
        q[k]=q[pj];
        q[pj]=tmp;
        for (i=0;i<n;i++){
            ftmp=a(i,k);
            a(i,k)=a(i,pj);
            a(i,pj)=ftmp;
        }
        //END PIVOT

        //check pivot size and decompose
        if ((fabs(a(k,k))>TINY)){
            for (i=k+1;i<n;i++){
                //Column normalisation
                ftmp=a(i,k)/=a(k,k);
                for (j=k+1;j<n;j++){
                    //a(ik)*a(kj) subtracted from lower right submatrix elements
                    a(i,j)-=(ftmp*a(k,j));
                }
            }
        }
        //END DECOMPOSE
    }
}


void h_solve(float *a, float *x, int *p, int *q){
    //forward substitution; see  Golub, Van Loan 96
    //And see http://www.cs.rutgers.edu/~richter/cs510/completePivoting.pdf
    int i,ii=0,ip,j,tmp;
    float ftmp;
    float xtmp[MAT1];
    //Swap rows (x=Px)
    puts("x=Px Stage");
    for (i=0; i<MAT1; i++){
        xtmp[i]=x[p[i]]; //value that should be here
        printf("x:%.1lf,q:%d\n",xtmp[i],q[i]);
    }
    //Lx=x
    puts("Lx=x Stage");
    for (i=0;i<MAT1;i++){
        ftmp=xtmp[i];
        if (ii != 0)
            for (j=ii-1;j<i;j++)
                ftmp-=a(i,j)*xtmp[j];
        else
            if (ftmp!=0.0)
                ii=i+1;
        xtmp[i]=ftmp;
        printf("x:%.1lf,q:%d\n",xtmp[i],q[i]);
    }
    puts("Ux=x");
    //backward substitution
    //partially taken from Sourcebook on Parallel Computing p577
    //solves Uy=z
    xtmp[MAT1-1]/=a(MAT1-1,MAT1-1);
    for (i=MAT1-2;i>=0;i--){
        ftmp=xtmp[i];
        for (j=i+1;j<MAT1;j++){
            ftmp-=a(i,j)*xtmp[j];
        }
        xtmp[i]=(ftmp)/a(i,i);
    }
    for (i=0;i<MAT1;i++)
        printf("x:%.1lf,q:%d\n",xtmp[i],q[i]);

    //Last bit
    //solves x=Qy
    puts("x=Qx Stage");
    for (i=0;i<MAT1;i++){
        x[i]=xtmp[q[i]];
        printf("x:%.1lf,q:%d\n",x[i],q[i]);
    }
}


void main(){
    //3x3 Matrix
    //float a[]={1,-2,3,2,-5,12,0,2,-10};
    float a[]={1,3,-2,3,5,6,2,4,3};
    float b[]={5,7,8};
    //float a[]={1,2,3,2,-1,1,3,4,-1};
    //float b[]={14,3,8};
    //float a[]={1,-2,1,0,2,2,-2,4,2};
    //float b[]={1,4,2};
    int sig;
    puts("Declared Stuff");

    //pivot array (not used currently)
    int* p_pivot = (int *)malloc(sizeof(int)*MAT1);
    int* q_pivot = (int *)malloc(sizeof(int)*MAT1);
    puts("Starting Stuff");
    for (unsigned int i=0; i<MAT1; i++){
        p_pivot[i]=i;
        q_pivot[i]=i;
        printf("%.1lf|",b[i]);
        for (unsigned int j=0;j<MAT1; j++){
            printf("%.1lf,",a(i,j));
        }
        printf("|%d,%d",p_pivot[i],q_pivot[i]);
        puts("");
    }

    h_pivot_decomp(&a[0],p_pivot,q_pivot);
    puts("After Pivot");
    for (unsigned int i=0; i<MAT1; i++){
        printf("%.1lf|",b[i]);
        for (unsigned int j=0;j<MAT1; j++){
            printf("%.1lf,",a(i,j));
        }
        printf("|%d,%d",p_pivot[i],q_pivot[i]);
        puts("");
    }

    h_solve(&a[0],&b[0],p_pivot,q_pivot);
    puts("Finished Solve");

    for (unsigned int i=0; i<MAT1; i++){
        printf("%.1lf|",b[i]);
        for (unsigned int j=0;j<MAT1; j++){
            printf("%.1lf,",a(i,j));
        }
        puts("");
    }
}