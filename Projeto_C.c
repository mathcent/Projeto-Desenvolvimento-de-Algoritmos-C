#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>

void inicializaMatriz(int N, double ***M);
void destroiMatriz(int N, double ***M);
void imprimeMatriz(int N, double **M);
void lerMatriz(int N, double ***M);

void inicializaMatriz(int N, double ***M){
    *M = (double**)malloc(sizeof(double*)*N);
    for(int i=0;i<N;i++){
        (*M)[i] = (double*)malloc(sizeof(double)*(N+1));
    }
    for(int i=0;i<N;i++){
        for(int j=0;j<=N;j++)
            (*M)[i][j] = 1.0;
    }

}

void destroiMatriz(int N, double ***M){
    for(int i=0;i<N;i++)
       free((*M)[i]);
    
    free(*M);
}

void imprimeMatriz(int N, double **M){
    for(int i=0;i<N;i++){
        for(int j=0;j<=N;j++)
            printf("%.1lf ", M[i][j]);
        printf("\n");
    }
    printf("\n");
}

void lerMatriz(int N, double ***M){
    double x;
    for(int i=0;i<N;i++){
        for(int j=0;j<=N;j++){
            scanf("%lf", &x);
            (*M)[i][j] = x;
        }
    }
}

void trocaLinhas(int N, double ***M, int row1, int row2){
    for(int i=0;i<=N;i++){
        double aux = (*M)[row1][i];
        (*M)[row1][i] = (*M)[row2][i];
        (*M)[row2][i] = aux; 
    }
}

void subtrLinha(int N, double ***M, int row1, int row2, double k){
    for(int i=0; i<=N;i++)
            (*M)[row1][i] -= (*M)[row2][i]*k;
}

int Gauss(int N, double ***M){
    for(int j=0;j<N-1;j++){
        double max = (*M)[j][j];
        double piv = j;
        for(int k=j+1;k<N;k++){
            if(abs((*M)[k][j]) > max){
                max = (*M)[k][j];
                piv = k;
            }
        }
        trocaLinhas(N, M, j, piv);
        
        for(int i=j+1;i<N;i++){
            double mi = (*M)[i][j]/(*M)[j][j];
            subtrLinha(N, M, i, j, mi);
        }
    }
    for(int i=0;i<N;i++){
        if((*M)[i][i] == 0)
            return 0;
    }
    return 1;
}

double* solveLin(int N, double **M){
    double *x = malloc(sizeof(double)*N);
    x[N-1] = M[N-1][N]/M[N-1][N-1];
    for(int i=N-2;i>=0;i--){
        double soma=0;
        for(int j=i+1;j<N;j++)
            soma += x[j]*M[i][j];
        x[i] = (M[i][N] - soma)/M[i][i];
    }
    return x;
}

int main(int argc, char* argv[]){
    int N;
    scanf("%d", &N);

    double **M;
    inicializaMatriz(N, &M);
    lerMatriz(N, &M);
    imprimeMatriz(N, M);
    if(Gauss(N, &M)){
        imprimeMatriz(N, M);
        double *result = solveLin(N, M);
        for(int i=0;i<N;i++)
            printf("%.1lf\n", result[i]);
        destroiMatriz(N, &M);
        free(result);   
    }
    else
        printf("Sistema sem solucao.\n");
    
    return 0;
}
