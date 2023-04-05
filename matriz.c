//
// Created by David on 16/09/2021.
//

#include "matriz.h"
/// Reserva memoria dinamicamente para una estructura matrix (definida en alglin/matriz.h)
/// \param matriz Puntero a la estructura matrix
/// \param m Numero de filas de la matriz
/// \param n Numero de columnas de la matriz
void matriz_inicializar(matrix *matriz, int m, int n) {

    double** puntero=NULL;
    int i;
    matriz->m = m, matriz->n =n;

    //Se asigna memoria dinamicamente a un array de punteros (double*), que representan las filas de la matriz
    puntero = (double**) malloc(m*sizeof (double*));
    if(puntero==NULL){
        fprintf(stderr,ERR_MEM);
        exit(-1);
    }
    //Se inicializa la matriz de punteros con el valor 0
    memset(puntero,0,m* sizeof(double*));
    matriz->datos = puntero;
    //Asignacion dinamica de memoria a cada fila
    for (i=0;i<m;i++){
        puntero[i]= (double*) malloc(n* sizeof(double));
        //Si la asignacion dinamica fracasa, liberar la memoria asignada a la matriz y salir del programa
        if (puntero[i]==NULL){
            matriz_liberar(*matriz);
            matriz->datos=NULL;
            fprintf(stderr,ERR_MEM);
            exit (-1);
        }
        //Finalmente, se inicializa la fila con ceros
        memset(puntero[i],0,n* sizeof(double ));
    }
}

void matriz_nula(matrix matriz) {
    int i, j;
    for (i = 0; i < matriz.m; i++) {
        for (j = 0; j < matriz.n; j++) {
            matriz.datos[i][j] = 0;
        }
    }
}

void matriz_imprimir(matrix matriz) {
    int i,j,k,l;
    int fils = ceil((double)matriz.n/MAX_COLS_IMPRIMIR);
    int j_max;
    int j_min;


    for (k=0;k<fils;k++){
        if (k==fils-1){
            j_max=matriz.n;
        }
        else{
            j_max = k*MAX_COLS_IMPRIMIR+MAX_COLS_IMPRIMIR;
        }
        j_min=MAX_COLS_IMPRIMIR*k;
        printf("%9s","");
        fflush(stdout);
        for (l=j_min;l<j_max;l++){
            printf("%.4d",l+1);
            printf("%13s","");
            fflush(stdout);
        }
        printf("\n");
        for (i=0;i<matriz.m;i++){
            printf("%.4d    ",i+1);
            for (j=k*MAX_COLS_IMPRIMIR;j<j_max;j++){
                if(matriz.datos[i][j] >= 0){
                    printf(" ");
                }
                if (fabsl(matriz.datos[i][j]) < TOLERANCIA_GRAL){
                    printf(ANSI_COLOR_CYAN"%.10lf     "ANSI_COLOR_RESET,matriz.datos[i][j]);

                }
                else{
                    printf("%.10lf     ",matriz.datos[i][j]);
                }

            }
            printf("\n");
        }
        printf("\n");
    }
    fflush(stdout);
}

void matriz_sumar(matrix *matrizA, matrix *matrizB, matrix *matrizS) {

    int i, j;

    if (matrizA->m != matrizB->m || matrizA->n != matrizB->n) {
        fprintf(stderr, ERR_ALG);
        exit(1);
    }

    if (matrizA->m != matrizS->m || matrizA->n != matrizS->n) {
        fprintf(stderr, ERR_ALG);
        exit(1);
    }

    for (i = 0; i < matrizA->m; i++) {
        for (j = 0; j < matrizA->n; j++) {
            matrizS->datos[i][j] = matrizA->datos[i][j] + matrizB->datos[i][j];
        }
    }

}
void matriz_restar(matrix *matrizA, matrix *matrizB, matrix *matrizS) {

    int i, j;

    if (matrizA->m != matrizB->m || matrizA->n != matrizB->n) {
        fprintf(stderr, ERR_ALG);
        exit(1);
    }

    if (matrizA->m != matrizS->m || matrizA->n != matrizS->n) {
        fprintf(stderr, ERR_ALG);
        exit(1);
    }

    for (i = 0; i < matrizA->m; i++) {
        for (j = 0; j < matrizA->n; j++) {
            matrizS->datos[i][j] = matrizA->datos[i][j] - matrizB->datos[i][j];
        }
    }

}

void matriz_copiar(matrix matrizA, matrix matrizC) {
    int i; //Número de filas de la matriz a copiar
    if (matrizA.m != matrizC.m || matrizA.n != matrizC.n){
        printf ("ERROR:LA MATRIZ DE ORIGEN Y LA DE DESTINO NO TIENEN LA MISMA DIMENSION\n");
        printf("Dimension A: %dx%d\n",matrizA.m,matrizA.n);
        printf("Dimension C: %dx%d\n",matrizC.m,matrizC.n);
        exit (-1);
    }
    for (i = 0; i < matrizA.m; i++) {
        memcpy(matrizC.datos[i], matrizA.datos[i], matrizA.n * sizeof(double));
    }
}

void matriz_trasponer(matrix matriz, matrix matriz_traspuesta) {
    int filas_M  =matriz.m;
    int columnas_M = matriz.n;
    int filas_Mt = matriz_traspuesta.m;
    int columnas_Mt = matriz_traspuesta.n;
    printf ("fM %d  cM %d   fMt %d  cMt %d\n",filas_M,columnas_M,filas_Mt,columnas_Mt);
    int i,j;
    if ((filas_M!=columnas_Mt)||(columnas_M!=filas_Mt)){
        fprintf(stderr,ERR_ALG);
        exit (1);
    }
    for (i=0;i<filas_Mt;i++){
        for (j=0;j<columnas_Mt;j++){
            matriz_traspuesta.datos[i][j]=matriz.datos[j][i];
        }
    }
}


void matriz_liberar(matrix matriz) {
    int i;

    for (i = 0; i < matriz.m; i++) {
        free(matriz.datos[i]);
    }
    free(matriz.datos);
}



void matriz_identidad(matrix matriz) {
    int i, j;
    if (matriz.m != matriz.n) {
        fprintf(stderr, ERR_ALG);
        exit(1);
    }
    for (i = 0; i < matriz.m; i++) {
        for (j = 0; j <= i; j++) {
            if (i != j) {
                matriz.datos[i][j] = 0;
                matriz.datos[j][i] = 0;
            } else {
                matriz.datos[i][j] = 1;
            }
        }
    }
}


void matriz_producto(matrix matrizA, matrix matrizB, matrix matrizC) {
    int m_matrizA = matrizA.m;
    int n_matrizA = matrizA.n;
    int m_matrizB = matrizB.m;
    int n_matrizB = matrizB.n;
    int i, j, k;
    matrix aux;

    if (n_matrizA!=m_matrizB){
        fprintf(stderr,ERR_ALG);
        exit(-1);
    }

    if(matrizA.datos==matrizC.datos || matrizB.datos==matrizC.datos)
        matriz_inicializar(&aux,m_matrizA,n_matrizB);
    else
        matriz_nula(matrizC);

    for (i = 0; i < m_matrizA; i++) {
        for (j = 0; j < n_matrizB; j++) {
            for (k = 0; k < n_matrizA; k++) {
                if(matrizA.datos==matrizC.datos || matrizB.datos==matrizC.datos)
                    aux.datos[i][j] += matrizA.datos[i][k] * matrizB.datos[k][j];
                else
                    matrizC.datos[i][j] += matrizA.datos[i][k] * matrizB.datos[k][j];
            }
        }
    }

    if (matrizA.datos==matrizC.datos || matrizB.datos==matrizC.datos){
        matriz_copiar(aux,matrizC);
        matriz_liberar(aux);
    }
}

double matriz_traza (matrix matrizA){
    int i;
    int N = matrizA.m;
    double traza = 0;
    for (i=0;i<N;i++){
        traza += matrizA.datos[i][i];
    }
    return traza;
}

void num_mat_prod (double num, matrix A, matrix B){
    int i,j;
    int M = A.m;
    int N = A.n;
    matriz_nula(B);
    for (i=0;i<M;i++){
        for (j=0;j<N;j++){
            B.datos[i][j]= num * A.datos[i][j];
        }
    }
}

double norma_euclidea (matrix mat){
    int i,j;
    int M = mat.m;
    int N = mat.n;
    double norma=0.0;
    for (i=0;i<M;i++){
        for (j=0;j<N;j++){
            norma += powl(mat.datos[i][j],2);
        }
    }
    norma = sqrtl(norma);
    return norma;
}

