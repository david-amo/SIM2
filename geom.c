//
// Created by David on 06/02/2023.
//

#include "geom.h"

void leer_geom(char* nom_fichero, mol* molecula){
    FILE* arch;
    int N_atoms;
    char coment[200];
    arch= fopen(nom_fichero, "rt");
    int i;
    char simb[2];
    double P_X, P_Y, P_Z;

    fscanf(arch,"%d",&N_atoms);
    molecula->N = N_atoms;
    molecula->arr_atoms = malloc(N_atoms* sizeof(atom));
    fscanf(arch,"%s",coment);
    for (i=0;i<N_atoms;i++){
        fscanf(arch, "%s %lf %lf %lf",simb, &P_X,&P_Y,&P_Z);
        printf("%lf %lf %lf\n",P_X,P_Y,P_Z);

        molecula->arr_atoms[i].pos[0]=P_X;
        molecula->arr_atoms[i].pos[1]=P_Y;
        molecula->arr_atoms[i].pos[2]=P_Z;




        molecula->arr_atoms[i].M = M_a(simb);
        molecula->arr_atoms[i].Z = Z(simb);
        molecula->arr_atoms[i].id_atom =i;
    }
    fclose(arch);
}

void impr_mol (mol molecula){
    int i;
    printf("El número de átomos es %d\n",molecula.N);
    printf("#Atomo   Z   Masa         x          y           z\n");
    printf("---------------------------------------------------------\n");
    for (i=0;i<molecula.N;i++){
        printf("   %lf   %lf    %lf\n",molecula.arr_atoms[i].pos[0],molecula.arr_atoms[i].pos[1],molecula.arr_atoms[i].pos[2]);
    }
    printf("\n");
}

int cons_sets_SEA(mol** array_SEAS, mol molecula){
    //Voy a modificar la dirección a la que apunta el puntero mol* array_SEAS. Por lo tanto, tengo que pasar un puntero a un puntero: mol** array_SEAS

    int N_atoms = molecula.N;
    int index_atomo, index_set;
    matrix MDIN;
    matriz_inicializar(&MDIN,N_atoms,N_atoms);
    cons_MDIN(MDIN,molecula);
    int num_sets=1;
    int pert_SEA_previo=0;

    //El primer átomo declarado en la molécula constituye el primer set de átomos simétricamente equivalentes.
    (*array_SEAS)[0].arr_atoms[0].id_atom=molecula.arr_atoms[0].id_atom;
    (*array_SEAS)[0].arr_atoms[0].Z=molecula.arr_atoms[0].Z;
    (*array_SEAS)[0].arr_atoms[0].M=molecula.arr_atoms[0].M;
    (*array_SEAS)[0].arr_atoms[0].pos[0]=molecula.arr_atoms[0].pos[0];
    (*array_SEAS)[0].arr_atoms[0].pos[1]=molecula.arr_atoms[0].pos[1];
    (*array_SEAS)[0].arr_atoms[0].pos[2]=molecula.arr_atoms[0].pos[2];
    int set_de_equivalencia;


    for (index_atomo = 1; index_atomo < N_atoms; index_atomo++) {
        pert_SEA_previo = 0;
        for (index_set = 0; index_set < num_sets; index_set++) {
            pert_SEA_previo = comp_SEA(index_atomo, (*array_SEAS)[index_set].arr_atoms[0].id_atom, MDIN); //Si

            if (pert_SEA_previo==1){
                set_de_equivalencia = index_set;
                break;
            }

        }

        if (pert_SEA_previo==0){
            num_sets++;
            *array_SEAS=realloc(*array_SEAS, (num_sets)* sizeof(mol));
            (*array_SEAS)[num_sets-1].arr_atoms = malloc(sizeof(atom));
            (*array_SEAS)[num_sets-1].N = 1;
            (*array_SEAS)[num_sets-1].arr_atoms[0].id_atom=molecula.arr_atoms[index_atomo].id_atom;
            (*array_SEAS)[num_sets-1].arr_atoms[0].Z=molecula.arr_atoms[index_atomo].Z;
            (*array_SEAS)[num_sets-1].arr_atoms[0].M=molecula.arr_atoms[index_atomo].M;
            (*array_SEAS)[num_sets-1].arr_atoms[0].pos[0]=molecula.arr_atoms[index_atomo].pos[0];
            (*array_SEAS)[num_sets-1].arr_atoms[0].pos[1]=molecula.arr_atoms[index_atomo].pos[1];
            (*array_SEAS)[num_sets-1].arr_atoms[0].pos[2]=molecula.arr_atoms[index_atomo].pos[2];
        }
        else{
            ad_atom_a_mol(&(*array_SEAS)[set_de_equivalencia],molecula.arr_atoms[index_atomo]);
        }
    }
    //Finalmente, calculo el vector posición del centro de masas para cada SEA y sitúo el origen de coordenadas de cada set en
    //el correspondiente CDM;

    double* pos_CDM = malloc(3* sizeof(double));
    /*for (index_set=0;index_set<num_sets;index_set++){
        calc_cdm((*array_SEAS)[index_set],pos_CDM);
        origen_cdm(pos_CDM,(*array_SEAS)[index_set]);
    }*/






    return num_sets;
}

void cons_tensor_inerc(matrix tensor, mol molecula){
    int N_atms = molecula.N;
    int atom;
    double m;
    double comp_x,comp_y,comp_z;
    matriz_nula(tensor);

    enum ccart {x,y,z};
    for(atom=0; atom < N_atms; atom++){
        m = molecula.arr_atoms[atom].M;
        comp_x = molecula.arr_atoms[atom].pos[x];
        comp_y = molecula.arr_atoms[atom].pos[y];
        comp_z = molecula.arr_atoms[atom].pos[z];

        tensor.datos[x][x]+= m*(comp_y*comp_y+comp_z*comp_z);
        tensor.datos[y][y]+= m*(comp_x*comp_x+comp_z*comp_z);
        tensor.datos[z][z]+= m*(comp_x*comp_x+comp_y*comp_y);

        tensor.datos[y][x]-= m*(comp_x*comp_y);
        tensor.datos[z][x]-= m*(comp_x*comp_z);
        tensor.datos[z][y]-= m*(comp_z*comp_y);
    }

    tensor.datos[x][y]=tensor.datos[y][x];
    tensor.datos[x][z]=tensor.datos[z][x];
    tensor.datos[y][z]=tensor.datos[z][y];



}


void ad_atom_a_mol (mol* molecula, atom atomo){ //Añade un nuevo átomo a una estructura molécula



    int N_atomos = molecula->N;

    N_atomos++;

    molecula->arr_atoms=realloc(molecula->arr_atoms,N_atomos* sizeof(atom));

    molecula->N = N_atomos;
    molecula->arr_atoms[N_atomos-1].id_atom=atomo.id_atom;
    molecula->arr_atoms[N_atomos-1].M=atomo.M;
    molecula->arr_atoms[N_atomos-1].Z=atomo.Z;
    molecula->arr_atoms[N_atomos-1].pos[0] = atomo.pos[0];
    molecula->arr_atoms[N_atomos-1].pos[1] = atomo.pos[1];
    molecula->arr_atoms[N_atomos-1].pos[2] = atomo.pos[2];
}

void cons_MDIN (matrix MDIN, mol molecula){
    int N_atom = molecula.N;
    size_t i,j;
    for (i=0;i<N_atom;i++){
        for (j=0;j<N_atom;j++){
            MDIN.datos[i][j]= dist(molecula.arr_atoms[i].pos,molecula.arr_atoms[j].pos);
        }
    }
}

int comp_SEA(int atom1, int atom2, matrix MDIN){
    int dist_i, dist_j;
    int i;
    int* arr_equiv = malloc(MDIN.m* sizeof(int));
    for (i=0;i<MDIN.m;i++){
        arr_equiv[i]=-1;
    }
    int k;
    size_t N = MDIN.m;
    bool asignado;
    int equiv_dist;
    for (dist_i = 0; dist_i < N; dist_i++) {
        for (dist_j = 0; dist_j < N; dist_j++) {
            asignado = false;
            for (k = 0; k < dist_i; k++) {
                if (arr_equiv[k] == dist_j) {
                    asignado = true;
                    break;
                }
            }
            if (asignado == false) {
                equiv_dist = comp_doubles((double)MDIN.datos[dist_i][atom1], (double)MDIN.datos[dist_j][atom2], TOL_GEOM);
                if (equiv_dist == 1) {
                    arr_equiv[dist_i] = dist_j;
                }
            }
        }
    }

    for (i = 0; i < N; i++) {
        if(arr_equiv[i]==-1){
            //free(arr_equiv);
            return 0; // No son átomos simétricamente equivalentes
        }
    }
    //free(arr_equiv);
    return 1; // Son átomos simétricamente equivalentes;
}

int comp_doubles (double val1, double val2, double toler){
    int equiv;
    if (fabs(val1-val2)<=toler){
        return 1;
    }
    else{
        return 0;
    }
}

double M_a (char* simbolo){
    if (strcmp(simbolo,"H")==0)
        return 1.007825;
    if (strcmp(simbolo,"C")==0)
        return 12.000000;
    if (strcmp(simbolo,"N")==0)
        return 14.003074;
    if(strcmp(simbolo,"O")==0)
        return 15.9994;
    if(strcmp(simbolo,"F")==0)
        return 18.998403;

    else
        return 0;
}

int Z (char* simbolo){
    if (strcmp(simbolo,"H")==0)
        return 1;
    if (strcmp(simbolo,"N")==0)
        return 7;
    if (strcmp(simbolo,"C")==0)
        return 6;
    if (strcmp(simbolo,"O")==0)
        return 8;
    if (strcmp(simbolo,"F")==0)
        return 9;

    else
        return 0;
}

double dist (double* vecA, double* vecB){
    int i;
    double dst=0;

    for (i=0;i<3;i++){
        dst += pow(vecA[i]-vecB[i],2);
    }

    dst = sqrt(dst);

    return dst;
}



void realloc_SEAS (mol** array_SEAS, size_t nuevo_tam){
    *array_SEAS = realloc(*array_SEAS,nuevo_tam*sizeof(mol));
    if (*array_SEAS==NULL){
        printf("Error critico en la asignacion de memoria\n");
        exit(-1);
    }
    //array_SEAS[num_sets-1].arr_atoms = malloc(sizeof(atom));

}

void calc_cdm(mol molecula, double* r_cm){
    double M = 0;
    double m;
    int atomo;
    int coord;
    r_cm[0]=0.0;
    r_cm[1]=0.0;
    r_cm[2]=0.0;
    for (atomo=0;atomo<molecula.N;atomo++){
        m = molecula.arr_atoms[atomo].M;
        for (coord=0;coord<3;coord++){
            r_cm[coord] += m*molecula.arr_atoms[atomo].pos[coord];
        }
        M += m;
    }


    for (coord=0;coord<3;coord++){
        r_cm[coord] /= M;
    }
}

void origen_cdm(double* cdm, mol molecula){ //Traslada una molecula al centro de masas
    int N=molecula.N;
    int coord;
    int atomo;
    for (atomo=0;atomo<N;atomo++){
        for (coord=0;coord<3;coord++){
            molecula.arr_atoms[atomo].pos[coord] -= cdm[coord];
        }
    }

}

bool test_equiv_mol (mol molA, mol molB){
    size_t atomA, atomB;
    size_t N = molA.N;
    size_t i,k;
    bool asignado;
    int test_x, test_y,test_z;
    int* arr_equiv = malloc(molA.N* sizeof(int));
    for (i=0;i<molA.N;i++){
        arr_equiv[i]=-1;
    }

    if (molA.N != molB.N){
        return false;
    }

    for (atomA = 0; atomA<N; atomA++) {
        for (atomB = 0; atomB < N; atomB++) {
            asignado = false;
            for (k = 0; k < atomA; k++) {
                if (arr_equiv[k] == atomB) {
                    asignado = true;
                    break;
                }
            }

            if (asignado == false) {
                test_x = comp_doubles(molA.arr_atoms[atomA].pos[0], molB.arr_atoms[atomB].pos[0], TOL_GEOM);
                test_y = comp_doubles(molA.arr_atoms[atomA].pos[1], molB.arr_atoms[atomB].pos[1], TOL_GEOM);
                test_z = comp_doubles(molA.arr_atoms[atomA].pos[2], molB.arr_atoms[atomB].pos[2], TOL_GEOM);

                if ((test_x==1) && (test_y==1) && (test_z==1)) {
                    arr_equiv[atomA] = (int)atomB;
                    break;
                }
            }

        }
    }

    for (i=0;i<N;i++){
        if (arr_equiv[i]==-1){
            free (arr_equiv);
            return false;
        }
    }
    free(arr_equiv);
    return true;
}

void rot_mol (mol molA, mol* molB, double* eje, double theta){

    int N = molA.N;
    double cos_theta = cos(theta);
    double sin_theta = sin(theta);
    double uno_m_cos_theta = 1-cos_theta;
    //Primero, normalizar el vector;
    double norm = eje[0]*eje[0]+eje[1]*eje[1]+eje[2]*eje[2];
    norm = sqrt(norm);
    double x = eje[0]/norm;
    double y = eje[1]/norm;
    double z = eje[2]/norm;
    size_t atom;

    matrix R;
    matriz_inicializar(&R,3,3);

    R.datos[0][0]=cos_theta+x*x*uno_m_cos_theta;
    R.datos[1][1]=cos_theta+y*y*uno_m_cos_theta;
    R.datos[2][2]=cos_theta+z*z*uno_m_cos_theta;

    R.datos[1][0]=y*x*uno_m_cos_theta+z*sin_theta;
    R.datos[2][0]=z*x*uno_m_cos_theta-y*sin_theta;

    R.datos[0][1]=x*y*uno_m_cos_theta-z*sin_theta;
    R.datos[2][1]=z*y*uno_m_cos_theta+x*sin_theta;

    R.datos[0][2]=x*z*uno_m_cos_theta+y*sin_theta;
    R.datos[1][2]=y*z*uno_m_cos_theta-x*sin_theta;

    int i;
    molB->N = N;
    for (i=0;i<N;i++){
        (*molB).arr_atoms[i]=molA.arr_atoms[i];
    }
    int coord;

    int j,b;
    //Transforma cada átomo en la molécula
    for (atom=0;atom<molA.N;atom++){
        molB->arr_atoms[atom].pos[0]=0;
        molB->arr_atoms[atom].pos[1]=0;
        molB->arr_atoms[atom].pos[2]=0;
        for (coord=0;coord<3;coord++){
                for (j=0;j<3;j++){
                    molB->arr_atoms[atom].pos[coord]+=R.datos[coord][j]*molA.arr_atoms[atom].pos[j];
                }

        }


    }
}

void reflex_mol(mol molA, mol* molB, double* v_normal){
    matrix Mreflex;
    matriz_inicializar(&Mreflex,3,3);
    enum ccart{x,y,z};
    int N = molA.N;
    int atom;
    int coord1,coord2;
    vector temp;
    vector_inic(&temp,3);


    Mreflex.datos[x][x]= 1-2*v_normal[x]*v_normal[x];
    Mreflex.datos[y][y]= 1-2*v_normal[y]*v_normal[y];
    Mreflex.datos[z][z]= 1-2*v_normal[z]*v_normal[z];
    for (coord1=y;coord1<=z;coord1++){
        for (coord2=x;coord2<coord1;coord2++){
            Mreflex.datos[coord1][coord2]= -2*v_normal[coord1]*v_normal[coord2];
            Mreflex.datos[coord2][coord1]=Mreflex.datos[coord1][coord2];
        }
    }
    int i,j,b;
    molB->N = N;
    for (i=0;i<N;i++){
        (*molB).arr_atoms[i]=molA.arr_atoms[i];
    }


    int coord;


    //Transforma cada átomo en la molécula
    for (atom=0;atom<molA.N;atom++){
        molB->arr_atoms[atom].pos[0]=0;
        molB->arr_atoms[atom].pos[1]=0;
        molB->arr_atoms[atom].pos[2]=0;
        for (coord=0;coord<3;coord++){
            for (j=0;j<3;j++){
                molB->arr_atoms[atom].pos[coord] += Mreflex.datos[coord][j] * molA.arr_atoms[atom].pos[j];
            }

        }
    }

}

int hallar_ejes_rot_prop(elem_sim** array_elems, mol set_SEAS, mol molecula, int* n_total_elems, double* CDM_SET_REF_CDM_MOL){
    matrix tensor_inercia, eigenvectores_inercia;
    int num_ejes = 0;
    elem_sim* arr_ejes_candidatos = NULL;
    //Cambiar código. Primero generar todos los posibles candidatos. A continuación, comprobar la existencia
    // de los elementos de simetría en la molécula.
    int index_set;
    int i;
    int N_atoms_set = set_SEAS.N;
    int N_atoms_molec = molecula.N;
    int orden;
    enum m_inercia{A,B,C};
    int eje;
    int num_ejes_candidatos = 0;
    mol mol_rotada;
    double prod_escalar;
    double cdm[3];

    mol_rotada.N = N_atoms_set;
    mol_rotada.arr_atoms = malloc(N_atoms_molec* sizeof(atom));

    mol set_ORG_CDMSET;
    set_ORG_CDMSET.N=set_SEAS.N;
    set_ORG_CDMSET.arr_atoms = malloc (N_atoms_set* sizeof(atom));

    memcpy(set_ORG_CDMSET.arr_atoms,set_SEAS.arr_atoms,set_ORG_CDMSET.N* sizeof(atom));
    calc_cdm(set_ORG_CDMSET,cdm);
    origen_cdm(cdm,set_ORG_CDMSET);
    double eje_candidato_ref_CDM_MOL[3];

    if (set_SEAS.N > 1) {

        matriz_inicializar(&tensor_inercia, 3, 3);
        matriz_inicializar(&eigenvectores_inercia, 3, 3);
        cons_tensor_inerc(tensor_inercia, set_ORG_CDMSET);
        double *eigenvalores_inerc = malloc(3 * sizeof(double));
        eigen_jacobi(tensor_inercia, eigenvalores_inerc, eigenvectores_inercia, false);
        //Distinguir casos de arreglo poliédrico y poligonal
        // ¿Los átomos constituyen un arreglo planar?
        double I_a, I_b, I_c;
        I_a = eigenvalores_inerc[A];
        I_b = eigenvalores_inerc[B];
        I_c = eigenvalores_inerc[C];
        printf("autovectores del tensor de inercia\n");
        matriz_imprimir(eigenvectores_inercia);
        double *eje_candidato = malloc(3 * sizeof(double));
        printf ("autovalores del tensor de inercia\n");
        printf("Ia:%lf\nIb:%lf\nIc:%lf\n",I_a,I_b,I_c);



        int elem_arr = 0;
        int j;
        bool preexistente=false;

        if ((comp_doubles(I_a, 0.0, TOL_DOUBLES) == 1) && (comp_doubles(I_b, I_c, TOL_DOUBLES))) {
            printf("El set describe un arreglo lineal\n");
            for (eje = A; eje <= C; eje++) {
                cop_col_a_vect(eigenvectores_inercia, eje_candidato, eje);
                rest_vector(CDM_SET_REF_CDM_MOL,eje_candidato,eje_candidato_ref_CDM_MOL);
                for (orden = 2; orden <= MAX_ORDER; orden++) {
                    printf("probando orden %d\n", orden);
                    rot_mol(molecula, &mol_rotada, eje_candidato_ref_CDM_MOL, 2 * M_PI / orden);
                    if (test_equiv_mol(molecula, mol_rotada) == true) {
                        (*n_total_elems)++;
                        *array_elems = realloc(*array_elems, *(n_total_elems) * sizeof(elem_sim));
                        (*array_elems)[*n_total_elems - 1].ref = 1;
                        (*array_elems)[*n_total_elems - 1].order = orden;
                        memcpy((*array_elems)[*n_total_elems - 1].vec_dic, eje_candidato_ref_CDM_MOL, 3 * sizeof(double));
                        printf("IDENTIFICADO EJE DE ROTACION DE ORDEN %d\n", orden);
                    }
                }
            }
        }
        else if (comp_doubles(I_a + I_b, I_c, TOL_DOUBLES) == 1) { //CASO I : Arreglo poligonal
            printf("El SEA conforma un arreglo poligonal");
            cop_col_a_vect(eigenvectores_inercia, eje_candidato, 2);
            rest_vector(CDM_SET_REF_CDM_MOL,eje_candidato,eje_candidato_ref_CDM_MOL);
            if (comp_doubles(I_a, I_b, TOL_DOUBLES) == 1) { //CASO Ia: Arreglo poligonal regular
                printf(" regular\n");
                //Los ejes de rotación posibles son perpendiculares (eigenvector Ic) al plano determinado por el polígono regular con
                //orden el número de vértices y sus divisores
                preexistente=false;
                for (i = -N_atoms_set; i <= N_atoms_set; i++) {
                    if (i != 0 && abs(i) != 1) {
                        if (((N_atoms_set % i) == 0)) {
                            printf("Probando en arreglo poligonal con orden %d\n",i);
                            rot_mol(molecula, &mol_rotada, eje_candidato_ref_CDM_MOL, 2 * M_PI / i);
                            if (test_equiv_mol(mol_rotada, mol_rotada) == true) {
                                for (j=1;j<*n_total_elems;j++){
                                    if ((*array_elems)[j].ref==1&&(*array_elems)[j].order==i){
                                        prod_escalar=prod_esc(eje_candidato_ref_CDM_MOL,(*array_elems)[j].vec_dic);
                                        if (comp_doubles(prod_escalar,1,TOL_DOUBLES)){
                                            preexistente=true;
                                        }
                                    }
                                }
                                if (preexistente==false){
                                    (*n_total_elems)++;
                                    *array_elems = realloc(*array_elems, *(n_total_elems) * sizeof(elem_sim));
                                    (*array_elems)[*n_total_elems - 1].ref = 1;
                                    (*array_elems)[*n_total_elems - 1].order = i;
                                    memcpy((*array_elems)[*n_total_elems - 1].vec_dic, eje_candidato_ref_CDM_MOL, 3 * sizeof(double));
                                    printf("IDENTIFICADO EJE DE ROTACION DE ORDEN %d\n", i);
                                }
                            }
                        }
                    }
                }
            }
        }
        else{//arreglo poliedrico

        }
    }

//            } else {
//                printf(" irregular\n"); //Caso Ib: Arreglo poligonal irregular
//                //Los ejes de rotación posibles son perpendiculares (eigenvector Ic) al plano determinado por el polígono regular con
//                //orden los divisores del número de vértices
//                for (i = 2; i < N_atoms_set; i++) {
//                    printf("Natoms:%d\n",N_atoms_set);
//                    if ((N_atoms_set % i) == 0) {
//                        arr_ejes_candidatos = realloc(arr_ejes_candidatos, num_ejes_candidatos * sizeof(elem_sim));
//                        arr_ejes_candidatos[num_ejes_candidatos - 1].ref = 1;
//                        arr_ejes_candidatos[num_ejes_candidatos - 1].order = i;
//                        cop_col_a_vect(tensor_inercia, arr_ejes_candidatos[num_ejes_candidatos - 1].vec_dic, C);
//                        printf("probando eje de rotacion propio de orden %d...\n",i);
//                        rot_mol(set_SEAS, &mol_rotada, cand_rot, 2 * M_PI / i);
//                        *//*if (test_equiv_mol(mol_rotada, set_SEAS) == true) {
//                            printf("CONFIRMADO EJE DE ROTACION DE ORDEN %d\n", i);
//                        }*//*
//                    }
//                }
//
//            }
//
//        }
//        else{
//            //De no cumplirse ninguna de las condiciones anteriores, el conjunto de átomos simétricamente equivalentes
//            //describe un poliedro
//            if ((I_a<I_b)&&(comp_doubles(I_b, I_c, TOL_DOUBLES) == 1)){
//                printf ("trompo simetrico prolate\n");
//                // El autovector asociado a I_A es un candidato a eje de rotacion propio;
//                cop_col_a_vect(eigenvectores_inercia, cand_rot, 0);
//                for (orden=-MAX_ORDER;orden<=MAX_ORDER;orden++){
//
//                    printf("Probando con orden %d\n",orden);
//                    rot_mol(molecula, &mol_rotada, cand_rot, 2 * M_PI / orden);
//                    if (test_equiv_mol(mol_rotada, molecula) == true) {
//                        printf("CONFIRMADO EJE DE ROTACION DE ORDEN %d\n", orden);
//                    }
//                }
//            }
//            else if ((I_b<I_c)&&(comp_doubles(I_a, I_b, TOL_DOUBLES))){
//                printf("trompo simetrico oblate\n");
//
//            }
//
//        }
//
//    }
//    else{
//        // Si el conjunto analizado consta de un único átomo, no proporciona ninguna información sobre los elementos
//        // de simetría presentes en la molécula.
//
//        return 0;
//    }
//        for (i = 0; i < num_ejes_candidatos; i++) {
//            rot_mol(molecula, &mol_rotada, arr_ejes_candidatos[i].vec_dic, 2 * M_PI / arr_ejes_candidatos[i].order);
//            if (test_equiv_mol(molecula, mol_rotada) == true) {
//                (*n_total_elems)++;
//                (*array_elems) = realloc(*array_elems, *n_total_elems*sizeof(atom));
//                (*array_elems)[*n_total_elems - 1] = arr_ejes_candidatos[i];
//            }
//        }
//        return 0;
//    }
    free_mol(&set_ORG_CDMSET);
    free_mol(&mol_rotada);
}
int hallar_c2_perps(elem_sim** array_elems, mol set_SEAS, mol molecula, int* n_total_elems){
    int atom1;
    int atom2;
    int N_atoms_set = set_SEAS.N;
    int N_atoms_molec = molecula.N;
    double cand_eje_A[3];
    enum ccart{x,y,z};
    enum ccart coord;
    mol mol_rotada;
    int j;
    double prod_escalar;
    bool preexistente;
    for (atom1=0;atom1<N_atoms_set;atom1++){
        for (atom2=atom1+1;atom2<N_atoms_set;atom2++){
            sum_vector(set_SEAS.arr_atoms[atom1].pos,set_SEAS.arr_atoms[atom2].pos,cand_eje_A);
            for (coord=x;coord<z;coord++){
                cand_eje_A[coord]/=2;
            }
            normaliza_vec(cand_eje_A);
            mol_rotada.N = N_atoms_set;
            mol_rotada.arr_atoms = malloc(N_atoms_molec* sizeof(atom));
            rot_mol(molecula, &mol_rotada, cand_eje_A, M_PI);
            if (test_equiv_mol(molecula, mol_rotada) == true){
                preexistente=false;

                for (j=1;j<*n_total_elems;j++){
                    if ((*array_elems)[j].ref==5&&(*array_elems)[j].order==2){
                        prod_escalar=prod_esc(cand_eje_A,(*array_elems)[j].vec_dic);
                        if (comp_doubles(prod_escalar,1,TOL_DOUBLES)==1||comp_doubles(prod_escalar,-1,TOL_DOUBLES)==1){
                            preexistente=true;
                        }
                    }
                }

                if (preexistente==false){
                    (*n_total_elems)++;
                    *array_elems = realloc(*array_elems, *(n_total_elems) * sizeof(elem_sim));
                    (*array_elems)[*n_total_elems - 1].ref = 5;
                    (*array_elems)[*n_total_elems - 1].order = 2;
                    memcpy((*array_elems)[*n_total_elems - 1].vec_dic, cand_eje_A, 3 * sizeof(double));
                }

            }

        }

    }


}

int hallar_inv (elem_sim** array_elems,mol molecula,int* n_total_elems){
    int N = molecula.N;
    enum ccart{x,y,z};
    enum ccart coord;
    int atomo;
    int tiene_inv=0;

    mol mol_invertida;
    mol_invertida.N = N;
    mol_invertida.arr_atoms = malloc(N* sizeof(atom));

    double c_inv[3];
    calc_cdm(molecula,c_inv);
    elem_sim* temp_arr_inv = NULL;


    memcpy(mol_invertida.arr_atoms,molecula.arr_atoms,N* sizeof(atom));
    for (atomo=0;atomo<N;atomo++){
        for (coord=x;coord<=z;coord++){
            mol_invertida.arr_atoms[atomo].pos[coord] *= -1;
        }
    }

    if(test_equiv_mol(molecula,mol_invertida)==true){
        (*n_total_elems)++;
        *array_elems = realloc(*array_elems, *(n_total_elems) * sizeof(elem_sim));
        (*array_elems)[*n_total_elems - 1].ref = 7;
        (*array_elems)[*n_total_elems - 1].order = 0;
        memcpy((*array_elems)[*n_total_elems - 1].vec_dic, c_inv, 3 * sizeof(double));
    }
    else{
        printf("La molécula no tiene centro de inversion\n");
        return 0;
    }
    return 1;







}

int hallar_planos_ref (elem_sim** array_elems, mol set_SEAS, mol molecula,int* n_total_elems){
    int N_atoms_SEA = set_SEAS.N;
    int N_atoms_molec = molecula.N;
    int atomA, atomB;
    int num_planos=0;

    double norm_plano[3];
    mol mol_reflejada;
    alloc_mol(N_atoms_molec,&mol_reflejada);
    elem_sim* temp_arr = NULL;
    int j;
    double prod_escalar;
    bool preexistente;
    int eje_princ = pos_arr_cmax(array_elems,*n_total_elems);
    double vec_eje_princ[3];
    memcpy(vec_eje_princ,(*array_elems)[eje_princ].vec_dic,3* sizeof(double));
    int ref=2; //2 para plano horizontal y 3 para plano vertical/diedro
    double test_horizontal;


    if (N_atoms_SEA > 1) {
        printf("IDENTIFICANDO PLANOS DE SIMETRIA....\n");
        for (atomA = 0; atomA < N_atoms_SEA; atomA++) {
            for (atomB = atomA + 1; atomB < N_atoms_SEA; atomB++) {
                rest_vector(set_SEAS.arr_atoms[atomA].pos, set_SEAS.arr_atoms[atomB].pos, norm_plano);
                normaliza_vec(norm_plano);
                reflex_mol(molecula, &mol_reflejada, norm_plano);
                preexistente = false;
                if (test_equiv_mol(mol_reflejada, molecula) == 1) {
                    for (j = 1; j < *n_total_elems; j++) {
                        if ((*array_elems)[j].ref == 2 || (*array_elems)[j].ref == 3) {
                            prod_escalar = prod_esc(norm_plano, (*array_elems)[j].vec_dic);
                            if (comp_doubles(fabs(prod_escalar), 1, TOL_DOUBLES)) {
                                preexistente = true;
                            }
                        }
                    }
                    if (preexistente == false) {
                        printf("Hallado plano de simetria\n");
                        num_planos++;
                        *n_total_elems = (*n_total_elems) + 1;
                        int NELS = *n_total_elems;
                        temp_arr = realloc(*array_elems, NELS * sizeof(elem_sim));
                        if (*array_elems == NULL) {
                            printf("ERROR DE MEMORIA\n");
                            exit(-1);
                        } else {
                            *array_elems = temp_arr;
                        }

                        test_horizontal = prod_esc(vec_eje_princ, norm_plano);
                        printf("Vector eje principal:%lf\t%lf\t%lf\n", vec_eje_princ[0], vec_eje_princ[1],
                               vec_eje_princ[2]);
                        printf("TEST HORIZONTAL:%lf\n", test_horizontal);
                        if (test_horizontal == 1) {
                            ref = 2; //plano horizontal
                        } else {
                            ref = 3; //plano vertical/diedro
                        }
                        (*array_elems)[NELS - 1].order = 0;
                        (*array_elems)[NELS - 1].ref = ref;
                        (*array_elems)[NELS - 1].vec_dic[0] = norm_plano[0];
                        (*array_elems)[NELS - 1].vec_dic[1] = norm_plano[1];
                        (*array_elems)[NELS - 1].vec_dic[2] = norm_plano[2];
                    }

                }
            }
        }


        printf("SE HAN ENCONTRADO %d PLANOS DE SIMETRIA EN TOTAL\n",num_planos);
        //free(norm_plano);
        free_mol(&mol_reflejada);
        int i;
        return num_planos;
    }

    else{
        printf("El conjunto solo incluye un atomo y no proporciona informacion sobre la simetria molecular\n");
        return 0;
    }
}


void vector_reflex (vector V, vector P){
    matrix Mreflex;
    matriz_inicializar(&Mreflex,3,3);
    enum ccart{x,y,z};
    vector temp;
    vector_inic(&temp,3);

    Mreflex.datos[x][x]= 1-2*x*x;
    Mreflex.datos[y][y]= 1-2*y*y;
    Mreflex.datos[z][z]= 1-2*z*z;

    Mreflex.datos[y][x]=-2*x*y;
    Mreflex.datos[x][y]=Mreflex.datos[y][x];

    Mreflex.datos[z][x]=-2*x*z;
    Mreflex.datos[x][z]=Mreflex.datos[z][x];

    Mreflex.datos[z][y]=-2*x*z;
    Mreflex.datos[y][z]=Mreflex.datos[z][y];

    printf("Matriz de reflexion\n");
    matriz_imprimir(Mreflex);
    vector_transf(Mreflex,V,temp);
    V.vector[x]=temp.vector[x];
    V.vector[y]=temp.vector[y];
    V.vector[z]=temp.vector[z];

}

void vector_transf (matrix T, vector A, vector A_T){
    size_t i,j;
    for (i=0;i<3;i++){
        for (j=0;j<3;j++){
            A_T.vector[i]+=T.datos[i][j]*A.vector[j];
        }
    }
}

void vector_inic(vector* Vec, size_t m) {
    double* ptr;


    ptr = malloc(m * sizeof(double));

    if (ptr == NULL) {
        free(ptr);
        exit(-1);

    }
    Vec->m = m;
    Vec->vector=ptr;
}
void vector_norm (vector V){
    double normal = 0;
    size_t i;

    for (i=0;i<V.m;i++)
        normal += pow(V.vector[i],2);

    normal = sqrt(normal);

    for (i=0;i<V.m;i++)
        V.vector[i] /= normal;
}

void cop_col_a_vect (matrix Mat, double* V, int col){

    int m = Mat.m; //Número de filas
    int i;
    for (i=0;i<m;i++)
        V[i]=Mat.datos[i][col];
}

void rest_vector (double* VecA, double* VecB, double* VecC){
    int coord;
    for (coord=0;coord<3;coord++){
        VecC[coord]=VecA[coord]-VecB[coord];
    }
}

void sum_vector (double* VecA, double* VecB, double* VecC){
    int coord;
    for (coord=0;coord<3;coord++){
        VecC[coord]=VecA[coord]+VecB[coord];
    }
}

void normaliza_vec (double* Vec){
    enum ccart {x,y,z};
    int coord;
    double cte_norm = Vec[x]*Vec[x] + Vec[y]*Vec[y] + Vec[z]*Vec[z];
    cte_norm = sqrt(cte_norm);
    for (coord=x;coord<=z;coord++){
        Vec[coord] /= cte_norm;
    }
}

void impr_elems_sim(elem_sim* arr_elems, int num_total_elems_sim){
    int i;
    printf("\n");
    printf("numero de elementos:%d\n",num_total_elems_sim);
    printf("************************************\n");
    printf("*      ELEMENTOS DE SIMETRIA       *\n");
    printf("************************************\n");
    printf("TIPO                    ORDEN                         X           Y             Z\n");
    printf("---------------------------------------------------------------------------------\n");
    int tipo;
    char descr[25];
    double x,y,z;
    int orden=0;
    for (i=0;i<num_total_elems_sim;i++){
        orden = arr_elems[i].order;
        tipo = arr_elems[i].ref;
        if (tipo == 0) {
            strcpy(descr, "IDENTIDAD              ");
        } else if (tipo == 1) {
            strcpy(descr, "EJE                    ");
        } else if (tipo == 2) {
            strcpy(descr, "\u03c3 (h)             ");
        } else if (tipo == 3) {
            strcpy(descr, "\u03c3 (v/d)           ");
        } else if (tipo == 5) {
            strcpy(descr, "C2                     ");
        }else if (tipo==7){
            strcpy(descr, "i                      ");
        }

        x = arr_elems[i].vec_dic[0];
        y = arr_elems[i].vec_dic[1];
        z = arr_elems[i].vec_dic[2];
        printf("%s\t\t%02.1d                %lf        %lf           %lf\n", descr,orden,x,y,z);
    }
}

double prod_esc (double* VecA, double* VecB){
    double prodesc = 0.0;
    enum ccart{x,y,z};
    int coord;

    for (coord=x;coord<=z;coord++)
        prodesc += VecA[coord]*VecB[coord];

    return prodesc;
}

int testSIMELS (elem_sim** array_elems, mol set_SEAS, mol molecula,int* n_total_elems){

    int incremento=10;
    (*n_total_elems)=incremento;
    *array_elems = realloc(*array_elems,(*n_total_elems)* sizeof(elem_sim));
    int i;
    for (i=1;i<incremento;i++){
        (*array_elems)[i].order =1;
        (*array_elems)[i].ref=2;
        (*array_elems)[i].vec_dic[0]=1.0;
        (*array_elems)[i].vec_dic[1]=2.0;
        (*array_elems)[i].vec_dic[2]=3.0;
    }
  return 0;
}

void free_mol(mol* molecula){
    free(molecula->arr_atoms);
    molecula->arr_atoms = NULL;
}

int alloc_mol(int N_at, mol* molecula){
    molecula->N = N_at;
    molecula->arr_atoms = malloc(N_at* sizeof(atom));
    if (molecula->arr_atoms==NULL)
        return 0;
    else
        return 1;
}

void copia_mol (mol destino, mol origen){
    destino.N = origen.N;
    int i;
    for (i=0;i<destino.N;i++){

    }
}

void rest_vect (double vecA[3], double vecB[3], double vecR[3]){
    int i;

    for (i=0;i<3;i++){
        vecR[i]=vecA[i]-vecB[i];
    }
}

int pos_arr_cmax(elem_sim** array_elems, int n_total_elems){
    int i;
    int index;
    int orden_max=0;
    for (i=0;i<n_total_elems;i++){
        if ((*array_elems)[i].ref==1 && (*array_elems)[i].order>orden_max){
            orden_max = (*array_elems)[i].order>orden_max;
            index = i;
        }
    }
return index;
}

void prodvect (double* vecA, double* vecB, double* vecP){
    vecP[0]=vecA[1]*vecB[2]-vecA[2]*vecB[1];
    vecP[1]=vecA[2]*vecB[0]-vecA[0]*vecB[2];
    vecP[2]=vecA[0]*vecB[1]-vecA[1]*vecB[0];
}

