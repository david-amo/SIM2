#include <stdio.h>
#include "geom.h"
#include "jacobi.h"

int main() {
    double array[3];
    double dec;
    mol molecula;
    leer_geom("benceno.xyz",&molecula);


    matrix MDIN;
    matriz_inicializar(&MDIN,molecula.N,molecula.N);
    cons_MDIN(MDIN,molecula);
    printf("matriz de distancias internucleares\n");
    matriz_imprimir(MDIN);






    mol* array_SEAS;
    array_SEAS= malloc(sizeof(mol));
    array_SEAS[0].N=1;
    array_SEAS[0].arr_atoms = malloc(sizeof(atom));



    double* cdm = malloc(3* sizeof(double));
    calc_cdm(molecula,cdm);
    origen_cdm(cdm,molecula);
    printf ("Molecula a analizar 1:\n");
    impr_mol(molecula);
    int num_sets=cons_sets_SEA(&array_SEAS,molecula);
    calc_cdm(molecula,cdm);
    printf("CDM MOLECULA------>    x:%lf\ty:%lf\tz:%lf\n",cdm[0],cdm[1],cdm[2]);


    int i;
    for (i=0;i<num_sets;i++){
        printf("CONJUNTO DE EQUIVALENCIA NUMERO %d\n",i);
        printf("-----------------------------------\n");
        impr_mol(array_SEAS[i]);
        calc_cdm(array_SEAS[i],cdm);
        printf("CDM SET %d------>    x:%lf\ty:%lf\tz:%lf\n",i,cdm[0],cdm[1],cdm[2]);
    }
    free(cdm);
    elem_sim* arr_elems = malloc (sizeof(elem_sim));

    double** arr_cdmsets ;
    arr_cdmsets = malloc(num_sets* sizeof(double*));
    for (i=0;i<num_sets;i++){
        arr_cdmsets[i]= malloc(3* sizeof(double));
    }

    for (i=0;i<num_sets;i++){
        calc_cdm(array_SEAS[i],arr_cdmsets[i]);
    }



    arr_elems[0].ref=0;
    arr_elems[0].order=0;
    arr_elems[0].vec_dic[0]=0.0;
    arr_elems[0].vec_dic[1]=0.0;
    arr_elems[0].vec_dic[2]=0.0;
    //impr_elems_sim(arr_elems,num_total_elementos);
    int ejes_rotp;
    printf("SIZE OF ATOM:%lu\n", sizeof(elem_sim));
    double* CDM_SET_REF_CDM_MOL = malloc(3* sizeof(double));


    int num_total_elementos = 1;
    for (i=0;i<num_sets;i++){
        printf("ESTUDIANDO EL CONJUNTO DE EQUIVALENCIA %d...\n",i);
        calc_cdm(array_SEAS[i],CDM_SET_REF_CDM_MOL);
        printf("EL CDM DEL SET NUMERO %d es:\nx:%lf\ty:%lf\tz:%lf\n",i,CDM_SET_REF_CDM_MOL[0],CDM_SET_REF_CDM_MOL[1],CDM_SET_REF_CDM_MOL[2]);
        hallar_ejes_rot_prop(&arr_elems,array_SEAS[i],molecula,&num_total_elementos,CDM_SET_REF_CDM_MOL);
        hallar_planos_ref(&arr_elems,array_SEAS[i],molecula,&num_total_elementos);
        if (array_SEAS[i].N>1)
        hallar_c2_perps(&arr_elems,array_SEAS[i],molecula,&num_total_elementos);
    }
    double cand_plano_horizontal[3];
    prodvect(molecula.arr_atoms[0].pos,molecula.arr_atoms[1].pos,cand_plano_horizontal);
    normaliza_vec(cand_plano_horizontal);



    hallar_inv(&arr_elems,molecula,&num_total_elementos);
    printf("TIPO\tORDEN\tx\ty\tz\n");
    printf("NUMERO DE ELEMENTOS DE SIMETRIA:%d\n",num_total_elementos);
    for (i=0;i<num_total_elementos;i++){
        printf("%d\t%d\t%lf\t%lf\t%lf\n",arr_elems[i].ref,arr_elems[i].order,arr_elems[i].vec_dic[0],arr_elems[i].vec_dic[1],arr_elems[i].vec_dic[2]);
    }
    impr_elems_sim(arr_elems,num_total_elementos);
    printf("Candidato a lano horizontal: X:%lf\tY:%lf\tZ:%lf\n",cand_plano_horizontal[0],cand_plano_horizontal[1],cand_plano_horizontal[2]);

    free(arr_elems);
    arr_elems = NULL;
    return 0;
}


