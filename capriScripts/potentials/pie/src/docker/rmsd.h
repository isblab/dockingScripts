/* PIEDOCK                                                            *
 * Authors: D. V. S. Ravikant                                         *
 * (C) 2011 D. V. S. Ravikant and Ron Elber                           */

#ifndef _rmsd_h_
#define _rmsd_h_

/*
 * calculate_rotation_rmsd()
 *
 *   given two lists of x,y,z coordinates, constructs
 *    - mov_com: the centre of mass of the mov list
 *    - mov_to_ref: vector between the com of mov and ref
 *    - U: the rotation matrix for least-squares
 *    - rmsd: measures similarity between the vectors
 * 
 * to transform a mov point to the ref frame translate by mov_to_ref
 * and rotate using U
 *           for (i=0; i<3; i++)
 *           {
 *             rotated_v[i] = 0.0; 
 *             for (j=0; j<3; j++)
 *               rotated_v[i] += U[i][j] * v[j];
 *           }
 */
void calculate_rotation_rmsd(float ref_xlist[][3], float mov_xlist[][3], int n_list, float mov_com[3], float mov_to_ref[3],
                             int *rank, float U[3][3], float* rmsd);

/*
 *
 * Fast calculation of rmsd w/o calculating a rotation matrix,
 * adapted from the BTK by Chris Saunders 11/2002.
 * code fixed by Ravikant 
 */
void cubicsolve_rmsd(float ref_xlist[][3], float mov_xlist[][3], int n_list, float* rmsd); 

void cubicsolve_rmsd_fromRE(double R[3][3], double Eo, int, float *rmsd);

#endif
