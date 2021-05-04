/* PIEDOCK                                                            *
 * Authors: D. V. S. Ravikant                                         *
 * (C) 2011 D. V. S. Ravikant and Ron Elber                           */

#include <stdio.h>
#include <math.h>
#include "rmsd.h"

#ifndef EPSILON
#define EPSILON 0.00000000001
#endif
#ifndef PI
#define PI 3.14159265358
//#define PI 3.14159265358979323846
#endif
//#define abs(x) (x > 0) ? x : -x;

/* vector functions using c arrays */

void normalize(float a[3]){
	float  b;
	
	b = sqrt((float)(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]));
	a[0] /= b;
	a[1] /= b;
	a[2] /= b;
}

static float dot(float a[3], float b[3]){
	return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
}

static void cross(float a[3], float b[3], float c[3]){
	a[0] = b[1]*c[2] - b[2]*c[1];
	a[1] = b[2]*c[0] - b[0]*c[2];
	a[2] = b[0]*c[1] - b[1]*c[0];
}

static float dot(double a[3], float b[3]){
	return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
}

static void cross(float a[3], double b[3], double c[3]){
	a[0] = b[1]*c[2] - b[2]*c[1];
	a[1] = b[2]*c[0] - b[0]*c[2];
	a[2] = b[0]*c[1] - b[1]*c[0];
}

/*
 * setup_rotation() 
 *
 * given two lists of x,y,z coordinates, constructs
 * the correlation R matrix and the E value needed to calculate the
 * least-squares rotation matrix.
 */
void setup_rotation(float ref_xlist[][3], float mov_xlist[][3], int n_list, float mov_com[3], float mov_to_ref[3], double R[3][3], double* E0) {
	int i, j, n;
	float ref_com[3];
	
	/* calculate the centre of mass */
	for (i=0; i<3; i++) { 
		mov_com[i] = 0.0;
		ref_com[i] = 0.0;
	}
	  
	/*for (n=0; n<n_list; n++){
		printf("ref\t(%f,%f,%f)\tmov\t(%f,%f,%f)\n",ref_xlist[n][0],ref_xlist[n][1],ref_xlist[n][2],mov_xlist[n][0],mov_xlist[n][1],mov_xlist[n][2]);  
	}*/
    
	for (n=0; n<n_list; n++) 
	    for (i=0; i<3; i++){ 
			mov_com[i] += mov_xlist[n][i];
			ref_com[i] += ref_xlist[n][i];
		}
	    
	for (i=0; i<3; i++){
		mov_com[i] /= n_list;
		ref_com[i] /= n_list;
		mov_to_ref[i] = ref_com[i] - mov_com[i];
	}
	
	/* shift mov_xlist and ref_xlist to centre of mass */
	for (n=0; n<n_list; n++) 
		for (i=0; i<3; i++){ 
			mov_xlist[n][i] -= mov_com[i];
			ref_xlist[n][i] -= ref_com[i];
		}
	
	/* initialize */
	for (i=0; i<3; i++)
		for (j=0; j<3; j++) 
	    	R[i][j] = 0.0;
	*E0 = 0.0;

	for (n=0; n<n_list; n++){
		/* E0 = 1/2 * sum(over n): y(n)*y(n) + x(n)*x(n) */
		for (i=0; i<3; i++)
			*E0 +=  mov_xlist[n][i] * mov_xlist[n][i] + ref_xlist[n][i] * ref_xlist[n][i];
		    
		    /*
		     * correlation matrix R:   
		     *   R[i,j] = sum(over n): y(n,i) * x(n,j)  
		     *   where x(n) and y(n) are two vector sets   
		     */
		    for (i=0; i<3; i++){
		      for (j=0; j<3; j++)
		        R[i][j] += mov_xlist[n][i] * ref_xlist[n][j];
			}
	}
	  
	/*printf("R");
	for (i=0; i<3; i++){
		printf("\t"); 
		for (j=0; j<3; j++)
	    	printf("%f ",R[i][j]);
	  	printf("\n");
	}*/
	  
	*E0 *= 0.5;
}

#define ROTATE(a,i,j,k,l) { g = a[i][j]; \
                            h = a[k][l]; \
                            a[i][j] = g-s*(h+g*tau); \
                            a[k][l] = h+s*(g-h*tau); } //g+s*(h-g*tau);
 
	                       
/*   
 * jacobi3
 *
 * computes eigenval and eigen_vec of a real 3x3
 * symmetric matrix. On output, elements of a that are above 
 * the diagonal are destroyed. d[1..3] returns the 
 * eigenval of a. v[1..3][1..3] is a matrix whose 
 * columns contain, on output, the normalized eigen_vec of
 * a. n_rot returns the number of Jacobi rotations that were required.
 */
int jacobi3(float a[3][3], float d[3], float v[3][3], int* n_rot){
	int count, k, i, j;
	float tresh, theta, tau, t, sum, s, h, g, c, b[3], z[3];

	/*Initialize v to the identity matrix.*/
	for (i=0; i<3; i++) { 
		for (j=0; j<3; j++) 
      		v[i][j] = 0.0;
    	v[i][i] = 1.0;
	}

	/* Initialize b and d to the diagonal of a */
  	for (i=0; i<3; i++) 
    	b[i] = d[i] = a[i][i];

  	/* z will accumulate terms */
  	for (i=0; i<3; i++) 
    	z[i] = 0.0; 
  
  	*n_rot = 0;
  	
  	for (count=0; count<50; count++) {
		/*float av[3][3],vtav[3][3],vtv[3][3];
	
		for (i=0; i<2; i++)
      	for (j=i+1; j<3; j++)
      		a[j][i] = a[i][j];
      	
    	for(i = 0 ; i < 3; i++)
    		a[i][i] = d[i];
    
		for(i = 0 ; i < 3; i++){
			for(j = 0 ; j < 3 ; j++){
				av[i][j] = 0;
				vtv[i][j] = 0;
				for(k = 0 ; k < 3 ; k++){
					av[i][j] += a[i][k]*v[k][j];
					vtv[i][j] += v[k][i]*v[k][j];
				}
			}
		}
	
		printf("v =");
		for(i = 0 ; i < 3 ; i++){
			printf("\t");
			for(j = 0 ; j < 3 ; j++){
				printf("%f ",v[i][j]);
			}
			printf("\n");
	  	}
	  	printf("\n");
	  	
	  	printf("a =");
		for(i = 0 ; i < 3 ; i++){
			printf("\t");
			for(j = 0 ; j < 3 ; j++){
				printf("%f ",a[i][j]);
			}
			printf("\n");
	  	}
	  	printf("\n");
	  	
		printf("a*v =");
		for(i = 0 ; i < 3 ; i++){
			printf("\t");
			for(j = 0 ; j < 3 ; j++){
				printf("%f ",av[i][j]);
			}
			printf("\n");
	  	}
	  	printf("\n");
  	
	  	printf("vt*v =");
		for(i = 0 ; i < 3 ; i++){
			printf("\t");
			for(j = 0 ; j < 3 ; j++){
				printf("%f ",vtv[i][j]);
			}
			printf("\n");
	  	}
	  	printf("\n");
	  	
		for(i = 0 ; i < 3; i++){
			for(j = 0 ; j < 3 ; j++){
				vtav[i][j] = 0;
				for(k = 0 ; k < 3 ; k++){
					vtav[i][j] += v[k][i]*av[k][j];
				}
			}
		}
		 
		printf("%d vt*a*v =",*n_rot);
		for(i = 0 ; i < 3 ; i++){
			printf("\t");
			for(j = 0 ; j < 3 ; j++){
				printf("%f ",vtav[i][j]);
			}
			printf("\n");
	  	}
	  	printf("\n");
		*/
		
	    /* sum off-diagonal elements */
		sum = 0.0;
	    for(i=0; i<2; i++){
	    	for(j=i+1; j<3; j++)
	         	sum += fabs(a[i][j]);
	    }

	    /* if converged to machine underflow */
	    if(sum == 0)	return(1);
	
	    /* on 1st three sweeps... */
	    if(count < 3)
	      	tresh = sum * 0.2 / 9.0;    
	    else
	      	tresh = 0.0;      
	
		for (i=0; i<2; i++) {
	      	for (j=i+1; j<3; j++) {
	        	g = 100.0 * fabs(a[i][j]);
	
	        	/*  after four sweeps, skip the rotation if
	         	*   the off-diagonal element is small */
	        	if ( count > 3  &&  fabs(d[i])+g == fabs(d[i]) &&  fabs(d[j])+g == fabs(d[j]) ) {
	          		a[i][j] = 0.0;
	        	} else if (fabs(a[i][j]) > tresh) {
	          		h = d[j] - d[i];
	          
	          		if (fabs(h)+g == fabs(h)) {
	            		t = a[i][j] / h;
	          		} else {
	            		theta = 0.5 * h / (a[i][j]);
	            		t = 1.0 / ( fabs(theta) + (float)sqrt(1.0 + theta*theta) );
	            		if (theta < 0.0) 
	              		t = -t;
	          		}
	          
	          		c = 1.0 / (float) sqrt(1 + t*t);
	          		s = t * c;
	          		tau = s / (1.0 + c);
	          		h = t * a[i][j];

	          		z[i] -= h;
	          		z[j] += h;
	          		d[i] -= h;
	          		d[j] += h;
	
	          		a[i][j] = 0.0;
	
	          		for (k=0; k<=i-1; k++) 
	            		ROTATE(a, k, i, k, j)
	
	          		for (k=i+1; k<=j-1; k++) 
	            		ROTATE(a, i, k, k, j)
	
	          		for (k=j+1; k<3; k++) 
	            		ROTATE(a, i, k, j, k)
	
	          		for (k=0; k<3; k++) 
	            		ROTATE(v, i, k, j, k)
	
	          		++(*n_rot);
	        	}
	      	}
		}
	
	    for(i=0; i<3; i++){
	      	b[i] += z[i];
	      	d[i] = b[i];
	      	z[i] = 0.0;
		}
	}

	printf("Too many iterations in jacobi3; sum %f\n",sum);
  	return (0);
}  


/* 
 * diagonalize_symmetric 
 *
 * Diagonalize a 3x3 matrix & sort eigenval by size
 */
int diagonalize_symmetric(double matrix[3][3], float eigen_vec[3][3], float eigenval[3]) {
  	int n_rot, i, j, k;
  	float vec[3][3];
  	float val; 

  	float m[3][3];
  	for(int i = 0 ; i < 3; i++)
  		for(int j = 0 ; j < 3 ; j++)
  			m[i][j] = matrix[i][j];
  		  			
  	if (!jacobi3(m, eigenval, vec, &n_rot)) {
    	printf("convergence failed\n");
    	return (0);
  	}

  	/* sort solutions by eigenval */
  	for (i=0; i<3; i++) {
    	k = i;
    	val = eigenval[i];
    
    	for (j=i+1; j<3; j++)
      	if (eigenval[j] >= val) { 
        	k = j;
        	val = eigenval[k];
      	}
       
    	if (k != i) {
	      	eigenval[k] = eigenval[i];
	      	eigenval[i] = val;
	      	for (j=0; j<3; j++) {
	        	val = vec[i][j];
	        	vec[i][j] = vec[k][j];
	        	vec[k][j] = val;
	      	}
    	}
	}

  	/* transpose such that first index refers to solution index */
  	for (i=0; i<3; i++)
   		for (j=0; j<3; j++){
      		eigen_vec[i][j] = vec[i][j];
    	}

  	/*printf("\neigenvalues\t%f\t%f\t%f\n",eigenval[0],eigenval[1],eigenval[2]);
  	printf("vectors\t%f %f %f\n\t%f %f %f\n\t%f %f %f\n",eigen_vec[0][0] , eigen_vec[1][0], eigen_vec[2][0]
  		,eigen_vec[0][1] , eigen_vec[1][1], eigen_vec[2][1] ,eigen_vec[0][2] , eigen_vec[1][2], eigen_vec[2][2]);
  	printf("M\t%f %f %f\n\t%f %f %f\n\t%f %f %f\n",matrix[0][0] , matrix[0][1], matrix[0][2]
  			,matrix[1][0] , matrix[1][1], matrix[1][2] ,matrix[2][0] , matrix[2][1], matrix[2][2]);
  	/*printf("M * ev\t%f %f %f\n\t%f %f %f\n\t%f %f %f\n"
  			, dot(&matrix[0][0],&eigen_vec[0][0]), dot(&matrix[0][0],&eigen_vec[1][0]), dot(&matrix[0][0],&eigen_vec[2][0])
  			, dot(&matrix[1][0],&eigen_vec[0][0]), dot(&matrix[1][0],&eigen_vec[1][0]), dot(&matrix[1][0],&eigen_vec[2][0])
  			, dot(&matrix[2][0],&eigen_vec[0][0]), dot(&matrix[2][0],&eigen_vec[1][0]), dot(&matrix[2][0],&eigen_vec[2][0]) );
  	*/
  	return (1);
}


/*
 * calculate_rotation_matrix() 
 *
 *   calculates the rotation matrix U and the
 * rmsd from the R matrix and E0:
 */
int calculate_rotation_matrix(double R[3][3], float U[3][3], int *rank, double E0, float* residual) {
  	int i, j, k;
  	double Rt[3][3], RtR[3][3],UtU[3][3];
  	float b[3][3], a[3][3], eigenval[3], v[3][3], vt[3][3];
  	float vtdv[3][3], vtd[3][3], Rvtdv[3][3], vtdvRt[3][3];
  	float vv[3];
  	float sigma;

  	/* build Rt, transpose of R  */
  	for (i=0; i<3; i++)
    	for (j=0; j<3; j++)
      		Rt[i][j] = R[j][i];

  	/* make symmetric RtR = Rt X R */
  	for (i=0; i<3; i++) 
    	for (j=0; j<3; j++) {
      		RtR[i][j] = 0.0;
      	for (k = 0; k<3; k++)
        	RtR[i][j] += Rt[i][k] * R[k][j];
    	}
  	
  	double m[3][3];
  	for(i = 0; i < 3 ;i++)
  		for(j = 0; j < 3 ; j++)
  			m[i][j] = RtR[i][j];
  		
  	if(!diagonalize_symmetric(m, a, eigenval))	return(0);
    
  	/* Let's force the eigen vectors of RtR into right-handed system.*/
  	cross(&a[2][0], &a[0][0], &a[1][0]);

  	/*printf("R\t%f %f %f\n\t%f %f %f\n\t%f %f %f\n",R[0][0] , R[0][1], R[0][2]
  			,R[1][0] , R[1][1], R[1][2] ,R[2][0] , R[2][1], R[2][2]);
  	printf("a\t%f %f %f\n\t%f %f %f\n\t%f %f %f\n",a[0][0] , a[1][0], a[2][0]
  			,a[0][1] , a[1][1], a[2][1] ,a[0][2] , a[1][2], a[2][2]);
  	printf("R * a\t%f %f %f\n\t%f %f %f\n\t%f %f %f\n" , dot(&R[0][0],&a[0][0]), dot(&R[0][0],&a[1][0]), dot(&R[0][0],&a[2][0])
  			, dot(&R[1][0],&a[0][0]), dot(&R[1][0],&a[1][0]), dot(&R[1][0],&a[2][0]) , dot(&R[2][0],&a[0][0]), dot(&R[2][0],&a[1][0]), dot(&R[2][0],&a[2][0]) );
  	printf("RtR * a\t%f %f %f\n\t%f %f %f\n\t%f %f %f\n" , dot(&RtR[0][0],&a[0][0]), dot(&RtR[0][0],&a[1][0]), dot(&RtR[0][0],&a[2][0])
  			, dot(&RtR[1][0],&a[0][0]), dot(&RtR[1][0],&a[1][0]), dot(&RtR[1][0],&a[2][0]) , dot(&RtR[2][0],&a[0][0]), dot(&RtR[2][0],&a[1][0]), dot(&RtR[2][0],&a[2][0]) );
  	*/

  	for (i=0; i<3; i++) 
    	for (j=0; j<3; j++){
    		v[i][j] = vt[j][i] = a[i][j];
    	} 
    
  	/*printf("v");
  	for(i = 0 ; i < 3 ; i++){
  		printf("\t");
		for(int j = 0 ; j < 3 ; j++)
			printf("%f ",v[i][j]);
		printf("\n");
  	}
  	printf("\n");
  	*/
	  
  	for (i=0; i<3; i++) 
    	for (j=0; j<3; j++) 
      		b[j][i] = dot(&R[i][0], &a[j][0]);
  			
  	*rank = 0;
  	for (i=0; i<3; i++){
  		if(fabs(eigenval[i]) > EPSILON)
  			(*rank)++;
  	}
  
  	if(*rank == 3){
	  	/* 
	   	* First check if the rotational matrices generated from the 
	   	* orthogonal eigenvectors are in a right-handed or left-handed
	   	* co-ordinate system - given by sigma. Sigma is needed to
	   	* resolve this ambiguity in calculating the RMSD.
	   	*/
	  	cross(vv, &b[0][0], &b[1][0]);
	  
	  	if (dot(vv, &b[2][0]) < 0.0)
	    	sigma = -1.0;
	  	else 
	    	sigma = 1.0;
  	} else
  		sigma = 1.0;
	  
	if(*rank == 3){
    	/*printf("%d b\t%f %f %f\n\t%f %f %f\n\t%f %f %f\n", *rank, b[0][0] , b[1][0], b[2][0]
  			,b[0][1] , b[1][1], b[2][1] ,b[0][2] , b[1][2], b[2][2]);
 
 	  	printf("\t%f\t%f\t%f\n", b[0][0]*b[0][0] + b[0][1]*b[0][1] + b[0][2]*b[0][2],
 	  			 b[1][0]*b[1][0] + b[1][1]*b[1][1] + b[1][2]*b[1][2], b[2][0]*b[2][0] + b[2][1]*b[2][1] + b[2][2]*b[2][2]);
  	  	*/
  	  
	  	/* calc optimal rotation matrix U that minimises residual */
	  	for (i=0;i<3; i++)
	    	for (j=0; j<3; j++)
	      		U[i][j] = 0.0;
	   
	  	float d[3];
	  	//compute u = Rvtdv
	  	for (k=0; k<3; k++)
	  		if(fabs(eigenval[k]) > EPSILON){
	    		float f = 1.0/sqrt(fabs(eigenval[k]));
	      		if(k == 2)
	      			f *= sigma;
	      		d[k] = f;
	      		for (i=0;i<3; i++)
	    			for (j=0; j<3; j++)
	    		  		/* need to apply U to the second set of coordinates,
				   		* U calculated here is for transforming the first
				   		* set of coordinates 
				  		*/
	  		      		U[j][i] += b[k][i] * a[k][j] *f;
	    	}
	  
	  	/*printf("U");
	  	for(i = 0 ; i < 3 ; i++){
	  		printf("\t");
			for(int j = 0 ; j < 3 ; j++)
				printf("%f ",U[i][j]);
			printf("\n");
	  	}
	 
	 	for(i = 0 ; i < 3; i++){
			for(j = 0 ; j < 3 ; j++){
				UtU[i][j] = 0;
				for(k = 0 ; k < 3 ; k++){
					UtU[i][j] += U[k][i]*U[k][j];
				}
			}
	  	}
	  	printf("UtU");
	  	for(i = 0 ; i < 3 ; i++){
	  		printf("\t");
			for(int j = 0 ; j < 3 ; j++)
				printf("%f ",UtU[i][j]);
			printf("\n");
	  	}
	  	printf("\n");
	  
	  	printf("R");
	  	for(i = 0 ; i < 3 ; i++){
	  		printf("\t");
			for(int j = 0 ; j < 3 ; j++)
				printf("%f ",R[i][j]);
			printf("\n");
	  	}
	  	printf("RtR");
	  	for(i = 0 ; i < 3 ; i++){
	  		printf("\t");
			for(j = 0 ; j < 3 ; j++){
				RtR[i][j] = 0;
				for(k = 0 ; k < 3; k++)
					RtR[i][j] += R[k][i]*R[k][j];
				printf("%f ",RtR[i][j]);
			}
			printf("\n");
	  	}
	  
	  	for(i=0;i<3;i++)
	  		for(j=0;j<3;j++)
	  			vtd[i][j] = vt[i][j]*d[j];
	  	for(i=0;i<3;i++)
	  		for(j=0;j<3;j++)
	   			vtdv[i][j] = dot(&vtd[i][0],&vt[j][0]);
	  
		float vtdvRtRvtdv[3][3];
	  	for(i=0;i<3;i++)
	  		for(j=0;j<3;j++){
	  			Rvtdv[i][j] = 0;
	  			vtdvRt[i][j] = 0;
	  			for(k=0;k<3;k++){
	  				Rvtdv[i][j] += R[i][k]*vtdv[k][j];
	  				vtdvRt[i][j] += vtdv[i][k]*Rt[k][j];
	  			}
	  		}
	  
		for(i=0;i<3;i++)
	  		for(j=0;j<3;j++){
	  			vtdvRtRvtdv[i][j] = 0;
	  			for(k=0;k<3;k++)
	  				vtdvRtRvtdv[i][j] += vtdvRt[i][k]*Rvtdv[k][j];
	  		}
	  
		printf("Rvtdv");
	  	for(i = 0 ; i < 3 ; i++){
	  		printf("\t");
			for(int j = 0 ; j < 3 ; j++)
				printf("%f ",Rvtdv[i][j]);
			printf("\n");
		}
	  
	  	printf("vtdvRtRvtdv");
	  	for(i = 0 ; i < 3 ; i++){
	  		printf("\t");
			for(int j = 0 ; j < 3 ; j++)
				printf("%f ",vtdvRtRvtdv[i][j]);
			printf("\n");
	  	}*/
	}
  	*residual = E0 - (float) sqrt(fabs(eigenval[0])) - (float) sqrt(fabs(eigenval[1])) - sigma * (float) sqrt(fabs(eigenval[2]));
  	//printf("\nregular %f\t%f %f %f %f %f\n",*residual , E0 , eigenval[0] , eigenval[1] ,eigenval[2], sigma);
  	
  	return (1);
}


void calculate_rotation_rmsd(float ref_xlist[][3], float mov_xlist[][3], int n_list, float mov_com[3], float mov_to_ref[3], int *rank,
	float U[3][3], float* rmsd) {
  	double Eo;
  	float residual;
  	double R[3][3];
  
  	setup_rotation(ref_xlist, mov_xlist, n_list, mov_com, mov_to_ref, R, &Eo);
  	int ret = calculate_rotation_matrix(R, U, rank, Eo, &residual);

  	if(ret == 0)
  		printf("no rotation #points %d\n",n_list);
  	/* check if the computed transformation is correct */
  	/*if(*rank ==3){
  		float rmsdc = 0;
  		for(int n = 0 ; n <n_list; n++){
  			float v[3];
  			for(int j = 0 ; j < 3 ; j++){
  				v[j] = 0;
  				for(int k = 0 ; k < 3 ; k++)
  					v[j] += U[j][k]*mov_xlist[n][k];
	  			rmsdc += (ref_xlist[n][j] - v[j])*(ref_xlist[n][j] - v[j]); 
  			}
  		}
  		printf("check rmsd %f\n",sqrt(rmsdc/n_list));
  	}*/
  
  	residual = fabs(residual); /* avoids the awkward case of -0.0 */
  	*rmsd = sqrt( fabs((float) (residual)*2.0/((float)n_list)) );
  	//float rmsd1;
  	//fast_rmsd(ref_xlist,mov_xlist,n_list,&rmsd1); 
}
 

/*
 * Fast calculation of rmsd w/o calculating a rotation matrix.
 *
 * Fast rmsd calculation by the method of  
 * Kabsch 1978, where the required eigenvalues are found by an 
 * analytical, rather than iterative, method to save time. 
 * The cubic factorization used to accomplish this only produces 
 * stable eigenvalues for the transpose(R]*R matrix of a typical 
 * protein after the whole matrix has been normalized.
 */
void cubicsolve_rmsd(float ref_xlist[][3], float mov_xlist[][3], int n_list, float* rmsd){ 
  	float mov_com[3];
  	float mov_to_ref[3];
  	double R[3][3], Eo;
    
  	setup_rotation(ref_xlist, mov_xlist, n_list, mov_com, mov_to_ref, R, &Eo);
  
  	cubicsolve_rmsd_fromRE(R, Eo, n_list, rmsd);
}

void cubicsolve_rmsd_fromRE(double R[3][3], double Eo, int n_list, float *rmsd){
    double d0,d1,d2,e0,e1,f0;
  	float omega;
	float v[3];
	float residual;
  	
  	/* cubic roots */
  	double r1,r2,r3;
  	double sr;
  
  	/* 
   	* check if the determinant is greater than 0 to
   	* see if R produces a right-handed or left-handed
   	* co-ordinate system.
   	*/
  	cross(v, &R[1][0], &R[2][0]);
  	if (dot(&R[0][0], v) > 0.0)
    	omega = 1.0;
  	else
    	omega = -1.0;

  	/*
   	* get elements we need from tran(R) x R 
   	*          matrix = d0 e0 f0
   	*                      d1 e1
   	*                         d2
   	*/
   
  	d0 =  R[0][0]*R[0][0] + R[1][0]*R[1][0] + R[2][0]*R[2][0];
  	d1 = (R[0][1]*R[0][1] + R[1][1]*R[1][1] + R[2][1]*R[2][1]);
  	d2 = (R[0][2]*R[0][2] + R[1][2]*R[1][2] + R[2][2]*R[2][2]);

  	e0 = (R[0][0]*R[0][1] + R[1][0]*R[1][1] + R[2][0]*R[2][1]);
  	e1 = (R[0][1]*R[0][2] + R[1][1]*R[1][2] + R[2][1]*R[2][2]);

  	f0 = (R[0][0]*R[0][2] + R[1][0]*R[1][2] + R[2][0]*R[2][2]);

  	/* cubic roots */
  	{
    	double B, C, D, q, q3, r, theta;
    	/*
     	* solving for eigenvalues as det(A-I*lambda) = 0
     	* yeilds the values below corresponding to:
     	* lambda**3 + B*lambda**2 + C*lambda + D = 0
     	*/
    	B = 0-(d0 + d1 + d2);
    	C = d0*d1 + d0*d2 + d1*d2 - (e0*e0 + f0*f0 + e1*e1);
    	D = (e0*e0*d2 + e1*e1*d0 + f0*f0*d1) - (d0*d1*d2 + 2*e0*f0*e1);

    	q = (B*B - 3.0*C) / 9.0;
    	r = (2.0*B*B*B - 9.0*B*C + 27.0*D) / 54.0;
    	// reduced to the form y^3 -3qy + 2r = 0, y = lambda + B/3
    	
    	//printf("x^3 + %f x^2 + %f x + %f = 0\n",B,C,D);
    	//printf("y^3 -3qy + 2r = 0\t q=%f\t r= %f\n",q,r);
    	
    	if( q == 0 && r == 0){
    		r1 = r2 = r3 = 0;
    	} else if( q != 0 && r == 0){
    		r1 = 0;
    		r2 = sqrt(3*q);
    		r3 = 0-r2;
    	} else if(r != 0 && q == 0){
    		r1 = r2 = r3 = cbrt(-2*r);
    	} else { //if(r != 0 && q != 0)
    		// use the information that the eigen values are real to reduce
    		q3 = q*q*q;
    		double costheta = r/sqrt(q3); 
    		if(costheta <= -1){
    			printf("WARNING : rmsd - solve cubic - had costheta = %f\n",costheta);
    			costheta = -1;
    		} else if(costheta >= 1){
    			printf("WARNING : rmsd - solve cubic - had costheta = %f\n",costheta);
    			costheta = 1;
    		}
    		theta = acos(costheta);
    		
    		r1 = r2 = r3 = -2.0*sqrt(q);
    		r1 *= cos(theta/3.0);
    		r2 *= cos((theta + 2.0*PI) / 3.0);
    		r3 *= cos((theta - 2.0*PI) / 3.0);
    	}
    	r1 -= B / 3.0;
    	r2 -= B / 3.0;
    	r3 -= B / 3.0;
  	}

	if(r1<0 || r2<0 || r3<0){
		printf("WARNING : rmsd - solve cubic - have negative roots %f,%f,%f\n",r1,r2,r3);
		printf("%f\t%f\t%f\n%f\t%f\t%f\n%f\t%f\t%f\n",R[0][0],R[0][1],R[0][2],R[1][0],R[1][1],R[1][2],R[2][0],R[2][1],R[2][2]);
		printf("%f\t%f\t%f\n%f\t%f\t%f\n%f\t%f\t%f\n",d0,e0,f0,e0,d1,e1,f0,e1,d2);
	}
		
  	/* expect positive eigenvalues */
  	r1 = fabs(r1);
  	r2 = fabs(r2);
  	r3 = fabs(r3);
  	
  	/* sort roots in increasing order */
  	if (r2<r1) {
    	sr = r2;
    	r2 = r1;
    	r1 = sr;
  	}
  	if (r3<r1){
  		sr = r1;
  		r1 = r3;
  		r3 = r2;
  		r2 = sr;
  	} else if (r3<r2){
  		sr = r2;
  		r2 = r3;
  		r3 = sr;
  	}

  	residual = Eo - sqrt(r3) - sqrt(r2) - omega*sqrt(r1);
  	residual = fabs(residual);
  	*rmsd = sqrt( (float) residual*2.0 / ((float) n_list) );
  	//printf("\nfast %f\t%f %f %f %f %f %f\n",residual , Eo , r3,r2,r1, omega, *rmsd);
  	//printf("r0 %f\tr1 %f\tr2 %f\tresidual %f\trmsd %f\n",, residual, *rmsd); 
}
