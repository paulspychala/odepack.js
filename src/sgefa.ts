
/*
 sgefa factors a real matrix by gaussian elimination.
 sgefa is usually called by sgeco, but it can be called
 directly with a saving in time if  rcond  is not needed.
 (time for sgeco) = (1 + 9/n)*(time for sgefa) .

 on entry
    a       real(lda, n)
            the matrix to be factored.
    lda     integer
            the leading dimension of the array  a .
    n       integer
            the order of the matrix  a 
            .
 on return
    a       an upper triangular matrix and the multipliers
            which were used to obtain it.
            the factorization can be written  a = l*u  where
            l  is a product of permutation and unit lower
            triangular matrices and  u  is upper triangular.
    ipvt    integer(n)
            an integer vector of pivot indices.
    info    integer
            = 0  normal value.
            = k  if  u(k,k) == 0.0 .  this is not an error
                 condition for this subroutine, but it does
                 indicate that sgesl or sgedi will divide by zero
                 if called.  use  rcond  in sgeco for a reliable
                 indication of singularity.

 linpack. this version dated 07/14/77 .
 cleve moler, university of new mexico, argonne national labs.
 translation to Typescript by Pawel Spychala 2020
*/
import {saxpy} from "./saxpy";
import {isamax} from "./isamax";
import {sscal} from "./sscal";

export function sgefa(a: number[][], n: number, ipvt: number[], info: number) : number
{
    let k: number
    let i: number
    let t: number
    let nm1: number
    let l: number

    // gaussian elimination with partial pivoting

    info = 0;
    nm1 = n - 1
    if (nm1 >= 1) {
        for (k = 1; k <= nm1; k++) {

            // find l = pivot index
            // l = isamax(n-k+1,a(k,k),1) + k - 1
            l = isamax(n - k + 1, a[k], k - 1, 1) + k - 1;
            ipvt[k] = l;

            // zero pivot implies this column already triangularized
            if (a[k][l] == 0.) {
                info = k;
                continue;
            }

            // interchange if necessary
            if (l != k) {
                t = a[k][l];
                a[k][l] = a[k][k];
                a[k][k] = t;
            }

            // compute multipliers
            t = -1. / a[k][k];
            // sscal2(n -k, t, a, k, k, 1);
            sscal(n - k, t, a[k], k, 1);

            // row elimination with column indexing
            for (i = k + 1; i <= n; i++) {
                t = a[i][l];
                if (l != k) {
                    a[i][l] = a[i][k];
                    a[i][k] = t;
                }
                saxpy(n - k, t, a[k], k, 1, a[i], k, 1);
            }
        }
    }

	ipvt[n] = n;
	if (a[n][n] == 0) info = n;
	return info
}