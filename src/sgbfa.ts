/*
 sgbfa factors a real band matrix by elimination.
 sgbfa is usually called by sgbco, but it can be called
 directly with a saving in time if  rcond  is not needed.

 on entry
    abd     real(lda, n)
            contains the matrix in band storage.  the columns
            of the matrix are stored in the columns of  abd  and
            the diagonals of the matrix are stored in rows
            ml+1 through 2*ml+mu+1 of  abd .
            see the comments below for details.
    lda     integer
            the leading dimension of the array  abd .
            lda must be .ge. 2*ml + mu + 1 .
    n       integer
            the order of the original matrix.
    ml      integer
            number of diagonals below the main diagonal.
            0 <= ml < n .
    mu      integer
            number of diagonals above the main diagonal.
            0 <= mu < n .
            more efficient if  ml <= mu .

 on return
    abd     an upper triangular matrix in band storage and
            the multipliers which were used to obtain it.
            the factorization can be written  a = l*u  where
            l  is a product of permutation and unit lower
            triangular matrices and  u  is upper triangular.
    ipvt    integer(n)
            an integer vector of pivot indices.
    info    integer
            = 0  normal value.
            = k  if  u(k,k) == 0.0 .  this is not an error
                 condition for this subroutine, but it does
                 indicate that sgbsl will divide by zero if
                 called.  use  rcond  in sgbco for a reliable
                 indication of singularity.
                 
 band storage
       if  a  is a band matrix, the following program segment
       will set up the input.
               ml = (band width below the diagonal)
               mu = (band width above the diagonal)
               m = ml + mu + 1
               do 20 j = 1, n
                  i1 = max0(1, j-mu)
                  i2 = min0(n, j+ml)
                  do 10 i = i1, i2
                     k = i - j + m
                     abd(k,j) = a(i,j)
            10    continue
            20 continue
       this uses rows  ml+1  through  2*ml+mu+1  of  abd .
       in addition, the first  ml  rows in  abd  are used for
       elements generated during the triangularization.
       the total number of rows needed in  abd  is  2*ml+mu+1 .
       the  ml+mu by ml+mu  upper left triangle and the
       ml by ml  lower right triangle are not referenced.

 linpack. this version dated 07/14/77 .
 cleve moler, university of new mexico, argonne national labs.
 translated to typescript by Pawel Spychala 2020
*/
import {isamax} from "./isamax";
import {sscal} from "./sscal";
import {saxpy} from "./saxpy";

export function sgbfa(abd: number[][],lda: number,n: number,ml: number,mu: number, ipvt: number[],info: number) : number {
    let t: number
    let i: number
    // let isamax: number
    let i0: number
    let j: number
    let ju: number
    let jz: number
    let j0: number
    let j1: number
    let k: number
    let kp1: number
    let l: number
    let lm: number
    let m: number
    let mm: number
    let nm1: number

    m = ml + mu + 1
    info = 0

    // zero initial fill-in columns
    j0 = mu + 2
    j1 = Math.min(n,m) - 1
    if (j1 >= j0) {
        for (jz = j0; jz <= j1; jz++) {
            i0 = m + 1 - jz
            for (i = i0; i <= ml; i++) {
                abd[i][jz] = 0
            }
        }
    }
    jz = j1
    ju = 0

    // gaussian elimination with partial pivoting
    nm1 = n - 1
    if (nm1 >= 1) {
        for (k = 1; k <= nm1; k++) {
            kp1 = k + 1
            // zero next fill-in column
            jz = jz + 1
            if (jz <= n) {
                if (ml >= 1) {
                    for (i = 1; i <= ml; i++) {
                        abd[i][jz] = 0
                    }
                }
            }
            // find l = pivot index
            lm = Math.min(ml,n-k)
            // l = isamax(lm+1,abd(m,k),1) + m - 1
            l = isamax(lm+1,abd[k], m-1, 1) + m - 1
            ipvt[k] = l + k - m

            // zero pivot implies this column already triangularized
            if (abd[l][k] != 0) {
                
                // interchange if necessary
                if (l != m) {
                    t = abd[l][k]
                    abd[l][k] = abd[m][k]
                    abd[m][k] = t
                }

                // compute multipliers
                t = -1/abd[m][k]
                // sscal(lm,t,abd(m+1,k),1)
                sscal(lm,t,abd[k],m,1)

                // row elimination with column indexing
                ju = Math.min(Math.max(ju,mu+ipvt[k]),n)
                mm = m
                if (ju >= kp1) {
                    for (j = kp1; j <= ju; j++) {
                        l = l - 1
                        mm = mm - 1
                        t = abd[j][l]
                        if (l != mm) {
                            abd[j][l] = abd[j][mm]
                            abd[j][mm] = t
                        }
                        saxpy(lm,t,abd[k],m,1,abd[j],mm,1)
                    }
                }

            } else { //  abd[k][l] == 0
                info = k
            }
        }
    }
    ipvt[n] = n
    if (abd[n][m] == 0) info = n
    return info
}