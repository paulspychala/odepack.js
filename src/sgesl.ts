/*
    sgesl solves the real system
    a * x = b  or  trans(a) * x = b
    using the factors computed by sgeco or sgefa.

    on entry

    a       real(lda, n)
            the output from sgeco or sgefa.

    lda     integer
            the leading dimension of the array  a .

    n       integer
            the order of the matrix  a .

    ipvt    integer(n)
            the pivot vector from sgeco or sgefa.

    b       real(n)
            the right hand side vector.

    job     integer
            = 0         to solve  a*x = b ,
            = nonzero   to solve  trans(a)*x = b  where
                        trans(a)  is the transpose.

    on return

    b       the solution vector  x .

    error condition

    a division by zero will occur if the input factor contains a
    zero on the diagonal.  technically this indicates singularity
    but it is often caused by improper arguments or improper
    setting of lda .  it will not occur if the subroutines are
    called correctly and if sgeco has set rcond .gt. 0.0
    or sgefa has set info .eq. 0 .

    to compute  inverse(a) * c  where  c  is a matrix
    with  p  columns
        call sgeco(a,lda,n,ipvt,rcond,z)
        if (rcond is too small) go to ...
        do 10 j = 1, p
            call sgesl(a,lda,n,ipvt,c(1,j),0)
    10 continue

    linpack. this version dated 07/14/77 .
    cleve moler, university of new mexico, argonne national labs.
    translated to Typrscript by Pawel Spychala 2020
*/
import {sdot} from "./sdot";
import {saxpy} from "./saxpy";

export function sgesl(a: number[][], n: number, ipvt: number[], b: number[], job: number)
{
    let nm1: number
    let k: number
    let t: number
    let l: number
    let kb: number

	nm1 = n - 1;
	if (job == 0) {

        // job = 0 , solve  a * x = b
        // first solve  l*y = b

        for (k = 1; k <= n; k++) {
            t = sdot(k-1, a[k], 0, 1, b, 0, 1)
            // t = sdot(k-1,a(1,k),1,b(1),1)
            b[k] = (b[k] - t) / a[k][k]
        }

        // now solve  u*x = y

        if (nm1 >= 1) {
            for (kb = 1; kb <= nm1; kb++) {
                k = n -kb
                b[k] = b[k] + sdot(n-k, a[k], k, 1, b, k, 1)
                // b(k) = b(k) + sdot(n-k,a(k+1,k),1,b(k+1),1)
                l = ipvt[k]
                if (l != k) {
                    t = b[l]
                    b[l] = b[k]
                    b[k] = t
                }
            }
        }

        

        return;
    }
    
    // job = nonzero, solve  trans(a) * x = b
    // first solve  trans(u)*y = b

    if(nm1 >= 1) {
        for (k = 1; k <= nm1; k++) {
            l = ipvt[k]
            t = b[l]
            if (l != k) {
                b[l] = b[k]
                b[k] = t
            }
            saxpy(n-k, t, a[k], k, 1, b, k, 1)
        }
    }

    // now solve trans(l)*x = y

    for (kb = 1; kb <= n; kb++) {
        k = n + 1 - kb
        b[k] = b[k]/a[k][k]
        t = -b[k]
        saxpy(k - 1, t, a[k], 0, 1, b, 0, 1);
    }
}