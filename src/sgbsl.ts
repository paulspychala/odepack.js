/*
    sgbsl solves the real band system
    a * x = b  or  trans(a) * x = b
    using the factors computed by sgbco or sgbfa.
    on entry
       abd     real(lda, n)
               the output from sgbco or sgbfa.
       lda     integer
               the leading dimension of the array  abd .
       n       integer
               the order of the original matrix.
       ml      integer
               number of diagonals below the main diagonal.
       mu      integer
               number of diagonals above the main diagonal.
       ipvt    integer(n)
               the pivot vector from sgbco or sgbfa.
       b       real(n)
               the right hand side vector.
       job     integer
               = 0         to solve  a*x = b ,
               = nonzero   to solve  trans(a)*x = b , where
                           trans(a)  is the transpose.
    on return
       b       the solution vector  x .
    error condition
       a division by zero will occur if the input factor contains a
       zero on the diagonal.  technically this indicates singularity
       but it is often caused by improper arguments or improper
       setting of lda .  it will not occur if the subroutines are
       called correctly and if sgbco has set rcond .gt. 0.0
       or sgbfa has set info .eq. 0 .
    to compute  inverse(a) * c  where  c  is a matrix
    with  p  columns
          call sgbco(abd,lda,n,ml,mu,ipvt,rcond,z)
          if (rcond is too small) go to ...
          do 10 j = 1, p
             call sgbsl(abd,lda,n,ml,mu,ipvt,c(1,j),0)
       10 continue
    linpack. this version dated 07/14/77 .
    cleve moler, university of new mexico, argonne national labs.
    translated to typescript by Pawel Spychala 2020
*/
import {sdot} from "./sdot";
import {saxpy} from "./saxpy";

export function sgbsl(abd: number[][],lda, n: number, ml: number, mu: number, ipvt: number[], b: number[], job: number) {
    let t: number
    let k: number
    let kb: number
    let l: number
    let la: number
    let lb: number
    let lm: number
    let m = mu + ml + 1
    let nm1 = n - 1
    if (job == 0) { 
        // job = 0 , solve  a * x = b
        // first solve l*y = b
        if (ml != 0) {
            if (nm1 >= 1) {
                for (k = 1; k <= nm1; k++) {
                    lm = Math.min(ml,n-k)
                    l = ipvt[k]
                    t = b[l]
                    if (l != k) {
                        b[l] = b[k]
                        b[k] = t
                    }
                    saxpy(lm,t,abd[k],m,1,b, k,1)
                }
            }
        }
        // now solve  u*x = y
        for (kb = 1; kb <= n; kb++) {
            k = n + 1 - kb
            b[k] = b[k]/abd[m][k]
            lm = Math.min(k,m) - 1
            la = m - lm
            lb = k - lm
            t = -b[k]
            saxpy(lm,t,abd[k], la-1,1,b, lb-1,1)
        }
    } else { // job != 0

        // job = nonzero, solve  trans(a) * x = b
        // first solve  trans(u)*y = b
        for (k = 1; k <= n; k++) {
            lm = Math.min(k,m) - 1
            la = m - lm
            lb = k - lm
            t = sdot(lm,abd[k],la-1,1,b,lb-1,1)
            b[k] = (b[k] - t)/abd[m][k]
        }
        // now solve trans(l)*x = y
        if (ml != 0) {
            if (nm1 >= 1) {
                for (kb = 1; kb <= nm1; kb++) {
                    k = n - kb
                    lm = Math.min(ml,n-k)
                    b[k] = b[k] + sdot(lm,abd[k],m,1,b,k,1)
                    l = ipvt[k]
                    if (l != k) {
                        t = b[l]
                        b[l] = b[k]
                        b[k] = t
                    }
                }
            }
        }
    }
}