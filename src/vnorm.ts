/*
-----------------------------------------------------------------------
 this function routine computes the weighted root-mean-square norm
 of the vector of length n contained in the array v, with weights
 contained in the array w of length n..
   vnorm = sqrt( (1/n) * sum( v(i)*w(i) )**2 )

    translated to Typescript by Pawel Spychala 2020
-----------------------------------------------------------------------
*/
export function vnorm (n: number, v: number[], w: number[]) {
    let sum = 0
    for (let i = 1; i <= n; i++) {
        sum += Math.pow((v[i] * w[i]),2)
    }
    return Math.sqrt(sum/n)
}