/*
-----------------------------------------------------------------------
 this subroutine sets the error weight vector ewt according to
     ewt(i) = rtol(i)*abs(ycur(i)) + atol(i),  i = 1,...,n,
 with the subscript on rtol and/or atol possibly replaced by 1 above,
 depending on the value of itol.

    translated to Typescript by Pawel Spychala 2020
-----------------------------------------------------------------------
*/
export function ewset(n: number, itol: number, rtol: number[], atol: number[], ycur: number[], ewt: number[]) {
    var i: number

    switch (itol) {
    case 1:
        for (i = 1; i <= n; i++)
            ewt[i] = rtol[1] * Math.abs(ycur[i]) + atol[1];
        break;
    case 2:
        for (i = 1; i <= n; i++)
            ewt[i] = rtol[1] * Math.abs(ycur[i]) + atol[i];
        break;
    case 3:
        for (i = 1; i <= n; i++)
            ewt[i] = rtol[i] * Math.abs(ycur[i]) + atol[1];
        break;
    case 4:
        for (i = 1; i <= n; i++)
            ewt[i] = rtol[i] * Math.abs(ycur[i]) + atol[i];
        break;
    }
}