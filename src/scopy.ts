/*
    Copies a vector, x, to a vector, y.

    Vector y should already be the same size vector as x

    the original copy allowed different incx and incy

    Pawel Spychala 2020
*/
export function scopy(sx: number[],sy: number[]) {
    for (let i = 1; i < sx.length; i++) {
        sy[i] = sx[i]
    }
}