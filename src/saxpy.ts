/*
    constant times a vector plus a vector.
    uses unrolled loop for increments equal to one.
    jack dongarra, linpack, 6/17/77.
    translated to Typescript by Pawel Spychala 2020
*/
export function saxpy(n: number, sa: number, sx: number[], sxoffset: number, incx, sy: number[], syoffset: number, incy: number)
{
    let ix: number
    let iy: number
    let i: number
    let m: number

	if (n < 0 || sa == 0) return;

    // code for unequal increments or equal increments not equal to 1

	if (incx != incy || incx < 1) {
		ix = 1;
		iy = 1;
		if (incx < 0) ix = (-n + 1) * incx + 1;
		if (incy < 0) iy = (-n + 1) * incy + 1;
		for (i = 1; i <= n; i++) {
			sy[iy+syoffset] = sy[iy+syoffset] + sa * sx[ix+sxoffset];
			ix = ix + incx;
			iy = iy + incy;
		}
		return;
	}
    // code for both increments equal to 1

    // clean-up loop

	if (incx == 1) {
		m = n % 4;
		if (m != 0) {
			for (i = 1; i <= m; i++) {
                sy[i+syoffset] = sy[i+syoffset] + sa * sx[i+sxoffset];
            }
			if (n < 4) return;
		}
		for (i = m + 1; i <= n; i = i + 4) {
			sy[i+syoffset] = sy[i+syoffset] + sa * sx[i+sxoffset];
			sy[i+syoffset + 1] = sy[i+syoffset + 1] + sa * sx[i+sxoffset + 1];
			sy[i+syoffset + 2] = sy[i+syoffset + 2] + sa * sx[i+sxoffset + 2];
			sy[i+syoffset + 3] = sy[i+syoffset + 3] + sa * sx[i+sxoffset + 3];
		}
		return;
	}
    // code for equal, positive, nonunit increments.

	for (i = 1; i <= n * incx; i = i + incx) {
        sy[i+syoffset] = sa * sx[i+sxoffset] + sy[i+syoffset];
    }
	return;
}