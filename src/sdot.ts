/*
    forms the dot product of a vector.
    uses unrolled loops for increments equal to one.
    jack dongarra, linpack, 6/17/77.
    translated by Pawel Spychala 2020
*/
export default function sdot(n: number, sx: number[], sxoffset: number, incx: number, sy: number[], syoffset: number, incy: number) : number
{
    let dotprod: number
    let ix: number
    let iy: number
    let i: number
    let m: number

	dotprod = 0;
	if (n <= 0)
		return dotprod;

    // code for unequal increments or equal increments not equal to 1

	if (incx != incy || incx < 1) {
		ix = 1;
		iy = 1;
		if (incx < 0) ix = (-n + 1) * incx + 1;
		if (incy < 0) iy = (-n + 1) * incy + 1;
		for (i = 1; i <= n; i++) {
			dotprod = dotprod + sx[ix+sxoffset] * sy[iy+syoffset];
			ix = ix + incx;
			iy = iy + incy;
		}
		return dotprod;
	}
    // code for both increments equal to 1
    // clean-up loop

	if (incx == 1) {
		m = n % 5;
		if (m != 0) {
			for (i = 1; i <= m; i++) {
                dotprod = dotprod + sx[i+sxoffset] * sy[i+syoffset];
            }
			if (n < 5) return dotprod;
		}
		for (i = m + 1; i <= n; i = i + 5) {
			dotprod = dotprod + sx[i+sxoffset] * sy[i+syoffset] + sx[i+sxoffset + 1] * sy[i+syoffset + 1] +
				sx[i+sxoffset + 2] * sy[i+syoffset + 2] + sx[i+sxoffset + 3] * sy[i+syoffset + 3] +
                sx[i+sxoffset + 4] * sy[i+syoffset + 4];
        }
		return dotprod;
	}
    //  code for positive equal nonunit increments.  

	for (i = 1; i <= n * incx; i = i + incx) {
        dotprod = dotprod + sx[i+sxoffset] * sy[i+syoffset];
    }
	return dotprod;

}