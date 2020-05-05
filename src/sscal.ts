/* 
    scales a vector by a constant.
    uses unrolled loops for increment equal to 1.
    jack dongarra, linpack, 6/17/77.
    translated to typescript by Pawel Spychala 2020
*/
export function sscal(n: number, sa: number, sx: number[], sxoffset: number, incx: number)
{
    let m: number
    let i: number

	if (n <= 0) return;

    // code for increment not equal to 1

	if (incx != 1) {
		for (i = 1; i <= n * incx; i = i + incx)
			sx[i+sxoffset] = sa * sx[i+sxoffset];
		return;
	}
    // code for increment equal to 1

    // clean-up loop

	m = n % 5;
	if (m != 0) {
		for (i = 1; i <= m; i++)
			sx[i+sxoffset] = sa * sx[i+sxoffset];
		if (n < 5)
			return;
	}
	for (i = m + 1; i <= n; i = i + 5) {
		sx[i+sxoffset] = sa * sx[i+sxoffset];
		sx[i+sxoffset + 1] = sa * sx[i+sxoffset + 1];
		sx[i+sxoffset + 2] = sa * sx[i+sxoffset + 2];
		sx[i+sxoffset + 3] = sa * sx[i+sxoffset + 3];
		sx[i+sxoffset + 4] = sa * sx[i+sxoffset + 4];
	}
	return;

}