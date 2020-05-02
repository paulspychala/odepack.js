/* 
    finds the index of element having max. absolute value.
    jack dongarra, linpack, 6/17/77.
    translated to typescript by Pawel Spychala 2020
*/
export default function isamax(n: number, sx: number[], sxoffset: number, incx: number) : number
{
    let smax: number
    let xmag: number

    let i: number
    let ii: number
    let indexMax: number

	indexMax = 0;
	if (n <= 0)
		return indexMax;
	indexMax = 1;
	if (n <= 1 || incx <= 0)
		return indexMax;

	// code for increment not equal to 1

	if (incx != 1) {
		smax = Math.abs(sx[1+sxoffset]);
		ii = 2;
		for (i = 1 + incx; i <= n * incx; i = i + incx) {
			xmag = Math.abs(sx[i+sxoffset]);
			if (xmag > smax) {
				indexMax = ii;
				smax = xmag;
			}
			ii++;
		}
		return indexMax;
	}
	// code for increment equal to 1

	smax = Math.abs(sx[1+sxoffset]);
	for (i = 2; i <= n; i++) {
		xmag = Math.abs(sx[i+sxoffset]);
		if (xmag > smax) {
			indexMax = i;
			smax = xmag;
		}
	}
	return indexMax;

}