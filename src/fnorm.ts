/*
-----------------------------------------------------------------------
 this function computes the norm of a full n by n matrix,
 stored in the array a, that is consistent with the weighted max-norm
 on vectors, with weights stored in the array w..
   fnorm = max(i=1,...,n) ( w(i) * sum(j=1,...,n) abs(a(i,j))/w(j) )

   translated to typescript by Pawel Spychala 2020
-----------------------------------------------------------------------
*/
export default function fnorm(n: number, a: number[][], w: number[]) : number {
	let i: number
	let j: number;
	let an: number
	let sum: number
	let ap1: number[];

	an = 0.;
	for (i = 1; i <= n; i++) {
		sum = 0.;
		ap1 = a[i];
		for (j = 1; j <= n; j++) { sum += Math.abs(ap1[j]) / w[j]; }
		an = Math.max(an, sum * w[i]);
	}
	return an;

}