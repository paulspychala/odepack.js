/*
-----------------------------------------------------------------------
 this function routine computes the weighted max-norm
 of the vector of length n contained in the array v, with weights
 contained in the array w of length n..
   vmnorm = max(i=1,...,n) abs(v(i))*w(i)

   translated to Typescript by Pawel Spychala 2020
-----------------------------------------------------------------------
*/
export function vmnorm(n: number, v: number[], w: number[]) : number {
	let i: number;
	let vm: number;

	vm = 0.;
	for (i = 1; i <= n; i++)
		vm = Math.max(vm, Math.abs(v[i]) * w[i]);
	return vm;

}