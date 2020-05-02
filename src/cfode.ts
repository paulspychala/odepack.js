/*
-----------------------------------------------------------------------
 cfode is called by the integrator routine to set coefficients
 needed there.  the coefficients for the current method, as
 given by the value of meth, are set for all orders and saved.
 the maximum order assumed here is 12 if meth = 1 and 5 if meth = 2.
 (a smaller value of the maximum order is also allowed.)
 cfode is called once at the beginning of the problem,
 and is not called again unless and until meth is changed.

 the elco array contains the basic method coefficients.
 the coefficients el(i), 1 .le. i .le. nq+1, for the method of
 order nq are stored in elco(i,nq).  they are given by a genetrating
 polynomial, i.e.,
     l(x) = el(1) + el(2)*x + ... + el(nq+1)*x**nq.
 for the implicit adams methods, l(x) is given by
     dl/dx = (x+1)*(x+2)*...*(x+nq-1)/factorial(nq-1),    l(-1) = 0.
 for the bdf methods, l(x) is given by
     l(x) = (x+1)*(x+2)* ... *(x+nq)/k,
 where         k = factorial(nq)*(1 + 1/2 + ... + 1/nq).

 the tesco array contains test constants used for the
 local error test and the selection of step size and/or order.
 at order nq, tesco(k,nq) is used for the selection of step
 size at order nq - 1 if k = 1, at order nq if k = 2, and at order
 nq + 1 if k = 3.

 translation to Typescript by Pawel Spychala 2020
-----------------------------------------------------------------------
*/
export default function cfode(meth: number, elco: number[][], tesco: number[][])
{
	let i: number
	let nq: number
	let nqm1: number
	let nqp1: number;
	let agamq: number
	let fnq: number
	let fnqm1: number
	let pc: number[] = new Array(13)
	let pint: number
	let ragq: number
	let rqfac: number
	let rq1fac: number
	let tsign: number
	let xpin: number;

	if (meth == 1) {
		elco[1][1] = 1.;
		elco[1][2] = 1.;
		tesco[1][1] = 0.;
		tesco[1][2] = 2.;
		tesco[2][1] = 1.;
		tesco[12][3] = 0.;
		pc[1] = 1.;
		rqfac = 1.;
		for (nq = 2; nq <= 12; nq++) {
			/*
			-----------------------------------------------------------------------
			 the pc array will contain the coefficients of the polynomial
			     p(x) = (x+1)*(x+2)*...*(x+nq-1).
			 initially, p(x) = 1.
			-----------------------------------------------------------------------
			*/
			rq1fac = rqfac;
			rqfac = rqfac / nq;
			nqm1 = nq - 1;
			fnqm1 = nqm1;
			nqp1 = nq + 1;

			// form coefficients of p(x)*(x+nq-1). ----------------------------------
			pc[nq] = 0.;
			for (i = nq; i >= 2; i--) pc[i] = pc[i - 1] + fnqm1 * pc[i];
			pc[1] = fnqm1 * pc[1];

			// compute integral, -1 to 0, of p(x) and x*p(x). -----------------------
			pint = pc[1];
			xpin = pc[1] / 2.;
			tsign = 1.;
			for (i = 2; i <= nq; i++) {
				tsign = -tsign;
				pint += tsign * pc[i] / i;
				xpin += tsign * pc[i] / (i + 1);
			}

			// store coefficients in elco and tesco. --------------------------------
			elco[nq][1] = pint * rq1fac;
			elco[nq][2] = 1.;
			for (i = 2; i <= nq; i++) elco[nq][i + 1] = rq1fac * pc[i] / i;
			agamq = rqfac * xpin;
			ragq = 1. / agamq;
			tesco[nq][2] = ragq;
			if (nq < 12) tesco[nqp1][1] = ragq * rqfac / nqp1;
			tesco[nqm1][3] = ragq;
		}		/* end for   */
		return;
	}			/* end if ( meth == 1 )   */
	/*
	   meth = 2.
	*/
	pc[1] = 1.;
	rq1fac = 1.;
	/*
	-----------------------------------------------------------------------
	 the pc array will contain the coefficients of the polynomial
	     p(x) = (x+1)*(x+2)*...*(x+nq).
	 initially, p(x) = 1.
	-----------------------------------------------------------------------
	*/
	for (nq = 1; nq <= 5; nq++) {
		fnq = nq;
		nqp1 = nq + 1;

		// form coefficients of p(x)*(x+nq). ------------------------------------
		pc[nqp1] = 0.;
		for (i = nq + 1; i >= 2; i--)
			pc[i] = pc[i - 1] + fnq * pc[i];
		pc[1] *= fnq;

		// store coefficients in elco and tesco. --------------------------------
		for (i = 1; i <= nqp1; i++)
			elco[nq][i] = pc[i] / pc[2];
		elco[nq][2] = 1.;
		tesco[nq][1] = rq1fac;
		tesco[nq][2] = (nqp1) / elco[nq][1];
		tesco[nq][3] = ((nq + 2)) / elco[nq][1];
		rq1fac /= fnq;
	}
	return;

} // ----------------------- end of subroutine cfode -----------------------