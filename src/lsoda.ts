/* 
  This is a typescript version of the LSODA library. The original is written in fortran and I also
  used Heng Li's C port as a guide. Since both fortran and c use pointer math I had to rewrite
  sections of the program.

  - Pawel Spychala <paulspychala@gmail.com>

  I kept most of the previous author's notes in here.

  Thanks to Heng Li's c port: 
  http://lh3lh3.users.sourceforge.net/download/this.c

*/
import {ewset} from "./ewset";
import {sgesl} from "./sgesl";
import {sgbsl} from "./sgbsl";
import {sgefa} from "./sgefa";
import {sgbfa} from "./sgbfa";
import {vmnorm} from "./vmnorm";
import {fnorm} from "./fnorm";
import {cfode} from "./cfode";

export type lsoda_func = (t: number, y: number[], ydot: number[], data: any) => any
export type jac_func = (t: number, y: number[], ml: number, mu: number, wm: number[][], data: any) => any

export class LSODA {

create2dArray(size1: number, size2: number) : number[][] {
	let ret = new Array(size1)
	for (let i = 0; i < size1; i++) {
		ret[i] = new Array(size2)
	}
	return ret
}

readonly ETA = 2.2204460492503131e-16

/* newly added static variables */

ml: number
mu: number
imxer: number
mord: number[] = [0, 12, 5]
sqrteta: number = Math.sqrt(this.ETA)
yp1: number[]
yp2: number[]
sm1: number[] = [0., 0.5, 0.575, 0.55, 0.45, 0.35, 0.25, 0.2, 0.15, 0.1, 0.075, 0.05, 0.025]

/* static variables for lsoda() */

ccmax: number
el0: number
h: number
hmin: number
hmxi: number
hu: number
rc: number
tn: number
illin: number = 0
init: number = 0
mxstep: number
mxhnil: number
nhnil: number
ntrep: number = 0
nslast: number
nyh: number
ierpj: number
iersl: number
jcur: number
jstart: number
kflag: number
l: number
meth: number
miter: number
maxord: number
maxcor: number
msbp: number
mxncf: number
n: number
nq: number
nst: number
nfe: number
nje: number
nqu: number
tsw: number
pdnorm: number
ixpr: number = 0
jtyp: number
mused: number
mxordn: number
mxords: number

/* no static variable for prja(), solsy() */
/* static variables for stoda() */

conit: number
crate: number
el: number[] = new Array(14)
elco: number[][] = this.create2dArray(13,14) 
hold: number
rmax: number
tesco: number[][] = this.create2dArray(13,4) 
ialth: number
ipup: number
lmax: number
nslp: number
pdest: number
pdlast: number
ratio: number
cm1: number[] = new Array(13)
cm2: number[] = new Array(6)
icount: number
irflag: number

/* static variables for various vectors and the Jacobian. */

yh: number[][]
wm: number[][]
ewt: number[]
savf: number[]
acor: number[]
ipvt: number[]

/*
   The following are useful statistics.
   hu,
   h,
   tn,
   tolsf,
   tsw,
   nst,
   nfe,
   nje,
   nqu,
   nq,
   imxer,
   mused,
   meth
*/


/* Terminate lsoda due to illegal input. */
terminate(istate: number) : number
{
	if (this.illin == 5) {
		console.error("[lsoda] repeated occurrence of illegal input. run aborted.. apparent infinite loop\n");
	} else {
		this.illin++;
		istate = -3;
    }
	return istate
}


/* Terminate lsoda due to various error conditions. */
terminate2(y: number[], t: number) : number
{
    let i: number
	this.yp1 = this.yh[1];
	for (i = 1; i <= this.n; i++)
		y[i] = this.yp1[i];
	t = this.tn;
	this.illin = 0;
	return t;
}

/*
-----------------------------------------------------------------------
 block g.
 the following block handles all successful returns from lsoda.
 if itask .ne. 1, y is loaded from yh and t is set accordingly.
 istate is set to 2, the illegal input counter is zeroed, and the
 optional outputs are loaded into the work arrays before returning.
 if istate = 1 and tout = t, there is a return with no action taken,
 except that if this has happened repeatedly, the run is terminated.
-----------------------------------------------------------------------
*/
successreturn(y: number[], t: number, itask: number, ihit: boolean, tcrit: number, istate: number) : [number,number]
{
	this.yp1 = this.yh[1];
	for (let i = 1; i <= this.n; i++)
		y[i] = this.yp1[i];
	t = this.tn;
	if (itask == 4 || itask == 5)
		if (ihit)
			t = tcrit;
	istate = 2;
	this.illin = 0;
	return [t, istate]
}

/*
-----------------------------------------------------------------------
 This is the March 30, 1987 version of
 lsoda.. Livermore solver for ordinary differential equations, with
         automatic method switching for stiff and nonstiff problems.

 This version is in double precision.

 lsoda solves the initial value problem for stiff or nonstiff
 systems of first order ode-s,
     dy/dt = f(t,y) ,  or, in component form,
     dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(neq)) (i = 1,...,neq).

 This a variant version of the lsode package.
 it switches automatically between stiff and nonstiff methods.
 This means that the user does not have to determine whether the
 problem is stiff or not, and the solver will automatically choose the
 appropriate method.  it always starts with the nonstiff method.

 authors..
                linda r. petzold  and  alan c. hindmarsh,
                computing and mathematics research division, l-316
                lawrence livermore national laboratory
                livermore, ca 94550.

 references..
 1.  alan c. hindmarsh,  odepack, a systematized collection of ode
     solvers, in scientific computing, r. s. stepleman et al. (eds.),
     north-holland, amsterdam, 1983, pp. 55-64.
 2.  linda r. petzold, automatic selection of methods for solving
     stiff and nonstiff systems of ordinary differential equations,
     siam j. sci. stat. comput. 4 (1983), pp. 136-148.
-----------------------------------------------------------------------
 summary of usage.

 communication between the user and the lsoda package, for normal
 situations, is summarized here.  this summary describes only a subset
 of the full set of options available.  see the full description for
 details, including alternative treatment of the jacobian matrix,
 optional inputs and outputs, nonstandard options, and
 instructions for special situations.  see also the example
 problem (with program and output) following this summary.

 a. first provide a subroutine of the form..
               subroutine f (neq, t, y, ydot)
               dimension y(neq), ydot(neq)
 which supplies the vector function f by loading ydot(i) with f(i).

 b. write a main program which calls subroutine lsoda once for
 each point at which answers are desired.  this should also provide
 for possible use of logical unit 6 for output of error messages
 by lsoda.  on the first call to lsoda, supply arguments as follows..
 f      = name of subroutine for right-hand side vector f.
          this name must be declared external in calling program.
 neq    = number of first order ode-s.
 y      = array of initial values, of length neq.
 t      = the initial value of the independent variable.
 tout   = first point where output is desired (.ne. t).
 itol   = 1 or 2 according as atol (below) is a scalar or array.
 rtol   = relative tolerance parameter (scalar).
 atol   = absolute tolerance parameter (scalar or array).
          the estimated local error in y(i) will be controlled so as
          to be less than
             ewt(i) = rtol*abs(y(i)) + atol     if itol = 1, or
             ewt(i) = rtol*abs(y(i)) + atol(i)  if itol = 2.
          thus the local error test passes if, in each component,
          either the absolute error is less than atol (or atol(i)),
          or the relative error is less than rtol.
          use rtol = 0.0 for pure absolute error control, and
          use atol = 0.0 (or atol(i) = 0.0) for pure relative error
          control.  caution.. actual (global) errors may exceed these
          local tolerances, so choose them conservatively.
 itask  = 1 for normal computation of output values of y at t = tout.
 istate = integer flag (input and output).  set istate = 1.
 iopt   = 0 to indicate no optional inputs used.
 rwork  = real work array of length at least..
             22 + neq * max(16, neq + 9).
          see also paragraph e below.
 lrw    = declared length of rwork (in user-s dimension).
 iwork  = integer work array of length at least  20 + neq.
 liw    = declared length of iwork (in user-s dimension).
 jac    = name of subroutine for jacobian matrix.
          use a dummy name.  see also paragraph e below.
 jt     = jacobian type indicator.  set jt = 2.
          see also paragraph e below.
 note that the main program must declare arrays y, rwork, iwork,
 and possibly atol.

 c. the output from the first call (or any call) is..
      y = array of computed values of y(t) vector.
      t = corresponding value of independent variable (normally tout).
 istate = 2  if lsoda was successful, negative otherwise.
          -1 means excess work done on this call (perhaps wrong jt).
          -2 means excess accuracy requested (tolerances too small).
          -3 means illegal input detected (see printed message).
          -4 means repeated error test failures (check all inputs).
          -5 means repeated convergence failures (perhaps bad jacobian
             supplied or wrong choice of jt or tolerances).
          -6 means error weight became zero during problem. (solution
             component i vanished, and atol or atol(i) = 0.)
          -7 means work space insufficient to finish (see messages).

 d. to continue the integration after a successful return, simply
 reset tout and call lsoda again.  no other parameters need be reset.

 e. note.. if and when lsoda regards the problem as stiff, and
 switches methods accordingly, it must make use of the neq by neq
 jacobian matrix, j = df/dy.  for the sake of simplicity, the
 inputs to lsoda recommended in paragraph b above cause lsoda to
 treat j as a full matrix, and to approximate it internally by
 difference quotients.  alternatively, j can be treated as a band
 matrix (with great potential reduction in the size of the rwork
 array).  also, in either the full or banded case, the user can supply
 j in closed form, with a routine whose name is passed as the jac
 argument.  these alternatives are described in the paragraphs on
 rwork, jac, and jt in the full description of the call sequence below.

-----------------------------------------------------------------------
*/

lsoda(f: lsoda_func, neq: number, y: number[], t: number, tout: number, 
			itol: number, rtol: number[], atol: number[],
			itask: number, istate: number, iopt: number, jac: jac_func, jt: number,
			iwork1: number, iwork2: number, iwork5: number, iwork6: number,
			iwork7: number, iwork8: number, iwork9: number,
			rwork1: number, rwork5: number, rwork6: number, rwork7: number, 
			_data: any) : [number, number]
/*
   If the user does not supply any of these values, the calling program
   should initialize those untouched working variables to zero.
   ml = iwork1
   mu = iwork2
   ixpr = iwork5
   mxstep = iwork6
   mxhnil = iwork7
   mxordn = iwork8
   mxords = iwork9
   tcrit = rwork1
   h0 = rwork5
   hmax = rwork6
   hmin = rwork7
*/


{
    let mxstp0: number = 500
    let mxhnl0: number = 10

	let i: number
	let iflag: number
	let lenyh: number
	let ihit: boolean
	let atoli: number
	let ayi: number
	let big: number
	let h0: number
	let hmax: number
	let hmx: number
	let rh: number
	let rtoli: number
	let tcrit: number
	let tdist: number
	let tnext: number
	let tol: number
	let tolsf: number
	let tp: number
	let size: number
	let sum: number
	let w0: number;

/*
	-----------------------------------------------------------------------
	 block a.
	 this code block is executed on every call.
	 it tests istate and itask for legality and branches appropriately.
	 if istate .gt. 1 but the flag init shows that initialization has
	 not yet been done, an error return occurs.
	 if istate = 1 and tout = t, jump to block g and return immediately.
	-----------------------------------------------------------------------
*/
	if (istate < 1 || istate > 3) {
		console.error("[lsoda] illegal istate = %d\n", istate);
		return [t,this.terminate(istate)];
	}
	if (itask < 1 || itask > 5) {
		console.error("[lsoda] illegal itask = %d\n", itask);
		return [t,this.terminate(istate)];
	}
	if (this.init == 0 && (istate == 2 || istate == 3)) {
		console.error("[lsoda] istate > 1 but lsoda not initialized\n");
		return [t,this.terminate(istate)];
	}
	if (istate == 1) {
		this.init = 0;
		if (tout == t) {
			this.ntrep++;
			if (this.ntrep < 5) return [t,istate];
			console.error("[lsoda] repeated calls with istate = 1 and tout = t. run aborted.. apparent infinite loop\n");
			return [t,istate];
		}
	}
/*
	-----------------------------------------------------------------------
	 block b.
	 the next code block is executed for the initial call (istate = 1),
	 or for a continuation call with parameter changes (istate = 3).
	 it contains checking of all inputs and various initializations.

	 first check legality of the non-optional inputs neq, itol, iopt,
	 jt, ml, and mu.
	-----------------------------------------------------------------------
*/

	if (istate == 1 || istate == 3) {
		this.ntrep = 0;
		if (neq <= 0) {
			console.error(`[lsoda] neq = ${neq} is less than 1`);
			return [t,this.terminate(istate)];
		}
		if (istate == 3 && neq > this.n) {
			console.error("[lsoda] istate = 3 and neq increased");
			return [t,this.terminate(istate)];
		}
		this.n = neq;
		if (itol < 1 || itol > 4) {
			console.error(`[lsoda] itol = ${itol} illegal`);
			return [t,this.terminate(istate)];
		}
		if (iopt < 0 || iopt > 1) {
			console.error(`[lsoda] iopt = ${iopt} illegal`);
			return [t,this.terminate(istate)];
		}
		if (jt == 3 || jt < 1 || jt > 5) {
			console.error(`[lsoda] jt = ${jt} illegal`);
			return [t,this.terminate(istate)];
		}
		this.jtyp = jt;
		if (jt > 2) {
			this.ml = iwork1;
			this.mu = iwork2;
			if (this.ml < 0 || this.ml >= this.n) {
				console.error(`[lsoda] ml = ${this.ml} not between 1 and neq.`);
				return [t,this.terminate(istate)];
			}
			if (this.mu < 0 || this.mu >= this.n) {
				console.error(`[lsoda] mu = ${this.mu} not between 1 and neq.`);
				return [t,this.terminate(istate)];
			}
		}
		// next process and check the optional inputs. --------------------------
		if (iopt == 0) {
			this.ixpr = 0;
			this.mxstep = mxstp0;
			this.mxhnil = mxhnl0;
			this.hmxi = 0.;
			this.hmin = 0.;
			if (istate == 1) {
				h0 = 0.;
				this.mxordn = this.mord[1];
				this.mxords = this.mord[2];
			}
		} else {		/* if ( iopt = 1 )  */
			this.ixpr = iwork5;
			if (this.ixpr < 0 || this.ixpr > 1) {
				console.error(`[lsoda] ixpr = ${this.ixpr} is illegal`);
				return [t,this.terminate(istate)];
			}
			this.mxstep = iwork6;
			if (this.mxstep < 0) {
				console.error("[lsoda] mxstep < 0\n");
				return [t,this.terminate(istate)];
			}
			if (this.mxstep == 0) this.mxstep = mxstp0;
			this.mxhnil = iwork7;
			if (this.mxhnil < 0) {
				console.error("[lsoda] mxhnil < 0\n");
				return [t,this.terminate(istate)];
			}
			if (istate == 1) {
				h0 = rwork5;
				this.mxordn = iwork8;
				if (this.mxordn < 0) {
					console.error(`[lsoda] mxordn = ${this.mxordn} is less than 0`);
					return [t,this.terminate(istate)];
				}
				if (this.mxordn == 0) this.mxordn = 100;
				this.mxordn = Math.min(this.mxordn, this.mord[1]);
				this.mxords = iwork9;
				if (this.mxords < 0) {
					console.error(`[lsoda] mxords = ${this.mxords} is less than 0`);
					return [t,this.terminate(istate)];
				}
				if (this.mxords == 0) this.mxords = 100;
				this.mxords = Math.min(this.mxords, this.mord[2]);
				if ((tout - t) * h0 < 0.) {
					console.error(`[lsoda] tout = ${tout} behind t = ${t}. integration direction is given by ${h0}`);
					return [t,this.terminate(istate)];
				}
			}	/* end if ( *istate == 1 )  */
			hmax = rwork6;
			if (hmax < 0.) {
				console.error("[lsoda] hmax < 0.\n");
				return [t,this.terminate(istate)];
			}
			this.hmxi = 0.;
			if (hmax > 0)
				this.hmxi = 1. / hmax;
			this.hmin = rwork7;
			if (this.hmin < 0.) {
				console.error("[lsoda] hmin < 0.\n");
				return [t,this.terminate(istate)];
			}
		}		/* end else   *//* end iopt = 1   */
	}			/* end if ( *istate == 1 || *istate == 3 )   */
	if (istate == 1) {
/*
		-----------------------------------------------------------------------
		 set work array pointers and check lengths lrw and liw.
		 if istate = 1, meth is initialized to 1 here to facilitate the
		 checking of work space lengths.
		 pointers to segments of rwork and iwork are named by prefixing l to
		 the name of the segment.  e.g., the segment yh starts at rwork(lyh).
		 segments of rwork (in order) are denoted  yh, wm, ewt, savf, acor.
		 if the lengths provided are insufficient for the current method,
		 an error return occurs.  this is treated as illegal input on the
		 first call, but as a problem interruption with istate = -7 on a
		 continuation call.  if the lengths are sufficient for the current
		 method but not for both methods, a warning message is sent.
		-----------------------------------------------------------------------
*/
		this.meth = 1;
		this.nyh = this.n;
		lenyh = 1 + Math.max(this.mxordn, this.mxords);

		// allocate memory for yh, wm, ewt, savf, acor, ipvt.
		this.yh = new Array(1 + lenyh)
		for (i = 1; i <= lenyh; i++){
			this.yh[i] = new Array(1 + this.nyh)
		}

		this.wm = new Array(1 + this.nyh)
		for (i = 1; i <= this.nyh; i++){
			this.wm[i] = new Array(1 + this.nyh)
		}

		this.ewt = new Array(1 + this.nyh)
		this.savf = new Array(1 + this.nyh)
		this.acor = new Array(1 + this.nyh)
		this.ipvt = new Array(1 + this.nyh)
	}
	// check rtol and atol for legality. ------------------------------------
	if (istate == 1 || istate == 3) {
		rtoli = rtol[1];
		atoli = atol[1];
		for (i = 1; i <= this.n; i++) {
			if (itol >= 3)
				rtoli = rtol[i];
			if (itol == 2 || itol == 4)
				atoli = atol[i];
			if (rtoli < 0.) {
				console.error(`[lsoda] rtol = ${rtoli} is less than 0.`);
				return [t,this.terminate(istate)];
			}
			if (atoli < 0.) {
				console.error(`[lsoda] atol = ${atoli} is less than 0.`);
				return [t,this.terminate(istate)];
			}
		}
	}
	// if istate = 3, set flag to signal parameter changes to stoda. --------
	if (istate == 3) {
		this.jstart = -1;
	}
/*
	-----------------------------------------------------------------------
	 block c.
	 the next block is for the initial call only (istate = 1).
	 it contains all remaining initializations, the initial call to f,
	 and the calculation of the initial step size.
	 the error weights in ewt are inverted after being loaded.
	-----------------------------------------------------------------------
*/
	if (istate == 1) {
		this.tn = t;
		this.tsw = t;
		this.maxord = this.mxordn;
		if (itask == 4 || itask == 5) {
			tcrit = rwork1;
			if ((tcrit - tout) * (tout - t) < 0.) {
				console.error("[lsoda] itask = 4 or 5 and tcrit behind tout\n");
				return [t,this.terminate(istate)];
			}
			if (h0 != 0. && (t + h0 - tcrit) * h0 > 0.)
				h0 = tcrit - t;
		}
		this.jstart = 0;
		this.nhnil = 0;
		this.nst = 0;
		this.nje = 0;
		this.nslast = 0;
		this.hu = 0.;
		this.nqu = 0;
		this.mused = 0;
		this.miter = 0;
		this.ccmax = 0.3;
		this.maxcor = 3;
		this.msbp = 20;
		this.mxncf = 10;
		// initial call to f.
		f(t,y, this.yh[2], _data)
		this.nfe = 1;
   		// load the initial value vector in yh.
		this.yp1 = this.yh[1];
		for (i = 1; i <= this.n; i++){
			this.yp1[i] = y[i];
		}
		// load and invert the ewt array.  (h is temporarily set to 1.0.) -------
		this.nq = 1;
		this.h = 1.;
		ewset(this.n, itol, rtol, atol, y, this.ewt);
		for (i = 1; i <= this.n; i++) {
			if (this.ewt[i] <= 0.) {
				console.error(`[lsoda] ewt[${i}] = ${this.ewt[i]} <= 0.`);
				t = this.terminate2(y, t);
				return [t,istate];
			}
			this.ewt[i] = 1. / this.ewt[i];
		}
/*
		-----------------------------------------------------------------------
		 the coding below computes the step size, h0, to be attempted on the
		 first step, unless the user has supplied a value for this.
		 first check that tout - t differs significantly from zero.
		 a scalar tolerance quantity tol is computed, as max(rtol(i))
		 if this is positive, or max(atol(i)/abs(y(i))) otherwise, adjusted
		 so as to be between 100*uround and 1.0e-3.
		 then the computed value h0 is given by..

		   h0**(-2)  =  1./(tol * w0**2)  +  tol * (norm(f))**2

		 where   w0     = max ( abs(t), abs(tout) ),
		         f      = the initial value of the vector f(t,y), and
		         norm() = the weighted vector norm used throughout, given by
		                  the vmnorm function routine, and weighted by the
		                  tolerances initially loaded into the ewt array.
		 the sign of h0 is inferred from the initial values of tout and t.
		 abs(h0) is made .le. abs(tout-t) in any case.
		-----------------------------------------------------------------------
*/
		if (h0 == 0.) {
			tdist = Math.abs(tout - t);
			w0 = Math.max(Math.abs(t), Math.abs(tout));
			if (tdist < 2. * this.ETA * w0) {
				console.error("[lsoda] tout too close to t to start integration\n ");
				return [t,this.terminate(istate)];
			}
			tol = rtol[1];
			if (itol > 2) {
				for (i = 2; i <= this.n; i++)
					tol = Math.max(tol, rtol[i]);
			}
			if (tol <= 0.) {
				atoli = atol[1];
				for (i = 1; i <= this.n; i++) {
					if (itol == 2 || itol == 4)
						atoli = atol[i];
					ayi = Math.abs(y[i]);
					if (ayi != 0.)
						tol = Math.max(tol, atoli / ayi);
				}
			}
			tol = Math.max(tol, 100. * this.ETA);
			tol = Math.min(tol, 0.001);
			sum = vmnorm(this.n, this.yh[2], this.ewt);
			sum = 1. / (tol * w0 * w0) + tol * sum * sum;
			h0 = 1. / Math.sqrt(sum);
			h0 = Math.min(h0, tdist);
			h0 = h0 * ((tout - t >= 0.) ? 1. : -1.);
		}		/* end if ( h0 == 0. )   */
		// adjust h0 if necessary to meet hmax bound. ---------------------------
		rh = Math.abs(h0) * this.hmxi;
		if (rh > 1.){
			h0 /= rh;
		}
		// load h with h0 and scale yh(*,2) by h0. ------------------------------
		this.h = h0;
		this.yp1 = this.yh[2];
		for (i = 1; i <= this.n; i++){
			this.yp1[i] *= h0;
		}
	}			/* if ( *istate == 1 )   */
/*
	-----------------------------------------------------------------------
	 block d.
	 the next code block is for continuation calls only (istate = 2 or 3)
	 and is to check stop conditions before taking a step.
	-----------------------------------------------------------------------
*/
	if (istate == 2 || istate == 3) {
		this.nslast = this.nst;
		switch (itask) {
		case 1:
			if ((this.tn - tout) * this.h >= 0.) {
				iflag = this.intdy(tout, 0, y, iflag);
				if (iflag != 0) {
					console.error("[lsoda] trouble from intdy, itask = %d, tout = %g\n", itask, tout);
					return [t,this.terminate(istate)];
				}
				t = tout;
				istate = 2;
				this.illin = 0;
				return [t,istate];
			}
			break;
		case 2:
			break;
		case 3:
			tp = this.tn - this.hu * (1. + 100. * this.ETA);
			if ((tp - tout) * this.h > 0.) {
				console.error("[lsoda] itask = %d and tout behind tcur - hu\n", itask);
				return [t,this.terminate(istate)];
			}
			if ((this.tn - tout) * this.h < 0.) break;
			return this.successreturn(y, t, itask, ihit, tcrit, istate);
		case 4:
			tcrit = rwork1;
			if ((this.tn - tcrit) * this.h > 0.) {
				console.error("[lsoda] itask = 4 or 5 and tcrit behind tcur\n");
				return [t,this.terminate(istate)];
			}
			if ((tcrit - tout) * this.h < 0.) {
				console.error("[lsoda] itask = 4 or 5 and tcrit behind tout\n");
				return [t,this.terminate(istate)];
			}
			if ((this.tn - tout) * this.h >= 0.) {
				iflag = this.intdy(tout, 0, y, iflag);
				if (iflag != 0) {
					console.error("[lsoda] trouble from intdy, itask = %d, tout = %g\n", itask, tout);
					return [t,this.terminate(istate)];
				}
				t = tout;
				istate = 2;
				this.illin = 0;
				return [t,istate];
			}
		case 5:
			if (itask == 5) {
				tcrit = rwork1;
				if ((this.tn - tcrit) * this.h > 0.) {
					console.error("[lsoda] itask = 4 or 5 and tcrit behind tcur\n");
					return [t,this.terminate(istate)];
				}
			}
			hmx = Math.abs(this.tn) + Math.abs(this.h);
			ihit = Math.abs(this.tn - tcrit) <= (100. * this.ETA * hmx);
			if (ihit) {
				t = tcrit;
				return this.successreturn(y, t, itask, ihit, tcrit, istate);
			}
			tnext = this.tn + this.h * (1. + 4. * this.ETA);
			if ((tnext - tcrit) * this.h <= 0.){
				break;
			}
			this.h = (tcrit - this.tn) * (1. - 4. * this.ETA);
			if (istate == 2)
			this.jstart = -2;
			break;
		}		/* end switch   */
	}			/* end if ( *istate == 2 || *istate == 3 )   */
	/*
	-----------------------------------------------------------------------
	 block e.
	 the next block is normally executed for all calls and contains
	 the call to the one-step core integrator stoda.

	 this is a looping point for the integration steps.

	 first check for too many steps being taken, update ewt (if not at
	 start of problem), check for too much accuracy being requested, and
	 check for h below the roundoff level in t.
	-----------------------------------------------------------------------
	*/
	while (1) {
		if (istate != 1 || this.nst != 0) {
			if ((this.nst - this.nslast) >= this.mxstep) {
				console.error("[lsoda] %d steps taken before reaching tout\n", this.mxstep);
				istate = -1;
				t = this.terminate2(y, t);
				return [t,istate];
			}
			ewset(this.n,itol, rtol, atol, this.yh[1], this.ewt);
			for (i = 1; i <= this.n; i++) {
				if (this.ewt[i] <= 0.) {
					console.error("[lsoda] ewt[%d] = %g <= 0.\n", i, this.ewt[i]);
					istate = -6;
					t = this.terminate2(y, t);
					return [t,istate];
				}
				this.ewt[i] = 1. / this.ewt[i];
			}
		}
		tolsf = this.ETA * vmnorm(this.n, this.yh[1], this.ewt);
		if (tolsf > 0.01) {
			tolsf = tolsf * 200.;
			if (this.nst == 0) {
				console.error("lsoda -- at start of problem, too much accuracy\n");
				console.error("         requested for precision of machine,\n");
				console.error("         suggested scaling factor = %g\n", tolsf);
				return [t,this.terminate(istate)];
			}
			console.error("lsoda -- at t = %g, too much accuracy requested\n", t);
			console.error("         for precision of machine, suggested\n");
			console.error("         scaling factor = %g\n", tolsf);
			istate = -2;
			t = this.terminate2(y, t);
			return [t,istate];
		}
		if ((this.tn + this.h) == this.tn) {
			this.nhnil++;
			if (this.nhnil <= this.mxhnil) {
				console.error("lsoda -- warning..internal t = %g and h = %g are\n", this.tn, this.h);
				console.error("         such that in the machine, t + h = t on the next step\n");
				console.error("         solver will continue anyway.\n");
				if (this.nhnil == this.mxhnil) {
					console.error("lsoda -- above warning has been issued %d times,\n", this.nhnil);
					console.error("         it will not be issued again for this problem\n");
				}
			}
		}
/*
		-----------------------------------------------------------------------
		     call stoda(neq,y,yh,nyh,yh,ewt,savf,acor,wm,iwm,f,jac,prja,solsy)
		-----------------------------------------------------------------------
*/
		this.stoda(neq, y, f, jac, _data);

		if (this.kflag == 0) {
/*
		-----------------------------------------------------------------------
		 block f.
		 the following block handles the case of a successful return from the
		 core integrator (kflag = 0).
		 if a method switch was just made, record tsw, reset maxord,
		 set jstart to -1 to signal stoda to complete the switch,
		 and do extra printing of data if ixpr = 1.
		 then, in any case, check for stop conditions.
		-----------------------------------------------------------------------
*/
			this.init = 1;
			if (this.meth != this.mused) {
				this.tsw = this.tn;
				this.maxord = this.mxordn;
				if (this.meth == 2)
				this.maxord = this.mxords;
				this.jstart = -1;
				if (this.ixpr) {
					if (this.meth == 2)
						console.error("[lsoda] a switch to the stiff method has occurred ");
					if (this.meth == 1)
						console.error("[lsoda] a switch to the nonstiff method has occurred");
					console.error("at t = %g, tentative step size h = %g, step nst = %d\n", this.tn, this.h, this.nst);
				}
			}	/* end if ( meth != mused )   */
			// itask = 1.  if tout has been reached, interpolate. -------------------
			if (itask == 1) {
				if ((this.tn - tout) * this.h < 0.)
					continue;
					iflag = this.intdy(tout, 0, y, iflag);
				t = tout;
				istate = 2;
				this.illin = 0;
				return [t,istate];
			}
   			// itask = 2.
			if (itask == 2) {
				return this.successreturn(y, t, itask, ihit, tcrit, istate);
			}
			// itask = 3.  jump to exit if tout was reached. ------------------------
			if (itask == 3) {
				if ((this.tn - tout) * this.h >= 0.) {
					return this.successreturn(y, t, itask, ihit, tcrit, istate);
				}
				continue;
			}
			// c itask = 4.  see if tout or tcrit was reached.  adjust h if necessary.
			if (itask == 4) {
				if ((this.tn - tout) * this.h >= 0.) {
					iflag = this.intdy(tout, 0, y, iflag);
					t = tout;
					istate = 2;
					this.illin = 0;
					return [t,istate];
				} else {
					hmx = Math.abs(this.tn) + Math.abs(this.h);
					ihit = Math.abs(this.tn - tcrit) <= (100. * this.ETA * hmx);
					if (ihit) {
						return this.successreturn(y, t, itask, ihit, tcrit, istate);
					}
					tnext = this.tn + this.h * (1. + 4. * this.ETA);
					if ((tnext - tcrit) * this.h <= 0.) continue;
					this.h = (tcrit - this.tn) * (1. - 4. * this.ETA);
					this.jstart = -2;
					continue;
				}
			}
			// itask = 5.  see if tcrit was reached and jump to exit. ---------------
			if (itask == 5) {
				hmx = Math.abs(this.tn) + Math.abs(this.h);
				ihit = Math.abs(this.tn - tcrit) <= (100. * this.ETA * hmx);
				return this.successreturn(y, t, itask, ihit, tcrit, istate);
			}
		}		/* end if ( kflag == 0 )   */
		/*
		   kflag = -1, error test failed repeatedly or with fabs(h) = hmin.
		   kflag = -2, convergence failed repeatedly or with fabs(h) = hmin.
		*/
		if (this.kflag == -1 || this.kflag == -2) {
			console.error(`lsoda -- at t = ${this.tn} and step size h = ${this.h}, the`);
			if (this.kflag == -1) {
				console.error("         error test failed repeatedly or");
				console.error("         with abs(h) = hmin");
				istate = -4;
			}
			if (this.kflag == -2) {
				console.error("         corrector convergence failed repeatedly or\n");
				console.error("         with abs(h) = hmin\n");
				istate = -5;
			}
			big = 0.;
			this.imxer = 1;
			for (i = 1; i <= this.n; i++) {
				size = Math.abs(this.acor[i]) * this.ewt[i];
				if (big < size) {
					big = size;
					this.imxer = i;
				}
			}
			t = this.terminate2(y, t);
			return [t,istate];
		}		/* end if ( kflag == -1 || kflag == -2 )   */
	}			/* end while   */
	return [t,istate];
} /* end lsoda   */


stoda(neq: number, y: number[], f: any, jac: any, _data: any[])
{
	let corflag: number
	let orderflag: number;
	let i: number
	let i1: number
	let j: number
	let m: number
	let ncf: number;
	let del: number
	let delp: number
	let dsm: number
	let dup: number
	let exup: number
	let r: number
	let rh: number
	let rhup: number
	let told: number;
	let pdh: number
	let pnorm: number;

/*
   stoda performs one step of the integration of an initial value
   problem for a system of ordinary differential equations.
   Note.. stoda is independent of the value of the iteration method
   indicator miter, when this is != 0, and hence is independent
   of the type of chord method used, or the Jacobian structure.
   Communication with stoda is done with the following variables:
   jstart = an integer used for input only, with the following
            values and meanings:
               0  perform the first step,
             > 0  take a new step continuing from the last,
              -1  take the next step with a new value of h,
                  n, meth, miter, and/or matrix parameters.
              -2  take the next step with a new value of h,
                  but with other inputs unchanged.
   kflag = a completion code with the following meanings:
             0  the step was successful,
            -1  the requested error could not be achieved,
            -2  corrector convergence could not be achieved,
            -3  fatal error in prja or solsy.
   miter = corrector iteration method:
             0  functional iteration,
            >0  a chord method corresponding to jacobian type jt.
*/
	this.kflag = 0;
	told = this.tn;
	ncf = 0;
	this.ierpj = 0;
	this.iersl = 0;
	this.jcur = 0;
	delp = 0.;

/*
	-----------------------------------------------------------------------
	 on the first call, the order is set to 1, and other variables are
	 initialized.  rmax is the maximum ratio by which h can be increased
	 in a single step.  it is initially 1.e4 to compensate for the small
	 initial h, but then is normally equal to 10.  if a failure
	 occurs (in corrector convergence or error test), rmax is set at 2
	 for the next increase.
	 cfode is called to get the needed coefficients for both methods.
	-----------------------------------------------------------------------
*/
	if (this.jstart == 0) {
		this.lmax = this.maxord + 1;
		this.nq = 1;
		this.l = 2;
		this.ialth = 2;
		this.rmax = 10000.;
		this.rc = 0.;
		this.el0 = 1.;
		this.crate = 0.7;
		this.hold = this.h;
		this.nslp = 0;
		this.ipup = this.miter;
		// initialize switching parameters.  meth = 1 is assumed initially. -----
		this.icount = 20;
		this.irflag = 0;
		this.pdest = 0.;
		this.pdlast = 0.;
		this.ratio = 5.;
		cfode(2,this.elco, this.tesco);
		for (i = 1; i <= 5; i++) this.cm2[i] = this.tesco[i][2] * this.elco[i][i + 1];

		cfode(1,this.elco,this.tesco);
		for (i = 1; i <= 12; i++) this.cm1[i] = this.tesco[i][2] * this.elco[i][i + 1];
		this.resetcoeff();
	}			/* end if ( jstart == 0 )   */
	/*
	-----------------------------------------------------------------------
	 the following block handles preliminaries needed when jstart = -1.
	 ipup is set to miter to force a matrix update.
	 if an order increase is about to be considered (ialth = 1),
	 ialth is reset to 2 to postpone consideration one more step.
	 if the caller has changed meth, cfode is called to reset
	 the coefficients of the method.
	 if h is to be changed, yh must be rescaled.
	 if h or meth is being changed, ialth is reset to l = nq + 1
	 to prevent further changes in h for that many steps.
	-----------------------------------------------------------------------
	*/
	if (this.jstart == -1) {
		this.ipup = this.miter;
		this.lmax = this.maxord + 1;
		if (this.ialth == 1) this.ialth = 2;
		if (this.meth != this.mused) {
			cfode(this.meth,this.elco,this.tesco);
			this.ialth = this.l;
			this.resetcoeff();
		}
		if (this.h != this.hold) {
			rh = this.h / this.hold;
			this.h = this.hold;
			[rh, pdh] = this.scaleh(rh, pdh);
		}
	}			/* if ( jstart == -1 )   */
	if (this.jstart == -2) {
		if (this.h != this.hold) {
			rh = this.h / this.hold;
			this.h = this.hold;
			[rh, pdh] = this.scaleh(rh, pdh);
		}
	}			/* if ( jstart == -2 )   */
	/*
	-----------------------------------------------------------------------
	 this section computes the predicted values by effectively
	 multiplying the yh array by the pascal triangle matrix.
	 rc is the ratio of new to old values of the coefficient  h*el(1).
	 when rc differs from 1 by more than ccmax, ipup is set to miter
	 to force pjac to be called, if a jacobian is involved.
	 in any case, pjac is called at least every msbp steps.
	-----------------------------------------------------------------------
	*/
	while (1) {
		while (1) {
			if (Math.abs(this.rc - 1.) > this.ccmax) this.ipup = this.miter;
			if (this.nst >= this.nslp + this.msbp) this.ipup = this.miter;
			this.tn += this.h;
			for (j = this.nq; j >= 1; j--){
				for (i1 = j; i1 <= this.nq; i1++) {
					this.yp1 = this.yh[i1];
					this.yp2 = this.yh[i1 + 1];
					for (i = 1; i <= this.n; i++)
					this.yp1[i] += this.yp2[i];
				}
			}
			pnorm = vmnorm(this.n, this.yh[1], this.ewt);

			[corflag, del, delp, told, ncf, rh, m] = this.correction(neq, y, f, corflag, pnorm, del, delp, told, ncf, rh, m, jac, _data);
			if (corflag == 0)
				break;
			if (corflag == 1) {
				rh = Math.max(rh, this.hmin / Math.abs(this.h));
				[rh, pdh] = this.scaleh(rh, pdh);
				continue;
			}
			if (corflag == 2) {
				this.kflag = -2;
				this.hold = this.h;
				this.jstart = 1;
				return;
			}
		}		/* end inner while ( corrector loop )   */
/*
   The corrector has converged.  jcur is set to 0
   to signal that the Jacobian involved may need updating later.
   The local error test is done now.
*/
		this.jcur = 0;
		if (m == 0)
			dsm = del / this.tesco[this.nq][2];
		if (m > 0)
			dsm = vmnorm(this.n, this.acor, this.ewt) / this.tesco[this.nq][2];
		if (dsm <= 1.) {
/*
   After a successful step, update the yh array.
   Decrease icount by 1, and if it is -1, consider switching methods.
   If a method switch is made, reset various parameters,
   rescale the yh array, and exit.  If there is no switch,
   consider changing h if ialth = 1.  Otherwise decrease ialth by 1.
   If ialth is then 1 and nq < maxord, then acor is saved for
   use in a possible order increase on the next step.
   If a change in h is considered, an increase or decrease in order
   by one is considered also.  A change in h is made only if it is by
   a factor of at least 1.1.  If not, ialth is set to 3 to prevent
   testing for that many steps.
*/
			this.kflag = 0;
			this.nst++;
			this.hu = this.h;
			this.nqu = this.nq;
			this.mused = this.meth;
			for (j = 1; j <= this.l; j++) {
				this.yp1 = this.yh[j];
				r = this.el[j];
				for (i = 1; i <= this.n; i++)
				this.yp1[i] += r * this.acor[i];
			}
			this.icount--;
			if (this.icount < 0) {
				[pdh,rh] = this.methodswitch(dsm, pnorm, pdh, rh);
				if (this.meth != this.mused) {
					rh = Math.max(rh, this.hmin / Math.abs(this.h));
					[rh, pdh] = this.scaleh(rh, pdh);
					this.rmax = 10.;
					this.endstoda();
					break;
				}
			}
/*
   No method switch is being made.  Do the usual step/order selection.
*/
			this.ialth--;
			if (this.ialth == 0) {
				rhup = 0.;
				if (this.l != this.lmax) {
					this.yp1 = this.yh[this.lmax];
					for (i = 1; i <= this.n; i++) this.savf[i] = this.acor[i] - this.yp1[i];
					dup = vmnorm(this.n, this.savf, this.ewt) / this.tesco[this.nq][3];
					exup = 1. / (this.l + 1);
					rhup = 1. / (1.4 * Math.pow(dup, exup) + 0.0000014);
				}
				[rhup, pdh, rh, orderflag] = this.orderswitch(rhup, dsm, pdh, rh, orderflag);
/*
   No change in h or nq.
*/
				if (orderflag == 0) {
					this.endstoda();
					break;
				}
/*
   h is changed, but not nq.
*/
				if (orderflag == 1) {
					rh = Math.max(rh, this.hmin / Math.abs(this.h));
					[rh, pdh] = this.scaleh(rh, pdh);
					this.rmax = 10.;
					this.endstoda();
					break;
				}
/*
   both nq and h are changed.
*/
				if (orderflag == 2) {
					this.resetcoeff();
					rh = Math.max(rh, this.hmin / Math.abs(this.h));
					[rh, pdh] = this.scaleh(rh, pdh);
					this.rmax = 10.;
					this.endstoda();
					break;
				}
			}	/* end if ( ialth == 0 )   */
			if (this.ialth > 1 || this.l == this.lmax) {
				this.endstoda();
				break;
			}
			this.yp1 = this.yh[this.lmax];
			for (i = 1; i <= this.n; i++)
			this.yp1[i] = this.acor[i];
			this.endstoda();
			break;
		}
		/* end if ( dsm <= 1. )   */
		/*
		   The error test failed.  kflag keeps track of multiple failures.
		   Restore tn and the yh array to their previous values, and prepare
		   to try the step again.  Compute the optimum step size for this or
		   one lower.  After 2 or more failures, h is forced to decrease
		   by a factor of 0.2 or less.
		 */ 
		else {
			this.kflag--;
			this.tn = told;
			for (j = this.nq; j >= 1; j--)
				for (i1 = j; i1 <= this.nq; i1++) {
					this.yp1 = this.yh[i1];
					this.yp2 = this.yh[i1 + 1];
					for (i = 1; i <= this.n; i++)
					this.yp1[i] -= this.yp2[i];
				}
				this.rmax = 2.;
			if (Math.abs(this.h) <= this.hmin * 1.00001) {
				this.kflag = -1;
				this.hold = this.h;
				this.jstart = 1;
				break;
			}
			if (this.kflag > -3) {
				rhup = 0.;
				[rhup, pdh, rh, orderflag] = this.orderswitch(rhup, dsm, pdh, rh, orderflag);
				if (orderflag == 1 || orderflag == 0) {
					if (orderflag == 0)
						rh = Math.min(rh, 0.2);
					rh = Math.max(rh, this.hmin / Math.abs(this.h));
					[rh, pdh] = this.scaleh(rh, pdh);
				}
				if (orderflag == 2) {
					this.resetcoeff();
					rh = Math.max(rh, this.hmin / Math.abs(this.h));
					[rh, pdh] = this.scaleh(rh, pdh);
				}
				continue;
			}
			/* if ( kflag > -3 )   */
			/*
			   Control reaches this section if 3 or more failures have occurred.
			   If 10 failures have occurred, exit with kflag = -1.
			   It is assumed that the derivatives that have accumulated in the
			   yh array have errors of the wrong order.  Hence the first
			   derivative is recomputed, and the order is set to 1.  Then
			   h is reduced by a factor of 10, and the step is retried,
			   until it succeeds or h reaches hmin.
			 */ 
			else {
				if (this.kflag == -10) {
					this.kflag = -1;
					this.hold = this.h;
					this.jstart = 1;
					break;
				} else {
					rh = 0.1;
					rh = Math.max(this.hmin / Math.abs(this.h), rh);
					this.h *= rh;
					this.yp1 = this.yh[1];
					for (i = 1; i <= this.n; i++)
						y[i] = this.yp1[i];
					// (*f) (this.tn, y + 1, this.savf + 1, _data); TODO
					f(this.tn, y, this.savf, _data)
					this.nfe++;
					this.yp1 = this.yh[2];
					for (i = 1; i <= this.n; i++)
					this.yp1[i] = this.h * this.savf[i];
					this.ipup = this.miter;
					this.ialth = 5;
					if (this.nq == 1)
						continue;
					this.nq = 1;
					this.l = 2;
					this.resetcoeff();
					continue;
				}
			}	/* end else -- kflag <= -3 */
		}		/* end error failure handling   */
	}			/* end outer while   */

}				/* end stoda   */

intdy(t: number, k: number, dky: number[], iflag: number) : number

/*
   Intdy computes interpolated values of the k-th derivative of the
   dependent variable vector y, and stores it in dky.  This routine
   is called within the package with k = 0 and *t = tout, but may
   also be called by the user for any k up to the current order.
   ( See detailed instructions in the usage documentation. )
   The computed values in dky are gotten by interpolation using the
   Nordsieck history array yh.  This array corresponds uniquely to a
   vector-valued polynomial of degree nqcur or less, and dky is set
   to the k-th derivative of this polynomial at t.
   The formula for dky is
             q
   dky[i] = sum c[k][j] * ( t - tn )^(j-k) * h^(-j) * yh[j+1][i]
            j=k
   where c[k][j] = j*(j-1)*...*(j-k+1), q = nqcur, tn = tcur, h = hcur.
   The quantities nq = nqcur, l = nq+1, n = neq, tn, and h are declared
   static globally.  The above sum is done in reverse order.
   *iflag is returned negative if either k or t is out of bounds.
*/

{
	let i: number
	let ic: number
	let j: number
	let jj: number
	let jp1: number;
	let c: number
	let r: number
	let s: number
	let tp: number;

	iflag = 0;
	if (k < 0 || k > this.nq) {
		console.error("[intdy] k = %d illegal\n", k);
		iflag = -1;
		return iflag;
	}
	tp = this.tn - this.hu - 100. * this.ETA * (this.tn + this.hu);
	if ((t - tp) * (t - this.tn) > 0.) {
		console.error(`intdy -- t = ${t} illegal. t not in interval tcur - hu to tcur\n`);
		iflag = -2;
		return iflag;
	}
	s = (t - this.tn) / this.h;
	ic = 1;
	for (jj = this.l - k; jj <= this.nq; jj++)
		ic *= jj;
	c = ic;
	this.yp1 = this.yh[this.l];
	for (i = 1; i <= this.n; i++)
		dky[i] = c * this.yp1[i];
	for (j = this.nq - 1; j >= k; j--) {
		jp1 = j + 1;
		ic = 1;
		for (jj = jp1 - k; jj <= j; jj++)
			ic *= jj;
		c = ic;
		this.yp1 = this.yh[jp1];
		for (i = 1; i <= this.n; i++)
			dky[i] = c * this.yp1[i] + s * dky[i];
	}
	if (k == 0)
		return iflag;
	r = Math.pow(this.h, (-k));
	for (i = 1; i <= this.n; i++)
		dky[i] *= r;
	return iflag
}				/* end intdy   */

scaleh(rh, pdh) : [number,number]//{rh: number, pdh: number}
{
	let r: number;
	let j: number
	let i: number;
/*
   If h is being changed, the h ratio rh is checked against rmax, hmin,
   and hmxi, and the yh array is rescaled.  ialth is set to l = nq + 1
   to prevent a change of h for that many steps, unless forced by a
   convergence or error test failure.
*/
	rh = Math.min(rh, this.rmax);
	rh = rh / Math.max(1., Math.abs(this.h) * this.hmxi * rh);
/*
   If meth = 1, also restrict the new step size by the stability region.
   If this reduces h, set irflag to 1 so that if there are roundoff
   problems later, we can assume that is the cause of the trouble.
*/
	if (this.meth == 1) {
		this.irflag = 0;
		pdh = Math.max(Math.abs(this.h) * this.pdlast, 0.000001);
		if ((rh * pdh * 1.00001) >= this.sm1[this.nq]) {
			rh = this.sm1[this.nq] / pdh;
			this.irflag = 1;
		}
	}
	r = 1.;
	for (j = 2; j <= this.l; j++) {
		r *= rh;
		this.yp1 = this.yh[j];
		for (i = 1; i <= this.n; i++) this.yp1[i] *= r;
	}
	this.h *= rh;
	this.rc *= rh;
	this.ialth = this.l;
	return [rh,pdh]//{rh, pdh}

} /* end scaleh   */

/*
-----------------------------------------------------------------------
 prja is called by stoda to compute and process the matrix
 p = i - h*el(1)*j , where j is an approximation to the jacobian.
 here j is computed by the user-supplied routine jac if
 miter = 1 or 4 or by finite differencing if miter = 2 or 5.
 j, scaled by -h*el(1), is stored in wm.  then the norm of j (the
 matrix norm consistent with the weighted max-norm on vectors given
 by vmnorm) is computed, and j is overwritten by p.  p is then
 subjected to lu decomposition in preparation for later solution
 of linear systems with p as coefficient matrix. this is done
 by sgefa if miter = 1 or 2, and by sgbfa if miter = 4 or 5.

 in addition to variables described previously, communication
 with prja uses the following..
 y     = array containing predicted values on entry.
 ftem  = work array of length n (acor in stoda).
 savf  = array containing f evaluated at predicted y.
 wm    = real work space for matrices.  on output it contains the
         lu decomposition of p.
         storage of matrix elements starts at wm(3).
         wm also contains the following matrix-related data..
         wm(1) = sqrt(uround), used in numerical jacobian increments.
 iwm   = integer work space containing pivot information, starting at
         iwm(21).   iwm also contains the band parameters
         ml = iwm(1) and mu = iwm(2) if miter is 4 or 5.
 el0   = el(1) (input).
 pdnorm= norm of jacobian matrix. (output).
 ierpj = output error flag,  = 0 if no trouble, .gt. 0 if
         p matrix found to be singular.
 jcur  = output flag = 1 to indicate that the jacobian matrix
         (or approximation) is now current.
 this routine also uses the common variables el0, h, tn, uround,
 miter, n, nfe, and nje.
-----------------------------------------------------------------------
*/
prja(neq: number, y: number[], f: any, jac: any, _data: number[])
{
	let i: number
	let ier: number
	let j: number;
	let fac: number
	let hl0: number
	let r: number
	let r0: number
	let yj: number
	let mband: number
	let meband: number

	this.nje++;
	this.ierpj = 0;
	this.jcur = 1;
	hl0 = this.h * this.el0;

	if (this.miter == 1) {
		// if miter = 1, call jac and multiply by scalar. -----------------------
		for (i = 1; i <= this.n; i++) {
			for(j = 1; j <= this.n; j++) {
				this.wm[i][j] = 0
			}
		}
		jac(this.tn, y, 0, 0, this.wm, _data)
		for (i = 1; i <= this.n; i++) {
			for(j = 1; j <= this.n; j++) {
				this.wm[i][j] *= -hl0
			}
		}

	} else if (this.miter == 2) {
		// if miter = 2, make n calls to f to approximate j. --------------------
		fac = vmnorm(this.n, this.savf, this.ewt);
		r0 = 1000. * Math.abs(this.h) * this.ETA * this.n * fac;
		if (r0 == 0.)
			r0 = 1.;
		for (j = 1; j <= this.n; j++) {
			yj = y[j];
			r = Math.max(this.sqrteta * Math.abs(yj), r0 / this.ewt[j]);
			y[j] += r;
			fac = -hl0 / r;
			f(this.tn, y, this.acor, _data)
			for (i = 1; i <= this.n; i++) this.wm[i][j] = (this.acor[i] - this.savf[i]) * fac;
			y[j] = yj;
		}
		this.nfe += this.n;
	} else if (this.miter == 3) {
		// dummy block only, since miter is never 3 in this routine. ------------
		return
	} else if (this.miter == 4) {
		// if miter = 4, call jac and multiply by scalar. -----------------------
		throw "Error in prja: We do not support miter=4";
		let ml3 = this.ml + 3
		mband = this.ml + this.mu + 1
		meband = mband + this.ml
		for (i = 1; i <= this.n; i++) {
			for(j = 1; j <= this.n; j++) {
				this.wm[i][j] = 0
			}
		}
		jac(this.tn, y, this.ml, this.mu, this.wm, _data)
		// pointerMath fixes ml3 = ml + 3 (1 + 2offset), so we shift it by ml
		for (i = meband; i >= 1; i--) {
			for (let j = 1; j <= this.n; j++) {
				this.wm[i+this.ml][j] = this.wm[i][j] * -hl0
				this.wm[i][j] = 0
			}
		}
	} else if(this.miter == 5) {
		// if miter = 5, make mband calls to f to approximate j. ----------------
		throw "Error in prja: We do not currently support miter=5";
		let yi: number
		let jj: number
		let yjj: number
		let i1: number
		let i2: number
		let ii: number
		mband = this.ml + this.mu + 1
     	let mba = Math.min(mband,this.n)
     	meband = mband + this.ml
     	let meb1 = meband - 1
     	fac = vmnorm (this.n, this.savf, this.ewt)
     	r0 = 1000.0e0*Math.abs(this.h)*this.ETA*this.n*fac
		if (r0 == 0.0e0) r0 = 1.0e0
		for (j = 1; j <= mba; j++) {
			for (i = j; i <= this.n; i = i + mband) {
				yi = y[i]
				r = Math.max(this.sqrteta*Math.abs(yi),r0/this.ewt[i])
				y[i] = y[i] + r
			}
        	f (this.tn, y, this.acor, _data)
			for (jj = j; jj <= this.n; jj = jj + mband) {
				y[jj] = this.yh[1][jj]
				yjj = y[jj]
				r = Math.max(this.sqrteta*Math.abs(yjj),r0/this.ewt[jj])
				fac = -hl0/r
				i1 = Math.max(jj-this.mu,1)
				i2 = Math.min(jj+this.ml,this.n)
				ii = jj*meb1 - this.ml// + 2
				for (i = i1; i <= i2; i++) {
					this.wm[i][jj] = (this.acor[i] - this.savf[i])*fac
				}
			}
		}
      	this.nfe = this.nfe + mba
		return
	}

	// compute norm of jacobian. --------------------------------------------
	this.pdnorm = fnorm(this.n, this.wm, this.ewt) / Math.abs(hl0);

	if (this.miter < 3) {
		// add identity matrix. -------------------------------------------------
		for (i = 1; i <= this.n; i++) this.wm[i][i] += 1.;

		// do lu decomposition on p. --------------------------------------------
		ier = sgefa(this.wm, this.n, this.ipvt, ier);
	} else {
		// add identity matrix. -------------------------------------------------
		for (i = 1; i <= this.n; i++) this.wm[mband][i] += 1.;
		
		// do lu decomposition on p. --------------------------------------------
		ier = sgbfa(this.wm, meband, this.n, this.ml, this.mu, this.ipvt, ier)
	}
	if (ier != 0) this.ierpj = 1;
	return;
} /* end prja   */

correction(neq: number, y: number[], f: any, corflag: number, pnorm: number, del: number, delp: number, told: number,
					   ncf: number, rh: number, m: number, jac: any, _data: number[]) : [number,number,number,number,number,number,number]//{corflag, del: number, delp: number, m, rh}
/*
   *corflag = 0 : corrector converged,
              1 : step size to be reduced, redo prediction,
              2 : corrector cannot converge, failure flag.
*/

{
	let i: number;
	let rm: number
	let rate: number
	let dcon: number;

/*
   Up to maxcor corrector iterations are taken.  A convergence test is
   made on the r.m.s. norm of each correction, weighted by the error
   weight vector ewt.  The sum of the corrections is accumulated in the
   vector acor[i].  The yh array is not altered in the corrector loop.
*/

	m = 0;
	corflag = 0;
	rate = 0.;
	del = 0.;
	this.yp1 = this.yh[1];
	for (i = 1; i <= this.n; i++){
		y[i] = this.yp1[i];
	}
	// (*f) (tn, y + 1, this.savf + 1, _data); TODO
	f(this.tn, y, this.savf, _data)
	this.nfe++;
/*
   If indicated, the matrix P = I - h * el[1] * J is reevaluated and
   preprocessed before starting the corrector iteration.  ipup is set
   to 0 as an indicator that this has been done.
*/
	while (1) {
		if (m == 0) {
			if (this.ipup > 0) {
				this.prja(neq, y, f, jac, _data);
				this.ipup = 0;
				this.rc = 1.;
				this.nslp = this.nst;
				this.crate = 0.7;
				if (this.ierpj != 0) {
					[rh, ncf, corflag] = this.corfailure(told, rh, ncf, corflag);
					return [corflag, del, delp, told, ncf, rh, m];
				}
			}
			for (i = 1; i <= this.n; i++)
			this.acor[i] = 0.;
		}		/* end if ( m == 0 )   */
		if (this.miter == 0) {
/*
   In case of functional iteration, update y directly from
   the result of the last function evaluation.
*/
			this.yp1 = this.yh[2];
			for (i = 1; i <= this.n; i++) {
				this.savf[i] = this.h * this.savf[i] - this.yp1[i];
				y[i] = this.savf[i] - this.acor[i];
			}
			del = vmnorm(this.n, y, this.ewt);
			this.yp1 = this.yh[1];
			for (i = 1; i <= this.n; i++) {
				y[i] = this.yp1[i] + this.el[1] * this.savf[i];
				this.acor[i] = this.savf[i];
			}
		}
		/* end functional iteration   */
		/*
		   In the case of the chord method, compute the corrector error,
		   and solve the linear system with that as right-hand side and
		   P as coefficient matrix.
		 */ 
		else {
			this.yp1 = this.yh[2];
			for (i = 1; i <= this.n; i++){ y[i] = this.h * this.savf[i] - (this.yp1[i] + this.acor[i]); }
			this.solsy(y);
			del = vmnorm(this.n, y, this.ewt);
			this.yp1 = this.yh[1];
			for (i = 1; i <= this.n; i++) {
				this.acor[i] += y[i];
				y[i] = this.yp1[i] + this.el[1] * this.acor[i];
			}
		}		/* end chord method   */
/*
   Test for convergence.  If m > 0, an estimate of the convergence
   rate constant is stored in crate, and this is used in the test.
   We first check for a change of iterates that is the size of
   roundoff error.  If this occurs, the iteration has converged, and a
   new rate estimate is not formed.
   In all other cases, force at least two iterations to estimate a
   local Lipschitz constant estimate for Adams method.
   On convergence, form pdest = local maximum Lipschitz constant
   estimate.  pdlast is the most recent nonzero estimate.
*/
		if (del <= 100. * pnorm * this.ETA)
			break;
		if (m != 0 || this.meth != 1) {
			if (m != 0) {
				rm = 1024.0;
				if (del <= (1024. * delp))
					rm = del / delp;
				rate = Math.max(rate, rm);
				this.crate = Math.max(0.2 * this.crate, rm);
			}
			dcon = del * Math.min(1., 1.5 * this.crate) / (this.tesco[this.nq][2] * this.conit);
			if (dcon <= 1.) {
				this.pdest = Math.max(this.pdest, rate / Math.abs(this.h * this.el[1]));
				if (this.pdest != 0.) this.pdlast = this.pdest;
				break;
			}
		}
/*
   The corrector iteration failed to converge.
   If miter != 0 and the Jacobian is out of date, prja is called for
   the next try.   Otherwise the yh array is retracted to its values
   before prediction, and h is reduced, if possible.  If h cannot be
   reduced or mxncf failures have occured, exit with corflag = 2.
*/
		(m)++;
		if (m == this.maxcor || (m >= 2 && del > 2. * delp)) {
			if (this.miter == 0 || this.jcur == 1) {
				[rh, ncf, corflag] = this.corfailure(told, rh, ncf, corflag);
				return [corflag, del, delp, told, ncf, rh, m];
			}
			this.ipup = this.miter;
/*
   Restart corrector if Jacobian is recomputed.
*/
			m = 0;
			rate = 0.;
			del = 0.;
			this.yp1 = this.yh[1];
			for (i = 1; i <= this.n; i++) y[i] = this.yp1[i];
			// (*f) (tn, y + 1, this.savf + 1, _data); TODO
			f(this.tn,y,this.savf, _data)
			this.nfe++;
		}
/*
   Iterate corrector.
*/
		else {
			delp = del;
			// (*f) (tn, y + 1, this.savf + 1, _data); TODO
			f(this.tn, y, this.savf, _data)
			this.nfe++;
		}
	}			/* end while   */
	return [corflag, del, delp, told, ncf, rh, m];
}				/* end correction   */

corfailure(told: number, rh: number, ncf: number, corflag: number) : [number,number,number]//{rh, ncf, corflag}
{
	let j: number
	let i1: number
	let i: number;

	ncf++;
	this.rmax = 2.;
	this.tn = told;
	for (j = this.nq; j >= 1; j--)
		for (i1 = j; i1 <= this.nq; i1++) {
			this.yp1 = this.yh[i1];
			this.yp2 = this.yh[i1 + 1];
			for (i = 1; i <= this.n; i++) this.yp1[i] -= this.yp2[i];
		}
	if (Math.abs(this.h) <= this.hmin * 1.00001 || ncf == this.mxncf) {
		corflag = 2;
		return [rh, ncf, corflag]
	}
	corflag = 1;
	rh = 0.25;
	this.ipup = this.miter;
	return [rh, ncf, corflag]
}

/*
-----------------------------------------------------------------------
 this routine manages the solution of the linear system arising from
 a chord iteration.  it is called if miter .ne. 0.
 if miter is 1 or 2, it calls sgesl to accomplish this.
 if miter = 3 it updates the coefficient h*el0 in the diagonal
 matrix, and then computes the solution.
 if miter is 4 or 5, it calls sgbsl.
 communication with solsy uses the following variables..
 wm    = real work space containing the inverse diagonal matrix if
         miter = 3 and the lu decomposition of the matrix otherwise.
         storage of matrix elements starts at wm(3).
         wm also contains the following matrix-related data..
         wm(1) = sqrt(uround) (not used here),
         wm(2) = hl0, the previous value of h*el0, used if miter = 3.
 iwm   = integer work space containing pivot information, starting at
         iwm(21), if miter is 1, 2, 4, or 5.  iwm also contains band
         parameters ml = iwm(1) and mu = iwm(2) if miter is 4 or 5.
 x     = the right-hand side vector on input, and the solution vector
         on output, of length n.
 tem   = vector of work space of length n, not used in this version.
 iersl = output flag (in common).  iersl = 0 if no trouble occurred.
         iersl = 1 if a singular matrix arose with miter = 3.
 this routine also uses the common variables el0, h, miter, and n.
-----------------------------------------------------------------------
*/
solsy(y: number[]) {
	this.iersl = 0;
	if (this.miter == 1 || this.miter == 2) {
		sgesl(this.wm, this.n, this.ipvt, y, 0);
	} else if(this.miter == 3) {
		throw `solsy with miter = ${this.miter}`
	} else if(this.miter == 4 || this.miter == 5) {
		let meband = 2*this.ml + this.mu + 1
		sgbsl(this.wm, meband, this.n, this.ml, this.mu, this.ipvt, y, 0)
	}
	return;

}

methodswitch(dsm: number, pnorm: number, pdh: number, rh: number) : [number,number]//{pdh, rh}
{
	let lm1: number
	let lm1p1: number
	let lm2: number
	let lm2p1: number
	let nqm1: number
	let nqm2: number;
	let rh1: number
	let rh2: number
	let rh1it: number
	let exm2: number
	let dm2: number
	let exm1: number
	let dm1: number
	let alpha: number
	let exsm: number;

/*
   We are current using an Adams method.  Consider switching to bdf.
   If the current order is greater than 5, assume the problem is
   not stiff, and skip this section.
   If the Lipschitz constant and error estimate are not polluted
   by roundoff, perform the usual test.
   Otherwise, switch to the bdf methods if the last step was
   restricted to insure stability ( irflag = 1 ), and stay with Adams
   method if not.  When switching to bdf with polluted error estimates,
   in the absence of other information, double the step size.
   When the estimates are ok, we make the usual test by computing
   the step size we could have (ideally) used on this step,
   with the current (Adams) method, and also that for the bdf.
   If nq > mxords, we consider changing to order mxords on switching.
   Compare the two step sizes to decide whether to switch.
   The step size advantage must be at least ratio = 5 to switch.
*/
	if (this.meth == 1) {
		if (this.nq > 5)
			return [pdh,rh];
		if (dsm <= (100. * pnorm * this.ETA) || this.pdest == 0.) {
			if (this.irflag == 0)
				return [pdh,rh];
			rh2 = 2.;
			nqm2 = Math.min(this.nq, this.mxords);
		} else {
			exsm = 1. / this.l;
			rh1 = 1. / (1.2 * Math.pow(dsm, exsm) + 0.0000012);
			rh1it = 2. * rh1;
			pdh = this.pdlast * Math.abs(this.h);
			if ((pdh * rh1) > 0.00001)
				rh1it = this.sm1[this.nq] / pdh;
			rh1 = Math.min(rh1, rh1it);
			if (this.nq > this.mxords) {
				nqm2 = this.mxords;
				lm2 = this.mxords + 1;
				exm2 = 1. / lm2;
				lm2p1 = lm2 + 1;
				dm2 = vmnorm(this.n, this.yh[lm2p1], this.ewt) / this.cm2[this.mxords];
				rh2 = 1. / (1.2 * Math.pow(dm2, exm2) + 0.0000012);
			} else {
				dm2 = dsm * (this.cm1[this.nq] / this.cm2[this.nq]);
				rh2 = 1. / (1.2 * Math.pow(dm2, exsm) + 0.0000012);
				nqm2 = this.nq;
			}
			if (rh2 < this.ratio * rh1)
				return [pdh,rh];
		}
/*
   The method switch test passed.  Reset relevant quantities for bdf.
*/
		rh = rh2;
		this.icount = 20;
		this.meth = 2;
		this.miter = this.jtyp;
		this.pdlast = 0.;
		this.nq = nqm2;
		this.l = this.nq + 1;
		return [pdh,rh];
	}			/* end if ( meth == 1 )   */
	/*
	   We are currently using a bdf method, considering switching to Adams.
	   Compute the step size we could have (ideally) used on this step,
	   with the current (bdf) method, and also that for the Adams.
	   If nq > mxordn, we consider changing to order mxordn on switching.
	   Compare the two step sizes to decide whether to switch.
	   The step size advantage must be at least 5/ratio = 1 to switch.
	   If the step size for Adams would be so small as to cause
	   roundoff pollution, we stay with bdf.
	*/
	exsm = 1. / this.l;
	if (this.mxordn < this.nq) {
		nqm1 = this.mxordn;
		lm1 = this.mxordn + 1;
		exm1 = 1. / lm1;
		lm1p1 = lm1 + 1;
		dm1 = vmnorm(this.n, this.yh[lm1p1], this.ewt) / this.cm1[this.mxordn];
		rh1 = 1. / (1.2 * Math.pow(dm1, exm1) + 0.0000012);
	} else {
		dm1 = dsm * (this.cm2[this.nq] / this.cm1[this.nq]);
		rh1 = 1. / (1.2 * Math.pow(dm1, exsm) + 0.0000012);
		nqm1 = this.nq;
		exm1 = exsm;
	}
	rh1it = 2. * rh1;
	pdh = this.pdnorm * Math.abs(this.h);
	if ((pdh * rh1) > 0.00001)
		rh1it = this.sm1[nqm1] / pdh;
	rh1 = Math.min(rh1, rh1it);
	rh2 = 1. / (1.2 * Math.pow(dsm, exsm) + 0.0000012);
	if ((rh1 * this.ratio) < (5. * rh2))
		return [pdh,rh];
	alpha = Math.max(0.001, rh1);
	dm1 *= Math.pow(alpha, exm1);
	if (dm1 <= 1000. * this.ETA * pnorm)
		return [pdh,rh];
/*
   The switch test passed.  Reset relevant quantities for Adams.
*/
	rh = rh1;
	this.icount = 20;
	this.meth = 1;
	this.miter = 0;
	this.pdlast = 0.;
	this.nq = nqm1;
	this.l = this.nq + 1;
	return [pdh,rh]
}				/* end methodswitch   */


/*
   This routine returns from stoda to lsoda.  Hence freevectors() is
   not executed.
*/
endstoda()
{
	let r: number;
	let i: number;

	r = 1. / this.tesco[this.nqu][2];
	for (i = 1; i <= this.n; i++) {this.acor[i] *= r;}
	this.hold = this.h;
	this.jstart = 1;

}

orderswitch(rhup: number, dsm: number, pdh: number, rh: number, orderflag: number) : [number,number,number,number]//{ rhup, pdh, rh, orderflag }

/*
   Regardless of the success or failure of the step, factors
   rhdn, rhsm, and rhup are computed, by which h could be multiplied
   at order nq - 1, order nq, or order nq + 1, respectively.
   In the case of a failure, rhup = 0. to avoid an order increase.
   The largest of these is determined and the new order chosen
   accordingly.  If the order is to be increased, we compute one
   additional scaled derivative.
   orderflag = 0  : no change in h or nq,
               1  : change in h but not nq,
               2  : change in both h and nq.
*/

{
	let newq: number
	let i: number;
	let exsm: number
	let rhdn: number
	let rhsm: number
	let ddn: number
	let exdn: number
	let r: number;

	orderflag = 0;

	exsm = 1. / this.l;
	rhsm = 1. / (1.2 * Math.pow(dsm, exsm) + 0.0000012);

	rhdn = 0.;
	if (this.nq != 1) {
		ddn = vmnorm(this.n, this.yh[this.l], this.ewt) / this.tesco[this.nq][1];
		exdn = 1. / this.nq;
		rhdn = 1. / (1.3 * Math.pow(ddn, exdn) + 0.0000013);
	}
/*
   If meth = 1, limit rh accordinfg to the stability region also.
*/
	if (this.meth == 1) {
		pdh = Math.max(Math.abs(this.h) * this.pdlast, 0.000001);
		if (this.l < this.lmax) { rhup = Math.min(rhup, this.sm1[this.l] / pdh); }
		rhsm = Math.min(rhsm, this.sm1[this.nq] / pdh);
		if (this.nq > 1) { rhdn = Math.min(rhdn, this.sm1[this.nq - 1] / pdh); }
		this.pdest = 0.;
	}
	if (rhsm >= rhup) {
		if (rhsm >= rhdn) {
			newq = this.nq;
			rh = rhsm;
		} else {
			newq = this.nq - 1;
			rh = rhdn;
			if (this.kflag < 0 && rh > 1.)
				rh = 1.;
		}
	} else {
		if (rhup <= rhdn) {
			newq = this.nq - 1;
			rh = rhdn;
			if (this.kflag < 0 && rh > 1.)
				rh = 1.;
		} else {
			rh = rhup;
			if (rh >= 1.1) {
				r = this.el[this.l] / this.l;
				this.nq = this.l;
				this.l = this.nq + 1;
				this.yp1 = this.yh[this.l];
				for (i = 1; i <= this.n; i++) { this.yp1[i] = this.acor[i] * r; }
				orderflag = 2;
				return [rhup, pdh,rh,orderflag]
			} else {
				this.ialth = 3;
				return [rhup, pdh,rh,orderflag]
			}
		}
	}
/*
   If meth = 1 and h is restricted by stability, bypass 10 percent test.
*/
	if (this.meth == 1) {
		if ((rh * pdh * 1.00001) < this.sm1[newq])
			if (this.kflag == 0 && rh < 1.1) {
				this.ialth = 3;
				return [rhup, pdh,rh,orderflag]
			}
	} else {
		if (this.kflag == 0 && rh < 1.1) {
			this.ialth = 3;
			return [rhup, pdh,rh,orderflag]
		}
	}
	if (this.kflag <= -2)
		rh = Math.min(rh, 0.2);
/*
   If there is a change of order, reset nq, l, and the coefficients.
   In any case h is reset according to rh and the yh array is rescaled.
   Then exit or redo the step.
*/
	if (newq == this.nq) {
		orderflag = 1;
		return [rhup, pdh,rh,orderflag]
	}
	this.nq = newq;
	this.l = this.nq + 1;
	orderflag = 2;
	return [rhup, pdh,rh,orderflag]
}				/* end orderswitch   */

/*
-----------------------------------------------------------------------
 the el vector and related constants are reset
 whenever the order nq is changed, or at the start of the problem.
-----------------------------------------------------------------------
*/
resetcoeff() {
	let i: number;
	let ep1: number[];

	ep1 = this.elco[this.nq];
	for (i = 1; i <= this.l; i++) { this.el[i] = ep1[i]; }
	this.rc = this.rc * this.el[1] / this.el0;
	this.el0 = this.el[1];
	this.conit = 0.5 / (this.nq + 2);

}

}