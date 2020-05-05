/*
    This is a typescript version of the LSODAR library. The original is written
    in fortran which contains pointer logic and needed to be rewritten since
    javascript and typescript do not have pointers.

    Most of the original author's notes are kept including directions on how to
    use the library.

    - Pawel Spychala <paulspychala@gmail.com>
*/

import { LSODA, jac_func, lsoda_func } from "./lsoda";
import {ewset} from "./ewset";
import {scopy} from "./scopy";
import {vmnorm} from "./vmnorm";

export type lsodar_root_func = (t: number, y: number[], gout: number[], data: any) => any

export class LSODAR extends LSODA {
    nge: number      // the number of g evaluation for the problem so far
    t0: number       // value of t at one endpoint of interval of interest.
    irfnd: number    // input flag showing whether the last step taken had a root. 
    tlast: number    // last value of t returned by solver
    itaskc: number   // copy of itask
    ngc: number      // copy of ng
    toutc: number    // copy of tout

    last: number     // for roots.ts
    imax: number
    alpha: number
    x2: number

    g0: number[]
    g1: number[]
    gx: number[]
    
/*
-----------------------------------------------------------------------
 this is the march 30, 1987 version of
 this.. Livermore Solver for Ordinary Differential Equations, with
          Automatic method switching for stiff and nonstiff problems,
          and with Root-finding.

 This typescript version is in double precision. (all numbers are double)

 lsodar solves the initial value problem for stiff or nonstiff
 systems of first order ode-s,
     dy/dt = f(t,y) ,  or, in component form,
     dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(neq)) (i = 1,...,neq).
 at the same time, it locates the roots of any of a set of functions
     g(i) = g(i,t,y(1),...,y(neq))  (i = 1,...,ng).

 This a variant version of the lsode package.  it differs from lsode
 in two ways..
 (a) it switches automatically between stiff and nonstiff methods.
 this means that the user does not have to determine whether the
 problem is stiff or not, and the solver will automatically choose the
 appropriate method.  it always starts with the nonstiff method.
 (b) it finds the root of at least one of a set of constraint
 functions g(i) of the independent and dependent variables.
 it finds only those roots for which some g(i), as a function
 of t, changes sign in the interval of integration.
 It then returns the solution at the root, if that occurs
 sooner than the specified stop condition, and otherwise returns
 the solution according the specified stop condition.

 authors..
                Linda R. Petzold and Alan C. Hindmarsh,
                computing and mathematics research division, l-316
                Lawrence Livermore National Laboratory
                Livermore, CA 94550.

 references..
 1.  Alan C. Hindmarsh,  Odepack, a systematized collection of ode
     solvers, in scientific computing, R. S. Stepleman et al. (eds.),
     north-holland, Amsterdam, 1983, pp. 55-64.
 2.  Linda R. Petzold, Automatic selection of methods for solving
     stiff and nonstiff systems of ordinary differential equations,
     siam j. sci. stat. comput. 4 (1983), pp. 136-148.
 3.  Kathie L. Hiebert and Lawrence F. Shampine, Implicitly defined
     output points for solutions of ode-s, sandia report sand80-0180,
     February, 1980.
-----------------------------------------------------------------------
 Summary of Usage.

 Communication between the user and the lsodar package, for normal
 situations, is summarized here.  This summary describes only a subset
 of the full set of options available.  See the full description for
 details, including alternative treatment of the jacobian matrix,
 optional inputs and outputs, nonstandard options, and
 instructions for special situations.  see also the example
 problem (with program and output) following this summary.

 a. first provide a subroutine of the form..
               subroutine f (neq, t, y, ydot)
               dimension y(neq), ydot(neq)
 which supplies the vector function f by loading ydot(i) with f(i).

 b. provide a subroutine of the form..
               subroutine g (neq, t, y, ng, gout)
               dimension y(neq), gout(ng)
 which supplies the vector function g by loading gout(i) with
 g(i), the i-th constraint function whose root is sought.

 c. write a main program which calls subroutine lsodar once for
 each point at which answers are desired.  this should also provide
 for possible use of logical unit 6 for output of error messages by
 this.  on the first call to lsodar, supply arguments as follows..
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
             22 + neq * max(16, neq + 9) + 3*ng.
          see also paragraph f below.
 lrw    = declared length of rwork (in user-s dimension).
 iwork  = integer work array of length at least  20 + neq.
 liw    = declared length of iwork (in user-s dimension).
 jac    = name of subroutine for jacobian matrix.
          use a dummy name.  see also paragraph f below.
 jt     = jacobian type indicator.  set jt = 2.
          see also paragraph f below.
 g      = name of subroutine for constraint functions, whose
          roots are desired during the integration.
          this name must be declared external in calling program.
 ng     = number of constraint functions g(i).  if there are none,
          set ng = 0, and pass a dummy name for g.
 jroot  = integer array of length ng for output of root information.
          see next paragraph.
 note that the main program must declare arrays y, rwork, iwork,
 jroot, and possibly atol.

 d. the output from the first call (or any call) is..
      y = array of computed values of y(t) vector.
      t = corresponding value of independent variable.  this is
          tout if istate = 2, or the root location if istate = 3,
          or the farthest point reached if lsodar was unsuccessful.
 istate = 2 or 3  if lsodar was successful, negative otherwise.
           2 means no root was found, and tout was reached as desired.
           3 means a root was found prior to reaching tout.
          -1 means excess work done on this call (perhaps wrong jt).
          -2 means excess accuracy requested (tolerances too small).
          -3 means illegal input detected (see printed message).
          -4 means repeated error test failures (check all inputs).
          -5 means repeated convergence failures (perhaps bad jacobian
             supplied or wrong choice of jt or tolerances).
          -6 means error weight became zero during problem. (solution
             component i vanished, and atol or atol(i) = 0.)
          -7 means work space insufficient to finish (see messages).
 jroot  = array showing roots found if istate = 3 on return.
          jroot(i) = 1 if g(i) has a root at t, or 0 otherwise.

 e. to continue the integration after a successful return, proceed
 as follows..
  (a) if istate = 2 on return, reset tout and call lsodar again.
  (b) if istate = 3 on return, reset istate to 2 and call lsodar again.
 in either case, no other parameters need be reset.

 f. note.. if and when lsodar regards the problem as stiff, and
 switches methods accordingly, it must make use of the neq by neq
 jacobian matrix, j = df/dy.  for the sake of simplicity, the
 inputs to lsodar recommended in paragraph c above cause lsodar to
 treat j as a full matrix, and to approximate it internally by
 difference quotients.  alternatively, j can be treated as a band
 matrix (with great potential reduction in the size of the rwork
 array).  also, in either the full or banded case, the user can supply
 j in closed form, with a routine whose name is passed as the jac
 argument.  these alternatives are described in the paragraphs on
 rwork, jac, and jt in the full description of the call sequence below.
*/

    lsodar(f: lsoda_func, neq: number, y: number[], t: number, tout: number,
        itol: number, rtol: number[], atol: number[],
        itask: number, istate: number, iopt: number, jac: jac_func, jt: number,
        iwork1: number, iwork2: number, iwork5: number, iwork6: number,
        iwork7: number, iwork8: number, iwork9: number,
        rwork1: number, rwork5: number, rwork6: number, rwork7: number,
        g: lsodar_root_func, ng: number, jroot: number[],
        _data: any): [number, number]
/*
optional inputs.
The following is a list of the optional inputs provided for in the
call sequence.  (see also part ii.)  For each such input variable,
this table lists its name as used in this documentation, its
location in the call sequence, its meaning, and the default value.
The use of any of these inputs requires iopt = 1, and in that
case all of these inputs are examined.  A value of zero for any
of these optional inputs will cause the default value to be used.
Thus to use a subset of the optional inputs, simply preload
locations 5 to 10 in rwork and iwork to 0.0 and 0 respectively, and
then set those of interest to nonzero values.

name    loc         meaning and default value

ml      iwork1  these are the lower
mu      iwork2  and upper half-bandwidths of the banded jacobian, excluding
                the main diagonal. the band is defined by the matrix locations
                (i,j) with i-ml <= j <= i+mu. ml and mu must satisfy 0 <= ml,mu
                <= neq-1. these are required if jt is 4 or 5, and ignored
                otherwise. ml and mu may in fact be the band parameters for
                a matrix to which df/dy is only approximately equal
ixpr    iwork5  flag to generate extra printing at method switches.
                ipxr = 0 means no extra printing (the default).
                ipxr = 1 means print data on each switch
                t, h, and nst will be printed on the same logical unit as used
                for error messages.
mxstep  iwork6  maximum number of (internally defined) steps allowed during one
                call to the solver.the default value is 500
mxhnil  iwork7  maximum number of messages printed (per problem) warning that
                t + h = t on a step (h = step size). this must be positive to
                result in a non-default value. the default value is 10
mxordn  iwork8  the maximum order to be allowed for the nonstiff (adams) 
                method. the default value is 12. if mxordn exceeds the default
                value, it will be reduced to the default value. mxordn is
                held constant during the problem.
mxords  iwork9  the maximum order to be allowed for the stiff (bdf) method.
                the default value is 5. if mxords exceeds the default value, 
                it will be reduced to the default value. mxords is held
                constant during the problem.

tcrit   rwork1
h0      rwork5  the step size to be attempted on the first step. 
                the default value is determined by the solver
hmax    rwork6  the maximum absolute step size allowed.
                the default value is infinite
hmin    rwork7  the minimum absolute step size allowed.
                the default value is 0. (this lower bound is not enforced)
                on the final step before reaching tcrit when itask = 4 or 5)

optional outputs.

As optional additional output from lsodar, the variables listed
below are quantities related to the performance of lsodar
which are available to the user.  these are communicated by way of
the work arrays, but also have internal mnemonic names as shown.
except where stated otherwise, all of these outputs are defined
on any successful return from lsodar, and on any return with
istate = -1, -2, -4, -5, or -6.  on an illegal input return
(istate = -3), they will be unchanged from their existing values
(if any), except possibly for tolsf, lenrw, and leniw.
on any error return, outputs relevant to the error will be defined,
as noted below.

hu      rwork11 the step size in t last used (successfully).
hcur    rwork12 the step size to be attempted on the next step.
tcur    rwork13 the current value of the independent variable
                which the solver has actually reached, i.e. the
                current internal mesh point in t.  on output, tcur
                will always be at least as far as the argument
                t, but may be farther (if interpolation was done).
tolsf   rwork14 a tolerance scale factor, greater than 1.0,
                computed when a request for too much accuracy was
                detected (istate = -3 if detected at the start of
                the problem, istate = -2 otherwise).  if itol is
                left unaltered but rtol and atol are uniformly
                scaled up by a factor of tolsf for the next call,
                then the solver is deemed likely to succeed.
                (the user may also ignore tolsf and alter the
                tolerance parameters in any other way appropriate.)
tsw     rwork15 the value of t at the time of the last method switch, if any.

nge     iwork10 the number of g evaluations for the problem so far.
nst     iwork11 the number of steps taken for the problem so far.
nfe     iwork12 the number of f evaluations for the problem so far.
nje     iwork13 the number of jacobian evaluations (and of matrix
                lu decompositions) for the problem so far.
nqu     iwork14 the method order last used (successfully).
nqcur   iwork15 the order to be attempted on the next step.
imxer   iwork16 the index of the component of largest magnitude in
                the weighted local error vector ( e(i)/ewt(i) ),
                on an error return with istate = -4 or -5.
lenrw   iwork17 the length of rwork actually required, assuming
                that the length of rwork is to be fixed for the
                rest of the problem, and that switching may occur.
                this is defined on normal returns and on an illegal
                input return for insufficient storage.
leniw   iwork18 the length of iwork actually required, assuming
                that the length of iwork is to be fixed for the
                rest of the problem, and that switching may occur.
                this is defined on normal returns and on an illegal
                input return for insufficient storage.
musde   iwork19 the method indicator for the last successful step..
                1 means adams (nonstiff), 2 means bdf (stiff).
mcur    iwork20 the current method indicator..
                1 means adams (nonstiff), 2 means bdf (stiff).
                this is the method to be attempted
                on the next step.  thus it differs from mused
                only if a method switch has just been made.
*/ {
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

        // LSODAR additions
        let irfp: number
        let irt: number

        /*
        -----------------------------------------------------------------------
         block a.
         this code block is executed on every call.
         it tests istate and itask for legality and branches appropriately.
         if istate > 1 but the flag init shows that initialization has
         not yet been done, an error return occurs.
         if istate = 1 and tout = t, jump to block g and return immediately.
        -----------------------------------------------------------------------
        */

        if (istate < 1 || istate > 3) {
            console.error("[lsodar] illegal istate = %d\n", istate);
            return [t,this.terminate(istate)];
        }
        if (itask < 1 || itask > 5) {
            console.error("[lsodar] illegal itask = %d\n", itask);
            return [t,this.terminate(istate)];
        }
        if (this.init == 0 && (istate == 2 || istate == 3)) {
            console.error("[lsodar] istate > 1 but lsoda not initialized\n");
            return [t,this.terminate(istate)];
        }
        if (istate == 1) {
            this.init = 0;
            if (tout == t) {
                this.ntrep++;
                if (this.ntrep < 5) return [t,istate];
                console.error("[lsodar] repeated calls with istate = 1 and tout = t. run aborted.. apparent infinite loop\n");
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
         jt, ml, mu, and ng.
        -----------------------------------------------------------------------
        */

        if (istate == 1 || istate == 3) {
            this.ntrep = 0;
            if (neq <= 0) {
                console.error(`[lsodar] neq = ${neq} is less than 1`);
                return [t,this.terminate(istate)];
            }
            if (istate == 3 && neq > this.n) {
                console.error("[lsodar] istate = 3 and neq increased");
                return [t,this.terminate(istate)];
            }
            this.n = neq;
            if (itol < 1 || itol > 4) {
                console.error(`[lsodar] itol = ${itol} illegal`);
                return [t,this.terminate(istate)];
            }
            if (iopt < 0 || iopt > 1) {
                console.error(`[lsodar] iopt = ${iopt} illegal`);
                return [t,this.terminate(istate)];
            }
            if (jt == 3 || jt < 1 || jt > 5) {
                console.error(`[lsodar] jt = ${jt} illegal`);
                return [t,this.terminate(istate)];
            }
            this.jtyp = jt;
            if (jt > 2) {
                this.ml = iwork1;
                this.mu = iwork2;
                if (this.ml < 0 || this.ml >= this.n) {
                    console.error(`[lsodar] ml = ${this.ml} not between 1 and neq.`);
                    return [t,this.terminate(istate)];
                }
                if (this.mu < 0 || this.mu >= this.n) {
                    console.error(`[lsodar] mu = ${this.mu} not between 1 and neq.`);
                    return [t,this.terminate(istate)];
                }
            }
            // LSODAR additions
            if (ng < 0) {
                console.error(`[lsodar] ng = ${ng} is less than 0`)
                return [t, this.terminate(istate)]
            }
            if (istate == 3 && this.irfnd == 0 && ng != this.ngc) {
                console.error(`[lsodar] ng changed (from ${this.ngc} to ${ng}) illegally.`)
                console.error(`        i.e not immediately after a root was found`)
                return [t, this.terminate(istate)]
            }
            this.ngc = ng
            // END LSODAR additions

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
                    console.error(`[lsodar] ixpr = ${this.ixpr} is illegal`);
                    return [t,this.terminate(istate)];
                }
                this.mxstep = iwork6;
                if (this.mxstep < 0) {
                    console.error("[lsodar] mxstep < 0");
                    return [t,this.terminate(istate)];
                }
                if (this.mxstep == 0) this.mxstep = mxstp0;
                this.mxhnil = iwork7;
                if (this.mxhnil < 0) {
                    console.error("[lsodar] mxhnil < 0");
                    return [t,this.terminate(istate)];
                }
                if (this.mxhnil == 0) this.mxhnil = mxhnl0 // not in LSODA c
                if (istate == 1) {
                    h0 = rwork5;
                    this.mxordn = iwork8;
                    if (this.mxordn < 0) {
                        console.error(`[lsodar] mxordn = ${this.mxordn} is less than 0`);
                        return [t, this.terminate(istate)];
                    }
                    if (this.mxordn == 0) this.mxordn = 100;
                    this.mxordn = Math.min(this.mxordn, this.mord[1]);
                    this.mxords = iwork9;
                    if (this.mxords < 0) {
                        console.error(`[lsodar] mxords = ${this.mxords} is less than 0`);
                        return [t, this.terminate(istate)];
                    }
                    if (this.mxords == 0) this.mxords = 100;
                    this.mxords = Math.min(this.mxords, this.mord[2]);
                    if ((tout - t) * h0 < 0.) {
                        console.error(`[lsodar] tout = ${tout} behind t = ${t}. integration direction is given by ${h0}`);
                        return [t, this.terminate(istate)];
                    }
                }	/* end if ( *istate == 1 )  */
                hmax = rwork6;
                if (hmax < 0.) {
                    console.error("[lsodar] hmax < 0.");
                    return [t,this.terminate(istate)];
                }
                this.hmxi = 0.;
                if (hmax > 0)
                    this.hmxi = 1. / hmax;
                this.hmin = rwork7;
                if (this.hmin < 0.) {
                    console.error("[lsodar] hmin < 0.");
                    return [t,this.terminate(istate)];
                }
            }		/* end else   *//* end iopt = 1   */
        }			/* end if ( *istate == 1 || *istate == 3 )   */
        /*
        -----------------------------------------------------------------------
         set work array pointers and check lengths lrw and liw.
         if istate = 1, meth is initialized to 1 here to facilitate the
         checking of work space lengths.
         pointers to segments of rwork and iwork are named by prefixing l to
         the name of the segment.  e.g., the segment yh starts at rwork(lyh).
         segments of rwork (in order) are denoted  g0, g1, gx, yh, wm,
         ewt, savf, acor.
         if the lengths provided are insufficient for the current method,
         an error return occurs.  this is treated as illegal input on the
         first call, but as a problem interruption with istate = -7 on a
         continuation call.  if the lengths are sufficient for the current
         method but not for both methods, a warning message is sent.
        -----------------------------------------------------------------------
        */
        if (istate == 1) {  // block# 60
            this.meth = 1;
            this.nyh = this.n;
            lenyh = 1 + Math.max(this.mxordn, this.mxords)
            // also allocate memory for g0, g1, gx, yh, wm, ewt, savf, acor, ipvt.
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
            
            // LSODAR additions
            this.g0 = new Array(1 + ng)
            this.g1 = new Array(1 + ng)
            this.gx = new Array(1 + ng)
        }
        // check rtol and atol for legality. ------------------------------------
        if (istate == 1 || istate == 3) {   // block# 70
            rtoli = rtol[1];
            atoli = atol[1];
            for (i = 1; i <= this.n; i++) {
                if (itol >= 3)
                    rtoli = rtol[i];
                if (itol == 2 || itol == 4)
                    atoli = atol[i];
                if (rtoli < 0.) {
                    console.error(`[lsodar] rtol = ${rtoli} is less than 0.`);
                    return [t,this.terminate(istate)];
                }
                if (atoli < 0.) {
                    console.error(`[lsodar] atol = ${atoli} is less than 0.`);
                    return [t,this.terminate(istate)];
                }
            }		/* end for   */
        }			/* end if ( *istate == 1 || *istate == 3 )   */
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
                    console.error("[lsodar] itask = 4 or 5 and tcrit behind tout\n");
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
            f(t, y, this.yh[2], _data)
            this.nfe = 1;
            // load the initial value vector in yh.
            this.yp1 = this.yh[1];
            for (i = 1; i <= this.n; i++) {
                this.yp1[i] = y[i];
            }
            // load and invert the ewt array.  ( h is temporarily set to 1. )
            this.nq = 1;
            this.h = 1.;
            ewset(this.n,itol, rtol, atol, y, this.ewt);
            for (i = 1; i <= this.n; i++) {
                if (this.ewt[i] <= 0.) {
                    console.error(`[lsodar] ewt[${i}] = ${this.ewt[i]} <= 0.`);
                    t = this.terminate2(y, t);
                    return [t, istate];
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
                    console.error(`[lsodar] tout: ${tout} too close to t: ${t} to start integration`);
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
            if (rh > 1.) {
                h0 /= rh;
            }
            // load h with h0 and scale yh(*,2) by h0. ------------------------------
            this.h = h0;
            this.yp1 = this.yh[2];
            for (i = 1; i <= this.n; i++) {
                this.yp1[i] *= h0;
            }
            // check for a zero of g at t. ------------------------------------------
            this.irfnd = 0
            this.toutc = tout
            if (this.ngc != 0) {
                irt = this.rchek(1, g, neq, y, this.yh, this.nyh, this.g0, this.g1, this.gx, jroot, irt, _data)
                if (irt != 0) {
                    console.error(`[lsodar] one or more components of g has a root too near to the initial point: ${t}`);
                    return [t,this.terminate(istate)];
                }
            }
        }			/* if ( *istate == 1 )   */
        /*
        -----------------------------------------------------------------------
         block d.
         the next code block is for continuation calls only (istate = 2 or 3)
         and is to check stop conditions before taking a step.
         first, rchek is called to check for a root within the last step
         taken, other than the last root found there, if any.
         if itask = 2 or 5, and y(tn) has not yet been returned to the user
         because of an intervening root, return through block g.
        -----------------------------------------------------------------------
        */
        if (istate == 2 || istate == 3) {
            this.nslast = this.nst;
            irfp = this.irfnd
            if (this.ngc != 0) {
                if (itask == 1 || itask == 4) this.toutc = tout
                irt = this.rchek (2, g, neq, y, this.yh, this.nyh, this.g0, this.g1, this.gx, jroot, irt, _data)
                if (irt == 1) {
                    this.irfnd = 1
                    istate = 3
                    t = this.t0;
                    return [t, istate]
                }
                this.irfnd = 0
                if (irfp == 1 && this.tlast != this.tn && itask == 2) {
                    return this.successreturn(y, t, itask, ihit, tcrit, istate)
                }
            }
            switch (itask) {
                case 1:
                    if ((this.tn - tout) * this.h >= 0.) {
                        iflag = this.intdy(tout, 0, y, iflag);
                        if (iflag != 0) {
                            console.error("[lsodar] trouble from intdy, itask = %d, tout = %g\n", itask, tout);
                            return [t, this.terminate(istate)];
                        }
                        t = tout;
                        istate = 2;
                        this.illin = 0;
                        this.tlast = t
                        return [t, istate];
                    }
                    break;
                case 2:
                    break;
                case 3:
                    tp = this.tn - this.hu * (1. + 100. * this.ETA);
                    if ((tp - tout) * this.h > 0.) {
                        console.error("[lsodar] itask = ${itask} and tout behind tcur - hu");
                        return [t, this.terminate(istate)];
                    }
                    if ((this.tn - tout) * this.h < 0.) break;
                    return this.successreturn(y, t, itask, ihit, tcrit, istate);
                case 4:
                    tcrit = rwork1;
                    if ((this.tn - tcrit) * this.h > 0.) {
                        console.error("[lsodar] itask = 4 or 5 and tcrit behind tcur\n");
                        return [t, this.terminate(istate)];
                    }
                    if ((tcrit - tout) * this.h < 0.) {
                        console.error("[lsodar] itask = 4 or 5 and tcrit behind tout\n");
                        return [t, this.terminate(istate)];
                    }
                    if ((this.tn - tout) * this.h >= 0.) {
                        iflag = this.intdy(tout, 0, y, iflag);
                        if (iflag != 0) {
                            console.error(`[lsodar] trouble from intdy, itask = ${itask}, tout = ${tout}`);
                            return [t, this.terminate(istate)];
                        }
                        t = tout;
                        istate = 2;
                        this.illin = 0;
                        this.tlast = t
                        return [t, istate];
                    }
                case 5:
                    if (itask == 5) {
                        tcrit = rwork1;
                        if ((this.tn - tcrit) * this.h > 0.) {
                            console.error("[lsodar] itask = 4 or 5 and tcrit behind tcur\n");
                            return [t, this.terminate(istate)];
                        }
                    }
                    hmx = Math.abs(this.tn) + Math.abs(this.h);
                    ihit = Math.abs(this.tn - tcrit) <= (100. * this.ETA * hmx);
                    if (ihit) {
                        t = tcrit;
                        return this.successreturn(y, t, itask, ihit, tcrit, istate);
                    }
                    if (irfp == 1 && this.tlast != this.tn && itask == 5) {
                        return this.successreturn(y, t, itask, ihit, tcrit, istate);
                    }
                    tnext = this.tn + this.h * (1. + 4. * this.ETA);
                    if ((tnext - tcrit) * this.h <= 0.) {
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
                    console.error("[lsodar] %d steps taken before reaching tout\n", this.mxstep);
                    istate = -1;
                    t = this.terminate2(y, t);
                    return [t, istate];
                }
                ewset(this.n, itol, rtol, atol, this.yh[1], this.ewt);
                for (i = 1; i <= this.n; i++) {
                    if (this.ewt[i] <= 0.) {
                        console.error("[lsodar] ewt[%d] = %g <= 0.\n", i, this.ewt[i]);
                        istate = -6;
                        t = this.terminate2(y, t);
                        return [t, istate];
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
                    console.error(`         suggested scaling factor = ${tolsf}\n`);
                    return [t,this.terminate(istate)];
                }
                console.error(`lsoda -- at t = ${t}, too much accuracy requested\n`);
                console.error("         for precision of machine, suggested\n");
                console.error(`         scaling factor = ${tolsf}\n`);
                istate = -2;
                t = this.terminate2(y, t);
                return [t, istate];
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
                 then call rchek to check for a root within the last step.
                 then, if no root was found, check for stop conditions.
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
                            console.error("[lsodar] a switch to the stiff method has occurred ");
                        if (this.meth == 1)
                            console.error("[lsodar] a switch to the nonstiff method has occurred");
                        console.error("at t = %g, tentative step size h = %g, step nst = %d\n", this.tn, this.h, this.nst);
                    }
                }

                /*
                    root check
                */
                if (this.ngc != 0){
                    irt = this.rchek (3, g, neq, y, this.yh, this.nyh, this.g0, this.g1, this.gx, jroot, irt, _data)
                    if (irt == 1) {
                        this.irfnd = 1
                        istate = 3
                        t = this.t0;
                        // go to 425
                        this.tlast = t
                        return [t, istate]
                    }
                }
                // itask = 1.  if tout has been reached, interpolate. -------------------
                if (itask == 1) {
                    if ((this.tn - tout) * this.h < 0.)
                        continue;
                    iflag = this.intdy(tout, 0, y, iflag);
                    t = tout;
                    istate = 2;
                    this.illin = 0;
                    this.tlast = t
                    return [t, istate];
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
                // itask = 4.  see if tout or tcrit was reached.  adjust h if necessary.
                if (itask == 4) {
                    if ((this.tn - tout) * this.h >= 0.) {
                        iflag = this.intdy(tout, 0, y, iflag);
                        t = tout;
                        istate = 2;
                        this.illin = 0;
                        this.tlast = t
                        return [t, istate];
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
               kflag = -1, error test failed repeatedly or with abs(h) = hmin.
               kflag = -2, convergence failed repeatedly or with abs(h) = hmin.
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
                return [t, istate];
            }		/* end if ( kflag == -1 || kflag == -2 )   */
        }			/* end while   */
        return [t, istate];
    }				/* end lsoda   */

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
    successreturn(y: number[], t: number, itask: number, ihit: boolean, tcrit: number, istate: number) : [number,number] {
        this.yp1 = this.yh[1];
        for (let i = 1; i <= this.n; i++)
            y[i] = this.yp1[i];
        t = this.tn;
        if (itask == 4 || itask == 5)
            if (ihit)
                t = tcrit;
        istate = 2;
        this.illin = 0;
        this.tlast = t
        return [t, istate]
    }

    /*
    -----------------------------------------------------------------------
     This routine checks for the presence of a root in the
     vicinity of the current t, in a manner depending on the
     input flag job.  It calls subroutine roots to locate the root
     as precisely as possible.

     In addition to variables described previously, rchek
     uses the following for communication..
     
     job    = integer flag indicating type of call..
             job = 1 means the problem is being initialized, and rchek
                     is to look for a root at or very near the initial t.
             job = 2 means a continuation call to the solver was just
                     made, and rchek is to check for a root in the
                     relevant part of the step last taken.
             job = 3 means a successful step was just taken, and rchek
                     is to look for a root in the interval of the step.
     g0     = array of length ng, containing the value of g at t = t0.
             g0 is input for job >= 2 and on output in all cases.
     g1,gx  = arrays of length ng for work space.
     irt    = completion flag..
             irt = 0  means no root was found.
             irt = -1 means job = 1 and a root was found too near to t.
             irt = 1  means a legitimate root was found (job = 2 or 3).
                     on return, t0 is the root location, and y is the
                     corresponding solution vector.
     t0     = value of t at one endpoint of interval of interest.  only
             roots beyond t0 in the direction of integration are sought.
             t0 is input if job >= 2, and output in all cases.
             t0 is updated by rchek, whether a root is found or not.
     tlast  = last value of t returned by the solver (input only).
     toutc  = copy of tout (input only).
     irfnd  = input flag showing whether the last step taken had a root.
             irfnd = 1 if it did, = 0 if not.
     itaskc = copy of itask (input only).
     ngc    = copy of ng (input only).
    -----------------------------------------------------------------------
    */
    rchek (job: number, g: lsodar_root_func, neq: number, y: number[], yh: number[][], nyh: number, g0: number[], g1: number[], gx: number[], jroot: number[], irt: number, data: any) {
        // local variables
        var t1: number
        var x: number
        var hming: number
        var temp1: number
        var temp2: number
        var jflag: number
        var iflag: number
        var zroot: boolean


        irt = 0
        for (let i = 1; i <= this.ngc; i++) {
            jroot[i] = 0
        }
        hming = (Math.abs(this.tn) + Math.abs(this.h))*this.ETA*100.0e0

        // go to (100, 200, 300), job
        switch (job) {
            case 1:
                // evaluate g at initial t, and check for zero values. ------------------
                this.t0 = this.tn
                g (this.t0, y, g0, data)
                this.nge = 1
                zroot = false
                for (let i = 1; i <= this.ngc; i++) {
                    if (Math.abs(g0[i]) <= 0.0e0) zroot = true
                }
                if (zroot) {
                    // g has a zero at t.  look at g at t + (small increment). --------------
                    temp1 = hming * ((this.h >= 0) ? 1 : -1)
                    this.t0 = this.t0 + temp1
                    temp2 = temp1/this.h
                    for (let i = 1; i <= this.n; i++) {
                        // y[i] = y[i] + temp2*this.yh[i][2]
                        y[i] = y[i] + temp2*this.yh[2][i]
                    }
                    g (this.t0, y, g0, data)
                    this.nge = this.nge + 1
                    zroot = false
                    for (let i = 1; i <= this.ngc; i++) {
                        if (Math.abs(g0[i]) <= 0.0e0) zroot = true
                        
                    }
                    if (zroot) {
                        // g has a zero at t and also close to t.  take error return. -----------
                        irt = -1
                    }
                }
                return irt
        
            case 2:
                if(this.irfnd != 0) {
                    // if a root was found on the previous step, evaluate g0 = g(t0). -------
                    iflag = this.intdy(this.t0, 0, y, iflag)
                    g(this.t0, y, g0, data)
                    this.nge = this.nge + 1
                    zroot = false
                    for (let i = 1; i <= this.ngc; i++) {
                        if (Math.abs(g0[i]) <= 0) zroot = true
                    }
                    if (zroot){
                        //g has a zero at t0.  look at g at t + (small increment). -------------
                        temp1 = hming * ((this.h >= 0) ? 1 : -1)
                        this.t0 = this.t0 + temp1
                        if ((this.t0 - this.tn)*this.h < 0.0e0) {
                            iflag = this.intdy (this.t0, 0, y, iflag)
                        } else {
                            temp2 = temp1/this.h
                            for (let i = 1; i <= this.n; i++){
                                // y[i] = y[i] + temp2*this.yh[i][2]
                                y[i] = y[i] + temp2*this.yh[2][i]
                            }
                        }
                        g (this.t0, y, g0, data)
                        this.nge = this.nge + 1
                        zroot = false
                        for (let i = 1; i <= this.ngc; i++) {
                            if (Math.abs(g0[i]) > 0.) break
                            jroot[i] = 1
                            zroot = true
                        }
                        if (zroot) {
                            //g has a zero at t0 and also close to t0.  return root. ---------------
                            irt = 1
                            return irt

                        }
                    }
                }
                // here, g0 does not have a root
                // g0 has no zero components.  proceed to check relevant interval. ------
                if (this.tn == this.tlast) return irt
                // no break, we fall through the switch statement

            case 3:
                // set t1 to tn or toutc, whichever comes first, and get g at t1. -------
                if ( (this.itaskc != 2 && this.itaskc != 3 && this.itaskc != 5) &&
                    ((this.toutc - this.tn)*this.h < 0.0e0)
                ) {
                    t1 = this.toutc
                    if ((t1 - this.t0)*this.h <= 0.0e0) return irt
                    iflag = this.intdy (t1, 0, y, iflag)
                } else {
                    t1 = this.tn
                    for (let i = 1; i <= this.n; i++) {
                        // y[i] = yh[i][1]
                        y[i] = yh[1][i]
                    }
                }
                g (t1, y, g1, data)
                this.nge = this.nge + 1
                // call roots to search for root in interval from t0 to t1. -------------
                jflag = 0
                while(true){
                    [jflag, this.t0, t1, x] = this.roots(this.ngc, hming, jflag, this.t0, t1, g0, g1, gx, x, jroot)
                    if (jflag > 1) break
                    iflag = this.intdy (x, 0, y, iflag)
                    g (x, y, gx, data)
                    this.nge = this.nge + 1
                }
                this.t0 = x
                scopy(gx,g0)
                if (jflag == 4) return irt
                // found a root.  interpolate to x and return. --------------------------
                iflag = this.intdy (x, 0, y, iflag)
                irt = 1
                return irt

            default:
                break;
        }
        return irt
    } /*----------------------- end of subroutine rchek ----------------------- */

    /*
    -----------------------------------------------------------------------
     This subroutine finds the leftmost root of a set of arbitrary
     functions gi(x) (i = 1,...,ng) in an interval (x0,x1).  Only roots
     of odd multiplicity (i.e. changes of sign of the gi) are found.
     Here the sign of x1 - x0 is arbitrary, but is constant for a given
     problem, and -leftmost- means nearest to x0.
     The values of the vector-valued function g(x) = (gi, i=1...ng)
     are communicated through the call sequence of roots.
     The method used is the illinois algorithm.
     reference..

     Kathie L. Hiebert and Lawrence F. Shampine, implicitly defined
     output points for solutions of ode-s, sandia report sand80-0180,
     February, 1980.

     Description of parameters.
     ng     = number of functions gi, or the number of components of
              the vector valued function g(x).  input only.
     hmin   = resolution parameter in x.  input only.  when a root is
              found, it is located only to within an error of hmin in x.
              typically, hmin should be set to something on the order of
                   100 * uround * max(abs(x0),abs(x1)),
              where uround is the unit roundoff of the machine.
     jflag  = integer flag for input and output communication.
              on input, set jflag = 0 on the first call for the problem,
              and leave it unchanged until the problem is completed.
              (the problem is completed when jflag >= 2 on return.)
              on output, jflag has the following values and meanings..
              jflag = 1 means roots needs a value of g(x).  set gx = g(x)
                        and call roots again.
              jflag = 2 means a root has been found.  the root is
                        at x, and gx contains g(x).  (actually, x is the
                        rightmost approximation to the root on an interval
                        (x0,x1) of size hmin or less.)
              jflag = 3 means x = x1 is a root, with one or more of the gi
                        being zero at x1 and no sign changes in (x0,x1).
                        gx contains g(x) on output.
              jflag = 4 means no roots (of odd multiplicity) were
                        found in (x0,x1) (no sign changes).
     x0,x1  = endpoints of the interval where roots are sought.
              x1 and x0 are input when jflag = 0 (first call), and
              must be left unchanged between calls until the problem is
              completed.  x0 and x1 must be distinct, but x1 - x0 may be
              of either sign.  however, the notion of -left- and -right-
              will be used to mean nearer to x0 or x1, respectively.
              when jflag >= 2 on return, x0 and x1 are output, and
              are the endpoints of the relevant interval.
     g0,g1  = arrays of length ng containing the vectors g(x0) and g(x1),
              respectively.  when jflag = 0, g0 and g1 are input and
              none of the g0(i) should be be zero.
              when jflag >= 2 on return, g0 and g1 are output.
     gx     = array of length ng containing g(x).  gx is input
              when jflag = 1, and output when jflag >= 2.
     x      = independent variable value.  output only.
              when jflag = 1 on output, x is the point at which g(x)
              is to be evaluated and loaded into gx.
              when jflag = 2 or 3, x is the root.
              when jflag = 4, x is the right endpoint of the interval, x1.
     jroot  = integer array of length ng.  output only.
              when jflag = 2 or 3, jroot indicates which components
              of g(x) have a root at x.  jroot(i) is 1 if the i-th
              component has a root, and jroot(i) = 0 otherwise.

     Note.. this routine uses the common block /lsr001/ to save
     the values of certain variables between calls (own variables).
    -----------------------------------------------------------------------
    */
    roots (ng: number, hmin, jflag, x0: number, x1: number, g0: number[], g1: number[], gx: number[], x: number, jroot: number[]) : [number, number, number, number] {
        // local variables
        var zroot: boolean
        var xroot: boolean
        var sgnchg: boolean
        var tmax: number
        var t2: number
        var imxold: number
        var nxlast: number

        if (jflag == 1) {
            // check to see in which interval g changes sign. -----------------------
            imxold = this.imax
            this.imax = 0
            tmax = 0
            zroot = false
            for (let i = 1; i <= ng; i++) {
                if (Math.abs(gx[i]) > 0) {
                    // neither g0[i] nor gx[i] can be zero at this point. -------------------
                    if (Math.sign(g0[i]) == Math.sign(gx[i])) continue
                    t2 = Math.abs(gx[i]/(gx[i] - g0[i]))
                    if (t2 <= tmax) continue
                    tmax = t2
                    this.imax = i
                } else {
                    zroot = true
                }
            }

            if (this.imax > 0) {
                sgnchg = true
            } else {
                sgnchg = false
                this.imax = imxold
            }
            
            nxlast = this.last
            if (!sgnchg) {
                if (!zroot) {
                    // no sign change between x0 and x2.  replace x0 with x2. ---------------
                    scopy(gx,g0)
                    x0 = this.x2
                    this.last = 0
                    xroot = false
                } else {
                    // zero value at x2 and no sign change in (x0,x2), so x2 is a root. -----
                    x1 = this.x2
                    scopy(gx,g1)
                    xroot = true
                }
            } else {
                // sign change between x0 and x2, so replace x1 with x2. ----------------
                x1 = this.x2
                scopy(gx,g1)
                this.last = 1
                xroot = false
            }
            if (Math.abs(x1-x0) <= hmin) xroot = true
        } else {
        // jflag != 1.  check for change in sign of g or zero at x1. ----------
            this.imax = 0
            tmax = 0 
            zroot = false
            for (let i = 1; i <= ng; i++) {
                if (Math.abs(g1[i]) > 0) { 
                    // at this point, g0(i) has been checked and cannot be zero. ------------
                    if (Math.sign(g0[i]) == Math.sign(g1[i])) continue
                    t2 = Math.abs(g1[i]/(g1[i]-g0[i]))
                    if (t2 <= tmax) continue
                    tmax = t2
                    this.imax = i
                } else {
                    zroot = true
                }
            }
            
            if (this.imax > 0) {
                sgnchg = true
            } else {
                sgnchg = false
            }

            if(!sgnchg) {
                // no sign change in the interval.  check for zero at right endpoint. ---
                if (zroot) {
                    // zero value at x1 and no sign change in (x0,x1).  return jflag = 3. ---
                    x = x1
                    scopy(g1,gx)
                    for (let i = 1; i <= ng; i++) {
                        jroot[i] = 0
                        if (Math.abs(g1[i]) <= 0) jroot[i] = 1
                    }
                    jflag = 3
                    return [jflag, x0, x1, x]
                }

                // no sign changes in this interval.  set x = x1, return jflag = 4. -----
                scopy(g1,gx)
                x = x1
                jflag = 4
                return [jflag, x0, x1, x]
            } // end of sgnchg

            // c there is a sign change.  find the first root in the interval. --------
            xroot = false
            nxlast = 0
            this.last = 1
        } // end of jflag != 1

        while(true) {
            // repeat until the first root in the interval is found.  loop point. ---
            if (xroot) {
                // return with x1 as the root.  set jroot.  set x = x1 and gx = g1. -----
                jflag = 2
                x = x1
                scopy(g1,gx)
                for (let i = 1; i <= ng; i++) {
                    jroot[i] = 0
                    if (Math.abs(g1[i]) > 0) {
                        if (Math.sign(g0[i]) != Math.sign(g1[i])) jroot[i] = 1
                    } else {
                        jroot[i] = 1
                    }
                }
                return [jflag, x0, x1, x]
            }
            if (nxlast == this.last) { 
                if (this.last == 0) {
                    this.alpha = 2.0e0*this.alpha
                } else {
                    this.alpha = 0.5e0*this.alpha
                }
            } else {
                this.alpha = 1.0e0
            }
            this.x2 = x1 - (x1-x0)*g1[this.imax]/(g1[this.imax] - this.alpha*g0[this.imax])
            if ((Math.abs(this.x2-x0) < hmin) && (Math.abs(x1-x0) > 10.0e0*hmin)) {
                this.x2 = x0 + 0.1e0*(x1-x0)
            }
            jflag = 1
            x = this.x2
            // return to the calling routine to get a value of gx = g(x). -----------
            return [jflag, x0, x1, x]

        } // while loop
    // c----------------------- end of subroutine roots -----------------------
    }
}