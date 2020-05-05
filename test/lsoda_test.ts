/*
-----------------------------------------------------------------------
 Demonstration program for the DLSODA package.
 This is the version of 14 June 2001.

 This version is in double precision.

 The package is used to solve two simple problems,
 one with a full Jacobian, the other with a banded Jacobian,
 with the 2 appropriate values of jt in each case.
 If the errors are too large, or other difficulty occurs,
 a warning message is printed.  All output is on unit lout = 6.
-----------------------------------------------------------------------
*/
import {LSODA} from "../src/lsoda";

export function run_lsoda_tests() {
     lsoda_test1();
     lsoda_test()
}

function lsoda_test1_f(t: number, y: number[], ydot: number[], data: number[])
{
	let offset = 1
	ydot[0+offset] = 1.0E4 * y[1+offset] * y[2+offset] - .04E0 * y[0+offset];
	ydot[2+offset] = 3.0E7 * y[1+offset] * y[1+offset];
	ydot[1+offset] = -1.0 * (ydot[0+offset] + ydot[2+offset]);
}

function lsoda_test1()
{
	let rwork1: number
	let rwork5: number
	let rwork6: number
	let rwork7: number;
	let atol: number[] = new Array(4)
	let rtol: number[] = new Array(4)
	let t: number
	let tout: number
	let y: number[] = new Array(4)
	let iwork1: number
	let iwork2: number
	let iwork5: number
	let iwork6: number
	let iwork7: number
	let iwork8: number
	let iwork9: number;
	let neq: number = 3;
	let itol: number
	let itask: number
	let istate: number
	let iopt: number
	let jt: number
     let iout: number;
     let expectedresults: number[][] = [ // t, y[1], y[2], y[3]
          [4.0000e-01,9.851712e-01,3.386380e-05,1.479493e-02],
          [4.0000e+00,9.055333e-01,2.240655e-05,9.444430e-02],
          [4.0000e+01,7.158403e-01,9.186334e-06,2.841505e-01],
          [4.0000e+02,4.505250e-01,3.222964e-06,5.494717e-01],
          [4.0000e+03,1.831976e-01,8.941773e-07,8.168015e-01],
          [4.0000e+04,3.898729e-02,1.621940e-07,9.610125e-01],
          [4.0000e+05,4.936362e-03,1.984221e-08,9.950636e-01],
          [4.0000e+06,5.161833e-04,2.065787e-09,9.994838e-01],
          [4.0000e+07,5.179804e-05,2.072027e-10,9.999482e-01],
          [4.0000e+08,5.283675e-06,2.113481e-11,9.999947e-01],
          [4.0000e+09,4.658667e-07,1.863468e-12,9.999995e-01],
          [4.0000e+10,1.431100e-08,5.724404e-14,1.000000e+00]]
     var myresults: number[]
     var diff: number
     var diffAtol = .01
     var diffRtol = .01

	iwork1 = iwork2 = iwork5 = iwork6 = iwork7 = iwork8 = iwork9 = 0;
	rwork1 = rwork5 = rwork6 = rwork7 = 0.0;
	y[1] = 1.0E0;
	y[2] = 0.0E0;
	y[3] = 0.0E0;
	t = 0.0E0;
	tout = 0.4E0;
	itol = 2;
	rtol[0] = 0.0;
	atol[0] = 0.0;
	rtol[1] = rtol[3] = 1.0E-4;
	rtol[2] = 1.0E-8;
	atol[1] = 1.0E-6;
	atol[2] = 1.0E-10;
	atol[3] = 1.0E-6;
	itask = 1;
	istate = 1;
	iopt = 0;
	jt = 2;

	let sim = new LSODA()

	for (iout = 1; iout <= 12; iout++) {
		[t,istate] = sim.lsoda(lsoda_test1_f, neq, y, t, tout, itol, rtol, atol, itask, istate, iopt, null, jt,
		      iwork1, iwork2, iwork5, iwork6, iwork7, iwork8, iwork9,
			  rwork1, rwork5, rwork6, rwork7, []);
          console.log(` at t= ${t.toExponential(4)} y= ${y[1].toExponential(6)} ${y[2].toExponential(6)} ${y[3].toExponential(6)}`);
          myresults = [t,y[1],y[2],y[3]]
          for (let i=0; i < myresults.length; i++) {
               diff = Math.abs(myresults[i] - expectedresults[iout-1][i])
               if (diff > diffAtol + Math.abs(diffRtol*expectedresults[iout-1][i])){
                    console.log(`Error: expected: ${expectedresults[iout-1]} got ${myresults[i]}`)
               }
          }
		if (istate <= 0) {
			console.log(`error istate = ${istate}`);
		}
		tout = tout * 10.0E0;
	}

	console.log(`
 The correct answer (up to certain precision):
 at t=   4.0000e-01 y=   9.851712e-01   3.386380e-05   1.479493e-02
 at t=   4.0000e+00 y=   9.055333e-01   2.240655e-05   9.444430e-02
 at t=   4.0000e+01 y=   7.158403e-01   9.186334e-06   2.841505e-01
 at t=   4.0000e+02 y=   4.505250e-01   3.222964e-06   5.494717e-01
 at t=   4.0000e+03 y=   1.831976e-01   8.941773e-07   8.168015e-01
 at t=   4.0000e+04 y=   3.898729e-02   1.621940e-07   9.610125e-01
 at t=   4.0000e+05 y=   4.936362e-03   1.984221e-08   9.950636e-01
 at t=   4.0000e+06 y=   5.161833e-04   2.065787e-09   9.994838e-01
 at t=   4.0000e+07 y=   5.179804e-05   2.072027e-10   9.999482e-01
 at t=   4.0000e+08 y=   5.283675e-06   2.113481e-11   9.999947e-01
 at t=   4.0000e+09 y=   4.658667e-07   1.863468e-12   9.999995e-01
 at t=   4.0000e+10 y=   1.431100e-08   5.724404e-14   1.000000e+00`)

    return 0;
}


function lsoda_test(){
     let testOtherMiter = false
     let tout1 = 16.921743e0
     let dtout = 17.341162e0
     let nerr = 0
     let itol = 1
     let rtol = [0,0,0]
     let atol = [1.0e-8,1.0e-8,1.0e-8]
     let lrw = 522
     let liw = 45
     let iopt = 0

     // First problem
     let expectedresults: number[][] = [[0.84609E+01,0.16731E+01,-0.464E-01,2,4,0.209E+00,0.311E+00],
          [0.16922E+02,-0.11574E-03,-0.141E+02,1,7,0.206E-02,0.158E+02],
          [0.25592E+02,-0.16828E+01,0.459E-01,2,4,0.240E+00,0.174E+02],
          [0.34263E+02,0.21448E-03,0.141E+02,1,8,0.293E-02,0.332E+02]];
     // the expected results are the same for both jt=1 and jt=2
     // only the number of f-s change
     let myresults: number[]
     let diff
     var diffAtol = .01
     var diffRtol = .01

     let neq = 2
     let nout = 4
     console.log(`
Demonstration program for DLSODA package

Problem 1:   Van der Pol oscillator:'/
             xdotdot - 20*(1 - x**2)*xdot + x = 0, x(0) = 2, xdot(0) = 0
neq =${neq} itol = ${itol}   rtol =${rtol}   atol = ${atol}
`)

    for (let jt = 1; jt <= 2; jt++) {
          console.log(`Solution with jt = ${jt}`)
          console.log(`  t               x               xdot       meth   nq     h           tsw`)
        
          let t = 0.0
          let y = [,2.0,0]
          let itask = 1
          let istate = 1
          let dtout0 = 0.5*tout1
          let dtout1 = 0.5*dtout
          let tout = dtout0
          let ero = 0.0
          
          var lsoda = new LSODA()
          for (let iout = 1; iout <= nout; iout++) {
               [t,istate] = lsoda.lsoda(f1,neq,y,t,tout,itol,rtol,atol,itask,istate,iopt,jac1,jt,0,0,0,0,0,0,0,0,0,0,0, [])
               console.log(`${t.toExponential(4)}    ${y[1].toExponential(4)}     ${y[2].toExponential(2)}    ${lsoda.mused}     ${lsoda.nqu}     ${lsoda.hu.toExponential(2)}    ${lsoda.tsw.toExponential(2)}`)
               if (istate < 0) break;
               // compare results here:
               myresults = [t,y[1],y[2],lsoda.mused,lsoda.nqu, lsoda.hu]
               for (let i = 0; i < expectedresults.length; i++) {
                    diff = Math.abs(myresults[i] - expectedresults[iout-1][i])
                    if (diff > diffAtol + Math.abs(diffRtol*expectedresults[iout-1][i])) {
                         console.log(` Error: t value problems Expected: ${expectedresults[iout-1][i]} and got ${t}`)
                         nerr++
                    }
               }
               if (iout == 1) tout = tout + dtout0
               if (iout > 1) tout = tout + dtout1
          }    
          if (istate < 0) nerr = nerr + 1
          let nfea = lsoda.nfe
          if (jt == 2) nfea = lsoda.nfe - neq*lsoda.nje
          console.log(`Final statistics for this run:
 rwork size =${lrw}   iwork size =${liw}
 number of steps = ${lsoda.nst}
 number of f-s   = ${lsoda.nfe}
 (excluding J-s) = ${nfea}
 number of J-s   = ${lsoda.nje}
 max. error at root = ${ero}`)
     }

     if(!testOtherMiter) return;

     // Second problem
     expectedresults = [[0.10000E-01,0.476E-06,1,2,0.714E-02,0.000E+00],
          [0.10000E+00,0.988E-06,1,4,0.343E-01,0.000E+00],
          [0.10000E+01,0.431E-06,1,5,0.724E-01,0.000E+00],
          [0.10000E+02,0.558E-07,1,3,0.323E+00,0.000E+00],
          [0.10000E+03,0.127E-11,2,1,0.239E+03,0.170E+02]]
     // again expected results are the same fore jt=4 and jt=5
     // only the run stats will change (number of f-s)

     neq = 25
     let ml = 5
     let mu = 0
     let mband = ml + mu + 1
     let y = new Array(neq+1)
     atol = [1.0e-6,1.0e-6]
     nout = 5
     console.log(`\n\nProblem 2: ydot = A * y , where A is a banded lower triangular matrix'
            derived from 2-D advection PDE'/
neq = ${neq}   ml = ${ml}   mu = ${mu}
itol = ${itol}   rtol = ${rtol}   atol = ${atol}`)
     for(let jt = 4; jt <= 5; jt++){
          console.log(` Solution with jt = ${jt}\n     t       max.err.    meth   nq      h            tsw`)
          let t = 0.0
          for (let i = 0; i <= neq; i++) {
               y[i] = 0
          }
          y[1] = 1
          let itask = 1
          let istate = 1
          let tout = 0.01
          let ero = 0.0
          let erm
          let lsoda = new LSODA()
          for (let iout = 1; iout <= nout; iout++) {
               [t,istate] = lsoda.lsoda(f2,neq,y,t,tout,itol,rtol,atol,itask,istate,iopt,jac2,jt,ml,mu,0,0,0,0,0,0,0,0,0, [])
               erm = edit2(y,t, erm)

               console.log(`${t.toExponential(4)}    ${erm.toExponential(2)}      ${lsoda.mused}    ${lsoda.nqu}    ${lsoda.hu.toExponential(2)}    ${lsoda.tsw.toExponential(2)}`)
               if (istate < 0) break;
               // compare results here:
               myresults = [t,erm,lsoda.mused,lsoda.nqu,lsoda.hu,lsoda.tsw]
               for (let i = 0; i < expectedresults.length; i++) {
                    diff = Math.abs(myresults[i] - expectedresults[iout-1][i])
                    if (diff > diffAtol + Math.abs(diffRtol*expectedresults[iout-1][i])) {
                         console.log(` Error: t value problems Expected: ${expectedresults[iout-1][i]} and got ${t}`)
                         nerr++
                    }
               }
               tout = tout*10.0
          }
          if (istate < 0) nerr = nerr + 1
          // nst = iwork(11)
          // nfe = iwork(12)
          // nje = iwork(13)
          // lenrw = iwork(17)
          // leniw = iwork(18)
          let nfea = lsoda.nfe
          if (jt == 5) nfea = lsoda.nfe - mband*lsoda.nje
          console.log(`Final statistics for this run:
 rwork size =${lrw}   iwork size =${liw}
 number of steps = ${lsoda.nst}
 number of f-s   = ${lsoda.nfe}
 (excluding J-s) = ${nfea}
 number of J-s   = ${lsoda.nje}
 max. error at root = ${ero}`)
     }
     console.log(`Number of errors encountered = ${nerr}`)
}

function f1 (t: number, y: number[], ydot: number[], _data) {
     //  integer neq
     //  double precision t, y, ydot
     //  dimension y(neq), ydot(neq)
      ydot[1] = y[2]
      ydot[2] = 20.0*(1.0 - y[1]*y[1])*y[2] - y[1]
      return
}

function jac1 (t, y, ml, mu, pd) {
     //  integer neq, ml, mu, nrowpd
     //  double precision t, y, pd
     //  dimension y(neq), pd(nrowpd,neq)
      pd[1][1] = 0.0
      pd[1][2] = 1.0
      pd[2][1] = -40.0*y[1]*y[2] - 1.0
      pd[2][2] = 20.0*(1.0 - y[1]*y[1])
      return
}

function f2 (t, y: number[], ydot: number[], _data) {
     //  integer neq, i, j, k, ng
     //  double precision t, y, ydot, alph1, alph2, d
     //  dimension y(neq), ydot(neq)
     //  data alph1/1.0/, alph2/1.0/, ng/5/
     let i: number
     let j: number
     let k: number
     let d: number
     let alph1 = 1
     let alph2 = 1
     let ng = 5
     for (j = 1; j <= ng; j++) {
          for (i = 1; i <= ng; i++) {
               k = i + (j - 1)*ng
               d = -2.0*y[k]
               if (i != 1) d = d + y[k-1]*alph1
               if (j != 1) d = d + y[k-ng]*alph2
               ydot[k] = d
          }
     }
}

function jac2 (t: number, y: number[], ml: number, mu: number, pd: number[][]) {
     let j: number
     //   integer neq, ml, mu, nrowpd, j, mband, mu1, mu2, ng
     //   double precision t, y, pd, alph1, alph2
     //   dimension y(neq), pd(nrowpd,neq)
     //   data alph1/1.0/, alph2/1.0/, ng/5/
     let mband = ml + mu + 1
     let mu1 = mu + 1
     let mu2 = mu + 2
     let alph1 = 1
     let alph2 = 1
     let ng = 5
     for (j = 1; j <= pd.length; j++) {
          pd[mu1][j] = -2.0
          pd[mu2][j] = alph1
          pd[mband][j] = alph2
     }
     for (j = ng; j < pd.length; j = j + ng) {
          // do 20 j = ng,neq,ng
          pd[mu2][j] = 0.0
     }
}

function edit2 (y: number[], t: number, erm: number) : number {
     //  integer i, j, k, ng
     //  double precision y, t, erm, alph1, alph2, a1, a2, er, ex, yt
     //  dimension y(*)
     //  data alph1/1.0/, alph2/1.0/, ng/5/
     let i: number
     let j: number
     let k: number
     let ng = 5
     let alph1 = 1
     let alph2 = 1
     erm = 0.0
     let a1: number
     let a2: number
     let er: number
     let ex: number
     let yt: number
     if (t == 0.0) return erm
     ex = 0.0
     if (t <= 30.0) ex = Math.exp(-2.0*t)
     a2 = 1.0
     for (j = 1; j <= ng; j++) {
          a1 = 1.0
          for (i = 1; i <= ng; i++) {
               k = i + (j - 1)*ng
               yt = t**(i+j-2)*ex*a1*a2
               er = Math.abs(y[k]-yt)
               erm = Math.max(erm,er)
               a1 = a1*alph1/i
          }
          a2 = a2*alph2/j
     }
     return erm
}