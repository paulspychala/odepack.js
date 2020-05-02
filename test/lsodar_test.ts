import LSODA from "../src/lsoda";
import LSODAR from "../src/lsodar";


function lsodar_test3_f(t: number, y: number[], ydot: number[], data: number[])
{
	ydot[1] = -.04*y[1] + 1.e4*y[2]*y[3]
	ydot[3] = 3.e7*y[2]*y[2]
	ydot[2] = -ydot[1] - ydot[3]
}

function lsodar_test3_g(t: number, y: number[], gout: number[]){
	gout[1] = y[1] - 1.e-4
	gout[2] = y[3] - 1.e-2
}

function lsodar_test3(){
	let neq = 3
	let y: number[] = [0,1,0,0]
	let t = 0.
	let tout = .4
	let itol = 2
	let rtol: number[] = [1.e-4,1.e-4,1.e-4,1.e-4]
	let atol: number[] = [1.0e-6,1.0e-6,1.0e-10,1.0e-6]
	let itask = 1
	let istate = 1
	let iopt = 0
	let lrw = 76
	let liw = 23
	let jt = 2
	let ng = 2
	let jroot: number[] = new Array(3)

	let lsodar = new LSODAR()
	console.log(`\n Running LSODAR test\n`)

	for (let i = 1; i < 13; i++) {
		[t, istate] = lsodar.lsodar(lsodar_test3_f, neq, y, t, tout, itol, rtol, 
			atol, itask, istate, iopt, null, jt, 0,0,0,0,0,0,0,0,0,0,0,
			lsodar_test3_g,ng,jroot,[])

		console.log(` at t = ${t.toExponential(4)} y = ${y[1].toExponential(6)} ${y[2].toExponential(6)} ${y[3].toExponential(6)}`);

		if(istate == 2){
			tout = tout * 10
		} else if(istate <= 0){
			console.log(`error istate = ${istate}`);
		} else {
			console.log(`   the above line is a root,  jroot =  ${jroot}`)
			istate = 2
			i--
		}
	}

	console.log(`\n no. steps = ${lsodar.nst}  no. f-s = ${lsodar.nfe}  no. j-s =  ${lsodar.nje}  no. g-s = ${lsodar.nge}`)
	console.log(` method last used = ${lsodar.mused}   last switch was at t =  ${lsodar.tsw.toExponential(4)}`)

	console.log(`
 Expected Output:
 at t =  2.6400e-01   y =  9.899653e-01  3.470563e-05  1.000000e-02
      the above line is a root,  jroot =    0    1
 at t =  4.0000e-01   y =  9.851712e-01  3.386380e-05  1.479493e-02
 at t =  4.0000e+00   y =  9.055333e-01  2.240655e-05  9.444430e-02
 at t =  4.0000e+01   y =  7.158403e-01  9.186334e-06  2.841505e-01
 at t =  4.0000e+02   y =  4.505250e-01  3.222964e-06  5.494717e-01
 at t =  4.0000e+03   y =  1.831975e-01  8.941774e-07  8.168016e-01
 at t =  4.0000e+04   y =  3.898730e-02  1.621940e-07  9.610125e-01
 at t =  4.0000e+05   y =  4.936363e-03  1.984221e-08  9.950636e-01
 at t =  4.0000e+06   y =  5.161831e-04  2.065786e-09  9.994838e-01
 at t =  2.0745e+07   y =  1.000000e-04  4.000395e-10  9.999000e-01
      the above line is a root,  jroot =    1    0
 at t =  4.0000e+07   y =  5.179817e-05  2.072032e-10  9.999482e-01
 at t =  4.0000e+08   y =  5.283401e-06  2.113371e-11  9.999947e-01
 at t =  4.0000e+09   y =  4.659031e-07  1.863613e-12  9.999995e-01
 at t =  4.0000e+10   y =  1.404280e-08  5.617126e-14  1.000000e+00
 
 no. steps = 361  no. f-s = 693  no. j-s =  64  no. g-s = 390
 method last used = 2   last switch was at t =  6.0092e-03
	`)	
}



/*
-----------------------------------------------------------------------
 Demonstration program for the DLSODAR package.
 This is the version of 14 June 2001.

 This version is in double precision.

 The DLSODAR package is used to solve two simple problems,
 one nonstiff and one intermittently stiff.
 If the errors are too large, or other difficulty occurs,
 a warning message is printed.  All output is on unit lout = 6.
-----------------------------------------------------------------------
*/
export default function run_lsodar_tests(){
	let nerr = 0

	console.log(`\nDemonstration program for DLSODAR package\n`)
	nerr += lsodar_test1()
	nerr += lsodar_test2() // can't change jac function
	lsodar_test3()
	console.log(`Total number of errors encountered = ${nerr}`)
}
/*
-----------------------------------------------------------------------
First problem.
The initial value problem is:
dy/dt = ((2*log(y) + 8)/t - 5)*y,  y(1) = 1,  1 .le. t .le. 6
The solution is  y(t) = exp(-t**2 + 5*t - 4)
The two root functions are:
g1 = ((2*log(y)+8)/t - 5)*y (= dy/dt)  (with root at t = 2.5),
g2 = log(y) - 2.2491  (with roots at t = 2.47 and 2.53)
-----------------------------------------------------------------------
*/
function lsodar_test1(){
	// Set all input parameters and print heading.
	let neq = 1
	let y: number[] = [0,1.0]
	let t = 1.0
	let tout = 2.0
	let itol = 1
	let rtol: number[] = [1.0e-6,1.0e-6]
	let atol: number[] = [1.0e-6,1.0e-6]
	let itask = 1
	let istate = 1
	let iopt = 0
	let lrw = 44
	let liw = 21
	let jt = 2
	let ng = 2
	let jroot: number[] = new Array(ng)

	let nerr = 0

	console.log(`\nFirst problem\n`)
	console.log(`Problem is  dy/dt = ((2*log(y)+8)/t - 5)*y,  y(1) = 1`)
	console.log(`Solution is  y(t) = exp(-t**2 + 5*t - 4)`)
	console.log(`Root functions are:`)
	console.log(`          g1 = dy/dt  (root at t = 2.5)`)
	console.log(`          g2 = log(y) - 2.2491  (roots at t = 2.47 and t = 2.53)`)
	console.log(`itol ='${itol}'   rtol ='${rtol[1]}'   atol ='${atol[1]}`)
	console.log(`jt ='${jt}\n`)

	let lsodar = new LSODAR()

	// Call DLSODAR in loop over tout values 2,3,4,5,6.
	let ero = 0.0
	while (tout < 7) {
		[t,istate] = lsodar.lsodar(f1,neq,y,t,tout,itol,rtol,atol,itask,istate,iopt,null,jt,
			0,0,0,0,0,0,0,0,0,0,0,gr1,ng,jroot,[])
		// Print y and error in y, and print warning if error too large.
		let yt = Math.exp(-t*t + 5.0*t - 4.0)
		let er = y[1] - yt
		console.log(`At t =  ${t.toExponential(7)}   y =  ${y[1].toExponential(7)}  error =  ${er.toExponential(4)}`)
		if (istate < 0) break
		er = Math.abs(er)/(rtol[1]*Math.abs(y[1]) + atol[1])
		ero = Math.max(ero,er)
		if (er > 1000) {
			console.log(`Warning: error exceeds 1000 * tolerance`)
			nerr = nerr + 1
		}
		if (istate == 3) {
			// c If a root was found, write results and check root location.
			// c Then reset istate to 2 and return to DLSODAR call.
			console.log(`\nRoot found at t = ${t}  jroot = ${jroot}`)
			let errt: number
			if (jroot[1] == 1) errt = t - 2.5
			if (jroot[2] == 1 && t < 2.5) errt = t - 2.47
			if (jroot[2] == 1 && t > 2.5) errt = t - 2.53
			console.log(`Error in t location of root is ${errt.toExponential(4)}\n`)

			if (Math.abs(errt) > 1.0e-3) {
				console.log(`Warning: root error exceeds 1.0d-3`)
				nerr = nerr + 1
			}
			istate = 2
		} else {
			// If no root found, increment tout and loop back.
			tout = tout + 1
		}
	}

	// Problem complete.  Print final statistics.
	if (istate < 0) nerr = nerr + 1
	let nfea = lsodar.nfe
	if (jt == 2) nfea = lsodar.nfe - neq*lsodar.nje
	console.log(`Final statistics for this run:`)
	console.log(`rwork size =${42}   iwork size =${21}`)
	console.log(`number of steps =  ${lsodar.nst}`)
	console.log(`number of f-s   =  ${lsodar.nfe}`)
	console.log(`(excluding j-s) =  ${nfea}`)
	console.log(`number of j-s   =  ${lsodar.nje}`)
	console.log(`number of g-s   =  ${lsodar.nge}`)
	console.log(`error overrun = ${ero}`)
	return nerr
}

/*
-----------------------------------------------------------------------
Second problem (Van der Pol oscillator).
The initial value problem is (after reduction of 2nd order ODE):
dy1/dt = y2,  dy2/dt = 100*(1 - y1**2)*y2 - y1,
y1(0) = 2,  y2(0) = 0,  0 .le. t .le. 200
The root function is  g = y1.
An analytic solution is not known, but the zeros of y1 are known
to 15 figures for purposes of checking the accuracy.
-----------------------------------------------------------------------
*/
function lsodar_test2(){
// Set tolerance parameters and print heading.
let itol = 2
let nerr = 0
let rtol: number[] = [1.0e-6,1.0e-6,1.0e-6]
let atol: number[] = [0,1.0e-6,1.0e-4]
console.log(`Second problem (Van der Pol oscillator)`)
console.log(`Problem is dy1/dt = y2,  dy2/dt = 100*(1-y1**2)*y2 - y1`)
console.log(`           y1(0) = 2,  y2(0) = 0`)
console.log(`Root function is  g = y1`)
console.log(`itol =${itol}   rtol =${rtol}  atol =${atol}`)

// Loop over jt = 1, 2.  Set remaining parameters and print jt.
for (let jt = 1; jt < 3; jt++) {
	let neq = 2
	let y: number[] = [,2.0,0]
	let t = 0.0
	let tout = 20.0
	let itask = 1
	let istate = 1
	let iopt = 0
	let lrw = 57
	let liw = 22
	let ng = 1

	let jroot = new Array(ng)
	let lsodar = new LSODAR()
	console.log(`Solution with jt = ${jt}`)
	// Call DLSODAR in loop over tout values 20,40,...,200.
	for (let iout = 1; iout < 11; iout++) {
		[t,istate] = lsodar.lsodar(f2,neq,y,t,tout,itol,rtol,atol,itask,istate,
			iopt,jac2,jt,0,0,0,0,0,0,0,0,0,0,0,gr2,ng,jroot,[])
		// Print y1 and y2.
		console.log(`At t = ${t}  y1 =  ${y[1]}  y2 =    ${y[2]}`)
		if (istate < 0) break
		if (istate == 3) {
			// If a root was found, write results and check root location.
			// Then reset istate to 2 and return to DLSODAR call.
			console.log(`Root found at t = ${t}`)
			let kroot = Math.floor(t/81.2 + 0.5)
			let tzero = 81.17237787055 + (kroot-1)*81.41853556212
			let errt = t - tzero
			console.log(`Error in t location of root is ${errt}`)
			if (Math.abs(errt) > 1.0e-1) {
				console.log(`Warning: root error exceeds 1.0d-1`)
				nerr = nerr + 1
			}
			istate = 2       
		} else {
			// If no root found, increment tout and loop back.
			tout = tout + 20
		}
	} // end for iout 1 < 11

	// Problem complete.  Print final statistics.
	if (istate < 0) nerr = nerr + 1
	let nfea = lsodar.nfe
	if (jt == 2) nfea = lsodar.nfe - neq*lsodar.nje
	console.log(`Final statistics for this run:`)
	console.log(`rwork size =${42}   iwork size =${21}`)
	console.log(`number of steps =  ${lsodar.nst}`)
	console.log(`number of f-s   =  ${lsodar.nfe}`)
	console.log(`(excluding j-s) =  ${nfea}`)
	console.log(`number of j-s   =  ${lsodar.nje}`)
	console.log(`number of g-s   =  ${lsodar.nge}`)
} // end of for 1 < jt < 2 

return nerr
}

function f1 (t: number, y: number[], ydot: number[], data: number[]){
ydot[1] = ((2.0*Math.log(y[1]) + 8.0)/t - 5.0)*y[1]
}

function gr1 (t: number, y: number[], groot: number[]){
groot[1] = ((2.0*Math.log(y[1]) + 8.0)/t - 5.0)*y[1]
groot[2] = Math.log(y[1]) - 2.2491
}

function f2 (t: number, y: number[], ydot: number[], data: number[]){
ydot[1] = y[2]
ydot[2] = 100.0*(1.0 - y[1]*y[1])*y[2] - y[1]
}

function jac2 (t: number, y: number[], ml, mu, pd: number[][]){
pd[1][1] = 0.0
pd[1][2] = 1.0
pd[2][1] = -200.0*y[1]*y[2] - 1.0
pd[2][2] = 100.0*(1.0 - y[1]*y[1])

}

function gr2 (t: number, y: number[], groot: number[]){
groot[1] = y[1]
}