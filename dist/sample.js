var odepack = require("./odepack.min");

var sim = new odepack.LSODA();

var neq = 3;
var y = [,1,0,0];
var t = 0
var tout = .4
var itol = 2
var rtol = [0,1e-4,1e-8,1e-4];
var atol = [0,1e-6,1e-10,1e-6];
var itask = 1;
var istate = 1;
var iopt = 0;
var jt = 2;

for (iout = 1; iout <= 12; iout++) {
    [t,istate] = sim.lsoda(lsoda_test1_f, neq, y, t, tout, itol, rtol, atol, itask, istate, iopt, null, jt, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, null);
    console.log(` at t= ${t.toExponential(4)} y= ${y[1].toExponential(6)} ${y[2].toExponential(6)} ${y[3].toExponential(6)}`);
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

function lsoda_test1_f(t, y, ydot, data)
{
	let offset = 1
	ydot[0+offset] = 1.0E4 * y[1+offset] * y[2+offset] - .04E0 * y[0+offset];
	ydot[2+offset] = 3.0E7 * y[1+offset] * y[1+offset];
	ydot[1+offset] = -1.0 * (ydot[0+offset] + ydot[2+offset]);
}

console.log("\nExample 2\n")

var neq = 3;
var y = [0,1,0,0];
var t = 0;
var tout = .4;
var itol = 2;
var rtol = [1.e-4,1.e-4,1.e-4,1.e-4];
var atol = [1.0e-6,1.0e-6,1.0e-10,1.0e-6];
var itask = 1;
var istate = 1;
var iopt = 0;
var jt = 2;
var ng = 2;
var jroot = new Array(3);

var lsodar = new odepack.LSODAR();

for (var i = 1; i < 13; i++) {
    [t, istate] = lsodar.lsodar(f2, neq, y, t, tout, itol, rtol, 
        atol, itask, istate, iopt, null, jt, 0,0,0,0,0,0,0,0,0,0,0,
        f2root,ng,jroot,null);

    console.log(` at t = ${t.toExponential(4)} y = ${y[1].toExponential(6)} ${y[2].toExponential(6)} ${y[3].toExponential(6)}`);

    if(istate == 2){
        tout = tout * 10;
    } else {
        console.log(`   the above line is a root,  jroot =  ${jroot}`);
        istate = 2;
        i--;
    }
}
function f2(t, y, ydot, data)
{
	ydot[1] = -.04*y[1] + 1.e4*y[2]*y[3];
	ydot[3] = 3e7*y[2]*y[2];
	ydot[2] = -ydot[1] - ydot[3];
}

function f2root(t, y, gout, data){
	gout[1] = y[1] - 1.e-4;
	gout[2] = y[3] - 1.e-2;
}
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
 at t =  4.0000e+10   y =  1.404280e-08  5.617126e-14  1.000000e+00`)