# odepack.js

ODEPACK is a collection of Fortran solvers for the initial value problem for ordinary differential equation systems. This is a javascript translation of the pack that can be used in your browser and also in NodeJS.

The original ODEPACK is a collection of 9 solvers, currently there are only 2 solvers translated in the pack: LSODA and LSODAR. You can view and download the original ODEPACK at https://computing.llnl.gov/casc/odepack/.

Since the original solvers are written in Fortran which has index of 1 for arrays, I have decided to maintain the same pattern here. This means that all arrays are expected to start at index of 1 (not 0 as is traditional in javascript).

## Quick Start

Using LSODA

```js
var solver = new odepack.LSODA();

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

[t,istate] = solver.lsoda(f, neq, y, t, tout, itol, rtol, atol, itask, istate, iopt, null, jt, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, null);

console.log(`at t = ${t} y(1) = ${y[1]} y(2) = ${y[2]} y(3) = ${y[3]}`);

function f(t, y, ydot, data)
{
	ydot[1] = 1e4 * y[2] * y[3] - .04e0 * y[1];
	ydot[3] = 3e7 * y[2] * y[2];
	ydot[2] = -1 * (ydot[1] + ydot[3]);
}
```

```
at t = 0.4 y(1) = 0.9851712049052344 y(2) = 0.00003386380058450371 y(3) = 0.014794931294181262
```

Using LSODAR

```js
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

    console.log(` at t = ${t.toExponential(4)} y = ${y[1]} ${y[2]} ${y[3]}`);

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

function f2root(t, y, gout){
	gout[1] = y[1] - 1.e-4;
	gout[2] = y[3] - 1.e-2;
}
```

```
 at t = 2.6400e-1 y = 0.9899652943680691 0.000034705631931064796 0.010000000000000004
     the above line is a root,  jroot =  ,0,1
 at t = 4.0000e-1 y = 0.9851712049052344 0.00003386380058450371 0.014794931294181262
 at t = 4.0000e+0 y = 0.9055332934421416 0.000022406552427739382 0.0944443000054307
 at t = 4.0000e+1 y = 0.7158403474138746 0.000009186333848221152 0.28415046625227736
 at t = 4.0000e+2 y = 0.4505250299535515 0.0000032229637737336674 0.5494717470826751
 at t = 4.0000e+3 y = 0.183197567007961 8.941773348573555e-7 0.8168015388147041
 at t = 4.0000e+4 y = 0.038987292756483966 1.6219395881389404e-7 0.9610125450495185
 at t = 4.0000e+5 y = 0.0049363622000517075 1.98422087011259e-8 0.9950636179577332
 at t = 4.0000e+6 y = 0.0005161832764573236 2.0657864796937233e-9 0.9994838146577499
 at t = 2.0745e+7 y = 0.0001 4.00039494400783e-10 0.9998999995999543
     the above line is a root,  jroot =  ,1,0
 at t = 4.0000e+7 y = 0.00005179804694124399 2.072027333913957e-10 0.9999482017458503
 at t = 4.0000e+8 y = 0.000005283663311688205 2.113476114412524e-11 0.9999947163155508
 at t = 4.0000e+9 y = 4.658688623047484e-7 1.863476268264901e-12 0.9999995341292733
 at t = 4.0000e+10 y = 1.430506778106123e-8 5.722032331499939e-14 0.9999999856948669
```

## How to use

In the browser

```js
<script src="odepack.min.js" type="text/javascript"></script>
```

In NodeJS

```js
var odepack = require("./odepack.min");
```

## Using LSODA

There are many inputs to the LSODA function.

```ts
lsoda(f: lsoda_func, neq: number, y: number[], t: number, tout: number, 
    itol: number, rtol: number[], atol: number[],
    itask: number, istate: number, iopt: number, jac: jac_func, jt: number,
    iwork1: number, iwork2: number, iwork5: number, iwork6: number,
    iwork7: number, iwork8: number, iwork9: number,
    rwork1: number, rwork5: number, rwork6: number, rwork7: number, 
    _data: any) : [number, number]
```

>* f      = name of subroutine for right-hand side vector f.
>          this name must be declared external in the calling program.
>* neq    = number of first order ode-s.
>* y      = array of initial values, of length neq.
>* t      = the initial value of the independent variable.
>* tout   = first point where output is desired (.ne. t).
>* itol   = 1 or 2 according as atol (below) is a scalar or array.
>* rtol   = relative tolerance parameter (array).
>* atol   = absolute tolerance parameter (scalar or array).
>* itask  = 1 for normal computation of output values of y at t = tout.
>* istate = integer flag (input and output).  set istate = 1.
>* iopt   = 0 to indicate no optional inputs used.
>* jac    = name of subroutine for jacobian matrix.
>          use a dummy name.  
>* jt     = jacobian type indicator.  set jt = 2. 
>* _data  = optional data passed to the f, jac, and g (see lsodar) subroutine 

The program also has several optional inputs to set and will be used if you set iopt = 1 as mentioned above.

>* iwork1 = ml = these are the lower
>* iwork2 = mu = and upper half-bandwidths of the banded jacobian, excluding
                the main diagonal. the band is defined by the matrix locations
                (i,j) with i-ml <= j <= i+mu. ml and mu must satisfy 0 <= ml,mu
                <= neq-1. these are required if jt is 4 or 5, and ignored
                otherwise. ml and mu may in fact be the band parameters for
                a matrix to which df/dy is only approximately equal
>* iwork5 = ixpr = flag to generate extra printing at method switches.
                ipxr = 0 means no extra printing (the default).
                ipxr = 1 means print data on each switch
                t, h, and nst will be printed on the same logical unit as used
                for error messages.
>* iwork6 = mxstep = maximum number of (internally defined) steps allowed during one
                call to the solver.the default value is 500
>* iwork7 = mxhnil =  maximum number of messages printed (per problem) warning that
                t + h = t on a step (h = step size). this must be positive to
                result in a non-default value. the default value is 10
>* iwork8 = mxordn = the maximum order to be allowed for the nonstiff (adams) 
                method. the default value is 12. if mxordn exceeds the default
                value, it will be reduced to the default value. mxordn is
                held constant during the problem.
>* iwork9 = mxords = the maximum order to be allowed for the stiff (bdf) method.
                the default value is 5. if mxords exceeds the default value, 
                it will be reduced to the default value. mxords is held
                constant during the problem.
>* rwork1 = tcrit = 
>* rwork5 = h0 = the step size to be attempted on the first step. 
                the default value is determined by the solver
>* rwork6 = hmax = the maximum absolute step size allowed.
                the default value is infinite
>* rwork7 = hmin = the minimum absolute step size allowed.
                the default value is 0. (this lower bound is not enforced)
                on the final step before reaching tcrit when itask = 4 or 5)

func f and jac will need to be of type:
```ts
type lsoda_func = (t: number, y: number[], ydot: number[], data: any) => any
type jac_func = (t: number, y: number[], ml: number, mu: number, wm: number[][], data: any) => any
```           

On call return the program will return: 
>* t = corresponding value of independent variable.  this is
          tout
>* istate = 2  if lsoda was successful, negative otherwise.
          -1 means excess work done on this call (perhaps wrong jt).
          -2 means excess accuracy requested (tolerances too small).
          -3 means illegal input detected (see printed message).
          -4 means repeated error test failures (check all inputs).
          -5 means repeated convergence failures (perhaps bad jacobian
             supplied or wrong choice of jt or tolerances).
          -6 means error weight became zero during problem. (solution
             component i vanished, and atol or atol(i) = 0.)
          -7 means work space insufficient to finish (see messages).

The program will also modify

>* y = array of computed values of y(t) vector.

## Using LSODAR                

```ts
lsodar(f: lsoda_func, neq: number, y: number[], t: number, tout: number,
    itol: number, rtol: number[], atol: number[],
    itask: number, istate: number, iopt: number, jac: jac_func, jt: number,
    iwork1: number, iwork2: number, iwork5: number, iwork6: number,
    iwork7: number, iwork8: number, iwork9: number,
    rwork1: number, rwork5: number, rwork6: number, rwork7: number,
    g: lsodar_root_func, ng: number, jroot: number[],
    _data: any): [number, number]
```

Using LSODAR will be nearly identical to LSODA except for a few added function parameters:

>* g = name of subroutine for constraint functions, whose
          roots are desired during the integration.
          this name must be declared external in calling program.
>* ng = number of constraint functions g(i).  if there are none,
          set ng = 0, and pass a dummy name for g.
>* jroot = array showing roots found if istate = 3 on return.
          jroot(i) = 1 if g(i) has a root at t, or 0 otherwise.

func g will need to be of type:
```ts
type lsodar_root_func = (t: number, y: number[], gout: number[], data: any) => any
```

On call return, the program will return:

>* t = corresponding value of independent variable.  this is
          tout if istate = 2, or the root location if istate = 3,
          or the farthest point reached if lsodar was unsuccessful.
>* istate = 2 or 3  if lsodar was successful, negative otherwise.
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

The program will also modify

>* y = array of computed values of y(t) vector.
>* jroot = array showing roots found if istate = 3 on return.
          jroot(i) = 1 if g(i) has a root at t, or 0 otherwise.