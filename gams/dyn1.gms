$TITLE This is a calibrated dynamic ramsey model


$ontext

Edward J. Balistreri
Senior Associate
Charles River Associates, Inc.
ejb@crai.com
http://www.tomco.net/~balistre

In this model the social accounts are taken at
face value and no assumption is made about the 
interest/discount rate.  

Original Base Year Data (rectangular SAM)


	   Y	  I	|	  RA
------------------------------------
PY	 200	-20	|	-180
RENT	 -50		|	  50
WAGE	-150		|	 150
Savings		 20	|	 -20

$offtext

Parameter
	Y0	Benchmark Gross Output,
	KS0	Benchmark Capital Supply,
	LS0	Benchmark Labor Supply,
	I0	Benchmark Savings,
	C0	Benchmark Consumption;

Y0	= 200;
KS0	=  50;
LS0	= 150;
I0	=  20;
C0	= 180;

*----------Setup the Dynamics----------*
Set	T		Time Periods	/2000*2100/
	TFIRST(T)	First Period	/2000/
	TLAST(T)	Terminal Period /2100/;

Scalars
	R	Rate of Interest
	G	Growth Rate		/0.02/,
	D	Depreciation Rate	/0.02/;

Parameter
	K0	Benchmark Capital Stock,
	RK0	Benchmark Return on a Unit of K
	PK0	Benchmark Price of a Unit of K
	QREF(T) Growth Path for Quantities
	PREF(T) Present Value Price Paths;

* 1) Solve for the value of the initial Capital Stock
* Investment must cover growth and depreciation.
* I0 = K0*(1+G) - K0*(1+D)	solve for K0

K0	= I0/(G+D);

* 2) Solve for the return to a unit of Capital

RK0	= KS0/K0;

* 3) Solve for the interest rate 

R	= (RK0 - D)/(1 + RK0 - D);

* 4) Find the price of a unit of Benchmark Capital

PK0	= 1/(1-r);

* 5) Set the steady-state Reference Paths

QREF(T)	=	(1+G)**(ORD(T)-1);
PREF(T) =	(1-r)**(ORD(T)-1);

Display K0,RK0,R,PK0,QREF,PREF;

Parameter TAX(T) Tax Rate on Investment;
TAX(T)	=  0;



$ONTEXT

$MODEL:DYN1

$Sectors:
	Y(T)	! Production
	I(T)	! Investment
	K(T)	! Capital Stock
	C(T)	! Consumption Index

$Commodities:
	PY(T)	! Price index on output
	PC(T)	! Price index on consumption
	RK(T)	! Present Value Return to capital
	PL(T)	! Present Value Wage
	PK(T)	! Price index on Capital
	PKT	! Price of Terminal Capital

$Consumers:
	RA	! Representative agent

$Auxiliary:
	TCAP	! Terminal Capital Demand

$PROD:Y(T) s:1
	O:PY(T)	Q:Y0
	I:RK(T)	Q:KS0
	I:PL(T)	Q:LS0

$PROD:I(T)
	O:PKT$TLAST(T)	Q:I0
	O:PK(T+1)	Q:I0
	I:PY(T)		Q:I0	A:RA T:TAX(T)

$PROD:K(T)
	O:PKT$TLAST(T)	Q:(K0*(1-D))
	O:PK(T+1)	Q:(K0*(1-D))
	O:RK(T)		Q:KS0
	I:PK(T)		Q:K0

$PROD:C(T) 
	O:PC(T)		Q:C0
	I:PY(T)		Q:C0

$Demand:RA s:0.5
	D:PC(T)		Q:(C0*QREF(T))	P:PREF(T)
	E:PK(TFIRST)	Q:K0
	E:PL(T)		Q:(LS0*QREF(T))
	E:PKT		Q:(-1)		R:TCAP

$Constraint:TCAP
	SUM(T$TLAST(T+1),C(T)*I(T+1) - I(T)*C(T+1)) =E= 0;

$offtext
$sysinclude mpsgeset dyn1

* Set Steady State Level Values:
Y.L(T)	= QREF(T);
I.L(T)	= QREF(T);
K.L(T)	= QREF(T);
C.L(T)	= QREF(T);

PY.L(T) = PREF(T);
PC.L(T) = PREF(T);
RK.L(T) = PREF(T);
PL.L(T) = PREF(T);
PK.L(T) = PREF(T)*PK0;


PKT.L	= (1-R)*SUM(TLAST, PK.L(TLAST));
TCAP.L  = SUM(TLAST, I0*QREF(TLAST)+K0*(1-D)*QREF(TLAST));

dyn1.ITERLIM = 0;
$include dyn1.gen
solve dyn1 using mcp;

* Run a simple experiment: 10% tax on investment from 2010 on.

SET 
	SHOCK(T) Time periods under the shock /2010*2100/;

TAX(SHOCK) = 0.1;

dyn1.ITERLIM = 100;
$include dyn1.gen
solve dyn1 using mcp;

Parameter MACRO(T,*)	"% Change in Macro Variables";
MACRO(T,"INVEST")	= 100*(I.L(T)/QREF(T)-1);
MACRO(T,"CONS")		= 100*(C.L(T)/QREF(T)-1);
MACRO(T,"CAPITAL")	= 100*(K.L(T)/QREF(T)-1);
MACRO(T,"OUTPUT")	= 100*(Y.L(T)/QREF(T)-1);


*	Plotting if you have gnuplot
*$libinclude gnuplot macro


Display macro;