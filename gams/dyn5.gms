$TITLE This is a ramsey model with multiple sectors


$ontext

This model extends \DYN4.GMS to 
incorporate sector specific growth 
projections. Non-balanced growth.

Assume that sector 1 grows at 2% over the horizon, but
sector 2 grows at 2.2% for the first 50yrs and at 1.6% 
for the last 50yrs.  The assumption is that the 
unbalanced growth is due to capital's productivity
changing over time. 

Original Base Year Data (rectangular SAM)


	   Y	X1	X2	  I	|	  RA
----------------------------------------------------
PY	 200			-20	|	-180
PX1	-120	120			|	
PX2	 -80		80		|	
RENT	        -20    -30		|	  50
WAGE	       -100    -50		|	 150
Savings			  	 20	|	 -20


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

Set J	Index on Production Sectors /1,2/;
Alias (J,JJ);

Parameter
	X0(J)	 Output by sector,
	LD(J)	 Labor Demand by Sector,
	KD(J)	 Capital Demand by Sector;

X0("1")	=	120;
X0("2")	=	 80;
LD("1")	=	100;
LD("2")	=	 50;
KD("1")	=	 20;
KD("2")	=	 30;

*----------Setup the Dynamics----------*
Set	T		Time Periods	/2000*2100/
	TFIRST(T)	First Period	/2000/
	TLAST(T)	Terminal Period /2100/
	FIRST50(T)	First 50 yrs	/2000*2050/
	LAST50(T)	Last 50 yrs	/2051*2100/;

*   Introduce the unbalanced growth paths

Parameter 
	PHI(J,T) Productivity Shift Parameter,
	QREFJ(J,T) Sectoral growth path,
	SECQREF(J,T) Sectoral growth path,
	AVGREF(T) Average reference path;

SECQREF("1",T) = (1.02)**(ORD(T)-1);
SECQREF("2",FIRST50) = (1.022)**(ORD(FIRST50)-1);
SECQREF("2",LAST50) = SECQREF("2","2050")*((1.016)**(ORD(LAST50)));

AVGREF(T) = SUM(J,SECQREF(J,T)*X0(J)/SUM(JJ,X0(J)));

*  If you want to see the paths
*Parameter graph(T,*);
*graph(T,J)=secqref(j,T);
*graph(T,"avg")=avgref(T);
*$libinclude gnuplot graph

*   Setup balanced growth as a starting point
*   This insures that we have the dynamics right

Scalars
	R	Rate of Interest	/0.05/
	G	Growth Rate		/0.02/,
	D	Depreciation Rate	/0.02/;

*   Adjust G to be the wtd. average rate from 
*   the sector specific projections

G = ((AVGREF("2100")/AVGREF("2000"))**(1/100))-1;

Display G;

Parameter
	K0	Benchmark Capital Stock,
	RK0	Benchmark Return on a Unit of K
	PK0	Benchmark Price of a Unit of K
	QREF(T) Growth Path for Quantities
	PREF(T) Present Value Price Paths
	KSCAL0   Calibrated level of Value to Capital Services;

* 1) Given the interest rate Calculate PK0 and RK0

PK0	= 1/(1-r);
RK0	= (R-R*D+D)/(1-R);

* 2) Solve for the Initial Capital Stock

K0	= I0/(G+D);

* 3) Find the Calibrated Flow to Capital

KSCAL0	=K0*RK0;

* 4) Adjust Social Accounts to match the new level
*     of capital services.  Here we employ a 
*     weighted-least-squares routine.  You may
*     want to use an alternative objective funct.?

POSITIVE VARIABLES 
	VK(J)	calibrated value of capital earnings
	VL(J)	calibrated value of labor earnings;

VARIABLE	OBJ	objective function;

EQUATIONS
	VABAL(J)	value-added consistency by sector
	VKBAL		capital income balance
	OBJDEF		define the objective function;

PARAMETER	VA0(J)		Base year value added;

*	sectoral value-added

VA0(J) = LD(J) + KD(J);

*  Value of capital earnings consistent with the steady-state
*  growth path:

VABAL(J)..	VK(J) + VL(J)  =E= VA0(J);

VKBAL..		SUM(J, VK(J)) =E= KSCAL0;

OBJDEF..	OBJ =E= SUM((J), (1/KD(J))*SQR(VK(J) - KD(J)));

MODEL KBAL / VABAL, VKBAL, OBJDEF /;

KBAL.ITERLIM = 1000;
SOLVE KBAL USING NLP MINIMIZING OBJ;

*  Show what we did in the listing file.
DISPLAY "++++++++++++++++Raw factor Demands++++++++++++++++",LD, KD;

*	New Calibrated labor and capital demands:
LD(J) = VL.L(J);  
KD(J) = VK.L(J);

DISPLAY "+++++++++++++Calibrated Factor Demands++++++++++++",LD, KD;


* 5) Set the steady-state Reference Paths

QREF(T)	=	(1+G)**(ORD(T)-1);
PREF(T) =	(1-r)**(ORD(T)-1);
QREFJ(J,T) =	QREF(T);
PHI(J,T) =	1;

Display K0,RK0,R,PK0,QREF,PREF;

Parameter TAX(T) Tax Rate on Capital Earnings;
TAX(T)	=  0;

$ONTEXT

$MODEL:DYN5

$Sectors:
	Y(T)	! Macro Output (transitory utility)
	X(J,T)	! Production
	KJ(J,T)	! Capital Input with productivity adj
	I(T)	! Investment
	K(T)	! Capital Stock
	C(T)	! Consumption Index

$Commodities:
	PY(T)	! Price index on macro output
	PX(J,T)	! Price index on sector output
	PKJ(J,T) ! Price index on Capital with prod adj
	PC(T)	! Price index on consumption
	RK(T)	! Present Value Return to capital
	PL(T)	! Present Value Wage
	PK(T)	! Price index on Capital
	PKT	! Price of Terminal Capital

$Consumers:
	RA	! Representative agent

$Auxiliary:
	TCAP	! Terminal Capital Demand
	SK(J,T)	! Endogenous Subsidy to find productivity shock
	ADJ_PKJ(J,T) ! Adjustment to negate income effects of SK

$PROD:Y(T) s:1
	O:PY(T)	  Q:Y0
	I:PX(J,T) Q:X0(J)

$PROD:X(J,T) s:1
	O:PX(J,T) Q:X0(J)
	I:PKJ(J,T) Q:KD(J)
	I:PL(T)	  Q:LD(J)

$PROD:KJ(J,T)
	O:PKJ(J,T)	Q:(PHI(J,T)*KD(J)) A:RA N:SK(J,T) M:-1
	I:RK(T)		Q:KD(J)

$PROD:I(T)
	O:PKT$TLAST(T)	Q:I0
	O:PK(T+1)	Q:I0
	I:PY(T)		Q:I0	A:RA T:TAX(T)

$PROD:K(T)
	O:PKT$TLAST(T)	Q:(K0*(1-D))
	O:PK(T+1)	Q:(K0*(1-D))
	O:RK(T)		Q:KSCAL0
	I:PK(T)		Q:K0

$PROD:C(T) 
	O:PC(T)		Q:C0
	I:PY(T)		Q:C0

$Demand:RA s:0.5
	D:PC(T)		Q:(C0*QREF(T))	P:PREF(T)
	E:PK(TFIRST)	Q:K0
	E:PL(T)		Q:(SUM(j,LD(j))*QREF(T))
	E:PKT		Q:(-1)	R:TCAP
	E:PKJ(J,T)	Q:1	R:ADJ_PKJ(J,T)

$Constraint:SK(J,T)
	X(J,T) - QREFJ(J,T) =E= 0;
$Constraint:ADJ_PKJ(J,T)
	ADJ_PKJ(J,T) =E= KD(J)*SK(J,T)*KJ(J,T);

$Constraint:TCAP
	SUM(T$TLAST(T+1),C(T)*I(T+1) - I(T)*C(T+1)) =E= 0;

$offtext
$sysinclude mpsgeset DYN5

* Set Steady State Level Values:
Y.L(T)	= QREF(T);
X.L(J,T)= QREF(T);
I.L(T)	= QREF(T);
KJ.L(J,T)= QREF(T);
K.L(T)	= QREF(T);
C.L(T)	= QREF(T);

PY.L(T) = PREF(T);
PX.L(J,T)= PREF(T);
PC.L(T) = PREF(T);
PKJ.L(J,T) = PREF(T);
RK.L(T) = PREF(T);
PL.L(T) = PREF(T);
PK.L(T) = PREF(T)*PK0;


PKT.L	= (1-R)*SUM(TLAST, PK.L(TLAST));
TCAP.L  = SUM(TLAST, I0*QREF(TLAST)+K0*(1-D)*QREF(TLAST));

*  Check the Balanced Growth Equilibrium
DYN5.ITERLIM = 0;
$include DYN5.gen
solve DYN5 using mcp;

*  Free bounds on adjustment constraints
ADJ_PKJ.LO(J,T) =-INF;
SK.LO(J,T)	=-0.99;

*  Use the solver to find the output coefficent
*  for sector specific growth paths
QREFJ(J,T) = SECQREF(J,T);

DYN5.ITERLIM = 200;
$include DYN5.gen
solve DYN5 using mcp;

*  Lock in the productivity assumption and check
*  That we have a unbalance equilibrium

PHI(J,T)	= 1+SK.L(J,T);
SK.FX(J,T)	= 0;
ADJ_PKJ.FX(J,T) = 0;

DYN5.ITERLIM = 0;
$include DYN5.gen
solve DYN5 using mcp;


* Store the Business as Usual level values
Parameter 
	BAUMACRO(T,*),
	BAUOUT(T,*);

BAUMACRO(T,"INVEST")	= I.L(T);
BAUMACRO(T,"CONS")	= C.L(T);
BAUMACRO(T,"CAPITAL")	= K.L(T);
BAUMACRO(T,"OUTPUT")	= Y.L(T);

BAUOUT(T,"Aggregate")	= Y.L(T);
BAUOUT(T,"Sector1")	= X.L("1",T);
BAUOUT(T,"Sector2")	= X.L("2",T);


* Run a simple experiment: 10% tax on investment from 2010 on.

SET 
	SHOCK(T) Time periods under the shock /2010*2100/;

TAX(SHOCK) = 0.1;

DYN5.ITERLIM = 100;
$include DYN5.gen
solve DYN5 using mcp;

Parameter MACRO(T,*)	"% Change in Macro Variables",
	  OUTPUT(T,*)	"% Change in output";
MACRO(T,"INVEST")	= 100*(I.L(T)/BAUMACRO(T,"INVEST")-1);
MACRO(T,"CONS")		= 100*(C.L(T)/BAUMACRO(T,"CONS")-1);
MACRO(T,"CAPITAL")	= 100*(K.L(T)/BAUMACRO(T,"CAPITAL")-1);
MACRO(T,"OUTPUT")	= 100*(Y.L(T)/BAUMACRO(T,"OUTPUT")-1);

OUTPUT(T,"Aggregate")	= 100*(Y.L(T)/BAUOUT(T,"AGGREGATE")-1);
OUTPUT(T,"Sector1")	= 100*(X.L("1",T)/BAUOUT(T,"Sector1")-1);
OUTPUT(T,"Sector2")	= 100*(X.L("2",T)/BAUOUT(T,"Sector2")-1);

$libinclude gnuplot macro
$libinclude gnuplot output
Display macro,output;