OBJ=QM.f90 chebft.o derivs.o dtVdiff1.o dtVdiff2.o eta.o expaPoint.o FRG.o fun.o gapEq.o GAULEG.o initial.o intLin.o newt.o nfdx.o odeint.o phypoint.o phypoint2.o PolyakovEq.o selfEQ.o vInf.o
F77=f95

exe: $(OBJ)
	$(F77) -o exe -O2 $(OBJ) $(LIBS)

chebft.o: chebft.f90
	$(F77)  -c -O2 chebft.f90

derivs.o: derivs.f90
	$(F77)  -c -O2 derivs.f90

dtVdiff1.o: dtVdiff1.f90
	$(F77)  -c -O2 dtVdiff1.f90

dtVdiff2.o: dtVdiff2.f90
	$(F77)  -c -O2 dtVdiff2.f90

eta.o: eta.f90
	$(F77)  -c -O2 eta.f90

expaPoint.o: expaPoint.f90
	$(F77)  -c -O2 expaPoint.f90

FRG.o: FRG.f90
	$(F77)  -c -O2 FRG.f90

fun.o: fun.f90
	$(F77)  -c -O2 fun.f90

gapEq.o: gapEq.f90
	$(F77)  -c -O2 gapEq.f90

GAULEG.o: GAULEG.f90
	$(F77)  -c -O2 GAULEG.f90

initial.o: initial.f90
	$(F77)  -c -O2 initial.f90

intLin.o: intLin.f90
	$(F77)  -c -O2 intLin.f90

newt.o: newt.f90
	$(F77)  -c -O2 newt.f90

nfdx.o: nfdx.f90
	$(F77)  -c -O2 nfdx.f90

odeint.o: odeint.f90
	$(F77)  -c -O2 odeint.f90

phypoint.o: phypoint.f90
	$(F77)  -c -O2 phypoint.f90

phypoint2.o: phypoint2.f90
	$(F77)  -c -O2 phypoint2.f90

PolyakovEq.o: PolyakovEq.f90
	$(F77)  -c -O2 PolyakovEq.f90

selfEQ.o: selfEQ.f90
	$(F77)  -c -O2 selfEQ.f90

vInf.o: vInf.f90
	$(F77)  -c -O2 vInf.f90
