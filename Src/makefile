all: gen_sheet dem mkUhyd RHmap

gen_sheet: gen_sheet.o
	g++ -o gen_sheet gen_sheet.cpp
gen_sheet.o: gen_sheet.cpp
	g++ -c gen_sheet.cpp

# SWELLING LAW DATA
mkUhyd: Uhyd.cpp
	g++ -o mkUhyd Uhyd.cpp
RHmap: RHmap.cpp
	g++ -o RHmap RHmap.cpp

# CG-MD PROGRAM
dem: dem.o force.o control.o prtcl2d.o sheet.o vec2.o spline.o uhyd_smec.o
	g++ -o dem dem.o force.o control.o prtcl2d.o sheet.o vec2.o spline.o uhyd_smec.o
	mv dem  ../
dem.o: dem.cpp dem.h
	g++ -c dem.cpp
force.o: force.cpp dem.h
	g++ -c force.cpp
control.o: control.cpp dem.h
	g++ -c control.cpp
prtcl2d.o: prtcl2d.cpp dem.h
	g++ -c prtcl2d.cpp
sheet.o: sheet.cpp dem.h
	g++ -c sheet.cpp
vec2.o: vec2.cpp dem.h
	g++ -c vec2.cpp
spline.o: spline.cpp dem.h
	g++ -c spline.cpp

uhyd_smec.o: uhyd_smec.cpp dem.h
	g++ -c uhyd_smec.cpp


