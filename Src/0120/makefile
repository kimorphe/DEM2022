all: gen_sheet  dem

gen_sheet: gen_sheet.o
	g++ -o gen_sheet gen_sheet.cpp
	mv gen_sheet  ../

gen_sheet.o: gen_sheet.cpp
	g++ -c gen_sheet.cpp

dem: dem.o force.o control.o prtcl2d.o sheet.o vec2.o spline.o uhyd_smec.o
	g++ -o dem dem.o force.o control.o prtcl2d.o sheet.o vec2.o spline.o uhyd_smec.o
	mv dem  ../../
dem.o: dem.cpp
	g++ -c dem.cpp
force.o: force.cpp
	g++ -c force.cpp
control.o: control.cpp
	g++ -c control.cpp
prtcl2d.o: prtcl2d.cpp
	g++ -c prtcl2d.cpp
sheet.o: sheet.cpp
	g++ -c sheet.cpp
vec2.o: vec2.cpp
	g++ -c vec2.cpp
spline.o: spline.cpp
	g++ -c spline.cpp

uhyd_smec.o: uhyd_smec.cpp
	g++ -c uhyd_smec.cpp


