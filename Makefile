all: build/thermo.o build/environment.o build/parcel.o build/pseudo.o build/RK_dynamic.o build/FD_dynamic.o | output
	g++ -O3 build/RK_dynamic.o build/FD_dynamic.o build/pseudo.o build/environment.o build/thermo.o build/parcel.o src/main.cpp -o simulator.exe
	rm -rf build

build/thermo.o: src/thermodynamic_calc.cpp src/thermodynamic_calc.h | build
	g++ -O3 -c src/thermodynamic_calc.cpp -o build/thermo.o

build/environment.o: src/environment.cpp src/environment.h | build
	g++ -O3 -c src/environment.cpp -o build/environment.o
	
build/parcel.o: src/parcel.cpp src/parcel.h | build
	g++ -O3 -c src/parcel.cpp -o build/parcel.o
	
build/pseudo.o: src/pseudoadiabatic_scheme.cpp src/pseudoadiabatic_scheme.h | build
	g++ -O3 -c src/pseudoadiabatic_scheme.cpp -o build/pseudo.o
	
build/RK_dynamic.o: src/runge_kutta_dynamics.cpp src/dynamic_scheme.h | build
	g++ -O3 -c src/runge_kutta_dynamics.cpp -o build/RK_dynamic.o

build/FD_dynamic.o: src/finite_difference_dynamics.cpp src/dynamic_scheme.h | build
	g++ -O3 -c src/finite_difference_dynamics.cpp -o build/FD_dynamic.o

build:
	mkdir build
	
output:
	mkdir output
	
clean:
	rm -rf build
	rm -f *.exe
