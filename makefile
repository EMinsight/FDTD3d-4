OBJS = calc_fdtd.o D_update.o D_update_pml.o E_update.o H_update.o H_update_pml.o \
		main.o memory_allocate3d.o pml_class.o PML_field_initialize.o PML_idx_initialize.o \
		sigma_calc.o
HEADERS = fdtd3d.h main.h pml.h
OPTS = -I/opt/include/eigen3 -std=c++1z -O3

main : $(OBJS)
	g++ -o $@ $(OBJS)
%.o : %.cpp $(HEADERS) $(OPTS)
	g++ -c $< $(HEADERS) $(OPTS)