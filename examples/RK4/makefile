#specify the directory of the vectorclass and the ParallelDDE libraries
vcldir = /home/user/VCL/version2-2.01.02/
ddedir = /home/user/Paralell-DDE-Solver/

#specify the compiler flags
flags = -O3 -std=c++17 -mavx

all: example1 example2 example3 example4
.PHONY: all

dir:
	@echo "VCL directory		: "$(vcldir)
	@echo "ParallelDDE directory	: "$(ddedir)
example1:
	@g++ $(flags) -I$(vcldir) -I$(ddedir) -o example1 example1_basics.cpp
	@echo "Compiled example1_basics.cpp to example1"
example2:
	@g++ $(flags) -I$(vcldir) -I$(ddedir) -o example2 example2_discontinuities.cpp
	@echo "Compiled example2_discontinuities.cpp to example2"
example3:
	@g++ $(flags) -I$(vcldir) -I$(ddedir) -o example3 example3_parametersweep.cpp
	@echo "Compiled example3_parametersweep.cpp to example3"
example4:
	@g++ $(flags) -I$(vcldir) -I$(ddedir) -o example4 example4.cpp
	@echo "Compiled example4.cpp to example4"
clean:
	@rm -f example1 example2 example3 example4