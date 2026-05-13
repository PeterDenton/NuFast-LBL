CXX=g++

CFlags=-c -Wall -MMD -std=c++17
CFlags+=-Ofast -ffast-math # can improve speed by 10s of percent over the default of -O0

Sources=$(wildcard src/*.cpp)
IncludeDir=-I./include
AllObjects=$(addprefix obj/,$(notdir $(Sources:.cpp=.o)))
Executables=Benchmark
Objects=$(filter-out $(addprefix obj/,$(Executables:=.o)),$(AllObjects))

all: $(Sources) $(Executables)

$(Executables): $(AllObjects)
	@mkdir -p obj data
	$(CXX) $(Objects) $(addprefix obj/,$@.o) -o $@

obj/%.o: src/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CFlags) $(IncludeDir) $< -o $@

-include $(AllObjects:.o=.d)

clean:
	rm -rf obj/*.o obj/*.d $(Executables)

