
optimization = -O3 
objects = ehhpop.o CmdLine.o 
libs = -L/opt/local/lib -I/usr/include -lgsl -lgslcblas -lm

xpehh: xpehh.o $(objects)
	g++ $(libs) -o xpehh xpehh.o $(objects) $(optimization)

hzy: hzy.o $(objects)
	g++ $(libs) -o hzy hzy.o $(objects) $(optimization)
 
test: test.o $(objects)
	g++ $(libs) -o test test.o $(objects) $(optimization)
 
%.o: %.cpp
	g++ -c $< -o $@ $(optimization)

clean:
	rm *.o mpop
