all: clean create run

create:
	clear
	g++ -o render main.cpp -lm -lncurses

run:
	./render ./models/cube.obj

clean:
	rm -f render