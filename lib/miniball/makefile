TARGET = libminiball.so

OBJS_CPP = src/Miniball.o src/Miniball_wrapper.o

all: $(TARGET)


$(TARGET): $(OBJS_CPP)
	ar rcs $@ $(OBJS_CPP)


%.o: %.cpp
	g++ -fPIC -o $@ -c $<

%.o: %.c
	gcc -fPIC -o $@ -c $<

install:
	@cp $(TARGET) ../../lib/

clean:
	@rm -f $(TARGET) $(OBJS_CPP)
