CC = gcc
TARGET := bin/main
SRCDIR := src
BUILDDIR := build

SRCEXT := c
SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
INCLUDE := -I include

ERROR_CFLAGS = -Wall
OPTI_FLAG = -O0
CFLAGS = -lm
FLAGS = $(ERROR_CFLAGS) $(OPTI_FLAG) 

$(TARGET): $(OBJECTS)
	@echo " Linking..."
	@echo " $(CC) $^ -o $(TARGET)"; $(CC) $^ -o  $(TARGET) $(INCLUDE) $(CFLAGS)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(BUILDDIR)
	@echo " $(CC) $(FLAGS) -c -o $@ $<"; $(CC) $(FLAGS) $(INCLUDE) -c -o $@ $< 
	@echo " $(CC) $^ -o $(TARGET)"; $(CC) $^ -o $(TARGET) $(INCLUDE) 

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(BUILDDIR)
	@echo " $(CC) $(FLAGS) -c -o $@ $<"; $(CC) $(FLAGS) $(INCLUDE) -c -o $@ $< 

run:
	./$(TARGET)

edit:
	subl ./src/*.c
	subl ./include/*.h

clean:
	@echo " Cleaning...";
	@echo " $(RM) -r $(BUILDDIR) $(TARGET)"; $(RM) -r $(BUILDDIR) $(TARGET)

.PHONY: clean
