TARGET=NonLocalStaticAnalysis

CC=icpc

HEADDIR = include
SRCDIR  = src
OBJDIR  = obj
BINDIR  = bin

SOURCES  := $(wildcard $(SRCDIR)/*.cpp)
INCLUDES := $(wildcard $(HEADDIR)/*.h)
OBJECTS  := $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)

CFLAGS=-Wall -I$(HEADDIR) -c -w -O2 -m64

$(BINDIR)/$(TARGET): $(OBJECTS)
	$(CC) -o $@ $(LFLAGS) $(OBJECTS)
	@echo "Linking complete!"

$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	@$(CC) $(CFLAGS) $< -o $@
	@echo "Compiled "$<" successfully!"

clean:
	rm -rf $(OBJDIR)/*.o $(BINDIR)/NonLocalStaticAnalysis
