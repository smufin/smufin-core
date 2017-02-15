TEST = $(sort $(wildcard */*.test))
RUN = $(TEST:.test=.run)

all: $(RUN)

%.run: %.test
	@./test $<
