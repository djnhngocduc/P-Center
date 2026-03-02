all:
	$(MAKE) -C "/home/g22028166/P-Center/solvers/kissat/build"
kissat:
	$(MAKE) -C "/home/g22028166/P-Center/solvers/kissat/build" kissat
tissat:
	$(MAKE) -C "/home/g22028166/P-Center/solvers/kissat/build" tissat
clean:
	rm -f "/home/g22028166/P-Center/solvers/kissat"/makefile
	rm -f "/home/g22028166/P-Center/solvers/kissat"/src/makefile
	-$(MAKE) -C "/home/g22028166/P-Center/solvers/kissat/build" clean
	rm -rf "/home/g22028166/P-Center/solvers/kissat/build"
coverage:
	$(MAKE) -C "/home/g22028166/P-Center/solvers/kissat/build" coverage
format:
	$(MAKE) -C "/home/g22028166/P-Center/solvers/kissat/build" format
test:
	$(MAKE) -C "/home/g22028166/P-Center/solvers/kissat/build" test
.PHONY: all clean coverage format kissat test tissat
