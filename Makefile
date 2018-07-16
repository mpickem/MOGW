all:
	cd src; make all
	python -m compileall gw.py

clean:
	cd src; make clean
	rm -rf *.pyc
