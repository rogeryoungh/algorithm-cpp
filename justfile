test:
	oj-verify run

debug:
	oj-verify --config-file .verify-helper/config-debug.toml run

expand src: pre-expand
	python3 scripts/expand.py {{src}}

docs:
	cd docs/ && hugo serve

pre-expand:
	python3 scripts/pre-expand.py

bench src:
	g++ {{src}} -lbenchmark -std=c++20 -o a.out -O3 -fno-split-paths -march=native -Wno-ignored-attributes
	./a.out

fmt:
	clang-format -i src/*.hpp src/*/*.hpp src/*/*/*.hpp test/*/*.cpp

# docs:
# 	oj-verify docs

# expand2 src:
# 	oj-bundle {{src}}
