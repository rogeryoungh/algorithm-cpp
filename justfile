test:
	oj-verify run

debug:
	oj-verify --config-file .verify-helper/config-debug.toml run

expand src: pre-expand
	python3 scripts/expand.py {{src}}

pre-expand:
	python3 scripts/generate-system-header.py

fmt:
	clang-format -i src/*.hpp src/*/*.hpp test/*/*.cpp
