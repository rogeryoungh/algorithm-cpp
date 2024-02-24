test:
	oj-verify run

debug:
	oj-verify --config-file .verify-helper/config-debug.toml run

expand src: pre-expand
	python3 scripts/expand.py {{src}}

pre-expand:
	python3 scripts/generate-system-header.py

fmt:
	clang-format -i include/algo/*.hpp include/algo/*/*.hpp test/*/*.cpp
