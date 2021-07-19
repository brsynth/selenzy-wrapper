include ../../extras/.env
include ../.env

SHELL := /bin/bash

# HELP
# This will output the help for each task
# thanks to https://marmelab.com/blog/2016/02/29/auto-documented-makefile.html
.PHONY: help test

help:
	@awk 'BEGIN {FS = ":.*?## "} /^[a-zA-Z_-]+:.*?## / {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}' $(MAKEFILE_LIST)

.DEFAULT_GOAL := all

MAKE_CMD = $(MAKE) -s --no-print-directory
ECHO = echo -n ">>>"


doc:
	@$(ECHO) "Building documentation...\n"
	cd ../../docs ; \
	sphinx-apidoc -f -o source ../${PACKAGE} \
	&& echo OK

doc-init:
	@$(ECHO) "Initialising documentation...\n"
	cd ../../docs ; \
	sphinx-init \
	&& echo OK
