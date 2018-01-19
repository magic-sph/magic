#!/bin/bash
## launch unit tests locally
## interactive mode
python -u magic_wizard.py --level=-1 2>&1 | tee log.test
## quiet mode
#python -u magic_wizard.py --level=-1 >log.test 2>&1
