#!/bin/bash

git ls-files -z '*.cpp' '*.cc' '*.cxx' '*.h' '*.hpp' '*.hh' | xargs -0 clang-format -i
