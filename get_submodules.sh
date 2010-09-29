#!/bin/sh
git submodule init
git submodule update
(cd barvinok; git submodule init polylib isl-polylib; git submodule update)
