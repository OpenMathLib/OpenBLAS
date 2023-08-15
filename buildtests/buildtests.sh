#!/bin/bash

#Remove test directory if present, then make a new one
rm -r ../../OpenBLAS-buildtests
mkdir ../../OpenBLAS-buildtests

#Store path to current directory for later use
startpath=$(pwd)

#First do a build using the default settings
mkdir ../../OpenBLAS-buildtests/default
cp -r ../* ../../OpenBLAS-buildtests/default/
cd ../../OpenBLAS-buildtests/default/
make
if [ $? -ne 0 ]; then
  echo "TEST ERROR: build failed"
  exit -127
fi
default_hash=$(shasum libopenblas.so)
cd "$startpath"

#Manual target should yield the same binary as the default
mkdir ../../OpenBLAS-buildtests/manual_target
cp -r ../* ../../OpenBLAS-buildtests/manual_target/
cp Makefile.rule_manual_target ../../OpenBLAS-buildtests/manual_target/Makefile.rule
cd ../../OpenBLAS-buildtests/manual_target/
make
if [ $? -ne 0 ]; then
  echo "TEST ERROR: build failed"
  exit -127
fi
manual_target_hash=$(shasum libopenblas.so)
cd "$startpath"

#DYNAMIC_ARCH = 0 should yield the same binary as the default
mkdir ../../OpenBLAS-buildtests/dynarch_disabled
cp -r ../* ../../OpenBLAS-buildtests/dynarch_disabled/
cp Makefile.rule_dynarch_disabled ../../OpenBLAS-buildtests/dynarch_disabled/Makefile.rule
cd ../../OpenBLAS-buildtests/dynarch_disabled/
make
if [ $? -ne 0 ]; then
  echo "TEST ERROR: build failed"
  exit -127
fi
dynarch_disabled_hash=$(shasum libopenblas.so)
cd "$startpath"

#DYNAMIC_ARCH = 1 should yield a different binary
mkdir ../../OpenBLAS-buildtests/dynarch_enabled
cp -r ../* ../../OpenBLAS-buildtests/dynarch_enabled/
cp Makefile.rule_dynarch_enabled ../../OpenBLAS-buildtests/dynarch_enabled/Makefile.rule
cd ../../OpenBLAS-buildtests/dynarch_enabled/
make
if [ $? -ne 0 ]; then
  echo "TEST ERROR: build failed"
  exit -127
fi
dynarch_enabled_hash=$(shasum libopenblas.so)
cd "$startpath"

#DYNAMIC_ARCH = 1 DYNAMIC_OLDER = 0 should be the same as DYNAMIC_ARCH = 1
mkdir ../../OpenBLAS-buildtests/dynarch_enabled_old_disabled
cp -r ../* ../../OpenBLAS-buildtests/dynarch_enabled_old_disabled/
cp Makefile.rule_dynarch_enabled_old_disabled ../../OpenBLAS-buildtests/dynarch_enabled_old_disabled/Makefile.rule
cd ../../OpenBLAS-buildtests/dynarch_enabled_old_disabled/
make
if [ $? -ne 0 ]; then
  echo "TEST ERROR: build failed"
  exit -127
fi
dynarch_enabled_old_disabled_hash=$(shasum libopenblas.so)
cd "$startpath"

#DYNAMIC_ARCH = 1 DYNAMIC_OLDER = 1 should be different
mkdir ../../OpenBLAS-buildtests/dynarch_enabled_old_enabled
cp -r ../* ../../OpenBLAS-buildtests/dynarch_enabled_old_enabled/
cp Makefile.rule_dynarch_enabled_old_enabled ../../OpenBLAS-buildtests/dynarch_enabled_old_enabled/Makefile.rule
cd ../../OpenBLAS-buildtests/dynarch_enabled_old_enabled/
make
if [ $? -ne 0 ]; then
  echo "TEST ERROR: build failed"
  exit -127
fi
dynarch_enabled_old_enabled_hash=$(shasum libopenblas.so)
cd "$startpath"

#DYNAMIC_OLDER = 1 alone should be ignored
mkdir ../../OpenBLAS-buildtests/old_enabled
cp -r ../* ../../OpenBLAS-buildtests/old_enabled/
cp Makefile.rule_old_enabled ../../OpenBLAS-buildtests/old_enabled/Makefile.rule
cd ../../OpenBLAS-buildtests/old_enabled/
make
if [ $? -ne 0 ]; then
  echo "TEST ERROR: build failed"
  exit -127
fi
old_enabled_hash=$(shasum libopenblas.so)
cd "$startpath"

#DYNAMIC_ARCH = 0 DYNAMIC_OLDER = 1 should be ignored
mkdir ../../OpenBLAS-buildtests/dynarch_disabled_old_enabled
cp -r ../* ../../OpenBLAS-buildtests/dynarch_disabled_old_enabled/
cp Makefile.rule_dynarch_disabled_old_enabled ../../OpenBLAS-buildtests/dynarch_disabled_old_enabled/Makefile.rule
cd ../../OpenBLAS-buildtests/dynarch_disabled_old_enabled/
make
if [ $? -ne 0 ]; then
  echo "TEST ERROR: build failed"
  exit -127
fi
dynarch_disabled_old_enabled_hash=$(shasum libopenblas.so)
cd "$startpath"

#BINARY=64 should yield the same binary as the default
mkdir ../../OpenBLAS-buildtests/bin64
cp -r ../* ../../OpenBLAS-buildtests/bin64/
cp Makefile.rule_bin64 ../../OpenBLAS-buildtests/bin64/Makefile.rule
cd ../../OpenBLAS-buildtests/bin64/
make
if [ $? -ne 0 ]; then
  echo "TEST ERROR: build failed"
  exit -127
fi
bin64_hash=$(shasum libopenblas.so)
cd "$startpath"

#BINARY=32 should be different
mkdir ../../OpenBLAS-buildtests/bin32
cp -r ../* ../../OpenBLAS-buildtests/bin32/
cp Makefile.rule_bin32 ../../OpenBLAS-buildtests/bin32/Makefile.rule
cd ../../OpenBLAS-buildtests/bin32/
make
if [ $? -ne 0 ]; then
  echo "TEST ERROR: build failed"
  exit -127
fi
bin32_hash=$(shasum libopenblas.so)
cd "$startpath"


echo "$default_hash"
echo "$manual_target_hash"
echo "$dynarch_disabled_hash"
echo "$dynarch_enabled_hash"
echo "$dynarch_enabled_old_disabled_hash"
echo "$dynarch_enabled_old_enabled_hash"
echo "$old_enabled_hash"
echo "$dynarch_disabled_old_enabled_hash"
echo "$bin64_hash"
echo "$bin32_hash"

if [ "$default_hash" != "$manual_target_hash" ]; then
  echo "TEST ERROR: manual target changes binary"
  exit -1
fi
if [ "$default_hash" != "$dynarch_disabled_hash" ]; then
  echo "TEST ERROR: DYNAMIC_ARCH = 0 changes binary"
  exit -2
fi
if [ "$default_hash" = "$dynarch_enabled_hash" ]; then
  echo "TEST ERROR: DYNAMIC_ARCH = 1 does not change binary"
  exit -3
fi
if [ "$dynarch_enabled_hash" != "$dynarch_enabled_old_disabled_hash" ]; then
  echo "TEST ERROR: DYNAMIC_ARCH = 1 is not the same as DYNAMIC_ARCH = 1 DYNAMIC_OLDER = 0"
  exit -4
fi
if [ "$default_hash" = "$dynarch_enabled_old_enabled_hash" ]; then
  echo "TEST ERROR: DYNAMIC_ARCH = 1 DYNAMIC_OLDER = 1 does not change binary"
  exit -5
fi
if [ "$dynarch_enabled_hash" = "$dynarch_enabled_old_enabled_hash" ]; then
  echo "TEST ERROR: DYNAMIC_ARCH = 1 is the same as DYNAMIC_ARCH = 1 DYNAMIC_OLDER = 1"
  exit -6
fi
if [ "$default_hash" != "$old_enabled_hash" ]; then
  echo "TEST ERROR: DYNAMIC_OLDER = 1 alone changes binary"
  exit -7
fi
if [ "$default_hash" != "$dynarch_disabled_old_enabled_hash" ]; then
  echo "TEST ERROR: DYNAMIC_ARCH = 0 DYNAMIC_OLDER = 1 changes binary"
  exit -8
fi
if [ "$default_hash" != "$bin64_hash" ]; then
  echo "TEST ERROR: BINARY=64 changes binary"
  exit -9
fi
if [ "$default_hash" = "$bin32_hash" ]; then
  echo "TEST ERROR: BINARY=32 does not change binary"
  exit -10
fi


echo "All build tests passed. Yay!"
exit 0