#!/bin/bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
SOURCE_DIR="$(dirname $(dirname $SCRIPT_DIR))/nextflow"

# Builds fake Java jars
cwd=${PWD}
cd ${SCRIPT_DIR}/java

javac FakeExtractionPipeline.java
jar cfe extraction.jar FakeExtractionPipeline FakeExtractionPipeline.class

javac FakeLoadingPipeline.java
jar cfe loading.jar FakeLoadingPipeline FakeLoadingPipeline.class

javac FakeClusteringPipeline.java
jar cfe clustering.jar FakeClusteringPipeline FakeClusteringPipeline.class

cd ${cwd}
