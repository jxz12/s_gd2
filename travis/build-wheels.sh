#!/bin/bash
set -e -x

# Auditwheel requirements
yum install -y atlas-devel

SOURCE_DIR="/io/cpp"
PKG_NAME="s_gd2"

OUT_DIR="${SOURCE_DIR}/dist/"
mkdir -p "${OUT_DIR}"
for PYBIN in /opt/python/*/bin; do
    # Select python version corresponding to this test
    if [ $("${PYBIN}/python" --version 2>&1 | grep -c "Python ${PYTHON}") -eq 0 ]; then
        continue
    fi

    # Setup
    TMP_DIR="${SOURCE_DIR}/wheelhouse_tmp/${PLAT}/${PYBIN}"
    REPAIR_DIR="${SOURCE_DIR}/wheelhouse_repair/${PLAT}/${PYBIN}"
    mkdir -p $TMP_DIR
    mkdir -p $REPAIR_DIR

    # Compile wheels
    cd $SOURCE_DIR
    "${PYBIN}/pip" install --upgrade pip
    "${PYBIN}/pip" install --prefer-binary numpy || continue
    "${PYBIN}/pip" wheel --prefer-binary . -w "${TMP_DIR}"

    # Bundle external shared libraries into the wheels
    ls -lrt $TMP_DIR
    for whl in $(ls -1 ${TMP_DIR}/s_gd2*.whl); do
      auditwheel repair --plat "$PLAT" -w "${REPAIR_DIR}" $whl 
    done

    cd $SOURCE_DIR
    # install deps
    "${PYBIN}/pip" install --prefer-binary .[test]
    # install s_gd2 directly from wheel
    "${PYBIN}/pip" uninstall -y $PKG_NAME
    "${PYBIN}/pip" install $PKG_NAME --no-index -f "${REPAIR_DIR}"
    # run tests
    "${PYBIN}/python" setup.py test
    cd ..

    # Clean up
    "${PYBIN}/pip" uninstall -y $PKG_NAME
    "${PYBIN}/pip" uninstall -y numpy
    rm -rf build
    
    # Move wheel to output directory
    for whl in $(ls -1 ${REPAIR_DIR} | grep $PKG_NAME); do
      mv ${REPAIR_DIR}/$whl "${OUT_DIR}"
    done
done

