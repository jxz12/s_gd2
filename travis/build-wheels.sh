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
    "${PYBIN}/pip" install build
    "${PYBIN}/python" -m build -w -o "${TMP_DIR}" -C--prefer-binary .

    # Bundle external shared libraries into the wheels
    ls -lrt $TMP_DIR
    for whl in $(ls -1 ${TMP_DIR}/*.whl); do
      if [[ "$whl" == *"s_gd2"* ]]; then
        auditwheel repair --plat "$PLAT" -w "${REPAIR_DIR}" $whl 
      else
        # don't need to repair dependency wheels
        mv $whl "${REPAIR_DIR}"
      fi
    done

    # Install and test
    "${PYBIN}/pip" install $PKG_NAME --no-index -f "${REPAIR_DIR}"
    cd $SOURCE_DIR
    "${PYBIN}/pip" install --prefer-binary .[test]
    "${PYBIN}/python" setup.py test
    cd ..

    # Clean up
    "${PYBIN}/pip" uninstall -y $PKG_NAME
    rm -rf build
    
    # Move wheel to output directory
    for whl in $(ls -1 ${REPAIR_DIR} | grep $PKG_NAME); do
      mv ${REPAIR_DIR}/$whl "${OUT_DIR}"
    done
done

