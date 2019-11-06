set -e
yum -y remove cmake
/opt/python/${PYTHON_VERSION}/bin/pip install -r requirements.txt
export PATH=$PATH:/opt/python/${PYTHON_VERSION}/bin/
/opt/python/${PYTHON_VERSION}/bin/pip wheel .
auditwheel repair eqlib-*.whl
