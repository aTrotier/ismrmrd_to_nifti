FROM python:3
# ismrmrd needs to be installed _after_ pyxb
RUN apt-get update && apt-get install -y libhdf5-dev python-numpy python-scipy && pip install nibabel h5py pyxb && pip install ismrmrd
WORKDIR /opt/src/python_version
COPY . .
ENTRYPOINT [ "/opt/src/python_version/main_conversion.py" ]
