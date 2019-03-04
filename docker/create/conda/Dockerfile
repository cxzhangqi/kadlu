# Use miniconda as parent image
FROM continuumio/miniconda3

# Set the working directory to /app
WORKDIR /app

# Copy kadlu package, setup and environment files to working directory
COPY kadlu/ /app/kadlu
COPY setup.py /app
COPY environment.yml /app/environment.yml

# Create conda environment
RUN conda env create -f /app/environment.yml

# Pull the environment name out of the environment.yml
RUN echo "source activate $(head -1 /tmp/environment.yml | cut -d' ' -f2)" > ~/.bashrc
ENV PATH /opt/conda/envs/$(head -1 /tmp/environment.yml | cut -d' ' -f2)/bin:$PATH

# Install the kadlu package using pip
RUN python setup.py sdist
RUN pip install dist/kadlu-0.0.1.tar.gz