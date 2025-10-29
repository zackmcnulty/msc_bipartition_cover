# Bipartition Covers in the Multispecies Coalescent Model

This repository contains notebooks and code for simulating and analyzing bipartition covers.  
It uses a **Docker-based JupyterLab environment** for reproducibility. The base image `upyter/datascience-notebook` includes the standard Python data science stack, and a few additional packages:

- `dendropy`
- `tqdm`
- `numba`

You can of course simply clone the repository and run the code on your local kernel if you prefer. 

## Prerequisites

1. Install [Docker](https://www.docker.com/get-started) for your platform.
2. Install [Docker Compose](https://docs.docker.com/compose/install/).

---

## Clone the Repository

First, we need to download the repository from GitHub. 

```bash
git clone https://github.com/zackmcnulty/msc_bipartition_cover.git
cd msc_bipartition_cover
```

## Starting the Container

Now we start the container and get the jupyterlab running. 

```bash
docker-compose up --build
```

## Opening Jupyter Lab

This should start a Jupyter server on your local device. Open any web-browser and navigate to the URL `http://localhost:8888/lab` or click [here](http://localhost:8888/lab). 

