# Bipartition Covers in the Multispecies Coalescent Model

This repository contains jupyter notebooks and code for reproducing figures in the paper [PAPER NAME](link). It focuses on exploring various upperbounds the paper develops on the number of gene trees required to achieve a bipartion cover of the species with high probability, comparing them to empirically estimated quantities under a variety of different species tree topologies. 

It uses a **Docker-based JupyterLab environment** for reproducibility, although you can of course simply clone the repository and run the code on your local kernel if you prefer. The former uses the base Docker image `jupyter/datascience-notebook` which includes the standard Python data science stack, and a few additional packages described in `requirements.txt`. Most notably, it uses `dendropy` for simulating gene trees under the multispecies coalescent model

## Prerequisites

If you would like to use Docker, you can install it following locations if you do not have it already.

1. Install [Docker](https://www.docker.com/get-started) for your platform.
2. Install [Docker Compose](https://docs.docker.com/compose/install/).

---

Then you can either directly download the Docker image [HERE](todo) or follow the below steps to build a local copy. 

## Building Docker Image Locally

First, we need to download the repository from GitHub. 

```bash
git clone https://github.com/zackmcnulty/msc_bipartition_cover.git
cd msc_bipartition_cover
```

Now we can build local Docker image. 

```bash
docker-compose up --build
```

This should start a Jupyter server on your local device. Open any web-browser and navigate to the URL `http://localhost:8888/lab` or simply click [here](http://localhost:8888/lab) to navigate to Jupyter lab. 

All the code for generating the relevant plots can be found in the notebook `notebooks/bipartitions_paper.ipynb`. Helper functions and their implementation details can be found in the `utility` folder.





