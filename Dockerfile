#Use the official Jupyter Data Science Notebook image with pinned Python version
FROM jupyter/datascience-notebook:latest

# Set working directory inside container
WORKDIR /usr/jovyan/bipartition_covers
ENV WORKDIR=/usr/jovyan/bipartition_covers

# Copy requirements and install dependencies
USER root
COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt
   
# Switch to non-root user (recommended by Jupyter)
USER jovyan

# Copy notebooks folder
COPY notebooks/ notebooks/
COPY README.md README.md

# Expose Jupyter port
EXPOSE 8888

# Start JupyterLab without token/password (local reproducibility)
CMD ["start-notebook.py", "--ServerApp.token=''", "--ServerApp.password=''"]
