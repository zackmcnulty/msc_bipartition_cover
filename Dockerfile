# Use the official Jupyter Data Science Notebook image with pinned Python version
FROM jupyter/datascience-notebook:python-3.12.1

# Set working directory inside container
WORKDIR /usr/jovyan/bipartition_covers

# Copy requirements and install dependencies
COPY requirements.txt ./
RUN pip install --upgrade pip && \
    pip install --no-cache-dir -r requirements.txt

# Copy notebooks folder
COPY notebooks/ notebooks/
COPY README.md README.md

# Switch to non-root user (recommended by Jupyter)
USER jovyan

# Expose Jupyter port
EXPOSE 8888

# Start JupyterLab without token/password (local reproducibility)
CMD ["start-notebook.sh", "--NotebookApp.token=''", "--NotebookApp.password=''"]
