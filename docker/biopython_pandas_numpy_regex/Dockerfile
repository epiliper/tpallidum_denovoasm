FROM python:3.11-slim
 
WORKDIR /app
 
RUN apt-get update && apt-get install -y --no-install-recommends \
    procps \
    && rm -rf /var/lib/apt/lists/*
 
RUN pip install --no-cache-dir \
    biopython==1.87 \
    pandas==3.0.2 \
    numpy==2.4.4 \
    regex==2026.4.4
 
CMD ["python3"]
 
