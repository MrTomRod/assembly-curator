FROM ghcr.io/astral-sh/uv:python3.12-bookworm-slim
# Inspiration: https://github.com/astral-sh/uv-docker-example/blob/main/Dockerfile

# Run QT without display
ENV QT_QPA_PLATFORM=offscreen
ENV XDG_RUNTIME_DIR=/tmp/runtime-root

# Install libqt5svg5 for gfaviz
RUN apt-get update  && \
    apt-get install --no-install-recommends -y libqt5svg5 && \
    apt-get clean

# Copy gfaviz-mrtomrod from existing Docker image
COPY --from=docker.io/troder/gfaviz:1.0.0-mrtomrod /opt/gfaviz/gfaviz-mrtomrod /usr/local/bin/

# Copy minimap2 from existing Docker image
COPY --from=docker.io/nanozoo/minimap2:2.28--9e3bd01 /opt/conda/bin/minimap2 /usr/local/bin/

WORKDIR /app

# Enable bytecode compilation
ENV UV_COMPILE_BYTECODE=1
# Copy from the cache instead of linking since it's a mounted volume
ENV UV_LINK_MODE=copy

# Install dependencies
RUN --mount=type=cache,target=/root/.cache/uv,Z \
    --mount=type=bind,source=uv.lock,target=uv.lock,Z \
    --mount=type=bind,source=pyproject.toml,target=pyproject.toml,Z \
    uv sync --frozen --no-dev --no-install-project

# Copy the current directory contents into the container at /app and install
COPY . /app
RUN --mount=type=cache,target=/root/.cache/uv,Z \
    uv sync --frozen --no-dev

# Set environment variables
ENV PATH="/app/.venv/bin:$PATH"
ENV PYTHONPATH="/plugins"

# Reset the entrypoint, don't invoke `uv`
ENTRYPOINT []

# Run the application
CMD ["python", "assembly_curator/main_flask.py", "--address=0.0.0.0", "--samples_dir=/data", "--plugin_dir=/plugins", "--n_workers=8"]

# podman build . --tag docker.io/troder/assembly-curator:0.0.1-alpha
# podman run -it --rm -v ./data-pb-share:/data:Z -v ./plugins-share:/plugins:Z -p 8080:8080 --name assembly-curator docker.io/troder/assembly-curator:0.0.1-alpha
