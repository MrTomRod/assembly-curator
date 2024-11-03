FROM ghcr.io/astral-sh/uv:python3.12-bookworm-slim
# Inspiration: https://github.com/astral-sh/uv-docker-example/blob/main/Dockerfile

# Run QT without display
ENV QT_QPA_PLATFORM=offscreen
ENV XDG_RUNTIME_DIR=/tmp/runtime-root

# Install curl to download rustup and libqt5svg5 for gfaviz
RUN apt-get update && \
    apt-get install --no-install-recommends -y libqt5svg5 && \
    apt-get clean

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

# Place executables in the environment at the front of the path
ENV PATH="/app/.venv/bin:/app/bin:$PATH"

# Add /plugins to PYTHONPATH=/plugins
ENV PYTHONPATH="/plugins"

# Reset the entrypoint, don't invoke `uv`
ENTRYPOINT []

CMD ["python", "assembly_curator/main_flask.py", "--address=0.0.0.0", "--samples_dir=/data", "--plugin_dir=/plugins", "--n_workers=8"]

# podman build . --tag docker.io/troder/assembly-curator:0.0.0-alpha
# podman run -it --rm -v ./data-pb-all:/data:Z -v ./plugins-pb:/plugins:Z -p 8080:8080 --name assembly-curator assembly-curator
