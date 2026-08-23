# syntax=docker/dockerfile:1.7
FROM debian:bookworm-slim AS build

ARG MINIMAP2_VERSION=2.31
ARG MINIMAP2_TAG=v2.31
ARG MINIMAP2_COMMIT=3c28777
ARG MINIMAP2_SOURCE_URL=https://github.com/lh3/minimap2/releases/download/v2.31/minimap2-2.31_x64-linux.tar.bz2
ARG MINIMAP2_SHA256
ARG BLAST_VERSION=2.17.0
ARG BLAST_SOURCE_URL=https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.17.0/ncbi-blast-2.17.0+-x64-linux.tar.gz
ARG BLAST_SHA256
ARG BLAST_MD5

RUN test -n "$MINIMAP2_SHA256" || (echo "MINIMAP2_SHA256 is required and must be independently verified" >&2; exit 64) \
    && test -n "$BLAST_SHA256" || (echo "BLAST_SHA256 is required and must be independently verified" >&2; exit 64) \
    && test -n "$BLAST_MD5" || (echo "BLAST_MD5 is required and must be independently verified" >&2; exit 64)

RUN apt-get update \
    && apt-get install --no-install-recommends -y ca-certificates curl bzip2 zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build
RUN curl --fail --location --silent --show-error "$MINIMAP2_SOURCE_URL" --output minimap2.tar.bz2 \
    && printf '%s  %s\n' "$MINIMAP2_SHA256" minimap2.tar.bz2 | sha256sum --check --strict - \
    && mkdir minimap2 \
    && tar --extract --bzip2 --file minimap2.tar.bz2 --directory minimap2 \
    && minimap2_bin="$(find minimap2 -type f -name minimap2 -perm -u+x -print -quit)" \
    && test -n "$minimap2_bin" \
    && install -D -m 0755 "$minimap2_bin" "/opt/comparators/minimap2-$MINIMAP2_VERSION"

RUN curl --fail --location --silent --show-error "$BLAST_SOURCE_URL" --output blast.tar.gz \
    && printf '%s  %s\n' "$BLAST_SHA256" blast.tar.gz | sha256sum --check --strict - \
    && printf '%s  %s\n' "$BLAST_MD5" blast.tar.gz | md5sum --check --strict - \
    && mkdir blast \
    && tar --extract --gzip --file blast.tar.gz --directory blast --strip-components=1 \
    && install -D -m 0755 blast/bin/blastn "/opt/comparators/blastn-$BLAST_VERSION"

FROM debian:bookworm-slim

ARG MINIMAP2_VERSION=2.31
ARG MINIMAP2_TAG=v2.31
ARG MINIMAP2_COMMIT=3c28777
ARG BLAST_VERSION=2.17.0

LABEL org.opencontainers.image.description="Pinned BLAST+ and minimap2 comparator environment" \
      org.jam-rs.minimap2.version="$MINIMAP2_VERSION" \
      org.jam-rs.minimap2.tag="$MINIMAP2_TAG" \
      org.jam-rs.minimap2.commit="$MINIMAP2_COMMIT" \
      org.jam-rs.blast.version="$BLAST_VERSION"

RUN apt-get update \
    && apt-get install --no-install-recommends -y ca-certificates libgomp1 zlib1g \
    && rm -rf /var/lib/apt/lists/* \
    && groupadd --system --gid 10001 comparator \
    && useradd --system --uid 10001 --gid comparator --create-home --home-dir /home/comparator comparator

COPY --from=build "/opt/comparators/minimap2-$MINIMAP2_VERSION" /usr/local/bin/minimap2
COPY --from=build "/opt/comparators/blastn-$BLAST_VERSION" /usr/local/bin/blastn

USER comparator
WORKDIR /work
ENTRYPOINT ["/bin/sh"]
